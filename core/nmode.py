"""
N-mode (toroidal mode number) compute core for PRISM.

Qt-free numerical pipeline for Mirnov-coil toroidal mode-number analysis,
extracted verbatim from ui/tabs/diagnostics/spectral/nmode_spectrum_tab.py so
it can be imported and run headlessly (batch / CLI) without pulling in PySide6.

Public pipeline:
    mirnov = load_mirnov_data(shot, tmin, tmax)          # MDS+ fetch + time align
    fft    = calculate_fft(mirnov, t_interval, ...)      # windowed FFT
    mode, freq, amp_evolution = calculate_mode_numbers(fft, fmin, fmax, ...)

The only behavioural change vs. the GUI version is that the MDS+ server is now a
parameter (`mds_server`) instead of the hardcoded NModeConfig.MDS_SERVER, with the
same value as the default so existing callers are unaffected.
"""

import json
import os
import getpass
from datetime import datetime
from pathlib import Path
import time as tclock
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.signal import detrend
from scipy.fft import fft
from scipy.ndimage import uniform_filter
# NOTE: MDSplus is imported lazily inside _load_single_channel so this module can be
# imported (and its FFT / mode-number math used) on a host without MDSplus.

# Load Mirnov coil configuration (core/ is one level below the version root)
_MIRNOV_CONFIG_PATH = Path(__file__).parents[1] / 'config' / 'mirnov_config.json'
with open(_MIRNOV_CONFIG_PATH, 'r') as f:
    _MIRNOV_CONFIG = json.load(f)


# =============================================================================
# Configuration
# =============================================================================

class NModeConfig:
    """N-Mode spectrum configuration"""
    MDS_SERVER = 'mdsr.kstar.kfe.re.kr:8005'

    # Local raw-Mirnov archive (PRISM NAS, mounted at /PRISM on nkstar & ukstar).
    # If a shot's archive exists here it is read instead of MDS+.
    ARCHIVE_DIR = '/PRISM/mirnov_archive'

    DEFAULT_SHOT = 0
    DEFAULT_TMIN = 0.0
    DEFAULT_TMAX = 10.0
    DEFAULT_FMIN = 0
    DEFAULT_FMAX = 100
    DEFAULT_TINTERVAL = 0.01
    DEFAULT_NMODES = 5
    DEFAULT_TOL = 0.8
    DEFAULT_FRAC = 1e-2
    DEFAULT_INTEGRATE = False
    DEFAULT_DETREND = True
    DEFAULT_MSIGN = 1
    DEFAULT_PLOT_TYPE = 'contour'
    DEFAULT_NUMC = 50

    FIGURE_SIZE = (10, 6)


# =============================================================================
# Data Structures
# =============================================================================

class MirnovData:
    """Container for Mirnov coil data"""
    def __init__(self):
        self.time = None
        self.signals = None
        self.positions = None
        self.channel_names = None
        self.shot = None


class FFTResult:
    """Container for FFT results"""
    def __init__(self):
        self.time = None
        self.freq = None
        self.amp = None
        self.phase = None
        self.positions = None
        self.shot = None
        self.array = 'TOROIDAL'
        self.p = None
        self.dt = None


# =============================================================================
# Mirnov Coil Configuration
# =============================================================================

def get_shot_year(shot):
    """Get year from shot number using config file"""
    for year_str, shot_range in _MIRNOV_CONFIG['shot_year_ranges'].items():
        if year_str.startswith('_'):
            continue
        if not isinstance(shot_range, list) or len(shot_range) != 2:
            continue
        shot_min, shot_max = shot_range
        if shot_max is None:
            shot_max = float('inf')
        if shot_min <= shot <= shot_max:
            return int(year_str)
    return 2025  # Default to latest year


def get_mirnov_config(year):
    """Get Mirnov coil configuration for given year from config file"""
    year_str = str(year)
    yearly_signs = _MIRNOV_CONFIG['yearly_signs']

    # Get signs for the year, fallback to default
    signs = yearly_signs.get(year_str, yearly_signs['default'])

    # Get base channel names and angles
    names = _MIRNOV_CONFIG['channel_names'].copy()
    angles = _MIRNOV_CONFIG['toroidal_angles']

    # Apply channel name overrides for 2017+
    if year >= 2017:
        names[9] = 'MC1P03'

    # Filter to valid channels only (signs != 0)
    valid_idx = [i for i, s in enumerate(signs) if s != 0]

    return {
        'names': [names[i] for i in valid_idx],
        'angles': [angles[i] for i in valid_idx],
        'signs': [signs[i] for i in valid_idx]
    }


# =============================================================================
# Data Loading
# =============================================================================

def _load_single_channel(args):
    """Load single Mirnov channel"""
    from MDSplus import Connection   # lazy: keep core.nmode importable without MDSplus
    shot, name, angle, sign, tmin, tmax, mds_server = args

    try:
        mds = Connection(mds_server)
        mds.openTree('kstar', shot)
        mds.get(f'SetTimeContext({tmin},{tmax},)').data()
        data = mds.get(f'\\{name}').data()
        time_arr = mds.get(f'dim_of(\\{name})').data()
        mds.get('SetTimeContext(,,)').data()
        mds.closeTree('kstar', shot)

        return {
            'name': name, 'angle': angle, 'sign': sign,
            'data': np.array(data, dtype=np.float32) * sign,
            'time': np.array(time_arr, dtype=np.float32),
            'n_points': len(data), 'error': None
        }
    except Exception as e:
        return {
            'name': name, 'angle': angle, 'sign': sign,
            'data': None, 'time': None, 'n_points': 0, 'error': str(e)
        }


# =============================================================================
# Raw-Mirnov NAS archive (HDF5)
# =============================================================================

def _archive_path(shot, archive_dir):
    """Path of a shot's archive file inside archive_dir."""
    return os.path.join(archive_dir, f"mirnov_{shot}.h5")


def _compression_kwargs(compression):
    """h5py create_dataset compression kwargs. Default is blosc-zstd (fastest to
    write AND smallest: ~1s / ~760MB per full shot) via hdf5plugin. If hdf5plugin
    is not installed it falls back to lzf (built into h5py, ~2s / ~1140MB) so a
    write never fails. gzip (~37s / ~940MB) and none are also selectable."""
    if compression in (None, 'none'):
        return {}
    if compression == 'gzip':
        return {'compression': 'gzip', 'compression_opts': 4}
    if compression == 'lzf':
        return {'compression': 'lzf'}
    # default / 'zstd' / 'blosc': blosc-zstd, falling back to lzf without the plugin
    try:
        import hdf5plugin
        return dict(hdf5plugin.Blosc(cname='zstd', clevel=5,
                                     shuffle=hdf5plugin.Blosc.BITSHUFFLE))
    except Exception:
        print("[n-Mode] hdf5plugin unavailable; using lzf instead of blosc-zstd")
        return {'compression': 'lzf'}


def _chmod_quiet(path, mode):
    """Best-effort chmod so a shared NAS archive stays writable by all users.
    Failures (e.g. not the owner) are ignored."""
    try:
        os.chmod(path, mode)
    except OSError:
        pass


def _ensure_archive_dir(archive_dir):
    """Create the archive dir and make it world-writable so any user can add /
    refresh shots in the shared NAS archive."""
    os.makedirs(archive_dir, exist_ok=True)
    _chmod_quiet(archive_dir, 0o777)


def _set_archive_attrs(f, shot, year, mds_server):
    """Write file-level provenance attributes."""
    f.attrs['shot'] = shot
    f.attrs['year'] = year
    f.attrs['mds_server'] = mds_server
    f.attrs['units'] = 'physical = int16 * scale (polarity applied at load)'
    f.attrs['created'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    f.attrs['created_by'] = getpass.getuser()


def _write_channel(f, name, raw, t, comp):
    """Write one raw (polarity-free) channel as int16 + scale/timebase attrs."""
    raw = np.asarray(raw)
    n = raw.size
    t0 = float(t[0])
    dt = float((t[100] - t[0]) / 100.0) if n > 100 else float((t[-1] - t[0]) / max(n - 1, 1))
    m = float(np.max(np.abs(raw))) or 1.0
    scale = m / 32767.0
    q = np.round(raw / scale).astype(np.int16)
    ds = f.create_dataset(name, data=q, chunks=(min(1_000_000, n),), **comp)
    ds.attrs.update({'scale': scale, 't0': t0, 'dt': dt, 'n': n})


def _auto_archive_from_results(shot, results, archive_dir, mds_server, compression='zstd'):
    """Archive an MDS+ download to the NAS so later loads skip MDS+. Best-effort:
    any failure (no h5py, unwritable NAS, partial data) is logged and ignored so
    it never breaks the analysis. Skips if an archive already exists.

    `results` come from the MDS+ path where data = raw*sign; the stored value is
    raw (sign removed) so the read path can re-apply the current config polarity.
    """
    if not archive_dir:
        return
    path = _archive_path(shot, archive_dir)
    if os.path.isfile(path):
        return
    try:
        import h5py
    except Exception:
        return
    tmp = path + '.tmp'
    try:
        _ensure_archive_dir(archive_dir)
        year = get_shot_year(shot)
        comp = _compression_kwargs(compression)
        with h5py.File(tmp, 'w') as f:
            _set_archive_attrs(f, shot, year, mds_server)
            for r in results:
                if r.get('error') or r.get('data') is None:
                    continue
                raw = np.asarray(r['data']) * np.float32(r['sign'])   # data was raw*sign
                _write_channel(f, r['name'], raw, r['time'], comp)
        os.replace(tmp, path)
        _chmod_quiet(path, 0o666)   # let other users overwrite/refresh
        size_mb = os.path.getsize(path) / 1024 / 1024
        print(f"[n-Mode] Auto-archived shot #{shot} -> {path} ({size_mb:.0f} MB)")
    except Exception as e:
        print(f"[n-Mode] Auto-archive skipped ({e})")
        try:
            if os.path.exists(tmp):
                os.remove(tmp)
        except Exception:
            pass


def archive_info(shot, archive_dir=None):
    """Return archive provenance for a shot ({created, created_by, shot, year,
    n_channels, size_mb, path}) or {} if not archived. Read-only, no data load."""
    if archive_dir is None:
        archive_dir = NModeConfig.ARCHIVE_DIR
    path = _archive_path(shot, archive_dir)
    if not os.path.isfile(path):
        return {}
    try:
        import h5py
        with h5py.File(path, 'r') as f:
            info = {k: f.attrs.get(k) for k in ('shot', 'year', 'created', 'created_by', 'mds_server')}
            info['n_channels'] = len(list(f.keys()))
    except Exception:
        return {}
    info['path'] = path
    info['size_mb'] = round(os.path.getsize(path) / 1024 / 1024, 1)
    return info


def _read_mirnov_archive(shot, tmin, tmax, archive_dir, config):
    """Read the [tmin, tmax] window for a shot from its HDF5 archive.

    Returns a results list compatible with the MDS+ loader (one dict per
    channel), or None if no usable archive is present. Stored data is raw
    int16*scale (no polarity); the current config sign/angle are applied here,
    exactly as the MDS+ path does.
    """
    if not archive_dir:
        return None
    path = _archive_path(shot, archive_dir)
    if not os.path.isfile(path):
        return None
    try:
        import h5py
    except Exception as e:
        print(f"[n-Mode] Archive present but h5py unavailable ({e}); using MDS+")
        return None
    # Register the blosc-zstd filter so zstd-compressed archives can be read.
    # Harmless (and unneeded) for lzf/gzip/uncompressed files if it's absent.
    try:
        import hdf5plugin  # noqa: F401
    except Exception:
        pass

    results = []
    try:
        with h5py.File(path, 'r') as f:
            # Normalize older ISO timestamps (2026-06-29T16:32:50) for display
            created = str(f.attrs.get('created', 'unknown')).replace('T', ' ')
            created_by = f.attrs.get('created_by', 'unknown')
            print(f"[n-Mode] Loaded shot #{shot} from archive: {path}")
            print(f"[n-Mode]   archived {created} by {created_by}")
            for name, angle, sign in zip(config['names'], config['angles'], config['signs']):
                if name not in f:
                    continue
                ds = f[name]
                t0 = float(ds.attrs['t0']); dt = float(ds.attrs['dt'])
                n = int(ds.attrs['n']); scale = float(ds.attrs['scale'])
                i0 = 0 if tmin is None else max(0, int(np.floor((tmin - t0) / dt)))
                i1 = n if tmax is None else min(n, int(np.ceil((tmax - t0) / dt)) + 1)
                if i1 <= i0:               # requested window outside data -> take all
                    i0, i1 = 0, n
                raw = ds[i0:i1].astype(np.float32)
                data = raw * np.float32(scale) * np.float32(sign)
                time_arr = (t0 + dt * np.arange(i0, i1)).astype(np.float32)
                results.append({'name': name, 'angle': angle, 'sign': sign,
                                'data': data, 'time': time_arr,
                                'n_points': len(data), 'error': None})
    except Exception as e:
        print(f"[n-Mode] Failed to read archive {path} ({e}); falling back to MDS+")
        return None
    return results if results else None


def archive_mirnov(shot, *, mds_server=None, archive_dir=None, compression='zstd', overwrite=False):
    """Download a full shot's toroidal Mirnov channels from MDS+ and archive them
    to archive_dir as one HDF5 file (int16 + per-channel scale, blosc-zstd by default).

    One-time operation; afterwards load_mirnov_data reads the shot from the
    archive instead of MDS+. If the archive already exists it is kept and the
    download is skipped unless overwrite=True (use that to re-download a shot).
    Returns the archive file path.
    """
    import h5py
    from MDSplus import Connection

    if mds_server is None:
        mds_server = NModeConfig.MDS_SERVER
    if archive_dir is None:
        archive_dir = NModeConfig.ARCHIVE_DIR
    _ensure_archive_dir(archive_dir)

    year = get_shot_year(shot)
    config = get_mirnov_config(year)
    names, angles, signs = config['names'], config['angles'], config['signs']
    path = _archive_path(shot, archive_dir)
    tmp = path + '.tmp'
    comp = _compression_kwargs(compression)

    if os.path.isfile(path) and not overwrite:
        info = archive_info(shot, archive_dir)
        print(f"[Archive] shot #{shot} already archived (archived {info.get('created','?')} "
              f"by {info.get('created_by','?')}); skipping. Use overwrite=True to re-download.")
        return path

    print(f"[Archive] shot #{shot} (year {year}, {len(names)} ch) from MDS+ -> {path}")
    t_start = tclock.time()
    mds = Connection(mds_server)
    mds.openTree('kstar', shot)
    try:
        with h5py.File(tmp, 'w') as f:
            _set_archive_attrs(f, shot, year, mds_server)
            for i, name in enumerate(names, 1):
                raw = np.asarray(mds.get('\\' + name).data())
                t = np.asarray(mds.get(f'dim_of(\\{name})').data())
                _write_channel(f, name, raw, t, comp)
                print(f"[Archive]   {i:>2d}/{len(names)} {name}  n={raw.size:,}", flush=True)
    finally:
        mds.closeTree('kstar', shot)
    os.replace(tmp, path)
    _chmod_quiet(path, 0o666)   # let other users overwrite/refresh
    size_mb = os.path.getsize(path) / 1024 / 1024
    print(f"[Archive] done: {path}  ({size_mb:.0f} MB, {tclock.time() - t_start:.0f}s)")
    return path


def load_mirnov_data(shot, tmin, tmax, *, mds_server=None, archive_dir=None, auto_archive=True):
    """Load Mirnov coil data, preferring the local NAS archive over MDS+.

    If an HDF5 archive for this shot exists under archive_dir (default
    NModeConfig.ARCHIVE_DIR, the PRISM NAS), the requested [tmin, tmax] window is
    read from it (partial read, no MDS+ traffic). Otherwise the data is fetched
    from MDS+ in parallel and, when auto_archive is True, written to the archive
    so subsequent loads of this shot skip MDS+. This applies equally to GUI and
    batch/API runs. Set auto_archive=False to disable the automatic save.
    """
    if mds_server is None:
        mds_server = NModeConfig.MDS_SERVER
    if archive_dir is None:
        archive_dir = NModeConfig.ARCHIVE_DIR

    year = get_shot_year(shot)
    config = get_mirnov_config(year)
    n_channels = len(config['names'])

    # --- Prefer the local archive (it prints its own provenance line) ---
    results = _read_mirnov_archive(shot, tmin, tmax, archive_dir, config)
    if results is None:
        # Explain why we are going to MDS+ (not in archive / unreadable)
        apath = _archive_path(shot, archive_dir) if archive_dir else None
        if apath and not os.path.isfile(apath):
            note = " (will auto-archive)" if auto_archive else ""
            print(f"[n-Mode] Shot #{shot} not in archive ({apath}); fetching from MDS+{note}")
        else:
            print(f"[n-Mode] Archive for shot #{shot} unavailable; fetching from MDS+")
        print(f"[n-Mode] Loading Mirnov data for shot #{shot} (year {year}) from MDS+ ...")
        print(f"[n-Mode]   Time range: {tmin:.2f} - {tmax:.2f} s")
        print(f"[n-Mode]   Channels ({n_channels}):")
        for i, (name, angle) in enumerate(zip(config['names'], config['angles']), 1):
            print(f"[n-Mode]     {i:>2d}. {name}  {angle:>6.2f}°")
        print(f"[n-Mode]   Threads: {n_channels}")

        args_list = [
            (shot, name, angle, sign, tmin, tmax, mds_server)
            for name, angle, sign in zip(config['names'], config['angles'], config['signs'])
        ]

        results = []
        completed = 0
        t0 = tclock.time()

        with ThreadPoolExecutor(max_workers=n_channels) as executor:
            futures = {executor.submit(_load_single_channel, args): args[1] for args in args_list}
            for future in as_completed(futures):
                results.append(future.result())
                completed += 1
                progress = completed / n_channels
                bar = '█' * int(40 * progress) + '░' * (40 - int(40 * progress))
                print(f'\r[n-Mode]   [{bar}] {completed}/{n_channels}', end='', flush=True)

        print(f'\n[n-Mode]   Completed in {tclock.time() - t0:.2f} sec')

        # Auto-archive the freshly downloaded shot so later loads skip MDS+
        if auto_archive:
            _auto_archive_from_results(shot, results, archive_dir, mds_server)

    valid_results = []
    results.sort(key=lambda x: x['angle'])

    # Collect valid results
    for result in results:
        if result['error'] is None and result['n_points'] > 0:
            valid_results.append(result)

    if len(valid_results) == 0:
        raise RuntimeError("No valid Mirnov channels found")

    # Find common time range across all channels
    common_tmin = max(r['time'][0] for r in valid_results)
    common_tmax = min(r['time'][-1] for r in valid_results)

    if common_tmin >= common_tmax:
        raise RuntimeError(f"No overlapping time range found across channels")

    print(f"[n-Mode]   Common time range: {common_tmin:.4f} - {common_tmax:.4f} s")

    # Use the channel with finest sampling as reference
    ref_result = min(valid_results, key=lambda r: (r['time'][-1] - r['time'][0]) / r['n_points'])
    ref_time = ref_result['time']

    # Trim reference time to common range
    ref_mask = (ref_time >= common_tmin) & (ref_time <= common_tmax)
    time_arr = ref_time[ref_mask]
    n_points = len(time_arr)

    print(f"[n-Mode]   Reference channel: {ref_result['name']} ({n_points} points)")

    # Interpolate all channels to common time grid
    signals, valid_names, valid_angles = [], [], []

    for result in valid_results:
        ch_time = result['time']
        ch_data = result['data']

        # Check if interpolation is needed
        if len(ch_time) == n_points and np.allclose(ch_time[ref_mask] if len(ch_time) == len(ref_time) else ch_time, time_arr, rtol=1e-6):
            # Same time grid, just trim
            ch_mask = (ch_time >= common_tmin) & (ch_time <= common_tmax)
            interp_data = ch_data[ch_mask]
        else:
            # Interpolate to common time grid
            interp_data = np.interp(time_arr, ch_time, ch_data).astype(np.float32)

        signals.append(interp_data)
        valid_names.append(result['name'])
        valid_angles.append(result['angle'])

    print(f"[n-Mode]   Aligned {len(signals)} channels to common time grid")

    mirnov = MirnovData()
    mirnov.time = np.array(time_arr, dtype=np.float32)
    mirnov.signals = np.array(signals, dtype=np.float32)
    mirnov.positions = np.array(valid_angles)
    mirnov.channel_names = valid_names
    mirnov.shot = shot

    print(f"[n-Mode]   Loaded {len(valid_names)} channels, {len(time_arr)} points")
    return mirnov


# =============================================================================
# FFT Calculation
# =============================================================================

def calculate_fft(mirnov, t_interval, integrate=False, run_detrend=True, tmin=None, tmax=None):
    """Calculate FFT for all channels"""
    t0 = tclock.time()

    time = mirnov.time
    signals = mirnov.signals.astype(np.float32)
    n_channels = signals.shape[0]

    if tmin is not None and tmax is not None:
        idx_start = np.searchsorted(time, tmin)
        idx_end = np.searchsorted(time, tmax)
        time = time[idx_start:idx_end]
        signals = signals[:, idx_start:idx_end]

    n_points = signals.shape[1]
    dt = (time[100] - time[0]) / 100.0
    fs = 1.0 / dt

    # Window the record into q intervals of p samples. p is fixed to the
    # requested t_interval (p = int(t_interval/dt)) -- NOT n_points/q -- so the
    # window length, and hence the frequency resolution (1/(p*dt)), always equals
    # the requested t_interval regardless of record length, matching MCspectrogram.
    # Trailing samples beyond q*p (< one window) are dropped, as in the reference.
    q = int(n_points * dt / t_interval)
    p = int(t_interval / dt)
    n_windows, points_per_window = q, p

    print(f"[n-Mode] FFT: fs={fs/1e6:.2f}MHz, windows={n_windows}, pts/win={points_per_window}")

    # Step 1: Preprocessing
    print('[n-Mode]   1/4 Preprocessing...')
    if integrate:
        for i in range(n_channels):
            signals[i] = np.cumsum(signals[i]) * dt * 1e4

    n_use = n_windows * points_per_window
    signals = signals[:, :n_use]
    time_trimmed = time[:n_use]
    signals_reshaped = signals.reshape(n_channels, n_windows, points_per_window)

    # Step 2: Detrend
    print('[n-Mode]   2/4 Detrending...')
    if run_detrend:
        n_segments = 8
        denom = points_per_window % n_segments
        main_size = points_per_window - denom
        seg_size = main_size // n_segments

        main_part = signals_reshaped[:, :, :main_size].reshape(
            n_channels, n_windows, n_segments, seg_size)
        main_part = detrend(main_part, axis=3, type='linear')
        signals_reshaped[:, :, :main_size] = main_part.reshape(
            n_channels, n_windows, main_size)

        # Handle remaining points (denom)
        if denom > 0:
            remainder = signals_reshaped[:, :, -denom:].reshape(
                n_channels, n_windows, denom, 1)
            remainder = detrend(remainder, axis=2, type='linear')
            signals_reshaped[:, :, -denom:] = remainder.reshape(
                n_channels, n_windows, denom)

    # Step 3: FFT calculation
    print('[n-Mode]   3/4 FFT...')
    n_freq_half = points_per_window // 2 + 1
    overp = 1.0 / float(points_per_window)

    amp = np.zeros((n_windows, n_freq_half, n_channels), dtype=np.float32)
    phase = np.zeros((n_windows, n_freq_half, n_channels), dtype=np.float32)

    fft_result_raw = fft(signals_reshaped, axis=2, workers=-1)
    fft_half = fft_result_raw[:, :, :n_freq_half]
    amp = np.transpose(np.abs(fft_half) * overp, (1, 2, 0)).astype(np.float32)
    phase = np.transpose(np.angle(fft_half), (1, 2, 0)).astype(np.float32)

    # Step 4: Post-processing
    print('[n-Mode]   4/4 Finalizing...')
    window_start_indices = np.arange(n_windows) * points_per_window
    time_centers = (time_trimmed[window_start_indices] +
                   time_trimmed[window_start_indices + points_per_window - 1]) / 2.0

    idt = 1.0 / (points_per_window * dt)
    freq = np.arange(n_freq_half, dtype='f') * idt * 1e-3

    result = FFTResult()
    result.time = time_centers.astype(np.float32)
    result.freq = freq
    result.amp = amp
    result.phase = phase
    result.positions = mirnov.positions
    result.shot = mirnov.shot
    result.p = points_per_window
    result.dt = dt

    print(f'\n[n-Mode]   Completed in {tclock.time() - t0:.2f} sec')
    return result


# =============================================================================
# Mode Number Calculation
# =============================================================================

def calculate_mode_numbers(fft_result, fmin, fmax, tol, nmodes, frac, msign,
                           integrate=False, amp_mode='max'):
    """Calculate toroidal mode numbers

    amp_mode: how to reduce per-mode amplitude across the [fmin, fmax] band into
              the amplitude-evolution trace -- 'max' (peak bin, default) or
              'sum' (band sum, matching MCspectrogram's active code path).
    """
    t0 = tclock.time()

    amp = fft_result.amp
    phase = fft_result.phase
    freq = fft_result.freq
    positions = fft_result.positions

    n_time, freq_range, n_channels = amp.shape
    nch = 1

    fmax_threshold = frac * np.max(amp[:, :, nch])
    mode = np.zeros((n_time, freq_range), dtype='int')

    wh = np.where(amp[:, :, nch] >= fmax_threshold)
    wh_i, wh_j = wh[0], wh[1]
    valid = wh_j != 0
    wh_i, wh_j = wh_i[valid], wh_j[valid]

    if len(wh_i) == 0:
        return mode, freq, np.zeros((2 * nmodes, n_time))

    pos_rad = np.array(positions) * np.pi / 180.0
    nmeas = n_channels
    nvec = np.arange(nmodes) + 1.0

    pharr = nvec * pos_rad[:, np.newaxis] * 1j
    exparr = np.exp(pharr) / float(nmeas)

    bamp = amp[wh_i, wh_j, :]
    barg = phase[wh_i, wh_j, :]
    bvec = bamp * np.exp(1j * barg)
    b2avg = np.sum(bamp ** 2, axis=1) / nmeas

    bmodep = np.dot(bvec, exparr)
    bmoden = np.dot(bvec, np.conj(exparr))

    if bmodep.ndim == 1:
        bmodep = bmodep.reshape(1, -1)
        bmoden = bmoden.reshape(1, -1)
        b2avg = np.array([b2avg]) if np.isscalar(b2avg) else b2avg

    b2avg_vec = np.tile(b2avg[:, np.newaxis], (1, nmodes))

    std_pos = np.sqrt(np.abs(1.0 - np.abs(bmodep) ** 2 / b2avg_vec))
    std_neg = np.sqrt(np.abs(1.0 - np.abs(bmoden) ** 2 / b2avg_vec))

    kmn_pos = std_pos.argmin(axis=1)
    mn_pos = std_pos[np.arange(len(kmn_pos)), kmn_pos]
    kmn_neg = std_neg.argmin(axis=1)
    mn_neg = std_neg[np.arange(len(kmn_neg)), kmn_neg]

    # Assign every bin within tolerance to its best positive mode, then let the
    # negative assignment overwrite any bin that is also within tolerance for a
    # negative mode. This reproduces MCspectrogram's behaviour exactly (its
    # `indices22 = indices20 and indices21` collapses to just `mn <= tol`, so the
    # negative pass overwrites the positive one). NOTE: this is intentionally the
    # reference's logic, not a min-std mutually-exclusive choice, so PRISM's
    # mode map / amplitude traces match MCspectrogram.
    pos_assign = mn_pos <= tol
    neg_assign = mn_neg <= tol

    mode[wh_i[pos_assign], wh_j[pos_assign]] = kmn_pos[pos_assign] + 1
    mode[wh_i[neg_assign], wh_j[neg_assign]] = -kmn_neg[neg_assign] - 1

    amp_evolution = compute_amp_evolution(
        amp, freq, mode, fmin, fmax, nmodes, integrate, amp_mode, nch
    )

    print(f"[n-Mode] Mode calculation completed in {tclock.time() - t0:.2f} sec")
    return mode, freq, amp_evolution


def compute_amp_evolution(amp, freq, mode, fmin, fmax, nmodes,
                          integrate=False, amp_mode='max', nch=1):
    """Build the per-mode amplitude-evolution trace from an existing mode map.

    Separated from calculate_mode_numbers so the GUI can switch between 'max'
    and 'sum' reductions (and re-plot) without re-running the FFT / mode fit.

    amp:      FFT amplitude array (n_time, n_freq, n_channels)
    mode:     toroidal mode-number map (n_time, n_freq)
    amp_mode: 'max' -> peak amplitude bin in [fmin, fmax] (default, raw)
              'sum' -> sum of amplitude over [fmin, fmax], then a size-3 moving
                       average over time (matches MCspectrogram's active path)
    nch:      representative coil channel index used for the amplitude trace
              (1, matching the threshold channel and MCspectrogram)
    """
    n_time = amp.shape[0]

    amp_for_evolution = np.abs(amp[:, :, nch]).copy()
    if not integrate:
        amp_for_evolution = 2.0 * amp_for_evolution
        freq_hz = freq * 1e3
        freq_hz_safe = np.where(freq_hz > 0, freq_hz, 1.0)
        amp_for_evolution = amp_for_evolution / (2.0 * np.pi * freq_hz_safe)
        amp_for_evolution[:, freq_hz == 0] = 0.0
        amp_for_evolution = amp_for_evolution * 1e4

    freq_idx = np.where((freq >= fmin) & (freq <= fmax))[0]
    amp_evolution = np.zeros((2 * nmodes, n_time))
    temp = np.zeros_like(amp_for_evolution)

    for i in range(nmodes):
        for sign_mult, offset in [(1, 0), (-1, nmodes)]:
            ind1 = np.where(mode == sign_mult * (i + 1))
            if len(ind1[0]) > 0:
                temp[ind1] = amp_for_evolution[ind1]
                if amp_mode == 'sum':
                    trace = temp[:, freq_idx].sum(axis=1)
                    amp_evolution[i + offset, :] = uniform_filter(trace, size=3)
                else:
                    amp_evolution[i + offset, :] = temp[:, freq_idx].max(axis=1)
                temp[ind1] = 0.0

    return amp_evolution
