"""
Typed result containers for the PRISM batch SDK, saved as .npz (one file per job).

The .npz layout matches what the GUI's n-mode tab already writes, so existing
NPZ-reading scripts (and the tab's bundled example reader) keep working unchanged:

    metadata       : 0-d object array holding the provenance dict
                     (load with allow_pickle=True, then .item())
    time           : window-center times [s]
    frequency      : frequency axis [kHz]
    mode_spectrum  : signed toroidal mode number, shape (n_time, n_freq)
    amplitude      : per-mode amplitude evolution (= amplitude_max), kept for
                     backward compatibility, shape (n_mode_rows, n_time)
    amplitude_max  : peak-bin amplitude evolution, shape (n_mode_rows, n_time)
    amplitude_sum  : band-sum amplitude evolution, shape (n_mode_rows, n_time)
"""

from dataclasses import dataclass, field

import numpy as np


def _shape(a):
    """Compact shape string for result reprs, e.g. '(244, 100, 100)' or 'None'."""
    shp = getattr(a, "shape", None)
    return str(tuple(shp)) if shp is not None else str(a)


@dataclass
class NModeResult:
    """N-mode spectrum result.

    Attributes:
        shot:          KSTAR shot number
        time:          window-center times [s], shape (n_time,)
        frequency:     frequency axis [kHz], shape (n_freq,)
        mode_spectrum: signed toroidal mode number, shape (n_time, n_freq)
        amplitude:     per-mode amplitude evolution (= amplitude_max), shape (n_mode_rows, n_time)
        amplitude_max: peak-bin amplitude evolution, shape (n_mode_rows, n_time)
        amplitude_sum: band-sum amplitude evolution, shape (n_mode_rows, n_time)
        meta:          provenance dict (shot, params, version, created_iso, ...)
    """

    shot: int
    time: np.ndarray
    frequency: np.ndarray
    mode_spectrum: np.ndarray
    amplitude: np.ndarray
    amplitude_max: np.ndarray = None
    amplitude_sum: np.ndarray = None
    meta: dict = field(default_factory=dict)

    SUBSYSTEM = "nmode"

    def __post_init__(self):
        # `amplitude` is the backward-compatible alias of the Max reduction.
        if self.amplitude_max is None:
            self.amplitude_max = self.amplitude

    def __repr__(self):
        return (
            f"NModeResult(shot={self.shot})\n"
            f"  time          {_shape(self.time)}   window-center times [s]\n"
            f"  frequency     {_shape(self.frequency)}   frequency axis [kHz]\n"
            f"  mode_spectrum {_shape(self.mode_spectrum)}   signed toroidal mode n (time, freq)\n"
            f"  amplitude_max {_shape(self.amplitude_max)}   peak-bin amplitude (n_mode_rows, time) [Gauss]\n"
            f"  amplitude_sum {_shape(self.amplitude_sum)}   band-sum amplitude (n_mode_rows, time) [Gauss]\n"
            f"  meta keys: {sorted(self.meta)}"
        )

    def save_npz(self, file):
        """Write this result as .npz to `file` (a path or open binary file object).

        Layout is identical to the GUI's saved n-mode NPZ (metadata stored as a
        0-d object array), so the same reader scripts work for both. Both
        amplitude reductions are written; `amplitude` aliases amplitude_max for
        older readers.
        """
        np.savez(
            file,
            metadata=np.array(self.meta),
            time=self.time,
            frequency=self.frequency,
            mode_spectrum=self.mode_spectrum,
            amplitude=self.amplitude_max,
            amplitude_max=self.amplitude_max,
            amplitude_sum=self.amplitude_sum,
        )

    @classmethod
    def load_npz(cls, path):
        """Rebuild an NModeResult from a .npz file written by save_npz.

        Falls back to the legacy single-`amplitude` layout for older files.
        """
        data = np.load(path, allow_pickle=True)
        try:
            meta = data["metadata"].item()
            if not isinstance(meta, dict):
                meta = {}
            files = data.files
            amp_max = data["amplitude_max"] if "amplitude_max" in files else data["amplitude"]
            amp_sum = data["amplitude_sum"] if "amplitude_sum" in files else None
            return cls(
                shot=int(meta.get("shot", 0)),
                time=data["time"],
                frequency=data["frequency"],
                mode_spectrum=data["mode_spectrum"],
                amplitude=amp_max,
                amplitude_max=amp_max,
                amplitude_sum=amp_sum,
                meta=dict(meta),
            )
        finally:
            data.close()


@dataclass
class IRVBResult:
    """IRVB result. Layout matches the GUI IRVB tab's saved NPZ exactly, so the
    tab's bundled example reader works on CLI-produced files too:

        metadata       : provenance dict (shot, efit_tree, psi_boundaries, ...)
        time           : IRVB frame times [s]
        R, Z           : IRVB reconstruction grid [m]
        prad_2d        : 2D radiated power density (time, Z, R) [MW/m^3]
        psi_n          : normalized psi on the IRVB (R,Z) grid (time, Z, R)
        region_prad    : per-region Prad (n_regions, time) [MW]
        ptot           : total radiated power (time,) [MW]
        efit_*         : EFIT grid, per-frame boundary/axis/psi, and limiter
    """

    shot: int
    time: np.ndarray
    R: np.ndarray
    Z: np.ndarray
    prad_2d: np.ndarray
    psi_n: np.ndarray
    region_prad: np.ndarray
    ptot: np.ndarray
    efit_r_grid: np.ndarray
    efit_z_grid: np.ndarray
    efit_psi_n: np.ndarray
    efit_bdry_r: np.ndarray
    efit_bdry_z: np.ndarray
    efit_nbdry: np.ndarray
    efit_rmaxis: np.ndarray
    efit_zmaxis: np.ndarray
    efit_limiter_r: np.ndarray
    efit_limiter_z: np.ndarray
    meta: dict = field(default_factory=dict)

    SUBSYSTEM = "irvb"

    def __repr__(self):
        return (
            f"IRVBResult(shot={self.shot})\n"
            f"  time        {_shape(self.time)}   IRVB frame times [s]\n"
            f"  R, Z        {_shape(self.R)}, {_shape(self.Z)}   IRVB grid [m]\n"
            f"  prad_2d     {_shape(self.prad_2d)}   radiated power density (time, Z, R) [MW/m^3]\n"
            f"  psi_n       {_shape(self.psi_n)}   psi_N on the IRVB grid (mask any region)\n"
            f"  region_prad {_shape(self.region_prad)}   per-region Prad (n_regions, time) [MW]\n"
            f"  ptot        {_shape(self.ptot)}   total radiated power [MW]\n"
            f"  efit_*      grid / per-frame boundary, axis, psi, limiter\n"
            f"  meta keys: {sorted(self.meta)}"
        )

    def save_npz(self, file):
        """Write this result as .npz (identical key layout to the GUI IRVB save)."""
        np.savez(
            file,
            metadata=np.array(self.meta),
            time=self.time, R=self.R, Z=self.Z,
            prad_2d=self.prad_2d, psi_n=self.psi_n,
            region_prad=self.region_prad, ptot=self.ptot,
            efit_r_grid=self.efit_r_grid, efit_z_grid=self.efit_z_grid,
            efit_psi_n=self.efit_psi_n,
            efit_bdry_r=self.efit_bdry_r, efit_bdry_z=self.efit_bdry_z,
            efit_nbdry=self.efit_nbdry,
            efit_rmaxis=self.efit_rmaxis, efit_zmaxis=self.efit_zmaxis,
            efit_limiter_r=self.efit_limiter_r, efit_limiter_z=self.efit_limiter_z,
        )

    @classmethod
    def load_npz(cls, path):
        """Rebuild an IRVBResult from a .npz file written by save_npz."""
        data = np.load(path, allow_pickle=True)
        try:
            meta = data["metadata"].item()
            if not isinstance(meta, dict):
                meta = {}
            return cls(
                shot=int(meta.get("shot", 0)),
                time=data["time"], R=data["R"], Z=data["Z"],
                prad_2d=data["prad_2d"], psi_n=data["psi_n"],
                region_prad=data["region_prad"], ptot=data["ptot"],
                efit_r_grid=data["efit_r_grid"], efit_z_grid=data["efit_z_grid"],
                efit_psi_n=data["efit_psi_n"],
                efit_bdry_r=data["efit_bdry_r"], efit_bdry_z=data["efit_bdry_z"],
                efit_nbdry=data["efit_nbdry"],
                efit_rmaxis=data["efit_rmaxis"], efit_zmaxis=data["efit_zmaxis"],
                efit_limiter_r=data["efit_limiter_r"], efit_limiter_z=data["efit_limiter_z"],
                meta=dict(meta),
            )
        finally:
            data.close()
