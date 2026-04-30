"""
TRANSP U-File directory loader.

Scans a directory for `<prefix><shot>.<ext>` files, parses each via
core.ufile_parser, and groups them into a "run" record:

    {
        'run_key':  '40848_F_dirname',
        'shot':     40848,
        'prefix':   'F',
        'label':    'dirname',          # display label
        'dir':      '/path/to/run',
        'mtime':    1713434820.0,
        'ufiles':   { 'CUR': <uf_dict>, 'NER': <uf_dict>, ... },
        'classes':  { 'CUR': 'time_trace', 'NER': 'profile_time', ... },
    }
"""

import os
import re
import time

from core.ufile_parser import parse_ufile, classify_ufile

_FNAME_RE = re.compile(r'^([A-Za-z]+)(\d+)\.([A-Za-z0-9]{2,4})$')


def _make_run_key(shot, prefix, label):
    return f"{shot}_{prefix}_{label}"


def _make_label(dir_path, mtime):
    """Derive a human-readable label for a run directory."""
    base = os.path.basename(os.path.normpath(dir_path))
    if base and not base.isdigit():
        return base
    # Fall back to mtime
    return time.strftime("%y%m%d-%H%M", time.localtime(mtime))


def scan_run_directory(dir_path, verbose=True):
    """Scan one directory and return a run record (or None if empty)."""
    if not os.path.isdir(dir_path):
        return None

    candidates = []
    for name in sorted(os.listdir(dir_path)):
        full = os.path.join(dir_path, name)
        if not os.path.isfile(full):
            continue
        m = _FNAME_RE.match(name)
        if not m:
            continue
        candidates.append((m.group(1), int(m.group(2)), m.group(3).upper(), full))

    if not candidates:
        return None

    # Pick dominant (prefix, shot) pair (most files share it)
    pairs = {}
    for prefix, shot, ext, _ in candidates:
        pairs[(prefix, shot)] = pairs.get((prefix, shot), 0) + 1
    (prefix, shot) = max(pairs, key=pairs.get)

    ufiles = {}
    classes = {}
    skipped = []
    for fprefix, fshot, ext, full in candidates:
        if fprefix != prefix or fshot != shot:
            continue
        if ext in ufiles:
            continue  # already loaded
        try:
            uf = parse_ufile(full)
            if uf is None:
                continue
            ufiles[ext] = uf
            classes[ext] = classify_ufile(uf)
        except Exception as e:
            skipped.append((ext, str(e)))
            if verbose:
                print(f"[UFILE] Failed to parse {os.path.basename(full)}: {e}")

    if not ufiles:
        return None

    mtime = max(os.path.getmtime(p) for p in
                (os.path.join(dir_path, f"{prefix}{shot}.{e.upper()}")
                 for e in ufiles)
                if os.path.isfile(p))

    label = _make_label(dir_path, mtime)
    run_key = _make_run_key(shot, prefix, label)

    if verbose:
        n_classes = {}
        for cls in classes.values():
            n_classes[cls] = n_classes.get(cls, 0) + 1
        cls_summary = ", ".join(f"{k}={v}" for k, v in n_classes.items())
        print(f"[UFILE] {prefix}{shot} · {label}: {len(ufiles)} files "
              f"({cls_summary})")
        for ext in sorted(ufiles):
            uf = ufiles[ext]
            x0 = uf.get('X0', {})
            x1 = uf.get('X1', {})
            shape = uf.get('F', {}).get('data')
            shape = shape.shape if hasattr(shape, 'shape') else ()
            print(f"  {ext:>4}  NDIM={uf.get('NDIM')}  shape={shape}  "
                  f"X0={x0.get('name','')!r}/{x0.get('units','')!r}  "
                  f"X1={x1.get('name','')!r}/{x1.get('units','')!r}  "
                  f"-> {classes[ext]}")

    return {
        'run_key': run_key,
        'shot': shot,
        'prefix': prefix,
        'label': label,
        'dir': dir_path,
        'mtime': mtime,
        'ufiles': ufiles,
        'classes': classes,
        'skipped': skipped,
    }


def list_extensions_by_class(run, target_classes):
    """Return sorted list of extensions in `run` matching given classes."""
    if isinstance(target_classes, str):
        target_classes = (target_classes,)
    target = set(target_classes)
    return sorted([ext for ext, cls in run['classes'].items() if cls in target])


# ----- Catalog (optional pretty labels) ---------------------------------------

_CATALOG_CACHE = None


def get_catalog():
    """Load (and cache) the optional `config/transp_ufiles.json` catalog."""
    global _CATALOG_CACHE
    if _CATALOG_CACHE is not None:
        return _CATALOG_CACHE
    import json
    cfg_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'config', 'transp_ufiles.json'
    )
    try:
        with open(cfg_path, 'r') as f:
            _CATALOG_CACHE = json.load(f)
    except Exception:
        _CATALOG_CACHE = {}
    return _CATALOG_CACHE


def pretty_label(ext, uf):
    """Build a display label for a UFILE entry: '{ext} {name} [{units}]'."""
    cat = get_catalog().get(ext.upper(), {})
    name = cat.get('label') or (uf.get('F', {}).get('name') or '').strip()
    units = (uf.get('F', {}).get('units') or '').strip()
    parts = [ext]
    if name and name != ext:
        parts.append(name)
    if units:
        parts.append(f"[{units}]")
    return '  '.join(parts)
