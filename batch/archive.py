"""
npz archive for the PRISM batch SDK.

A result is stored as one .npz file per shot, named simply by shot:

    <archive_root>/<subsystem>/<subsystem>_<shot>.npz   (e.g. nmode/nmode_40848.npz)

One file per shot (no parameter hash in the name): nmode_<shot>.npz. A run reuses an
existing file unless overwrite=True; the CLI checks for an existing file per shot and
asks whether to overwrite or skip. The file's metadata records the parameters it was
computed with. No version segment, so a result persists across PRISM version updates.

Writes are atomic (write a .tmp.npz sibling, then os.replace) so a concurrent GUI
or batch session never sees a torn file.
"""

import os
import uuid
import zipfile
from pathlib import Path


def cache_key(spec) -> str:
    return f"{spec.SUBSYSTEM}_{spec.shot}"


def archive_path(archive_root, spec) -> Path:
    # Flat under the subsystem dir, one file per shot (e.g. nmode/nmode_40848.npz).
    return Path(archive_root) / spec.SUBSYSTEM / (cache_key(spec) + ".npz")


def save_record(archive_root, spec, result) -> Path:
    """Atomically write `result` to its archive path. Returns the final path."""
    path = archive_path(archive_root, spec)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Unique per-writer tmp name (pid + uuid) so a concurrent GUI + batch write of
    # the same key never shares or clobbers a sibling tmp; still matches *.tmp.*.
    tmp = path.with_name(f"{path.stem}.{os.getpid()}.{uuid.uuid4().hex}.tmp.npz")
    try:
        with open(tmp, "wb") as fh:
            result.save_npz(fh)          # write into the tmp handle, no extension games
        os.replace(str(tmp), str(path))  # atomic rename on the same filesystem
    finally:
        if tmp.exists():
            tmp.unlink()                 # honor the *.tmp.* cleanup rule
    return path


def load_record(archive_root, spec, result_cls):
    """Return the cached result, or None if absent / unreadable (-> recompute)."""
    path = archive_path(archive_root, spec)
    if not path.exists():
        return None
    try:
        return result_cls.load_npz(str(path))
    except (OSError, ValueError, EOFError, KeyError, zipfile.BadZipFile):
        # torn / unreadable / structurally-valid-but-incomplete (a missing array
        # raises KeyError) -> treat as a miss so get_or_compute self-heals.
        return None


def get_or_compute(archive_root, spec, compute_fn, result_cls, overwrite=False):
    """Return (result, status) where status is 'cached' or 'computed'.

    With overwrite=False, an existing readable file for this shot is reused as-is
    ('cached'). With overwrite=True, the shot is recomputed and the file overwritten.
    The CLI decides per shot (prompting on collisions) and passes overwrite=True for
    the shots the user chose to (re)compute.
    """
    if not overwrite:
        hit = load_record(archive_root, spec, result_cls)
        if hit is not None:
            return hit, "cached"
    result = compute_fn(spec)
    save_record(archive_root, spec, result)
    return result, "computed"
