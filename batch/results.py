"""
Typed result containers for the PRISM batch SDK, saved as .npz (one file per job).

The .npz layout matches what the GUI's n-mode tab already writes, so existing
NPZ-reading scripts (and the tab's bundled example reader) keep working unchanged:

    metadata       : 0-d object array holding the provenance dict
                     (load with allow_pickle=True, then .item())
    time           : window-center times [s]
    frequency      : frequency axis [kHz]
    mode_spectrum  : signed toroidal mode number, shape (n_time, n_freq)
    amplitude      : per-mode amplitude evolution, shape (n_mode_rows, n_time)
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class NModeResult:
    """N-mode spectrum result.

    Attributes:
        shot:          KSTAR shot number
        time:          window-center times [s], shape (n_time,)
        frequency:     frequency axis [kHz], shape (n_freq,)
        mode_spectrum: signed toroidal mode number, shape (n_time, n_freq)
        amplitude:     per-mode amplitude evolution, shape (n_mode_rows, n_time)
        meta:          provenance dict (shot, params, version, created_iso, ...)
    """

    shot: int
    time: np.ndarray
    frequency: np.ndarray
    mode_spectrum: np.ndarray
    amplitude: np.ndarray
    meta: dict = field(default_factory=dict)

    SUBSYSTEM = "nmode"

    def save_npz(self, file):
        """Write this result as .npz to `file` (a path or open binary file object).

        Layout is identical to the GUI's saved n-mode NPZ (metadata stored as a
        0-d object array), so the same reader scripts work for both.
        """
        np.savez(
            file,
            metadata=np.array(self.meta),
            time=self.time,
            frequency=self.frequency,
            mode_spectrum=self.mode_spectrum,
            amplitude=self.amplitude,
        )

    @classmethod
    def load_npz(cls, path):
        """Rebuild an NModeResult from a .npz file written by save_npz."""
        data = np.load(path, allow_pickle=True)
        try:
            meta = data["metadata"].item()
            if not isinstance(meta, dict):
                meta = {}
            return cls(
                shot=int(meta.get("shot", 0)),
                time=data["time"],
                frequency=data["frequency"],
                mode_spectrum=data["mode_spectrum"],
                amplitude=data["amplitude"],
                meta=dict(meta),
            )
        finally:
            data.close()
