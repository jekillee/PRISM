"""
Job specifications for the PRISM batch SDK.

A JobSpec is a frozen, hashable description of *what* to compute (never *how*).
Its canonical form drives the archive cache key, so two specs with identical
parameters resolve to the same cached record. Floats are rounded before hashing
so that text like "2.0" vs "2.0000001" from a manifest never thrashes the cache.

Note: msign is part of the spec because the stored arrays are sign-filtered to
match the GUI's saved NPZ (pos/neg slice the amplitude, abs takes |mode|). If a
future revision archives the full signed array instead, drop msign from the key.
"""

import hashlib
import json
import math
from dataclasses import dataclass, asdict
from typing import Optional


def _canon_value(value, ndigits=6):
    """Normalize a value for hashing so equivalent spellings collide.

    Numbers are rounded, and integer-valued numbers are collapsed to int, so
    tmin=2, tmin=2.0 and tmin=2.0000001 all canonicalize identically (avoiding
    cache thrash from int/float spelling or float-text noise). bool is left
    as-is (it is technically an int but means something different).
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        rounded = round(float(value), ndigits)
        if not math.isfinite(rounded):
            return rounded   # NaN/inf: keep as-is; int() would raise, json/hash handle it
        return int(rounded) if rounded == int(rounded) else rounded
    if isinstance(value, (list, tuple)):
        return [_canon_value(v, ndigits) for v in value]
    if isinstance(value, dict):
        return {k: _canon_value(v, ndigits) for k, v in value.items()}
    return value


def _canonical_json(params):
    return json.dumps(_canon_value(params), sort_keys=True, default=str)


@dataclass(frozen=True)
class JobSpec:
    """Base job spec. Subclasses set SUBSYSTEM and declare their parameter fields."""

    shot: int

    # Subclasses override this class attribute (NOT a dataclass field).
    SUBSYSTEM = "base"

    def params(self) -> dict:
        """All spec fields as a plain dict."""
        return asdict(self)

    def hash_params(self) -> dict:
        """Params that define the analysis, EXCLUDING shot. The shot is already in
        the filename/dir, so dropping it here means a fixed parameter set yields the
        same spec_hash for every shot (a --shot-range batch shares one suffix)."""
        return {k: v for k, v in self.params().items() if k != "shot"}

    def canonical(self) -> str:
        """Deterministic JSON of the analysis params (sorted keys, rounded floats)."""
        return _canonical_json(self.hash_params())

    @property
    def spec_hash(self) -> str:
        """Short stable hash of the analysis params (shot-independent)."""
        return hashlib.sha1(self.canonical().encode("utf-8")).hexdigest()[:12]


@dataclass(frozen=True)
class NModeJobSpec(JobSpec):
    """N-mode (toroidal mode number) spectrum job.

    Mirror of NModeSpectrumTab's inputs. The Mirnov coil set is not specified
    here: it is auto-derived from the shot via config/mirnov_config.json.

    tmin/tmax default to None = "full shot": compute_nmode then loads the whole
    shot and uses its actual data range (tmin clamped up to 0.0, like the GUI's
    "Use full shot length"). A None bound is part of the cache key, so a full-shot
    request caches separately from any explicit window.
    """

    tmin: Optional[float] = None
    tmax: Optional[float] = None
    t_interval: float = 0.01
    fmin: float = 0.0
    fmax: float = 100.0
    tol: float = 0.8
    nmodes: int = 5
    frac: float = 1e-2
    msign: int = 1            # 0=abs, 1=pos, -1=neg, 2=all
    integrate: bool = False
    detrend: bool = True

    SUBSYSTEM = "nmode"
