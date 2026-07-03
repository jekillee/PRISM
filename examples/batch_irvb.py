#!/usr/bin/env python3
"""Example: compute an IRVB result with the PRISM batch API and inspect it.

Run on nkstar/ukstar (inside the KSTAR network) with the `prism python` helper, which
supplies PRISM's bundled Python + import path (no env setup needed):

    prism python examples/batch_irvb.py        # from the PRISM install directory
    prism python /path/to/your/copy.py         # or a copy anywhere

Nothing else to configure: results cache under ~/prism_results/irvb/ (override with
the $PRISM_ARCHIVE_ROOT env var). Re-running is instant (cache hit). To compute
in-memory with no file, use `from batch.compute import compute_irvb` instead of `run`.
"""
import numpy as np

from batch import IRVBJobSpec, run

SHOT = 40085

# 1) Describe WHAT to compute (not how): full shot, EFIT tree, psi_N region bounds.
spec = IRVBJobSpec(shot=SHOT, efit_tree="efit01", psi_boundaries=(0.7, 1.0))

# 2) Run it. status is "computed" the first time, "cached" afterwards.
result, status = run(spec)
print(f"status = {status}\n")

# 3) What comes out: printing the result summarizes every field and its shape.
print(result)

# 4) Access fields directly.
print(f"\n{len(result.time)} frames; regions = {result.meta['region_labels']}")
print(f"peak core-region Prad = {result.region_prad[0].max():.3f} MW")

# 5) Recompute Prad for ANY psi_N range straight from the arrays (psi_n shares the
#    prad_2d grid), not just the boundaries above:
R, Z = result.R, result.Z
RR, _ = np.meshgrid(R, Z)
vol = 2 * np.pi * RR * (R[1] - R[0]) * (Z[1] - Z[0])          # 2*pi*R*dR*dZ [m^3]
mask = (result.psi_n >= 0.3) & (result.psi_n < 0.85)
prad_custom = np.sum(result.prad_2d * mask * vol, axis=(1, 2))   # [MW] per frame
print(f"custom region 0.3<=psi_N<0.85: "
      f"{prad_custom.min():.3f} .. {prad_custom.max():.3f} MW")
