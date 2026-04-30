"""
TV (Visible Camera) data loader.

Resolves TV ZIP file paths under /Diag_TV/{campaign}/ for a given shot.

Lookup strategy:
    1. JSON config (config/tv_campaigns.json) — fast O(1) lookup by shot range
    2. Fallback: scan /Diag_TV for campaign-pattern directories and check if
       the shot's ZIP exists. Logs a warning asking the user to notify
       Jekil Lee so the config can be updated.
"""

import json
import os
import re
from pathlib import Path
from typing import List, Optional


# TV camera parameters
TV_FPS = 210.0
TV_OFFSET = 0.1  # 100ms trigger offset


_CONFIG_PATH = Path(__file__).parent.parent / 'config' / 'tv_campaigns.json'
with open(_CONFIG_PATH, 'r') as f:
    _CONFIG = json.load(f)

_BASE_PATH = _CONFIG.get('base_path', '/Diag_TV')
_CAMPAIGNS = _CONFIG['campaigns']

_CAMPAIGN_DIR_PATTERN = re.compile(r'^\d{4}C\d{2}$')

# Session cache for fallback lookups: {shot: campaign}
_fallback_cache = {}


def _campaign_from_config(shot_number: int) -> Optional[str]:
    """Look up campaign in JSON config by shot range."""
    for entry in _CAMPAIGNS:
        smax = entry['shot_max']
        if smax is None:
            smax = float('inf')
        if entry['shot_min'] <= shot_number <= smax:
            return entry['name']
    return None


def _campaign_from_scan(shot_number: int, tv_name: str,
                       startup: bool = False) -> Optional[str]:
    """Fallback: scan /Diag_TV directory for matching campaign folder.

    Returns campaign name if a ZIP file is found for this shot, else None.
    """
    if shot_number in _fallback_cache:
        return _fallback_cache[shot_number]

    if not os.path.exists(_BASE_PATH):
        return None

    try:
        entries = sorted(os.listdir(_BASE_PATH), reverse=True)
    except Exception as e:
        print(f"[TV] Scan error: {e}")
        return None

    shot_str = f'{shot_number:06d}'
    tv_num = tv_name.replace('TV', '').lower()

    for entry in entries:
        if not _CAMPAIGN_DIR_PATTERN.match(entry):
            continue
        if startup:
            zip_path = os.path.join(_BASE_PATH, entry, f'{tv_name}_startup',
                                    f'{shot_str}_{tv_name.lower()}_startup.zip')
        else:
            zip_path = os.path.join(_BASE_PATH, entry, tv_name,
                                    f'{shot_str}_tv{tv_num}.zip')
        if os.path.exists(zip_path):
            _fallback_cache[shot_number] = entry
            print(f"[TV] Shot #{shot_number} not in tv_campaigns.json, "
                  f"found in {entry} via directory scan.")
            print(f"[TV] Please notify Jekil Lee (jklee@kfe.re.kr) "
                  f"to update config/tv_campaigns.json.")
            return entry

    return None


def get_campaign_for_shot(shot_number: int, tv_name: str = 'TV01',
                          startup: bool = False) -> Optional[str]:
    """Resolve campaign directory name for a given shot.

    Args:
        shot_number: KSTAR shot number.
        tv_name: 'TV01' or 'TV02'. Used only for fallback scan.
        startup: Whether to look for startup ZIP files.

    Returns:
        Campaign directory name (e.g. '2026C19') or None if not found.
    """
    campaign = _campaign_from_config(shot_number)
    if campaign:
        return campaign
    return _campaign_from_scan(shot_number, tv_name, startup=startup)


def get_tv_zip_path(shot_number: int, tv_name: str) -> Optional[str]:
    """Full path to TV ZIP file (regular plasma view)."""
    campaign = get_campaign_for_shot(shot_number, tv_name, startup=False)
    if campaign is None:
        return None
    shot_str = f'{shot_number:06d}'
    tv_num = tv_name.replace('TV', '').lower()
    return f'{_BASE_PATH}/{campaign}/{tv_name}/{shot_str}_tv{tv_num}.zip'


def get_tv_startup_zip_path(shot_number: int, tv_name: str) -> Optional[str]:
    """Full path to TV startup ZIP file."""
    campaign = get_campaign_for_shot(shot_number, tv_name, startup=True)
    if campaign is None:
        return None
    return (f'{_BASE_PATH}/{campaign}/{tv_name}_startup/'
            f'{shot_number:06d}_{tv_name.lower()}_startup.zip')


def find_available_tvs(shot_number: int, startup: bool = False) -> List[str]:
    """Return list of TV cameras that have ZIP data for this shot."""
    available = []
    for tv_name in ['TV01', 'TV02']:
        if startup:
            zip_path = get_tv_startup_zip_path(shot_number, tv_name)
        else:
            zip_path = get_tv_zip_path(shot_number, tv_name)
        if zip_path and os.path.exists(zip_path):
            available.append(tv_name)
            print(f"[TV] Found {zip_path}")
    return available
