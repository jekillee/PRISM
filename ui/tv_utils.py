#!/usr/bin/python3.8

"""
TV (Visible Camera) shared utilities
Common constants and functions used by TV and TV Startup tabs
"""

import os


# TV camera parameters
TV_FPS = 210.0
TV_OFFSET = 0.1  # 100ms trigger offset


def get_year_from_shot(shot_number):
    """Get year from shot number (KSTAR shot ranges)"""
    if shot_number >= 40464:
        return 2026
    elif shot_number > 37741:
        return 2025
    elif shot_number > 34836:
        return 2024
    elif shot_number > 32768:
        return 2023
    elif shot_number > 30445:
        return 2022
    elif shot_number > 27400:
        return 2021
    elif shot_number > 24081:
        return 2020
    elif shot_number > 21758:
        return 2019
    elif shot_number > 19396:
        return 2018
    elif shot_number > 17376:
        return 2017
    elif shot_number > 14407:
        return 2016
    elif shot_number > 11724:
        return 2015
    elif shot_number > 9427:
        return 2014
    elif shot_number > 8354:
        return 2013
    elif shot_number > 6470:
        return 2012
    elif shot_number > 4468:
        return 2011
    elif shot_number > 2342:
        return 2010
    elif shot_number > 1283:
        return 2009
    else:
        return 2008


def get_campaign_from_shot(shot_number):
    """Get campaign folder name from shot number"""
    # Campaign boundaries (shot > threshold means this campaign)
    # Format: (threshold, campaign_name)
    campaigns = [
        (40464, "2026C19"),
        (38000, "2025C18"),
        (35000, "2025C17"),
        (34000, "2024C16"),
        (33000, "2023C16"),
        (27805, "2021C14"),
        (24330, "2020C13"),
        (0, "2019C12"),
    ]

    for threshold, campaign in campaigns:
        if shot_number > threshold:
            return campaign

    return "2019C12"


def get_campaign_from_year(year):
    """Get campaign folder name from year by scanning directory"""
    base_path = '/Diag_TV'

    if not os.path.exists(base_path):
        print(f"[TV] Base path {base_path} not found")
        return None

    try:
        entries = os.listdir(base_path)
        campaign_dirs = [d for d in entries if d.startswith(f'{year}C')]

        if campaign_dirs:
            campaign_dirs.sort()
            return campaign_dirs[-1]
        return None
    except Exception as e:
        print(f"[TV] Error listing campaigns: {str(e)}")
        return None


def get_tv_zip_path(shot_number, tv_name):
    """
    Get full path to TV ZIP file

    Args:
        shot_number: Shot number
        tv_name: 'TV01' or 'TV02'

    Returns:
        Full path to ZIP file or None if campaign not found
    """
    year = get_year_from_shot(shot_number)
    campaign = get_campaign_from_year(year)

    if campaign is None:
        return None

    shot_str = f'{shot_number:06d}'
    tv_num = tv_name.replace('TV', '').lower()
    return f'/Diag_TV/{campaign}/{tv_name}/{shot_str}_tv{tv_num}.zip'


def get_tv_startup_zip_path(shot_number, tv_name):
    """
    Get full path to TV startup ZIP file

    Args:
        shot_number: Shot number
        tv_name: 'TV01' or 'TV02'

    Returns:
        Full path to startup ZIP file
    """
    campaign = get_campaign_from_shot(shot_number)
    tv_num = tv_name.lower()  # tv01 or tv02
    return f"/Diag_TV/{campaign}/{tv_name}_startup/{shot_number:06d}_{tv_num}_startup.zip"


def find_available_tvs(shot_number, startup=False):
    """
    Find available TV01/TV02 ZIP files for given shot

    Args:
        shot_number: Shot number
        startup: If True, look for startup ZIP files

    Returns:
        List of available TVs (e.g., ['TV01', 'TV02'])
    """
    available_tvs = []

    for tv_name in ['TV01', 'TV02']:
        if startup:
            zip_path = get_tv_startup_zip_path(shot_number, tv_name)
        else:
            zip_path = get_tv_zip_path(shot_number, tv_name)

        if zip_path and os.path.exists(zip_path):
            available_tvs.append(tv_name)
            print(f"[TV] Found {zip_path}")

    return available_tvs


def frame_to_time(frame_idx, fps=TV_FPS, offset=TV_OFFSET):
    """
    Convert frame index to time in seconds

    Args:
        frame_idx: Frame index (0-based)
        fps: Frames per second
        offset: Time offset in seconds

    Returns:
        Time in seconds
    """
    return (frame_idx + 1) / fps - offset


def time_to_frame(time_sec, total_frames, fps=TV_FPS, offset=TV_OFFSET):
    """
    Convert time in seconds to nearest frame index

    Args:
        time_sec: Time in seconds
        total_frames: Total number of frames
        fps: Frames per second
        offset: Time offset in seconds

    Returns:
        Frame index (0-based, clamped to valid range)
    """
    frame_idx = int(round((time_sec + offset) * fps - 1))
    return max(0, min(frame_idx, total_frames - 1))


def frame_to_time_ms(frame_idx, fps=TV_FPS):
    """
    Convert frame index (1-based) to time in milliseconds
    Used for startup comparison where frame 1 = 0ms

    Args:
        frame_idx: Frame index (1-based)
        fps: Frames per second

    Returns:
        Time in milliseconds
    """
    return (frame_idx - 1) / fps * 1000
