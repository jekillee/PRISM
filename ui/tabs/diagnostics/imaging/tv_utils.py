"""
TV (Visible Camera) UI helpers — frame/time conversions.

Path resolution and ZIP file lookup are in data_loaders/tv_loader.py.
TV_FPS / TV_OFFSET are re-exported from tv_loader for backward compatibility.
"""

from data_loaders.tv_loader import TV_FPS, TV_OFFSET


def frame_to_time(frame_idx, fps=TV_FPS, offset=TV_OFFSET):
    """Convert frame index (0-based) to time in seconds."""
    return (frame_idx + 1) / fps - offset


def time_to_frame(time_sec, total_frames, fps=TV_FPS, offset=TV_OFFSET):
    """Convert time in seconds to nearest frame index (clamped)."""
    frame_idx = int(round((time_sec + offset) * fps - 1))
    return max(0, min(frame_idx, total_frames - 1))


def frame_to_time_ms(frame_idx, fps=TV_FPS):
    """Convert frame index (1-based) to time in milliseconds.
    Used for startup comparison where frame 1 = 0ms.
    """
    return (frame_idx - 1) / fps * 1000
