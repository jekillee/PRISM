"""
Animated toggle switch widget for PySide6
"""

from PySide6.QtWidgets import QWidget
from PySide6.QtCore import Qt, QPropertyAnimation, QEasingCurve, Property, Signal, QSize
from PySide6.QtGui import QPainter, QColor


class ToggleSwitch(QWidget):
    """iOS-style animated toggle switch"""

    toggled = Signal(bool)

    def __init__(self, parent=None, checked=False):
        super().__init__(parent)
        self._checked = checked
        self._handle_pos = 1.0 if not checked else self._max_handle_pos()
        self.setFixedSize(36, 20)
        self.setCursor(Qt.PointingHandCursor)

        self._animation = QPropertyAnimation(self, b"handle_pos", self)
        self._animation.setDuration(150)
        self._animation.setEasingCurve(QEasingCurve.InOutCubic)

    def _max_handle_pos(self):
        return self.width() - self.height() + 1

    def sizeHint(self):
        return QSize(36, 20)

    def _get_handle_pos(self):
        return self._handle_pos

    def _set_handle_pos(self, pos):
        self._handle_pos = pos
        self.update()

    handle_pos = Property(float, _get_handle_pos, _set_handle_pos)

    def isChecked(self):
        return self._checked

    def setChecked(self, checked, animate=True):
        if self._checked == checked:
            return
        self._checked = checked
        target = self._max_handle_pos() if checked else 1.0
        if animate:
            self._animation.setStartValue(self._handle_pos)
            self._animation.setEndValue(target)
            self._animation.start()
        else:
            self._handle_pos = target
            self.update()
        self.toggled.emit(checked)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.setChecked(not self._checked)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)

        w, h = self.width(), self.height()
        r = h / 2

        # Track
        if self._checked:
            track_color = QColor("#4CAF50")
        else:
            track_color = QColor("#888888")
        p.setBrush(track_color)
        p.setPen(Qt.NoPen)
        p.drawRoundedRect(0, 0, w, h, r, r)

        # Handle
        handle_r = h - 4
        p.setBrush(QColor("#ffffff"))
        p.drawEllipse(int(self._handle_pos), 2, handle_r, handle_r)

        p.end()
