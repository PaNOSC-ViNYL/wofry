"""
Represents an ideal lens.
"""
from syned.beamline.optical_elements.ideal_elements.screen import Screen
from wofry.beamline.decorators import WOOpticalElementDecorator

class WOScreen(Screen, WOOpticalElementDecorator):
    def __init__(self, name="Undefined"):
        Screen.__init__(self, name=name)

    def applyOpticalElement(self, wavefront, parameters=None):
        return wavefront