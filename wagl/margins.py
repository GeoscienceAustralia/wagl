"""
Image Margins
"""
from __future__ import absolute_import, print_function

class ImageMargins(object):

    """
    Holds some value for each side of an image. This was initially
    created to hold buffer widths (in pixels) for a scene, but does
    not care about the type of the values passed to the constructer
    and can hence be used for any type.
    """

    def __init__(self, left, right=None, top=None, bottom=None):
        """
        Constructor.

        The arguments are copied directly to members of the same name
        with the restrictions that:

        - if right is None then top and bottom must also be None, and

        - if right is not None then top and bottom must not be none
          either.
        """
        if right is None:
            msg = "if right is None then top and bottom must also be None"
            assert top is None and bottom is None, msg
            self.left = self.right = self.top = self.bottom = left
        else:
            msg = ("if right is not None then top and bottom must "
                   "also not be None")
            assert top is not None and bottom is not None, msg
            self.left = left
            self.right = right
            self.top = top
            self.bottom = bottom

    def __str__(self):
        msg = "ImageMargins({left}, {right}, {top}, {bottom})"
        msg = msg.format(left=self.left, right=self.right, top=self.top,
                         bottom=self.bottom)
        return msg
