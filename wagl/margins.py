"""
Defining and buffering image extents/margins.
"""
from __future__ import absolute_import, print_function
from math import ceil


class ImageMargins:

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
            msg = "if right is not None then top and bottom must " "also not be None"
            assert top is not None and bottom is not None, msg
            self.left = left
            self.right = right
            self.top = top
            self.bottom = bottom

    def __str__(self):
        msg = "ImageMargins({left}, {right}, {top}, {bottom})"
        msg = msg.format(
            left=self.left, right=self.right, top=self.top, bottom=self.bottom
        )
        return msg


def pixel_buffer(acquisition, distance=8000):
    """
    Determine a buffer in pixel units given an `Acquisition` and
    distance.
    If the acquistion's pixel units are in metres, then a distance
    of 8000 would equate to 8000 metres.
    The result of the number of pixels to buffer is rounded up to
    the next whole integer.
    For determining the approproate distance to use as a buffer
    within your region of interest. You need to take into account
    not just the highest elevation, but also steepest solar angle.
    For Australia, this was roughly 6.25km, and in order to be
    extra conservative, a default value of 8km was selected.

    :param acquisition:
        An instance of an `Acquistion` object.

    :param distance:
        A number representing the desired distance (in the same
        units as the acquisition) in which to calculate the extra
        number of pixels required to buffer an image.
        Default is 8000.

    :return:
        An instance of an `ImageMargins` object with each of:

        * ImageMargins.left
        * ImageMargins.right
        * ImageMargins.top
        * ImageMargins.bottom

        set to buffer in pixel units equivalent to that given by
        distance.
    """
    pixels = [ceil(distance / i) for i in acquisition.resolution]
    margins = ImageMargins(pixels[0], pixels[0], pixels[1], pixels[1])
    return margins
