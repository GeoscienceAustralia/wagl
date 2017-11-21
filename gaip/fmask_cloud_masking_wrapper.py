from __future__ import absolute_import, print_function
from gaip.acca_cloud_masking import majority_filter

from . import fmask_cloud_masking as _fmask

def fmask_cloud_mask(mtl, null_mask=None, cloud_prob=None, wclr_max=None,
                     sat_tag=None, aux_data=None):

    Lnum = int(sat_tag[-1:])
    (_, _, _, _, _,
     _, _, _, fmask_byte,
     _, _, _, _, _,
     _, _) = _fmask.plcloud(filename=mtl, mask=null_mask, num_Lst=Lnum,
                            aux_data=aux_data or {})

    # Convert to bool, True = Cloud, False not Cloud
    fmask_byte = fmask_byte == 1
    # Use a majority filter to fill holes, 2 iterations works well to smoothe
    # things over
    fmask_byte = majority_filter(array=fmask_byte, iterations=2)

    return ~fmask_byte # Invert
