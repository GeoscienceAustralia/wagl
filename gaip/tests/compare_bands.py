from __future__ import absolute_import
from __future__ import print_function
import rasterio as rio
import sys
import gaip

f1 = sys.argv[1]
f2 = sys.argv[2]
f3 = sys.argv[3]

print("files are %s and %s" % (f1, f2))

with rio.open(f1) as ds1:
    with rio.open(f2) as ds2:
        assert(ds1.count == ds2.count)

        for band in range(1, ds1.count+1):
            d1 = ds1.read(band)
            d2 = ds2.read(band)

            diff = sum(sum(d1-d2))
            print("band=%d, diff=%d" % (band, diff))

            gaip.write_img((d1-d2), f3, fmt='GTiff')
            
