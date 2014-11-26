import rasterio as rio
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]

print "files are %s and %s" % (f1, f2)

with rio.open(f1) as ds1:
    with rio.open(f2) as ds2:
        assert(ds1.count == ds2.count)

        for band in range(1, ds1.count+1):
            d1 = ds1.read_band(band)
            d2 = ds2.read_band(band)

            diff = sum(sum(d1-d2))
            print "band=%d, diff=%d" % (band, diff)
