#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function

from gaip.modtran_profiles import (
    MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO,
    MIDLAT_SUMMER_TRANSMITTANCE, TROPICAL_TRANSMITTANCE
)


def main():
    with open('data/test-midlat-summer.tp5', 'w') as src:
        src.write(MIDLAT_SUMMER_ALBEDO.format(albedo=0.0,
                                              water=1.07000122070313,
                                              ozone=0.28499999642372,
                                              filter_function='landsat7_vsir.flt',
                                              visibility=-0.02264800000000,
                                              elevation=0.70900000000000,
                                              sat_height=705.0,
                                              sat_view=171.000748,
                                              doy=212,
                                              lat=-29.33856209871443,
                                              lon=209.88857485506449,
                                              time=23.73920805027778,
                                              sat_azimuth=279.408417))
    
    
    # TODO: get a tropical dataset
    with open('data/test-tropical.tp5', 'w') as src:
        src.write(TROPICAL_ALBEDO.format(albedo=0.0,
                                         water=1.07000122070313,
                                         ozone=0.28499999642372,
                                         filter_function='landsat7_vsir.flt',
                                         visibility=-0.02264800000000,
                                         elevation=0.70900000000000,
                                         sat_height=705.0,
                                         sat_view=171.000748,
                                         doy=212,
                                         lat=-29.33856209871443,
                                         lon=209.88857485506449,
                                         time=23.73920805027778,
                                         sat_azimuth=279.408417))


    with open('data/test-midlat-summer-trans.tp5', 'w') as src:
        src.write(MIDLAT_SUMMER_TRANSMITTANCE.format(albedo=0.0,
                                                     water=1.07000122070313,
                                                     ozone=0.28499999642372,
                                                     filter_function='landsat7_vsir.flt',
                                                     visibility=-0.02264800000000,
                                                     elevation=0.70900000000000,
                                                     sat_height=705.0,
                                                     sat_view=171.000748,
                                                     doy=212,
                                                     sat_view_offset=180.0 - 171.000748))
    
    
    # TODO: get a tropical dataset
    with open('data/test-tropical-trans.tp5', 'w') as src:
        src.write(TROPICAL_TRANSMITTANCE.format(albedo=0.0,
                                                water=1.07000122070313,
                                                ozone=0.28499999642372,
                                                filter_function='landsat7_vsir.flt',
                                                visibility=-0.02264800000000,
                                                elevation=0.70900000000000,
                                                sat_height=705.0,
                                                sat_view=171.000748,
                                                doy=212,
                                                sat_view_offset=180.0 - 171.000748))


    with open('data/TL_alb_0.tp5', 'r') as src:
        ref_albedo = src.readlines()

    with open('data/test-midlat-summer.tp5', 'r') as src:
        test_albedo = src.readlines()

    print("Testing mid lat albedo")
    for i in range(len(ref_albedo)):
        if not ref_albedo[i] == test_albedo[i]:
            print("Line {} not equivilent".format(i))

    with open('data/TL_alb_t.tp5', 'r') as src:
        ref_trans = src.readlines()

    with open('data/test-midlat-summer-trans.tp5', 'r') as src:
        test_trans = src.readlines()

    print("Testing mid lat transmittance")
    for i in range(len(ref_trans)):
        if not ref_trans[i] == test_trans[i]:
            print("Line {} not equivilent".format(i))


if __name__ == '__main__':
    main()
