#! /usr/bin/env python

import os, sys, re, datetime

class PQALogExtractor(object):

    def __init__(self, pqa_tif_path):

        logfile_list = [os.path.join(os.path.dirname(pqa_tif_path), logfile) for logfile in['ACCA_LOGFILE.txt',
                        'FMASK_LOGFILE.txt',
                        'ACCA_CLOUD_SHADOW_LOGFILE.txt',
                        'FMASK_CLOUD_SHADOW_LOGFILE.txt']]

        # Returns the required line from a list of strings
        def linefinder(array, string = ""):
            '''Searches a list for the specified string.

               Args:
                   array: A list containing searchable strings.
                   string: User input containing the string to search.

               Returns:
                   The line containing the found sting.
            '''

            for line in array:
                if string in str(line):
                    return line

        #=======================================================================
        # for pqa_tif_path in logfile_list:
        #     f = open(pqa_tif_path, 'r')
        #     log = f.readlines()
        #     f.close()
        #=======================================================================

        # get ACCA percent
        f = open(logfile_list[0], 'r')
        acca_log = f.readlines()
        f.close()
        find  = linefinder(acca_log,'Final Cloud Layer Percent')
        self.acca_percent = float(find.split()[-1])

        #get Fmask percent
        f = open(logfile_list[1], 'r')
        fmask_log = f.readlines()
        f.close()
        find  = linefinder(fmask_log,'Final Cloud Layer Percent')
        self.fmask_percent = float(find.split()[-1])

        #get ACCA cloud shadow
        f = open(logfile_list[2], 'r')
        acca_cloud_shadow_log = f.readlines()
        f.close()
        find  = linefinder(acca_cloud_shadow_log,'Cloud Shadow Percent')
        self.acca_cloud_shadow_percent = float(find.split()[-1])

        #get Fmask cloud shadow
        f = open(logfile_list[3], 'r')
        fmask_cloud_shadow_log = f.readlines()
        f.close()
        find  = linefinder(fmask_cloud_shadow_log,'Cloud Shadow Percent')
        self.fmask_cloud_shadow_percent = float(find.split()[-1])


        #extract which tests have been run from TIF file name and return the list
        self.tests_run = list((os.path.splitext(os.path.basename(pqa_tif_path))[0]).split('_')[-1])

    #===========================================================================
    #     run_not_run = {'0' : 'Not Run',
    #                    '1' : 'Run'
    #                   }
    #     sat = os.path.basename(pqa_tif_path).split('_')[0]
    #
    #     LS5_band6_run_not_run = {'0' : 'Not Run',
    #                              '1' : 'Duplicated Band61'
    #                             }
    #
    #     sensor_test = {'LS5': LS5_band6_run_not_run[tests_run[6]],
    #                  'LS7': run_not_run[tests_run[6]]
    #                 }
    #     self.sat_band1  = run_not_run[tests_run[0]]
    #     self.sat_band2  = run_not_run[tests_run[1]]
    #     self.sat_band3  = run_not_run[tests_run[2]]
    #     self.sat_band4  = run_not_run[tests_run[3]]
    #     self.sat_band5  = run_not_run[tests_run[4]]
    #     self.sat_band61 = run_not_run[tests_run[5]]
    #     self.sat_band62 = sensor_test[sat]
    #     self.sat_band7  = run_not_run[tests_run[7]]
    #     self.contiguity = run_not_run[tests_run[8]]
    #     self.land_sea   = run_not_run[tests_run[9]]
    #     acca       = run_not_run[tests_run[10]]
    #     fmask      = run_not_run[tests_run[11]]
    #     acca_shad  = run_not_run[tests_run[12]]
    #     fmask_shad = run_not_run[tests_run[13]]
    #     empty_1    = run_not_run[tests_run[14]]
    #     empty_2    = run_not_run[tests_run[15]]
    #
    #===========================================================================

