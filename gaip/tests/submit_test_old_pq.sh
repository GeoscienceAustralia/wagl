#!/bin/bash
#
# Submit test runs for PQA using old image_processor code
#
# Set of 10 L1T datasets for PQ process testing
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS5_TM_OTH_P51_GALPGS01-002_100_081_20100228,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS5_TM_NBAR_P54_GANBAR01-002_100_081_20100228 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS5_TM_OTH_P51_GALPGS01-002_105_080_20100215,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS5_TM_NBAR_P54_GANBAR01-002_105_080_20100215 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS5_TM_OTH_P51_GALPGS01-002_107_076_20100213,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS5_TM_NBAR_P54_GANBAR01-002_107_076_20100213 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS5_TM_OTH_P51_GALPGS01-002_110_073_20100218,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS5_TM_NBAR_P54_GANBAR01-002_110_073_20100218 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS7_ETM_OTH_P51_GALPGS01-002_107_074_20140710,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS7_ETM_NBAR_P54_GANBAR01-002_107_074_20140710 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS7_ETM_OTH_P51_GALPGS01-002_110_081_20140731,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS7_ETM_NBAR_P54_GANBAR01-002_110_081_20140731 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS7_ETM_OTH_P51_GALPGS01-002_113_079_20140720,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS7_ETM_NBAR_P54_GANBAR01-002_113_079_20140720 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS7_ETM_OTH_P51_GALPGS01-002_115_076_20140718,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS7_ETM_NBAR_P54_GANBAR01-002_115_076_20140718 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS8_OLITIRS_OTH_P51_GALPGS01-002_115_075_20141014,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS8_OLI_TIRS_NBAR_P54_GANBAR01-002_115_075_20141014 test_old_pq.pbs
qsub -v L1T_PATH=/g/data1/v10/projects/Luigi_work_flow_test/L1T/LS8_OLITIRS_OTH_P51_GALPGS01-015_101_078_20141012,NBAR_PATH=/g/data1/v10/projects/Luigi_work_flow_test/NBAR/LS8_OLI_TIRS_NBAR_P54_GANBAR01-015_101_078_20141012 test_old_pq.pbs

