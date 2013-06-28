#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define _GNU_SOURCE
#include <getopt.h>

#define FORTRAN_DATA

#define TRUE 1
#define FALSE 0

#define DEBUG 0
#define STRSZ 256

#define FLOAT_ATOF(x) ((float) atof(x))



typedef struct _cmdlinedata {
    float west, east, south, north;
    char pix_file[STRSZ], cmp_file[STRSZ], date[STRSZ], time[STRSZ];
    int status;
} CMDLINE_DATA;


typedef struct _pixdata {
    long int nx, ny;
    long int nrecs;
    float *xgrid, *ygrid;        // [nx], [ny]
    short int *idxlon, *idxlat;  // [nrecs]
    short int *dt[3], *tm[3];    // [nrecs]
    float *lat, *lon, *aot;      // [nrecs]
} PIXDATA;


typedef struct _cmpdata {
    int nx, ny;
    float *x;    // [nx]
    float *y;    // [ny]
    float *aot;  // [nx * ny]
} CMPDATA;


//
// Globals
//

double X24 = 1.0 / 24;
double X60 = 1.0 / 60;

double SCENE_DT_TOLERANCE = 0.5;

const char *exe_name = "aot_loader";

static int print_values;


//
// Function decs
//


PIXDATA* load_pix(char *fname);

CMPDATA* load_cmp(char *fname);

void print_pix_entry(PIXDATA *pix, int k);

double decimal_date_time(short int *year, short int *month, short int *day,
                         short int *hour, short int *min, short int *sec);

int point_in_region(float *lon, float *lat,
                    float *lon_min, float *lon_max,
                    float *lat_min, float *lat_max);

float eval_pix_aot(PIXDATA* pix, float lon_min, float lon_max, float lat_min, float lat_max,
                                 short int year, short int month, short int day,
                                 short int hour, short int min, short int sec);

CMDLINE_DATA* parse_args(int argc, char *argv[]);

void usage(void);

int test_dt_diff(void);

int test_pix(void);

#ifdef FORTRAN_DATA
int read_frecord(void *buffer, int n, size_t element_size, FILE *fp);
#endif




#ifdef FORTRAN_DATA
int read_frecord(void *buffer, int n, size_t element_size, FILE *fp)
{
    //
    // Reads a FORTRAN data record.
    // Not required for second AATSR data format.
    //

    int head_value, tail_value;
    int nread;
    int nchk;

    if (DEBUG) {
        printf("\nread_frecord:\n");
        printf("\tn             = %d\n", n);
        printf("\telement_size  = %d bytes\n", (int)element_size);
    }

    assert(buffer);
    assert(fp);

    // Record head

    nchk = fread(&head_value, sizeof(int), 1, fp);
    assert(nchk == 1);
    if (DEBUG)
        printf("\thead_value    = %d bytes\n", head_value);
    assert(head_value == element_size * n);

    // Data values

    nread = fread(buffer, element_size, n, fp);
    assert(nread == n);

    // Record tail

    nchk = fread(&tail_value, sizeof(int), 1, fp);
    assert(nchk == 1);
    if (DEBUG)
        printf("\ttail_value    = %d bytes\n", tail_value);
    assert(head_value == tail_value);

    return nread;
}
#endif



PIXDATA* load_pix(char *fname)
{
    int i;

    FILE *fp = fopen(fname, "rb");
    assert(fp);

    // allocate data structure

    PIXDATA *data = (PIXDATA*) malloc(sizeof(PIXDATA));
    assert(data);

    // dimensions

    fread(&(data->nx), 1, sizeof(int), fp);
    fread(&(data->ny), 1, sizeof(int), fp);
    fread(&(data->nrecs), 1, sizeof(int), fp);

    // lon and lat grid coordinates

    data->xgrid = (float*) malloc(data->nx * sizeof(float));
    data->ygrid = (float*) malloc(data->ny * sizeof(float));

    fread(data->xgrid, data->nx, sizeof(float), fp);
    fread(data->ygrid, data->ny, sizeof(float), fp);

    // lon and lat indices

    data->idxlon = (short int*) malloc(data->nrecs * sizeof(short int));
    data->idxlat = (short int*) malloc(data->nrecs * sizeof(short int));

    fread(data->idxlon, data->nrecs, sizeof(short int), fp);
    fread(data->idxlat, data->nrecs, sizeof(short int), fp);

    // dt arrays: year, month, day

    for (i=0; i<3; i++) {
        data->dt[i] = (short int*) malloc(data->nrecs * sizeof(short int));
        fread(data->dt[i], data->nrecs, sizeof(short int), fp);
    }

    // tm arrays: hour. minute, second

    for (i=0; i<3; i++) {
        data->tm[i] = (short int*) malloc(data->nrecs * sizeof(short int));
        fread(data->tm[i], data->nrecs, sizeof(short int), fp);
    }

    // lat, lon, aot

    data->lat = (float*) malloc(data->nrecs*sizeof(float));
    data->lon = (float*) malloc(data->nrecs*sizeof(float));
    data->aot = (float*) malloc(data->nrecs*sizeof(float));

    fread(data->lat, data->nrecs, sizeof(float), fp);
    fread(data->lon, data->nrecs, sizeof(float), fp);
    fread(data->aot, data->nrecs, sizeof(float), fp);

    printf("LOADED PIX: %s\n", fname);

    return data;
}


CMPDATA* load_cmp(char *fname)
{
    FILE *fp = fopen(fname, "rb");
    assert(fp);

    // allocate data structure

    CMPDATA *data = (CMPDATA*) malloc(sizeof(CMPDATA));

    data->x = (float*) malloc(data->nx*sizeof(float));
    data->y = (float*) malloc(data->ny*sizeof(float));
    data->aot = (float*) malloc(data->nx*data->ny*sizeof(float));

    // dimensions

    fread(&(data->nx), 1, sizeof(int), fp);
    fread(&(data->ny), 1, sizeof(int), fp);

    // x, y, aot

    data->x = (float*) malloc(data->nx*sizeof(float));
    data->y = (float*) malloc(data->ny*sizeof(float));
    data->aot = (float*) malloc(data->nx*data->ny*sizeof(float));

    fread(data->x, data->nx, sizeof(float), fp);
    fread(data->y, data->ny, sizeof(float), fp);
    fread(data->aot, data->nx*data->ny, sizeof(float), fp);

    fclose(fp);

    printf("LOADED CMP: %s\n", fname);

    return data;
}


void print_pix_entry(PIXDATA *pix, int k)
{
    char *pix_fmt = "[%8d] %3d %3d  (%4d %2d %2d %2d %2d %2d)  %f %f %f\n";

    printf(pix_fmt, k, pix->idxlon[k], pix->idxlat[k],
                       pix->dt[0][k], pix->dt[1][k], pix->dt[2][k],
                       pix->tm[0][k], pix->tm[1][k], pix->tm[2][k],
                       pix->lon[k], pix->lat[k], pix->aot[k]);
}


double decimal_date_time(short int *year, short int *month, short int *day,
                         short int *hour, short int *min, short int *sec)
{
    return ((double) (*year) * 10000 + (double) (*month) * 100 + (double) (*day)) +
           ((double) (*hour) + ((double) (*min) + (double) (*sec) * X60) * X60) * X24;
}


int point_in_region(float *lon, float *lat,
                    float *lon_min, float *lon_max,
                    float *lat_min, float *lat_max)
{
    return (*lon >= *lon_min && *lon <= *lon_max && *lat >= *lat_min && *lat <= *lat_max);
}


float eval_pix_aot(PIXDATA* pix, float lon_min, float lon_max, float lat_min, float lat_max,
                                 short int year, short int month, short int day,
                                 short int hour, short int min, short int sec)
{
    int k;
    float aot_sum = 0;
    int aot_count = 0;
    int in_scene_count = 0;

    assert(pix);

    double scene_dt = decimal_date_time(&year, &month, &day, &hour, &min, &sec);

    for (k=0; k<pix->nrecs; k++) {

        // From Yi:
        // "	lat/lon	- The mean latitude and longitude of the native pixels
        //                within the cells. They might be slightly different
        //                from the grid latitude/longitude. No particular use,
        //                can be ignored."
        //
        //float lon = pix->lon[k];
        //float lat = pix->lat[k];

        float lon = pix->xgrid[ (int) pix->idxlon[k] ];
        float lat = pix->ygrid[ (int) pix->idxlat[k] ];

        if (point_in_region(&lon, &lat, &lon_min, &lon_max, &lat_min, &lat_max)) {

            //printf("\n");
            //printf("DATE %2d %2d %2d\n", pix->dt[0][k], pix->dt[1][k], pix->dt[2][k]);
            //printf("TIME %2d %2d %2d\n", pix->tm[0][k], pix->tm[1][k], pix->tm[2][k]);

            double dtk = decimal_date_time(&pix->dt[0][k], &pix->dt[1][k], &pix->dt[2][k],
                                           &pix->tm[0][k], &pix->tm[1][k], &pix->tm[2][k]);

            // Accumulate aot values in the range [0.0 - 1.0] obtained within
            // +-12 hours of the scene datetime.

            float aot_value_k = pix->aot[k];

            if (aot_value_k > 0.0 && aot_value_k <= 1.0 && fabs(dtk - scene_dt) < SCENE_DT_TOLERANCE) {
                aot_sum += aot_value_k;
                ++aot_count;
                if (print_values)
                    //printf("@ %12.6f (%12.6f, %12.6f) (%4d-%02.2d-%02.2d %02.2d:%02.2d:%02.2d) %f %f %f\n",
                    //       aot_value_k, lon, lat, year, month, day, hour, min, sec, scene_dt, dtk, fabs(dtk - scene_dt));
                    printf("@  %12.6f (%12.6f, %12.6f) %f %f %f\n",
                           aot_value_k, lon, lat, scene_dt, dtk, fabs(dtk - scene_dt));
            } else {
                if (print_values)
                    printf("  %12.6f (%12.6f, %12.6f) %f %f %f\n",
                           aot_value_k, lon, lat, scene_dt, dtk, fabs(dtk - scene_dt));
            }

            ++in_scene_count;
        }
    }

    float aot_mean = 0;

    if (aot_count)
        aot_mean = aot_sum / aot_count;

    if (0) {
        printf("\n");
        printf("in_scene_count %d\n", in_scene_count);
        printf("aot_count %d\n", aot_count);
        printf("aot_mean %f\n", aot_mean);
    }

    return aot_mean;
}


float eval_cmp_aot(CMPDATA *cmp, float lon_min, float lon_max, float lat_min, float lat_max)
{
    int i, j;
    float aot_sum = 0;
    int aot_count = 0;

    assert(cmp);

    for (j=0; j<cmp->ny; j++) {
        //float yj = cmp->y[j];
        float lat = cmp->y[j];
        for (i=0; i<cmp->nx; i++) {
            float lon = cmp->x[i];
            int k = j * cmp->nx + i;
            assert(lat >= -90.0 && lat <= 90.0);
            float aot_value_k = cmp->aot[k];

            if (aot_value_k > 0.0 && aot_value_k <= 1.0 &&
                point_in_region(&lon, &lat, &lon_min, &lon_max, &lat_min, &lat_max)) {
                aot_sum += aot_value_k;
                ++aot_count;
                if (print_values)
                    printf("@ %12.6f (%12.6f, %12.6f) (%5d, %5d, %5d)\n",
                           aot_value_k, lon, lat, i, j, k);
            } else {
                if (print_values)
                    printf("  %12.6f (%12.6f, %12.6f) (%5d, %5d, %5d)\n",
                           aot_value_k, lon, lat, i, j, k);
            }
        }
    }

    float aot_mean = 0;

    if (aot_count)
        aot_mean = aot_sum / aot_count;

    if (0) {
        printf("\n");
        printf("aot_count %d\n", aot_count);
        printf("aot_mean %f\n", aot_mean);
    }

    return aot_mean;
}



int test_dt_diff(void)
{
    double x = 20020727.008680;
    double y = 20020727.018680;

    printf("%f\n", x - y);
    printf("%f\n", fabs(x - y));

    return 0;
}


int test_pix(void)
{
    int i, j, k;
    int debug_print_cmp = TRUE;
    int debug_print_pix = FALSE;

    // ...todo...

    //PIX_DATA *pix = load_pix("/data/nbar-data/ancillary/aerosol/AASTR/ATSR_LF_200208.pix");
    //CMP_DATA *cmp = load_cmp("/data/nbar-data/ancillary/aerosol/AASTR/aot_mean_2002_08.cmp");

    PIXDATA *pix = load_pix("ATSR_LF_200207.pix");

    if (debug_print_pix) {
        for (k=0; k<pix->nrecs; k++)
            print_pix_entry(pix, k);
    }

    CMPDATA *cmp = load_cmp("aot_mean_Apr2003_AllAerosols.cmp");

    if (debug_print_cmp) {
        char *cmp_fmt = "%8d %8d  %f %f %f\n";
        for (j=0; j<cmp->ny; j++) {
            float yj = cmp->y[j];
            for (i=0; i<cmp->nx; i++) {
                int k = j * cmp->ny + i;
                printf(cmp_fmt, i, j, cmp->x[i], yj, cmp->aot[k]);
            }
        }
    }

    // dummy scene extents and datetime for testing

    float roi_lon[2] = { 129.5, 130.5 };
    float roi_lat[2] = { -16.5, -15.5 };

    short int dt_foo_year = 2002;
    short int dt_foo_month = 7;
    short int dt_foo_day = 27;
    short int dt_foo_hour = 0;
    short int dt_foo_min = 12;
    short int dt_foo_sec = 30;

    float test_aot_mean = eval_pix_aot(
                              pix,
                              roi_lon[0], roi_lon[1], roi_lat[0], roi_lat[1],
                              dt_foo_year, dt_foo_month, dt_foo_day, dt_foo_hour, dt_foo_min, dt_foo_sec
                          );

    printf("MEAN AOT VALUE: %f\n", test_aot_mean);

    return !(test_aot_mean > 0);
}


void usage(void)
{
    printf("usage: %s [--pix|--cmp] <data_file.pix|.cmp> \\ \n", exe_name);
    printf("                  --west <west_extent> \\ \n");
    printf("                  --east <east_extent> \\ \n");
    printf("                  --south <south_extent> \\ \n");
    printf("                  --north <north_extent> \\ \n");
    printf("                  [--date <scene_date: 'YYYY-MM-DD'>] \\ \n");
    printf("                  [--time <scene_time: 'HH:MM:SS.dddd'>] \n");
    printf("                  [--print-values] \n\n");
    printf("Argument '--pix' requires '--date <date>' and '--time <time>'.\n\n");
    printf("Argument '--print-values' prints all values loaded from the data file.\n");
    printf("Entries starting with '@' lie within the specified extents. These are\n");
    printf("used to calculate the mean.\n\n");
}


CMDLINE_DATA* parse_args(int argc, char *argv[])
{
    static float undefined_extent = 99999;

    static char *opt_spec = "w:e:n:s:c:p:ht:d:";
    static struct option long_options[] = {
        {"west",  required_argument, 0, 'w'},
        {"east",  required_argument, 0, 'e'},
        {"south", required_argument, 0, 's'},
        {"north", required_argument, 0, 'n'},
        {"cmp",   required_argument, 0, 'c'},
        {"pix",   required_argument, 0, 'p'},
        {"date",  required_argument, 0, 'd'},
        {"time",  required_argument, 0, 't'},
        {"usage", no_argument,       0, 'h'},
        {"print-values", no_argument, &print_values, 1},
        {0, 0, 0, 0}
    };

    CMDLINE_DATA *cmdline = (CMDLINE_DATA*) malloc(sizeof(CMDLINE_DATA));
    cmdline->east = cmdline->west = cmdline->south = cmdline->north = undefined_extent;
    cmdline->pix_file[0] = '\0';
    cmdline->cmp_file[0] = '\0';
    cmdline->date[0] = '\0';
    cmdline->time[0] = '\0';

    while (1) {
        int option_index = 0;
        int c = getopt_long(argc, argv, opt_spec, long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case 'p':
                strcpy(cmdline->pix_file, optarg);
                break;
            case 'c':
                strcpy(cmdline->cmp_file, optarg);
                break;
            case 'w':
                cmdline->west = FLOAT_ATOF(optarg);
                break;
            case 'e':
                cmdline->east = FLOAT_ATOF(optarg);
                break;
            case 's':
                cmdline->south = FLOAT_ATOF(optarg);
                break;
            case 'n':
                cmdline->north = FLOAT_ATOF(optarg);
                break;
            case 'd':
                strcpy(cmdline->date, optarg);
                break;
            case 't':
                strcpy(cmdline->time, optarg);
                break;
            case 'h':
                usage();
                exit(1);
            default:
                //printf ("?? getopt returned character code 0%o ??\n", c);
                break;
        }
    }

    // Error: --pix and --cmp.

    if (strlen(cmdline->pix_file) && strlen(cmdline->cmp_file)) {
        printf("ERROR: PIX (-p,--pix) and CMP (-c,--cmp) options cannot be used together\n");
        exit(1);
    }

    // Error: no --pix or --cmp.

    if (!strlen(cmdline->pix_file) && !strlen(cmdline->cmp_file)) {
        printf("ERROR: PIX (-p,--pix) or CMP (-c,--cmp) option must be used\n");
        exit(1);
    }

    // Error: undefined extent.

    if (cmdline->west  == undefined_extent || cmdline->east  == undefined_extent ||
        cmdline->south == undefined_extent || cmdline->north == undefined_extent) {
        printf("ERROR: extents are not fully defined\n");
        printf("    north = %f\n", cmdline->north);
        printf("    south = %f\n", cmdline->south);
        printf("    east  = %f\n", cmdline->east);
        printf("    west  = %f\n", cmdline->west);
        printf("\n");
        usage();
        exit(1);
    }

    // Error: --pix but no date and/or time.

    if ( strlen(cmdline->pix_file) && !(strlen(cmdline->date) && strlen(cmdline->time))) {
        printf("ERROR: PIX (--pix) file requires '--date <date>' and '--time <time>'\n");
        printf("\n");
        usage();
        exit(1);
    }

    // Additional args are ignored.

    if (optind < argc) {
        printf ("WARNING: ignored additional arguments: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
    }

    return cmdline;
}




int main(int argc, char *argv[])
{
    // static flag
    print_values = 0;

    float result = 0;
    void *data = 0;


    CMDLINE_DATA *cmdline = parse_args(argc, argv);

    if (1) {
        printf("\n");
        printf("WEST  %f\n", cmdline->west);
        printf("EAST  %f\n", cmdline->east);
        printf("SOUTH %f\n", cmdline->south);
        printf("NORTH %f\n", cmdline->north);
        printf("PIXFILE %s\n", cmdline->pix_file);
        printf("CMPFILE %s\n", cmdline->cmp_file);
        printf("DATE    %s\n", cmdline->date);
        printf("TIME    %s\n", cmdline->time);
    }

    if (strlen(cmdline->cmp_file)) {

        data = load_cmp(cmdline->cmp_file);
        result = eval_cmp_aot(data, cmdline->west, cmdline->east,
                                    cmdline->south, cmdline->north);

    } else if (strlen(cmdline->pix_file)) {

        // Parse date and time strings of the form:
        // DATE    "2011-08-30"
        // TIME    "10:31:57.123"

        char s_year[5], s_month[3], s_day[3];
        char s_hour[3], s_min[3], s_sec[3];

        int year = 0, month = 0, day = 0;
        int hour = 0, min = 0, sec = 0;

        // Date
        // C string ops are such fun :)

        strncpy(s_year, (cmdline->date + 0), 4);
        s_year[4] = '\0';
        strncpy(s_month, (cmdline->date + 5), 2);
        s_month[2] = '\0';
        strncpy(s_day, (cmdline->date + 8), 2);
        s_day[2] = '\0';

        sscanf(s_year, "%d", &year);
        sscanf(s_month, "%d", &month);
        sscanf(s_day, "%d", &day);

        //printf("[%s] [%s] [%s]\n", s_year, s_month, s_day);
        //printf("(%4d) (%2d) (%2d)\n", year, month, day);

        assert(year >= 1970);
        assert(month > 0 && month <= 12);
        assert(day > 0 && day <= 31);

        // Time string

        strncpy(s_hour, (cmdline->time + 0), 2);
        s_hour[2] = '\0';
        strncpy(s_min, (cmdline->time + 3), 2);
        s_min[2] = '\0';
        strncpy(s_sec, (cmdline->time + 6), 2);
        s_sec[2] = '\0';

        sscanf(s_hour, "%d", &hour);
        sscanf(s_min, "%d", &min);
        sscanf(s_sec, "%d", &sec);

        //printf("[%s] [%s] [%s]\n", s_hour, s_min, s_sec);
        //printf("(%2d) (%2d) (%2d)\n", hour, min, sec);

        assert(hour >= 0 && hour < 24);
        assert(min >= 0 && min < 60);
        assert(sec >= 0 && sec < 60);

        // Load the data file and extract value(s).

        data = load_pix(cmdline->pix_file);
        result = eval_pix_aot(data, cmdline->west, cmdline->east,
                                    cmdline->south, cmdline->north,
                                    (short int) year,
                                    (short int) month,
                                    (short int) day,
                                    (short int) hour,
                                    (short int) min,
                                    (short int) sec);

    }

    printf("AOT AATSR value: %f\n", result);

    return 0;
}

