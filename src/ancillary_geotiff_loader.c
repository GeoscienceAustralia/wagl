//
// Generic DEM, ozone and water vapour loader.
//
// Opens GeoTIFF and reads data value(s) for specified location (lon, lat)
// and date/time. Code contains hardwired extent and scale values for each
// data type, as defined in the original versions of extractDem, etc.
//
// Dependencies: This code utilizes libgeotiff.
//
// Note that geolocation is achieved by indexing arrays/scanlines of known
// size. Future versions should use geotiff dimension and extent information.
//
// Roger Edberg / NEO Science and Strategy Group (Aug 2011)
//

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define _GNU_SOURCE
#include <getopt.h>

#include "geo_tiffp.h"


#define DEBUG 0
#define STRSZ 256


//
// Globals
//

static char global_errmsg[256];



//
// Extent and scale
//

enum { DATA_TYPE_DEM=0, DATA_TYPE_O3, DATA_TYPE_H2O };

typedef struct _extent {
    // UL[0] = LON (decimal degrees)
    // UL[1] = LAT (decimal degrees)
    // UL[2] = ??

    float UL[3], LL[3], UR[3], LR[3];
} EXTENT;


typedef struct _scale {
    float x, y;
} SCALE;


EXTENT extents[3] = {
    //
    // Geobounds for data files.
    //

    { // dem
        { 108.00,   -8.00, 0.0 },       // UL
        { 108.00,  -48.00, 0.0 },       // LL
        { 157.999999999,  -8.00, 0.0 }, // UR
        { 157.999999999, -48.00, 0.0 }, // LR
    },

    { // ozone
        { 109.00,   -7.00, 0.0 }, // UL
        { 109.00,  -47.00, 0.0 }, // LL
        { 159.00,   -7.00, 0.0 }, // UR
        { 159.00,  -47.00, 0.0 }, // LR
    },

    { // water vapour
        {  -1.25,   91.25, 0.0 }, // UL
        {  -1.25,  -91.25, 0.0 }, // LL
        { 358.75,   91.25, 0.0 }, // UR
        { 358.75,  -91.25, 0.0 }, // LR
    }
};

SCALE scales[3] = {
    //
    // Scale factors for data types.
    //

    { 1.0, 1.0 }, // dem
    { 2.0, 2.0 }, // o3
    { 2.5, 2.5 }, // h2o
};


//
// GeoTIFF data
//

typedef struct _gtiff_struct {
    TIFF *tif;
    int byte_order;
    int width, height;
    int window_size;
    int scanline_size;
    double *tie_points;
    double *pixel_scale;
    double UL[3], UR[3], LL[3], LR[3];
    unsigned int samples_per_pixel;
    unsigned int bits_per_sample;
    int dtype;
} GTIFF_STRUCT;


GTIFF_STRUCT *load_image(char *file_name);

void close_image(GTIFF_STRUCT *p);


//
// Command line args
//

typedef struct _ARGS {
    int data_type;
    char file_name[STRSZ];
    char date_time[STRSZ];
    double lat, lon;
} ARGS;




void error_exit(char *msg, int action)
{
    // Simple error handler

    char *fmt = "ERROR: %s\n";

    if (msg && strlen(msg))
        fprintf(stderr, fmt, msg);
    else
        fprintf(stderr, fmt, "<unspecified>");

    if (action)
        exit(action);
}


void print_usage(void)
{
    // Print usage instructions

    printf(
        "USAGE: \n\n"
        "    ancillary_geotiff_loader [--elevation|--ozone|--water] \\ \n"
        "        --file <datafile.tif> \\ \n"
        "        --lat <latitude (DD.xxxx)> \\ \n"
        "        --lon <longiitude (DD.xxxx)> \\ \n"
        "        [--datetime <datetime (DOY:HH.xxxx)>]\n\n"

        "Option '--datetime' must be used with '--water'. \n"
    );

}


ARGS* parse_args(int argc, char *argv[])
{
    // Parse command line arguments

    static char *opt_spec = "eowf:x:y:d:";

    static struct option long_options[] = {
        {"elevation", no_argument, 0, 'e'},
        {"ozone",     no_argument, 0, 'o'},
        {"water",     no_argument, 0, 'w'},
        {"file",      required_argument, 0, 'f'},
        {"lon",       required_argument, 0, 'x'},
        {"lat",       required_argument, 0, 'y'},
        {"datetime",  required_argument, 0, 'd'},
        {"usage",     no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    ARGS *args = (ARGS*) malloc(sizeof(ARGS));

    float dummy_lat_lon = -999.999;

    args->data_type = -1;
    args->lat = dummy_lat_lon;
    args->lon = dummy_lat_lon;
    args->file_name[0] = '\0';
    args->date_time[0] = '\0';

    while (1) {
        int option_index = 0;
        int c = getopt_long(argc, argv, opt_spec, long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {
            case 'e':
                args->data_type = DATA_TYPE_DEM;
                break;
            case 'o':
                args->data_type = DATA_TYPE_O3;
                break;
            case 'w':
                args->data_type = DATA_TYPE_H2O;
                break;
            case 'f':
                strcpy(args->file_name, optarg);
                break;
            case 'd':
                strcpy(args->date_time, optarg);
                break;
            case 'x':
                args->lon = atof(optarg);
                break;
            case 'y':
                args->lat = atof(optarg);
                break;
            case 'h':
                print_usage();
                exit(1);
            default:
                //printf ("?? getopt returned character code 0%o ??\n", c);
                ;
        }
    }

    // Calculation type must be [0 (elevation), 1 (ozone), 2 (water vapour)].

    if (args->data_type < 0 || args->data_type > 2) {
        //sprintf(global_errmsg, "ERROR: Data type not specified. Use '--elevation', '--ozone' or '--water'\n", args->data_type);
        print_usage();
        error_exit(NULL, 1);
    }

    // File, latitide and longitude are required.

    if (!strlen(args->file_name)) {
        sprintf(global_errmsg, "Option '--file <filename>' is required\n");
        error_exit(global_errmsg, 1);
    }

    if (args->lon == dummy_lat_lon) {
        sprintf(global_errmsg, "Option '--lon <value>' is required\n");
        error_exit(global_errmsg, 1);
    }

    if (args->lat == dummy_lat_lon) {
        sprintf(global_errmsg, "Option '--lat <value>' is required\n");
        error_exit(global_errmsg, 1);
    }

    // Option '--water' requires '--datetime'.

    if (args->data_type == DATA_TYPE_H2O && !strlen(args->date_time) ) {
        sprintf(global_errmsg, "Option '--water' requires '--datetime'\n");
        error_exit(global_errmsg, 1);
    }

    //
    // TODO check additional bad arg combinations
    //

    return args;
}


GTIFF_STRUCT* load_image(char *filename)
{
    // Load GeoTIFF image/data file

    int i, n;
    unsigned int samples_per_pixel = 0;
    unsigned int bits_per_sample = 0;
    int scanline_size = 0;
    int line_number = -1;
    int window_size= 0;
    int byte_order = 0;
    int image_width = 0;
    int image_height = 0;
    int n_pixel_scale = 0;
    int n_tie_points = 0;
    double *pixel_scale;
    double *tie_points;

    assert(filename);

    TIFF *ptiff = TIFFOpen(filename, "r");

    if (!ptiff) {
        sprintf(global_errmsg, "TIFFOpen(%s) FAILED\n", filename);
        error_exit(global_errmsg, 1);
    }

    TIFFGetField(ptiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
    TIFFGetField(ptiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    TIFFGetField(ptiff, TIFFTAG_IMAGEWIDTH, &image_width);
    TIFFGetField(ptiff, TIFFTAG_IMAGELENGTH, &image_height);
    TIFFGetField(ptiff, GTIFF_PIXELSCALE, &n_pixel_scale, &pixel_scale);

    // TODO tie points are not used, remove?
    //TIFFGetField(ptiff, GTIFF_TIEPOINTS, &n_tie_points, &tie_points);

    // From original: still required???
    bits_per_sample = (bits_per_sample > 8 && bits_per_sample <= 16)
                          ? 16 : bits_per_sample;

    //printf("BITSPERSAMPLE %d\n", bits_per_sample);

    scanline_size = image_width * (bits_per_sample / 8);

    byte_order = LITTLE_ENDIAN;

    // Populate data struct.

    GTIFF_STRUCT *p = (GTIFF_STRUCT*) malloc(sizeof(GTIFF_STRUCT));

    p->tif = ptiff;
    p->byte_order = byte_order;
    p->samples_per_pixel = samples_per_pixel;
    p->bits_per_sample = bits_per_sample;
    p->width = image_width;
    p->height = image_height;
    p->scanline_size = scanline_size;
    p->window_size = window_size;

    p->pixel_scale = pixel_scale;
    p->tie_points = tie_points;

    return p;
}


void close_image(GTIFF_STRUCT *p)
{
    // Close image file

    if (p && p->tif)
        TIFFClose(p->tif);
}


uint16 get_dem_value(GTIFF_STRUCT *p, int line, int pixel)
{
    // Extract DEM data value from GeoTIFF

    uint16 value = 0;

    assert(p);
    assert(p->scanline_size > 0);

    uint16 *buffer = (uint16*) malloc(p->scanline_size);
    assert(buffer);

    if (TIFFReadScanline(p->tif, buffer, (uint32)line, 0) > 0)
        value = buffer[pixel];
    else {
        sprintf(global_errmsg, "TIFFReadScanline failed (DEM)\n");
        error_exit(global_errmsg, 1);
    }

    free(buffer);
    return value;
}


float get_ozone_value(GTIFF_STRUCT *p, int line, int pixel)
{
    // Extract ozone data value from GeoTIFF

    float value = 0;

    assert(p);
    assert(p->scanline_size > 0);

    float *buffer = (float*) malloc(p->scanline_size);
    assert(buffer);

    if (TIFFReadScanline(p->tif, buffer, (uint32)line, 0) > 0)
        value = buffer[pixel];
    else {
        sprintf(global_errmsg, "TIFFReadScanline failed (O3)\n");
        error_exit(global_errmsg, 1);
    }

    free(buffer);
    return value;
}


int16 get_watervapour_value(GTIFF_STRUCT *p,
                            int line, int pixel,
                            char *datetime, float scaley)
{
    // Extract water vapour data value from GeoTIFF
    // NOTE int16 instead of uint16 ???

    int err = 0;
    int16 value = 0;

    int day_of_year;
    float hour;

    assert(p);
    assert(p->scanline_size > 0);

    sscanf(datetime, "%d:%f", &day_of_year, &hour);

    int part_day = (int) ((hour + 3.0) / 6.0);
    int band = (day_of_year - 1) * 4 + part_day;

    //

    double sx = (double) p->width / 1000;
    double sy = (double) p->height / 1000;
    int wsize = (int) sx + 1;

    int scanline_size = TIFFScanlineSize(p->tif);
    int window_size = (int) scaley + 1;
    p->byte_order = LITTLE_ENDIAN;

    if (DEBUG) {
        printf("scanline_size %d\n", scanline_size);
        printf("window_size   %d\n", window_size);
        printf("wsize         %d\n", wsize);
    }

    int16 *buffer = (int16*) malloc(scanline_size);
    assert(buffer);

    int index = p->samples_per_pixel*pixel + band;

    if (DEBUG) {
        printf("malloc'ed buffer %d\n", scanline_size);
        printf("line %d\n", line);
        printf("pixel %d\n", pixel);
        printf("band %d\n", band);
        printf("buffer index %d\n", index);
    }

    if (TIFFReadScanline(p->tif, buffer, (uint32)line, 0) > 0)
        value = buffer[index];
    else {
        sprintf(global_errmsg, "TIFFReadScanline failed (H2O)\n");
        error_exit(global_errmsg, 1);
    }

    free(buffer);
    return value;
}


int main(int argc, char** argv)
{
    ARGS *args = parse_args(argc, argv);

    if (DEBUG) {
        printf("CMDLINE.data_type %d\n", args->data_type);
        printf("CMDLINE.file_name %s\n", args->file_name);
        printf("CMDLINE.lon %f\n", args->lon);
        printf("CMDLINE.lat %f\n", args->lat);
        printf("CMDLINE.date_time %s\n", args->date_time);
    }

    // Extent and scale factor for the data type.

    int data_type = args->data_type;
    assert(data_type >= 0 && data_type < 3);

    EXTENT ext = extents[data_type];

    if (args->lat < ext.LR[1] || args->lat > ext.UL[1] ||
        args->lon < ext.UL[0] || args->lon > ext.LR[0]) {
        sprintf(global_errmsg, "Point (%f, %f) is out of range",
                         args->lon, args->lat);
        error_exit(global_errmsg, 1);
    }

    SCALE scale = scales[data_type];
    assert(abs(scale.x) > 0 && abs(scale.y) > 0);

    // Open the file...

    GTIFF_STRUCT *p = load_image(args->file_name);

    // Get indices: line and pixel.

    int line  = (int) ((ext.UL[1] - args->lat) / scale.y);
    int pixel = (int) ((args->lon - ext.LL[0]) / scale.x);

    if (DEBUG) {
        printf("SCALE.X %f\n", scale.x);
        printf("SCALE.Y %f\n", scale.y);
        printf("LINE %d\n", line);
        printf("PIXEL %d\n", pixel);
        printf("WIDTH %d\n", p->width);
        printf("HEIGHT %d\n", p->height);
    }

    if (line < 0 || line > p->height || pixel< 0 || pixel > p->width) {
        sprintf(global_errmsg, "LINE, PIXEL (%d, %d) out of range {[0,%d], [0,%d]}\n",
                         line, pixel, p->height, p->width);
        error_exit(global_errmsg, 1);
    }

    // Extract the data value.

    int status = 0;

    switch (args->data_type) {
        case DATA_TYPE_DEM:
            fprintf(stdout, "Line: %d, Pixel: %d, DEM value: %d\n",
                            line, pixel, get_dem_value(p, line, pixel));
            break;
        case DATA_TYPE_O3:
            fprintf(stdout, "Line: %d, Pixel: %d, OZONE value: %f\n",
                            line, pixel, get_ozone_value(p, line, pixel));
            break;
        case DATA_TYPE_H2O:
            fprintf(stdout, "Line: %d, Pixel: %d, WATER VAPOUR value: %d\n",
                            line, pixel, get_watervapour_value(p, line, pixel, args->date_time, scale.y));
            break;
        default:
            assert(! "switch(data_type) fall-through case...should never happen");
    }

    close_image(p);
    return status;
}


