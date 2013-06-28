/****************************************************************
 * Contains all definitions used in the GtifMosaic Program
 *
 * @author: Frank Q. Fu,  Friday, 15/11/2008
 ****************************************************************
 */

#ifndef _STDIO_H
#include <stdio.h>
#endif

#include "jpeglib.h"

// libgeotiff
#include "geo_tiffp.h"
#include "geo_keyp.h"


/* --------------------------------------------------------
    Constants used by the system
   -------------------------------------------------------
*/

// boolean values for convenience
//typedef short int boolean;

#ifndef TRUE
    #define TRUE  1
#endif

#ifndef FALSE
    #define FALSE 0
#endif


// byte order of pixe value
// They may be already defined in "/usr/include/endian.h"
// TIFF_LITTLEENDIAN=18761
#ifndef LITTLE_ENDIAN
    #define LITTLE_ENDIAN TIFF_LITTLEENDIAN
#endif

// TIFF_BIGENDIAN=19789
#ifndef BIG_ENDIAN
    #define BIG_ENDIAN TIFF_BIGENDIAN
#endif

//
#ifndef PDP_ENDIAN
    #define PDP_ENDIAN 18765
#endif


#define NULL_VALUE_INT32 262144   // should be 4294967296L
#define NULL_VALUE  0  // for int16 -32768

// This is the number of input image files that can be supported for mosaic
// Currently it is set to 4, but can be changed to any number
#define MAX_IN_FILES 4

/* --------------------------------------------------------
    Data Structures used by the system
	-------------------------------------------------------
*/

typedef struct {
	int x, y, z;  // coordinates of the a pixel
} Pixel;


// struct for double 2D point
typedef struct {
	double x, y, z;   // geographic location in meter
} Point3d;


// An object of geotiff files
typedef struct
{
	TIFF *tif;		// TIFF file pointer
    char *filename;
	int byteOrder;		// byte order of big/little endian
	// int type;		// =1 for 1km, =2 for 500m, =4 for 500m
	boolean isValid;		// isValid
	int width, height;    // number of pixel in each direction

	// scanline used to buffer pixel values
	uint8 *scanline;
	int scanlineSize;
	int lineNo;   // the number of  scanline buffered

	double *pixelSize;   // pixel scale in meter: 500m, 500m or 1km
	double *tiePoints;    // Tie Points of the image in meter
	double UL[3];          // upper left corner of the image

	uint16 planarConf;
	uint16 samplesPerPixel;
	uint16 bitsPerSample;
} GeoTIFF;


/*---------------------------------------------------------
   Global variables to be used by the system
   --------------------------------------------------------
*/

boolean verbose;		// boolean for whether to print run-time info
char *outFileName;     // output folder for storing resultant iamges
char *inFileNames[MAX_IN_FILES]; // input folders containing input geotiff files
int numOfInFiles;	// number of input geotiff files in each folder

char *quadCode;

// maximum size of output image
int outWidth, outHeight;

double outTiePoint[6];
double outPixelSize;

// Variables for the CPU time when running the program
// CPU time (in seconds) at the beginning of a program
time_t startTimer;

/* --------------------------------------------------------
    All functions declared in the header file
   --------------------------------------------------------
*/

// Convert the time in milli-seconds into a formatted string
// of "hours:minutes:seconds" which is stored in timeString
char *toTimeString (long int time);

// Parse command line input
void parseCommandLine (int , char**);

// parse input folder string into folder names delimited by comma ','
void parseInputFileList (char *folderString);

void do2QuadrantsMosaic ();

// Do the composite of 4 input GeoTIFF images
void do4QuadrantsMosaic ();

// Return true if the filename contains ".tif"
boolean isGeoTIFF (char *filename);

// Remove extra white spaces in the front and/or end of the string
extern void trim (char *string);

// Return true if the full_str ends with the sub_str
boolean endswith (char *full_str, char *sub_str);

//Return true if file exists
boolean isFileExist(char *filename);

// ------------- signed 16/32/8 bits convert ---------------
// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order
extern int16 uint8ToInt16 (int byteOrder, uint8 a, uint8 b);

// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order. In GeoTIFF, all bytes (char) are unsigned 8-bit integers
void int16ToUint8 (int byteOrder, uint8 *bytes, int16 data);

// convert two unsigned 8-bit integers into an unsigned 16-bit integer in big endian
extern uint16 uint8ToUint16 (int byteOrder, uint8 a, uint8 b);

// convert four unsigned 8-bit integers into a signed 32-bit
// integer using given byte order
extern int32 uint8ToInt32 (int byteOrder, uint8 *bytes);

// convert a signed 32-bit integers into four unsigned 8-bit
// integers using given byte order
extern void int32ToUint8 (int byteOrder, uint8 *bytes, int32 data);

// set up all input files of the given file order
extern boolean setUpAllInputGtifs (GeoTIFF *inFiles);

extern void setOutputRange(GeoTIFF *inGtifList);

extern void setOutputRange_2quads(GeoTIFF *inGtifList);

// set up a new unsigned GeoTIFF file with given file name, type and bit number for output
extern boolean setUpOutputGtif (GeoTIFF  *outGtif, char *outFileName, GeoTIFF *inGtif);

// release all memory allocated to the mask file before closing it
extern void closeGeoTIFF (GeoTIFF *gtif);

// Return the Upper Left position of the pixel (x, y) in the output unsigned GeoTIFF
// Allowing an error of half pixel size
Point3d getOutPixelUL (GeoTIFF outGtif, int x, int y);

// Return a pixel defined by the range of outUL and outLR in the input file
Pixel getInPixel (GeoTIFF *inGtif, Point3d outUL);

// return one pixel value (int16) from a 1000m-pixelsize file
int16 getInt16PixelValue (GeoTIFF *gtif, int x, int y);

// Return the pixelNo-th pixel value of 16-bit integer from the scanline
// int32 getInt32PixelValue (GeoTIFF *gtif, int x, int y);
