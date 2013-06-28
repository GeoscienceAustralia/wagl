#ifndef TO_JPEG
#define TO_JPEG


#include "geo_tiffp.h"
#include "geo_keyp.h"
#include "JpegUtil.h"


/* --------------------------------------------------------
    Constants used by the system
   -------------------------------------------------------
*/

// boolean values for convenience
// typedef short int boolean;

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

// null/Minimum and maximum pixel values of a 16-bit GeoTIFF images
#define GTIF_MAX_VALUE  1023
#define GTIF_NULL_VALUE 0    // -1024

// Minimum and Maximum pixel values of a JPEG image
#define JPEG_MIN_VALUE 0
#define JPEG_MAX_VALUE  255

// maximum size of a GeoTIFF image with pixel size of 500m
// #define WIDTH  9000
// #define HEIGHT 7000

#define JPEG_QUALITY 90

/*
    R -120 4500 - band7
    G -63 4432 - band2
    B -68 2400 - band4

    For example for BAND7
    Anything smaller than 0 should set to 0,
    Anything greater than 4500 should set to 255
    Anything between 0~4500 should scale to 0~255
*/
#define NULL_VALUE 0
#define GTIF_MAX_VALUE 1023
//int band_ranges[] = {1023, 1023, 1023};
#define jpgWidth  1000 //inGtifs[0].width/buffLineNo;
#define jpgHeight 1000 //inGtifs[0].height/buffLineNo;


#endif
/* --------------------------------------------------------
    Structures used by the system
	-------------------------------------------------------
*/

typedef struct {
	int x, y, z;  // coordinates of the a pixel
} Pixel;


// struct for double 2D point
typedef struct {
	int16 min, max;   // geographic location in meter
} RangeInt16;


// A simplified GeoTIFF files object
typedef struct
{
  TIFF *tif;		// TIFF file pointer
  int byteOrder;		// big or little endian
  int width, height;    // size of the GeoTIFF

  // scanline used to buffer pixel values
  uint16 **scanlines;
  int windowSize;
  int scanlineSize;
  double scaleX, scaleY;
  double *tiePoints;    // Tie Points of the image in meters
  double *pixelSize;   // pixel scale in meter: 500m, 500m or 1km
  double UL[3];        // upper left corner of the image
  double UR[3];        // upper right corner of the image
  double LL[3];        // lower left corner of the image
  double LR[3];        // lower right corner of the image
  int lineNo;   // the number of  scanline buffered

  uint16 samplesPerPixel;
  uint16 bitsPerSample;

  int numBands;

} Gtif;

/*
// A simplified GeoTIFF files object
typedef struct
{
	FILE *fp;		// Fast file pointer
	int byteOrder;		// big or little endian
	int width, height;    // size of the GeoTIFF

	// scanline used to buffer pixel values
	uint16 **scanlines;
	int scanlineSize;
	int windowSize;
	float scaleX, scaleY;
	int bitsPerSample;
} Fast;
*/

/*---------------------------------------------------------
   Global variables to be used by the system
   --------------------------------------------------------
*/

boolean verbose;		// boolean for whether to print run-time info
char *outJpgFilename;     // output folder for storing resultant iamges
char *inFileNames[3]; // input folders containing input GeoTIFF files

int numOfInFiles;	// number of input GeoTIFF files in each folder

// Variables for the CPU time when running the program
// CPU time (in seconds) at the beginning of a program
time_t startTimer;

// A timer recording previously temporary CPU time
time_t prevTimer;


/* --------------------------------------------------------
    All extern functions declared in each header file
   --------------------------------------------------------
*/
// print the usage information about the program and then exit
void usage();

// set the current CPU timer
void setCPUTimer ();

// Convert the time in milli-seconds into a formatted string
// of "hours:minutes:seconds" which is stored in timeString
extern char *toTimeString (long int time);

// parse input folder string into folder names delimited by comma ','
void parseInputFileList (char *fileList);

// return the value of unsigned 8-bit integer pixel from a gtif GeoTIFF after scaling it to 0-255
uint8 getScaledUint8Value (Gtif *gtif, int x, int y);

// Find out the range of the input Gtif
// void findDataValueRanges ();

// Run the one-day composite algorithm and output composite products
extern void compositeRgbJpeg ();

// Return true if the filename contains ".tif"
boolean isGeoTIFF (char *filename);

// Remove extra white spaces in the front and/or end of the string
extern void trim (char *string);

// ------------- signed 16-bit vs unsigned 8-bit ---------------
// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order
extern int16 uint8ToInt16 (int byteOrder, uint8 a, uint8 b);

// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order. In Gtif, all bytes (char) are unsigned
// 8-bit integers
void int16ToUint8 (int byteOrder, uint8 *bytes, int16 data);

// ------------- signed 16 bits vs signed 8 bits ---------------
// convert 2 signed 8-bit integers into a signed 16-bit integer
// using given byte order
extern int16 int8ToInt16 (int byteOrder, int8 a, int8 b);

// convert 2 bytes into a int16 integer assuming big-endian
extern void int16ToInt8 (int byteOrder, int8 *bytes, int16 data);

// ------------- unsigned 16 bits vs unsigned 8 bits -------------
// convert an unsigned 16-bit integer into two unsigned 8-bit integers using given byte order
extern void uint16ToUint8 (int byteOrder, uint8 *bytes, uint16 data);

// convert two unsigned 8-bit integers into an unsigned 16-bit integer in big endian
uint16 uint8ToUint16 (int byteOrder, uint8 a, uint8 b);

// Returns a GeoTIFF file from all files contained in the folder
//  based on keywords
extern char *lookupFileName (char **file_list, int numOfFiles, const char **keywords);

// set up all fileNo-th GeoTIFF files as inputs
boolean setupAllGtifs (Gtif *inFiles);

// Set up an input Gtif file with the file name in the folder of folderName
boolean setupGtif (Gtif *gtif, char *fileName);

// release all memory allocated to the mask file before closing it
void closeGtif (Gtif *gtif);

// Find out the valid ranges of each input GeoTIFF
void findDataValueRanges ();

// return the pixel value of an unsigned 16-bit integer from a 1km file
extern uint16 getUint16PixelValue (Gtif *gtif, int x, int y);

// return one pixel value (int16) from a 1000m-pixelsize file
extern int16 getInt16PixelValue (Gtif *gtif, int x, int y);

/*
void fastToRgbJpeg ();
void readLinesToProcessingWindow(Fast *inFileObjList);
uint8 getAverageValue(Fast fast, int inPixelNo);
boolean setupAllInputFiles (Fast *inFileObjList);
boolean setupFile (Fast *fastFile, char *fileName);
void closeFile (Fast *fastFile);
*/


int16 readVapValue(Gtif *inGtif, int LineNo, int PixelNo, int BandNo);

uint8 getAverageValue(Gtif gtif, int winSizeY, int startPixelNo, int stopPixelNo);
