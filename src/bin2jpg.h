
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "JpegUtil.h"

// 07889_05
// ./bin/bin2jpg  ./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31_b7.bin,./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31_b4.bin,./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31_b1.bin 9200 9200 ./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31.jpg

/* --------------------------------------------------------
    Constants used by the system
   -------------------------------------------------------
*/

#ifndef int8
    #define int8 signed char
#endif

#ifndef uint8
    #define uint8 unsigned char
#endif

#ifndef int16
    #define int16 signed short
#endif


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


// Minimum and Maximum pixel values of a JPEG image
#define JPEG_MIN_VALUE 0
#define JPEG_MAX_VALUE  255

#define JPEG_QUALITY 90

#define INT16_NULL_VALUE -999
#define INT16_MAX_VALUE 32768
#define jpgWidth  1000
#define jpgHeight 1000


const char *usageStr = "Usage: bin2jpg <Input File List> <Width> <Height> <Output File Name> ";

/* --------------------------------------------------------
    Structures used by the system
	-------------------------------------------------------
*/

typedef struct {
	int x, y, z;  // coordinates of the a pixel
} Pixel;


// struct for double 2D point
typedef struct {
	int8 min, max;   // geographic location in meter
} RangeInt8;


// A simplified GeoTIFF files object
typedef struct
{
	FILE *fp;		// INFILE file pointer
	int byteOrder;		// big or little endian
	int width, height;    // size of the GeoTIFF

	// scanline used to buffer pixel values
	int16 **scanlines;
    int scanlineSize;
	int windowSize;
	double scaleX, scaleY;
	int bitsPerSample;
} INFILE;


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
char *toTimeString (long int time);

// parse input folder string into folder names delimited by comma ','
void parseInputFileList (char *fileList);

// Find out the range of the input Gtif
// void findDataValueRanges ();

// Remove extra white spaces in the front and/or end of the string
void trim (char *string);

// Returns a GeoTIFF file from all files contained in the folder
//  based on keywords
extern char *lookupFileName (char **file_list, int numOfFiles, const char **keywords);

// Find out the valid ranges of each input GeoTIFF
void findDataValueRanges ();

// Composite three int16 binary files into a RGB JPEG image of scaled size
void int16ToRgbJpeg ();

void readLinesToProcessingWindow(INFILE *inFileObjList, int winSizeY);

uint8 getAverageValue(INFILE infile, int winSizeY, int start, int end);

boolean setupAllInputFiles (INFILE *inFileObjList);
boolean setupFile (INFILE *int16File, char *fileName);
void closeFile (INFILE *int16File);

void printImageRange(uint8 *image_buffer, long numOfPixels);

// stretch three bands of bufffer image to enhance the contrast
// @param image_buffer: a array of pixels buffered for the 3 bands image
// @param totalPixels: total number of pixels
// @param threshold: cut-off threshold, half of which will be cut off over each side
void stretchThreeBands(uint8 *image_buffer, long totalPixels, double threshold);

