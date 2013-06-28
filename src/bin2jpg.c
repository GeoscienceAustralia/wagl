/*****************************************************************
 * bin2jpg: composite three 16 bits binary files to a 24 bits RGB JPEG image
 *
 * @author: Frank Fu  Created on: 16/02/2010
 *
 *****************************************************************
 */

// ./bin/bin2jpg  ./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31_b7.bin,./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31_b4.bin,./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31_b1.bin 9200 9200 ./tmp/LS5_TM_OTH_P54_NBAR_092_084_2009_01_31.jpg


#include <time.h>
#include "bin2jpg.h"

int inWidth, inHeight;
uint8 maxValue = 0, minValue = 255;

int main (int argc, char** argv)
{
	verbose = FALSE; // TRUE;

    // Variables for the CPU time when running the program
    setCPUTimer ();    // set the first stage time
    startTimer = prevTimer;

    // Parse command line inputs including arguments and options
	if (argc != 5)
		usage ();


    parseInputFileList (argv[1]);
	if (numOfInFiles != 3)
	{
		fprintf (stdout, "Error: number of input INFILE-Format files must be 3.\n");
		usage ();
	}

    inWidth = atoi(argv[2]);
    inHeight = atoi(argv[3]);

    outJpgFilename = (char *)malloc(strlen(argv[4])*sizeof(char)+2);
    strcpy (outJpgFilename, argv[4]);

	// print command line inputs and options
	int i;
	if (verbose)
	{
		fprintf (stdout, "    Input Files: \n");
		for (i = 0; i < numOfInFiles; i ++)
			fprintf (stdout, "      (%d) \"%s\"\n", (i+1), inFileNames[i]);

		fprintf (stdout, "    Input Files size (%d*%d)\n", inWidth, inHeight);
		fprintf (stdout, "    Output File Name: \"%s\"\n", outJpgFilename);
	}

    // Find out the valid ranges of each input GeoTIFF
    // findDataValueRanges ();
    int16ToRgbJpeg ();

    // count CPU time for the whole processing
    fprintf (stdout, "\n    Total CPU time for the whole composite processing: %s.\n",
        toTimeString(difftime(time(NULL), startTimer)));

    // Set exit status to be 0 on successful execution
    exit (0);
}

// Composite three int16 binary files into a RGB JPEG image of scaled size
void int16ToRgbJpeg ()
{
    int jpgDepth = numOfInFiles;  // = 3 for RGB color images

	uint8 *image_buffer = malloc(sizeof(uint8)*jpgWidth*jpgHeight*jpgDepth);
	if (image_buffer == (uint8 *)NULL)
	{
		fprintf (stderr, "!!! Error: cannot allocate memory to image_buffer.\n");
		exit(13);
	}

    // set up all input inGtifs
    INFILE inFileObjList[numOfInFiles];
	setupAllInputFiles (inFileObjList);

	register int inLineNo, inPixelNo, fileNo;
    int winSizeY = 0;
	long numOfPixels = 0L;
	for (inLineNo = 0; inLineNo < jpgHeight; inLineNo ++)
	{
        winSizeY = (int)(inFileObjList[0].scaleY*(inLineNo+1)+0.5) - (int)(inFileObjList[0].scaleY*inLineNo+0.5);
        readLinesToProcessingWindow(inFileObjList, winSizeY);
        int startNo = 0, stopNo = 0;
		for (inPixelNo = 0; inPixelNo < jpgWidth; inPixelNo ++)
		{
            startNo = stopNo; // (int)(inPixelNo*inFileObjList[0].scaleX+0.5);
            stopNo = (int)((inPixelNo+1)*inFileObjList[0].scaleX+0.5);
            // read RGB values from the three input files
			for (fileNo = 0; fileNo < numOfInFiles; fileNo ++)
			{
                image_buffer[numOfPixels] = getAverageValue(inFileObjList[fileNo], winSizeY, startNo, stopNo);
                if (image_buffer[numOfPixels] > maxValue)
                    maxValue = image_buffer[numOfPixels];
                if (image_buffer[numOfPixels] < minValue && image_buffer[numOfPixels] != 0)
                    minValue = image_buffer[numOfPixels];
                numOfPixels++;
			}
		}  // end of inPixelNo
	}  // end of inLineNo

    // Close all open GeoTIFF files before output JPEG image
	for (fileNo = 0; fileNo < numOfInFiles; fileNo ++)
        closeFile (&(inFileObjList[fileNo]));

    // stretch image to enhance the contrast
    double threshold = 0.02;  // 2%
    stretchThreeBands(image_buffer, numOfPixels, threshold);

    if (verbose)
        printImageRange(image_buffer, numOfPixels);

    //write the buffered line into JPG image file
    if (createJpegImage (outJpgFilename, image_buffer, JPEG_QUALITY, jpgWidth, jpgHeight, jpgDepth))
    {
		fprintf(stderr, "!!! Error in writing JPEG image: \"%s\"\n", outJpgFilename);
		free (image_buffer);
		exit(10);
	}

    free (image_buffer);
    exit(0);
}

// Print out the range of buffered image
void printImageRange(uint8 *image_buffer, long numOfPixels){

    int maxValue = 0, minValue = 256;
    int i;
    for (i = 0; i < numOfPixels; i ++){
        if (image_buffer[i] > maxValue)
            maxValue = image_buffer[i];
        if (image_buffer[i] < minValue && image_buffer[i] != 0)
            minValue = image_buffer[i];
    }

    printf("numOfPixels=%ld; minValue=%d; maxValue=%d\n",numOfPixels, minValue, maxValue);
}

// Read a line of pixels into the buffered widnow
void readLinesToProcessingWindow(INFILE *inFileObjList, int winSizeY)
{
    int fileNo, winLineNo;
    int i, numRead = 0;
    for (fileNo = 0; fileNo < numOfInFiles; fileNo ++)
    {
        INFILE infile = inFileObjList[fileNo];
        for (winLineNo = 0; winLineNo < winSizeY; winLineNo ++)
        {
            numRead = fread (infile.scanlines[winLineNo], sizeof(int16), infile.scanlineSize, infile.fp);
            if (numRead < infile.scanlineSize)
            {
                for (i = numRead; i < infile.scanlineSize; i ++)
                    infile.scanlines[winLineNo][i] = 0;
            }
        }
    }
}

uint8 getAverageValue(INFILE infile, int winSizeY, int startNo, int stopNo)
{
    int sum = 0;
    int validPixelNo = 0;
    register int winLineNo, winPixelNo;
    for (winLineNo = 0; winLineNo < winSizeY; winLineNo ++)
    {
        for (winPixelNo = startNo; winPixelNo < stopNo; winPixelNo ++)
        {
            int16 tmp = (int16)infile.scanlines[winLineNo][winPixelNo];
            if (tmp > 0)  // != JPEG_MIN_VALUE)
            {
                sum += tmp;
                validPixelNo ++;
            }
        }
    }

    if (validPixelNo == 0)
    {
        return 0;
    }
    // Return null if majority is null
    else if (validPixelNo < infile.scaleX*infile.windowSize/2.0)
    {
        return 0;
    }
    // otherwise return the valid average value
    else
    {
    	// convert it to 0-255
        double tmpValue = (double)sum/validPixelNo;
        uint8 avgValue = (uint8)(tmpValue/10000*255 + 0.5);
        return avgValue;
    }
}

// set up all GeoTIFF files with given scale window size
boolean setupAllInputFiles (INFILE *inFileObjList)
{
	char *fileName;

	int numOfFilesOpened = 0;
	register int fileNo;
	for (fileNo = 0; fileNo < numOfInFiles; fileNo ++)
	{
		fileName = inFileNames[fileNo];

		if (setupFile (&inFileObjList[fileNo], fileName))
			numOfFilesOpened ++;
	}

	if (verbose && numOfFilesOpened < numOfInFiles)
	{
        fprintf(stderr, "    Warning: only %d/%d of input files have been setup!\n",
			numOfFilesOpened,numOfInFiles);
    }

	return (numOfFilesOpened == numOfInFiles);
}

// Return true if having successfully opened a GeoTIFF file and extract relavent metadata.
boolean setupFile (INFILE *int16File, char *fileName)
{
	if (fileName == (char *)NULL)
		return FALSE;

	// Open the file, read the Gtif information, and print to stdout.
    int16File->fp = (FILE *)NULL;  /* TIFF-level descriptor */
	int16File->fp = fopen(fileName, "rb");
	if (!int16File->fp)
	{
		fprintf(stderr, "    Error: cannot open file \"%s\".\n", fileName);
		return FALSE;
	}

	// size of the GeoTIFF image
	int16File->width = inWidth;
	int16File->height = inHeight;

	int16File->bitsPerSample = 16;
	int16File->scaleX = (double)inWidth/jpgWidth;
	int16File->scaleY = (double)inHeight/jpgHeight;

    // window size should be the larger one between scaleX and scaleY
    int16File->windowSize = (int)((int16File->scaleX>int16File->scaleX ? int16File->scaleX :int16File->scaleY) + 1);

	int16File->scanlineSize = inWidth;
	int16File->byteOrder = LITTLE_ENDIAN;

	// allocate memory to all scanline buffers
    int16File->scanlines = (int16 **)malloc(int16File->scanlineSize*int16File->windowSize*sizeof(int16));
    int i;
    for (i = 0; i < int16File->windowSize; i ++)
	{
        int16File->scanlines[i] = (int16 *)malloc(int16File->scanlineSize*sizeof(int16));
        if (int16File->scanlines[i] == (int16 *)NULL)
        {
            fprintf (stderr, "    Error: cannot allocate memory to %d-th scanline.\n", (i+1));
            exit(13);
        }
    }

    if (verbose){
        printf("windowSize=%d\n", int16File->windowSize);
        printf("scaleX=%f; scaleY=%f\n", int16File->scaleX, int16File->scaleY);
        printf("scanlineSize=%d; bitsPerSample=%d\n", int16File->scanlineSize, int16File->bitsPerSample);
    }

	// successfully extracted all metadata
	return TRUE;
}

// close the opened getiff file and release its memory
void closeFile (INFILE *int16File)
{
	if (int16File == (INFILE *)NULL)
		return;

    int i;
	for (i = 0; i < int16File->windowSize; i ++)
    {
        if (int16File->scanlines[i] != (int16 *)NULL)
            free (int16File->scanlines[i]);
    }

    free (int16File->scanlines);

	if (int16File->fp != (FILE *)NULL)
		fclose (int16File->fp);
}


// print the usage information about the program and then exit
void usage()
{
	fprintf(stdout, "%s\n", usageStr);
	fprintf(stdout, "Where:\n");
	fprintf(stdout, "      <Input File List>: three input int16 binary files separated by comma \",\"\n");
    fprintf(stdout, "                         e.g. \"file1,file2,file3\"\n");
    fprintf(stdout, "      <width>: width of input images\n");
    fprintf(stdout, "      <height>: height of input images\n");
	fprintf(stdout, "      <Output File Name>: file name of the output color JPEG image\n");
	exit(12);
}


// parse input folder string into folder names delimited by comma ','
void parseInputFileList (char *fileList)
{
	// parse input folder names from command line
	char *word = strtok (fileList, ", ");
	while (word != (char *)NULL)
	{
		if (inFileNames[numOfInFiles] == (char *)NULL)
			inFileNames[numOfInFiles] = (char *)malloc (strlen(word)+1);
		if (inFileNames[numOfInFiles] == (char *)NULL)
		{
			fprintf(stderr, "    Error: cannot allocate memory to \"inputFileNames\".\n");
			exit (10);
		}

		trim (word);
		if (word == (char *)NULL)
			continue;

		strcpy (inFileNames[numOfInFiles++], word);
		word = strtok (NULL, ",");
	}
}

// set the CPU timer
void setCPUTimer ()
{
    prevTimer = time(NULL); // start to count CPU time
}

// get the CPU time duration in a specific stage
long getStageTime ()
{
    time_t thisTimer = time(NULL);
    double duration = difftime (thisTimer, prevTimer);

    prevTimer = thisTimer;

    return (long)duration;
}

// Convert the time in milli-seconds into a formatted string
// of "hours:minutes:seconds.milli-seconds"
// e.g., "001:38:24.567" and "000:00:00.000"
char *toTimeString (long int duration)
{
    int hours = duration / 3600;

    // total seconds less than a hour
    long total_seconds = duration % 3600;

    int minutes = total_seconds/60;
    int seconds = total_seconds%60;

    char *timeString = (char *)malloc(10*sizeof(char));
    if (hours == 0)
    {
        if (minutes == 0 && seconds == 0)
                sprintf (timeString, "< %d second", 1);
        else
            sprintf (timeString, "%dm:%ds", minutes, seconds);
    }
    else
        sprintf (timeString, "%dh:%dm:%ds", hours, minutes, seconds);

    return timeString;
}


// remove spare white spaces in the beginning or back of the string
void trim (char *string)
{
	if (string == (char *)NULL)
		return;

	int length = strlen (string);
	int i, first = 0, last = length;

	// count how many spaces in the front of the string
	for (i = 0; i < strlen (string); i ++)
	{
		if (string[i]==' ' || string[i]=='\t' || string[i]=='\n' || string[i]=='\r')
			first ++;
		else
			break;
	}

	// count how many spaces at the back of the string
	for (i = strlen(string)-1; i >= 0; i --)
	{
		if (string[i]==' ' || string[i]=='\t' || string[i]=='\n' || string[i]=='\r')
			last --;
		else
			break;
	}

	if (first == 0 && last == length) // no any extra spaces
		return ;
	else if (first >= length)   // empty/null string
	{
		string = (char *)NULL;
		return ;
	}

	char newString[last-first+1];

	for (i = first; i < last; i ++)
		newString[i-first] = string[i];
	newString[last] = '\0';  // add null char

	strcpy (string, newString);  // copy new string
}

// stretch three bands of bufffer image to enhance the contrast
// @param image_buffer: a array of pixels buffered for the 3 bands image
// @param totalPixels: total number of pixels
// @param threshold: cut-off threshold, half of which will be cut off over each side
void stretchThreeBands(uint8 *image_buffer, long totalPixels, double threshold){
    const int BAND_NO = 3;
    const int MAX_SIZE = 256;
    long histogram[MAX_SIZE];
    int start = 0, end = MAX_SIZE;
    long offset = 0L, sum = 0L;
    int bandNo, pixelNo, i, j;

    int numOfPixels = totalPixels/BAND_NO;

    // number of pixels to be skipped on both sides
    offset = (long)(numOfPixels*threshold/2.0);

    // stretch over each band based on the histogram
    for (bandNo=0; bandNo < BAND_NO; bandNo ++){

        // Build the histogram
        memset(histogram, 0, MAX_SIZE * sizeof(long));
        for( i = 0; i < numOfPixels; i ++ ){
            pixelNo = i*BAND_NO+bandNo;
            histogram[image_buffer[pixelNo]]++;
        }

        // find the start position in the histogram
        sum = 0L;
        for (i = 1; i < MAX_SIZE-1; i ++ ){
            sum += histogram[i];
            if (sum >= offset){
                start = i;
                break;
            }
        }

        // find the end position in the histogram
        sum = 0L;
        for (j = MAX_SIZE-2; j > 0; j --){
            sum += histogram[j];
            if (sum >= offset){
                end = j;
                break;
            }
        }

        if (verbose)
            printf("offset=%ld; start=%d; end=%d\n", offset, start, end);

        if (end <= start)
            break;

        // stretch to the new range (start, end)
        double newValue;
        uint8 oldValue;
        for( i = 0; i < numOfPixels; i ++ ){
            pixelNo = i*BAND_NO + bandNo;
            oldValue = image_buffer[pixelNo];
            if (oldValue <= start){
                image_buffer[pixelNo] = 0;
            } else if (oldValue >= end) {
                image_buffer[pixelNo] = 255;
            } else {
                newValue = (double)(oldValue-start)/(end-start)*255;
                image_buffer[pixelNo] = (uint8)(newValue+0.5);
            }
        }
    }
}


void stretch_off(uint8 *image_buffer, long numOfPixels){
    long histogram[256];
    double lut[256];
    long isize;
    int i, j, high_pos, low_pos;
    long low = -1, high = 256, flags, val, lower;
    double width;
    int min_flags = 30;

    /* Build the histogram */
    memset(histogram, 0, 256 * sizeof(long));
    for( i = 0; i < numOfPixels; i ++ )
        histogram[image_buffer[i]]++;

    //for( i = 0; i < 256; i++ )
    //    printf("histogram[%d]=%ld ", i, histogram[i]);

    /* Calculate the initial step size */
    isize = numOfPixels / 256;

    /* Keep iterating, using smaller and smaller step sizes until we get a
    usable number of flag points */
    flags = 0;
    while( flags < min_flags && isize > 1 ) {

        /* Preset the look up table */
        for( i = 0; i < 256; i++ )
            lut[i] = (double)i;

        /* Find the bottom of the histogram */
        low = -1;
        while( histogram[++low] == 0.0 )
            lut[low] = 0.0;

        /* Find the top of the histogram */
        high = 256;
        while( histogram[--high] == 0.0 )
            lut[high] = 255.0;
        high++;

        /* Check for a single spike or no data */
        if( high <= low )
            break;

        /* Flag the end of each piece with a -1, working from the ends into
        the middle */
        flags = 0;
        val = 0;
        low_pos = low - 1;
        high_pos = high;
        while( high_pos > low_pos ) {
            while( (val < isize) && (low_pos < high_pos) )
                val += histogram[++low_pos];
            lut[low_pos] = -1.0;
            flags++;
            val = 0;

            while( (val < isize) && (high_pos > low_pos) )
                val += histogram[--high_pos];
            lut[high_pos] = -1.0;
            flags++;
            val = 0;
        }

        isize /= 2;
    }

    printf("%3ld    %6ld\n", flags, isize);

    /* If we reach zero isize then the data is almost a single spike.
    Don't stretch it */
    if( isize == 0 ) {
        printf("Couldn't generate enough reference points. Not stretching!\n");
        return;
    }

    width = 255.0 / flags;

    for( i = 0; i < flags; i++ ) {
        lower = low;

        for( j = low + 1; j < high; j++ ) {
            low = j;

            if( lut[j] < 0.0 )
            break;
        }

        for( j = lower; j < low; j++ ) {
            lut[j] = width * (i + ((double)(j - lower) / (double)(low - lower)));
        }
    }
    lut[high - 1] = 255.0;

    /* Apply the look up table back to the buffer */
    for( i = 0; i < numOfPixels; i++ )
        image_buffer[i] = lut[image_buffer[i]] + 0.5;

    return ;
}

