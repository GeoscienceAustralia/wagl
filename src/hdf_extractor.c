/*******************************************************************
 * HdfExtractor.c: this program extracts a subset of image from
 * a input HDF file and then calculate the average value over the subset
 *
 * @author: Frank Q. Fu        03 February May 2010
 *******************************************************************/

// bin/hdfExtractor 139.083625, -29.183625 141.283375 -31.383375 inputs/MCD43A1.2009.049.aust.005.b01.500m_0620_0670nm_brdf_par_fiso.hdf

// bin/hdfExtractor 149.97167435, -35.47110925 150.06562035 -35.57914715 inputs/MCD43A1.2009.049.aust.005.b01.500m_0620_0670nm_brdf_par_fiso.hdf


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PLATFORM_VAYU)
#include "hdf4_netcdf.h"
#endif

#include "hdf.h"
#include "mfhdf.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

// Fixed Upper-left corner of all input HDF files
const double ULX = 110.0;
const double ULY = -10.0;
const double PIXEL_SIZE = 0.0046973;


// range of the subset
double ulX = 0.0, ulY = 0.0, lrX = 0.0, lrY = 0.0;

// corresponding pixel and line numbers of the subset
int x0 = 0, y0 = 0, x1 = 0, y1 = 0;

int16 fillValue;

int verbose;

const char *getDataTypeName (int type);
double getAvgValue(char *hdf);
void trim (char *string);  // copied from bin2jpg.c


int main(int argc, char *argv[])
{
    void usage();
    char *inHdf = NULL;
    verbose = FALSE;

    if ( argc >= 5 ) {
        ulX = atof(argv[1]);
        ulY = atof(argv[2]);
        lrX = atof(argv[3]);
        lrY = atof(argv[4]);
        inHdf = argv[5];
    } else {
        usage();
    }

    // Check UL coordinates against BRDF HDF UL coordinates.

    if (ulY > ULY || ulX < ULX) {
        fprintf(stdout,
                "\nhdfExtractor ERROR: Scene UL (%f, %f), BRDF UL (%f, %f)\n\n",
                ulX, ulY, ULX, ULY);
        fprintf(stderr,
                "\nhdfExtractor ERROR: Scene UL (%f, %f), BRDF UL (%f, %f)\n\n",
                ulX, ulY, ULX, ULY);
        exit(99);
    }

    double avgValue = getAvgValue (inHdf);
    printf("Avg Value=%f\n", avgValue); // pass value to python
    return 0;
}


void usage(){
    fprintf (stdout, "Usage: hdfExtractor <ulX> <ulY> <lrX> <lrY> <input_HDF>\n");
    fprintf (stdout, "       <ulX> <ulY>: coodinates of Upper-left corner\n");
    fprintf (stdout, "       <lrX> <lrY>: coodinates of Lower-right corner\n");
    fprintf (stdout, "       <input_HDF>: input HDF file name\n");
    //fprintf (stdout, "       [dataset_no]: number of dataset to be extracted\n");
    exit(11);
}

// -L/hdf_c_4.2r1-linux/lib -lmfhdf -ldf -ljpeg -lsz -lz -L/usr/lib -lnsl
// Extract metadata and metadata from the input HDF file
double getAvgValue(char *hdf)
{
    int contains (const char *s, const char *w);
    int isInSdsList(char *, int, char **);
    double parse2dSds(int32 hdfID, int index);

    int32 hdfID, numOfDataSets, numOfFileAttrs, index;
    int32 istat;

    //int32 sdsID, numOfAttrs, dataType, numOfDims;
    //int32 dims[MAX_VAR_DIMS];
    //char name[64];

    /* Open the file and initiate the SD interface. */
    hdfID = SDstart(hdf, DFACC_READ);
    if (hdfID == FAIL)
    {
        fprintf (stderr, " Error: cannot open input HDF file: %s\n", hdf);
        return 0;
    }

    if (verbose)
        printf ("Open input HDF file: \"%s\" Done.\n", hdf);

    /* Determine the contents of the file. */
    istat = SDfileinfo(hdfID, &numOfDataSets, &numOfFileAttrs);
    if (verbose)
        fprintf (stdout, "Number of Datasets=%ld, Number of File Attributes=%ld\n", numOfDataSets, numOfFileAttrs);
/*
    // Access and print the name of every dataset in the file.
    for (index = 0; index < numOfDataSets; index++)
    {
        sdsID = SDselect(hdfID, index);        // data set ID
        istat = SDgetinfo(sdsID, name, &numOfDims, dims, &dataType, &numOfAttrs);

        if (verbose)
            fprintf (stdout, "  Extracting Dataset %ld: \"%s\" done.\n", (index+1), name);

        // we need to extract data from 1D dataset first
        if (numOfDims == 1)
        {
            //fprintf (stdout, "%s %ld %ld \n", getDataTypeName(dataType), numOfDims, dims[0]);
            if (verbose)
                fprintf (stdout, "     1D Dataset(%s): Width=%ld\n",
                getDataTypeName(dataType), dims[0]);

            float32 float32_data[dims[0]*sizeof(float32)];
            int32 start[1];
            start[0] = 0;
            int32 edges[1];
            edges[0] = dims[0];
            istat = SDreaddata (sdsID, start, NULL, edges, float32_data); //(VOIDP)

            int i;
            for (i = 0; i < dims[0]; i ++) {
                //fprintf (stdout, "(%d)%0.5f ", i+1, float32_data[i]);
                if (contains(name, "latitude")){
                    if (abs(ulY-float32_data[i]) < 0.00001)
                        y0 = i;
                    if (abs(lrY-float32_data[i]) < 0.00001)
                        y1 = i;
                } else if (contains(name, "longitude")){
                    if (abs(ulX-float32_data[i]) < 0.00001)
                        x0 = i;
                    if (abs(lrX-float32_data[i]) < 0.00001)
                        x1 = i;
                }
            }
            fprintf (stdout, "\n");
            fprintf (stdout, "UL(%d, %d) LR(%d, %d)\n", x0, y0, x1, y1);
        }

        istat = SDendaccess(sdsID);
        if (istat == FAIL)
            fprintf (stderr, "     Error: cannot end access to sdsID(%ld): %ld\n", sdsID, istat);
    }
*/

    x0 = (int)((ulX - ULX)/PIXEL_SIZE+0.5);
    y0 = (int)((ULY - ulY)/PIXEL_SIZE+0.5);
    x1 = (int)((lrX - ULX)/PIXEL_SIZE+0.5);
    y1 = (int)((ULY - lrY)/PIXEL_SIZE+0.5);

    index = 0;  // the first SDS is real data
    double avg = parse2dSds(hdfID, index);

    // Terminate access to the SD interface and close the file
    istat = SDend(hdfID);
    if (istat == FAIL)
    {
        fprintf (stderr, "Error: cannot close input HDF file\n");
        return 0;
    }
    return avg;

}

double parse2dSds(int32 hdfID, int index){

    double avgValue = 0.0;
    char name[64];
    int32 dims[MAX_VAR_DIMS];
    int32 numOfDims, dataType, numOfAttrs, istat;

    int32 sdsID = SDselect(hdfID, index);        // data set ID
    istat = SDgetinfo(sdsID, name, &numOfDims, dims, &dataType, &numOfAttrs);

    if (numOfDims != 2)
        return 0.0;

    float64 scale_factor [1];        // 0.004789327298677885
    float64 scale_factor_err [1];
    float64 add_offset [1];            // -53833.56628700022
    float64 add_offset_err [1];
    int32 nt [1];
    SDgetcal (sdsID, scale_factor, scale_factor_err, add_offset, add_offset_err, nt);

    SDgetfillvalue(sdsID, &fillValue);

    // data_type dimension width height ...
    /*
    if (verbose)
    {
        fprintf (stdout, "%s %ld %ld %ld %1.12f %f %6.12f %f\n",
        getDataTypeName(dataType), numOfDims, dims[0], dims[1],
        scale_factor[0], scale_factor_err[0], add_offset[0], add_offset_err[0]);
    }
    */

    if (verbose)
    {
        fprintf (stdout, "     2D Data Set(%s): Height=%ld; Width=%ld\n",
            getDataTypeName(dataType), dims[0], dims[1]);
        printf ("     Fill Value=%d\n", fillValue);  // 0 or -32768
        fprintf (stdout, "     Scale_Factor=%1.21f, Scale_Factor_Error=%f\n",
            scale_factor[0], scale_factor_err[0]);
        fprintf (stdout, "     Add_Offset=%6.21f, Add_Offset_Error=%f\n",
            add_offset[0], add_offset_err[0]);
    }

    // copy all data
    int16 *int16_data = (int16 *)malloc(dims[0]*dims[1]*sizeof(int16));
    if (int16_data == (int16 *)NULL)
    {
        fprintf (stderr, "     Error: cannot allocate memory to int16 data\n");
        return 0;
    }

    int32 start[2];
    start[0] = 0;
    start[1] = 0;
    int32 edges[2];
    edges[0] = dims[0];
    edges[1] = dims[1];
    istat = SDreaddata (sdsID, start, NULL, edges, int16_data);


    if (verbose)
    {
        fprintf (stdout, "input window UL(%f, %f) LR(%f, %f)\n", ulX, ulY, lrX, lrY);
        fprintf (stdout, "x0=%d, y0=%d; x1=%d, y1=%d\n", x0, y0, x1, y1);
    }

    // calculate
    long numOfValidPixels = 0, numOfNullPixels = 0;
    double totalValue = 0.0;
    int x, y;
    float32 value = 0;
    for (y = y0; y <= y1; y ++)
    {
        for (x = x0; x <= x1; x ++)
        {
            int pixelNo = y*dims[1]+x;

            /*if ((x==6174 && y==5345) || (x==6599 && y==4470))
                fprintf (stdout, "Value(%d, %d)=%d (%f)\n",
                x, y,int16_data[pixelNo],
                (int16_data[pixelNo]-add_offset[0])*scale_factor[0]);
            */

            // ignore fillValue (-32768)
            if (int16_data[pixelNo] == fillValue)
            {
                numOfNullPixels ++;
                continue;
            }
            else
            {
                // convert 16-bit integers into float values
                // Equation: floatValue = (int16 - add_offset)*scale_factor
                value = (int16_data[pixelNo]-add_offset[0])*scale_factor[0];
                totalValue += value;
                numOfValidPixels ++;
            }
        }
    }

    if (numOfValidPixels != 0)
        avgValue = totalValue/numOfValidPixels;

    if (verbose) {
        fprintf (stdout, " Total Value=%f\n Number Of Valid Pixels=%ld\n Avg Value=%f\n",
        totalValue, numOfValidPixels, avgValue);
        fprintf (stdout, " Number Of Null Pixels(%d)=%ld\n", fillValue, numOfNullPixels);
    }

    if (int16_data != (int16 *)NULL)
        free (int16_data);

    return avgValue;
}


int isInSdsList (char *sdsName, int numOfDatasets, char ** sdsList)
{
    int contains (const char *s, const char *w);

    int i;
    for (i = 0; i < numOfDatasets; i ++)
    {
        if (contains(sdsName, sdsList[i]) || contains(sdsList[i], sdsName))
            return 1;
    }
    return 0;
}

// Return a name string representing the data type
const char *getDataTypeName (int type)
{
    switch (type)
    {
        case DFNT_UCHAR:            // 3
            return "UCHAR";
        case DFNT_CHAR:                // 4
            return "CHAR";
        case DFNT_FLOAT32:            // 5
            return "FLOAT32";
        case DFNT_FLOAT64:            // 6
            return "FLOAT64";
        case DFNT_FLOAT128:            // 7
            return "FLOAT128";
        case DFNT_INT8:                // 20
            return "INT8";
        case DFNT_UINT8:            // 21
            return "UINT8";
        case DFNT_INT16:            // 22
            return "INT16";
        case DFNT_UINT16:            // 23
            return "UINT16";
        case DFNT_INT32:            // 24
            return "INT32";
        case DFNT_UINT32:            // 25
            return "UINT32";
        case DFNT_INT64:            // 26
            return "INT64";
        case DFNT_UINT64:            // 27
            return "UINT64";
        case DFNT_INT128:            // 28
            return "INT128";
        case DFNT_UINT128:            // 30
            return "UINT128";
        default:
            return "UNKNOWN_DATA_TYPE"; // -1 or 0 ???
    }
}


// return TRUE if string s contains string w ignoring case
int contains (const char *s, const char *w)
{
    char chs, chw;
    // TRUE for null string case
    if (s == NULL || strlen(s) == 0)
        return FALSE;
    else if (w == NULL || strlen(w) == 0)
        return TRUE;
    else if (strlen(s) < strlen(w))
        return FALSE;

    int i, j;
    for (i = 0; i < strlen(s)-strlen(w)+1; i ++)
    {
        chs = tolower (s[i]);
        chw = tolower (w[0]);
        if (chs == chw) // match the first letter
        {
            for (j = 1; j < strlen(w); j ++)
            {
                if (i+j >= strlen(s))
                    break;
                chs = tolower (s[i+j]);
                chw = tolower (w[j]);
                if (chs != chw)
                    break;
            }
            if (j == strlen(w))
                return TRUE;
        }
    }
    return FALSE;  // no match
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





