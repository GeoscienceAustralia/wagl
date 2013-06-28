/* -----------------------------------------------------------------
 * MODULE:
 *	High level routines to output image in JPEG format
 * RCS:
 *	$Id: JpegUtil.h
 * DESTINATION:
 *	dpfpr4
 * DESCRIPTION:
 *	Routines to convert
 *      - createJpegImage()
 *      - image_to_jpeg_memory()
 * ---------------------------------------------------------------- */

#ifndef _STDIO_H
#include <stdio.h>
#endif

#ifndef _JPEG_H
#define _JPEG_H

    #define MAX_IMAGE_HEIGHT 10240

    #include "jpeglib.h"
    GLOBAL(void) jpeg_memory_dest (j_compress_ptr cinfo, void * outbuf);
    GLOBAL(size_t) jpeg_memory_size (j_compress_ptr cinfo);

    #if 0

        /* open input file and allocate image construction space */
        int jpeg_source(char *fname, int width, int height, int maximum_depth);
        int jpeg_source_mem(char *buffer, int width, int height, int maximum_depth);

        /* close input file and free image construction space */
        void jpeg_close();

        /* black fill image construction space */
        void jpeg_new();

        /* add "count" bytes starting at "which" in file to "where" line of image */
        /* returns 0 on failure, 1 on success */
        int jpeg_add(int which, int where, int count);
        int jpeg_add_mem(int which, int where, int count);

        /* write constructed jpeg*/
        int jpeg_write(char *outfilename, int jpeg_depth);

    #endif

    /* One step mode: to write image (24bits RGB) data to a jpeg file */
    int createJpegImage(char *outfilename, char *img_buf, int quality,
                   int width, int height, int depth);

    /* One step mode: compress image (24bits RGB) data and put into memory */
    /* return the compressed size, (-1) means failed */
    int image_to_jpeg_memory(char *dest, char *img_buf, int quality,
                 int width, int height, int depth);

#endif
