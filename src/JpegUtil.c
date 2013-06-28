/* ----------------------------------------------------------------
 * MODULE:
 *	High level routines to perform JPEG data I/O
 * RCS:
 *	$Id: JpegUtil.c
 * DESTINATION:
 *
 * DESCRIPTION:
 *	Routines to read/write jpegs
 * --------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jerror.h>
#include <setjmp.h>
#include "JpegUtil.h"

/******************** JPEG COMPRESSION SAMPLE INTERFACE *******************/

/*
 * We present a minimal version that does not worry about refinements such
 * as error recovery (the JPEG code will just exit() if it gets an error).
 */

/*
 * IMAGE DATA FORMATS:
 *
 * The standard input image format is a rectangular array of pixels, with
 * each pixel having the same number of "component" values (color channels).
 * Each pixel row is an array of JSAMPLEs (which typically are unsigned chars).
 * If you are working with color data, then the color values for each pixel
 * must be adjacent in the row; for example, R,G,B,R,G,B,R,G,B,... for 24-bit
 * RGB color.
 *
 */

/***************************************************************************
 * ERROR HANDLING:
 *
 * The JPEG library's standard error handler (jerror.c) is divided into
 * several "methods" which you can override individually.  This lets you
 * adjust the behavior without duplicating a lot of code, which you might
 * have to update with each future release.
 *
 * Our example here shows how to override the "error_exit" method so that
 * control is returned to the library's caller when a fatal error occurs,
 * rather than calling exit() as the standard error_exit method does.
 *
 * We use C's setjmp/longjmp facility to return control.  This means that the
 * routine which calls the JPEG library must first execute a setjmp() call to
 * establish the return point.  We want the replacement error_exit to do a
 * longjmp().  But we need to make the setjmp buffer accessible to the
 * error_exit routine.  To do this, we make a private extension of the
 * standard JPEG error handler object.  (If we were using C++, we'd say we
 * were making a subclass of the regular error handler.)
 *
 * Here's the extended error handler struct:
 */

struct write_error_mgr {
	struct jpeg_error_mgr pub;	/* "public" fields */
	jmp_buf setjmp_buffer;		/* for return to caller */
};

typedef struct write_error_mgr * write_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

/***************************************************************************/
static void
write_error_exit (j_common_ptr cinfo)
{
  /* Do all the error processing when we return to the setjmp point... */
  /*
   * cinfo->err really points to a my_error_mgr struct,
   * so coerce pointer
   */
  write_error_ptr write_err = (write_error_ptr) cinfo->err;

  /* here is where the message is outputed
   * ...we may want to do something else..*/
  /*
    (*cinfo->err->output_message) (cinfo);
  */

  /* Return control to the setjmp point */
  longjmp(write_err->setjmp_buffer, 1);
}

/***************************************************************************
 * Routine for JPEG compression.
 * Creates a JPEG compressed image of the pixels stored in memory
 *
 * Parameters:
 *	image_buf - pointer to a 2-d array of pixel values
 *	filename  - target file name for JPEG image
 *	quality	  - JPEG compression quality of output image
 *			(ranges from 0 - 100)
 *	image_width  - width of the image in pixels
 *	image_height - height of the image in pixels
 *
 * Returns:
 *	 0 - succesfully wrote a JPEG file
 *	-1 - an error occured and errno
 */

static int
write_JPEG_file (JSAMPARRAY *image_buf, char * filename,
		 int quality, int image_width, int image_height,
		 int input_components, J_COLOR_SPACE in_color_space)
{
  /* This struct contains the JPEG compression parameters and pointers to
   * working space (which is allocated as needed by the JPEG library).
   * It is possible to have several such structures, representing multiple
   * compression/decompression processes, in existence at once.  We refer
   * to any one struct (and its associated working data) as a
   * "JPEG object".
   */
  struct jpeg_compress_struct cinfo;

  /* This struct represents a JPEG error handler.
   * It is declared separately because applications often want to
   * supply a specialized error handler (see the second half of
   * this file for an example).  But here we just take the easy way
   * out and use the standard error handler, which will
   * print a message on stderr and call exit() if compression fails.
   * Note that this struct must live as long as the main JPEG parameter
   * struct, to avoid dangling-pointer problems.
   */
  struct write_error_mgr jerr;

  /* More stuff */
  FILE * outfile;		/* target file */

  /* Step 1: allocate and initialize JPEG compression object */

  /* We have to set up the error handler first, in case the
   * initialization step fails.  (Unlikely, but it could happen
   * if you are out of memory.)
   * This routine fills in the contents of struct jerr, and returns jerr's
   * address which we place into the link field in cinfo.
   */

  /*
    cinfo.err = jpeg_std_error(&jerr);
  */

  /* Here we use the library-supplied code to send compressed data to a
   * stdio stream.  You can also write your own code to do something else.
   * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
   * requires it in order to write binary files.
   */
  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    exit(1);
  }

  /* We set up the normal JPEG error routines, then override error_exit */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = write_error_exit;

  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input
     * file, and return.
     */

    /* Need to do something with the error */
    char buffer[JMSG_LENGTH_MAX];
    (cinfo.err->format_message) ((j_common_ptr)&cinfo, buffer);
    fprintf(stderr, "JPEG Error=>%s\n", buffer);

    /* Clean up the JPEG object */
    jpeg_destroy_compress(&cinfo);
    fclose(outfile);
    return -1;
  }


  /* Now we can initialize the JPEG compression object. */
  jpeg_create_compress(&cinfo);

  /* Step 2: specify data destination (eg, a file) */
  /* Note: steps 2 and 3 can be done in either order. */

  jpeg_stdio_dest(&cinfo, outfile);

  /* Step 3: set parameters for compression */

  /* First we supply a description of the input image.
   * Four fields of the cinfo struct must be filled in:
   */

  /* Set image width and height, in pixels */
  cinfo.image_width = image_width;
  cinfo.image_height = image_height;

  /* Set # of color components per pixel */
  cinfo.input_components = input_components;

  /* colorspace of input image */
  cinfo.in_color_space = in_color_space;


  /* Now use the library's routine to set default compression parameters.
   * (You must set at least cinfo.in_color_space before calling this,
   * since the defaults depend on the source color space.)
   */
  jpeg_set_defaults(&cinfo);

  /* Now you can set any non-default parameters you wish to.
   * Here we just illustrate the use of quality (quantization table) scaling:
   */
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

  /* Step 4: Start compressor */

  /* TRUE ensures that we will write a complete interchange-JPEG file.
   * Pass TRUE unless you are very sure of what you're doing.
   */
  jpeg_start_compress(&cinfo, TRUE);

  /* Step 5: while (scan lines remain to be written) */
  /*           jpeg_write_scanlines(...); */

  /* Here we use the library's state variable cinfo.next_scanline as the
   * loop counter, so that we don't have to keep track ourselves.
   * To keep things simple, we pass one scanline per call; you can pass
   * more if you wish, though.
   */
#if 0
  while (cinfo.next_scanline < cinfo.image_height) {
    /* jpeg_write_scanlines expects an array of pointers
     * to scanlines.
     * Here the array is only one element long, but you could pass
     * more than one scanline at a time if that's more convenient.
     */

    JSAMPLE *rpoint;	/* pointer to JSAMPLE row[s] */
    rpoint = &(*(*image_buf)[cinfo.next_scanline]);

    (void) jpeg_write_scanlines(&cinfo, &rpoint, 1);
  }
#endif
  jpeg_write_scanlines(&cinfo, *image_buf, cinfo.image_height);

  /* Step 6: Finish compression */

  jpeg_finish_compress(&cinfo);
  /* After finish_compress, we can close the output file. */
  fclose(outfile);

  /* Step 7: release JPEG compression object */

  /* This is an important step since it will release a good
   * deal of memory.
   */
  jpeg_destroy_compress(&cinfo);

  /* And we're done! */
  return(1);
}

/***************************************************************************
 * SOME FINE POINTS:
 *
 * In the above loop, we ignored the return value of jpeg_write_scanlines,
 * which is the number of scanlines actually written.  We could get away
 * with this because we were only relying on the value of cinfo.next_scanline,
 * which will be incremented correctly.  If you maintain additional loop
 * variables then you should be careful to increment them properly.
 * Actually, for output to a stdio stream you needn't worry, because
 * then jpeg_write_scanlines will write all the lines passed (or else exit
 * with a fatal error).  Partial writes can only occur if you use a data
 * destination module that can demand suspension of the compressor.
 * (If you don't know what that's for, you don't need it.)
 *
 * If the compressor requires full-image buffers (for entropy-coding
 * optimization or a multi-scan JPEG file), it will create temporary
 * files for anything that doesn't fit within the maximum-memory setting.
 * (Note that temp files are NOT needed if you use the default parameters.)
 * On some systems you may need to set up a signal handler to ensure that
 * temporary files are deleted if the program is interrupted.  See libjpeg.doc.
 *
 * Scanlines MUST be supplied in top-to-bottom order if you want your JPEG
 * files to be compatible with everyone else's.  If you cannot readily read
 * your data in that order, you'll need an intermediate array to hold the
 * image.  See rdtarga.c or rdbmp.c for examples of handling bottom-to-top
 * source data using the JPEG code's internal virtual-array mechanisms.
 */

static FILE	    *imagefile;
static char         *imagebuffer;
static JSAMPARRAY   jpeg_buffer;	/* Pointer to the whole image */
static int	    jpeg_height;	/* Number of rows in image */
static int	    jpeg_width;	/* Number of columns in image */
static int	    jpeg_depth;	/* Number of bytes per element(maximum) */

int jpeg_source(char *fname, int width, int height, int maximum_depth) {
  int i;

  if ((imagefile=fopen(fname, "rb")) == NULL) {
    fprintf(stderr, "cannot open subsample file: %s\n", fname);
    return(-1);
  }
  imagebuffer = NULL;

  jpeg_width = width;
  jpeg_height = height;
  jpeg_depth = maximum_depth;

  /* Allocate memory for a pointer to each row */
  if ((jpeg_buffer = (JSAMPROW *)malloc(sizeof(JSAMPROW)*jpeg_height)) == NULL) {
    fprintf(stderr, "Unable to alloc memory for image\n");
    return(-1);
  }

  for(i = 0; i < jpeg_height; i++) {
    /* Allocate memory for all the pixels in this row */
    if ((jpeg_buffer[i] = (JSAMPLE *)malloc(jpeg_depth*jpeg_width*sizeof(JSAMPLE)))
		== NULL) {
      fprintf(stderr, "Unable to malloc a row\n");
      return(-1);
    }
  }
  return(0);
}

int jpeg_source_mem(char *buffer, int width, int height, int maximum_depth) {
  int i;

  imagefile = NULL;
  imagebuffer = buffer;

  jpeg_width = width;
  jpeg_height = height;
  jpeg_depth = maximum_depth;

  /* Allocate memory for a pointer to each row */
  if ((jpeg_buffer = (JSAMPROW *)malloc(sizeof(JSAMPROW)*jpeg_height)) == NULL) {
    fprintf(stderr, "Unable to alloc memory for image\n");
    return(-1);
  }

  for(i = 0; i < jpeg_height; i++) {
    /* Allocate memory for all the pixels in this row */
    if ((jpeg_buffer[i] = (JSAMPLE *)malloc(jpeg_depth*jpeg_width*sizeof(JSAMPLE)))
		== NULL) {
      fprintf(stderr, "Unable to malloc a row\n");
      return(-1);
    }
  }
  return(0);
}

/****************************/
void jpeg_close()
{
  int i;

  if (imagefile) {
    fclose(imagefile);
    imagefile = NULL;
  }
  if (jpeg_buffer) {
    for(i = 0; i < jpeg_height; i++) {
      if (jpeg_buffer[i]) {
        free(jpeg_buffer[i]);
        jpeg_buffer[i] = NULL;
      }
    }
    free(jpeg_buffer);
    jpeg_buffer = NULL;
  }
}

/****************************/
void jpeg_new()
{
  int		i;

  if (jpeg_buffer)
    for (i = 0; i < jpeg_height; i++)
      if (jpeg_buffer[i])
        memset(jpeg_buffer[i], 0, jpeg_width * jpeg_depth);
}

/****************************/
int jpeg_add(int which, int where, int count)
{
  if (jpeg_buffer && (where >= 0) && (where < jpeg_height)) {
#ifdef DEBUG
    printf("putting offset %d into scene at line %d\n", which, where);
#endif
    fseek(imagefile, which, SEEK_SET);
    if (jpeg_buffer[where])
      return(fread(jpeg_buffer[where], count, 1, imagefile) == 1);
  } else
    printf("error: putting offset %d into scene at line %d\n", which, where);

  return(0);
}

int jpeg_add_mem(int which, int where, int count)
{
  if (jpeg_buffer && (where >= 0) && (where < jpeg_height)) {
    if (jpeg_buffer[where]) {
      memcpy(jpeg_buffer[where], &imagebuffer[which], count);
      return 1;
    }
  } else
    printf("error: putting offset %d into scene at line %d\n", which, where);

  return(0);
}

// MAX_IMAGE_HEIGHT=10240; JCS_RGB=2; JCS_GRAYSCALE=1
/****************************/
int jpeg_write(char *outfilename, int jpeg_depth)
{
    J_COLOR_SPACE	in_color_space = JCS_RGB;	/* colorspace of input image */
    if (jpeg_depth == 1)
        in_color_space = JCS_GRAYSCALE;

    write_JPEG_file (&jpeg_buffer, outfilename, 75, jpeg_width, jpeg_height, jpeg_depth, in_color_space);
    return(0);
}

int createJpegImage(char *outfilename, char *img_buf, int quality, int width, int height, int depth)
{
    JSAMPROW jbuf[MAX_IMAGE_HEIGHT];
    JSAMPARRAY jarray = jbuf;

    J_COLOR_SPACE	in_color_space = JCS_RGB;	/* colorspace of input image */
    if (depth == 1)
        in_color_space = JCS_GRAYSCALE;

    if (height > MAX_IMAGE_HEIGHT)
        height = MAX_IMAGE_HEIGHT;

    int i;
    for (i = 0; i < height; i ++)
        jbuf[i] = (JSAMPROW)(img_buf + i * width * depth);

    write_JPEG_file (&jarray, outfilename, quality, width, height, depth, in_color_space);
    return(0);
}

static int jpeg_compress(void *dest, JSAMPARRAY image_buf,
			 int quality, int image_width, int image_height,
			 int input_components,
			 J_COLOR_SPACE in_color_space)
{
  struct jpeg_compress_struct cinfo;
  struct write_error_mgr jerr;
  size_t outsize;

  /* Step 1: allocate and initialize JPEG compression object */

  /* We set up the normal JPEG error routines, then override error_exit */
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = write_error_exit;

  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(jerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object, close the input
     * file, and return.
     */

    /* Need to do something with the error */
    char buffer[JMSG_LENGTH_MAX];
    (cinfo.err->format_message) ((j_common_ptr)&cinfo, buffer);
    fprintf(stderr, "JPEG Error=>%s\n", buffer);

    /* Clean up the JPEG object */
    jpeg_destroy_compress(&cinfo);
    return -1;
  }

  /* Now we can initialize the JPEG compression object. */
  jpeg_create_compress(&cinfo);

  /* Step 2: specify data destination (eg, a file) */
  /* Note: steps 2 and 3 can be done in either order. */

  jpeg_memory_dest(&cinfo, dest);

  /* Step 3: set parameters for compression */

  /* Set image width and height, in pixels */
  cinfo.image_width = image_width;
  cinfo.image_height = image_height;

  /* Set # of color components per pixel */
  cinfo.input_components = input_components;

  /* colorspace of input image */
  cinfo.in_color_space = in_color_space;

  jpeg_set_defaults(&cinfo);

  /* Now you can set any non-default parameters you wish to.
   * Here we just illustrate the use of quality (quantization table) scaling:
   */
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

  /* Step 4: Start compressor */

  jpeg_start_compress(&cinfo, TRUE);

  /* Step 5: while (scan lines remain to be written) */
  /*           jpeg_write_scanlines(...); */

  jpeg_write_scanlines(&cinfo, image_buf, cinfo.image_height);

  /* Step 6: Finish compression */

  jpeg_finish_compress(&cinfo);

  outsize = jpeg_memory_size(&cinfo);

  /* Step 7: release JPEG compression object */

  /* This is an important step since it will release a good
   * deal of memory.
   */
  jpeg_destroy_compress(&cinfo);

  /* And we're done! */
  return(outsize);
}

int image_to_jpeg_memory(char *dest, char *img_buf, int quality,
			 int width, int height, int depth)
{
  JSAMPROW jbuf[MAX_IMAGE_HEIGHT];
  JSAMPARRAY jarray = jbuf;
  int i;

  J_COLOR_SPACE	in_color_space = JCS_RGB;	/* colorspace of input image */
  if (depth == 1)
    in_color_space = JCS_GRAYSCALE;

  if (height > MAX_IMAGE_HEIGHT) height =MAX_IMAGE_HEIGHT;
  for(i = 0; i < height; i++) {
    jbuf[i] = (JSAMPROW)(img_buf + i * width * depth);
  }
  return jpeg_compress(dest, jarray, quality, width, height, depth,
		       in_color_space);
}


