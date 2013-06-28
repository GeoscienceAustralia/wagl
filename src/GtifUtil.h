/* *************************************************************
 * IntConvert.h: definitions for the convertion of integers of different bits
 *
 *
 * @author: Frank Q. Fu  10 April ,2006
 * Copyright (c) 2006 ACRES GEMD, Geoscience Australia
 *
 * **************************************************************
 */

#include <stdio.h>
#include <stdlib.h>

#include "GtifMosaic.h"

// ------------- signed 16-bit vs unsigned 8-bit ---------------
// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order
int16 uint8ToInt16 (int byteOrder, uint8 a, uint8 b);

// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order. In GeoTIFF, all bytes (char) are unsigned
// 8-bit integers
void int16ToUint8 (int byteOrder, uint8 *bytes, int16 data);

// ------------- signed 16 bits vs signed 8 bits ---------------
// convert 2 signed 8-bit integers into a signed 16-bit integer
// using given byte order
int16 int8ToInt16 (int byteOrder, int8 a, int8 b);

// convert 2 bytes into a int16 integer assuming big-endian
void int16ToInt8 (int byteOrder, int8 *bytes, int16 data);

// ------------- unsigned 16 bits vs unsigned 8 bits -------------
// convert an unsigned 16-bit integer into two unsigned 8-bit integers using given byte order
void uint16ToUint8 (int byteOrder, uint8 *bytes, uint16 data);

// convert two unsigned 8-bit integers into an unsigned 16-bit integer in big endian
uint16 uint8ToUint16 (int byteOrder, uint8 a, uint8 b);

// --------------- signed 32-bit vs unsigned 8-bit -------------

// convert four unsigned 8-bit integers into a signed 32-bit
// integer using given byte order
int32 uint8ToInt32 (int byteOrder, uint8 *bytes);

// convert a signed 32-bit integers into four unsigned 8-bit
// integers using given byte order
void int32ToUint8 (int byteOrder, uint8 *bytes, int32 data);

// --------------- unsigned 32-bit vs unsigned 8-bit -------------

// convert four unsigned 8-bit integers into a signed 32-bit
// integer using given byte order
uint32 uint8ToUint32 (int byteOrder, uint8 *bytes);

// convert a signed 32-bit integers into four unsigned 8-bit
// integers using given byte order
void uint32ToUint8 (int byteOrder, uint8 *bytes, uint32 data);

// --------------- signed 32-bit vs signed 8-bit -------------

// convert four 8-bit signed integers into a 32-bit signed
// integer using given byte order
int32 int8ToInt32 (int byteOrder, int8 *bytes);

// convert a 32-bit signed integers into four 8-bit unsigned
// integers using given byte order
void int32ToInt8 (int byteOrder, int8 *bytes, int32 data);

// remove spare white spaces in the beginning or back of the string
void trim (char *string);

// Return true if the filename contains ".tif"
boolean isGeoTIFF (char *filename);

