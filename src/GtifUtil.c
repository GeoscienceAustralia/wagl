/* *************************************************************
 * GtifUtil.c: A set of functions for converting integers of
 * different bits, either signed or unsigned
 *
 * Note: in GeoTIFF, all 8-bit integers are unsigned. That is,
 *       the scanline should be the data type of unsigned char.
 *       The convert of unsigned 8-bit integers into signed 16-bit
 *       integers is using the overflow of integer.
 *       For example, two unsigned 8-bit integers are: a=0, b=128.
 *       Using formula: b*256+a=128*256+0=32768. If it is for
 *       unsigned 16-bit, the converted value is 32768. But if it
 *       is for signed 16-bit, this value is overflow because the
 *       maximum signed 16-bit integer is 32767. So the result
 *       value should be: 32768-65536=32768. 65536 is the maimum
 *       signed 16-bit integer value.
 *
 * @author: Frank Q. Fu  17 April ,2006
 * Copyright (c) 2006-? ACRES GEMD, Geoscience Australia
 *
 * **************************************************************
 */

#include <stdio.h>
#include <stdlib.h>

#include "GtifUtil.h"

// ------------- signed 16-bit vs unsigned 8-bit ---------------
// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order
int16 uint8ToInt16 (int byteOrder, uint8 a, uint8 b)
{
	if (byteOrder == BIG_ENDIAN)
		return a*256 + b;
	else
		return b*256 + a;
}

// convert 2 unsigned 8-bit integers into a signed 16-bit integer
// using given byte order. In GeoTIFF, all bytes (char) are unsigned
// 8-bit integers
void int16ToUint8 (int byteOrder, uint8 *bytes, int16 data)
{
	// first cast signed 16-bit integer into unsigned 16-bit integer
	uint16 data_uint16 = (uint16)data;

	if (byteOrder == BIG_ENDIAN)
	{
		bytes[0] = data_uint16 / 256;
		bytes[1] = data_uint16 % 256;
	}
	else // litlte endian
	{
		bytes[1] = data_uint16 / 256;
		bytes[0] = data_uint16 % 256;
	}
}

// ------------- signed 16 bits vs signed 8 bits ---------------
// convert 2 signed 8-bit integers into a signed 16-bit integer
// using given byte order
// It seems geotif uses little-endian. E.g., the NULL pixel
// is -32767 and corresponding 2 bytes are 1 and -128
// Using the following formular: -128*256+1=-32768+1=-32767
int16 int8ToInt16 (int byteOrder, int8 a, int8 b)
{
	if (byteOrder == BIG_ENDIAN)
		return a*256 + b;
	else
		return b*256 + a;
}

// convert 2 bytes into a int16 integer assuming big-endian
void int16ToInt8 (int byteOrder, int8 *bytes, int16 data)
{
	int16 tmp = data % 256; // -255 ~ 255
	int8 b = tmp;
	int8 a = data/256;

	// singed 8-bit: -128 ~ 127
	if (tmp > 127) // 128~255
	{
		a += 1;
		b = tmp - 256;
	}
	else if (tmp < -128) // -255~-129
	{
		a -= 1;
		b = tmp + 256;
	}

	if (byteOrder == BIG_ENDIAN)
	{
		bytes[0] = a;
		bytes[1] = b;
	}
	else
	{
		bytes[0] = b;
		bytes[1] = a;
	}
}

// ------------- unsigned 16 bits vs unsigned 8 bits -------------
// convert an unsigned 16-bit integer into two unsigned 8-bit
// integers using given byte order
void uint16ToUint8 (int byteOrder, uint8 *bytes, uint16 data)
{
	if (byteOrder == BIG_ENDIAN)
	{
		bytes[0] = data/256;    // 2^8
		bytes[1] = data % 256;
	}
	else // if (byteOrder == LITTLE_ENDIAN)
	{
		bytes[1] = data/256;      // 2^8
		bytes[0] = data % 256;
	}
}

// convert two unsigned 8-bit integers into an unsigned 16-bit integer in big endian
uint16 uint8ToUint16 (int byteOrder, uint8 a, uint8 b)
{
	if (byteOrder == BIG_ENDIAN)
		return a*256 + b;
	else
		return b*256 + a;
}

// --------------- signed 32-bit vs unsigned 8-bit -------------

// convert four unsigned 8-bit integers into a signed 32-bit
// integer using given byte order
int32 uint8ToInt32 (int byteOrder, uint8 *bytes)
{
	// 256*256*256=16777216 and 256*256=65536
	if (byteOrder == BIG_ENDIAN)
		return bytes[0]*16777216+bytes[1]*65536+bytes[2]*256+bytes[3];
	else
		return bytes[3]*16777216+bytes[2]*65536+bytes[1]*256+bytes[0];
}

// convert a signed 32-bit integers into four unsigned 8-bit
// integers using given byte order
void int32ToUint8 (int byteOrder, uint8 *bytes, int32 data_int32)
{
	// cast signed 32-bit integer into unsigned 32-bit integer
	uint32 data = (uint32)data_int32;

	uint8 a = data/16777216; //256*256*256=16777216
	data = data%16777216;

	uint8 b = data/65536;    //256*256=65536=2^16
	data = data%65536;

	uint8 c = data/256;      // 2^8
	uint8 d = data % 256;

	if (byteOrder == BIG_ENDIAN)
	{
		bytes[0] = a;
		bytes[1] = b;
		bytes[2] = c;
		bytes[3] = d;
	}
	else
	{
		bytes[0] = d;
		bytes[1] = c;
		bytes[2] = b;
		bytes[3] = a;
	}
}
// --------------- unsigned 32-bit vs unsigned 8-bit -------------

// convert four unsigned 8-bit integers into a signed 32-bit
// integer using given byte order
uint32 uint8ToUint32 (int byteOrder, uint8 *bytes)
{
	// 256*256*256=16777216 and 256*256=65536
	if (byteOrder == BIG_ENDIAN)
		return bytes[0]*16777216+bytes[1]*65536+bytes[2]*256+bytes[3];
	else
		return bytes[3]*16777216+bytes[2]*65536+bytes[1]*256+bytes[0];
}

// convert a signed 32-bit integers into four unsigned 8-bit
// integers using given byte order
void uint32ToUint8 (int byteOrder, uint8 *bytes, uint32 data)
{
	uint32 a = data/16777216; //256*256*256=16777216
	data = data%16777216;

	uint32 b = data/65536;    //256*256=65536=2^16
	data = data%65536;

	uint32 c = data/256;      // 2^8
	uint32 d = data % 256;

	if (byteOrder == BIG_ENDIAN)
	{
		bytes[0] = a;
		bytes[1] = b;
		bytes[2] = c;
		bytes[3] = d;
	}
	else
	{
		bytes[0] = d;
		bytes[1] = c;
		bytes[2] = b;
		bytes[3] = a;
	}
}

// --------------- signed 32-bit vs signed 8-bit -------------

// convert four 8-bit signed integers into a 32-bit signed
// integer using given byte order
int32 int8ToInt32 (int byteOrder, int8 *bytes)
{
	// 256*256*256=16777216 and 256*256=65536
	if (byteOrder == BIG_ENDIAN)
		return bytes[0]*16777216+bytes[1]*65536+bytes[2]*256+bytes[3];
	else
		return bytes[3]*16777216+bytes[2]*65536+bytes[1]*256+bytes[0];
}

// convert a 32-bit signed integers into four 8-bit unsigned
// integers using given byte order
void int32ToInt8 (int byteOrder, int8 *bytes, int32 data)
{
	int32 a = data/16777216; //256*256*256=16777216=2^24
	data = data%16777216;

	int32 b = data/65536;    //256*256=65536=2^16
	data = data%65536;

	int32 c = data/256;      // 2^8
	int32 d = data%256;

	if (d >= 128 || d < -128)   // 2^7
	{
		c += d/128;
		int sign = 0;
		if (d != 0)
			sign = d/abs(d); // =1, 0 or -1
		d -= sign*256;
	}

	if (c >= 256 || c <= -256)
	{
		b += c/256;
		c -= (c/abs(c))*256;
	}

	if (c >= 128 || c < -128)  // 2^15
	{
		b += c/128;
		int sign = 0;
		if (c != 0)
			sign = c/abs(c); // =1, 0 or -1
		c -= sign*256;
	}

	if (b >= 256 || b <= -256)
	{
		a += b/256;
		b -= (b/abs(b))*256;
	}

	if (b >= 128 || b < -128)  // 2^23
	{
		a += b/128;
		int sign = 0;
		if (b != 0)
			sign = b/abs(b); // =1, 0 or -1
		b -= sign*256;
	}

	if (byteOrder == BIG_ENDIAN)
	{
		bytes[0] = a;
		bytes[1] = b;
		bytes[2] = c;
		bytes[3] = d;
	}
	else
	{
		bytes[0] = d;
		bytes[1] = c;
		bytes[2] = b;
		bytes[3] = a;
	}
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

// Return true if the filename contains ".tif"
boolean isGeoTIFF (char *filename)
{
    if (filename == (char *)NULL)
        return FALSE;
    char *dot_chr = strrchr (filename, '.');
    if (dot_chr == (char *)NULL)
        return FALSE;
    else if (strcmp(dot_chr, ".tif") == 0)
        return TRUE;
    else
        return FALSE;
}

