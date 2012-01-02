/*--------------------------------------------------------------------
 *	$Id:$
 *
 *	Copyright (c) 2004-2012 by J. Luis
 *
 * 	This program is part of Mirone and is free software; you can redistribute
 * 	it and/or modify it under the terms of the GNU Lesser General Public
 * 	License as published by the Free Software Foundation; either
 * 	version 2.1 of the License, or any later version.
 * 
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * 	Lesser General Public License for more details.
 *
 *	Contact info: w3.ualg.pt/~jluis/mirone
 *--------------------------------------------------------------------*/

/* Program:	cvcolor_mex.c
 * Purpose:	matlab callable routine to do color space conversions
 * 		using the OpenCV library
 *
 *		Usage: B = cvcolor_mex(A,'TRF');
 *		       where A is a uint8 MxNx3 rgb image OR a Mx3 colormap (double):
 *		       TRF is a string controling the transformation. Possibilities are:
 *		       rgb2lab,lab2rgb, rgb2luv,luv2rgb, rgb2xyz,xyz2rgb
 *		       rgb2yiq,yiq2rgb, rgb2hsv,luv2hsv, rgb2hsl,hsl2rgb, rgb2YCrCb,YCrCb2rgb
 *
 *	NOTE1:	The transformation of ColorMaps is done by transforming the [0 1] interval
 *		into [0 255]. This is obviously much less precise than Matlab direct computations
 *
 *	NOTE2:	I left in the code the option to do a GRAY2RGB - which does not work - and because
 *		of it the code is a bit confusing in some parts.
 *
 *	NOTE3:	Only UInt8 are allowed because the OpenCV cvCvtColor function give either errors
 *		or non-sense results. Oddly, there is no documentation for this function in the docs.
 *
 *	NOTE4:	This code is amazingly fast. For the rgb2lab tansformation it is more than 60 times
 *		faster than Matlab and with only one copy of input as a memory overload.
 *
 * Revision 1.0  12/6/2006 Joaquim Luis
 *
 */

#include <math.h>
#include "mex.h"
#include <cv.h>

void interleave(unsigned char in[], unsigned char out[], int nx, int ny, int nBands, int dir);
void cmap2colorvec(double *in8_img, unsigned char *in_img, double *out8_img, unsigned char *out_img, int ny, int dir);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	unsigned char *out_img, *in_img, *tmp_img;
	int nx, ny, n_bands_in, n_bands_out, out_dims[3], n_xy, m, n, k;
	int n_dims, cv_code, mx_CLASS, nBytes;
	int src_step = 0, dst_step;
	const int *dim_array;
	char	*argv;
	double *in8_img, *out8_img;
	IplImage* src_img = 0;
	IplImage* dst_img = 0;
	mxArray *ptr, *ptr_u8, *ptr_d;

	/* ---- Check for errors in user's call to function.  -------------------- */
	if (nrhs == 0) {
		mexPrintf("Usage: B = cvcolor_mex(A,'TRF');\n");
		mexPrintf("       where A is a uint8 MxNx3 rgb image OR a Mx3 colormap (double):\n");
		mexPrintf("       TRF is a string controling the transformation. Possibilities are:\n");
		mexPrintf("       rgb2lab,lab2rgb, rgb2luv,luv2rgb, rgb2xyz,xyz2rgb\n");
		mexPrintf("       rgb2yiq,yiq2rgb, rgb2hsv,luv2hsv, rgb2hsl,hsl2rgb, rgb2YCrCb,YCrCb2rgb\n");
		return;
	}
	if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]))
		mexErrMsgTxt("CVCOLOR_MEX: Input image must be a real matrix");
	else if (nlhs > 2)
		mexErrMsgTxt("CVCOLOR_MEX returns only one output argument!");


	if (nrhs != 2)
		mexErrMsgTxt("CVCOLOR_MEX requires 2 input arguments!");

	if(!mxIsChar(prhs[1]))
		mexErrMsgTxt("RGB2XYZ Second argument must be a valid string!");
	else {
		argv = (char *)mxArrayToString(prhs[1]);
		//<X>/<Y> = RGB, BGR, GRAY, HSV, YCrCb, XYZ, Lab, Luv, HLS
		if (!strcmp(argv,"rgb2lab"))
			cv_code = CV_RGB2Lab;
		else if (!strcmp(argv,"lab2rgb"))
			cv_code = CV_Lab2RGB;

		else if (!strcmp(argv,"rgb2luv"))
			cv_code = CV_RGB2Luv;
		else if (!strcmp(argv,"luv2rgb"))
			cv_code = CV_Luv2RGB;

		else if (!strcmp(argv,"rgb2xyz"))
			cv_code = CV_RGB2XYZ;
		else if (!strcmp(argv,"xyz2rgb"))
			cv_code = CV_XYZ2RGB;

		else if (!strcmp(argv,"rgb2yiq") || !strcmp(argv,"rgb2gray"))
			cv_code = CV_RGB2GRAY;
		else if (!strcmp(argv,"yiq2rgb") || !strcmp(argv,"gray2rgb"))
			cv_code = CV_GRAY2RGB;

		else if (!strcmp(argv,"rgb2hsv"))
			cv_code = CV_RGB2HSV;
		else if (!strcmp(argv,"hsv2rgb"))
			cv_code = CV_HSV2RGB;

		else if (!strcmp(argv,"rgb2hls"))
			cv_code = CV_RGB2HLS;
		else if (!strcmp(argv,"hls2rgb"))
			cv_code = CV_HLS2RGB;

		else if (!strcmp(argv,"rgb2YCrCb"))
			cv_code = CV_RGB2YCrCb;
		else if (!strcmp(argv,"YCrCb2rgb"))
			cv_code = CV_YCrCb2RGB;

		else {
			mexPrintf("CVCOLOR_MEX ERROR: Unknown conversion type string.\n");
			mexPrintf("Valid types: rgb2lab,lab2rgb, rgb2luv,luv2rgb, rgb2xyz,xyz2rgb\n");
			mexPrintf("             rgb2yiq,yiq2rgb, rgb2hsv,luv2hsv, rgb2gray,gray2rgb\n");
			mexErrMsgTxt("             rgb2hsl,hsl2rgb, rgb2YCrCb,YCrCb2rgb.\n");
		}
	}

	/* Find out in which data type was given the input array */
	if (mxIsUint8(prhs[0])) {
		mx_CLASS  = mxUINT8_CLASS;
		nBytes    = 1;
	}
	else if (mxIsDouble(prhs[0])) {
		mx_CLASS  = mxDOUBLE_CLASS;
		nBytes    = 8;
	}
	else {
		mexPrintf("CVCOLOR_MEX ERROR: Invalid input data type.\n");
		mexErrMsgTxt("Valid types are: UInt8 or double.\n");
	}

	n_dims = mxGetNumberOfDimensions(prhs[0]);
	dim_array=mxGetDimensions(prhs[0]);
	ny = dim_array[0];
	nx = dim_array[1];
	n_bands_in = dim_array[2];

	if (nBytes == 8 && n_dims != 2)	/* Doubles only in ML colormap */
		mexErrMsgTxt("CVCOLOR_MEX ERROR: input double only for Mx3 colormaps.\n");

	if (nBytes == 1 && n_dims == 2)
		mexErrMsgTxt("CVCOLOR_MEX ERROR: input must be MxNx3 UInt8 array.\n");

	if (cv_code == CV_GRAY2RGB && n_dims != 2)	/* It's not working anyway */
		mexErrMsgTxt("CVCOLOR_MEX ERROR: 'gray2rgb' option implies that input data is MxN.\n");
	/* -------------------- End of input error tests ------------------------------------- */

	if (n_dims == 2 && nBytes == 1)	/* Otherwise it would stay undefined */
		n_bands_in = 1;
	else if (nBytes == 8)		/* ML cmap, so lets pretend it is a 3D (changed below into one) */
		n_bands_in = 3;

	/* ------ GET POINTER TO INPUT IMAGE DATA --------- */ 
	if (nBytes == 8) {
		in8_img = mxGetPr(prhs[0]);
		in_img = (unsigned char *)mxCalloc (3*ny, sizeof (unsigned char));
		cmap2colorvec(in8_img,in_img,out8_img,out_img,ny,1);/* Convert the ML colormap into a RGB vector */
	}
	else
		in_img = mxGetData(prhs[0]);
	/* ------------------------------------------------ */ 

	out_dims[0] = ny;
	out_dims[1] = nx;
	if (nBytes == 8) {	/* ML cmap. We need to fake that the output is 1xMx3 */
		n_bands_out = 3;
		plhs[0] = mxCreateNumericMatrix(ny, 3, mx_CLASS, mxREAL);
		out_dims[2] = 3;	/* May seam silly but it is used for intermediary array */
	}
	else if (cv_code == CV_RGB2GRAY) {	/* rgb2gray */
		n_bands_out = 1;
		out_dims[2] = 3;	/* May seam silly but it is used for intermediary array */
		plhs[0] = mxCreateNumericMatrix(ny, nx, mx_CLASS, mxREAL);
	}
	else if (cv_code == CV_GRAY2RGB) {	/* Here we must force a 3D output array */
		n_bands_out = 3;
		out_dims[2] = n_bands_out;
		plhs[0] = mxCreateNumericArray(3, out_dims, mx_CLASS, mxREAL);
	}
	else {
		n_bands_out = n_bands_in;
		out_dims[2] = n_bands_out;
		plhs[0] = mxCreateNumericArray(3, out_dims, mx_CLASS, mxREAL);
	}

	/* ------ GET POINTER TO OUTPUT IMAGE DATA --------- */ 
	if (nBytes == 1)
		out_img = mxGetData(plhs[0]);
	else {
		out_img = (unsigned char *)mxCalloc (3*ny, sizeof (unsigned char));
		out8_img = mxGetPr(plhs[0]);	/* Remember that output is double */
	}
	/* ------------------------------------------------- */ 

	if (cv_code != CV_GRAY2RGB) 	/* tmp_img will be 3D */
		ptr = mxCreateNumericArray(3, out_dims, mxUINT8_CLASS, mxREAL);
	else				/* tmp_img will be 2D */ 
		ptr = mxCreateNumericMatrix(ny, nx, mxUINT8_CLASS, mxREAL);
	tmp_img = mxGetData(ptr);

	if (nBytes == 8) {	/* From here on nx & ny must reflect the fact its a row vector... */
		nx = ny;	ny = 1;
		if (cv_code == CV_RGB2GRAY) 	/* Case of a colormap to gray */
			n_bands_out = 1; 
	}

	if (cv_code == CV_GRAY2RGB) 	/* We need to copy input to tmp array */
		for (m = 0; m < nx*ny; m++) tmp_img[m] = in_img[m];
	else				/* We need to change interleaving from ML order to BIP */
		interleave (in_img, tmp_img, nx, ny, n_bands_in, 1);
    
	src_step = nx * n_bands_in;
	dst_step = nx * n_bands_out;
    
	src_img = cvCreateImageHeader( cvSize(nx, ny), IPL_DEPTH_8U, n_bands_in );
	cvSetImageData( src_img, (void *)tmp_img, src_step );

	dst_img = cvCreateImageHeader(cvGetSize(src_img), IPL_DEPTH_8U, n_bands_out );
	cvSetImageData( dst_img, (void *)out_img, dst_step );

	cvCvtColor(src_img,dst_img,cv_code); // src -> dst


	if (n_bands_out == 3) {		/* NOT rgb2gray or ML cmap */
		if (cv_code == CV_GRAY2RGB) {
			mxDestroyArray(ptr);		/* Destroy it because it was 2D and we need a 3D */
			ptr = mxCreateNumericArray(3, out_dims, mx_CLASS, mxREAL);
			tmp_img = mxGetData(ptr);
		}
		/* Copy the pixel interleaved result into the temporary array */
		for (m = 0; m < nx*ny*3; m++) tmp_img[m] = out_img[m];
		interleave (tmp_img, out_img, nx, ny, 3, -1); /* Change from BIP to ML order */

	}
	if (nBytes == 8) { 	/* Convert the RGB row vector into a Mx3 ML colormap */
		if (cv_code == CV_RGB2GRAY)
			cmap2colorvec(in8_img,in_img,out8_img,out_img,nx,-1);	/* nx because we swapped it above */
		else
			cmap2colorvec(in8_img,in_img,out8_img,out_img,nx,0);	/* nx because we swapped it above */

		mxFree ((void *)in_img);
		mxFree ((void *)out_img);
	}

	cvReleaseImageHeader( &src_img );
	cvReleaseImageHeader( &dst_img );
	mxDestroyArray(ptr);
}

/* -------------------------------------------------------------------------------------------- */
void interleave(unsigned char in[], unsigned char out[], int nx, int ny, int nBands, int dir) {
	int n_xy, m, n, k, c = 0;

	n_xy = nx*ny;

	if (dir == 1) {		/* Matlab order to r0,g0,b0, r1,g1,b1, etc ... */
		for (n = 0; n < nx; n++) {
			for (m = 0; m < ny; m++) {
				for (k = 0; k < nBands; k++)
					out[c++] = in[m+n*ny+k*n_xy];
			}
		}
	}
	else {			/* r0,g0,b0, r1,g1,b1, etc ...  to Matlab order */
		for (k = 0; k < nBands; k++) {		/* nBands MUST BE == 3 */
			for (n = 0; n < n_xy; n++)
				out[c++] = in[n*nBands + k];
		}
	}
}

/* -------------------------------------------------------------------------------------------- */
void cmap2colorvec(double *in8_img, unsigned char *in_img, double *out8_img, unsigned char *out_img, int ny, int dir) {
	/* dir = 1 -> Convert the Mx3 ML colormap with range [0,1] into a pixel interleaved RGB row vector
	   dir = 0 -> convert the row color vector into a Mx3 double ML cmap
	   dir = -1 -> convert the row color vector into a gray Mx3 double ML cmap */
	int m, n, k, c = 0;

	if (dir == 1) {
		for (m = 0; m < 3*ny; m++)
			in_img[m] = (unsigned char)(in8_img[m]*255);
	}
	else if (dir == 0) {
		for (m = 0; m < 3*ny; m++)
			out8_img[m] = (double)(out_img[m]/255.);
	}
	else {		/* 2GRAY, so repeat out_img 3 times */
		for (k = 0; k < 3; k++) {
			for (m = 0; m < ny; m++)
				out8_img[c++] = (double)(out_img[m]/255.);
		}
	}
}
