#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

/**********************************************************************************************************************\
File: ferroic_core.c

Implementation of the homogenized free energy hysteresis model, with and without thermal relaxation.  Only the
core of the model is is implemented here, with all setup assumed to come from outside.  The paper
"Efficient Implementation of Homogenized Free Energy Model with Emphasis on Thermal Relaxation" by Braun and Smith
(presented at 2005 SPIE smart materials conference) should be consulted for reference, as the faster models
discussed in that paper are presented here.  

This file is implemented as a MATLAB mex function to do the core processing.  All setup is in the matlab file
ferroic_hyst.m.  This c file has been tested to run on Windows, Solaris, and Linux platforms.  Other platforms
should also be able to implement it without difficulty -- it is all ansi.

LICENSE: This work is distributed under BSD license, as copied below
    Copyright (c) 2005, Tom Braun / North Carolina State University / Department of Mathematics and Center for 
        Research in Scientific Computation
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
    following conditions are met:
        * Redistributionsc of source code must retain the above copyright notice, this list of conditions and the 
          following disclaimer.
        * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the 
          following disclaimer in the documentation and/or other materials provided with the distribution.
        * Neither the name of the Tom Braun / North Carolina State University nor the names of its contributors may 
          be used to endorse or promote products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
    INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
    WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
    USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

NOTE: This code is derived from a MATLAB homogenized free energy model written by Andrew Hatch, North Carolina State 
University, Department of Mathematics.  It is used with permission, and much thanks.  

Written by:
Tom Braun
North Carolina State University
Department of Mathematics/Center for Research in Scientific Computation
Raleigh NC  27695
tbraun@pobox.com
January 2005
\**********************************************************************************************************************/

/* NOTE: Pr is assumed to be 1 throughout this code, and is not mentioned or otherwise included */

typedef struct {
  double *points, *weights;
  unsigned int num;
} dist_type;

typedef struct {
  double step, *weights;
  unsigned int num;
} fixeddist_type;

typedef struct {
  double *in, *out;
  unsigned int num;
} field_type;

typedef struct {
  double *p_negpos, *p_posneg, *P_neg, *P_pos;
  unsigned int num;
} likes_type;

extern double erfc(double x);

/**********************************************************************************************************************\
Function: Neg_Relax_Hyst

This represents the negligable relaxation model. 

Input Paramters:
    field.in -- input electric or magnetic field points
    field.num -- number of values in field.in array.  This will also be the number of values placed into field.out
    coercive.points -- quadrature points for the coercive field integration.
    coercive.weights -- quadrature and distribution weights for the coercive field distribution.
    coercive.num -- number of points in the coercive.points and coercive.weights arrays
    effective.points -- quadrature points for the effective field integration
    effective.weights -- quadrature and distribution weights for the effective field integration
    effective.num -- number of points in the effective.points and effecctive.weights arrays
    weightsum -- (sum_i sum_j coercive.weights[i] * effective.weights[j]) / eta
    addit -- (sum_i sum_j coercive.weights[i] * effective.weights[j] * effective.points[j]) / eta

In/Out Parameter:
    xpos -- array of coercive.num*effective.num [indexed coercive_index * effective.num + effective_index].  This is
        assumed to contain the state of the dipoles corresponding to each pair of quadrature points.  It is returned
        back to MATLAB, in case the state needs to be stored for future operations.

Output Parameter:
    field.out -- output polarizations/magnetizations
\**********************************************************************************************************************/

void Neg_Relax_Hyst(field_type field, dist_type coercive, dist_type effective, double weightsum, double addit, double *xpos) {
  unsigned int i, j, k;
  double out;

  for (k = 0; k < field.num; ++k) {
    field.out[k] = addit + weightsum * field.in[k];
    for (i = 0; i < coercive.num; ++i) {
      out = 0;
      for (j = 0; j < effective.num; ++j) {
	if (effective.points[j] + field.in[k] + coercive.points[i]*xpos[i*effective.num + j] > 0) {
	  xpos[i*effective.num + j] = 1;
	  out += effective.weights[j];
	} else {
	  xpos[i*effective.num + j] = -1;
	  out -= effective.weights[j];
	}
      }
      field.out[k] += out * coercive.weights[i];
    }
  }
} /* end void Neg_Relax_Hyst */

/**********************************************************************************************************************\
Function: Relax_Hyst

A speed optimized version of the model including thermal relaxation.  The effective distribution is assumed to be 
use a Newton-Coates formula for integration (fixed, equally spaced quadrature points).

Input Paramters:
    field.in -- input electric or magnetic field points
    field.num -- number of values in field.in array.  This will also be the number of values placed into field.out
    coercive.points -- quadrature points for the coercive field integration.
    coercive.weights -- quadrature and distribution weights for the coercive field distribution.
    coercive.num -- number of points in the coercive.points and coercive.weights arrays
    effective.step -- stepsize for the effective field integration
    effective.weights -- quadrature and distribution weights for the effective field integration
    effective.num -- number of points in the effecctive.weights array
    likes.p_negpos -- values for the likelihood of switching from negative to positive for all emu in range
    likes.p_posneg -- values for the likelihood of switching from positive to negative for all emu in range
    likes.P_neg -- a component of the average polarization of negatively oriented dipoles for all emu in range
    likes.P_pos -- a component of the average polarization of positively oriented dipoles for all emu in range
    likes.num -- number of points in the likes.p_negpos, likes.p_posneg, likes.P_neg, and likes.P_pos arrays
    weightsum -- (sum_i sum_j coercive.weights[i] * effective.weights[j]) / eta
    addit -- (sum_i sum_j coercive.weights[i] * effective.weights[j] * effective.points[j]) / eta
    deltat -- stepsize between successive input field values
    increase -- the number of points added to each side of the original Ee_pts spread when finding emu
    resolution -- an integer divisor of effective.step, used to improve accuracy when rounding to the nearest emu

In/Out Parameter:
    xpos -- array of coercive.num*effective.num [indexed coercive_index * effective.num + effective_index].  This is
        assumed to contain the state of the dipoles corresponding to each pair of quadrature points.  It is returned
        back to MATLAB, in case the state needs to be stored for future operations.

Output Parameter:
    field.out -- output polarizations/magnetizations
\**********************************************************************************************************************/

void Relax_Hyst(field_type field, dist_type coercive, fixeddist_type effective, likes_type likes, double weightsum, 
		double addit, double deltat, unsigned int increase, unsigned int resolution, double *xpos) {
  unsigned int i, j, k, num_res;
  int j_neg, j_pos;
  double k_shift, i_shift, out;

  num_res = effective.num * resolution;

  for (k = 0; k < field.num; ++k) {
    field.out[k] = addit + weightsum*field.in[k];
    k_shift = field.in[k] / effective.step;

    for (i = 0; i < coercive.num; ++i) {
      i_shift = coercive.points[i] / effective.step;

      j_neg = increase + (int)(k_shift - i_shift + 0.5);
      if (j_neg < 0)
	j_neg = 0;
      if (j_neg + num_res >= likes.num)
	j_neg = likes.num - num_res - 1;
      j_pos = increase + (int)(k_shift + i_shift + 0.5);
      if (j_pos < 0)
	j_pos = 0;
      if (j_pos + num_res >= likes.num)
	j_pos = likes.num - num_res - 1;
      
      out = 0;
      for (j = 0; j < effective.num; ++j, j_neg += resolution, j_pos += resolution) {
	xpos[i*effective.num + j] = (xpos[i*effective.num + j] + likes.p_negpos[j_neg] * deltat) / 
          (1 + deltat * (likes.p_negpos[j_neg] + likes.p_posneg[j_pos]));
	out += effective.weights[j] * (xpos[i*effective.num + j] * (likes.P_pos[j_pos] - likes.P_neg[j_neg] + 2) + likes.P_neg[j_neg]);
      }
      field.out[k] += out * coercive.weights[i];
    }
  }
} /* end void Relax_Hyst */

/**********************************************************************************************************************\
Function: mexFunction

Standard MATLAB interface function.  See the MATLAB documentation for more information.  

Inputs:
The mex function should be called with either 8 (for negligable relaxation) or 15 (with relaxation) inputs, in the 
following order.  Note, the 1st 8 are always given.
    field.in -- the input electric or magnetic field
    coercive.points -- coercive quadrature integration points
    coercive.weights -- coercive quadrature integeration and distribution weights
    effective.points (negligable relaxation) or effective.step (with relaxation) -- quadrature points or stepsize
        for effective field integration
    effective.weights -- effective quadrature integration and distribution weights
    weightsum -- see functions above for value
    addit -- see functions above for value
    xpos -- state of the material at each point in the mesh formed by coercive.points and effective.points
The parameters should only be given if the model with relaxation is desired.  In fact, that's how the code tells
which model you want -- whether or not these parameters are present.
    likes.p_negpos -- likelihood of switching from negative to positive
    likes.p_posneg -- likelihood of switching from positive to negative
    likes.P_neg -- a component of the average polarizations for negatively oriented dipoles
    likes.P_pos -- a component of the average polarizations for positively oriented dipoles
    deltat -- stepsize between successive field inputs
    increase -- number of points added to each side of the range of the effective points when forming emu
    resolution -- an integer divisor used to increase accuracy of calcualtions.

Outputs:
    field.out -- resulting polarization or magnetization
    xpos -- the final state of the system
\**********************************************************************************************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  field_type field;
  dist_type effective, coercive;
  fixeddist_type f_effective;
  likes_type likes;
  unsigned int i, increase, resolution;
  double weightsum, addit, *xpos, *tmp, deltat;

  /* Check number of inputs/outputs */
  if (nlhs > 2)
    mexErrMsgTxt("There are only two outputs");
  if ((nrhs != 8) && (nrhs != 15))
    mexErrMsgTxt("Wrong number of input arguments -- must be 8 or 15");
  
  /* Check for errors and get the common parameters */
  if ((! mxIsDouble(prhs[0])) || mxIsComplex(prhs[0]) || ((mxGetN(prhs[0]) != 1) && (mxGetM(prhs[0]) != 1)))
    mexErrMsgTxt("Input field must be a real vector");
  field.num = mxGetM(prhs[0])*mxGetN(prhs[0]); /* allow both row and column vectors */
  field.in = mxGetPr(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL); /* make output same orientation (row or column) as input */
  field.out = mxGetPr(plhs[0]);  

  if ((! mxIsDouble(prhs[1])) || mxIsComplex(prhs[1]) || ((mxGetN(prhs[1]) != 1) && (mxGetM(prhs[1]) != 1)))
    mexErrMsgTxt("coercive distribution quadrature points must be a real vector");
  if ((! mxIsDouble(prhs[2])) || mxIsComplex(prhs[2]) || ((mxGetN(prhs[2]) != 1) && (mxGetM(prhs[2]) != 1)))
    mexErrMsgTxt("coercive distribution quadrature weights must be a real vector");
  coercive.num = mxGetM(prhs[1])*mxGetN(prhs[1]);  /* allow both row and column vectors */
  if (coercive.num != mxGetM(prhs[2])*mxGetN(prhs[2])) 
    mexErrMsgTxt("There must be the same number of coercive distribution quadrature points and weights");
  if (coercive.num < 2)
    mexErrMsgTxt("At least 2 coercive distribution quadrature points are required");
  coercive.points = mxGetPr(prhs[1]);
  coercive.weights = mxGetPr(prhs[2]);

  /* skip 3 for now, since it depends on the whether we're using the relaxation model or not */
  if ((! mxIsDouble(prhs[4])) || mxIsComplex(prhs[4]) || ((mxGetN(prhs[4]) != 1) && (mxGetM(prhs[4]) != 1)))
    mexErrMsgTxt("effective distribution quadrature weights must be a real vector");
  effective.num = mxGetM(prhs[4])*mxGetN(prhs[4]);  /* allow both row and column vectors */
  if (effective.num < 2)
    mexErrMsgTxt("At least 2 effective distribution quadrature points/weights are required");
  if (effective.num % 2)
    mexErrMsgTxt("There must be an even number of effective distribution quadrature points.  Please change the quadrature method or number of intervals");
  effective.weights = mxGetPr(prhs[4]);  

  if ((! mxIsDouble(prhs[5])) || mxIsComplex(prhs[5]) || (mxGetN(prhs[5])*mxGetM(prhs[5]) != 1))
    mexErrMsgTxt("weightssum must be a real scalar");
  weightsum = mxGetScalar(prhs[5]);
  if ((! mxIsDouble(prhs[6])) || mxIsComplex(prhs[6]) || (mxGetN(prhs[6])*mxGetM(prhs[6]) != 1))
    mexErrMsgTxt("addit must be a real scalar");
  addit = mxGetScalar(prhs[6]);

  if ((! mxIsDouble(prhs[7])) || mxIsComplex(prhs[7]) || ((mxGetN(prhs[7]) != coercive.num) && (mxGetM(prhs[7]) != effective.num)))
    mexErrMsgTxt("xpos must be a real matrix of effective.num by coercive.num");
  xpos = mxGetPr(prhs[7]);
  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(effective.num, coercive.num, mxREAL);
    tmp = xpos;
    xpos = mxGetPr(plhs[1]);
    /* we copy xpos to the output matrix, and then work of that one */
    for (i = 0; i < coercive.num*effective.num; ++i) {
      xpos[i] = tmp[i];
    }
  } 

  if (nrhs == 8) {
    /* We still have to resolve the effective.points array */
    if ((! mxIsDouble(prhs[3])) || mxIsComplex(prhs[3]) || ((mxGetN(prhs[3]) != 1) && (mxGetM(prhs[3]) != 1)))
      mexErrMsgTxt("effective distribution quadrature points must be a real vector");
    if (effective.num != mxGetM(prhs[3])*mxGetN(prhs[3])) 
      mexErrMsgTxt("There must be the same number of effective distribution quadrature points and weights");
    effective.points = mxGetPr(prhs[3]);
    
    Neg_Relax_Hyst(field, coercive, effective, weightsum, addit, xpos);
  } else {
    /* Get the remaining parameters */
    if ((! mxIsDouble(prhs[3])) || mxIsComplex(prhs[3]) || (mxGetN(prhs[3])*mxGetM(prhs[3]) != 1))
      mexErrMsgTxt("effective distribution stepsize must be a real scalar");
    f_effective.num = effective.num;
    f_effective.weights = effective.weights;
    f_effective.step = mxGetScalar(prhs[3]);

    if ((! mxIsDouble(prhs[8])) || mxIsComplex(prhs[8]) || ((mxGetN(prhs[8]) != 1) && (mxGetM(prhs[8]) != 1)))
      mexErrMsgTxt("p_negpos must be a real vector");
    if ((! mxIsDouble(prhs[9])) || mxIsComplex(prhs[9]) || ((mxGetN(prhs[9]) != 1) && (mxGetM(prhs[9]) != 1)))
      mexErrMsgTxt("p_posneg must be a real vector");
    if ((! mxIsDouble(prhs[10])) || mxIsComplex(prhs[10]) || ((mxGetN(prhs[10]) != 1) && (mxGetM(prhs[10]) != 1)))
      mexErrMsgTxt("P_neg must be a real vector");
    if ((! mxIsDouble(prhs[11])) || mxIsComplex(prhs[11]) || ((mxGetN(prhs[11]) != 1) && (mxGetM(prhs[11]) != 1)))
      mexErrMsgTxt("P_pos must be a real vector");
    likes.num = mxGetM(prhs[8])*mxGetN(prhs[8]);  /* allow both row and column vectors */
    if ((likes.num != mxGetM(prhs[9])*mxGetN(prhs[9])) || (likes.num != mxGetM(prhs[10])*mxGetN(prhs[10])) ||
	(likes.num != mxGetM(prhs[11])*mxGetN(prhs[11])))
      mexErrMsgTxt("p_negpos, p_posneg, P_neg, and P_pos must all be the same length");
    likes.p_negpos = mxGetPr(prhs[8]);
    likes.p_posneg = mxGetPr(prhs[9]);
    likes.P_neg = mxGetPr(prhs[10]);
    likes.P_pos = mxGetPr(prhs[11]);
    
    if ((! mxIsDouble(prhs[12])) || mxIsComplex(prhs[12]) || (mxGetN(prhs[12])*mxGetM(prhs[12]) != 1))
      mexErrMsgTxt("delta t must be a real scalar");
    deltat = mxGetScalar(prhs[12]);
    if ((! mxIsDouble(prhs[13])) || mxIsComplex(prhs[13]) || (mxGetN(prhs[13])*mxGetM(prhs[13]) != 1))
      mexErrMsgTxt("increase must be a real integer scalar > 0");
    increase = (unsigned int)mxGetScalar(prhs[13]);
    if ((! mxIsDouble(prhs[14])) || mxIsComplex(prhs[14]) || (mxGetN(prhs[14])*mxGetM(prhs[14]) != 1))
      mexErrMsgTxt("increase must be a real integer scalar > 1");
    resolution = (unsigned int)mxGetScalar(prhs[14]);

    Relax_Hyst(field, coercive, f_effective, likes, weightsum, addit, deltat, increase, resolution, xpos);
  }
} /* end void mexFunction */
