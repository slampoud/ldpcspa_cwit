/*********************************************************************
Copyright (c) 2007, Regents of the University of California

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.  Neither the name of the
University of California, Santa Barbara nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.  In addition, permission to use, copy,
modify, and distribute this software and its documentation for educational,
research and non-profit purposes, without fee, and without a written agreement
is hereby granted.  Permission to incorporate this software into commercial
products may be obtained by contacting the University of California, Santa
Barbara's Sponsored Projects Office through

SHERYLLE MILLS ENGLANDER
Director
Office of Technology & Industry Alliances
552 University Ave.
Trailer #342
University of California, Santa Barbara
Santa Barbara, CA 93106
Mail Code: 2055
englander@research.ucsb.edu

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/
/*********************************************************************
 $Id: SPA.h,v 1.7 2010-05-10 23:31:17 slampoud Exp $
*********************************************************************/
/*********************************************************************
Author: Sotiria Lampoudi
*********************************************************************/


#ifndef SPA_H
#define SPA_H

struct spa_stc
{
  double zero_threshold;
  double inf_threshold;
  int max_spa_iters;
};

typedef struct spa_stc SPAOpts;

void InitSPAOptions(SPAOpts * SPAoptions);

double sfunc(double i);
int Bit2Check(Node *bit);
int Check2Bit(Node *check);
double CalcNewEstimate(Node * bit);
int SPASingle(Node ** bit_array, int bit_count,
	      Node ** check_array, int check_count);
int SPAMultiple(SPAOpts * SPAoptions, 
		 Node ** bit_array, int bit_count,
		 Node ** check_array, int check_count);
int CountBitErrors(Node ** bit_array, int bit_count, double shouldbe);
#endif

