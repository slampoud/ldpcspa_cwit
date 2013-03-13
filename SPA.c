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
 $Id: SPA.c,v 1.15 2010-05-11 01:47:29 slampoud Exp $
*********************************************************************/
/*********************************************************************
Author: Sotiria Lampoudi
*********************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "elements.h"
#include "SPA.h"

#define NEWS 1

void InitSPAOptions(SPAOpts * SPAoptions){
  SPAoptions->zero_threshold = 0.001;
  SPAoptions->inf_threshold = 1000.0;
  SPAoptions->max_spa_iters = 500;
}

double checkmult(double i, double j){
 if( (i == 0) && isinf(j) ){
  return 1.0;
 } else if( isinf(i) && (j==0)){
  return 1.0;
 } else {
  return i*j;
 } 
}

double sfunc(double i){
#if NEWS
  if(isinf(i)){
    return -1;
  } else {
    return ((1.0-i)/(1.0+i));
  }
#else
return ((1.0-i)/(1.0+i));
#endif
}

int Bit2Check(Node *bit)
{
  int i,j;
  double prod;
  double rho_p;
  rho_p = bit->p/(1.0-bit->p);

  for(j=0; j < bit->max_edge; j++){
    prod=1.0;
    for(i=0; i<bit->max_edge; i++){
      if( i == j ) continue;
      prod = checkmult(prod, sfunc(bit->edges[i]->d_nu));
      //prod *= sfunc(bit->edges[i]->d_nu);
    }
    //bit->edges[j]->d_mu = sfunc(rho_p*prod);
    bit->edges[j]->d_mu = sfunc(checkmult(rho_p, prod));
  } 
  return(-1);
}

int Check2Bit(Node *check)
{
  int i,j;
  double prod;
  
  for(j=0; j < check->max_edge; j++){
    prod=1.0;
    for(i=0; i < check->max_edge; i++){
      if( i == j) continue;
      prod *= check->edges[i]->d_mu;
    }
    check->edges[j]->d_nu = prod;
  }
  return(-1);
}

double CalcNewEstimate(Node * bit){
  int i;
  double prod = 1.0;
  for(i=0; i < bit->max_edge; i++){
	//prod *= sfunc(bit->edges[i]->d_nu);
        prod = checkmult(prod, sfunc(bit->edges[i]->d_nu));
  }
  //prod *= bit->p/(1.0-bit->p);
  prod = checkmult(prod, bit->p/(1.0-bit->p));
  bit->u_l_hat=prod;
  return (bit->u_l_hat);
}

int SPASingle(Node ** bit_array, int bit_count,
	      Node ** check_array, int check_count){
  int i;
  int ret;
  /* do check -> bit */
  for(i=0; i<check_count; i++){
    ret = Check2Bit(check_array[i]);
    if (ret >= 0){
      fprintf(stderr, "SPA halted due to exhausted numerical precision in check %d edge %d.\n", i, ret);
      return(0);
    }
  }
  /* do bit -> check */
  for(i=0; i<bit_count; i++){
    ret = Bit2Check(bit_array[i]);
    if (ret >= 0){
      fprintf(stderr, "SPA halted due to exhausted numerical precision in bit %d edge %d.\n", i, ret);
      return(0);
    }
  }
  return(1);
}

int SPAMultiple(SPAOpts * SPAoptions,
                 Node ** bit_array, int bit_count,
                 Node ** check_array, int check_count){
  int do_another_iter;
  int cur_spa_iter;
  int ret, i;
  double u_l_hat;

  do_another_iter=1;
  for(cur_spa_iter=0; cur_spa_iter < SPAoptions->max_spa_iters; cur_spa_iter++)\
    {
      if(do_another_iter){
	do_another_iter=0;
	ret = SPASingle(bit_array, bit_count,
			check_array, check_count);
	if(ret == 0){
	  /*                                                                      
           error condition: a d_nu is 1 (p(0)=1), or -1 (p(1)=1)                
           if we continue past here we will end up with NaN's                   
           so we are forced to stop.                                            
	  */
	  break;
	}
	/*                                                                        
         check if bits' opinions of themselves have all converged               
	*/
	for(i = 0; i<bit_count; i++){
	  u_l_hat = CalcNewEstimate(bit_array[i]);
	  if( (u_l_hat > SPAoptions->zero_threshold) &&
	      (u_l_hat < SPAoptions->inf_threshold) ){
	    do_another_iter++;
	  }
	} /* checked all bits */
      } else {
	break; /* we've decoded, get out */
      }
    }/* spa iteration check */
return cur_spa_iter;
}

int CountBitErrors(Node ** bit_array, int bit_count, double shouldbe){
     int bit_error_count = 0;
     int i;
     if (shouldbe == 0.0){
       for(i = 0; i<bit_count; i++){
	 if(bit_array[i]->u_l_hat > 1.0 ){ /* p(1) > p(0) */
	   bit_error_count++;
	 }
       }
     } else {
       for(i = 0; i<bit_count; i++){
	 if(bit_array[i]->u_l_hat < 1.0 ){ /* p(1) < p(0) */
	   bit_error_count++;
	 }
       }
     }
     return bit_error_count;
}
