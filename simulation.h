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
 $Id: simulation.h,v 1.8 2010-05-10 19:49:34 slampoud Exp $
*********************************************************************/
/*********************************************************************
Author: Sotiria Lampoudi
*********************************************************************/


#ifndef SIMULATION_H
#define SIMULATION_H

#include "SPA.h"

#define ARAND() ((double)random() / ((double)pow(2.0,31.0)-1.0))
#define RAND() (drand48())
#define SEED(s) (srand48(s))

struct simulation_stc
{
double SD;
double SNR;
double rate;
int max_frames;
int max_error_frames;
int quiet;
char fname[255];
char msg_fname[255];
int test_symmetric;
};

typedef struct simulation_stc SimulationOpts;

void InitSimulationOptions(SimulationOpts * SimOptions);

struct simulationstats_stc
{
int frame_total;
int frame_error_total;
int bit_error_total;
int symmetric_fail_total;
};

typedef struct simulationstats_stc SimulationStats;

void InitSimulationStats(SimulationStats * SimStats);

void ParseCmdline(int argc, char ** argv, SimulationOpts * SimOptions, 
                  SPAOpts * SPAoptions);

void FinalReport(SimulationOpts * SimOptions, SimulationStats * SimStats,
                 int bit_count, int check_count);

#endif
