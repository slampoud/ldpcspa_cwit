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
 $Id: simulation.c,v 1.16 2010-05-11 01:47:29 slampoud Exp $
*********************************************************************/
/*********************************************************************
Author: Sotiria Lampoudi
*********************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "elements.h"
#include "SPA.h"
#include "norm.h"
#include "voutput.h"
#include "simulation.h"

#define ARGS "g:s:I:M:E:F:S:Z:R:Q1"
char *Usage = "simulation [options]\n\
Code options:\n\
\t-g [filename], graph file (required)\n\
\t-S [double], standard deviation/ -s [double], SNR (one required)\n\
\t-R [double], rate of the code (required)\n\
SPA options:\n\
\t-Z [double], zero threshold (default=0.001)\n\
\t-I [int], number of max SPA iterations (default=500)\n\
Simulation options:\n\
\t-F [int], terminate after this many total frames (required)\n\
\t-E [int], terminate after this many frame errors, unless -F is exceeded (optional)\n\
\t-M [filename], file containing messages to decode (optional)\n\
\t-1 , test for symmetric decoding of the 1-vector\n\
\t-Q , quiet mode (default=off)\n";

int main(int argc, char *argv[])
{
  /* Internal */
  Node **check_array; 
  int check_count; 
  Node **bit_array;
  int bit_count;
  int frame_index=0;
  int bit_error_count=0;
  int bit_error_count_sym=0;
  struct timeval tm;
  FILE * msg_file = 0;

  SPAOpts SPAoptions;
  SimulationStats SimStats;
  SimulationOpts SimOptions;

  InitSPAOptions(&SPAoptions);
  InitSimulationOptions(&SimOptions);
  InitSimulationStats(&SimStats);

  gettimeofday(&tm,NULL);
  SEED(tm.tv_sec+tm.tv_usec);
  
  /* get the cmdline args and make sure they are sane */
  ParseCmdline(argc, argv, &SimOptions, &SPAoptions);

  /* read in graph */
  ParseInput(SimOptions.fname, &bit_array, &bit_count, 
	     &check_array, &check_count);

  /* if messages are coming from a file, open it */
  if( SimOptions.msg_fname[0] != 0){
    msg_file = fopen(SimOptions.msg_fname, "r");
  }
  
  /* Start simulation */
  for(frame_index=0; frame_index < SimOptions.max_frames; frame_index++){

    if(SimStats.frame_error_total < SimOptions.max_error_frames){
      /* initialize probabilities */
      /* randomly */
      if( SimOptions.msg_fname[0] == 0){ 
	ResetBitProbs(bit_array, bit_count, SimOptions.SD);
      } else { /* or from file */
	ReadBitProbs(bit_array, bit_count, msg_file);
      }
      
      SPAMultiple(&SPAoptions, bit_array, bit_count, 
		  check_array, check_count);
      
      bit_error_count = CountBitErrors(bit_array, bit_count, 0.0);

      if(bit_error_count){
	SimStats.frame_error_total++;
	SimStats.bit_error_total+=bit_error_count;
	if(! SimOptions.quiet){
	  /* print out error */
	  print_undecoded(bit_array, bit_count, SimOptions.SD);
	}
      }
      /* Symmetry test: check that the all ones vector would have the 
         same number of bit errors
       */ 
      if(SimOptions.test_symmetric){
	/* flip probabilities */
	FlipBitProbs(bit_array, bit_count);
	/* do SPA */
	SPAMultiple(&SPAoptions, bit_array, bit_count, 
		    check_array, check_count);
	/* count # of bit errors */      
	bit_error_count_sym = CountBitErrors(bit_array, bit_count, 1.0);
	/* check if it agrees with the inverse */	
	if(bit_error_count_sym != bit_error_count){
          SimStats.symmetric_fail_total++;
	  fprintf(stderr, "Failed the symmetry test: \
                  zero-vector had %d errors, one-vector had %d errors\n", 
                  bit_error_count, bit_error_count_sym);
	  fflush(stderr);
          if(! SimOptions.quiet){
	    print_undecoded(bit_array, bit_count, SimOptions.SD);
          }
	}
      }
      SimStats.frame_total=frame_index+1;
    }  /* SimStats.frame_error_total check */
    else { /* we have reached max_error_frames */
      break;
    }
  } /* total frame count check */
  /* Generate report */
  FinalReport(&SimOptions, &SimStats, bit_count, check_count);

  return(0);
}

void FinalReport(SimulationOpts * SimOptions, SimulationStats * SimStats, 
                 int bit_count, int check_count){
  fprintf(stdout, "SNR: %f SD: %f total_frames: %d \
          BER: %e bit_errors: %d frame_errors: %d\n",
          SimOptions->SNR,
          SimOptions->SD,
          SimStats->frame_total,
          ((double) SimStats->bit_error_total)/
              ((double) SimStats->frame_total * (double)bit_count),
          SimStats->bit_error_total,
          SimStats->frame_error_total);
  if(SimOptions->test_symmetric){
    fprintf(stdout, "Symmetric failure count: %d\n", SimStats->symmetric_fail_total);
  }
}

void ParseCmdline(int argc, char ** argv, 
		  SimulationOpts * SimOptions, 
		  SPAOpts * SPAoptions){
  int c;
  while((c = getopt(argc,argv,ARGS)) != EOF)
    {
      switch(c)
	{
        case '1': 
          SimOptions->test_symmetric = 1;
          break;
	case 'g':
	  strncpy(SimOptions->fname,optarg,sizeof(SimOptions->fname));
	  break;
	case 'S':
	  SimOptions->SD = atof(optarg);
	  break;
	case 'R':
	  SimOptions->rate = atof(optarg);
	  break;
	case 's':
	  SimOptions->SNR = atof(optarg);
	  break;
	case 'Z':
	  SPAoptions->zero_threshold = atof(optarg);
	  SPAoptions->inf_threshold = 1.0 / SPAoptions->zero_threshold;
	  break;
	case 'I':
	  SPAoptions->max_spa_iters = atoi(optarg);
	  break;
	case 'M':
	  strncpy(SimOptions->msg_fname,optarg,sizeof(SimOptions->msg_fname));
	  break;
	case 'F':
	  SimOptions->max_frames = atoi(optarg);
	  break;
	case 'E':
	  SimOptions->max_error_frames = atoi(optarg);
	  break;
	case 'Q':
	  SimOptions->quiet = 1;
	  break;
	default:
	  fprintf(stderr,"unrecognized option %c",
		  (char)c);
	  fprintf(stderr,"%s",Usage);
	  fflush(stderr);
	  exit(1);
	}
    }
  
  /* 
     Check that all required args have been specified 
     and options make sense
  */
  if(SimOptions->fname[0] == 0){
    fprintf(stderr,"You must specify a graph file\n");
    fprintf(stderr,"%s",Usage);
    fflush(stderr);
    exit(1);
  }
  if((SimOptions->rate == -1) && (SimOptions->SNR > 0)){
    fprintf(stderr, "You must specify a rate to use SNR\n");
    fprintf(stderr, "%s", Usage);
    fflush(stderr);
    exit(1);
  }
  if(SimOptions->SNR > 0){
    SimOptions->SD = ConvertSNR(SimOptions->SNR, SimOptions->rate);
  }
  if(SimOptions->SD == -1){
    fprintf(stderr, "You must specify a SD or SNR & rate\n");
    fprintf(stderr,"%s",Usage);
    fflush(stderr);
    exit(1);
  } 
  if(SimOptions->max_error_frames == -1){
    SimOptions->max_error_frames=SimOptions->max_frames;
  }
  if(SimOptions->max_error_frames > SimOptions->max_frames){
    fprintf(stderr, "The number of total frames must equal or exceed the number of error frames\n");
    fprintf(stderr, "%s", Usage);
    fflush(stderr);
    exit(1);
  } 
}

void InitSimulationOptions(SimulationOpts * SimOptions){
  SimOptions->SD=-1;
  SimOptions->SNR=-1;
  SimOptions->rate=-1;
  SimOptions->max_frames=-1;
  SimOptions->max_error_frames=-1;
  SimOptions->quiet=0;
  SimOptions->test_symmetric=0;
  memset(SimOptions->fname,0,sizeof(SimOptions->fname));
  memset(SimOptions->msg_fname,0,sizeof(SimOptions->msg_fname));
}

void InitSimulationStats(SimulationStats * SimStats){
  SimStats->frame_total=0;
  SimStats->frame_error_total = 0;
  SimStats->bit_error_total = 0;
  SimStats->symmetric_fail_total = 0;
}

