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
 $Id: elements.c,v 1.5 2010-05-10 16:46:51 slampoud Exp $
*********************************************************************/

/*********************************************************************
Authors: Rich Wolski, Sotiria Lampoudi
*********************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

#include "elements.h"
#include "voutput.h"
#include "SPA.h"
#include "norm.h"

#define ARAND() ((double)random() / ((double)pow(2.0,31.0)-1.0))
#define RAND() (drand48())
#define SEED(s) (srand48(s))

double ConvertSNR(double snr, double rate)
{
  double sd;
  sd = (1.0 / sqrt(2*rate)) * pow(10.0,-1.0*snr/20.0);
  return(sd);
}

void print_undecoded(Node ** bit_array, int bit_count, double SD){
  double * out;
  double * in;
  int i;
  
  out = (double *)malloc(bit_count*sizeof(double));
  if(out == NULL)
      exit(1);

  in = (double *)malloc(bit_count*sizeof(double));
  if(in == NULL)
    exit(1);

  for(i=0; i < bit_count; i++){
    /* p(1) opinion */
    if(isinf(bit_array[i]->u_l_hat)){
      out[i] = 1.0;
    } else {
      out[i] = bit_array[i]->u_l_hat / (1.0 + bit_array[i]->u_l_hat); 
    }
    /* initial p(1) */
    in[i] = bit_array[i]->p;
  }
  fprintf(stdout, "%f | ", SD);
  PrintVector(stdout,in,bit_count);
  fprintf(stdout, "| ");
  PrintVector(stdout,out,bit_count);
  fprintf(stdout, "\n");

}

double BayesRatio(double p, double sigma);

double X_0(double p_0, double sigma)
{
	double r;

	r = ((sigma*sigma)/2.0)*log((1.0-p_0)/p_0);


	return(r);
}
	

double GaussOfOne(double x, double sd)
{
	double p;

	p = InvNormal(x,1.0,sd);
	return(p);
}

/*
 * computes the cond probability of a 1 given p < q
 */
double CondRANDGT(double sd, double q)
{
        double r = RAND();
        double upper;
        double c;

        c = ((sd * sd) / 2.0) * log((1.0 - q) / q);

        upper = NormalCDF(c,1.0,sd);

        /*
         * U(0,c)
         */
        return(r * upper);
}


double CondRANDLT(double sd, double q)
{
        double r = RAND();
        double upper;
        double c;

        c = ((sd * sd) / 2.0) * log((1.0 - q) / q);

        upper = NormalCDF(c,1.0,sd);

        /*
         * U(0,c)
         */
        return(r * (1.0 - upper) + upper);
}

double NChooseK(double n, double k)
{
        double *b;
        int i;
        int j;
        double r;

        b = (double *)malloc((n+1)*sizeof(double));
        if(b == NULL)
        {
                exit(1);
        }

        b[0] = 1.0;
        for(i=1; i <= n; ++i)
        {
                b[i] = 1;
                for(j=i-1; j > 0; --j)
                {
                        b[j] += b[j-1];
                }
        }

        r = b[(int)k];
        free(b);
        return(r);
}

double BayesRatio(double x, double sd)
{
	double b;
	double f_0;
	double f_1;

	f_1 = Normal(x,-1.0,sd);
	f_0 = Normal(x,1.0,sd);

	b = f_1 / (f_1+f_0);

	return(b);
}

Node *MakeNode(int id)
{
	Node *n;

	n = (Node *)malloc(sizeof(Node));
	if(n == NULL)
		exit(1);
	memset(n,0,sizeof(Node));

	n->id = id;
    	n->p = 0.0;
        n->u_l_hat = 0.0;

	return(n);
}

void FreeNode(Node *n)
{
	free(n);
}

void SetNodeProb(Node *bit, double p)
{
  bit->p = p;
  bit->u_l_hat=p/(1.0-p);
  return;
}

void SetEdgeProb(Edge *e, double p)
{
  e->d_mu = (1.0-p)-p;
  e->d_nu = 0;
  return;
}

Edge *MakeEdge()
{
	Edge *e;

	e = (Edge *)malloc(sizeof(Edge));
	if(e == NULL)
		exit(1);

	memset(e,0,sizeof(Edge));

	return(e);
}

void FreeEdge(Edge *e)
{
	free(e);
}

void LinkNodes(Node *check, Node *bit)
{
	Edge *e;
	int i;

	e = MakeEdge();
	if(e == NULL)
		exit(1);
	e->bit = bit;
	e->check = check;

	i = bit->max_edge;
	bit->edges[i] = e;
	bit->max_edge = i+1;
	if(i >= MAXEDGES)
	{
		fprintf(stderr,"LinkNodes: too many src edges\n");
		fflush(stderr);
		exit(1);
	}

	i = check->max_edge;
	check->edges[i] = e;
	check->max_edge = i+1;
	if(i >= MAXEDGES)
	{
		fprintf(stderr,"LinkNodes: too many dst edges\n");
		fflush(stderr);
		exit(1);
	}

	return;
}

void PrintGraph(Node **check_array, int max_checks)
{
	int i;
	int j;

	for(i=0; i < max_checks; i++)
	{
		for(j=0; j < check_array[i]->max_edge; j++)
		{
			printf(" %d:%d",
				check_array[i]->id,
				check_array[i]->edges[j]->bit->id);
		}
	}

	return;
}


int GetIntStr(FILE *fd, char *str, int ssize)
{
	int i = 0;
	int c;

	c = fgetc(fd);
	if(c == EOF)
	{
		return(EOF);
	}

	memset(str,0,ssize);

	while((c != EOF) && !isdigit(c))
	{
		c = fgetc(fd);
		if(c == EOF)
			return(EOF);
	} 

	while(isdigit(c) && (i < ssize))
	{
		str[i] = c;
		i++;
		c = fgetc(fd);
		if(c == EOF)
		{
			fprintf(stderr,"GetIntStr EOF\n");
			fflush(stderr);
			exit(1);
		}
	}

	return(1);
}


int ParseIntStr(char *s, char *str, int ssize, char **next_str)
{
	int i = 0;
	char *c;

	c = s;
	if(*c == 0)
	{
		return(EOF);
	}

	memset(str,0,ssize);

	while((*c != 0) && !isdigit(*c))
	{
		c++;
		if(*c == 0)
			return(EOF);
	} 

	while(isdigit(*c) && (i < ssize))
	{
		str[i] = *c;
		i++;
		c++;
		if(*c == 0)
		{
			break;
		}
	}

	*next_str = c;

	return(1);
}


void ResetBitProbs(Node **bit_array, int bit_count, double sd)
{
	double prob = sd;
	int i;
	int j;

	for(i=0; i < bit_count; i++){
	  prob = GaussOfOne(RAND(),sd);
	  SetNodeProb(bit_array[i], BayesRatio(prob,sd));
	  for(j=0; j < bit_array[i]->max_edge; j++){
	    SetEdgeProb(bit_array[i]->edges[j],bit_array[i]->p);
	  }
	}
	return;
}

void FlipBitProbs(Node **bit_array, int bit_count)
{
  int i;
  int j;
  double p;
  for(i=0; i < bit_count; i++){
    p = 1.0-bit_array[i]->p;
    SetNodeProb(bit_array[i], p);
    for(j=0; j < bit_array[i]->max_edge; j++){
      SetEdgeProb(bit_array[i]->edges[j],bit_array[i]->p);
    }
  }
  return;
}

void ReadBitProbs(Node **bit_array, int bit_count, FILE * FH){
  double prob = -1;
  int i;
  int j;
  fscanf(FH, " ");
  for(i=0; i < bit_count; i++){
    fscanf(FH, "%lf ", &prob);
    SetNodeProb(bit_array[i], prob);
    for(j=0; j < bit_array[i]->max_edge; j++){
      SetEdgeProb(bit_array[i]->edges[j],bit_array[i]->p);
    }
  }
  fscanf(FH, "\n");
  return;
}


/*
 * columns are bits
 * rows are checks
 *
 * first two numbers are row and col dimenstions
 * row:col indicates that (row,col) element is connected
 */
void ParseInput(char *fname, 
		Node ***bits, int *bit_count, 
	     	Node ***checks, int *check_count)
{
	FILE *fd;
	int max_bits;
	int max_checks;
	Node **bit_array;
	Node **check_array;
	char num_string[1024];
	int c;
	int row;
	int col;
	int i;

	fd = fopen(fname,"r");
	if(fd == NULL)
	{
		fprintf(stderr,"can't open file %s\n",fname);
		fflush(stderr);
		exit(1);
	}

	/*
	 * get the row number
	 */
	
	GetIntStr(fd,num_string,sizeof(num_string)); 
	max_checks = atoi(num_string);

	/*
	 * get the col number
	 */
	GetIntStr(fd,num_string,sizeof(num_string)); 
	max_bits = atoi(num_string);


	bit_array = (Node **)malloc(max_bits*sizeof(Node *));
	if(bit_array == NULL)
		exit(1);
	check_array = (Node **)malloc(max_checks*sizeof(Node *));
	if(check_array == NULL)
		exit(1);

	for(i=0; i < max_checks; i++)
		check_array[i] = NULL;
	for(i=0; i < max_bits; i++)
		bit_array[i] = NULL;

	while(!feof(fd))
	{
		c = GetIntStr(fd,num_string,sizeof(num_string));
		if(c == EOF)
			break;
		row = atoi(num_string);
		c = GetIntStr(fd,num_string,sizeof(num_string));
		if(c == EOF)
		{
			fprintf(stderr,"EOF between row:col (%d)\n",
				row);
			fflush(stderr);
			exit(1);
		}
		col = atoi(num_string);
		if((row < 0) || (col < 0))
		{
			fprintf(stderr,"bad row:col -- %d:%d\n",row,col);
			fflush(stderr);
			continue;
		}
		if((row >= max_checks) || (col >= max_bits))
		{
			fprintf(stderr,"bad row:col -- %d:%d\n",row,col);
			fflush(stderr);
			continue;
		}
		if(check_array[row] == NULL)
		{
			check_array[row] = MakeNode(row);
		}
			
		if(bit_array[col] == NULL)
		{
			bit_array[col] = MakeNode(col);
		}

		LinkNodes(check_array[row],bit_array[col]);
	}


	*bits = bit_array;
	*checks = check_array;
	*bit_count = max_bits;
	*check_count = max_checks;

	return;

}


void MakeGraph(char *graph_string, 
		Node ***bits, int *bit_count, 
	     	Node ***checks, int *check_count)
{
	int max_bits;
	int max_checks;
	Node **bit_array;
	Node **check_array;
	char num_string[1024];
	int c;
	int row;
	int col;
	int i;
	char *next_str;


	/*
	 * get the row number
	 */
	
	next_str = graph_string;
	ParseIntStr(next_str,num_string,sizeof(num_string),&next_str); 
	max_checks = atoi(num_string);

	/*
	 * get the col number
	 */
	ParseIntStr(next_str,num_string,sizeof(num_string),&next_str); 
	max_bits = atoi(num_string);


	bit_array = (Node **)malloc(max_bits*sizeof(Node *));
	if(bit_array == NULL)
		exit(1);
	check_array = (Node **)malloc(max_checks*sizeof(Node *));
	if(check_array == NULL)
		exit(1);

	for(i=0; i < max_checks; i++)
		check_array[i] = NULL;
	for(i=0; i < max_bits; i++)
		bit_array[i] = NULL;

	while(*next_str != 0)
	{
		c = ParseIntStr(next_str,num_string,sizeof(num_string),
				&next_str);
		if(c == EOF)
			break;
		row = atoi(num_string);
		c = ParseIntStr(next_str,num_string,sizeof(num_string),
			&next_str);
		if(c == EOF)
		{
			fprintf(stderr,"EOF between row:col (%d)\n",
				row);
			fflush(stderr);
			exit(1);
		}
		col = atoi(num_string);
		if((row < 0) || (col < 0))
		{
			fprintf(stderr,"bad row:col -- %d:%d\n",row,col);
			fflush(stderr);
			continue;
		}
		if((row >= max_checks) || (col >= max_bits))
		{
			fprintf(stderr,"bad row:col -- %d:%d\n",row,col);
			fflush(stderr);
			continue;
		}
		if(check_array[row] == NULL)
		{
			check_array[row] = MakeNode(row);
		}
			
		if(bit_array[col] == NULL)
		{
			bit_array[col] = MakeNode(col);
		}

		LinkNodes(check_array[row],bit_array[col]);
	}


	*bits = bit_array;
	*checks = check_array;
	*bit_count = max_bits;
	*check_count = max_checks;

	return;

}
void FreeGraph(Node **check_array, int max_checks, Node **bit_array, int
	       max_bits)
{
	int i;
	int j;

	for(i=0; i < max_bits; i++)
	{
		for(j=0; j < bit_array[i]->max_edge; j++)
		{
			FreeEdge(bit_array[i]->edges[j]);
		}
		FreeNode(bit_array[i]);
	}

	for(i=0; i < max_checks; i++)
	{
		FreeNode(check_array[i]);
	}

	free(check_array);
	free(bit_array);

	return;
}

