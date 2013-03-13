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
 $Id: voutput.c,v 1.2 2010-05-05 05:25:32 slampoud Exp $
*********************************************************************/

/*********************************************************************
Authors: Rich Wolski, Sotiria Lampoudi
*********************************************************************/


#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

void PrintVector(FILE *fd, double *vec, int size)
{
	int i;
	for(i=0; i < size; i++)
	{
		fprintf(fd,"%f ",vec[i]);
	}
	fflush(fd);

	return;
}
void PrintVectorE(FILE *fd, double *vec, int size)
{
	int i;
	for(i=0; i < size; i++)
	{
		fprintf(fd,"%20.16e ",vec[i]);
	}
	fflush(fd);

	return;
}

void PrintVectorJustify(FILE *fd, double *vec, int size, int just)
{
	int i;
	int j;

	j = 0;
	for(i=0; i < size; i++)
	{
		if((i % just) == 0)
		{
			fprintf(fd,"(%d) ",i);
		}
		fprintf(fd,"%f ",vec[i]);
		if(j == (just - 1))
		{
			fprintf(fd,"\n");
			j = 0;
		}
		else
		{
			j++;
		}
	}
	fprintf(fd,"\n");
	fflush(fd);

	return;
}

void PrintVectorThresh(FILE *fd, double *vec, int size, double low, double
high)
{
	int i;
	for(i=0; i < size; i++)
	{
		if(vec[i] < low)
			fprintf(fd,"0 ");
		else if(vec[i] > high)
			fprintf(fd,"1 ");
		else
			fprintf(fd,"e ");
	}
	fflush(fd);

	return;
}

void Centroid(double *a, int dim, int *x, int *y)
{
	int i;
	int j;
	double x_acc = 0.0;
	double y_acc = 0.0;
	double xcount = 0.0;
	double ycount = 0.0;

	for(i=0; i < dim; i++)
	{
		for(j=0; j < dim; j++)
		{
			if(a[i*dim+j] > 0)
			{
				x_acc += (double)i;
				xcount += 1.0;
				y_acc += (double)j;
				ycount += 1.0;
			}
		}
	}

	*x = (int)(x_acc / xcount);
	*y = (int)(y_acc / ycount);

	return;
}


#ifdef PRINT_VEC
void PrintBinCounts(FILE *fd, int dim, double low, double high)
{
	double *a;
	int i;
	int j;
	int status;
	double snr;
	double *in;
	double *out;
	int size;
	int ones;
	int ees;
	int x;
	int y;

	status = ReadLogRecord(fd,&snr,&in,&out,&size);
	if(status == 0)
		return;

	rewind(fd);

	a = (double *)malloc(size*size*sizeof(double));
	if(a == NULL)
	{
		exit(1);
	}

	memset(a,0,size*size*sizeof(double));

	while(status = ReadLogRecord(fd,&snr,&in,&out,&size))
	{
		PrintVector(stdout,in,size);
		printf("\n");
		ones = 0;
		ees = 0;
		for(i=0; i < size; i++)
		{
			if(in[i] > high)
				ones++;
			else if(in[i] > low)
				ees++;
		}
		a[ones*dim+ees]++;
// printf("a[%d,%d]: %d\n",ones,ees,a[ones*50+ees]);
//   		if(ones > 50)
//  			PrintVector(stdout,in,size);
	}

	for(i=0; i < dim; i++)
	{
		fprintf(stdout,"(%d) ",i);
		for(j=0; j < dim; j++)
		{
			fprintf(stdout,"%d ",(int)a[i*dim+j]);
		}
		fprintf(stdout,"\n");
	}

	Centroid(a,dim,&x,&y);

	fprintf(stdout,"centroid: (%d, %d)\n",x,y);

	free(a);
	return;
}
			



#define ARGS "f:c:i:"
char *Usage = "print_vectors -f filename\n\
\t-i input_vec\n\
\t-c cut_off\n";

char Fname[255];
double CutOff;
int Justify;
int UseVec;

int main(int argc, char *argv[])
{
	int c;
	FILE *fd;
	int status;
	double *in;
	double *out;
	double snr;
	int size;
	int i;
	int biggies;

	CutOff = 0.7;

	while((c = getopt(argc,argv,ARGS)) != EOF)
	{
		switch(c)
		{
			case 'f':
				strncpy(Fname,optarg,sizeof(Fname));
				break;
			case 'c':
				CutOff = atof(optarg);
				break;
			default:
				fprintf(stderr,"bad command %c\n",(char)c);
				fprintf(stderr,"usage: %s",Usage);
				fflush(stderr);
				exit(1);
		}
	}

	if(Fname[0] == 0)
	{
		fprintf(stderr,"must specify filename\n");
		fprintf(stderr,"usage: %s",Usage);
		fflush(stderr);
		exit(1);
	}

	fd = fopen(Fname,"r");
	if(fd == NULL)
	{
		fprintf(stderr,"error opening %s\n",Fname);
		fflush(stderr);
		exit(1);
	}

	while(status = ReadLogRecord(fd,&snr,&in,&out,&size))
	{
		biggies = 0;
		for(i=0; i < size; i++)
		{
			if(out[i] >= CutOff)
			{
				fprintf(stdout,"%d ",i);
			}
		}
		fprintf(stdout,"\n");
		free(in);
		free(out);
	}

	fclose(fd);

	return(0);
}

#endif
