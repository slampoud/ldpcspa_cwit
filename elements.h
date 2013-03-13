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
 $Id: elements.h,v 1.4 2010-05-05 05:25:32 slampoud Exp $
*********************************************************************/

/*********************************************************************
Authors: Rich Wolski, Sotiria Lampoudi
*********************************************************************/

#ifndef ELEMENTS_H
#define ELEMENTS_H

#define MAXEDGES 50

struct edge_stc
{
	double d_mu;
	double d_nu;
	struct node_stc *bit;
	struct node_stc *check;
};

typedef struct edge_stc Edge;

struct node_stc
{
	int id;
	int max_edge;
	double p;
        double u_l_hat; /* opinion of value: odds of 1, p(1)/p(0)  */
	Edge *edges[MAXEDGES];
};

typedef struct node_stc Node;

double ConvertSNR(double snr, double rate);
void print_undecoded(Node ** bit_array, int bit_count, double SD);

Node *MakeNode(int id);
void LinkNodes(Node *check, Node *bit);
void PrintGraph(Node **check_array, int max_checks);
char *MakeGraphString(Node **check_array, int max_checks,
		Node **bit_array, int bit_count, int size);

void ResetBitProbs(Node **bit_array, int bit_count, double sd);
void ReadBitProbs(Node **bit_array, int bit_count, FILE * FH);
void FlipBitProbs(Node **bit_array, int bit_count);

void PrintGraph(Node **check_array, int max_checks);

void ParseInput(char *fname, 
                Node ***bits, int *bit_count,
                Node ***checks, int *check_count);

void MakeGraph(char *graph_string, 
                Node ***bits, int *bit_count,
                Node ***checks, int *check_count);

void FreeGraph(Node **checks, int check_count, Node **bits, int bit_count);

#endif
