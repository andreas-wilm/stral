#include <stdio.h>
#include "mhash.h"

/* NAME and/or VERSION already defined in rnafold::config.h */
#define STRAL_NAME "stral"
#define STRAL_VERSION "0.3.4"


struct seqdata {
	float weight;
	float dummy;
	int format;
	int seqlength;
	int seqID;
	char *sequence;
	char *ID;
	char *orisequence;
};
/* size is 32 */

struct dist {
	float distance;
	float dummy1;
	float dummy2;
	float dummy3;
	int seqID1;
	int seqID2;
	struct dist *next_ptr;		/*pointer to next vector */
	struct dist *previous_ptr;
};

float **distmatrix;                     /* distance matrix derived from pairwise alignment*/
float **array;
struct profile *prof;
struct seqdata *sdata;

float alpha;
float gapOpenM;
float gapExtM ;
float gapOpen ;
float gapExt  ;
int   njoin   ;
int   positive;
int   ribosum ;
int   weighted;
int   rooted  ;
int   verbose ;
int	 ppair   ;
int   bion    ;
int   weigh   ;
int   adjust  ;
int   upgma   ;
int   nofile  ;
int   gapmodel;
int 	 noterm  ;
int 	 palign  ;
int 	 fft     ;
int   spread  ;
int   submat  ;

/*double temperature;*/
int nseqs; /* no of sequences in input-file */

float msascore;

struct vector **lastcpyptr;
struct vector **firstcpyptr;

char workingdir[200];
char *seqfile;

MHASH td;

int ct; 	
