#include <stdio.h>

#include "bionj.h"

/* defining data structure for probabilty records*/
/* double-linked-list*/
struct vector {
  double p0;				/* probabilty unpaired */					
  double p1;				/* downstream pairing probabilty */
  double p2;				/* upstream pairing probabilty */
  double dotbracket;
  double dummy2;
  int seqID;				/* sequence ID from *.vie file */
  int replbase;			/* 0 1 2 3 4 5 */
  int pos;				/* position in sequence*/
  int oripos;				/* nt position in dealigned sequence */
  char base;				/* A C G U N - */
  char oribase;			/* original characters from input file*/
  struct vector *next_ptr;	/*pointer to next vector */
  struct vector *previous_ptr;
};


/* storing clustering information from guide trees*/
struct profile {
	int *seqprofA;
	int nA;
	int *seqprofB;
	int nB;
	int typeA;
	int typeB;
	int pnoA;
	int pnoB;
};
	


/*double linked list to hold the score values of one group */
struct lscore {
	float score;
	struct lscore *next_ptr;		
	struct lscore *previous_ptr;
	int pos;
};

/* FIXME: remove ????? */
struct profilearray {
	int a;
	int c;
	int g;
	int u;
	int n;
};

/* pointers to first resp. last element of sequences to be aligned */
struct vector **last_ptr;
struct vector **first_ptr;

/* pointers to first resp. last element of one group already aligned */
struct lscore **lps;
struct lscore **fps;

/*storing position specific profiles */
float **arrayI;
float **arrayJ;



/* functions */
void align(int nseqs, char *seqfile);
void malign(int n);
void method(struct vector *first1_ptr, struct vector *last1_ptr, struct vector *first2_ptr, struct vector *last2_ptr, int nseqs, char *seqfile);
float methodmult(struct profile p, int lenA, int lenB, int nseqs, int cp);
void levelpos(int cluster);
void print_data_content2(struct vector *current_ptr);
float **data_content(struct vector *current_ptr, float **arrayX);
