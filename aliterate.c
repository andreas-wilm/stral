#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>

/* StrAl */
#include "align.h"
#include "alfuncs.h"
#include "aliterate.h"
#include "create_matrix.h"
#include "stral.h"
#include "nrutils.h"
#include "iofuncs.h"




/* Function:		methodmultcpy
 * 
 * Purpose:		do the multiple sequence alignment
 *           
 * Args:			p: contains the groups to be aligned
 *				lenA, lenB: length of the sequences in group I, J
 *				nseqs: number of sequences in alignment
 *
 * Returns:		float 
 */
float methodmultcpy(struct profile p, int lenA, int lenB, int nseqs){
	int i = 0;
	int j = 0;
	int k = 0;
	int m = 0;
	int n = 0;
	int gapa=0;
	int gapb=0;
	int maxRow, maxCol;
	struct vector **ptrarrayA;
	struct vector **ptrarrayB;
	float **v;    
	float **f;
	float **e;
	float **g;
	float **scorematrix;
	float va,ga,ea,fa,maxxV=0;;
	char *backtrack;
	float sumWa =0;
	float sumWb =0;
	struct lscore *spA;
	struct lscore *spB;
	
#ifdef DEBUG
	if (!verbose){
		printf("\nlenA =  %d, lenB  =   %d\n", lenA,lenB);
	}
#endif
	
	ct = 0;
	
	ptrarrayA = (struct vector **)malloc(p.nA*sizeof(struct vector *));
	ptrarrayB = (struct vector **)malloc(p.nB*sizeof(struct vector *));
	
	scorematrix=(float **)malloc((lenA+1)*sizeof(float *));
	
	arrayI = (float **)calloc(lenA+1,sizeof(float *));
	arrayJ = (float **)calloc(lenB+1,sizeof(float *));
	
	
	v = (float **)malloc((lenA+1)*sizeof(float*));
	e = (float **)malloc((lenA+1)*sizeof(float*));
	f = (float **)malloc((lenA+1)*sizeof(float*));
	g = (float **)malloc((lenA+1)*sizeof(float*));
	
	
	
	for (k = 0; k <= lenA; ++k){
		arrayI[k] = (float *)calloc(9,sizeof(float));
		scorematrix[k] = (float *)malloc((lenB+1)*sizeof(float));
		v[k] = (float *)malloc((lenB+1)*sizeof(float));
		e[k] = (float *)malloc((lenB+1)*sizeof(float));
		f[k] = (float *)malloc((lenB+1)*sizeof(float));
		g[k] = (float *)malloc((lenB+1)*sizeof(float));
	
		
	}
	for (k = 0; k <= lenB; ++k){
		arrayJ[k] = (float *)calloc(9,sizeof(float));
	}
	for (k = 0; k < p.nA;++k ){
		ptrarrayA[k] = lastcpyptr[p.seqprofA[k]-1];
		data_content(lastcpyptr[p.seqprofA[k]-1],arrayI);
		sumWa += sdata[ptrarrayA[k]->seqID].weight; 
	}
	
	for (k = 0; k < p.nB;++k ){
		ptrarrayB[k] = lastcpyptr[p.seqprofB[k]-1];
		data_content(lastcpyptr[p.seqprofB[k]-1],arrayJ);
		sumWb += sdata[ptrarrayB[k]->seqID].weight;
	}
	
	spA=NULL;
	spB=NULL;
	
	if (lenA < lenB){
		for (i=1;i<=lenA;++i){	
			v[i][0]=-i*gapOpenM;
			f[i][0]=-i*gapOpenM;
			e[i][0]=-i*gapOpenM;
			g[i][0]=-i*gapOpenM;
		}
	
		for (j=0;j<=lenB;++j){
			v[0][j]=0;
			f[0][j]=0;
			g[0][j]=0;
			e[0][j]=0;
		}
	} else if	(lenA > lenB){
		for (i=0;i<=lenA;++i){	
			v[i][0]=0;
			f[i][0]=0;
			e[i][0]=0;
			g[i][0]=0;
		}
	
		for (j=1;j<=lenB;++j){
			v[0][j]=-j*gapOpenM;
			f[0][j]=-j*gapOpenM;
			e[0][j]=-j*gapOpenM;
			g[0][j]=-j*gapOpenM;
		}
	
	} else {
		for (i=0;i<=lenA;++i){	
			v[i][0]=0;
			f[i][0]=0;
			e[i][0]=0;
			g[i][0]=0;
		}
		
		for (j=1;j<=lenB;++j){
			v[0][j]=0;
			f[0][j]=0;
			g[0][j]=0;
			e[0][j]=0;
		}		
	}
	
		
	
	gapa = 0;
	gapb = 0;
	
	
	
	i=1;
	while(ptrarrayA[0]!=NULL){
		j=1;
			while(ptrarrayB[0]!=NULL){
			scorePP(ptrarrayA,ptrarrayB,p,scorematrix,spA,spB);
			g[i][j] = v[i-1][j-1] + scorematrix[i][j];
			e[i][j] = FMAX(e[i][j-1],v[i][j-1]-(p.nA-arrayI[i][5])*gapOpenM) - gapExtM * p.nA;
			f[i][j] = FMAX(f[i-1][j],v[i-1][j]-(p.nB-arrayJ[j][5])*gapOpenM) - gapExtM * p.nB;
			v[i][j] = FMAX(FMAX(g[i][j],e[i][j]), f[i][j] );
			maxxV = (maxxV > v[i][j] ? maxxV : v[i][j]);
			
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->previous_ptr;
			
			}
			++j;
		}
		
		
		for (k = 0; k < p.nA;++k ){
			ptrarrayA[k] = ptrarrayA[k]->previous_ptr;
		}
	
		
		for (k = 0; k < p.nB;++k ){
			ptrarrayB[k] = lastcpyptr[p.seqprofB[k]-1];
		}
		++i;
		 
	}

	
	
#ifdef DEBUG	
	if (!verbose) {
		printf ("maxV = %f\n",v[lenA][lenB]);
	}
#endif

	backtrack = (char*)calloc(1024,sizeof(char));
	
	for (k = 0; k < p.nA;++k ){
		ptrarrayA[k] = firstcpyptr[p.seqprofA[k]-1];
	}
	
	for (k = 0; k < p.nB;++k ){
		ptrarrayB[k] = firstcpyptr[p.seqprofB[k]-1];
	}
	
	i = lenA;
	j = lenB;
	
	
	maxRow=0;
	maxCol=0;
	for (k = i ; k >0 ; --k){
		maxRow = (maxRow < v[k][j] ? v[k][j] : maxRow);
	}
		
	for (k = j; k > 0; --k)	{
		maxCol = (maxCol < v[i][k] ? v[i][k] : maxCol);	
	}
	
	/*
	printf("maxRow: %f, maxCol %f\n",maxRow, maxCol);
	*/
	
	if (maxRow || maxCol > v[i-1][j-1]){
		if (maxRow > maxCol) {
			while (maxRow > v[i][j]) {
			/*fa*/	
				
			
			/* insert a gap in profile B */
			for (k = 0; k < p.nB; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = ptrarrayB[k];
				new_item_ptr->previous_ptr = ptrarrayB[k]->previous_ptr;
				new_item_ptr->seqID=ptrarrayB[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				if (ptrarrayB[k]->previous_ptr ==NULL){
					firstcpyptr[p.seqprofB[k]-1] = new_item_ptr;
					ptrarrayB[k]->previous_ptr = new_item_ptr;
				} else {
					ptrarrayB[k]->previous_ptr->next_ptr = new_item_ptr;
					ptrarrayB[k]->previous_ptr = new_item_ptr;

				}
		
			}
			strncat(backtrack ,"-" , 1);
			--i;
			for (k = 0; k < p.nA; ++k){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
			
				
					
			}
		} else {
			while (maxCol > v[i][j]) {
			/*ea*/
				
			
						
			/* insert a gap in profile A */
			for (k = 0; k < p.nA; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = ptrarrayA[k];
				new_item_ptr->previous_ptr = ptrarrayA[k]->previous_ptr;
				new_item_ptr->seqID = ptrarrayA[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				if (ptrarrayA[k]->previous_ptr ==NULL){
					firstcpyptr[p.seqprofA[k]-1] = new_item_ptr;
					ptrarrayA[k]->previous_ptr = new_item_ptr;
				} else {
					ptrarrayA[k]->previous_ptr->next_ptr = new_item_ptr;
					ptrarrayA[k]->previous_ptr = new_item_ptr;
				}
			}
			strncat(backtrack ,"-" , 1);	
			--j;
			for (k = 0; k < p.nB; ++k){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
			
				
				
			}
		
		}
	} 
	
	
	
	while (i > 0 && j > 0 ){
		va = v[i][j];
		ga = g[i][j];
		ea = e[i][j];
		fa = f[i][j];

		
		if (va == ga){
			strncat(backtrack ,symbolsPPcpy(ptrarrayA,ptrarrayB,p.nA,p.nB) , 1);
			for (k = 0; k < p.nA;++k ){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
			--i;
			--j;
			continue;
		} else if (va==ea) {
			
			/* insert a gap in profile A */
			for (k = 0; k < p.nA; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = ptrarrayA[k];
				new_item_ptr->previous_ptr = ptrarrayA[k]->previous_ptr;
				new_item_ptr->seqID = ptrarrayA[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				if (ptrarrayA[k]->previous_ptr ==NULL){
					firstcpyptr[p.seqprofA[k]-1] = new_item_ptr;
					ptrarrayA[k]->previous_ptr = new_item_ptr;
				} else {
					ptrarrayA[k]->previous_ptr->next_ptr = new_item_ptr;
					ptrarrayA[k]->previous_ptr = new_item_ptr;
				}
			}
			
			strncat(backtrack ,"-" , 1);	
			--j;
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
			continue;
		} else if (va==fa){
			
			/* insert a gap in profile B */
			for (k = 0; k < p.nB; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = ptrarrayB[k];
				new_item_ptr->previous_ptr = ptrarrayB[k]->previous_ptr;
				new_item_ptr->seqID=ptrarrayB[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				if (ptrarrayB[k]->previous_ptr ==NULL){
					firstcpyptr[p.seqprofB[k]-1] = new_item_ptr;
					ptrarrayB[k]->previous_ptr = new_item_ptr;
				} else {
					ptrarrayB[k]->previous_ptr->next_ptr = new_item_ptr;
					ptrarrayB[k]->previous_ptr = new_item_ptr;

				}
		
			}
			
			strncat(backtrack ,"-" , 1);
			--i;
			for (k = 0; k < p.nA;++k ){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
			continue;
		} else {
			printf("Should not be reached\n");
			exit(1);
		}
	}

	if (j>0 && i==0 ) {
		for (k = 0; k < p.nA; ++k){
			ptrarrayA[k] = lastcpyptr[p.seqprofA[k]-1];
		}
		/* insert gap(s) on 5' in profile A*/
		for ( m = j ; m > 0 ; --m){
			for (k = 0; k < p.nA; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = NULL;
				new_item_ptr->previous_ptr = ptrarrayA[k];
				new_item_ptr->seqID = ptrarrayA[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;

				ptrarrayA[k]->next_ptr = new_item_ptr;
				ptrarrayA[k] = new_item_ptr;
				lastcpyptr[p.seqprofA[k]-1] = new_item_ptr;
			}
			strncat(backtrack,"-",1);
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
		}
	} else if ( i > 0 &&  j==0 ){
		
		for (k = 0; k < p.nB; ++k){
			ptrarrayB[k] = lastcpyptr[p.seqprofB[k]-1];
		
		}
		/* insert gap(s) on 5' in profile B*/
		for ( n = i ; n > 0 ; --n) { 
			for (k = 0; k < p.nB; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = NULL;
				new_item_ptr->previous_ptr = ptrarrayB[k];
				new_item_ptr->seqID=ptrarrayB[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				ptrarrayB[k]->next_ptr = new_item_ptr;
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
				lastcpyptr[p.seqprofB[k]-1] = new_item_ptr;
			}
			strncat(backtrack,"-",1);
			for (k = 0; k < p.nA;++k ){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
		}
	}
#ifdef DEBUG
	if (!verbose){
		printf("Aligning ");
		for (k = 0; k < p.nA; ++k){
			printf(" %d  ",p.seqprofA[k]);
		}
		printf("\n");
	}
#endif	
	
	for (k = 0; k < p.nA; ++k){
		print_data_content2(lastcpyptr[p.seqprofA[k]-1]);
#ifdef DEBUG
		if (!verbose){
			printf("\n");
		}
#endif		
	}	

#ifdef DEBUG
	if (!verbose){		
		for (k = IMAX(lenA+gapa,lenB+gapb); k >= 0;--k){
			if (backtrack[k]!='\0')
			printf("%c",backtrack[k]);
		}	 
		printf("\n");
	}
#endif
	
	for (k = 0; k < p.nB; ++k){
		print_data_content2(lastcpyptr[p.seqprofB[k]-1]);
#ifdef DEBUG
		if (!verbose){	
			printf("\n");
		}
#endif		
	}	

#ifdef DEBUG
	if (!verbose){	
		printf(" with ");
		for (k = 0; k < p.nB; ++k){
			printf(" %d  ",p.seqprofB[k]);
		}
		printf("\n");
	
		printf("\nScore is:\n");
		printf("maxx V: %f\n", maxxV);
		printf("v[%d][%d]: %f\n",lenA,lenB, v[lenA][lenB]);
	}
#endif		
	
	
	for (k = 0; k <= lenA; ++k){
		free(arrayI[k]);
		free(scorematrix[k]);
	}
	for (k = 0; k <= lenB; ++k){
		free(arrayJ[k]);
	}
	
	
	free(ptrarrayA);
	free(ptrarrayB);
	free(backtrack);
	free(arrayI);	
	free(arrayJ);
	free(scorematrix);
	return(maxxV);
}


/* Function:		alignpieces
 * 
 * Purpose:		do the multiple sequence alignment with anchored sequences
 *           
 * Args:			p: contains the groups to be aligned
 *				lenA, lenB: length of the sequences in group I, J
 *				nseqs: number of sequences in alignment
 *
 * Returns:		float 
 */
float alignpieces(struct profile p, int part, int splits, struct vector ***segptr){
	int i = 0;
	int j = 0;
	int k = 0;
	int m = 0;
	int n = 0;
	int gapa=0;
	int gapb=0;
	
	struct vector **ptrarrayA;
	struct vector **ptrarrayB;
	float **v;    
	float **f;
	float **e;
	float **g;
	float **scorematrix;
	float va,ga,ea,fa,maxxV=0;
	char *backtrack;
	float sumWa =0;
	float sumWb =0;
	struct lscore *spA;
	struct lscore *spB;
	int lenA,lenB,beginA,beginB;
	
	/* fill arrays with pointers to vectors */
	ptrarrayA = (struct vector **)malloc(p.nA*sizeof(struct vector *));
	ptrarrayB = (struct vector **)malloc(p.nB*sizeof(struct vector *));
	
	
	
	/* get lenght of the parts to be aligned */

	for (k = 0; k < p.nA;++k ){
		ptrarrayA[k] = segptr[part][p.seqprofA[k]-1];
	}
	for (k = 0; k < p.nB;++k ){
		ptrarrayB[k] = segptr[part][p.seqprofB[k]-1];
	}

	
	
	printf("splits %d , part %d\n",splits, part);
	fflush(stdout);
	if (part != splits && part !=0) {
		fflush(stdout);
		lenA = segptr[part+1][p.seqprofA[0]-1]->pos-segptr[part][p.seqprofA[0]-1]->pos-1;
		lenB = segptr[part+1][p.seqprofB[0]-1]->pos-segptr[part][p.seqprofB[0]-1]->pos-1; 
		beginA = segptr[part][p.seqprofA[0]-1]->pos+1;
		beginB = segptr[part][p.seqprofB[0]-1]->pos+1;
		printf("\nbeginA =  %d, beginB  =   %d\n", beginA,beginB);
	} else if (part != splits || part ==0){
		fflush(stdout);
		lenA = segptr[part+1][p.seqprofA[0]-1]->pos-segptr[part][p.seqprofA[0]-1]->pos;
		lenB = segptr[part+1][p.seqprofB[0]-1]->pos-segptr[part][p.seqprofB[0]-1]->pos; 
		beginA = segptr[part][p.seqprofA[0]-1]->pos;
		beginB = segptr[part][p.seqprofB[0]-1]->pos;
		printf("\nbeginA =  %d, beginB  =   %d\n", beginA,beginB);
	} else {
		fflush(stdout);
		lenA = firstcpyptr[p.seqprofA[0]-1]->pos-segptr[part][p.seqprofA[0]-1]->pos;
		lenB = firstcpyptr[p.seqprofB[0]-1]->pos-segptr[part][p.seqprofB[0]-1]->pos; 
		beginA = segptr[part][p.seqprofA[0]-1]->pos+1;
		beginB = segptr[part][p.seqprofB[0]-1]->pos+1;
		printf("\nbeginA =  %d, beginB  =   %d\n", beginA,beginB);
	}
	
	
		
#ifdef DEBUG
	if (!verbose){
		printf("\nlenA =  %d, lenB  =   %d\n", lenA,lenB);
	}
#endif
	
	/* don't use precomputed scores in iteration */
	ct = 0;
	

		
	/* matrix storing the scores computed by scorePPcpy*/
	scorematrix=(float **)malloc((lenA+1)*sizeof(float *));
	
	/* arrays holding sequence information(nucleotide, gap frequencies) */
	arrayI = (float **)calloc(lenA+1,sizeof(float *));
	arrayJ = (float **)calloc(lenB+1,sizeof(float *));

	/* allocate matrices for the dynamic programming part */
	v = (float **)malloc((lenA+1)*sizeof(float*));
	e = (float **)malloc((lenA+1)*sizeof(float*));
	f = (float **)malloc((lenA+1)*sizeof(float*));
	g = (float **)malloc((lenA+1)*sizeof(float*));

	/* in 2 dimensions */
	for (k = 0; k <= lenA; ++k){
		arrayI[k] = (float *)calloc(9,sizeof(float));
		scorematrix[k] = (float *)malloc((lenB+1)*sizeof(float));
		v[k] = (float *)malloc((lenB+1)*sizeof(float));
		e[k] = (float *)malloc((lenB+1)*sizeof(float));
		f[k] = (float *)malloc((lenB+1)*sizeof(float));
		g[k] = (float *)malloc((lenB+1)*sizeof(float));
	}

	for (k = 0; k <= lenB; ++k){
		arrayJ[k] = (float *)calloc(9,sizeof(float));
	}
	
	
	for (k = 0; k < p.nA;++k ){
		ptrarrayA[k] = segptr[part][p.seqprofA[k]-1];
		data_content_cpy(ptrarrayA[k],lenA,arrayI);
		sumWa += sdata[ptrarrayA[k]->seqID].weight; 
	}

	for (k = 0; k < p.nB;++k ){
		ptrarrayB[k] = segptr[part][p.seqprofB[k]-1];
		data_content_cpy(ptrarrayB[k],lenB,arrayJ);
		sumWb += sdata[ptrarrayB[k]->seqID].weight; 
	}
	

	
	
	spA=NULL;
	spB=NULL;
	
	
	for (i=0;i<=lenA;++i){	
		v[i][0]=0;
		f[i][0]=0;
		e[i][0]=0;
		g[i][0]=0;
	}
		
	for (j=1;j<=lenB;++j){
		v[0][j]=0;
		f[0][j]=0;
		g[0][j]=0;
		e[0][j]=0;
	}
		
	
	gapa = 0;
	gapb = 0;
	
	
	
	for (k = 0; k < p.nA;++k ){
		if (part==0){
			ptrarrayA[k] = segptr[part][p.seqprofA[k]-1];
		} else {
			ptrarrayA[k] = segptr[part][p.seqprofA[k]-1]->previous_ptr;
		}
		/*fprintf(stdout,"ptrarrayA[k]->base: %c\n",ptrarrayA[k]->base);*/
	}

	for (k = 0; k < p.nB;++k ){
		if (part == 0){
			ptrarrayB[k] = segptr[part][p.seqprofB[k]-1];	
		} else {
			ptrarrayB[k] = segptr[part][p.seqprofB[k]-1]->previous_ptr;
		}
		/*fprintf(stdout,"ptrarrayB[k]->base: %c\n",ptrarrayB[k]->base);*/
	}
	
	

	
	/*i=1; */
	for (i=1; i <= lenA; i++){
		/*	while(ptrarrayA[0]!=NULL){*/
		/*j=1;*/
		for (j=1; j <= lenB; j++){
			/*printf("v[%2d][%2d]=%8.2f  ",i,j,v[i][j]);*/
			/*while(ptrarrayB[0]!=NULL){*/
			scorePPcpy(ptrarrayA,ptrarrayB,p,scorematrix,spA,spB,i,j);
			g[i][j] = v[i-1][j-1] + scorematrix[i][j];
			/*printf("scorematrix[%2d][%2d]=%8.2f  ",i,j,scorematrix[i][j]);*/
			e[i][j] = FMAX(e[i][j-1],v[i][j-1]-(p.nA-arrayI[i][5])*gapOpenM) - gapExtM * p.nA;
			f[i][j] = FMAX(f[i-1][j],v[i-1][j]-(p.nB-arrayJ[j][5])*gapOpenM) - gapExtM * p.nB;
			v[i][j] = FMAX(FMAX(g[i][j],e[i][j]), f[i][j] );
			maxxV = (maxxV > v[i][j] ? maxxV : v[i][j]);
			/*printf("v[%2d][%2d]=%8.2f    ",i,j,v[i][j]);*/
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->previous_ptr;
			
			}

		}
		/*printf("\n");*/
		
		for (k = 0; k < p.nA;++k ){
			ptrarrayA[k] = ptrarrayA[k]->previous_ptr;
		}
	
		for (k = 0; k < p.nB;++k ){
			if (part == 0){
				ptrarrayB[k] = segptr[part][p.seqprofB[k]-1];	
			} else {
				ptrarrayB[k] = segptr[part][p.seqprofB[k]-1]->previous_ptr;
			}
			/*fprintf(stdout,"ptrarrayB[k]->base: %c\n",ptrarrayB[k]->base);*/
		}
	}
/*
	printf("mat g\n");
	for (i=0; i <= lenA; i++){
		for (j=0; j <= lenB; j++){
			fprintf(stdout,"%8.2f ",g[i][j]);
		}
		printf("\n");
	}

	printf("mat e\n");
	for (i=0; i <= lenA; i++){
		for (j=0; j <= lenB; j++){
			fprintf(stdout,"%8.2f ",e[i][j]);
		}
		printf("\n");
	}

	printf("mat f\n");
	for (i=0; i <= lenA; i++){
		for (j=0; j <= lenB; j++){
			fprintf(stdout,"%8.2f ",f[i][j]);
		}
		printf("\n");
	}
	
	printf("mat v\n");
	for (i=0; i <= lenA; i++){
		for (j=0; j <= lenB; j++){
			fprintf(stdout,"%8.2f ",v[i][j]);
		}
		printf("\n");
	}
	
*/	
#ifdef DEBUG	
	if (!verbose) {
		printf ("maxV = %f\n",v[lenA][lenB]);
	}
#endif

	backtrack = (char*)calloc(1024,sizeof(char));
	
	if (part < splits){
		for (k = 0; k < p.nA;++k ){
			ptrarrayA[k] = segptr[part+1][p.seqprofA[k]-1]->next_ptr;
			/*fprintf(stdout,"X ptrarrayA[k]->base: %c\n",ptrarrayA[k]->base);	*/
		}	
	} else {
		for (k = 0; k < p.nA;++k ){
			ptrarrayA[k] = firstcpyptr[p.seqprofA[k]-1];
			/*fprintf(stdout,"X ptrarrayA[k]->base: %c\n",ptrarrayA[k]->base);	*/
		}
	}
		
	if (part < splits){
		for (k = 0; k < p.nB;++k ){
			ptrarrayB[k] = segptr[part+1][p.seqprofB[k]-1]->next_ptr;
			/*fprintf(stdout,"X ptrarrayB[k]->base: %c\n",ptrarrayB[k]->base);*/
		}
	} else {
		for (k = 0; k < p.nB;++k ){
			ptrarrayB[k] = firstcpyptr[p.seqprofB[k]-1];
			/*fprintf(stdout,"X ptrarrayB[k]->base: %c\n",ptrarrayB[k]->base);	*/
		}
	}
	
	i = lenA;
	j = lenB;
	while (i > 0 && j > 0 ){
		va = v[i][j];
		ga = g[i][j];
		ea = e[i][j];
		fa = f[i][j];

		printf("v[%d][%d] %f\n",i,j,v[i-1][j-1]);
		if (va == ga){
			strncat(backtrack ,symbolsPPcpy(ptrarrayA,ptrarrayB,p.nA,p.nB) , 1);
			for (k = 0; k < p.nA;++k ){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
			--i;
			--j;
			continue;
		} else if (va==ea) {
			
			/* insert a gap in profile A */
			for (k = 0; k < p.nA; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = ptrarrayA[k];
				new_item_ptr->previous_ptr = ptrarrayA[k]->previous_ptr;
				new_item_ptr->seqID = ptrarrayA[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				if (ptrarrayA[k]->previous_ptr ==NULL){
					firstcpyptr[p.seqprofA[k]-1] = new_item_ptr;
					ptrarrayA[k]->previous_ptr = new_item_ptr;
				} else {
					ptrarrayA[k]->previous_ptr->next_ptr = new_item_ptr;
					ptrarrayA[k]->previous_ptr = new_item_ptr;
				}
			}
			
			strncat(backtrack ,"-" , 1);	
			--j;
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
			continue;
		} else if (va==fa){
			
			/* insert a gap in profile B */
			for (k = 0; k < p.nB; ++k){
				struct vector *new_item_ptr;
				printf("working on 5' end\n");
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = ptrarrayB[k];
				new_item_ptr->previous_ptr = ptrarrayB[k]->previous_ptr;
				new_item_ptr->seqID=ptrarrayB[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				if (ptrarrayB[k]->previous_ptr ==NULL){
					firstcpyptr[p.seqprofB[k]-1] = new_item_ptr;
					ptrarrayB[k]->previous_ptr = new_item_ptr;
				} else {
					ptrarrayB[k]->previous_ptr->next_ptr = new_item_ptr;
					ptrarrayB[k]->previous_ptr = new_item_ptr;

				}
		
			}
			
			strncat(backtrack ,"-" , 1);
			--i;
			for (k = 0; k < p.nA;++k ){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
			continue;
		} else {
			printf("i %d , j %d",i,j);
			printf("Should not be reached\n");
			exit(1);
		}
	}

	if (j>0 && i==0 ) {
		for (k = 0; k < p.nA; ++k){
			ptrarrayA[k] = lastcpyptr[p.seqprofA[k]-1];
		}
		/* insert gap(s) on 5' in profile A*/
		for ( m = j ; m > 0 ; --m){
			for (k = 0; k < p.nA; ++k){
				struct vector *new_item_ptr;
				printf("working on 5' end\n");				
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = NULL;
				new_item_ptr->previous_ptr = ptrarrayA[k];
				new_item_ptr->seqID = ptrarrayA[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;

				ptrarrayA[k]->next_ptr = new_item_ptr;
				ptrarrayA[k] = new_item_ptr;
				lastcpyptr[p.seqprofA[k]-1] = new_item_ptr;
			}
			strncat(backtrack,"-",1);
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
		}
	} else if ( i > 0 &&  j==0 ){
		
		for (k = 0; k < p.nB; ++k){
			ptrarrayB[k] = lastcpyptr[p.seqprofB[k]-1];
		
		}
		/* insert gap(s) on 5' in profile B*/
		for ( n = i ; n > 0 ; --n) { 
			for (k = 0; k < p.nB; ++k){
				struct vector *new_item_ptr;
				new_item_ptr = malloc(sizeof (struct vector));
				new_item_ptr->next_ptr = NULL;
				new_item_ptr->previous_ptr = ptrarrayB[k];
				new_item_ptr->seqID=ptrarrayB[k]->seqID;
				new_item_ptr->base='-';
				new_item_ptr->oribase='-';
				new_item_ptr->replbase=5;
				new_item_ptr->p0=0;
				new_item_ptr->p1=0;
				new_item_ptr->p2=0;
				ptrarrayB[k]->next_ptr = new_item_ptr;
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
				lastcpyptr[p.seqprofB[k]-1] = new_item_ptr;
			}
			strncat(backtrack,"-",1);
			for (k = 0; k < p.nA;++k ){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
		}
	}
#ifdef DEBUG
	if (!verbose){
		printf("Aligning ");
		for (k = 0; k < p.nA; ++k){
			printf(" %d  ",p.seqprofA[k]);
		}
		printf("\n");
	}
#endif	
	
	for (k = 0; k < p.nA; ++k){
		print_data_content2(lastcpyptr[p.seqprofA[k]-1]);
#ifdef DEBUG
		if (!verbose){
			printf("\n");
		}
#endif		
	}	

#ifdef DEBUG
	if (!verbose){		
		for (k = 1; k < beginA; k++){
			printf(" ");
		}
		for (k = IMAX(lenA+gapa,lenB+gapb); k >= 0;--k){
			if (backtrack[k]!='\0')
			printf("%c",backtrack[k]);
		}	 
		printf("\n");
	}
#endif
	
	for (k = 0; k < p.nB; ++k){
		print_data_content2(lastcpyptr[p.seqprofB[k]-1]);
#ifdef DEBUG
		if (!verbose){	
			printf("\n");
		}
#endif		
	}	

#ifdef DEBUG
	if (!verbose){	
		printf(" with ");
		for (k = 0; k < p.nB; ++k){
			printf(" %d  ",p.seqprofB[k]);
		}
		printf("\n");
	
		printf("\nScore is:\n");
		printf("maxx V: %f\n", maxxV);
		printf("v[%d][%d]: %f\n",lenA,lenB, v[lenA][lenB]);
	}
#endif		
	
	
	for (k = 0; k <= lenA; ++k){
		free(arrayI[k]);
		free(scorematrix[k]);
	}
	for (k = 0; k <= lenB; ++k){
		free(arrayJ[k]);
	}
	
	
	free(ptrarrayA);
	free(ptrarrayB);
	free(backtrack);
	free(arrayI);	
	free(arrayJ);
	free(scorematrix);
	return(maxxV);
}



float **retree(int n){
	int i, j;
	float **remat;
	float scorev;
	struct vector *i_ptr, *j_ptr;
	remat = (float **)calloc(n,sizeof(float *));

	
	for (i = 0; i < n; i++){
		remat[i] = (float *)calloc(n,sizeof(float));		
		for (j = i+1; j < n; j++){
			i_ptr = last_ptr[i];
			j_ptr = last_ptr[j];
			scorev=0;
			fprintf(stdout,"%d %d\n",i,j);
			while (i_ptr!=NULL){
				scorev += score(i_ptr,j_ptr); 
				i_ptr=i_ptr->previous_ptr;
				j_ptr=j_ptr->previous_ptr;
			}
			fprintf(stdout,"score %f\n",scorev);
			remat[i][j]=scorev;
		}
	}	


	for (i = 0; i < n-1; i++){
		remat[i][i]=0;
		for (j = i+1; j < n; j++){
			remat[j][i]=remat[i][j];		
		}
	}

	fprintf(stdout,"\n ");
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			fprintf(stdout,"%8.2f ",remat[i][j]);
		}
		fprintf(stdout,"\n ");
	}
	
	return remat;
}












/* Function:		data_content
 * 
 * Purpose:		count occurence of each symbol in group X
 *           
 * Args:			current_ptr: points to first residue start from 
 *				arrayX: contains information about counts of symbols and sums up probabilities in a group
 *
 * Returns:		float arrayX
 */
float **data_content_cpy(struct vector *current_ptr, int len, float **arrayX){
	int i = 0;
	for (i = 0; i <=len; ++i) {
		if (current_ptr->base=='-')
			arrayX[i][5] +=sdata[current_ptr->seqID].weight;
		else if (current_ptr->base=='A')
			arrayX[i][0] +=1;
		else if (current_ptr->base=='C')
			arrayX[i][1] +=1;
		else if (current_ptr->base=='G')
			arrayX[i][2] +=1;
		else if (current_ptr->base=='U')
			arrayX[i][3] +=1;
		else if (current_ptr->base=='N')
			arrayX[i][4] +=1;
		arrayX[i][6] +=current_ptr->p0;
		arrayX[i][7] +=current_ptr->p1;
		arrayX[i][8] +=current_ptr->p2;
		current_ptr = current_ptr->previous_ptr;
	}
	return arrayX;
}













