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
#include "create_matrix.h"
#include "stral.h"
#include "nrutils.h"
#include "iofuncs.h"
#include "fft.h"



/* Function:		align
 * 
 * Purpose:		do the paiwise alignment
 *           
 * Args:			nseqs: number of sequences in alignment
 *				seqfile: name of input file
 *
 * Returns:		void
 */
void align(int nseqs, char *seqfile){
	int i, k;
	mode_t mode = S_IRWXU;
	array =(float **)malloc(nseqs*sizeof(float *));
	distmatrix =(float **)malloc(nseqs*sizeof(float *));
	zscore =(float **)malloc(nseqs*sizeof(float *));
	if (!nofile){
		chdir(workingdir);
		if (opendir("./pairwiseDIR")==NULL){
			mkdir("./pairwiseDIR",mode | S_IRGRP | S_IXGRP| S_IXOTH | S_IROTH);
		}
		chdir("./pairwiseDIR");
	}
	

	for ( i = 0; i < nseqs; ++i){
		*(array+i) = (float *)malloc(nseqs* sizeof(float)) ;
		*(zscore+i) = (float *)malloc(nseqs* sizeof(float)) ;
		*(distmatrix+i) = (float *)malloc(nseqs* sizeof(float)) ;
		for ( k=i+1; k < nseqs; ++k){
			
			
			if (fft){
				fftalign(last_ptr[i], last_ptr[k]);
			}
			
			
			method(first_ptr[i], last_ptr[i], first_ptr[k], last_ptr[k],nseqs, seqfile);

			
#ifdef DEBUG
			if (!verbose){
				printf("Aligning seq %d with %d\n", i+1,k+1);
				printf("Similarity %f\n", array[i][k]);
			}
#endif			
		}
	}
	if (spread){
		calczscore(nseqs);
	}
}


/* Function:		malign
 * 
 * Purpose:		do the multiple alignment
 *           
 * Args:			n: number of sequences in alignment
 *
 * Returns:		void
 */
void malign(int n){
	int i = 0;
	int j = 0;
	int k = 0;
	
	int sizeA = 0;
	int sizeB = 0;
	struct vector **ptrarray;
	
	
	ct=1;
	
	

	/* allocate memory for scoring vectors */
	fps = (struct lscore **)malloc(n * sizeof(struct lscore *));
	lps = (struct lscore **)malloc(n * sizeof(struct lscore *));	
	ptrarray = (struct vector **)malloc(n*sizeof(struct vector *));
	
	
	for (k = 0; k < n ; ++k){
		fps[k] = (struct lscore *)malloc(sizeof(struct lscore));
		lps[k] = (struct lscore *)malloc(sizeof(struct lscore));
	}
	
	
	
	for (i=n-2; i>=0;--i){
		
		for(j=i+1; j <= n-2 ;++j){
			if (prof[i].typeA==0 && prof[j].seqprofA[0] == prof[i].seqprofA[0]){
				prof[i].typeA = 1;
			}
		}
		
		if (prof[i].typeA==1){
			for (j=i+1; j <= n-2 ;++j){
				if (prof[j].seqprofA[0]==prof[i].seqprofA[0]){
					break;
				}
			}
			sizeA = prof[j].nA;
			sizeB = prof[j].nB;
			free(prof[i].seqprofA);
			prof[i].seqprofA=(int *)calloc((sizeA+sizeB), sizeof(int));

			/* kopieren aus knoten in seqprofA aus Teil A*/
			for (k=0;k < prof[j].nA; ++k){
				
				prof[i].seqprofA[k]=prof[j].seqprofA[k];
			}
			/* kopieren aus knoten in seqprofA aus Teil B*/
			for (k=0;k < prof[j].nB; ++k){
				
				prof[i].seqprofA[prof[j].nA+k]=prof[j].seqprofB[k];
			}
			/*aktualisieren der Sequenzzahl in A */
			prof[i].pnoA=j;
			prof[i].nA = prof[j].nA + prof[j].nB;
			prof[i].typeA = 1; /* bei neighbor gibt es hier (manchmal) einen fehler!! */
		}
		
		for(j=i+1; j <= n-2 ;++j){
			if (prof[i].typeB==0 && prof[j].seqprofA[0] == prof[i].seqprofB[0]){
				prof[i].typeB = 1;
			}
		}
		
		if (prof[i].typeB==1){
			for (j=i+1; j <= n-2 ;++j){
				if (prof[j].seqprofA[0]==prof[i].seqprofB[0]){
					break;
				}
			}

			
			sizeA = prof[j].nA;
			sizeB = prof[j].nB;
			free(prof[i].seqprofB);
			
			prof[i].seqprofB=(int *)calloc((sizeA+sizeB) ,sizeof(int));
			
			
			/* kopieren aus knoten in seqprofA aus Teil A*/
			for (k=0;k < prof[j].nA; ++k){
				
				prof[i].seqprofB[k]=prof[j].seqprofA[k];
			}
			/* kopieren aus knoten in seqprofA aus Teil B*/
			for (k=0;k < prof[j].nB; ++k){
				
				prof[i].seqprofB[prof[j].nA+k]=prof[j].seqprofB[k];
			}
			/*aktualisieren der Sequenzzahl in A */
			prof[i].pnoB = j;
			prof[i].nB = prof[j].nA + prof[j].nB;
			prof[i].typeA = 1;
			
	
		}
#ifdef DEBUG		
		if (!verbose){
				printf("\nCycle %d: %s %d(%d) joins %s %d(%d)\n",i,(prof[i].typeA==0 ? "Sequence" : "Cluster" ), prof[i].seqprofA[0],prof[i].pnoA, (prof[i].typeB==0 ? "Sequence" : "Cluster" ),
				prof[i].seqprofB[0], prof[i].pnoB);
		}
#endif		
		if (fft){	
			fftmultalign(prof[i] , sdata[prof[i].seqprofA[0]-1].seqlength , sdata[prof[i].seqprofB[0]-1].seqlength, n);
		}
		methodmult(prof[i] , sdata[prof[i].seqprofA[0]-1].seqlength , sdata[prof[i].seqprofB[0]-1].seqlength,n,i);
	}
}







/* pairwise sequence alignment */


/* Function:		method
 * 
 * Purpose:		do the paiwise alignment
 *           
 * Args:			first1_ptr, first2_ptr: points to the last residue of the sequence(first element of list)
 *				last1_ptr, last2_ptr: points to the first residue of the sequence(last element of list) 
 *				nseqs: number of sequences in alignment
 *				seqfile: input file
 *
 * Returns:		void
 */
void method(struct vector *first1_ptr, struct vector *last1_ptr, struct vector *first2_ptr, struct vector *last2_ptr, int nseqs, char *seqfile){
	
	struct vector *i_ptr, *j_ptr;
	float v[first1_ptr->pos+1][first2_ptr->pos+1];    /* first_ptr->pos equivalent to seq length*/
	float f[first1_ptr->pos+1][first2_ptr->pos+1];
	float e[first1_ptr->pos+1][first2_ptr->pos+1];
	float g[first1_ptr->pos+1][first2_ptr->pos+1];
	float endV, va, fa, ea, ga, maxxV = 0;
	int k, l, num, m, n, i, j;

	
	float diff, old,maxRow,maxCol, sum=0;
	char del='-';
	char *backtrack, *seqB, *seqA; 

	if (noterm){
		for (k=0;k<=first1_ptr->pos;++k){	
			v[k][0]=0;
			f[k][0]=0;
			e[k][0]=-gapOpen-k*gapExt;
			g[k][0]=0;
		}
			
		for (l=0;l<=first2_ptr->pos;++l){
			v[0][l]=0;
			f[0][l]=-gapOpen-l*gapExt;
			g[0][l]=0;
			e[0][l]=0;
		}
	} else {
		for (k=0;k<=first1_ptr->pos;++k){	
			v[k][0]=-gapOpen-k*gapExt;
			f[k][0]=0;
			e[k][0]=-gapOpen-k*gapExt;
			g[k][0]=0;
		}
			
		for (l=0;l<=first2_ptr->pos;++l){
			v[0][l]=-gapOpen-l*gapExt;
			f[0][l]=-gapOpen-l*gapExt;
			g[0][l]=0;
			e[0][l]=0;
		}
	}
	
	i_ptr=last1_ptr;
	j_ptr=last2_ptr;
	
	k=1;
	while(i_ptr!=NULL){
	l=1;
		while(j_ptr!=NULL){
				
	
			g[k][l] = v[k-1][l-1] + score(i_ptr, j_ptr);
			e[k][l] = FMAX(e[k][l-1],v[k][l-1]-gapOpen) - gapExt;
			f[k][l] = FMAX(f[k-1][l],v[k-1][l]-gapOpen) - gapExt;
			v[k][l] = FMAX( FMAX(g[k][l],e[k][l]), f[k][l] );
			
			
			maxxV = FMAX(maxxV,v[k][l]);
			++l;
			j_ptr=j_ptr->previous_ptr;
		}
		++k;
		i_ptr=i_ptr->previous_ptr;
		j_ptr=last2_ptr;
				
	}
	
	/*
		
	
	printf("matrix v\n");
	for (k=0; k <= first1_ptr->pos; ++k){
		for (l=0; l <= first2_ptr->pos; ++l){
			printf("%13.5f ",v[k][l]);
		}
		printf("\n");
	}
	printf("matrix e\n");
	for (k=0; k <= first1_ptr->pos; ++k){
		for (l=0; l <= first2_ptr->pos; ++l){
			printf("%13.5f ",e[k][l]);
		}
		printf("\n");
	}
	printf("matrix f\n");
	for (k=0; k <= first1_ptr->pos; ++k){
		for (l=0; l <= first2_ptr->pos; ++l){
			printf("%13.5f ",f[k][l]);
		}
		printf("\n");
	}
	printf("matrix g\n");
	for (k=0; k <= first1_ptr->pos; ++k){
		for (l=0; l <= first2_ptr->pos; ++l){
			printf("%13.5f ",g[k][l]);
		}
		printf("\n");
	}
	*/
	
	
	
	endV = v[k-1][l-1];
	
	
	
	
	backtrack = (char*)calloc(1024,sizeof(char));
	seqB = (char*)calloc(1024,sizeof(char));
	seqA = (char*)calloc(1024,sizeof(char));

	
	k = (first1_ptr->pos);
	l = (first2_ptr->pos);
	i_ptr = first1_ptr;
	j_ptr = first2_ptr;
	old = v[k-1][l-1];
	/* search for maximum in last row / last column */
	
	
	
	maxRow=0;
	maxCol=0;
	if (noterm){
		for (i = k ; i >0 ; --i){
			maxRow = (maxRow < v[i][l] ? v[i][l] : maxRow);
		}
		
		for (j = l; j > 0; --j)	{
			maxCol = (maxCol < v[k][j] ? v[k][j] : maxCol);	
		}
	}
	/*
	printf("maxRow: %f, maxCol %f\n",maxRow, maxCol);
	*/
	
	if (maxRow  > v[k-1][l-1] || maxCol > v[k-1][l-1]){
		if (maxRow > maxCol) {
			while (maxRow > v[k][l]) {
			/*fa*/	
				
				
				seqA[strlen(seqA)]=i_ptr->base;
				seqA[strlen(seqA)+1]='\0';
				seqB[strlen(seqB)]=del;
				seqB[strlen(seqB)+1]='\0';
				backtrack[strlen(backtrack)]=del;
				backtrack[strlen(backtrack)+1]='\0';
				if (i_ptr==NULL)	{
					fprintf(stderr, "FIXME: i_ptr invalid (va==fa)\n");
					exit(1);
				}
				i_ptr = i_ptr->next_ptr;
				--k;	
			}
		} else {
			while (maxCol > v[k][l]) {
			/*ea*/
				
				seqA[strlen(seqA)]=del;
				seqA[strlen(seqA)+1]='\0';
				seqB[strlen(seqB)]=j_ptr->base;
				seqB[strlen(seqB)+1]='\0';
				backtrack[strlen(backtrack)]=del;
				backtrack[strlen(backtrack)+1]='\0';
				if (j_ptr==NULL)	{
					fprintf(stderr, "FIXME: j_ptr invalid (va==ea)\n");
					exit(1);
				}
				j_ptr = j_ptr->next_ptr;
				
				
				
				--l;
			}
		
		}
	} 
	
	
	
	/*
	printf("k: %d, l: %d\n", k ,l );
	printf("%s\n",seqA);
	printf("%s\n",seqB);*/
	
	while (k > 0 && l > 0 && i_ptr!=NULL && j_ptr !=NULL){ /* paranoia test */

		ga = g[k][l];
		ea = e[k][l];
		fa = f[k][l];
		va = v[k][l];
		if (va == ga){
			
		/*	printf("match, va = %13.5f, ea = %13.5f, fa = %13.5f, ga = %13.5f i: %d, j %d v[k][l]= %13.5f\n",va,ea,fa,ga,k,l,v[k][l]);*/
					
			strncat(seqA,&i_ptr->base,1);
			strncat(seqB,&j_ptr->base,1);
			strncat(backtrack ,symbols(i_ptr,j_ptr) , 1);
			
			if (i_ptr==NULL)	{
				fprintf(stderr, "FIXME: i_ptr invalid (va == ga)\n");
				exit(1);
			}
			if (j_ptr==NULL)	{
				fprintf(stderr, "FIXME: j_ptr invalid (va == ga)\n");
				exit(1);
			}
			i_ptr = i_ptr->next_ptr;
			j_ptr = j_ptr->next_ptr;
			
			--k;
			--l;
		} else if (va==ea) {
/*			printf("ea gap, va = %13.5f, ea = %13.5f, fa = %13.5f, ga = %13.5f i: %d, j %d v[k][l]= %13.5f,v[k][l-1]= %13.5f \n",va,ea,fa,ga,k,l,v[k][l],v[k][l-1]);*/
			
			
			seqA[strlen(seqA)]=del;
			seqA[strlen(seqA)+1]='\0';
			seqB[strlen(seqB)]=j_ptr->base;
			seqB[strlen(seqB)+1]='\0';
			backtrack[strlen(backtrack)]=del;
			backtrack[strlen(backtrack)+1]='\0';
			if (j_ptr==NULL)	{
				fprintf(stderr, "FIXME: j_ptr invalid (va==ea)\n");
				exit(1);
			}
			j_ptr = j_ptr->next_ptr;
			
			--l;
			/*printf("i: %d, j %d\n",k,l);*/
		} else if (va==fa){
			
			
/*			printf("fa gap, va = %13.5f, ea = %13.5f, fa = %13.5f, ga = %13.5f i: %d, j %d v[k][l]= %13.5f\n",va,ea,fa,ga,k,l,v[k][l]);*/
			
			seqA[strlen(seqA)]=i_ptr->base;
			seqA[strlen(seqA)+1]='\0';
			seqB[strlen(seqB)]=del;
			seqB[strlen(seqB)+1]='\0';
			backtrack[strlen(backtrack)]=del;
			backtrack[strlen(backtrack)+1]='\0';
			if (i_ptr==NULL)	{
				fprintf(stderr, "FIXME: i_ptr invalid (va==fa)\n");
				exit(1);
			}
			i_ptr = i_ptr->next_ptr;
			--k;
			
		} else {
			printf("Should not be reached\n");
			exit(1);
		}
		diff = old-v[k][l];
		sum = sum +diff;
		++num;
		old = v[k][l];
	}
	
	if (l>0 && k==0 ) {
	    for ( m = l ; m > 0 ; --m){
		 	
			seqA[strlen(seqA)]=del;
			seqA[strlen(seqA)+1]='\0';
			seqB[strlen(seqB)]=j_ptr->base;
			seqB[strlen(seqB)+1]='\0';
			backtrack[strlen(backtrack)]=del;
			backtrack[strlen(backtrack)+1]='\0';
			if (j_ptr==NULL)	{
				fprintf(stderr, "FIXME: j_ptr invalid (l>0 && k==0)\n");
				exit(1);
			}
			j_ptr=j_ptr->next_ptr;
			
		}
	} else if ( k > 0 &&  l==0 ){
	    for ( n = k ; n > 0 ; --n) {
			
			seqA[strlen(seqA)]=i_ptr->base;
			seqA[strlen(seqA)+1]='\0';
			seqB[strlen(seqB)]=del;
			seqB[strlen(seqB)+1]='\0';
			backtrack[strlen(backtrack)]=del;
			backtrack[strlen(backtrack)+1]='\0';
			if (i_ptr==NULL)	{
				fprintf(stderr, "FIXME: i_ptr invalid ( k > 0 &&  l==0)\n");
				exit(1);
			}
			i_ptr=i_ptr->next_ptr;
			
	    }
    }
	
	fprintf(stdout,"aln->length: %d\n",strlen(seqA));
	zscore[first1_ptr->seqID][first2_ptr->seqID] = v[first1_ptr->pos][first2_ptr->pos]/strlen(seqA);
	fprintf(stdout,"zscore: %f\n",zscore[first1_ptr->seqID][first2_ptr->seqID]);
	
#ifdef DEBUG	
	if (!verbose){
		printf("\n");
		for (k = strlen(seqA); k >= 0;--k){
			if (seqA[k]!='\0')
			printf("%c",seqA[k]);
		}
		printf("\n");
		for (k = strlen(backtrack); k >= 0;--k){
			if (backtrack[k]!='\0')
			printf("%c",backtrack[k]);
		}
		printf("\n");
		for (k = strlen(seqB); k >= 0;--k){
			if (seqB[k]!='\0')	
			printf("%c",seqB[k]);
		}
		printf("\n");
	}
#endif
		
	if (ppair && nseqs == 2){
		printpairwise(first1_ptr->seqID, first2_ptr->seqID, seqA, seqB,1,seqfile);
	} else if (ppair){
		printpairwise(first1_ptr->seqID, first2_ptr->seqID, seqA, seqB,0,"");
	}
	
	if (backtrack!=NULL){
		free(backtrack);
		backtrack = NULL;
	}
	if (seqA!=NULL){
		free(seqA);
		seqA = NULL;
	}
	if (seqB!=NULL){
		free(seqB);
		seqB = NULL;		
	}
		
	array[first1_ptr->seqID][first2_ptr->seqID]=maxxV;
}





/* Function:		methodmult
 * 
 * Purpose:		do the multiple sequence alignment
 *           
 * Args:			p: contains the groups to be aligned
 *				lenA, lenB: length of the sequences in group I, J
 *				nseqs: number of sequences in alignment
 *				cp: current profile number
 *
 * Returns:		float 
 */
float methodmult(struct profile p, int lenA, int lenB, int nseqs, int cp){
	int i,j,k,m,n;
	int numseq;
	struct vector **ptrarrayA;
	struct vector **ptrarrayB;
	struct lscore *lptr, *spA, *spB;
	float **v;    
	float **f;
	float **e;
	float **g;
	float **scorematrix;	
	float va,ga,ea,fa;
	char *backtrack;
	float sumWa =0;
	float sumWb =0;
	float maxRow, maxCol;

#ifdef DEBUG
	if (!verbose){
		printf("\nlenA =  %d, lenB  =   %d\n", lenA,lenB);
	}
#endif

	ptrarrayA = (struct vector **)malloc(p.nA*sizeof(struct vector *));
	ptrarrayB = (struct vector **)malloc(p.nB*sizeof(struct vector *));
	
	
	numseq = (p.nA+p.nB) * (p.nA+p.nB-1) / 2;
		
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

	spA=NULL;
	spB=NULL;
	
	
	for (k = 0; k < p.nA;++k ){
		ptrarrayA[k] = last_ptr[p.seqprofA[k]-1];
		data_content(last_ptr[p.seqprofA[k]-1],arrayI);
		sumWa += sdata[ptrarrayA[k]->seqID].weight; 
	}
	

	for (k = 0; k < p.nB;++k ){
		ptrarrayB[k] = last_ptr[p.seqprofB[k]-1];
		data_content(last_ptr[p.seqprofB[k]-1],arrayJ);
		sumWb += sdata[ptrarrayB[k]->seqID].weight;
	}



	if (noterm) {
		for (i=0;i<=lenA;++i){	
			v[i][0]=0;
			f[i][0]=0;
			e[i][0]=(-gapOpenM-i*gapExtM)*numseq;
			g[i][0]=0;
		}
	
		for (j=1;j<=lenB;++j){
			v[0][j]=0;
			f[0][j]=(-gapOpenM-j*gapExtM)*numseq;
			e[0][j]=0;
			g[0][j]=0;
		}
	
	} else {
		for (i=0;i<=lenA;++i){	
			v[i][0]=(-gapOpenM-i*gapExtM)*numseq;
			f[i][0]=0;
			e[i][0]=(-gapOpenM-i*gapExtM)*numseq;
			g[i][0]=0;
		}
	
		for (j=1;j<=lenB;++j){
			v[0][j]=(-gapOpenM-j*gapExtM)*numseq;
			f[0][j]=(-gapOpenM-j*gapExtM)*numseq;
			g[0][j]=0;
			e[0][j]=0;
		}		
	}
	
	
	i=1;
	if(ct) { 
		spA = fps[p.pnoA];
		spB = fps[p.pnoB];
	}
	
	while(ptrarrayA[0]!=NULL){
		j=1;
	
		if(ct) { 
			spB = fps[p.pnoB];
		}
		while(ptrarrayB[0]!=NULL){
			scorePP(ptrarrayA,ptrarrayB,p,scorematrix,spA,spB);
			g[i][j] = v[i-1][j-1] + scorematrix[i][j];
			e[i][j] = FMAX(e[i][j-1],v[i][j-1]-(p.nA-arrayI[i][5])*gapOpenM) - gapExtM * p.nA;
			f[i][j] = FMAX(f[i-1][j],v[i-1][j]-(p.nB-arrayJ[j][5])*gapOpenM) - gapExtM * p.nB;
			v[i][j] = FMAX(FMAX(g[i][j],e[i][j]),f[i][j]);
		
			
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->previous_ptr;
			}
			++j;
			if(ct){
				if (p.nB>1)
					spB=spB->next_ptr;
			}
		}
		
		
		for (k = 0; k < p.nA;++k ){
			ptrarrayA[k] = ptrarrayA[k]->previous_ptr;
		}
	
		
		for (k = 0; k < p.nB;++k ){
			ptrarrayB[k] = last_ptr[p.seqprofB[k]-1];
		}
		++i;
		if(ct) { 
			if (p.nA>1)
				spA=spA->next_ptr;
		}
	}
	 
	
	
	
#ifdef DEBUG	
	if (!verbose) {
		printf ("maxV = %f\n",v[lenA][lenB]);
/*	
	printf("matrix v\n");
	for (i=0; i <= lenA; ++i){
		for (j=0; j <= lenB; ++j){
			printf("%13.5f ",v[i][j]);
		}
		printf("\n");
	}
	

	printf("matrix e\n");
	for (i=0; i <= lenA; ++i){
		for (j=0; j <= lenB; ++j){
			printf("%13.5f ",e[i][j]);
		}
		printf("\n");
	}

	
	printf("matrix f\n");
	for (i=0; i <= lenA; ++i){
		for (j=0; j <= lenB; ++j){
			printf("%13.5f ",f[i][j]);
		}
		printf("\n");
	}


	printf("matrix g\n");
	for (i=0; i <= lenA; ++i){
		for (j=0; j <= lenB; ++j){
			printf("%13.5f ",g[i][j]);
		}
		printf("\n");
	}*/
}
	
#endif
	
	backtrack = (char*)calloc(1024,sizeof(char));
		
	for (k = 0; k < p.nA;++k ){
		ptrarrayA[k] = first_ptr[p.seqprofA[k]-1];
	}
	
	for (k = 0; k < p.nB;++k ){
		ptrarrayB[k] = first_ptr[p.seqprofB[k]-1];
	}
	i = lenA;
	j = lenB;

	
	lptr = malloc(sizeof(struct lscore));

	
	
	maxRow=0;
	maxCol=0;
	if (noterm){
	for (k = i ; k >0 ; --k){
		maxRow = (maxRow < v[k][j] ? v[k][j] : maxRow);
	}
		
	for (k = j; k > 0; --k)	{
		maxCol = (maxCol < v[i][k] ? v[i][k] : maxCol);	
	}
	
	/*
	printf("maxRow: %f, maxCol %f\n",maxRow, maxCol);
	*/
	}
	if (maxRow > v[i-1][j-1] || maxCol > v[i-1][j-1]){
		if (maxRow > maxCol) {
			while (maxRow > v[i][j]) {
			/*fa*/	
				if (ct){
				if (i==1 && j ==1){
		     		lptr = lptr->next_ptr;
		     		free(lptr->previous_ptr);
		     		fps[cp]->score = recalculatescorePP(ptrarrayA,p.nA);
		     		fps[cp]->previous_ptr = NULL;
		     		fps[cp]->next_ptr = lptr;
		     		lptr->previous_ptr = fps[cp];
				} else if (i < lenA && j < lenB){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayA,p.nA);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		     	} else if ((i == lenA && j < lenB) || (j == lenB && i < lenA)){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayA,p.nA);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		     	} else if (i == lenA && j == lenB){
		     		lps[cp]->score = recalculatescorePP(ptrarrayA,p.nA);
		     		lps[cp]->next_ptr = NULL;
		     		lps[cp]->previous_ptr = lptr;
		     		lptr->next_ptr = lps[cp];
		     	} else {
		     		printf("lenA: %d, lenB: %d\ni: %d, j: %d\n",lenA,lenB,i,j);
		     		printf("exception(va=fa)\nProgram terminates\n");
		     		exit(1);
		     	}
			}
			
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
					first_ptr[p.seqprofB[k]-1] = new_item_ptr;
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
				if (ct){
				if (i==1 && j ==1){
		     		lptr = lptr->next_ptr;
		     		free(lptr->previous_ptr);
					fps[cp]->score = recalculatescorePP(ptrarrayB,p.nB);
		     		fps[cp]->previous_ptr = NULL;
		     		fps[cp]->next_ptr = lptr;
		     		lptr->previous_ptr = fps[cp];
			} else if (i < lenA && j < lenB){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayB,p.nB);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if ((i == lenA && j < lenB) || (j == lenB && i < lenA)){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayB,p.nB);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if (i == lenA && j == lenB){
		     		lps[cp]->score = recalculatescorePP(ptrarrayB,p.nB);
		     		lps[cp]->next_ptr = NULL;
		     		lps[cp]->previous_ptr = lptr;
		     		lptr->next_ptr = lps[cp];
		     	} else {
		     		printf("lenA: %d, lenB: %d\ni: %d, j: %d\n",lenA,lenB,i,j);
		     		printf("exception(va=ea)\nProgram terminates\n");
		     		exit(1);
		     	}
		}
			
						
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
					first_ptr[p.seqprofA[k]-1] = new_item_ptr;
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
	} else {
		
	}
	



	
	while (i > 0 && j > 0 ){
		ga = g[i][j];
		ea = e[i][j];
		fa = f[i][j];
		va = v[i][j];

		if (va == ga){
			if (ct){
				if (i==1 && j ==1){
		     		lptr = lptr->next_ptr;
		     		free(lptr->previous_ptr);
		     		fps[cp]->score = scorematrix[i][j];
		     		fps[cp]->previous_ptr = NULL;
		     		fps[cp]->next_ptr = lptr;
		     		lptr->previous_ptr = fps[cp];
				} else if (i < lenA && j < lenB){
					struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score=scorematrix[i][j];
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if ((i == lenA && j < lenB) || (j == lenB && i < lenA)){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score=scorematrix[i][j];
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if (i == lenA && j == lenB){
		     		lps[cp]->score = scorematrix[i][j];
		     		lps[cp]->next_ptr = NULL;
		     		lps[cp]->previous_ptr = lptr;
		     		lptr->next_ptr = lps[cp];
		     	} else {
		     		printf("lenA: %d, lenB: %d\ni: %d, j: %d\n",lenA,lenB,i,j);
		     		printf("exception(va=ga)\nProgram terminates\n");
		     		exit(1);
		     	}
			}
			strncat(backtrack ,symbolsPP(ptrarrayA,ptrarrayB,p.nA,p.nB) , 1);
			for (k = 0; k < p.nA;++k ){
				ptrarrayA[k] = ptrarrayA[k]->next_ptr;
			}
			for (k = 0; k < p.nB;++k ){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
			--i;
			--j;
			continue;
		
		
		} else if (va == ea) {
			if (ct){
				if (i==1 && j ==1){
		     		lptr = lptr->next_ptr;
		     		free(lptr->previous_ptr);
					fps[cp]->score = recalculatescorePP(ptrarrayB,p.nB);
		     		fps[cp]->previous_ptr = NULL;
		     		fps[cp]->next_ptr = lptr;
		     		lptr->previous_ptr = fps[cp];
			} else if (i < lenA && j < lenB){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayB,p.nB);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if ((i == lenA && j < lenB) || (j == lenB && i < lenA)){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayB,p.nB);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if (i == lenA && j == lenB){
		     		lps[cp]->score = recalculatescorePP(ptrarrayB,p.nB);
		     		lps[cp]->next_ptr = NULL;
		     		lps[cp]->previous_ptr = lptr;
		     		lptr->next_ptr = lps[cp];
		     	} else {
		     		printf("lenA: %d, lenB: %d\ni: %d, j: %d\n",lenA,lenB,i,j);
		     		printf("exception(va=ea)\nProgram terminates\n");
		     		exit(1);
		     	}
		}
			
						
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
					first_ptr[p.seqprofA[k]-1] = new_item_ptr;
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
			continue;
		} else if (va == fa){
			if (ct){
				if (i==1 && j ==1){
		     		lptr = lptr->next_ptr;
		     		free(lptr->previous_ptr);
		     		fps[cp]->score = recalculatescorePP(ptrarrayA,p.nA);
		     		fps[cp]->previous_ptr = NULL;
		     		fps[cp]->next_ptr = lptr;
		     		lptr->previous_ptr = fps[cp];
			} else if (i < lenA && j < lenB){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayA,p.nA);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if ((i == lenA && j < lenB) || (j == lenB && i < lenA)){
		     		struct lscore *nip;
		     		nip = malloc(sizeof(struct lscore));
		     		lptr->score = recalculatescorePP(ptrarrayA,p.nA);
		     		nip->next_ptr=lptr;
		     		lptr->previous_ptr = nip;
		     		lptr=nip;
		
		     	} else if (i == lenA && j == lenB){
		     		lps[cp]->score = recalculatescorePP(ptrarrayA,p.nA);
		     		lps[cp]->next_ptr = NULL;
		     		lps[cp]->previous_ptr = lptr;
		     		lptr->next_ptr = lps[cp];
		     	} else {
		     		printf("lenA: %d, lenB: %d\ni: %d, j: %d\n",lenA,lenB,i,j);
		     		printf("exception(va=fa)\nProgram terminates\n");
		     		exit(1);
		     	}
			}
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
					first_ptr[p.seqprofB[k]-1] = new_item_ptr;
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
			continue;
		} else {
			printf("Should not be reached\n");
			exit(1);
		}
	}

	if (j>0 && i==0 ) {
		for (k = 0; k < p.nA; ++k){
			ptrarrayA[k] = last_ptr[p.seqprofA[k]-1];
		}
		
		/* insert gap(s) on 5' in profile A*/
		for ( m = j ; m > 0 ; --m){
			if (ct){
				if (m == 1){
					lptr = lptr->next_ptr;
		     		free(lptr->previous_ptr);
					fps[cp]->score = recalculatescorePP(ptrarrayB,p.nB);
					fps[cp]->previous_ptr = NULL;
					fps[cp]->next_ptr = lptr;
					lptr->previous_ptr = fps[cp];
				} else if (m < lenB){
					struct lscore *nip;
					nip = malloc(sizeof(struct lscore));
					lptr->score = recalculatescorePP(ptrarrayB,p.nB);
					nip->next_ptr=lptr;
					lptr->previous_ptr = nip;
					lptr=nip;
				} else {
					printf("exception(A)\nProgram terminates\n");
					exit(1);
				}
			}		
		
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
				last_ptr[p.seqprofA[k]-1] = new_item_ptr;
			}
			strncat(backtrack,"-",1);
			for (k = 0; k < p.nB; ++k){
				ptrarrayB[k] = ptrarrayB[k]->next_ptr;
			}
		}
	
	} else if ( i > 0 &&  j==0 ){
		for (k = 0; k < p.nB; ++k){
			ptrarrayB[k] = last_ptr[p.seqprofB[k]-1];
		}
		
				
		
		/* insert gap(s) on 5' in profile B*/
		for ( n = i ; n > 0 ; --n) { 
			if(ct){
				if (n == 1){
					lptr = lptr->next_ptr;
		     		free(lptr->previous_ptr);					
					fps[cp]->score = recalculatescorePP(ptrarrayA,p.nA);
					fps[cp]->previous_ptr = NULL;
					fps[cp]->next_ptr = lptr;
					lptr->previous_ptr = fps[cp];
				} else if (n < lenA){
					struct lscore *nip;
					nip = malloc(sizeof(struct lscore));
					lptr->score=recalculatescorePP(ptrarrayA,p.nA);
					nip->next_ptr=lptr;
					lptr->previous_ptr = nip;
					lptr=nip;
				} else {
					printf("exception(B)\nProgram terminates\n");
					exit(1);
				}
			}
			
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
				last_ptr[p.seqprofB[k]-1] = new_item_ptr;
			}
			strncat(backtrack,"-",1);
			for (k = 0; k < p.nA; ++k){
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
		print_data_content2(last_ptr[p.seqprofA[k]-1]);
#ifdef DEBUG			
		if (!verbose){
			printf("\n");
		}
#endif		
	}	
	
#ifdef DEBUG	
	if (!verbose){		
		for (k = IMAX(lenA,lenB); k >= 0;--k){
			if (backtrack[k]!='\0')
			printf("%c",backtrack[k]);
		}	 
		printf("\n");
	}
#endif
	
	for (k = 0; k < p.nB; ++k){
		print_data_content2(last_ptr[p.seqprofB[k]-1]);
#ifdef DEBUG	
		if (!verbose){	
			printf("\n");
		}
#endif
	}	


/*	levelpos(cp);*/
	

#ifdef DEBUG	
	if (!verbose){	
		printf(" with ");
		for (k = 0; k < p.nB; ++k){
			printf(" %d  ",p.seqprofB[k]);
		}
		printf("\n");
	}
#endif	
	if(!nofile){
		if( (p.nA + p.nB) == nseqs){
			msascore= v[lenA][lenB];
			printf("\nScore is:\n");
			printf("%f\n", v[lenA][lenB]);
		}
	}
	
	
	
	
	
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
	
	return(v[lenA][lenB]);
}



void levelpos(int cluster){
	struct lscore *lptr;
	int pstn = 1;
	lptr = fps[cluster];
	
	while(lptr!=NULL){
		lptr->pos = pstn;
		lptr = lptr->next_ptr;
		pstn++;	
	}
}




/* Function:		print_data_content2
 * 
 * Purpose:		printout aligned sequence and update residue positions in alignment
 *           
 * Args:			current_ptr: pointer to starting residue
 *
 * Returns:		void
 */
void print_data_content2(struct vector *current_ptr){
	int pstn=1;
	int seqIDTY;
	seqIDTY = current_ptr->seqID;
	while(current_ptr!=NULL) {
#ifdef DEBUG
		if (!verbose){	
			printf("%c", current_ptr->base);
		}
#endif	
	
		current_ptr->pos = pstn;
		current_ptr = current_ptr->previous_ptr;
		++pstn;
	}
	sdata[seqIDTY].seqlength = pstn-1;
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
float **data_content(struct vector *current_ptr, float **arrayX){
	struct vector *old_ptr;
	while(current_ptr!=NULL) {
		if (current_ptr->base=='-')
			arrayX[current_ptr->pos][5] +=sdata[current_ptr->seqID].weight;
		else if (current_ptr->base=='A')
			arrayX[current_ptr->pos][0] +=1;
		else if (current_ptr->base=='C')
			arrayX[current_ptr->pos][1] +=1;
		else if (current_ptr->base=='G')
			arrayX[current_ptr->pos][2] +=1;
		else if (current_ptr->base=='U')
			arrayX[current_ptr->pos][3] +=1;
		else if (current_ptr->base=='N')
			arrayX[current_ptr->pos][4] +=1;
		arrayX[current_ptr->pos][6] +=current_ptr->p0;
		arrayX[current_ptr->pos][7] +=current_ptr->p1;
		arrayX[current_ptr->pos][8] +=current_ptr->p2;
		old_ptr= current_ptr;
		current_ptr = current_ptr->previous_ptr;

	}
	return arrayX;
}



void calczscore(int n){
	float meanscore;	/* Mittelwert */
	float zmin,zmax,ztmp;
	float *zi;
	int combine;	/* Anzahl der Stichproben */
	float sdsquare;	/* Standardabeweichung ^2 */
	float sdeviation;	/* Standardabweichung */
	/* zscore ist der beobachtete Wert der Zufallsvariablen */
	int i, j;
	int boolean  = 0;
	zmin = 10.e6;
	zmax = -10.e6;
	zi = (float *)malloc(n*sizeof(float));
	combine = (n * (n -1)) / 2;
	for (i = 0; i < n-1; i++){
		for (j = i + 1; j < n; j++){
			meanscore += zscore[i][j];
		}
	}
	meanscore /=combine;
	fprintf(stdout, "Mean Score: %f\n", meanscore);

	for (i = 0; i < n-1; i++){
		for (j = i + 1; j < n; j++){
			sdsquare += pow((zscore[i][j]-meanscore),2);
		}
	}
	sdeviation = sqrt(sdsquare/combine);
	fprintf(stdout, "Standard deviation^2: %f\n", sdsquare);
	fprintf(stdout, "Standard deviation: %f\n", sdeviation);
	
	for (i = 0; i < n-1; i++){
		for (j = i + 1; j < n; j++){
			ztmp = (zscore[i][j]-meanscore)/sdeviation;
			fprintf(stdout, "Z[%d][%d]: %f\n",i,j,ztmp);
			zmin = (zmin < ztmp ? zmin : ztmp);
			zmax = (zmax > ztmp ? zmax : ztmp);
		}
	}

	fprintf(stdout, "Zmax: %f\n", zmax);
	fprintf(stdout, "Zmin: %f\n", zmin);
	
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			if (i < j){
				zi[i] += (zscore[i][j]-meanscore)/sdeviation;
			} else if (j < i) {
				zi[i] += (zscore[j][i]-meanscore)/sdeviation;
			}
		}
		zi[i] = zi[i]/( n-1);
		/* ??? is a z-value of one a proper value ??? */
		if (zi[i] < -1 || zi	[i] > 1){
			fprintf(stdout, "Scores for sequence # %d (%s) deviate: Avg. z-value is %f\n", i+1,  sdata[i].ID, zi[i]);
			boolean = 1;
		}
	}
	if (boolean){
		fprintf(stdout, "Some sequences do not seem to fit to the others. Think of removing them. Calculation stops here!\n");
		exit(1);
	}
}
