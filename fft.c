#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrutils.h"
#include "stral.h"
#include "align.h"
#include "alfuncs.h"
#include "fft.h"

#define SWAP( a, b ) tempr=( a ); ( a ) = ( b ); ( b ) = tempr;
#define PI 3.14159265358979323846
#define min(a,b) ( (a) < (b) ? a : b)

struct segment **smatrix;
struct listPath **pathContainer;

int stepSize = 6;
int maxk,maxl;
int MAXI, MAXJ;

/* Function:		fftalign
 * 
 * Purpose:		calculates correlation between two sequences
 *           
 * Args:			last1_ptr, last2_ptr: pointer to the first sequence item
 *
 * Returns:		void
 */
void fftalign(struct vector *last1_ptr, struct vector *last2_ptr){
	float *data1;
	float *data2;
	float *ansp0, *ansp1, *ansp2, *ansA, *ansC,*ansG,*ansU, *anssump1, *anssump2;
	struct vector *item_ptr1,*item_ptr2, *old1_ptr, *old2_ptr;
	int i,j;
	int lengthSEQ1,lengthSEQ2;
	int winsize;
	int sizeA,sizeB;
	int elmax;
	float max;
	int *array;
	
		
	old1_ptr=last1_ptr;
	old2_ptr=last2_ptr;
	
	lengthSEQ1 = sdata[last1_ptr->seqID].seqlength;
	lengthSEQ2 = sdata[last2_ptr->seqID].seqlength;
	
	
	for (i=1; i < 12; ++i){
		if (lengthSEQ1 % (int)pow(2.0,(double)i)==lengthSEQ1){
			sizeA=(int)pow(2.0,(double)i);
			break;
		}
	}
	
	for (i=1; i < 12; ++i){
		if (lengthSEQ2 % (int)pow(2.0,(double)i)==lengthSEQ2){
			sizeB=(int)pow(2.0,(double)i);
			break;
		}
	}
	
	winsize=1*(sizeA > sizeB ? sizeA : sizeB);
	
	printf("SeqID-1:%s, length: %d, sizeA: %d\n",sdata[last1_ptr->seqID].ID, lengthSEQ1, sizeA);
	printf("SeqID-2:%s, length: %d, sizeB: %d\n",sdata[last2_ptr->seqID].ID, lengthSEQ2, sizeB);	
	printf("winsize: %d\n",winsize);
	
	
		
	data1=(float*)calloc(winsize,sizeof(float));
	data2=(float*)calloc(winsize,sizeof(float));
	ansp0=(float*)calloc(winsize*2,sizeof(float));
	ansp1=(float*)calloc(winsize*2,sizeof(float));
	ansp2=(float*)calloc(winsize*2,sizeof(float));
	ansA=(float*)calloc(winsize*2,sizeof(float));
	ansC=(float*)calloc(winsize*2,sizeof(float));
	ansG=(float*)calloc(winsize*2,sizeof(float));
	ansU=(float*)calloc(winsize*2,sizeof(float));
	anssump1=(float*)calloc(winsize*2,sizeof(float));
	anssump2=(float*)calloc(winsize*2,sizeof(float));
	
	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;

	
	
	/* beginning of fft part */
	/* correlation c(p0) */
	
	i = 0; 
	while(item_ptr1!=NULL){
		data1[i]=item_ptr1->p0;
		item_ptr1=item_ptr1->previous_ptr;
		++i;
	}

	i = 0; 
	while(item_ptr2!=NULL){
		data2[i]=item_ptr2->p0;
		item_ptr2=item_ptr2->previous_ptr;	
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansp0-1);
	
	printf("\n");
	
	
	/* correlation c(p1) */	
	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;
		i = 0; 
	while(item_ptr1!=NULL){
		data1[i]=item_ptr1->p1;
		item_ptr1=item_ptr1->previous_ptr;
		++i;
	}

	i = 0; 
	while(item_ptr2!=NULL){
		data2[i]=item_ptr2->p1;
		item_ptr2=item_ptr2->previous_ptr;	
		++i;
	}
	

	
	correl(data1-1,data2-1,winsize,ansp1-1);
	
	printf("\n");
	
	/* correlation c(p2) */		
	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;
	i = 0; 
	while(item_ptr1!=NULL){
		data1[i]=item_ptr1->p2;
		item_ptr1=item_ptr1->previous_ptr;
		++i;
	}

	i = 0; 
	while(item_ptr2!=NULL){
		data2[i]=item_ptr2->p2;
		item_ptr2=item_ptr2->previous_ptr;	
		++i;
	}
	

	correl(data1-1,data2-1,winsize,ansp2-1);
	
	printf("\n");
	
	
	/* correlation c(A) */			
	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;
	i = 0; 
	while(item_ptr1!=NULL){
		data1[i]=(item_ptr1->replbase==1 ? 1 : 0);
		item_ptr1=item_ptr1->previous_ptr;
		++i;
	}

	i = 0; 
	while(item_ptr2!=NULL){
		data2[i]=(item_ptr2->replbase==1 ? 1 : 0);
		item_ptr2=item_ptr2->previous_ptr;	
		++i;
	}
	

	correl(data1-1,data2-1,winsize,ansA-1);
	
	printf("\n");


	/* correlation c(C) */	
	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;
	i = 0; 
	while(item_ptr1!=NULL){
		data1[i]=(item_ptr1->replbase==1 ? 1 : 0);
		item_ptr1=item_ptr1->previous_ptr;
		++i;
	}

	i = 0; 
	while(item_ptr2!=NULL){
		data2[i]=(item_ptr2->replbase==1 ? 1 : 0);
		item_ptr2=item_ptr2->previous_ptr;	
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansC-1);

	printf("\n");



	/* correlation c(G) */	
	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;
	i = 0; 
	while(item_ptr1!=NULL){
		data1[i]=(item_ptr1->replbase==1 ? 1 : 0);
		item_ptr1=item_ptr1->previous_ptr;
		++i;
	}

	i = 0; 
	while(item_ptr2!=NULL){
		data2[i]=(item_ptr2->replbase==1 ? 1 : 0);
		item_ptr2=item_ptr2->previous_ptr;	
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansG-1);
	
	printf("\n");



	/* correlation c(U) */		
	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;
	i = 0; 
	while(item_ptr1!=NULL){
		data1[i]=(item_ptr1->replbase==1 ? 1 : 0);
		item_ptr1=item_ptr1->previous_ptr;
		++i;
	}

	i = 0; 
	while(item_ptr2!=NULL){
		data2[i]=(item_ptr2->replbase==1 ? 1 : 0);
		item_ptr2=item_ptr2->previous_ptr;	
		++i;
	}
	
		
	correl(data1-1,data2-1,winsize,ansU-1);

	/* end of fft part */


	item_ptr1=old1_ptr;
	item_ptr2=old2_ptr;
	
	for (i=0; i < winsize; ++i){
		if (item_ptr1){
			printf("%c",item_ptr1->base);
			item_ptr1=item_ptr1->previous_ptr;
		} else {
			printf("0");
		}
	}	

	printf("\n");
	for (i=0; i < winsize; ++i){
		if (item_ptr2){
			printf("%c",item_ptr2->base);
			item_ptr2=item_ptr2->previous_ptr;
		} else {
			printf("0");
		}
	}	
	printf("\n");



	for (i=0; i < winsize; ++i){
		/*char *array1,*array2;
		int j;*/
		
		anssump1[i]=ansA[i]+ansC[i]+ansG[i]+ansU[i]+alpha*ansp1[i];
		anssump2[i]=ansA[i]+ansC[i]+ansG[i]+ansU[i]+alpha*ansp2[i];
	}
	
	
/*********** ein verbesserungswürdiger sortieralgorihmus */	
	
	
	
	
	
	
	array=(int *)calloc(20,sizeof(int));

	for (i = 0; i < 20; ++i){
		int j=0;
		max = -1;
		elmax = -1; 
		for (j = 0; j < winsize; ++j){
			elmax = (anssump1[j] > max ?           j  : elmax);
			max   = (anssump1[j] > max ? 	anssump1[j] :   max);
		}
		
		anssump1[elmax]=0;
		if (elmax > winsize/2){
			elmax-=winsize;	
		}
		array[i]=(-elmax)%winsize;
		
	}
	
	for (i = 0; i < 20; ++i){
		printf("arrayp1[%d]:%d\n",i,array[i]);	
	}

/*********************************************************/

	
	smatrix=(struct segment **)calloc(sdata[last1_ptr->seqID].seqlength,sizeof(struct segment *));
	for (i=0; i < 	sdata[last1_ptr->seqID].seqlength; ++i){
		smatrix[i]=(struct segment *)calloc(sdata[last2_ptr->seqID].seqlength,sizeof(struct segment ));
	}
	
	
	/* do the sliding wondow analysis */
	segmentscorepair(array,last1_ptr,last2_ptr);

	
	

	

	/*********** sortieren zum zweiten **************/

	for (i = 0; i < 20; ++i){
		int j=0;
		max = -1;
		elmax = -1; 
		for (j = 0; j < winsize; ++j){
			elmax = (anssump2[j] > max ?           j  : elmax);
			max   = (anssump2[j] > max ? 	anssump2[j] : max  );
		}
		anssump2[elmax]=0;
		if (elmax > winsize/2){
			elmax-=winsize;	
		}
		array[i]=(-elmax)%winsize;
	}

	for (i = 0; i < 20; ++i){
		printf("arrayp2[%2d]:%d\n",i,array[i]);	
	}	

	/************************************************/


	
	
	segmentscorepair(array,last1_ptr,last2_ptr);
	




	printf("\n    ");
	
	for (i=0; i < sdata[last1_ptr->seqID].seqlength; ++i){
		printf("%7d ",i+1);
	}
	printf("\n"); 

	

	
	for (j=0; j < sdata[last2_ptr->seqID].seqlength; ++j){
		printf("%2d: ",j+1);
		for (i=0; i < sdata[last1_ptr->seqID].seqlength; ++i){
			printf("%7.1f ",smatrix[i][j].score);
		}
		printf("\n");
	}
	

	
	/* build a segment score matrix */
	findmaximizingpaths(sdata[last1_ptr->seqID].seqlength,sdata[last2_ptr->seqID].seqlength);
	
	free(smatrix);
	free(array);
	free(ansp0);
	free(ansp1);
	free(ansp2);
	free(ansA);
	free(ansC);
	free(ansG);
	free(ansU);
	free(anssump1);
	free(anssump2);
	
}

	
/* Function:		fftmultalign
 * 
 * Purpose:		calculates correlation between two sequences
 *           
 * Args:			last1_ptr, last2_ptr: pointer to the first sequence item
 *
 * Returns:		void
 */
void fftmultalign(struct profile p, int lengthSEQ1, int lengthSEQ2, int nseqs){
	float *data1, *data2;
	float *ansp0, *ansp1, *ansp2, *ansA, *ansC, *ansG, *ansU, *anssump1, *anssump2;
	int i,j,k;
	int winsize;
	int sizeA,sizeB;
	int elmax;
	float max;
	int *array;
	int paths;
	
	struct vector **ptrarrayA;
	struct vector **ptrarrayB;
	
	float **arrayI;
	float **arrayJ;
	
	ct = 0;	
	ptrarrayA = (struct vector **)malloc(p.nA*sizeof(struct vector *));
	ptrarrayB = (struct vector **)malloc(p.nB*sizeof(struct vector *));
	
	arrayI = (float **)calloc(lengthSEQ1+1,sizeof(float *));
	arrayJ = (float **)calloc(lengthSEQ2+1,sizeof(float *));
	
	
	for (k = 0; k <= lengthSEQ1; ++k){
		arrayI[k] = (float *)calloc(9,sizeof(float));
	}
	for (k = 0; k <= lengthSEQ2; ++k){
		arrayJ[k] = (float *)calloc(9,sizeof(float));
	}
	

	for (k = 0; k < p.nA;++k ){
		ptrarrayA[k] = last_ptr[p.seqprofA[k]-1];
	
		data_content(last_ptr[p.seqprofA[k]-1],arrayI);
	
		/*sumWa += sdata[ptrarrayA[k]->seqID].weight; */
	

	}
	

	for (k = 0; k < p.nB;++k ){
		ptrarrayB[k] = last_ptr[p.seqprofB[k]-1];
	
		data_content(last_ptr[p.seqprofB[k]-1],arrayJ);
	
		/*sumWb += sdata[ptrarrayB[k]->seqID].weight;*/
	}
	
	
	
	
	for (i=1; i < 12; ++i){
		if (lengthSEQ1 % (int)pow(2.0,(double)i)==lengthSEQ1){
			sizeA=(int)pow(2.0,(double)i);
			break;
		}
	}
	
	for (i=1; i < 12; ++i){
		if (lengthSEQ2 % (int)pow(2.0,(double)i)==lengthSEQ2){
			sizeB=(int)pow(2.0,(double)i);
			break;
		}
	}
	
	winsize=2*(sizeA > sizeB ? sizeA : sizeB);
	
	
	printf("SeqID-1:%s, length: %d, sizeA: %d\n",sdata[ptrarrayA[0]->seqID].ID, lengthSEQ1, sizeA);
	printf("SeqID-2:%s, length: %d, sizeB: %d\n",sdata[ptrarrayB[0]->seqID].ID, lengthSEQ2, sizeB);	
	printf("winsize: %d\n",winsize);
	
	
	data1=(float*)calloc(winsize,sizeof(float));
	data2=(float*)calloc(winsize,sizeof(float));
	ansp0=(float*)calloc(winsize*2,sizeof(float));
	ansp1=(float*)calloc(winsize*2,sizeof(float));
	ansp2=(float*)calloc(winsize*2,sizeof(float));
	ansA=(float*)calloc(winsize*2,sizeof(float));
	ansC=(float*)calloc(winsize*2,sizeof(float));
	ansG=(float*)calloc(winsize*2,sizeof(float));
	ansU=(float*)calloc(winsize*2,sizeof(float));
	anssump1=(float*)calloc(winsize*2,sizeof(float));
	anssump2=(float*)calloc(winsize*2,sizeof(float));
	
	
	
	/* correlation c(p0) */
	
	i = 0; 
	while(i < lengthSEQ1){
		data1[i]=arrayI[i][6]/p.nA;
		++i;
	}

	i = 0; 
	while(i < lengthSEQ2){
		data2[i]=arrayJ[i][6]/p.nB;
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansp0-1);
	
	printf("\n");
	
	
	/* correlation c(p1) */
	
	i = 0; 
	while(i < lengthSEQ1){
		data1[i]=arrayI[i][7]/p.nA;
		++i;
	}

	i = 0; 
	while(i < lengthSEQ2){
		data2[i]=arrayJ[i][7]/p.nB;
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansp1-1);
	
	printf("\n");
	
	
	/* correlation c(p2) */
	
	i = 0; 
	while(i < lengthSEQ1){
		data1[i]=arrayI[i][8]/p.nA;
		++i;
	}

	i = 0; 
	while(i < lengthSEQ2){
		data2[i]=arrayJ[i][8]/p.nB;
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansp2-1);
	
	printf("\n");


	/* correlation c(A) */
	
	i = 0; 
	while(i < lengthSEQ1){
		data1[i]=arrayI[i][0]/p.nA;
		++i;
	}

	i = 0; 
	while(i < lengthSEQ2){
		data2[i]=arrayJ[i][0]/p.nB;
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansA-1);
	
	printf("\n");


	/* correlation c(C) */
	
	i = 0; 
	while(i < lengthSEQ1){
		data1[i]=arrayI[i][0]/p.nA;
		++i;
	}

	i = 0; 
	while(i < lengthSEQ2){
		data2[i]=arrayJ[i][0]/p.nB;
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansC-1);
	
	printf("\n");


	/* correlation c(G) */
	
	i = 0; 
	while(i < lengthSEQ1){
		data1[i]=arrayI[i][2]/p.nA;
		++i;
	}

	i = 0; 
	while(i < lengthSEQ2){
		data2[i]=arrayJ[i][2]/p.nB;
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansG-1);
	
	printf("\n");


	/* correlation c(U) */
	
	i = 0; 
	while(i < lengthSEQ1){
		data1[i]=arrayI[i][3]/p.nA;
		++i;
	}

	i = 0; 
	while(i < lengthSEQ2){
		data2[i]=arrayJ[i][3]/p.nB;
		++i;
	}
	
	correl(data1-1,data2-1,winsize,ansU-1);
	
	printf("\n");




	for (i=0; i < winsize; ++i){
		anssump1[i]=ansA[i]+ansC[i]+ansG[i]+ansU[i]+alpha*ansp1[i];
		anssump2[i]=ansA[i]+ansC[i]+ansG[i]+ansU[i]+alpha*ansp2[i];
		printf("anssump1[%d]: %f\n",i, anssump1[i]);
		printf("anssump2[%d]: %f\n",i, anssump2[i]);
		fflush(stdout);
	}


/****************** sortieren der correlation ************/
	array=(int *)calloc(20,sizeof(int));
	for (i = 0; i < 20; ++i){
		int j=0;
		max = -1;
		elmax = -1; 
		for (j = 0; j < winsize; ++j){
			elmax = (anssump1[j] > max ? j : elmax);
			max = (anssump1[j] > max ? 	anssump1[j] : max);
		}
		
		anssump1[elmax]=0;
		if (elmax > winsize/2){
			elmax-=winsize;	
		}
		array[i]=(-elmax)%winsize;
		
	}
	
	for (i = 0; i < 20; ++i){
		printf("arrayp1[%d]:%d\n",i,array[i]);
		fflush(stdout);
	}
	
/****************** sortieren der correlation ************/	

	smatrix=(struct segment **)calloc(lengthSEQ1,sizeof(struct segment *));
	for (i=0; i < 	lengthSEQ1; ++i){
		smatrix[i]=(struct segment *)calloc(lengthSEQ2,sizeof(struct segment ));
	}


	segmentscoremult(array,ptrarrayA,ptrarrayB,p.nA,p.nB,p);



	for (i = 0; i < 20; ++i){
		int j=0;
		max = -1;
		elmax = -1; 
		for (j = 0; j < winsize; ++j){
			elmax = (anssump2[j] > max ? j : elmax);
			max = (anssump2[j] > max ? 	anssump2[j] : max);
		}
		anssump2[elmax]=0;
		if (elmax > winsize/2){
			elmax-=winsize;	
		}
		array[i]=(-elmax)%winsize;
	}
	
	for (i = 0; i < 20; ++i){
		fflush(stdout);
		printf("arrayp2[%d]:%d\n",i,array[i]);
			
	}
	segmentscoremult(array,ptrarrayA,ptrarrayB,p.nA,p.nB,p);	
	

	printf("\n    ");
	
	for (i=0; i < lengthSEQ1; ++i){
		printf("%5d ",i+1);
		fflush(stdout);
	}
	printf("\n"); 

	

	
	for (j=0; j < lengthSEQ2; ++j){
		printf("%2d: ",j+1);
		fflush(stdout);
		for (i=0; i < lengthSEQ1; ++i){
			fflush(stdout);
			printf("%5.2f ",smatrix[i][j].score);
		}
		printf("\n");
	}
	



	/* build a segment score matrix */
	paths = findmaximizingpaths(lengthSEQ1,lengthSEQ1);


	for (i=paths-1; i >= 0; i--){
		/*retrievePath(pathContainer[i]->startI,pathContainer[i]->startJ);*/
	}
	
	printf("\npaths: %d\n",paths);
	
	
	
	
	retrievePath(pathContainer[paths-1]->startI,pathContainer[paths-1]->startJ);


	
	
	/*
	p die Sequenzen die alignieren sind
	blmerge.count Anzahl der Blöcke
	i welcher Block
	segptr startzeiger
	
	alignpieces(p,i,blmerge.count,segptr);
	*/
	
	
	
	free(pathContainer);
	free(array);
	free(ansp0);
	free(ansp1);
	free(ansp2);
	free(ansA);
	free(ansC);
	free(ansG);
	free(ansU);
	free(anssump1);
	free(anssump2);
	ct = 1;

}



/* Function:		segmentscorepair
 * 
 * Purpose:		sliding window analysis to retrieve score between two sequences of a segment with lag k
 *           
 * Args:			array: contents lags of the 20 highest correlations
 *				
 * Returns:		void
 */
void segmentscorepair(int *array, struct vector *last1_ptr, struct vector *last2_ptr){
	int i,j;
	int lengthSEQ1, lengthSEQ2;
	struct vector *item_ptr1,*item_ptr2;
	
		
	

	for (i = 0; i < 20; ++i){
		int lag;
		lag = array[i];
		lengthSEQ1 = sdata[last1_ptr->seqID].seqlength;
		lengthSEQ2 = sdata[last2_ptr->seqID].seqlength;
		
		lengthSEQ1=(array[i] < 0 ? lengthSEQ1 + lag:lengthSEQ1);
		lengthSEQ2=(array[i] > 0 ? lengthSEQ2 - lag:lengthSEQ2);
		
		item_ptr1 = last1_ptr;
		item_ptr2 = last2_ptr;
		
		if (lag > 0){
			while(lag!=0 && item_ptr2!=NULL){
				item_ptr2=item_ptr2->previous_ptr;
				lag--;
			}
		
		}
		if (lag < 0){
			while(lag!=0 && item_ptr1!=NULL){
				item_ptr1=item_ptr1->previous_ptr;
				lag++;
			}
		
		}
		for(j=0; j < min(lengthSEQ1 ,lengthSEQ2); j+=stepSize){

			int k=0;
			float fftscore=0;
			char *segmentA;
			char *segmentB;

			segmentA = (char *)calloc(8,sizeof(char));
			segmentB = (char *)calloc(8,sizeof(char));
			while(k < stepSize){
				segmentA[k]=item_ptr1->base;
				segmentB[k]=item_ptr2->base;			
				fftscore+=score(item_ptr1,item_ptr2);		

				item_ptr1=item_ptr1->previous_ptr;
				if(item_ptr1==NULL){
					++k;
					break;
				}
				
				item_ptr2=item_ptr2->previous_ptr;
				if(item_ptr2==NULL){
					++k;
					break;
				}
				++k;
			}
			
			int beginA,beginB;
			beginA = (array[i] < 0 ? -array[i] : 0);
			beginB = (array[i] > 0 ?  array[i] : 0);
				
			if ((fftscore/k) > (0.7 *stepSize) /*&& k == stepSize*/){
				smatrix[j+beginA][j+beginB].score=fftscore/stepSize;
				smatrix[j+beginA][j+beginB].beginA=j+beginA;
				smatrix[j+beginA][j+beginB].endA=j+beginA+k-1;
				smatrix[j+beginA][j+beginB].beginB=j+beginB;
				smatrix[j+beginA][j+beginB].endB=j+beginB+k-1;
				smatrix[j+beginA][j+beginB].blocks=k;
				printf("beginA: %2d, endA: %2d   ",j+beginA, j+k+beginA-1);	
				printf("beginB: %2d, endB: %2d   ",j+beginB, j+k+beginB-1);	
				printf("score: %f ",fftscore);	
				printf("score: %f\n",fftscore/k);	
			}
			free(segmentA);
			free(segmentB);
		}
	}	
}




/* Function:		segmentscoremult
 * 
 * Purpose:		sliding window analysis to retrieve score between two groups of sequences of a segment with lag k
 *           
 * Args:			array: contents lags of the 20 highest correlations
 *				
 * Returns:		void
 */
void segmentscoremult(int *array,struct vector **ptrarrayA,struct vector **ptrarrayB, int a, int b, struct profile p){
	int i,j,l;
	int stepSize=6;
	int lengthSEQ1, lengthSEQ2;

	for (i = 0; i < 20; ++i){
		int lag;
		lag = array[i];
		
		
		for (l = 0; l < p.nA;++l ){
			ptrarrayA[l] = last_ptr[p.seqprofA[l]-1];
		}
	

		for (l = 0; l < p.nB;++l ){
			ptrarrayB[l] = last_ptr[p.seqprofB[l]-1];
		}
				
		
		lengthSEQ1 = sdata[ptrarrayA[0]->seqID].seqlength;
		lengthSEQ2 = sdata[ptrarrayB[0]->seqID].seqlength;
		
		lengthSEQ1=(array[i] < 0 ? lengthSEQ1 + lag:lengthSEQ1);
		lengthSEQ2=(array[i] > 0 ? lengthSEQ2 - lag:lengthSEQ2);
		
		printf("\nl1: %d, l2: %d\n",lengthSEQ1,lengthSEQ2);
		
		
		
		
		if (lag > 0){
			while(lag!=0){
				for (l=0; l < b; ++l){
					ptrarrayB[l]=ptrarrayB[l]->previous_ptr;
				}
				lag--;
			}
		}

		if (lag < 0){
			while(lag!=0){
				for (l=0; l < a; ++l){
					ptrarrayA[l]=ptrarrayA[l]->previous_ptr;
				}
				lag++;
			}
		}
	
		for(j=0; j < min(lengthSEQ1 ,lengthSEQ2); j+=stepSize){
			int k=0;
			float fftmultscore=0;
			char **segmentA;
			char **segmentB;
			
			segmentA = (char **)calloc(p.nA,sizeof(char*));
			segmentB = (char **)calloc(p.nB,sizeof(char*));
			
			printf("\n\narray[%d]:lag %d\n",i,array[i]);	

			for (l=0; l < a; ++l){
				segmentA[l] = (char *)calloc(stepSize,sizeof(char));					
			}
			
			for (l=0; l < b; ++l){
				segmentB[l] = (char *)calloc(stepSize,sizeof(char));					
			}
			
			while(k < stepSize){
				for (l=0; l < a; ++l){
					segmentA[l][k]=ptrarrayA[l]->base;
				}
			
				for (l=0; l < b; ++l){
					segmentB[l][k]=ptrarrayB[l]->base;
				}
				
				
/*				scorePP(struct vector **ptrarrayA,struct vector **ptrarrayB,struct profile p, float **scorematrix, struct lscore *spA, struct lscore *spB)*/
				fftmultscore+=scorePPfft(ptrarrayA,ptrarrayB, a, b);
				
				for (l=0; l < a; ++l){
					ptrarrayA[l]=ptrarrayA[l]->previous_ptr;
				}
				
				if(ptrarrayA[0]==NULL){
					++k;
					break;
				}
				
				
				for (l=0; l < b; ++l){
					ptrarrayB[l]=ptrarrayB[l]->previous_ptr;
				}

				if(ptrarrayB[0]==NULL){
					++k;
					break;
				}
				++k;
			}

			int beginA,beginB;
			beginA = (array[i] < 0 ? -array[i] : 0);
			beginB = (array[i] > 0 ?  array[i] : 0);			
			
			if ((fftmultscore/k) > (0.7 * stepSize)) {
				smatrix[j+beginA][j+beginB].score=fftmultscore/stepSize;
				smatrix[j+beginA][j+beginB].beginA=j+beginA;
				smatrix[j+beginA][j+beginB].endA=j+beginA+k-1;
				smatrix[j+beginA][j+beginB].beginB=j+beginB;
				smatrix[j+beginA][j+beginB].endB=j+beginB+k-1;
				smatrix[j+beginA][j+beginB].blocks=k;
				printf("beginA: %2d, endA: %2d   ",j+beginA, j+k+beginA-1);	
				printf("beginB: %2d, endB: %2d   ",j+beginB, j+k+beginB-1);	
				printf("score: %f ",fftmultscore);	
				printf("score: %f\n",fftmultscore/k);	
			}
			
			for (l=0; l < a; ++l){
				printf("segA: %s\n",segmentA[l]);
				free(segmentA[l]);
			}
			
			for (l=0; l < b; ++l){
				printf("segB: %s\n",segmentB[l]);
				free(segmentB[l]);
			}
			
			printf("score: %f",fftmultscore/(k*(p.nA+p.nB)));	
			free(segmentA);
			free(segmentB);
		}
	}
}


	
int findmaximizingpaths(int lenA,int lenB) {
	int i,j;
	float max;
	int count = 0;

	/* begin last row last column */
	


	for (j = lenB-1; j >=0; j--){
		
		for (i = lenA-1; i >=0; i--){
			if (smatrix[i][j].score!=0) {
				count++;
				smatrix[i][j].father=NULL;
				if ((i==(lenA-1)) || (j==(lenB-1)) ){
					smatrix[i][j].linkScore=smatrix[i][j].score;
					smatrix[i][j].child=NULL;
					
				} else if (i < lenA - stepSize && j < lenB -stepSize){
					lookForChildren(i,j,lenA,lenB);
				} else {
					smatrix[i][j].linkScore=smatrix[i][j].score;
					
				}
			} else {
				
			}
		}
		
	}
	
	printf("\n        ");
	for (i = 0 ; i < lenA; i++){
		printf("%2d      ",i+1);
	}

	printf("\n"); 

	for (j = 0; j < lenB; j++){
		printf("%2d ",j+1);
		for (i = 0 ; i < lenA; i++){
			printf("%7.1f ",smatrix[i][j].linkScore);
		}
		printf("\n");
	}
	
	return getAllPaths(count,lenA,lenB);
}

int getAllPaths(int counts,int lenA, int lenB){
	struct listPath *startItem;
	struct listPath *newItem;
	struct listPath *iterator;
	
	
	
	int i,j,paths=0;
	
	startItem=(struct listPath *)calloc(counts,sizeof(struct listPath ));
	
	startItem->score =0;
	startItem->father=NULL;
	startItem->child =NULL;
	
	
	for (i = lenA-1; i >=0; i--){
		for (j = lenB-1; j >=0; j--){
			if (smatrix[i][j].linkScore){
				
				/* iterate over list */						
				newItem=(struct listPath *)calloc(counts,sizeof(struct listPath ));			
				iterator = startItem;
				while(iterator->child!=NULL && iterator->score < smatrix[i][j].linkScore) {
					iterator= iterator->child;
				} 

				newItem->score = smatrix[i][j].linkScore;
				newItem->startpoint = &smatrix[i][j];
				newItem->startI = i;
				newItem->startJ = j;				
				newItem->father = iterator->father;
				if (iterator->father!=NULL){
					iterator->father->child = newItem; 
				}
				iterator->father = newItem;
				newItem->child = iterator;

				startItem = moveToStart(newItem);
				paths++;
			}
		}
	}
	
	pathContainer = (struct listPath **)calloc(paths, sizeof(struct listPath *));

	paths=0;
	iterator = startItem;
	while (iterator->child!=NULL){
		pathContainer[paths]=iterator;
		iterator = iterator->child;
		paths++;
	}
	return paths;
}



struct listPath *moveToStart(struct listPath *newItem){
	struct listPath *iterator;
	iterator = newItem;
	while(iterator->father!=NULL){
		iterator = iterator->father;	
	}
	return iterator;
}



void retrievePath(int i, int j){	
	int count = 1;
	struct segment *firstItem;
	firstItem = &smatrix[i][j];
	printf("\n\nfirstItem: %p and his child %p, firstItem.linkScore=%5.2f, startA %d, startB %d  ",firstItem,firstItem->child,firstItem->linkScore,firstItem->beginA,firstItem->beginB);
	getNext(firstItem, count);
	
}

struct segment *getNext(struct segment *treeItem , int count) {
	if (treeItem->child!=NULL) {
		printf("\nnextItem : %p and his child %p,",treeItem->child,treeItem->child->child);		
		printf("  nextItem.linkScore=%5.2f, startA %d, startB %d  ",treeItem->child->linkScore,treeItem->child->beginA,treeItem->child->beginB);			
		count++;
		printf("count %d",count);			
		getNext(treeItem->child,count);
	}
	return NULL;
}

void lookForChildren(int i, int j,int lenA, int lenB){
	int k,l;
	int boolean = 0;
	int linkK, linkL;
	int counter = 0;
	float maxscore = 0;
	struct segment *maxTreeItem;

	for ( l = j+stepSize; l < lenB ; l++){
		for ( k = i+stepSize; k < lenA ; k++){
			fflush(stdout);
		
			counter++;
			if(checkIsValidChildren(i,j,k,l)){
				boolean = 1;
				if (maxscore < smatrix[i][j].score + smatrix[k][l].linkScore ){
					fflush(stdout);
					smatrix[i][j].child = &smatrix[k][l];
					smatrix[i][j].linkScore = smatrix[i][j].score + smatrix[k][l].linkScore;
					fflush(stdout);
					if (getValidHasChildren(k,l)) { 
						maxTreeItem= smatrix[k][l].child; 
					}
				}	
				smatrix[i][j].child = (maxscore > smatrix[i][j].score + smatrix[k][l].linkScore ? smatrix[i][j].child :  &smatrix[k][l]);
				smatrix[i][j].linkScore = (maxscore > smatrix[i][j].score + smatrix[k][l].linkScore ? smatrix[i][j].linkScore : smatrix[i][j].score + smatrix[k][l].linkScore);
				maxscore = (maxscore > smatrix[i][j].score + smatrix[k][l].linkScore ? maxscore : smatrix[i][j].score + smatrix[k][l].linkScore);
			}
		}
	}

	if (!boolean) {
		smatrix[i][j].linkScore=smatrix[i][j].score;
	}
	printf("%7.1f ",smatrix[i][j].linkScore);
}


int checkIsValidChildren(int i, int j, int k, int l){ 

	if (smatrix[i][j].endA < smatrix[k][l].beginA && smatrix[i][j].endB < smatrix[k][l].beginB)
		return 1;
	return 0;
}



int getValidHasChildren(int k, int l) {
	if(smatrix[k][l].child!=NULL)
		return 1;
	return 0;
}



/* Function:		four1
 * 
 * Purpose:		replaces the contents in data array by its (inverse) discrete Fourier transform
 *           
 * Args:			data: array, which contents are alternating real resp imaginary floats
 *				nn: number of complex numbers (must be an integer power of 2)
 *				isign: if isign=1, the discrete Fourier transform is calculated, else isign=-1, the inverse 
 *				discrete Fourier transform is calculated
 *
 * Returns:		void
 */
void four1( float *data, int nn, int isign ){
	unsigned long n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	float tempr, tempi;

	n = nn << 1;
	j = 1;
	for( i=1; i<n; i+=2 )
	{
		if( j > i )
		{
			SWAP( data[j], data[i] );
			SWAP( data[j+1], data[i+1] );
		}
		m = n >> 1;
		while( m >= 2 && j > m )
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while( n > mmax )
	{
		istep = mmax << 1;
		theta = isign * ( 2 * PI / mmax );
		wtemp = sin( 0.5 * theta );
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin( theta );
		wr = 1.0;
		wi = 0.0;
		for( m=1; m<mmax; m+=2 )
		{
			for( i=m; i<=n; i+=istep )
			{
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr * data[j+1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = ( wtemp = wr ) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}






/* Function:		realft
 * 
 * Purpose:		real data arrays become discrete Fourier transformed
 *           
 * Args:			data1, data2: real data array 
 *				fft1,fft2: complex data array		
 *				n: number of numbers n *data (must be an integer power of 2)
 *
 * Returns:		void
 */
void realft(float data[], int n, int isign){
	int i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;
	
	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2; i <=(n>>2); i++){
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4]= -h1i+wr*h2i-wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign==1){
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1] = c1*((h1r=data[1])+data[2]);
		data[2] = c1*(h1r-data[2]);
		four1(data,n>>1,-1);	
	}
}




/* Function:		twofft
 * 
 * Purpose:		real data arrays become discrete Fourier transformed
 *           
 * Args:			data1, data2: real data array 
 *				fft1,fft2: complex data array		
 *				n: number of numbers n *data (must be an integer power of 2)
 *
 * Returns:		void
 */
void twofft(float data1[], float data2[],float fft1[], float fft2[],int n){
	int nn3, nn2, jj, j;
	float rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	
	for (j=1, jj=2; j<=n; j++,jj+=2){
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	
	four1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	
	for (j=3; j <=n+1; j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;	
		fft1[nn3-j]=-aim;
		fft2[j]=aip;
		fft2[j+1]=-rem;
		fft2[nn2-j]=aip;	
		fft2[nn3-j]=rem;
	}
}




/* Function:		correl
 * 
 * Purpose:		computes the correlation of the two data sets
 *           
 * Args:			data1,data2: array, which contents are alternating real resp imaginary floats
 *				n: number of complex numbers (must be an integer power of 2)
 *				ans: array which holds the lags
 *
 * Returns:		void
 */
void correl(float data1[],float data2[], int n, float ans[]){
	
	
	int no2,i; 
	float dum, *fft;
	
	fft=nrvector(1,n<<1);
	twofft(data1, data2, fft, ans, n);
	no2=n>>1;
	for (i=2; i <=n+2; i+=2){
		ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
		ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;	
	}	
	ans[2]=ans[n+1];
	realft(ans,n,-1);
	free_nrvector(fft,1,n<<1);
}
