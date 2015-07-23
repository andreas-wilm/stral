#include <stdio.h>
#include <math.h>
#include "align.h"
#include "alfuncs.h"

#include "stral.h"
#include "create_matrix.h"

float recalculatescorePP(struct vector **ptrarray, int a){
	float value=0.;
	int k,l;
	struct vector *i_ptr, *j_ptr;
	float weight=1;

	if (a > 1){
		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				i_ptr = ptrarray[k];
				j_ptr = ptrarray[l];	
				weight = sqrt(sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight);
				value = value + weight * ( d[i_ptr->replbase][j_ptr->replbase]  + alpha * (sqrt(i_ptr->p1 * j_ptr->p1))	 + alpha * (sqrt(i_ptr->p2 * j_ptr->p2))); 
			}
		}
		return value;
	}
	return value;
}


/* Function:		symbols
 * 
 * Purpose:		calculate symbol presented in backtrack according to RNAfold's backtrack
 *           
 * Args:			i_ptr, j_ptr: pointers to the vectors of residue i in sequence s
 *
 * Returns:		char
 */
char *symbols(struct vector *i_ptr,struct vector *j_ptr){
	float x,y,z;
	x = sqrt(i_ptr->p0 * j_ptr->p0);
	y = sqrt(i_ptr->p1 * j_ptr->p1);
	z = sqrt(i_ptr->p2 * j_ptr->p2);	
	
	if ( x > 0.667)
		return ".";
	if ( y > 0.667) 
		return "(";
	if ( z > 0.667) 
		return ")";
	if ( (y + z) > x) {
		if ((y/(y+z)) > 0.667) 
			return "{";
		if ((z/(y+z)) > 0.667) 
			return "}";
		return "|";
	}
	if ( x > (y + z) ) 
		return ",";
	return ":";
}


/* Function:		symbolsPP
 * 
 * Purpose:		calculate symbol presented in backtrack (in msa) accordinig to RNAfold's backtrack
 *           
 * Args:			ptrarrayA, ptrarrayB: array of pointers to the vectors of residue i in sequence s in group g
 *				a, b: number of sequences in one group(array)
 *
 * Returns:		char
 */
char *symbolsPP(struct vector **ptrarrayA,struct vector **ptrarrayB,int a, int b){
	float x=0;
	float y=0;
	float z=0;
	int k,l,n;
	struct vector *i_ptr, *j_ptr;
	
#ifdef DEBUG	
	if (verbose){
		return (".");
	} else {	
		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];		    
				x = x + (sqrt(i_ptr->p0 * j_ptr->p0)) ;
			}
			for (l = 0; l < b; ++l)  {
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayB[l];
				x = x + (sqrt(i_ptr->p0 * j_ptr->p0));
			}					    
		}
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];
				x = x + (sqrt(i_ptr->p0 * j_ptr->p0));
			}
		}

		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];		    
				y = y + (sqrt(i_ptr->p1 * j_ptr->p1)) ;
			}
			for (l = 0; l < b; ++l)  {
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayB[l];
				y = y + (sqrt(i_ptr->p1 * j_ptr->p1));
			}					    
		}
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];
				y = y + (sqrt(i_ptr->p1 * j_ptr->p1));
			}
		}

		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];		    
				z = z + (sqrt(i_ptr->p2 * j_ptr->p2)) ;
			}
			for (l = 0; l < b; ++l)  {
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayB[l];
				z = z + (sqrt(i_ptr->p2 * j_ptr->p2));
			}					    
		}
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];
				z = z + (sqrt(i_ptr->p2 * j_ptr->p2));
			}
		}
		n = (a+b);
		x = x / ((n*(n-1))/2);
		y = y / ((n*(n-1))/2);
		z = z / ((n*(n-1))/2);

		if ( x > 0.667) return ".";
		if ( y > 0.667) return "(";
		if ( z > 0.667) return ")";
		if ( (y + z) > x) {
		if ((y/(y+z)) > 0.667) return "{";
		if ((z/(y+z)) > 0.667) return "}";
			return "|";
		}
		if ( x > (y + z) ) return ",";
		return ":";
	 }
#endif

#ifndef DEBUG
	return ".";
#endif
}

/* Function:		symbolsPPcpy
 * 
 * Purpose:		calculate symbol presented in backtrack (in msa) accordinig to RNAfold's backtrack (during iterative phase)
 *           
 * Args:			ptrarrayA, ptrarrayB: array of pointers to the vectors of residue i in sequence s in group g
 *				a, b: number of sequences in one group(array)
 *
 * Returns:		void
 */
char *symbolsPPcpy(struct vector **ptrarrayA,struct vector **ptrarrayB,int a, int b){
	float x=0;
	float y=0;
	float z=0;
	int k,l,n;
	struct vector *i_ptr, *j_ptr;

#ifdef DEBUG	
	if (verbose){
		return (".");
	} else {
		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];
				x = x + (sqrt(i_ptr->p0 * j_ptr->p0)) ;
			}
			for (l = 0; l < b; ++l)	{
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayB[l];
				x = x + (sqrt(i_ptr->p0 * j_ptr->p0));
			}						
		}
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];
		
				x = x + (sqrt(i_ptr->p0 * j_ptr->p0));
			}
		}
		
		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];			
				y = y + (sqrt(i_ptr->p1 * j_ptr->p1)) ;
			}
			for (l = 0; l < b; ++l)	{
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayB[l];
				
				y = y + (sqrt(i_ptr->p1 * j_ptr->p1));
			}						
		}
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];
				
				y = y + (sqrt(i_ptr->p1 * j_ptr->p1));
			}
		}
		
		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];			
				z = z + (sqrt(i_ptr->p2 * j_ptr->p2)) ;
			}
			for (l = 0; l < b; ++l)	{
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayB[l];
			
				z = z + (sqrt(i_ptr->p2 * j_ptr->p2));
			}						
		}
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];
				
				z = z + (sqrt(i_ptr->p2 * j_ptr->p2));
			}
		}
		
		n = (a+b);
		x = x / ((n*(n-1))/2);
		y = y / ((n*(n-1))/2);
		z = z / ((n*(n-1))/2);
	
	
	
		if ( x > 0.667) return ".";
		if ( y > 0.667) return "(";
		if ( z > 0.667) return ")";
		if ( (y + z) > x) {
			if ((y/(y+z)) > 0.667) return "{";
			if ((z/(y+z)) > 0.667) return "}";
      		return "|";
   		}
		if ( x > (y + z) ) return ",";
			return ":";
	}
#endif

#ifndef DEBUG
	return ".";
#endif
}


/* Function:		score
 * 
 * Purpose:		calculate score in pairwise alignment		
 *           
 * Args:			s1_ptr,s2_ptr: pointers to the vectors of residue i in sequence s
 *
 * Returns:		float
 */
float score(struct vector *s1_ptr, struct vector *s2_ptr){
	float value;
	value=  d[s1_ptr->replbase][s2_ptr->replbase] + alpha * ((sqrt( s1_ptr->p1 * s2_ptr->p1)) +  (sqrt(	s1_ptr->p2 *	s2_ptr->p2)));
	return value;
}



/* Function:		scorePP
 * 
 * Purpose:		calculate score in msa
 *           
 * Args:			ptrarrayA, ptrarrayB: array of pointers to the vectors of residue i in sequence s in group g
 *				a, b: number of sequences in one group(array)
 *				scorematrix: 2D array storing all calculated scores
 *
 * Returns:		void
 */
void scorePP(struct vector **ptrarrayA,struct vector **ptrarrayB,struct profile p, float **scorematrix, struct lscore *spA, struct lscore *spB){
	float value=0.;
	int k,l;
	struct vector *i_ptr, *j_ptr;
	float weight=1;
	int a = p.nA;
	int b = p.nB;
	
	for (k = 0; k < a ; ++k){	
		for (l = 0; l < b; ++l){
			
			i_ptr = ptrarrayA[k];
			j_ptr = ptrarrayB[l];		
			if (p.nA == 1 && p.nB == 1) {
				weight =sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight;
			}
			value = value + weight * (d[i_ptr->replbase][j_ptr->replbase] + alpha * (sqrt(i_ptr->p1 * j_ptr->p1))+ alpha * (sqrt(i_ptr->p2 * j_ptr->p2)));
		}						
	}
	
	if (ct){
		if (p.nA > 1){
			value += spA->score;
		}
	
		if (p.nB > 1){
			value += spB->score;
		}
	} else {
		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];			
				
				if (p.nA == 1 && p.nB == 1) {
					weight =sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight;
				}
				
					
				value = value + weight * (d[i_ptr->replbase][j_ptr->replbase] + alpha * (sqrt(i_ptr->p1 * j_ptr->p1)) + alpha * (sqrt(i_ptr->p2 * j_ptr->p2)));
			
			}
		}	
		
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
			
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];	
							
								
				if (p.nA == 1 && p.nB == 1) {
					weight =sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight;
				}
			
				value = value + weight * ( d[i_ptr->replbase][j_ptr->replbase]  + alpha * (sqrt(i_ptr->p1 * j_ptr->p1))	 + alpha * (sqrt(i_ptr->p2 * j_ptr->p2))); 
			}
		}
	}

	scorematrix[ptrarrayA[0]->pos][ptrarrayB[0]->pos]=value;

}


/* Function:		scorePP
 * 
 * Purpose:		calculate score in msa
 *           
 * Args:			ptrarrayA, ptrarrayB: array of pointers to the vectors of residue i in sequence s in group g
 *				a, b: number of sequences in one group(array)
 *				scorematrix: 2D array storing all calculated scores
 *
 * Returns:		void
 */
void scorePPcpy(struct vector **ptrarrayA,struct vector **ptrarrayB,struct profile p, float **scorematrix, struct lscore *spA, struct lscore *spB,int i, int j){
	float value=0.;
	int k,l;
	struct vector *i_ptr, *j_ptr;
	float weight=1;
	int a = p.nA;
	int b = p.nB;
	
	for (k = 0; k < a ; ++k){	
		for (l = 0; l < b; ++l){
			
			i_ptr = ptrarrayA[k];
			j_ptr = ptrarrayB[l];		
			if (p.nA == 1 && p.nB == 1) {
				weight =sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight;
			}
			fprintf(stdout,"d[%d][%d],d[%c][%c],d[%d][%d]\n",	i_ptr->pos, j_ptr->pos,i_ptr->base,j_ptr->base,i_ptr->seqID,j_ptr->seqID);
			value = value + weight * (d[i_ptr->replbase][j_ptr->replbase] + alpha * (sqrt(i_ptr->p1 * j_ptr->p1))+ alpha * (sqrt(i_ptr->p2 * j_ptr->p2)));
		}						
	}
	
	if (ct){
		if (p.nA > 1){
			value += spA->score;
		}
	
		if (p.nB > 1){
			value += spB->score;
		}
	} else {
		for (k = 0; k < a ; ++k){
			for (l = k + 1; l < a; ++l){
				
				i_ptr = ptrarrayA[k];
				j_ptr = ptrarrayA[l];			
				
				if (p.nA == 1 && p.nB == 1) {
					weight =sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight;
				}
				
					
				value = value + weight * (d[i_ptr->replbase][j_ptr->replbase] + alpha * (sqrt(i_ptr->p1 * j_ptr->p1)) + alpha * (sqrt(i_ptr->p2 * j_ptr->p2)));
			
			}
		}	
		
		for (k = 0; k < b ; ++k){
			for (l = k + 1; l < b; ++l){
			
				i_ptr = ptrarrayB[k];
				j_ptr = ptrarrayB[l];	
							
								
				if (p.nA == 1 && p.nB == 1) {
					weight =sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight;
				}
			
				value = value + weight * ( d[i_ptr->replbase][j_ptr->replbase]  + alpha * (sqrt(i_ptr->p1 * j_ptr->p1))	 + alpha * (sqrt(i_ptr->p2 * j_ptr->p2))); 
			}
		}
	}

	scorematrix[i][j]=value;

}



/* Function:		scoreSP
 * 
 * Purpose:		calculate score in msa
 *           
 * Args:			ptrarrayA, ptrarrayB: array of pointers to the vectors of residue i in sequence s in group g
 *				a, b: number of sequences in one group(array)
 *				scorematrix: 2D array storing all calculated scores
 *
 * Returns:		void
 */
float scoreSP(struct vector **ptrarray,int spp){
	float value=0.;
	int k,l;
	struct vector *i_ptr, *j_ptr;
	float weight=1;
	
	for (k = 0; k < spp ; ++k){	
		for (l = k+1; l < spp; ++l){
			
			i_ptr = ptrarray[k];
			j_ptr = ptrarray[l];		
			
			weight =sdata[i_ptr->seqID].weight*sdata[j_ptr->seqID].weight;
			
			value = value + weight * (d[i_ptr->replbase][j_ptr->replbase] );
		}						
	}
	

	return value;

}


struct blocks phelix(int spp){
	int alen,i,l,k=0;
	float *H,threshold;
	int *L,*R,*M,flag=0;
	struct blocks blhel;
	/*struct lscore *lptr;*/
	
	struct vector **ptrarray;

	ptrarray = (struct vector **)malloc(spp*sizeof(struct vector *));
	for (l = 0; l < spp;++l ){
		ptrarray[l] = last_ptr[l];
	}
	
	alen = sdata[1].seqlength;
	H = (float *)calloc(alen+1,sizeof(float));
	L = (int *)calloc(1,(sizeof (int)));
	R = (int *)calloc(1,(sizeof (int)));
	M = (int *)calloc(1,(sizeof (int)));
	H[0]=0;
	threshold = (alpha + 3.5)*0.7*((spp * (spp-1))/2);
	/*
	printf("threshold %f\n",threshold);
	printf("alen %d\n",alen);
	lptr = fps[1];*/


	for ( i = 1; i <= alen; ++i){
		H[i]=recalculatescorePP(ptrarray,spp);
		
				
		if (H[i] > threshold && H[i-1]==0) {
			L = (int *)realloc(L,(k+2)*(sizeof (int)));
			/*printf("Anfang\n");*/
			H[i-1]=0;
			L[k]=i;
			flag = 1;
		}
		/*fprintf(stdout,"H[%2d] %8.2f %8.2f\n",i,H[i],H[i-1]);*/
		if (H[i] < threshold){
			/*printf("%d-4=%d > L[%d]= %d\n",i,i-4,k,L[k]);*/
			if (flag && i-4 > L[k]){
				R = (int *)realloc(R,(k+2)*(sizeof (int)));
				M = (int *)realloc(M,(k+2)*(sizeof (int)));			
				R[k]=i-1;						
				M[k]=(R[k]+L[k])/2;
				flag=0;
				++k;
				/*printf("Ende\n");*/
			}	
			flag=0;
			H[i]=0;
		}
		for (l = 0; l < spp;++l ){
			ptrarray[l] = ptrarray[l]->previous_ptr;
		}

		
		/*lptr=lptr->next_ptr;*/
	}
	blhel.L = L;
	blhel.R = R;
	blhel.M = M;
	blhel.count = k;
		
	return blhel;
}

struct blocks findsegment(int spp){
	int alen,i,l,k=0;
	float *H,threshold;
	int *L,*R,*M,flag=0;
	struct blocks blseg;
	struct vector **ptrarray;
	ptrarray = (struct vector **)malloc(spp*sizeof(struct vector *));
	for (l = 0; l < spp;++l ){
		ptrarray[l] = last_ptr[l];
	}
	
	
	alen = sdata[1].seqlength;
	H = (float *)calloc(alen+1,sizeof(float));
	L = (int *)calloc(1,(sizeof (int)));
	R = (int *)calloc(1,(sizeof (int)));	
	M = (int *)calloc(1,(sizeof (int)));
	H[0]=0;
	threshold = 3.5*0.7*((spp * (spp-1))/2);
	/*printf("threshold %f\n",threshold);
	printf("alen %d\n",alen);*/
	

	for ( i = 1; i <= alen; ++i){
		H[i]=scoreSP(ptrarray,spp);
		
		if (H[i] > threshold && H[i-1]==0) {
			L = (int *)realloc(L,(k+2)*(sizeof (int)));
			/*printf("Anfang\n");*/
			H[i-1]=0;
			L[k]=i;
			flag = 1;
		}
		/*fprintf(stdout,"H[%d] %f %f\n",i,H[i],H[i-1]);*/
		if (H[i] < threshold){
			if (flag && ((i-4) > L[k])){
				R = (int *)realloc(R,(k+2)*(sizeof (int)));
				M = (int *)realloc(M,(k+2)*(sizeof (int)));							
				R[k]=i-1;
				M[k]=(R[k]+L[k])/2;						
				flag=0;
				++k;
				/*printf("Ende\n");*/
			}
			flag=0;
			H[i]=0;
		}
		for (l = 0; l < spp;++l ){
			ptrarray[l] = ptrarray[l]->previous_ptr;
		}

		
		
	}
	
	blseg.L = L;
	blseg.R = R;
	blseg.M = M;
	blseg.count = k;
	
	return blseg;
}


