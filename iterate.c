#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "stral.h"
#include "align.h"
#include "iterate.h"
#include "alfuncs.h"
#include "aliterate.h"
#include "iofuncs.h"
#include "create_matrix.h"
/* Function:		doiteration
 * 
 * Purpose:		remove one sequence randomly from alignment and realign her
 *           
 * Args:			n: number of sequences in alignment
 *				score: aligment score
 *				rems: number of sequences to be removed at once
 *
 * Returns:		score
 */
float doiteration(int n, int rems, float score){
	
	int i,j,k;
	struct vector ***segptr;
	float mxs, **remat;
	
	int *removes;
	
	struct profile p;
	struct vector *curr_ptr, *old_ptr;
	struct blocks blhel,blseg,blmerge;




	/* pick number rems sequences randomly from all n sequences */
	rems = ((rand())% (n-1))+1;
	removes = (int*)malloc(rems*sizeof(int));
	/*
	remat = retree(n);
	
	distmatrix =(float **)malloc(nseqs*sizeof(float *));
	
	
	for (i = 0; i < n; i++){
		distmatrix[i] = (float *)malloc(nseqs* sizeof(float)) ;
	}
	
	
	createDistanceMatrix(remat, n);
	createTree(n);
	exit (1);
	*/
	
	for (i=0; i < rems; ++i){
		int rn;
		int j, boolean;
		rn = (rand()% n);
		boolean =1;
		j=0;
		while(j < i){
			if (rn==removes[j]){
				boolean=0;
				j=i;
			}
			j++;
		}
		if (boolean){
			removes[i]=rn;
		} else {
			--i;
		}
	}

		
	
	
	/* construct a "guide" profile */
	p.nA = n-rems;
	p.nB = rems;
	p.typeA=(n-rems ==1 ? 0 : 1);
	p.typeB=(rems ==1 ? 0 : 1);
	p.seqprofA = (int *)malloc((n-rems)*sizeof(int));
	p.seqprofB = (int *)malloc(rems*sizeof(int));
	p.pnoA=0;
	p.pnoB=0;
	
	
	for (j = 0; j < rems ;++j){
		p.seqprofB[j]=removes[j]+1;
	}
	k=0;	
	for (i = 0; i < n; ++i){
		int boolean = 1;
		for (j = 0; j < rems ;++j){ 
			if (i+1==p.seqprofB[j]){
				boolean = 0;
				break;
			}
		}
		if (boolean){
			p.seqprofA[k]=i+1;
			++k;
	 	}
	}
	
	copyvectors(n);
	
	
	/*****************************************************
	* Test the achored alignment                         *
	*****************************************************/
	if (palign){
		blhel = phelix(n);
		blseg = findsegment(n);
		blmerge = mergeblocks(blhel, blseg);
		segptr = (struct vector ***)calloc(blmerge.count+1,sizeof(struct vector **));
		segptr[0] = (struct vector **)calloc(n,sizeof(struct vector *));
		for ( i = 0; i < blmerge.count; ++i){
			segptr[i+1] = (struct vector **)calloc(n,sizeof(struct vector *));
			for (k = 0; k < n;  ++k){
				segptr[0][k] = lastcpyptr[k];						
				segptr[i+1][k] = lastcpyptr[k];						
				for (j = 1; j < blmerge.M[i]; j++){
					segptr[i+1][k] = segptr[i+1][k]->previous_ptr;						
				}
			}	
		}
		/*
		for ( i = 0; i <= blmerge.count; ++i){
			for (k = 0; k < n;  ++k){
				fprintf(stdout,"pos: %d, char %c\n",segptr[i][k]->pos,segptr[i][k]->base);
			}
			fprintf(stdout,"\n");
		}
		*/
	}
	
	remgapcolsfromcpy(p.seqprofB,rems);
	remgapcolsfromcpy(p.seqprofA,n-rems);
	
	
	
	for (i =0; i < n; ++i){
		print_data_content2(lastcpyptr[i]);
#ifdef DEBUG
		if (!verbose){
			printf("\n");
		}
#endif
	}

	
	
	
	if (palign && blmerge.count > 0){
		for ( i = 0; i <= blmerge.count; i+=2){
			alignpieces(p,i,blmerge.count,segptr);
		}

	for (j =0; j < n; ++j){
		print_data_content2(last_ptr[j]);
#ifdef DEBUG
			if (!verbose){
			printf("\n");
		}
#endif
	}
			
	} else {	
		mxs = methodmultcpy(p , firstcpyptr[p.seqprofA[0]-1]->pos , firstcpyptr[p.seqprofB[0]-1]->pos,n);
		for (i = 0; i < n; ++i){
			curr_ptr = lastcpyptr[i];
			while(curr_ptr!=NULL){
				old_ptr=curr_ptr->previous_ptr;
				free(curr_ptr);
				curr_ptr = old_ptr;
			}
		}
		free(lastcpyptr);
		free(firstcpyptr);
	}
	
	if (score < mxs){
		remgapcols(p.seqprofB, rems);
		remgapcols(p.seqprofA, n - rems);
		
		for (i =0; i < n; ++i){
			print_data_content2(last_ptr[i]);
#ifdef DEBUG
			if (!verbose){
				printf("\n");
			}
#endif
		}
		ct=0;
		score = methodmult(p , first_ptr[p.seqprofA[0]-1]->pos , first_ptr[p.seqprofB[0]-1]->pos,n,i);
	}
	free(removes);
	free(p.seqprofA);
	free(p.seqprofB);
	return(score);
}



/* Function:	adjustalpha
 * 
 * Purpose:		adjust alpha value according to average pairwise identity after aligning the sequences
 *           
 * Args:			pwident:	pairwise identity					
 *
 * Returns:		void
 */
void adjustalpha(float pwident, int n){
	
	int i;
	
	if (pwident > 0.5){
		alpha = -16 * pwident +16;
	} else {
		alpha = 8;	
	}
	
	for (i = 0; i < n; ++i){
		removegaps(i, 1);
	}
	printf("alpha: %f\n", alpha);
	
	ct=0;
	for (i= n-2; i>=0;--i){
		methodmult(prof[i] , sdata[prof[i].seqprofA[0]-1].seqlength , sdata[prof[i].seqprofB[0]-1].seqlength,n,i);
	}
}



/* Function:	bestfirst
 * 
 * Purpose:		removes each sequence once from alignment
 *					the vectors of the sequences are copied and 
 *					the removed sequence is going to be realigned
 *					the alignment scoring best is kept 
 *           
 * Args:			n: number of sequences in alignment
 *          
 * Returns:		void
 */
void bestfirst(int n, int rems){
	int i = 0;
	int j;
	int bf=0;
	float mxs;
	float	mscore = 0;
	struct profile p;
	struct vector *curr_ptr, *old_ptr;
	struct blocks blhel,blseg;
	
	
	
	blhel = phelix(n);
	blseg = findsegment(n);
	
	for (i = 0; i < blseg.count; ++i){
		fprintf(stdout,"segment b %d, e %d\n",blseg.L[i],blseg.R[i]);
		fprintf(stdout,"segment middle %d\n",blseg.M[i]);
	}
	for (i = 0; i < blhel.count; ++i){
		fprintf(stdout,"helix b %d, e %d\n",blhel.L[i],blhel.R[i]);
		fprintf(stdout,"helix middle %d\n",blhel.M[i]);
	}
	
	for (i=0; i < n; ++i){
	
		copyvectors(n);
		
		/* construct an input profile */
		p.nA = n-1;
		p.nB = 1;
		p.typeA=1;
		p.typeB=0;
		p.seqprofA = (int *)malloc((n-1)*sizeof(int));
		p.seqprofB = (int *)malloc(sizeof(int));
		p.seqprofB[0] = i+1;
		p.pnoA=0;
		p.pnoB=0;
		for (j = 0; j < n; ++j){
			if (j<i) {
				p.seqprofA[j] = j+1;
			} else if (j > i){
				p.seqprofA[j-1]=j+1;			
			}
		}		
		remgapcolsfromcpy(p.seqprofB, 1);
		remgapcolsfromcpy(p.seqprofA, n-1);
		
		for (j =0; j < n; ++j){
		print_data_content2(lastcpyptr[j]);
#ifdef DEBUG
			if (!verbose){
			printf("\n");
			}
#endif
		}
		
		
		mxs = methodmultcpy(p , firstcpyptr[(i+1)%n]->pos , firstcpyptr[i]->pos,n);
		bf = (mxs > mscore ? i : bf);
		mscore = (mxs > mscore ? mxs : mscore);
		for (j = 0; j < n; ++j){
			curr_ptr = lastcpyptr[j];
			while(curr_ptr!=NULL){
				old_ptr=curr_ptr->previous_ptr;
				free(curr_ptr);
				curr_ptr = old_ptr;
			}
		}
		free(lastcpyptr);
		free(firstcpyptr);
		free(p.seqprofA);
		free(p.seqprofB);
	}
	p.nA = n-1;
	p.nB = 1;
	p.typeA=1;
	p.typeB=0;
	p.seqprofA = (int *)malloc((n-1)*sizeof(int));
	p.seqprofB = (int *)malloc(sizeof(int));
	p.seqprofB[0] = bf+1;
	p.pnoA=0;
	p.pnoB=0;
	for (j = 0; j < n; ++j){
		if (j<bf) {
			p.seqprofA[j] = j+1;
		} else if (j > bf){
			p.seqprofA[j-1]=j+1;			
		}
	}		
	remgapcols(p.seqprofB, 1);
	remgapcols(p.seqprofA, n-1);
	
	for (j =0; j < n; ++j){
		print_data_content2(last_ptr[j]);
#ifdef DEBUG
			if (!verbose){
			printf("\n");
		}
#endif
	}
	printf("\nbf %d\n", bf+1);
	ct=0;
	methodmult(p , first_ptr[(bf+1)%n]->pos , first_ptr[bf]->pos,n,0);
	free(p.seqprofA);
	free(p.seqprofB);
}


/* Function:	removegapsfromcopy
 * 
 * Purpose:		removes gaps from the removed sequence
 *					instead of working on the original vectors
 *					they are to be copied
 *           
 * Args:			rn: number of sequence removed
 *          	rems: number of sequences to be removed at once
 *
 * Returns:		void
 */
void removegapsfromcpy(int *removes,int rems){
	struct vector *current_ptr;
	struct vector *item_ptr;
	int i; 
	
	for (i=0; i < rems; ++i){
		current_ptr = lastcpyptr[removes[i]];
		while(current_ptr!=NULL) {
			if (current_ptr->base=='-'){
				if (current_ptr==lastcpyptr[removes[i]]){
					lastcpyptr[removes[i]] = lastcpyptr[removes[i]]->previous_ptr;
					lastcpyptr[removes[i]]->next_ptr = NULL;
				} else if (current_ptr != firstcpyptr[removes[i]]){
					item_ptr = current_ptr->next_ptr; 
					item_ptr->previous_ptr = current_ptr->previous_ptr;
					item_ptr->previous_ptr->next_ptr = item_ptr;
				} else {
					firstcpyptr[removes[i]] = current_ptr->next_ptr;
					firstcpyptr[removes[i]]->previous_ptr = NULL;
					current_ptr = current_ptr->next_ptr;			
				}	
			}
			current_ptr = current_ptr->previous_ptr;
		}
		print_data_content2(lastcpyptr[removes[i]]);
#ifdef DEBUG
		if (!verbose){
			printf("\n");
		}
#endif
	}
	
}



/* Function:	removegapcolsfromcopy
 * 
 * Purpose:		removes gap only columns from alignment
 *					remaining after removing sequence(rm)
 *					instead of working on the original vectors
 *					they are to be copied
 *           
 * Args:			profileA:	sequences remainig after random pick
 *					n:		number of sequuences in profileA
 *          	
 * Returns:		void
 */
void remgapcolsfromcpy(int *profileA, int n){
	struct vector **ptrarray;
	struct vector *item_ptr;
	int i;
	ptrarray = (struct vector **)malloc(n*sizeof(struct vector *));
	
	for (i=0; i < n; ++i){
		ptrarray[i] = lastcpyptr[profileA[i]-1];
		
	}
	
	while(ptrarray[0]!=NULL){
		int boolean =1;
		for (i = 0; i < n; ++i){
			if (ptrarray[i]->replbase!=5){
				boolean = 0;
				i = n;
			}
		}
	
		if (boolean){
			for (i=0; i < n; ++i){
				if (ptrarray[i]==lastcpyptr[profileA[i]-1]){
					
					item_ptr = ptrarray[i];
					lastcpyptr[profileA[i]-1] = lastcpyptr[profileA[i]-1]->previous_ptr;
					lastcpyptr[profileA[i]-1]->next_ptr = NULL;
					/*free(item_ptr);*/
				} else if (ptrarray[i] != firstcpyptr[profileA[i]-1]){
					item_ptr = ptrarray[i]->next_ptr; 
					item_ptr->previous_ptr = ptrarray[i]->previous_ptr;
					item_ptr->previous_ptr->next_ptr = item_ptr;
					item_ptr = ptrarray[i];
					/*free(item_ptr);*/
				} else {
					firstcpyptr[profileA[i]-1] = ptrarray[i]->next_ptr;
					firstcpyptr[profileA[i]-1]->previous_ptr = NULL;
					item_ptr = ptrarray[i];
					ptrarray[i] = ptrarray[i]->next_ptr;			
					/*free(item_ptr);*/
				}
			}	
		}
		
		for (i =0; i < n; ++i){
			ptrarray[i]= ptrarray[i]->previous_ptr;
		}	
	}
	free(ptrarray);
}



/* Function:	removegaps
 * 
 * Purpose:		removes gaps from the removed sequence
 *           
 * Args:			rn: number of sequence removed
 *          	rems: number of sequences to be removed at once
 *
 * Returns:		void
 */
void removegaps(int rn, int rems){
	struct vector *current_ptr;
	struct vector *item_ptr;
	current_ptr = last_ptr[rn];
	while(current_ptr!=NULL) {
		if (current_ptr->base=='-'){
			if (current_ptr==last_ptr[rn]){
				last_ptr[rn] = last_ptr[rn]->previous_ptr;
				last_ptr[rn]->next_ptr = NULL;
			}
			else if (current_ptr != first_ptr[rn]){
				item_ptr = current_ptr->next_ptr; 
				item_ptr->previous_ptr = current_ptr->previous_ptr;
				item_ptr->previous_ptr->next_ptr = item_ptr;
			} else {
				first_ptr[rn] = current_ptr->next_ptr;
				first_ptr[rn]->previous_ptr = NULL;
				current_ptr = current_ptr->next_ptr;			
			}
		}
	
		current_ptr = current_ptr->previous_ptr;
	}
	
	print_data_content2(last_ptr[rn]);
#ifdef DEBUG
	if (!verbose){
		printf("\n");
	}
#endif
}

/* Function:	removegapcols
 * 
 * Purpose:		removes gap only columns from alignment
 *					remaining after removing sequence(rm)
 *         
 * Args:			profileA:	sequences remainig after random pick
 *					n:		number of sequuences in profileA
 *
 * Returns:		void
 */
void remgapcols(int *profile, int n){
	struct vector **ptrarray;
	struct vector *item_ptr;
	int i;
	ptrarray = (struct vector **)malloc(n*sizeof(struct vector *));
	
	
	for (i=0; i < n; ++i){
		ptrarray[i] = last_ptr[profile[i]-1];
	}

	
	while(ptrarray[0]!=NULL){
		int boolean = 1;
		for (i =0; i < n; ++i){
			/*fprintf(stderr, "i=%d ptrarray=%p\n", i, ptrarray); fflush(stderr);
   		fprintf(stderr, "ptrarray[i=%d]=%p\n", i, ptrarray[i]); fflush(stderr);
			fprintf(stderr, "ptrarray[%d]=%p ->replbase=%d\n", i, ptrarray[i], ptrarray[i]->replbase); fflush(stderr);*/
			if (ptrarray[i]->replbase!=5){
				boolean = 0;
				break;
			}
		}
		if (boolean){
			for (i=0; i < n; ++i){
				/*printf("ptrarray[%d]->seqID:%d,last_ptr[%d]->seqID=%d\n",i,ptrarray[i]->seqID,profile[i]-1,last_ptr[profile[i]-1]->seqID);*/
				if (ptrarray[i]==last_ptr[profile[i]-1]){
					last_ptr[profile[i]-1] = last_ptr[profile[i]-1]->previous_ptr;
					last_ptr[profile[i]-1]->next_ptr = NULL;
				}
				else if (ptrarray[i] != first_ptr[profile[i]-1]){
					item_ptr = ptrarray[i]->next_ptr; 
					item_ptr->previous_ptr = ptrarray[i]->previous_ptr;
					item_ptr->previous_ptr->next_ptr = item_ptr;
				} else {
					first_ptr[profile[i]-1] = ptrarray[i]->next_ptr;
					first_ptr[profile[i]-1]->previous_ptr = NULL;
					ptrarray[i] = ptrarray[i]->next_ptr;			
				}
			}	
		}
		
		for (i =0; i < n; ++i){
			ptrarray[i]= ptrarray[i]->previous_ptr;
		}	
	}
	free(ptrarray);
}

/* Function:	copyvectors
 * 
 * Purpose:		clones each vector("deep" copy) 
 *         
 * Args:			n:		number of sequences in alignment
 *          
 * Returns:		void
 */
void copyvectors(int n){
	int i=0;
	struct vector *curr_ptr, *new_item_ptr, *item_ptr, *old_ptr;
	lastcpyptr  = (struct vector **)calloc(n, sizeof(struct vector *));
	firstcpyptr = (struct vector **)calloc(n, sizeof(struct vector *));
	for (i=0; i < n; ++i){
		curr_ptr = last_ptr[i];
		new_item_ptr = malloc(sizeof(struct vector));
		new_item_ptr->p0 = curr_ptr->p0;
		new_item_ptr->p1 = curr_ptr->p1;
		new_item_ptr->p2 = curr_ptr->p2;
		new_item_ptr->seqID = curr_ptr->seqID;
		new_item_ptr->replbase = curr_ptr->replbase;
		new_item_ptr->pos = curr_ptr->pos;
		new_item_ptr->base = curr_ptr->base;
		new_item_ptr->oribase = curr_ptr->oribase;
		new_item_ptr->next_ptr = NULL;
		item_ptr = malloc(sizeof(struct vector));
		new_item_ptr->previous_ptr = item_ptr;
		lastcpyptr[i] = new_item_ptr;
		new_item_ptr = item_ptr;
		old_ptr = new_item_ptr;
		curr_ptr = curr_ptr->previous_ptr;			
		while (curr_ptr->previous_ptr!=NULL){
			new_item_ptr->p0 = curr_ptr->p0;
			new_item_ptr->p1 = curr_ptr->p1;
			new_item_ptr->p2 = curr_ptr->p2;
			new_item_ptr->seqID = curr_ptr->seqID;
			new_item_ptr->replbase = curr_ptr->replbase;
			new_item_ptr->pos = curr_ptr->pos;
			new_item_ptr->base = curr_ptr->base;
			new_item_ptr->oribase = curr_ptr->oribase;
			new_item_ptr->next_ptr = old_ptr;
			item_ptr = malloc(sizeof(struct vector));
			old_ptr = new_item_ptr;
			new_item_ptr->previous_ptr = item_ptr;
			new_item_ptr=item_ptr;					
			curr_ptr = curr_ptr->previous_ptr;
		}
		new_item_ptr->p0 = first_ptr[i]->p0;
		new_item_ptr->p1 = first_ptr[i]->p1;
		new_item_ptr->p2 = first_ptr[i]->p2;
		new_item_ptr->seqID = first_ptr[i]->seqID;
		new_item_ptr->replbase = first_ptr[i]->replbase;
		new_item_ptr->pos = first_ptr[i]->pos;
		new_item_ptr->base = first_ptr[i]->base;
		new_item_ptr->oribase = first_ptr[i]->oribase;
		new_item_ptr->next_ptr = old_ptr;
		new_item_ptr->previous_ptr = NULL;
		firstcpyptr[i]=new_item_ptr;
	}
}


struct blocks mergeblocks(struct blocks blhel, struct blocks blseg){
	int i,j,k,boolean,*L,*R,*Lm,*Rm,*Mm;
	struct blocks blmerge,blblock;
	
	L = (int *)calloc(1,(sizeof (int)));
	R = (int *)calloc(1,(sizeof (int)));	
	
	
	/*
	for (i = 0; i < blseg.count; ++i){
		fprintf(stdout,"segment[%2d] begins at %3d and ends on %3d\n",i+1,blseg.L[i],blseg.R[i]);
	}
	for (i = 0; i < blhel.count; ++i){
		fprintf(stdout,"  helix[%2d] begins at %3d and ends on %3d\n",i+1,blhel.L[i],blhel.R[i]);
	}
	*/

	/* merge possible helix parts and high homology parts */
	boolean=1;
	
	i=0;j=0;k=0;
	while(boolean){
		int begseg,endseg;
		int beghel,endhel;
		
		
		if (j == blhel.count ){
			while ( i < blseg.count){
				begseg = blseg.L[i];
				endseg = blseg.R[i];
				L = (int *)realloc(L,(k+2)*(sizeof (int)));
				L[k] = begseg;
				R = (int *)realloc(R,(k+2)*(sizeof (int)));
				R[k] = endseg;
				++i;
				++k;		
			}
		} else if (i == blseg.count){
			while ( j < blhel.count){
				beghel=blhel.L[j];
				endhel=blhel.R[j];
				L = (int *)realloc(L,(k+2)*(sizeof (int)));
				L[k] = beghel;
				R = (int *)realloc(R,(k+2)*(sizeof (int)));
				R[k] = endhel;
				++j;
				++k;	
			}		
		}
		
		begseg = blseg.L[i];beghel=blhel.L[j];
		endseg = blseg.R[i];endhel=blhel.R[j];
		
		if (j == blhel.count && i == blseg.count){
			boolean=0;
			continue;
		}
		if (begseg < beghel){
			/*fprintf(stdout,"begseg %d < %d beghel\n",begseg,beghel);*/
			L = (int *)realloc(L,(k+2)*(sizeof (int)));
			L[k] = begseg;
			/*fprintf(stdout,"L[%d]=%d\n",k,L[k]);*/
			R = (int *)realloc(R,(k+2)*(sizeof (int)));
			R[k] = endseg;
			++i;
			++k;		
		}
		if (j == blhel.count && i == blseg.count){
			boolean=0;
			continue;
		}
		if (begseg > beghel){
			/*fprintf(stdout,"begseg %d > %d beghel\n",begseg,beghel);*/
			L = (int *)realloc(L,(k+2)*(sizeof (int)));
			L[k] = beghel;
			/*fprintf(stdout,"L[%d]=%d\n",k,L[k]);*/
			R = (int *)realloc(R,(k+2)*(sizeof (int)));
			R[k] = endhel;
			++j;
			++k;		
		}
		if (j == blhel.count && i == blseg.count){
			boolean=0;
			continue;
		}
		if (begseg == beghel){
			if (endseg < endhel){
				L = (int *)realloc(L,(k+2)*(sizeof (int)));
				L[k] = begseg;
				R = (int *)realloc(R,(k+2)*(sizeof (int)));
				R[k] = endseg;
				++i;
				++k;		
			} else if (endseg > endhel){
				L = (int *)realloc(L,(k+2)*(sizeof (int)));
				L[k] = beghel;
				R = (int *)realloc(R,(k+2)*(sizeof (int)));
				R[k] = endhel;
				++i;
				++k;		
			} else if (endseg == endhel){
				L = (int *)realloc(L,(k+2)*(sizeof (int)));
				L[k] = beghel;
				R = (int *)realloc(R,(k+2)*(sizeof (int)));
				R[k] = endhel;
				++i;
				++j;
				++k;		
			}
		}
	}
	
	blblock.L = L;
	blblock.R = R;
	blblock.count = k;

	/*
	for (i = 0; i < blblock.count; ++i){
		fprintf(stdout,"  block[%2d] begins at %3d and ends on %3d\n",i+1,blblock.L[i],blblock.R[i]);
	}
	*/




	k = 0;
	
	fprintf(stdout,"\n");
	
	Lm = (int *)calloc(1,(sizeof (int)));
	Rm = (int *)calloc(1,(sizeof (int)));	
	Mm = (int *)calloc(1,(sizeof (int)));
	/* j is used for temporary storing R[i] */
	for (i = 0; i < blblock.count; ++i){
		/*fprintf(stdout,"  block[%2d] begins at %3d and ends on %3d\n",i+1,blblock.L[i],blblock.R[i]);*/
		
		Lm = (int *)realloc(Lm,(k+2)*(sizeof (int)));
		Rm = (int *)realloc(Rm,(k+2)*(sizeof (int)));	
		Mm = (int *)realloc(Mm,(k+2)*(sizeof (int)));	
		
		Lm[k] = blblock.L[i];
		Rm[k] = blblock.R[i];
		/*fprintf(stdout,"  merge[%2d] begins at %3d and ends on %3d\n",k,Lm[k],Rm[k]);*/
		
		label_1:
		if (blblock.L[i+1] >= Lm[k] && blblock.L[i+1] < Rm[k]){
			if (Rm[k] >= blblock.R[i+1]){
				++i;
				goto label_1;
			} else if (Rm[k] < blblock.R[i+1]){
				Rm[k] = blblock.R[i+1];
				++i;
				goto label_1;
			}
		}
		++k;
	}
	
	
	
	blmerge.L=Lm;
	blmerge.R=Rm;
	blmerge.count = k;
	for (i = 0; i < blmerge.count; ++i){
		fprintf(stdout,"  merge[%2d] begins at %3d and ends on %3d\n",i+1,blmerge.L[i],blmerge.R[i]);
		Mm[i]=(blmerge.L[i]+blmerge.R[i])/2;
		fprintf(stdout,"  merge[%2d]'s middle is %3d\n",i+1,Mm[i]);
	}

	blmerge.M=Mm;
		
		
	
	free(blblock.L);
	free(blblock.R);
	

	return blmerge;

}
