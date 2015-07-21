#include <stdio.h>
#include "phylip.h"
#include "dist.h"
#include "stral.h"
#include "align.h"
#include "weighted_nj.h"



/* Function:		wjoin
 * 
 * Purpose:				
 *           
 * Args:			wtree:
 *				spp: number of sequences in alignment
 *
 * Returns:		void
 */
void wjoin(tree wtree, int spp){
	int i;
	int odd;
	int nonodes2;
	float brlgth;
	float left=0;
	float right=0;
	float y=0;
	float x=0;
	float value=0;
	float value1, minvalue1;
	float value2, minvalue2;
	float value3, minvalue3;
	float mindeltap = 10.e6;
	float maxW = 0;
	float minW = 10.e6;
	int divnode;
	int j;
	node *p, *q, *start1 , *start2, *start3, *wroot, *next1, *next2, *divide=NULL;

	/* setting negative and tiny branch lengths to zero */
	
	nonodes2 = spp * 2 - 2;


	
	for(i=0; i < nonodes2; ++i){
  		if (wtree.nodep[i]->v < 0.0001){
			wtree.nodep[i]->v = 0.0;
		}
		
#ifdef DEBUG
		if(!verbose){
			printf("node %d , ptr %p  ",i, wtree.nodep[i]);
			printf("backptr %p,  nextptr %p  ",wtree.nodep[i]->back, wtree.nodep[i]->next);
			printf("index %ld,  tip %d  iter %d  v %f ",wtree.nodep[i]->index,wtree.nodep[i]->tip,wtree.nodep[i]->iter, wtree.nodep[i]->v );
		}
#endif
		if (i>= spp){
		if (wtree.nodep[i]->next->v < 0.0001){
			wtree.nodep[i]->next->v = 0.0;
		}
		if (wtree.nodep[i]->next->next->v < 0.0001){
			wtree.nodep[i]->next->next->v = 0.0;
		}
#ifdef DEBUG
		if(!verbose){	
			printf("\nnextptr  %p,     backptr %p,  nextptr %p, index %ld,  tip %d, iter %d  v %f\n", wtree.nodep[i]->next,wtree.nodep[i]->next->back, wtree.nodep[i]->next->next,
			wtree.nodep[i]->next->index, wtree.nodep[i]->next->tip, wtree.nodep[i]->next->iter, wtree.nodep[i]->next->v);
			printf("nextptr->nextptr %p,backptr %p, nextptr %p, index %ld,  tip %d, iter %d  v %f",wtree.nodep[i]->next->next, wtree.nodep[i]->next->next->back, wtree.nodep[i]->next->next->next,
			wtree.nodep[i]->next->next->index, wtree.nodep[i]->next->next->tip, wtree.nodep[i]->next->next->iter,
			wtree.nodep[i]->next->next->v);
		}
#endif			
	  	}
#ifdef DEBUG
		if (!verbose){
			printf("\n\n");
		}
#endif		
	}

	
	
	
	
	brlgth = branchtravel(wtree.nodep[0], wtree );
	
	
	
	divnode = spp;
	for (i = spp; i < nonodes2; ++i){
		
		start1 = wtree.nodep[i];
		start2 = wtree.nodep[i]->next;
		start3 = wtree.nodep[i]->next->next;
		p = wtree.nodep[i];
	
	/* branch lenghts beginning from start1 travelling the tree */
		j = 0;
		odd = 1;
		
		value1 = 0;
		value2 = 0;
		value3 = 0;
		
		while(j!=1){
		
		if (odd){
			q = p->back;
			if(p == q->back){
				value1 = value1 + p->v;
			}
			
		} else if (!odd && (p->tip)) {
			q = p->back;
			odd = !odd;
			if(p == q->back){
				value1 = value1 + p->v;
			}
			
		} else {	
			q = p->next;
		}
		odd = !odd;
		p = q;
		if(q == start1)
			j = 1;
		}
#if DEBUG
		
		printf("value 1: %f\n",value1);
#endif
		minvalue1 = fabs(2*value1-brlgth);
		
		
	/* branch lenghts beginning from start2 travelling the tree */	
	
		p = start2;
		j = 0;
		odd = 1;
		
		while(j!=1){
		if (odd){
			q = p->back;
			if(p == q->back){
				value2 = value2 + p->v;
			}
			
		} else if (!odd && (p->tip)) {
			q = p->back;
			odd = !odd;
			if(p == q->back){
				value2 = value2 + p->v;
			}
			
		} else {	
			q = p->next;
		}
		odd = !odd;
		p = q;
		if(q == start2)
			j = 1;
	
		}
	minvalue2 = fabs(2*value2-brlgth);
	
	/* branch lenghts beginning from start3 travelling the tree */	
		
		p = start3;
		
		j = 0;
		odd = 1;
		
		while(j!=1){
		if (odd){
			q = p->back;
			if(p == q->back){
				value3 = value3 + p->v;
			}
			
		} else if (!odd && (p->tip)) {
			q = p->back;
			odd = !odd;
			if(p == q->back){
				value3 = value3 + p->v;
			}
			
		} else {	
			q = p->next;
		}
		odd = !odd;
		p = q;
		if(q == start3)
			j = 1;
	
		}

		minvalue3 = fabs(2*value3-brlgth);
	
		
		
		if (minvalue1 == minvalue2){
			int bool1 = 0;
			int bool2 = 0;
#ifdef DEBUG
			printf("one equals two: don't know what to do!!!!!!!!!!");
#endif			
			divnode = 1;
			value = value1;
			divide = start1;
			left  = brlgth - value,
			right = value - 2*divide->v,
			x = (right-left)/4+(divide->v)/2,
			y = (divide->v)-x;
#ifdef DEBUG
			printf("divnode 1: left %f, right %f, x %f, y %f\n",left,right,x,y);
#endif
			if (x > 0 && y > 0){
				bool1 = 1;
				printf("bool1 %d\n",bool1);
				minvalue3+=10000;
			}
			divnode = 3;
			value = value3;
			divide = start3;
			
			left  = value - 2*divide->v,
			right = brlgth - value,
			x = fabs(right-left)/4+(divide->v)/2,
			y = (divide->v)-x;
#ifdef DEBUG
			printf("divnode 3: left %f, right %f, x %f, y %f\n",left,right,x,y);	
#endif
			if (x > 0 && y > 0){
				bool2 = 1;
				printf("bool2 %d\n",bool2);
				minvalue1+=10000;
			}
			minvalue2 = 10e7;
			
			
			/* TODO */
		}
		if (minvalue1 == minvalue3){
			int bool1 = 0;
			int bool2 = 0;
#ifdef DEBUG						
			printf("one equals three: don't know what to do!!!!!!!!!!\n");
			fflush(stdout);
#endif			
			divnode = 1;
			value = value1;
			divide = start1;
			left  = brlgth - value,
			right = value - 2*divide->v,
			x = (right-left)/4+(divide->v)/2,
			y = (divide->v)-x;
#ifdef DEBUG
			printf("divnode 1: left %f, right %f, x %f, y %f\n",left,right,x,y);
#endif
			if (x > 0 && y > 0){
				bool1 = 1;
/*				printf("bool1 %d\n",bool1);*/
				minvalue3+=10000;
			}
			divnode = 3;
			value = value3;
			divide = start3;
			
			left  = value - 2*divide->v,
			right = brlgth - value,
			x = fabs(right-left)/4+(divide->v)/2,
			y = (divide->v)-x;
#ifdef DEBUG
			printf("divnode 3: left %f, right %f, x %f, y %f\n",left,right,x,y);	
#endif
			if (x > 0 && y > 0){
				bool2 = 1;
/*				printf("bool2 %d\n",bool2);*/
				minvalue1+=10000;
			}
			minvalue2 = 10e7;
			
			/* TODO */
		}
		if (minvalue2 == minvalue3){
			int bool1 = 0;
			int bool2 = 0;
#ifdef DEBUG			
			printf("two equals three: don't know what to do!!!!!!!!!!");
			fflush(stdout);
#endif			
			divnode = 2;
			value = value1;
			divide = start1;
			left  = brlgth - value,
			right = value - 2*divide->v,
			x = (right-left)/4+(divide->v)/2,
			y = (divide->v)-x;
#ifdef DEBUG
			printf("divnode 1: left %f, right %f, x %f, y %f\n",left,right,x,y);
#endif
			if (x > 0 && y > 0){
				bool1 = 1;
/*				printf("bool1 %d\n",bool1);*/
				minvalue3+=10000;
			}
			divnode = 3;
			value = value3;
			divide = start3;
			
			left  = value - 2*divide->v,
			right = brlgth - value,
			x = fabs(right-left)/4+(divide->v)/2,
			y = (divide->v)-x;
#ifdef DEBUG
			printf("divnode 3: left %f, right %f, x %f, y %f\n",left,right,x,y);	
#endif
			if (x > 0 && y > 0){
				bool2 = 1;
/*				printf("bool2 %d\n",bool2);*/
				minvalue2+=10000;
			}
			minvalue1 = 10e7;
			
			/* TODO */
		}
		
		
		
		if (mindeltap > minvalue1){
#ifdef DEBUG
			printf("mindeltap %f,minvalue1 %f,minvalue3 %f\n",mindeltap,minvalue1,minvalue3);
#endif
			mindeltap = minvalue1;
#ifdef DEBUG
			printf("A mindeltap %f,minvalue1 %f,minvalue3 %f\n",mindeltap,minvalue1,minvalue3);
#endif
			divide = start1;
			value = value1;
			divnode = 1;
			if (fabs(minvalue1 -minvalue2) <= 0.00001 && (fabs(value2 - value1) > 0.00001) && value2 > value1){    /* due to inaccuracy of floating-point arithmetic */
				divide = start2;
				value = value2;
				divnode = 2;
			}
		} 
		if (mindeltap > minvalue2 ) {
#ifdef DEBUG
			printf("mindeltap %f,minvalue1 %f,minvalue3 %f\n",mindeltap,minvalue1,minvalue3);
#endif
			mindeltap = minvalue2;
			divide = start2;
			value = value2;
			divnode = 2;			
		} 
		if (mindeltap > minvalue3){
#ifdef DEBUG
			printf("mindeltap %f,minvalue1 %f,minvalue3 %f\n",mindeltap,minvalue1,minvalue3);
#endif
			mindeltap = minvalue3;
			divide = start3;
			value = value3;		
			divnode = 3;
		} 
	
#ifdef DEBUG		
		printf( "node %d\n1  %f minp   %f\n2  %f minp   %f\n3  %f minp   %f\n\n",i, value1, minvalue1, value2, minvalue2 , value3, minvalue3);
#endif	
	}
	
	

	
	
#ifdef DEBUG
	if (!verbose){
		printf("divide %p mindeltap %f divnode %d", divide, mindeltap, divnode);
	}
#endif	
	/* calculate branch lenghts after dividing */
	
	if (divnode == 2)
		left  = value - 2*divide->v,
		right = brlgth - value,
		x = fabs(right-left)/4+(divide->v)/2,
		y = (divide->v)-x;
	if (divnode == 1)
		left  = brlgth - value,
		right = value - 2*divide->v,
		x = (right-left)/4+(divide->v)/2,
		y = (divide->v)-x;
	if (divnode == 3)
		left  = value - 2*divide->v,
		right = brlgth - value,
		x = fabs(right-left)/4+(divide->v)/2,
		y = (divide->v)-x;

#ifdef DEBUG
	if (!verbose){
		printf("\nleft: %f   right: %f  divide->v: %f   x: %f    y: %f\n",left, right, divide->v, x, y);		
	}
#endif
	
	/* allocating new root node */
	wroot = (node*)malloc(sizeof(node)); 
	next1 = (node*)malloc(sizeof(node));
	next2 = (node*)malloc(sizeof(node));
	
	wroot->tip = 0;
	wroot->next = next1;
	wroot->back = NULL;
	wroot->index = nonodes2+1;                    
	wroot->next->next = next2;             /* next1->next = next2 */
	next1->index = nonodes2+1;
	wroot->next->next->next = wroot;			 /* next2->next = wroot */	
	next2->index = nonodes2+1;
	if(x < 0 || y < 0) {
		if (divnode==1){
			float brl = divide->v;
			p = divide->back;
			divide->v = 0;
			divide->back = next1;
			next1->back = divide;
			next1->v = 0;
			next1->tip = 0;
			next2->back = p;
			next2->v = brl;
			next2->tip = 0;
			p->v = brl;
			p->back = next2;
		} else {
			float brl = divide->v;
			p = divide->back;
			p->v = brl ;
			p->back = next1;
			next1->back = p;
			next1->v = brl;
			next1->tip = 0;
			next2->back = divide;
			next2->v = 0;
			next2->tip = 0;
			divide->v = 0;
			divide->back = next2;
		}
	} else {
		if (divnode==1){
			p = divide->back;
			divide->v = x;
			divide->back = next1;
			next1->back = divide;
			next1->v = x;
			next1->tip = 0;
			next2->back = p;
			next2->v = y;
			next2->tip = 0;
			p->v = y;
			p->back = next2;
		} else {
			p = divide->back;
			p->v = x;
			p->back = next1;
			next1->back = p;
			next1->v = x;
			next1->tip = 0;
			next2->back = divide;
			next2->v = y;
			next2->tip = 0;
			divide->v = y;
			divide->back = next2;
		}
	}
	
	branchtravel(wroot, wtree);
	/*treeoutr(wroot,&col,&wtree);
	printree(wroot, 1	, !njoin , 1); 	*/
	
	/* weighting the sequences */
	if(!rooted){
		for (i = 0; i < spp ; ++i){
			p = wtree.nodep[i];
			sdata[i].weight = lookforleaf(wroot, wtree, p);
#ifdef DEBUG
			if (!verbose){
				printf ("seq %d with distance %f to root\n",i+1, sdata[i].weight);		
			}
#endif
	
			maxW = ( (maxW > sdata[i].weight) ? maxW : sdata[i].weight);
			minW = ( (minW < sdata[i].weight) ? minW : sdata[i].weight);
		}	
#ifdef DEBUG
		if (!verbose){
			printf ("maxW %f\n",maxW);		
		}
#endif		
		/* normalize weights */
	
		for (i = 0; i < spp ; ++i){
			sdata[i].weight = (sdata[i].weight/maxW);
			
#ifdef DEBUG			
			if (!verbose){
				printf ("seq %d with weight %f\n",i+1, sdata[i].weight);		
			}
#endif			
		}
	}
	
	setprofile(wtree, wroot); 	
}
	

/* Function:		branchtravel
 * 
 * Purpose:		calculates the length of the eularian path
 *           
 * Args:			start:
 *				wtree:
 *
 * Returns:		float brlgth
 */
float branchtravel( node *start , tree wtree){
	node *p, *q;
	int i, odd;
	float brlgth = 0;
	p = start;
	i = 0;
	odd = 0;
	while(i!=1){
		if (odd){
			q = p->back;
			if(p == q->back){
				brlgth = brlgth + p->v;
			}
			
		} else if (!odd && (p->tip)) {
			q = p->back;
			odd = !odd;
			if(p == q->back){
				brlgth = brlgth + p->v;
			}
			
		} else {	
			q = p->next;
		}
		odd = !odd;
		p = q;
		
		if(q == start)
			i = 1;
	
	}
#ifdef DEBUG
	if (!verbose){
		printf("\noverall branch length  %f \n",brlgth);
	}
#endif
	return brlgth;
}


/* Function:		
 * 
 * Purpose:		calculates the path from root to leaf		
 *           
 * Args:			wroot
 *				wtree:
 *				p
 *
 * Returns:		void
 */
float lookforleaf(node *wroot,tree wtree, node *p){
	int i=0;
	int k;
	node **nodearray;
	int *branchorder;
	int *numseq;
	float lengthtoleaf =0;
	node *q;
	int odd=1;
	
	q=wroot->next;
	nodearray = (struct node **)malloc(sizeof(struct node *));
	*(nodearray+i) = (struct node *)malloc(sizeof(struct node));
	nodearray[i] = wroot->next;
	branchorder = (int *)calloc(1,sizeof(int));
	numseq = (int *)calloc(1,sizeof(int));
	while(q!=p){
		
		if (odd && !(q->back->tip)) {
			int j;
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			
			/* avoid use of realloc */
			
			/*nodearray = (struct node **)calloc(i+2,sizeof(struct node *));
			
			free(nodearray);*/

			branchorder = (int *)realloc(branchorder,(i+2)*(sizeof (int)));
			numseq = (int *)realloc(numseq,(i+2)*(sizeof (int)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = nodearray[i]->back;
			branchorder[i+1] = branchorder[i]+1;
			++i;
			odd = !odd;
			q = q->back;
			j = findnode(nodearray, q, i);
						
			if (j < i){
			
				nodearray = (struct node **)realloc(nodearray,(j+2)*sizeof(struct node *));
				branchorder = (int *)realloc(branchorder,(j+2)*(sizeof (int)));
				numseq =  (int *)realloc(numseq,(j+2)*(sizeof (int)));
				nodearray[j] = q;
				i = j;
			} else {
				numseq[i] = countLeafs(q);
			}
		} 
		
		
		else if (odd && (q->back->tip) && (q->back!=p)){
			int j;
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			branchorder = (int *)realloc(branchorder,(i+2)*(sizeof(int)));
			numseq = (int *)realloc(numseq,(i+2)*(sizeof (int)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = nodearray[i]->next;
			branchorder[i+1] = branchorder[i];
			
			++i;
			q = q->next;
			j = findnode(nodearray, q, i);
			if (j < i){
				nodearray = (struct node **)realloc(nodearray,(j+2)*sizeof(struct node *));
				branchorder = (int *)realloc(branchorder,(j+2)*(sizeof (int)));
				numseq =  (int *)realloc(numseq,(j+2)*(sizeof (int)));
				nodearray[j] = q;
				i = j;
			}
		} 
		
		
		else if (odd && (q->back->tip) && (q->back==p)){
			
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			branchorder = (int *)realloc(branchorder,(i+2)*(sizeof(int)));
			numseq = (int *)realloc(numseq,(i+2)*(sizeof (int)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = nodearray[i]->back;
			branchorder[i+1] = branchorder[i]+1;
			numseq[i+1]=1;
			++i;
			q = q->back;

			
		} 
		
		
		else if (!odd){
			int j;
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			branchorder = (int *)realloc(branchorder,(i+2)*(sizeof(int)));
			numseq = (int *)realloc(numseq,(i+2)*(sizeof (int)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = nodearray[i]->next;
			branchorder[i+1] = branchorder[i];
			++i;
			odd = !odd;
			q = q->next;
			j = findnode(nodearray, q, i);
			if (j < i){
				nodearray = (struct node **)realloc(nodearray,(j+2)*sizeof(struct node *));
				branchorder = (int *)realloc(branchorder,(j+2)*(sizeof (int)));
				numseq =  (int *)realloc(numseq,(j+2)*(sizeof (int)));
				nodearray[j] = q;
				i = j;
			}
			
		}
		
	
	}
	
	
	for (k=0; k < i; ++k){

		if (branchorder[k]!=branchorder[k+1]){
			lengthtoleaf = lengthtoleaf + nodearray[k]->v/(numseq[k+1]);
		}
	}
	
	free(nodearray);
	free(branchorder);
	return lengthtoleaf;
}

/* Function:		findnode
 * 
 * Purpose:				
 *           
 * Args:			nodearray:
 *				p:
 *				i:
 *
 * Returns:		void
 */
int findnode(node **nodearray, node *p, int i){
	int j=0;
	for (j=0 ; j < i; ++j){
		if	(nodearray[j]==p)
			return j; 
	}
	return i;
}


/* Function:		setprofile
 * 
 * Purpose:				
 *           
 * Args:			wtree:
 *				wroot:
 *
 * Returns:		void
 */
void setprofile (tree wtree, node *wroot){
	int i,m, j = spp-2;
	int nonodes2;
	int k = 0;
	node **narray;
	node *p;
	narray = (struct node**)calloc(spp,sizeof(struct node *));
	prof = (struct profile *)calloc(spp,sizeof(struct profile));
	for (i = 0; i < spp-1; ++i) {
		narray[i] = (struct node *)calloc(1,sizeof(struct node));
		narray[i] = NULL;
	}
	
	nonodes2 = spp * 2 - 2;
		for(i=0; i < nonodes2; ++i){
  		if (wtree.nodep[i]->v < 0.0001){
			wtree.nodep[i]->v = 0.0;
		}
		
		
		
		
#ifdef DEBUG
		if(!verbose){
			printf("node %d , ptr %p  ",i, wtree.nodep[i]);
			printf("backptr %p,  nextptr %p  ",wtree.nodep[i]->back, wtree.nodep[i]->next);
			printf("index %ld,  tip %d  iter %d  v %f ",wtree.nodep[i]->index,wtree.nodep[i]->tip,wtree.nodep[i]->iter, wtree.nodep[i]->v );
		}
#endif
		if (i>= spp){
		if (wtree.nodep[i]->next->v < 0.0001){
			wtree.nodep[i]->next->v = 0.0;
		}
		if (wtree.nodep[i]->next->next->v < 0.0001){
			wtree.nodep[i]->next->next->v = 0.0;
		}
#ifdef DEBUG
		if(!verbose){	
			printf("\nnextptr  %p,     backptr %p,  nextptr %p, index %ld,  tip %d, iter %d  v %f\n", wtree.nodep[i]->next,wtree.nodep[i]->next->back, wtree.nodep[i]->next->next,
			wtree.nodep[i]->next->index, wtree.nodep[i]->next->tip, wtree.nodep[i]->next->iter, wtree.nodep[i]->next->v);
			printf("nextptr->nextptr %p,backptr %p, nextptr %p, index %ld,  tip %d, iter %d  v %f",wtree.nodep[i]->next->next, wtree.nodep[i]->next->next->back, wtree.nodep[i]->next->next->next,
			wtree.nodep[i]->next->next->index, wtree.nodep[i]->next->next->tip, wtree.nodep[i]->next->next->iter,
			wtree.nodep[i]->next->next->v);
		}
#endif			
	  	}
#ifdef DEBUG
		if (!verbose){
			printf("\n\n");
		}
#endif		
	}
	
	
#ifdef DEBUG
	if(!verbose){
		printf("wroot  %p , wroot->next %p ,  wroot->next->next %p ",wroot,wroot->next,wroot->next->next);
	}
#endif	
	
	
	
	
	
	
	
	/* find leafs clustering with a neighbored leaf */
	
	for (i = 0; i < spp; ++i){
	
		p = wtree.nodep[i];	
		if (!(p->back->next == wroot)){
			if (p->back->next->back->tip){
				/* cluster found */

				prof[j].seqprofA = (int *)calloc(1,sizeof(int));
				prof[j].seqprofA[0] = ( (i+1 < p->back->next->back->index) ? i+1 : p->back->next->back->index);
				prof[j].typeA = 0;
				prof[j].nA=1;
					
				prof[j].seqprofB = (int *)calloc(1,sizeof(int));
				prof[j].seqprofB[0] = ( (i+1 > p->back->next->back->index) ? i+1 : p->back->next->back->index);
				prof[j].typeB = 0;
				prof[j].nB=1;
			
				j--;
				narray[j] = p->back->next->next;
				k +=2;		
			}
		}
	}
	
	k = 0;	
	while(j >=0){
		k+=1;
		
		/* find leafs clustering with yet a clustered node */
		for (i = 0; i < spp; ++i){
			p = wtree.nodep[i];
			for (m = spp-3; m >= j ; --m){
				if (m < 0)
					break;
				if (p->back->next->back == narray[m] && narray[m]!=NULL){
					
					if (j<0)
						break;
					free(prof[j].seqprofA);
					prof[j].seqprofA = (int *)calloc(1,sizeof(int));
					prof[j].seqprofA[0] = ( (i+1 < prof[m+1].seqprofA[0]) ? i+1 : prof[m+1].seqprofA[0]);
					prof[j].typeA = ((i+1 < prof[m+1].seqprofA[0]) ? 0 : 1);
	   				prof[j].nA=1;
						
					prof[j].seqprofB = (int *)calloc(1,sizeof(int));
					prof[j].seqprofB[0] = ( (i+1 > prof[m+1].seqprofA[0]) ? i+1 : prof[m+1].seqprofA[0]);
					prof[j].typeB = ((i+1 < prof[m+1].seqprofA[0]) ? 1 : 0);
   					prof[j].nB=1;
					--j;
					if (j<0)
						break;
					narray[j] = p->back->next->next;
					narray[m]=NULL;
					
				} else if (p->back->next->next->back == narray[m] && narray[m]!=NULL){
		 			if (j<0)
						break;
					prof[j].seqprofA = (int *)calloc(1,sizeof(int));
					prof[j].seqprofA[0] = ( (i+1 < prof[m+1].seqprofA[0]) ? i+1 : prof[m+1].seqprofA[0]);
					prof[j].typeA = ((i+1 < prof[m+1].seqprofA[0]) ? 0 : 1);
   					prof[j].nA=1;
							
					prof[j].seqprofB = (int *)calloc(1,sizeof(int));
					prof[j].seqprofB[0] = ( (i+1 > prof[m+1].seqprofA[0]) ? i+1 : prof[m+1].seqprofA[0]);
					prof[j].typeB = ((i+1 < prof[m+1].seqprofA[0]) ? 1 : 0);
  		 			prof[j].nB=1;
					
					--j;
					if (j<0)
						break;
					narray[j] = p->back->next;
					narray[m]=NULL;
					
				}
			}
		
		
		/* find nodes clustering with a yet clustered node */
		
		

			if (j<0)
				break;
			
			for (m = spp-3; m >= j ; --m){
				int n;
				if (m < 0)
					break;				
				
				p = narray[m];
				if (p==NULL) {
					continue;
				} else if (p->back==NULL) {
					continue;
				} else if (p->back->next==NULL) {
				
				}
				for (n = spp-3; n >= j ; --n){
					if (n<0)
						break;
					
					if (!(p->back->tip)&& p->back->next->back == narray[n] && narray[n]!=NULL && p->back->next->back!=NULL){
						prof[j].seqprofA = (int *)calloc(1,sizeof(int));
						prof[j].seqprofA[0] = ( (prof[n+1].seqprofA[0] < prof[m+1].seqprofA[0]) ? prof[n+1].seqprofA[0] : prof[m+1].seqprofA[0]);
						prof[j].typeA = 1;
	   					prof[j].nA=1;
						prof[j].seqprofB = (int *)calloc(1,sizeof(int));
						prof[j].seqprofB[0] = ( (prof[n+1].seqprofA[0] > prof[m+1].seqprofA[0]) ? prof[n+1].seqprofA[0] : prof[m+1].seqprofA[0]);
						prof[j].typeB = 1;
   						prof[j].nB=1;
						--j;
						if (j<0)
							break;
						narray[j] = p->back->next->next;
		
						narray[n]=NULL;
						narray[m]=NULL;
						k +=1;
						
	 				}	
	 			}
			}
		}
	}	

	/********* ???????????????? *******************/
	for (i = 0; i < j; ++i)
		free(narray[i]);
	free(narray);
}


/* Function:		countLeafs
 * 
 * Purpose:				
 *           
 * Args:			r:
 *
 * Returns:		int count
 */
int countLeafs(node *r){
	int count=0;
	int odd = 1;
	node *q,*p;
		q = r->back;
	p = r->next;
	while (p!=r){
		if (p->back->tip)
			count +=1,
			p = p->next,
			odd=1;
		else if(odd && !(p->back->tip))
			p = p->back,
			odd=0;
		else
			p = p->next,
			odd=1;
	}
	return count;
}
