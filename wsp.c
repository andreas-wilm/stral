#include <stdio.h>

#include "algebraic.h"
#include "phylip.h"
#include "dist.h"
#include "weighted_nj.h"
#include "wsp.h"
#include "align.h"
#include "nrutils.h"

node ****brancharray;
float ***blgth;
int **arraysize;
float **matl1;
float **matl2;
int **spparray;



/* Function:		mush
 * 
 * Purpose:		calculates the path from root to leaf		
 *           
 * Args:			p,q:
 *				wtree:
 *
 * Returns:		void
 */
void mush(tree wtree, int spp){
	float res; 
	int i,j,k,a;
	int count = 1;
	float **u,**ut;
	float **c,**d,**e;
	brancharray = (node ****)malloc(spp*sizeof(node ***));
	arraysize = (int **)calloc(spp,sizeof(int *));
	blgth = (float ***)calloc(spp,sizeof(float **));
	matl1 = (float **)calloc(spp,sizeof(float *));
	spparray = (int **)calloc(spp,sizeof(int *));
	a = (spp*(spp-1)/2);
	matl2 = (float **)calloc(a+1,sizeof(float *));
	u = (float **)calloc(1+1,sizeof(float*));
	ut = (float **)calloc(a+1,sizeof(float*));
	ut[0] = (float *)calloc(1+1,sizeof(float));
	ut[1] = (float *)calloc(1+1,sizeof(float));
	for (i = 0 ; i <= a ; ++i){
		u[i] = (float *)calloc(a+1,sizeof(float));
		u[i][1] = 1;
		ut[1][i] = 1;
		matl2[i] = (float *)calloc(a+1,sizeof(float ));	
	}
	for (i = 0; i < spp; i++){
		brancharray[i] = (node ***)malloc(spp*sizeof(node **));
		arraysize[i] = (int *)calloc(spp,sizeof(int));
		blgth[i] = (float **)calloc(spp,sizeof(float*));
		matl1[i] = (float *)calloc(spp,sizeof(float ));
		spparray[i] = (int*)calloc(spp,sizeof(int));
	}
	
	
	
	
	for (i = 0; i < spp; i++){
		for (j=i+1; j < spp; j++){
			spparray[i][j]=count;
			spparray[j][i]=spparray[i][j];
			res = leaf2leaf(wtree.nodep[i],wtree,wtree.nodep[j],i,j);
			matl1[i][j]=res;
			printf("path %d->%d: %f\n",i,j,res);
			count++;
		}
	}
	
	
	
	/* a very bad solution */	
	for (i = 0; i < spp; i++){
		for (j = i+1; j < spp; j++){
			int l,o;
			for (o=i; o < spp; ++o){
				for (l=(j > o  ?  j + 1 : o + 1); l < spp; l++){
					float lgth=0;
					for (k= 0; k < arraysize[i][j]; ++k) {
						int n;
						for (n=0; n < arraysize[o][l]; ++n){
							if(brancharray[i][j][k]==brancharray[o][l][n]){
									lgth += blgth[i][j][k];
									/*printf("%p %p",brancharray[i][j][k],brancharray[i][l][n]);*/
									n=100000;
							}
						}
					}
					printf("p %d-%d q %d-%d,length %f\n",i,j,o,l,lgth);
					
					
					matl2[spparray[i][j]][spparray[o][l]] = pow(lgth,2)/((matl1[i][j]==0 ? 1 : matl1[i][j]) *( matl1[o][l]==0 ? 1 :matl1[o][l]));
					printf("p %d-%d q %d-%d,mat %f\n",i,j,o,l,matl2[spparray[i][j]][spparray[o][l]]);
				}
			}	
		}
	}
	
	for (i = 0; i <= a; i++){
		matl2[i][i]=1;
		for (j = i+1; j <= a; j++){
			matl2[j][i]=matl2[i][j];
		}
	}
	
	for (i = 0; i <= a; i++){
		for (j =0; j <= a; j++){
			printf("%f ",matl2[i][j]);
		}
		printf("\n");
	}
	
	matl2 = matinverse(matl2,a);
	
	
	
printf("\n");	
	for (i = 0; i <= a; i++){
		for (j =0; j <= a; j++){
			printf("%f ",matl2[i][j]);
		}
		printf("\n");
	}
	
printf("\n");
	c = matmult(ut,matl2,c,1,a,a);

	for (i = 0; i <= 1; i++){
		for (j =0; j <= a; j++){
			printf("%8.3f ",c[i][j]);
		}
		printf("\n");
	}
printf("\n");	
	d = matmult(matl2,u,d,a,a,1);
	for (i = 0; i <= a; i++){
		for (j =0; j <= 1; j++){
			printf("%8.3f ",d[i][j]);
		}
		printf("\n");
	}

	e = matmult(c,u,e,1,a,1);
printf("\n");
	for (i = 0; i <= 1; i++){
		for (j =0; j <= 1; j++){
			printf("%8.3f ",e[i][j]);
		}
		printf("\n");
	}







}



/* Function:		
 * 
 * Purpose:		calculates the path from root to leaf		
 *           
 * Args:			p,q:
 *				wtree:
 *
 * Returns:		void
 */
float leaf2leaf(node *q,tree wtree, node *p,int k , int l){
	int i=0;
	int m;
	node **nodearray;

	float *branchorder;
	float lengthtoleaf =0;
	int odd=1;
	
	struct node *n;
	n = (struct node *)malloc(sizeof(struct node));
	n->tip=1;
	nodearray = (struct node **)malloc(sizeof(struct node *));
	*(nodearray+i) = (struct node *)malloc(sizeof(struct node));
	
	

	
	
	nodearray[i] = q;
	branchorder = (float *)calloc(1,sizeof(float));
	branchorder[i]=0;
	while(q!=p){
		
		if (q->back == NULL)
			q->back=n;
		if (odd && !(q->back->tip)) {
			int j;
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			branchorder = (float *)realloc(branchorder,(i+2)*(sizeof (float)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = q->back;
			branchorder[i+1]=q->v;
			++i;
/*			printf("nodearray[%d]->next: %p\n",i, nodearray[i]);
			printf("i: %d\n",i);
			fflush(stdout);*/
			odd = !odd;
			q = q->back;
			j = findnode(nodearray, q, i);
/*			printf("j: %d\n",j);*/
			fflush(stdout);	
			if (j < i){
				nodearray = (struct node **)realloc(nodearray,(j+2)*sizeof(struct node *));
				branchorder = (float *)realloc(branchorder,(j+2)*(sizeof (float)));
				nodearray[j+1]=q->next;

				branchorder[j+1]=0;
				j++;
/*				printf("nodearray[%d]: %p\n",j, nodearray[j]);				*/
				i=j;

/*				printf("i: %d, j: %d\n",i,j);*/
				fflush(stdout);
			} 
		} 
		
		
		else if (odd && (q->back->tip) && (q->back!=p)){
			int j;
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			branchorder = (float *)realloc(branchorder,(i+2)*(sizeof(float)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = nodearray[i]->next;
			branchorder[i+1] = 0;
			
			++i;
			q = q->next;
			j = findnode(nodearray, q, i);
/*			printf("nodearray[%d]->next: %p\n",i, nodearray[i]);
			printf("odd,tip: i: %d, j: %d\n",i,j);
			fflush(stdout);*/
			if (j < i){
				nodearray = (struct node **)realloc(nodearray,(j)*sizeof(struct node *));
				branchorder = (float *)realloc(branchorder,(j)*(sizeof (float)));
				i = j;
/*				printf("i: %d, j: %d\n",i,j);
				fflush(stdout);*/
			} 
		} 
		
		
		else if (odd && (q->back->tip) && (q->back==p)){
/*			int j = 0;
			printf("odd,q=p: i: %d, j: %d\n",i,j);*/
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			branchorder = (float *)realloc(branchorder,(i+2)*(sizeof(float)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = nodearray[i]->back;
			branchorder[i+1] = q->v;
			++i;
/*			printf("nodearray[%d]->next: %p\n",i, nodearray[i]);
			printf("odd,q=p: i: %d, j: %d\n",i,j);
			fflush(stdout);*/
			q = q->back;
		} 
		
		
		else if (!odd){
			int j;
			nodearray = (struct node **)realloc(nodearray,(i+2)*sizeof(struct node *));
			branchorder = (float *)realloc(branchorder,(i+2)*(sizeof(float)));
			*(nodearray+(i+1)) = (struct node *)malloc(sizeof(struct node));
			nodearray[i+1] = nodearray[i]->next;
			branchorder[i+1] = 0;
			++i;
/*			printf("nodearray[%d]->next: %p\n",i, nodearray[i]);
			fflush(stdout);*/
			odd = !odd;
			q = q->next;
			j = findnode(nodearray, q, i);
/*			printf("!odd: i: %d, j: %d\n",i,j);
			fflush(stdout);*/
			if (j < i){
				nodearray = (struct node **)realloc(nodearray,(j+2)*sizeof(struct node *));
				branchorder = (float *)realloc(branchorder,(j+2)*(sizeof (float)));
				nodearray[j] = q;
				i = j;
/*				printf("i: %d, j: %d\n",i,j);
				fflush(stdout);*/
			}
		}
	}
	
/*	printf("i: %d\n",i);*/
	brancharray[k][l] = (struct node **)malloc((i+1)*sizeof(struct node *));
	blgth[k][l] = (float *)calloc((i+1),sizeof(float));
	for (m=0; m <= i; ++m){
		lengthtoleaf = lengthtoleaf + branchorder[m];
		arraysize[k][l]=i+1; 
		brancharray[k][l][m]=nodearray[m];
		blgth[k][l][m]=branchorder[m];
		/*printf("nodearray[%d]->v: %f; p: %p\n",k,branchorder[m], nodearray[m]);*/
	}
	
	
	free(nodearray);
	free(branchorder);
	return lengthtoleaf;
}
