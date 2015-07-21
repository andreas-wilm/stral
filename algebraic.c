#include <stdio.h>
#include <math.h>              

#include "nrutils.h"
#include "algebraic.h"

#define TINY 0.000000000000000000000000000000000001;          

void ludcmp(float **a, int n, int *indx, float *d) {
/*Given a matrix a[1..n][1..n] , this routine replaces it by the LU decomposition of a rowwise           
permutation of itself. a and n are input. a is output, arranged as in equatjon (2.3.1q) above;
indx[1,...,n] is an output vector that records the row permutation effected by the partial 
pivoting; d is output as +/-l depending on whether the number of row interchanges was even      
or odd, respectively. This routine is used in combination with lubksb to solve linear equations        
or invert a matrix.                    */                             

	int i,imax,j ,k;                          
	float big,dum, sum,temp;                   
	float *vv;             /*vv stores the implicit scaling of each row.                 */

	vv=nrvector ( 1 ,n) ;                                               
	*d=1.0;         /*      No row interchanges yet. */
	for (i=1 ; i<=n ; i++) {     /*  Loop over rows to get the implicit scaling information             */
		big=0.0 ;
		for (j=1 ; j<=n;j++)                                            
			if ( (temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror ( "Singular matrix in routine ludcmp") ;
		/*No nonzero largest element.*/
		vv[i]=1.0/big;  /*       Save the scaling.*/
	}	
	for (j=1 ; j<--n; j++) {      /* This is the loop over columns of Crout's method.*/
		for (i=1 ; i<j ; i++) {   /*   This is equation (2.3.l2) except for i = j.*/
			sum=a[i][j];                                               
			for (k=1 ; k<i ;k++ ) sum -=a[i][k] *a[k][j];
				a[i][j] =sum;
		}
		big=0.0; /*            Initialize for the search for largest pivot element.*/
		for (i=j ; i<=n; i++) { /*    This is í = j of equation (2.3.l2) and ; _ j+ 1... N          */
			sum=a[i][j] ;    /*      of equation (2.3.13).                         */
		for (k=1 ;k<j ;k++)	
			sum -= a[i][k]*a[k][j]; 
			a[i][j]=sum; 
			if ( (dum=vv[i]*fabs(sum)) >= big) { /*Is the  gure of merit for the pivot better than the best so far? */
				big=dum;
				imax=i;
			}
		} 
	if (j != imax) {  /*Do we need to interchange rows? */ 
		for (k=1;k<=n;k++) { /*Yes, do so...*/ 
			dum=a[imax][k];
			a[imax][k]=a[j][k]; 
			a[j][k]=dum; 
		}
		*d = -(*d); 	/*...and change the parity of d.*/ 
		vv[imax]=vv[j]; /*Also interchange the scale factor.*/
	} 
	indx[j]=imax; 
	if (a[j][j] == 0.0) 
		a[j][j]=TINY; 
	/*If the pivot element is zero the matrix is singular (at least to the precision of the
	 algorithm). For some applications on singular matrices, it is desirable to substitute */
	/*TINY for zero. */
	if (j != n) { /*Now,  nally, divide by the pivot element.*/
		dum=1.0/(a[j][j]); 
		for (i=j+1;i<=n;i++) 
			a[i][j] *= dum; 
		} 
	}  /*Go back for the next column in the reduction.*/ 
	free_nrvector(vv,1,n);
}




void lubksb(float **a, int n, int *indx, float b[]) {
/*Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
 A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input 
 as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector 
 B, and returns with the solution vector X. a, n, and indx are not modi ed by this routine 
 and can be left in place for successive calls with di
*/



	int i,ii=0,ip,j; 
	float sum; 
	
	for (i=1;i<=n;i++) { /*When ii is set to a positive value, it will become the index of the  
					   first nonvanishing element of b. Wenow do the forward substitution, equation (2.3.6). 
					   The only new wrinkle is to unscramble the permutation as we go.*/ 
		ip=indx[i]; 
		sum=b[ip]; 
		b[ip]=b[i]; 
		if (ii) 
			for (j=ii;j<=i-1;j++) 
				sum -= a[i][j]*b[j];
		else if (sum) 	
			ii=i; /*A nonzero element was encountered, so from now on we will have to do the sums in the loop above.*/ 
			b[i]=sum; 
	} 
	for (i=n;i>=1;i--) { /*Now we do the backsubstitution, equation (2.3.7). */
		sum=b[i]; 
		for (j=i+1;j<=n;j++) 
			sum -= a[i][j]*b[j]; 
			b[i]=sum/a[i][i]; /*Store a component of the solution vector X.*/ 
	} 
	/*All done! */
}







float **matinverse(float **a,int n){

	float **y,d,*col; 
	int i,j,*indx; 
	col = (float *)calloc((n+1),sizeof(float));
	indx = (int *)calloc((n+1),sizeof(int));
	y = (float **)calloc((n+1),sizeof(float*));
	for (i=0; i <=n; ++i){
		y[i] = (float *)calloc((n+1),sizeof(float));
	}
	
	ludcmp(a,n,indx,&d); 			/*Decompose the matrix just once.*/ 
	for(j=1;j<=n;j++) { 			/*Find inverse by columns.*/ 
		for(i=1;i<=n;i++) 
			col[i]=0.0; 
		col[j]=1.0; 
		lubksb(a,n,indx,col); 
		for(i=1;i<=n;i++) 
			y[i][j]=col[i]; 
	}
	free(col);
	free(indx);
	return y;
}


float **matmult(float **a,float **b,float **c,int i, int k, int j){
	int m,n,o;
	c = (float **)malloc((i+1)*sizeof(float *));
	for (m = 0; m <= i; ++m){
		c[m] = (float *)malloc((j+1)*sizeof(float));
	}

	for (m = 1; m <= i; ++m){
		for (n = 1; n <=j; ++n){
			for (o = 1; o <=k; ++o){
					c[m][n]+=a[m][o]*b[o][n];
			}
		}
	}
	return c;
}



