#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "nrutils.h"
#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) > (dminarg2) ? (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) > (minarg2) ? (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) > (lminarg2) ? (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) > (iminarg2) ? (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define ANSI

#if defined(__STRICT_ANSI__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void numrecerror(char error_text[]);
float *nrvector(int nl, int nh);
int *inrvector(long nl, long nh);
unsigned char *cnrvector(long nl, long nh);
unsigned long *lnrvector(long nl, long nh);
double *dnrvector(long nl, long nh);

float **nrmatrix(long nrl, long nrh, long ncl, long nch);
double **dnrmatrix(long nrl, long nrh, long ncl, long nch);
int **inrmatrix(long nrl, long nrh, long ncl, long nch);

float **subnrmatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl);
float **convert_nrmatrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void free_nrvector(float *v, long nl, long nh);
void ifree_nrvector(int *v, long nl, long nh);
void cfree_nrvector(unsigned char *v, long nl, long nh);
void lfree_nrvector(unsigned long *v, long nl, long nh);
void dfree_nrvector(double *v, long nl, long nh);

void free_nrmatrix(float **m, long nrl, long nrh, long ncl, long nch);
void ifree_nrmatrix(int **m, long nrl, long nrh, long ncl, long nch);
void dfree_nrmatrix(double **m, long nrl, long nrh, long ncl, long nch);

void free_subnrmatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_nrmatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);

#else /* not ANSI */

void numrecerror();
float *nrvector();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */

#define NR_END 1
#define FREE_ARG char *



/* Function:		numrecerror
 * 
 * Purpose:		numerical recipes standard error handler 
 *           
 * Args:			error_text:
 *
 * Returns:		void
 */
void numrecerror(char error_text[]){
	fprintf(stderr,"Runtime error...\n");
	fprintf(stderr,"%s\n", error_text);
	fprintf(stderr,"program does not run proper...\n");
	exit(1);
}

/* Function:		nrvector
 * 
 * Purpose:		allocate a float nrvector with subscript range v[nl...nh] 
 *           
 * Args:			nl, nh:
 *
 * Returns:		*nrvector
 */
float *nrvector(int nl, int nh){
	float *v;
	
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) numrecerror("allocation failure in nrvector()");
	return v-nl+NR_END;
}

/* Function:		inrvector
 * 
 * Purpose:		allocate an int nrvector with subscript range v[nl...nh] 
 *           
 * Args:			nl, nh:
 *
 * Returns:		*inrvector
 */
int *inrvector(long nl, long nh){
	int *v;
	
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) numrecerror("allocation failure in inrvector()");
	return v-nl+NR_END;
}



/* Function:		cnrvector
 * 
 * Purpose:		allocate an unsigned char nrvector with subscript range v[nl...nh] 
 *           
 * Args:			nl, nh:
 *
 * Returns:		*cnrvector
 */
unsigned char *cnrvector(long nl, long nh){
	unsigned char *v;
	
	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) numrecerror("allocation failure in cnrvector()");
	return v-nl+NR_END;
}


/* Function:		linrvector
 * 
 * Purpose:		allocate a long nrvector with subscript range v[nl...nh] 
 *           
 * Args:			nl, nh:
 *
 * Returns:		*lnrvector
 */
unsigned long *lnrvector(long nl, long nh){
	unsigned long *v;
	
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned long)));
	if (!v) numrecerror("allocation failure in lnrvector()");
	return v-nl+NR_END;
}


/* Function:		dnrvector
 * 
 * Purpose:		allocate a double nrvector with subscript range v[nl...nh] 
 *           
 * Args:			nl, nh:
 *
 * Returns:		*dnrvector
 */
double *dnrvector(long nl, long nh){
	double *v;
	
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) numrecerror("allocation failure in dnrvector()");
	return v-nl+NR_END;
}

/* Function:		nrmatrix
 * 
 * Purpose:		allocate a float matirx with subscript range m[nrl...nrh][ncl...nch] 
 *           
 * Args:			nrl, nrh , ncl, nch:
 *
 * Returns:		*nrmatrix
 */
float **nrmatrix(long nrl, long nrh, long ncl, long nch){
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) numrecerror("allocation failure 1 in nrmatrix()");
	m += NR_END;
	m -= nrl;
	
	m[nrl]=(float *)malloc((size_t)((nrow+ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) numrecerror("allocation failure 2 in nrmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1; i <=nrh; i++) 
		m[i]=m[i-1]+ncol;
	
	return m;
}



/* Function:		dnrmatrix
 * 
 * Purpose:		allocate a double nrmatrix with subscript range m[nrl...nrh][ncl...nch] 
 *           
 * Args:			nrl, nrh , ncl, nch:
 *
 * Returns:		*dnrmatrix
 */
double **dnrmatrix(long nrl, long nrh, long ncl, long nch){
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	m=(double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) numrecerror("allocation failure 1 in dnrmatrix()");
	m += NR_END;
	m -= nrl;
	
	m[nrl]=(double *)malloc((size_t)((nrow+ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) numrecerror("allocation failure 2 in dnrmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1; i <=nrh; i++) 
		m[i]=m[i-1]+ncol;
	
	return m;
}

/* Function:		inrmatrix
 * 
 * Purpose:		allocate a int nrmatrix with subscript range m[nrl...nrh][ncl...nch] 
 *           
 * Args:			nrl, nrh , ncl, nch:
 *
 * Returns:		*inrmatrix
 */
int **inrmatrix(long nrl, long nrh, long ncl, long nch){
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) numrecerror("allocation failure 1 in inrmatrix()");
	m += NR_END;
	m -= nrl;
	
	m[nrl]=(int *)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) numrecerror("allocation failure 2 in inrmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1; i <=nrh; i++) 
		m[i]=m[i-1]+ncol;
	
	return m;
}


float **subnrmatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl){
	long i,j, nrow=oldrh-oldrl+1,ncol=oldch-newcl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) numrecerror("allocation failure in subnrmatrix()");
	m += NR_END;
	m -= newrl;
	
	for (i=oldrl+1, j =newrl; i <=oldrh; i++, j++) 
		m[j]=a[i]+ncol;
	
	return m;

}
float **convert_nrmatrix(float *a, long nrl, long nrh, long ncl, long nch){
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) numrecerror("allocation failure in subnrmatrix()");
	m += NR_END;
	m -= nrl;
	
	
	m[nrl]=a-ncl;
	for (i=1, j =nrl+1; i <=nrow; i++, j++) 
		m[j]=m[j-1]+ncol;
	
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh){
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1, ndep=ndh-ndl+1;
	float ***t;
	
	t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) numrecerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	
	t[nrl]=(float **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) numrecerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) numrecerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	

	for (j=ncl+1; j <=nch; j++) {
		t[nrl][j]=t[nrl][j-1]+ndep;
	}
	for (i=nrl+1; i <=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for (j=ncl+1; j <=nch;j++) {
			t[i][j]=t[i][j-1]+ndep;		
		}
	}
	return t;
}



/* Function:		free_nrvector
 * 
 * Purpose:		free a float nrvector allocated with nrvector
 *           
 * Args:			v:
 *				nl:
 *				nh:
 *
 * Returns:		void
 */
void free_nrvector(float *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}

/* Function:		free_nrvector
 * 
 * Purpose:		free a float nrvector allocated with nrvector
 *           
 * Args:			v:
 *				nl:
 *				nh:
 *
 * Returns:		void
 */
void ifree_nrvector(int *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}

/* Function:		free_nrvector
 * 
 * Purpose:		free a float nrvector allocated with nrvector
 *           
 * Args:			v:
 *				nl:
 *				nh:
 *
 * Returns:		void
 */
void cfree_nrvector(unsigned char *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}

/* Function:		free_nrvector
 * 
 * Purpose:		free a float nrvector allocated with nrvector
 *           
 * Args:			v:
 *				nl:
 *				nh:
 *
 * Returns:		void
 */
void lfree_nrvector(unsigned long *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}

/* Function:		dfree_nrvector
 * 
 * Purpose:		free a double nrvector allocated with nrvector
 *           
 * Args:			v:
 *				nl:
 *				nh:
 *
 * Returns:		void
 */
void dfree_nrvector(double *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}


/* Function:		free_nrmatrix
 * 
 * Purpose:		free a float nrmatrix allocated by nrmatrix
 *           
 * Args:			m:
 *				nrl,nrh,ncl,nch:
 *
 * Returns:		void
 */ 
void free_nrmatrix(float **m, long nrl, long nrh, long ncl, long nch){
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* Function:		ifree_nrmatrix
 * 
 * Purpose:		free an int nrmatrix allocated by inrmatrix
 *           
 * Args:			m:
 *				nrl,nrh,ncl,nch:
 *
 * Returns:		void
 */
void ifree_nrmatrix(int **m, long nrl, long nrh, long ncl, long nch){
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* Function:		dfree_nrmatrix
 * 
 * Purpose:		free a double nrmatrix allocated by dnrmatrix
 *           
 * Args:			m:
 *				nrl,nrh,ncl,nch:
 *
 * Returns:		void
 */
void dfree_nrmatrix(double **m, long nrl, long nrh, long ncl, long nch){
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* Function:		free_subnrmatrix
 * 
 * Purpose:		free a subnrmatrix allocated by subnrmatrix
 *           
 * Args:			m:
 *				nrl,nrh,ncl,nch:
 *
 * Returns:		void
 */
void free_subnrmatrix(float **b, long nrl, long nrh, long ncl, long nch){
	free((FREE_ARG) (b+nrl-NR_END));
}

/* Function:		free_convert_nrmatrix
 * 
 * Purpose:		free a convert_nrmatrix allocated by convert_nrmatrix
 *           
 * Args:			b:
 *				nrl,nrh,ncl,nch:
 *
 * Returns:		void
 */
void free_convert_nrmatrix(float **b, long nrl, long nrh, long ncl, long nch){
	free((FREE_ARG) (b+nrl-NR_END));
}

/* Function:		free_f3tensor
 * 
 * Purpose:		free a float f3tensor allocated by f3tensor
 *           
 * Args:			
 *
 * Returns:		void
 */
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh){
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
