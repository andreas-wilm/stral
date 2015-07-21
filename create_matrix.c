#include <stdio.h>
#include <string.h>
#include <math.h>

#include "create_matrix.h"
#include "stral.h"



/* Function:		create_sub_matrix
 * 
 * Purpose:		Build a substitution matrix                        
 *				according to option -R: ribosum or gotoh             
 *				and make it positive (according to                 
 *				Thompson, Higgins & Gibson (1994) CABIOS 10, 19-29
 *           
 * Args:			none
 *
 * Returns:		void
 */
void create_sub_matrix () {
	int i, j; 		   /*values in {A,C,G,U,N,-} */

	memset (d, 0, sizeof (d));

	if (ribosum){
		/*	 Substitution matrix according to  				   *
		 * Klein, R.J. & Eddy, S.R. (2003) BMC Bioinf. 4, 44-59	   *
		 * RSEARCH: Finding homologs of single structured RNA sequences*/

		d[0][0] = +2.22;
		d[0][1] = -1.86;
		d[0][2] = -1.46;
		d[0][3] = -1.39;

		d[1][0] = -1.86;
		d[1][1] = +1.16;
		d[1][2] = -2.48;
		d[1][3] = -1.05;

		d[2][0] = -1.46;
		d[2][1] = -2.48;
		d[2][2] = +1.03;
		d[2][3] = -1.74;

		d[3][0] = -1.39;
		d[3][1] = -1.05;
		d[3][2] = -1.74;
		d[3][3] = +1.65;
		

	} else {
		/* Substitution matrix according to					    * 
		 * Gotoh, O. (1999) Adv. Biophys. 36, 159-206.			    *
		 * Multiple sequence alignment: algorithms and applications     */

		
		d[0][0] = +2;
		d[0][1] = -2;
		d[0][2] = -1;
		d[0][3] = -2;

		d[1][0] = -2;
		d[1][1] = +2;
		d[1][2] = -2;
		d[1][3] = -1;

		d[2][0] = -1;
		d[2][1] = -2;
		d[2][2] = +2;
		d[2][3] = -2;

		d[3][0] = -2;
		d[3][1] = -1;
		d[3][2] = -2;
		d[3][3] = +2;
		
		
	}

		
		
		
		
/**** setting values to positive floats*****/
	if (positive == 1) {
		float min = 100;
		for (i = 0; i < 4; ++i)
			for (j = 0; j < 4; ++j){
				min = ((d[i][j] < min) ? d[i][j] : min);
			}

		for (i = 0; i < 4; ++i)
			for (j = 0; j < 4; ++j){
				d[i][j] = d[i][j] - min;
			}
	}
}


/* Function:		
 * 
 * Purpose:		calculate distance matrix as input for tree construiction
 *           
 * Args:			none
 *
 * Returns:		void
 */
void createDistanceMatrix(float **array, int n){

	float max=-1;
	float min=10.e7;
	int i,j;
	
			
	for (i=0; i < n; ++i){
		for (j=i+1; j < n; ++j){
			max = (max > array[i][j] ? max : array[i][j]);
			min = (min < array[i][j] ? min : array[i][j]);
		}
	}

	
	

	for (i=0; i < n; ++i){
		distmatrix[i][i] = 0;
		for (j=i+1; j < n; ++j){
			if( array[i][j]== max) {
				distmatrix[i][j]=0;	
				distmatrix[j][i]=0;
				continue;
			}
			distmatrix[i][j] = -1 * log((array[i][j]-min+1)/(max-min));
			distmatrix[j][i] = distmatrix[i][j];
		}
	}
	
#ifdef DEBUG	
	if (!verbose) {
		for (i=0; i < n; ++i){
			for (j=0; j < n; ++j){
				printf("%f   ",distmatrix[i][j]);
			}
			printf("\n");
		}
	}
#endif
}
	
