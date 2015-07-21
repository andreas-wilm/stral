#define _GNU_SOURCE

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h> 
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mhash.h"

#include "stral.h"
#include "align.h"
#include "iofuncs.h"
#include "create_matrix.h"


int pc;
int ac=-1;
int cn=0;
/* Function:		print2probfile
 * 
 * Purpose:		print the RNAfold propabilities to file named like sequence id
 *           
 * Args:			ptr: pointer to first residue to start from the printout
 *
 * Returns:		void
 */
void print2probfile(struct vector *ptr){
	int fd, j, i;
	char buffer[16];
	char c='\n';
	char e=' ';
	unsigned char *hash;
	char *hashfilename;
	
		
	td = mhash_init(MHASH_MD5);
    	if (td == MHASH_FAILED) exit(1);
    	i = 0;	
	while (i < strlen(sdata[ptr->seqID].sequence)) {
		mhash(td, &sdata[ptr->seqID].sequence[i], 1);
		++i;
	}

	hash = mhash_end(td);
	
	hashfilename = (char *)malloc(mhash_get_block_size(MHASH_MD5)*2*sizeof(char));
	
	for (i = 0; i < mhash_get_block_size(MHASH_MD5); i++) {
     	 sprintf(&hashfilename[i*2],"%.2x",hash[i]);
	}
	
    	
	fd = open(hashfilename,O_RDWR|O_CREAT|O_TRUNC, S_IRWXU|S_IRWXG|S_IRWXO);
	j = 0;
	sprintf(buffer, "%s\n",STRAL_VERSION);
	write(fd,&buffer,strlen(buffer));
	while(sdata[ptr->seqID].ID[j]!='\0'){ 
		write(fd,&sdata[ptr->seqID].ID[j],1);
		++j;
	}
	
	sprintf(buffer,"%d",sdata[ptr->seqID].seqlength);
	write(fd,&c,1);
	write(fd, &buffer, strlen(buffer));
	write(fd,&c,1);
	while (ptr!=NULL){
		write(fd,&ptr->base,1);
		write(fd,&e,1);
		sprintf(buffer,"%d",ptr->replbase);
		write(fd, &buffer, strlen(buffer));
		write(fd,&e,1);
		sprintf(buffer,"%f",ptr->p0);
		write(fd, &buffer, strlen(buffer));
		write(fd,&e,1);
		sprintf(buffer,"%f",ptr->p1);
		write(fd, &buffer, strlen(buffer));
		write(fd,&e,1);
		sprintf(buffer,"%f",ptr->p2);
		write(fd, &buffer, strlen(buffer));
		write(fd,&e,1);
		write(fd,&ptr->oribase,1);
		write(fd,&c,1);
		ptr = ptr->previous_ptr;
	}
	close(fd);
}

/* Function:		readfromprobfile
 * 
 * Purpose:		read the RNAfold probabilities from file
 *           
 * Args:			filename: name of the file(the sequence id)
 *				n: programm intern sequence id 
 *
 * Returns:		void
 */
int readfromprobfile(char *filename, int n){
	char buffer[34];
	int seqlen;
	FILE *fz;
	if ((fz = fopen(filename, "r"))==NULL){
		return 0;
	} else {
		int count=0;
		last_ptr[n] = first_ptr[n] = NULL;
		while(fgets(buffer, 34, fz)){
			if (count == 0){
				if (strncmp(buffer,STRAL_VERSION,5)){
					printf("version conflict in file %s\n", filename);
					printf("please remove file from probDIR\n");
					exit(1);				
				}				
			}
			if (count == 2){
				
				sscanf(buffer,"%d",&seqlen);
				if (seqlen!=sdata[n].seqlength){
					printf("input sequence and prob-file sequence differ (in length %d vs %d) %s\n", seqlen,sdata[n].seqlength, filename);
					printf("please remove file from probDIR\n");
					exit(1);				
				}				
			}
			if (count > 2){
				struct vector *new_item_ptr;
				double p0 = -1;
				double p1 = -1;
				double p2 = -1;
				int rb = -1;
				char p[8];
				int i, j;
				
				new_item_ptr = malloc (sizeof (struct vector));
				new_item_ptr->seqID = sdata[n].seqID;
				for (i = 0; i < strlen(buffer); ++i){
					if (i==0){
						new_item_ptr->base = buffer[i];
					} else if (i==2){
						sscanf(&buffer[2],"%d", &rb);
						new_item_ptr->replbase = rb;
					} else if (i==4){
						for (j=0; j<8; ++j){
							p[j]=buffer[i];
							++i;
						}
						sscanf(p,"%lf", &p0);
						new_item_ptr->p0 = p0;
					} else if (i==13){
						for (j=0; j<8; ++j){
							p[j]=buffer[i];
							++i;
						}
						sscanf(p,"%lf", &p1);
						new_item_ptr->p1 = p1;
					} else if (i==22){
						for (j=0; j<8; ++j){
							p[j]=buffer[i];
							++i;
						}
						sscanf(p,"%lf", &p2);
						new_item_ptr->p2 = p2;
					} else if (i==31){
						new_item_ptr->oribase = buffer[i];
					}
				}
				new_item_ptr->next_ptr = first_ptr[n];
				if (last_ptr[n] == NULL){
					last_ptr[n] = new_item_ptr;
					last_ptr[n]->previous_ptr = NULL;
				} else {
					new_item_ptr->next_ptr->previous_ptr = new_item_ptr;
					new_item_ptr->previous_ptr = NULL;
				}
				first_ptr[n] = new_item_ptr;
		
			}
			++count;
		}
		return 1;
	}
}

/* Function:		printdistmatrix
 * 
 * Purpose:		print out the distance matrix to file
 *           
 * Args:			n:
 *
 * Returns:		void
 */
void printdistmatrix(int n){
	int i, fd, j;
	char c='\n';
	char e=' ';
	char *fname;
	char buffer[128];
	fname = strdup("dist.phy");
	fd = open(fname,O_RDWR|O_CREAT|O_TRUNC, S_IRWXU|S_IRWXG|S_IRWXO);
	sprintf(buffer, "%d\n", n);
	write(fd, &buffer, strlen(buffer));
	for (i=0; i < n; ++i){
		j=0;
		sprintf(buffer,"%d",i);
		write(fd,&buffer,strlen(buffer));
		write(fd, &c, 1);
		for (j=0; j < n; ++j){
			sprintf(buffer, "%f",distmatrix[i][j]);
			write(fd, &buffer, strlen(buffer));					
			write(fd, &e, 1);
		}
		write(fd, &c, 1);
	}
	free(fname);
}

/* Function:		print2file
 * 
 * Purpose:		print out the alignment to file
 *           
 * Args:			nseqs: number of sequnces in alignment
 *				filename: name of zhe output file
 *
 * Returns:		void
 */
void print2file(int nseqs, char *filename){
	int i,fd,j;
	struct vector *current_ptr;
	char c='\n';
	char b='>';
	char e=' ';
	int len;
	remgapcols1(nseqs);
	for (len=strlen(filename)-1; len >=0 ; --len){
		if(filename[len]=='/')
			break;
	}
	wname = calloc(strlen(filename)-len+1, sizeof(char));
	for (i=0; i < strlen(filename)-len+1; ++i){
		wname[i]=filename[len+i+1];
	}

	
	fd = open(wname,O_RDWR|O_CREAT|O_TRUNC, S_IRWXU|S_IRWXG|S_IRWXO);

	for (i=0; i < nseqs; ++i){
	write(fd,&b,sizeof(b));
	write(fd,&e,sizeof(e));
	j = 0;
	while(sdata[i].ID[j]!='\0'){ 
		write(fd,&sdata[i].ID[j],1);
		++j;
	}
	write(fd,&c,1);
	current_ptr = last_ptr[i];
	while(current_ptr!=NULL) {
		write(fd,&current_ptr->oribase,1);
		current_ptr = current_ptr->previous_ptr;
		}
	write(fd,&c,1);
	write(fd,&c,1);
	}
	close(fd);
	
}

/* Function:		print2file
 * 
 * Purpose:		print out the alignment to file
 *           
 * Args:			nseqs: number of sequnces in alignment
 *				filename: name of zhe output file
 *
 * Returns:		void
 */
void print2stdout(int nseqs){
	int i;
	struct vector *current_ptr;
	
	for (i=0; i < nseqs; ++i){
		printf("> %s\n",sdata[i].ID);
		current_ptr = last_ptr[i];
		while(current_ptr!=NULL) {
			printf("%c",current_ptr->oribase);
			current_ptr = current_ptr->previous_ptr;
		}
		printf("\n\n");
	}
}


/* Function:		remgaps
 * 
 * Purpose:		remove gaps from one seqeunce
 *           
 * Args:			seq:
 *				len:
 *
 * Returns:		char
 */
char *remgaps(char *seq, int len){
	int i;
	int count=0;
	char *string;

	for (i = 0; i < len ; ++i){
		if (seq[i]=='-')
			++count;
	}
	string = (char *)calloc(len - count+1, sizeof(char));
	
	for (i = 0; i < len ; ++i){
		if (seq[i]!='-'){
			strncat(string, &seq[i],1);	
		}
	}
	seq = realloc(string, (len -count+1) *sizeof(char));
	return seq;
}

/* Function:		remgapcols1
 * 
 * Purpose:		remove gaps from columns only containing gaps in alignment
 *           
 * Args:			n:
 *
 * Returns:		void
 */
void remgapcols1(int n){
	struct vector **ptrarray;
	struct vector *item_ptr;
	int i;
	ptrarray = (struct vector **)malloc(n*sizeof(struct vector *));
	
	
	for (i=0; i < n; ++i){
		ptrarray[i] = last_ptr[i];
		
	}

	
	while(ptrarray[0]!=NULL){
		int boolean =1;
		for (i =0; i < n; ++i){
			if (ptrarray[i]->replbase!=5){
					boolean = 0;
					i = n;
			}
		}
		if (boolean){
			for (i=0; i < n; ++i){
				if (ptrarray[i]==last_ptr[i]){
					last_ptr[i] = last_ptr[i]->previous_ptr;
					last_ptr[i]->next_ptr = NULL;
				} else if (ptrarray[i] != first_ptr[i]){
					item_ptr = ptrarray[i]->next_ptr; 
					item_ptr->previous_ptr = ptrarray[i]->previous_ptr;
					item_ptr->previous_ptr->next_ptr = item_ptr;
				} else {
					first_ptr[i] = ptrarray[i]->next_ptr;
					first_ptr[i]->previous_ptr = NULL;
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


/* Function:		ToACGUN
 * 
 * Purpose:		convert other symbols than "ACGUNacgun" to N's
 *           
 * Args:			seq: the sequence to be checked
 *
 * Returns:		void
 */
void ToACGUN(char *seq){
	for (; *seq != '\0'; seq++){
		if (strchr("ACGUNacgun", *seq) == NULL) *seq = 'N';
	}
}




/* Function:		parsenewicktree
 * 
 * Purpose:		construct guide tree from bionj output
 *           
 * Args:			nseqs: number of sequences in alignment
 *
 * Returns:		void
 */
void parsenewicktree(int nseqs){
	char **nwtree;
	char buffer[16];
	int count=0;
	int i;
	int ***array;
	int cast1;
	FILE *fz;

	nwtree = (char **)calloc(1,sizeof(char *));
	prof = (struct profile *)calloc(nseqs-1,sizeof(struct profile));
	array = (int ***)malloc(20*sizeof(int**));
	array[0] = (int **)malloc(3*sizeof(int*));
	for(i=0; i < 20; ++i){
		int j=0;
		array[i]=(int **)calloc(2,sizeof(int*));	
		for (j=0; j< 2; ++j){
			array[i][j] = (int *)calloc(2,sizeof(int));	
		}
	}
	array[0][2]= (int *)calloc(2,sizeof(int));
	
	
	
	if ((fz = fopen("bionjout", "r"))==NULL){
		/*exit(1);*/
	} 
	
	while(fgets(buffer, 11, fz)){
		if (strncmp(buffer,"\n",1)!=0){
			nwtree = realloc(nwtree,(count+1)*sizeof(char *));
			nwtree[count] = (char *)calloc(8,sizeof(char));
			nwtree[count] = strdup(buffer);
			count++;
		}
	}
	
	
	pc=nseqs-2;
	
	if (strncmp(nwtree[1],"(",1)!=0){
		sscanf(nwtree[1],"%d",&cast1);
		array[0][0][0]=(cast1+1);
		array[0][0][1]=0;
		cn=1;		
		if (strncmp(nwtree[5],"(",1)!=0){
			sscanf(nwtree[5],"%d",&cast1);
			array[0][1][0]=(cast1+1);
			array[0][1][1]=0;
			cn=2;	
		}
	}
	recursion(0,nwtree,count,0,array,0);
	free(array);
}



/* Function:		recursion
 * 
 * Purpose:		construct guide tree from bionj output
 *           
 * Args:			begin:
 *				nwtree:
 *				count:
 *				rn:
 *				array:
 *				p
 *
 * Returns:		void
 */
int recursion(int begin,char **nwtree,int count,int rn, int ***array,int p){
	int i;
	int n;
	int boolean=1;
	int cast1, cast2;

	
	
	for (i = begin+1; i < count; ++i){
		if (strncmp(nwtree[i],"(",1)==0){
			if (rn==0){
				p=pc;			
				/*printf("***************************************\n                 p is %d               				\n****************************************\n",p);*/
			}
			i = recursion(i,nwtree,count,rn+1,array,p)+1;		
			boolean = 0;
		}
		
		/* case 1 "(x,y)" */
		
		if ((strncmp(nwtree[i],")",1)==0) && boolean){
			sscanf(nwtree[i-7],"%d",&cast1);
			sscanf(nwtree[i-3],"%d",&cast2);

			cast1++;
			cast2++;			

			prof[pc].nA=1;
			prof[pc].nB=1;
			prof[pc].typeA=0;
			prof[pc].typeB=0;
			prof[pc].seqprofA = (int*)malloc(sizeof(int));
			prof[pc].seqprofB = (int*)malloc(sizeof(int));
			prof[pc].seqprofA[0] = (cast1 < cast2 ? cast1 : cast2);
			prof[pc].seqprofB[0] = (cast1 > cast2 ? cast1 : cast2);

			

			if (rn>1){
				if(array[rn-1][0][0]==0){
					array[rn-1][0][0]=prof[pc].seqprofA[0];
					array[rn-1][0][1]=1;
				} else {
					array[rn-1][1][0]=prof[pc].seqprofA[0];
					array[rn-1][1][1]=1;
				}
			}

			if (rn==1){
				array[0][cn][0]=prof[pc].seqprofA[0];
				array[0][cn][1]=1;
				array[rn][0][0]  =0;
				array[rn][0][1]  =0;
				array[rn][1][0]  =0;
				array[rn][1][1]  =0;
				cn++;
				
			}
			if (rn==1 && ((i+4) < count-1) && (strncmp(nwtree[i+4],"(",1)!=0)){
				sscanf(nwtree[i+4],"%d",&cast1);
				array[0][cn][0]=cast1+1;
				array[0][cn][1]=0;
				cn++;		
			}
			
			--pc;			
			return i;
		
		/* case 2.1 "((x,y),z)" */
			
		}	else if ((strncmp(nwtree[i],")",1)==0) && (strncmp(nwtree[i-4],",",1)==0) && (rn !=0)){

			cast1 = prof[pc+1].seqprofA[0];
			sscanf(nwtree[i-3],"%d",&cast2);

			cast2++;			

			prof[pc].nA=1;
			prof[pc].nB=1;
			prof[pc].seqprofA = (int*)malloc(sizeof(int));
			prof[pc].seqprofB = (int*)malloc(sizeof(int));
			prof[pc].seqprofA[0] = (cast1 < cast2 ? cast1 : cast2);
			prof[pc].seqprofB[0] = (cast1 > cast2 ? cast1 : cast2);
			prof[pc].typeA=(cast1 < cast2 ? 1 : 0);
			prof[pc].typeB=(cast1 > cast2 ? 1 : 0);

			

			if (rn>1){
				if(array[rn-1][0][0]==0){
					array[rn-1][0][0]=(array[rn][0][0] < cast2 ? array[rn][0][0] : cast2);
					array[rn-1][0][1]=1;
					array[rn][0][0]  =0;
					array[rn][0][1]  =0;					
				} else {
					array[rn-1][1][0]=(array[rn][0][0] < cast2 ? array[rn][0][0] : cast2);
					array[rn-1][1][1]=1;
					array[rn][0][0]  =0;
					array[rn][0][1]  =0;
				}
			}
			
			
			if (rn==1){
				array[0][cn][0]=prof[pc].seqprofA[0];
				array[0][cn][1]=1;
				array[rn][0][0]  =0;
				array[rn][0][1]  =0;
				array[rn][1][0]  =0;
				array[rn][1][1]  =0;
				cn++;
			}
			if (rn==1 && ((i+4) < count-1) && (strncmp(nwtree[i+4],"(",1)!=0)){
				sscanf(nwtree[i+4],"%d",&cast1);
				array[0][cn][0]=cast1+1;
				array[0][cn][1]=0;
				cn++;		
			}
			--pc;
			return i;
			
		/* case 2.2 "(z,(x,y))" */
		
		}	else if ((strncmp(nwtree[i],")",1)==0) && (strncmp(nwtree[begin+4],",",1)==0) && (rn !=0)){

			cast1 = prof[pc+1].seqprofA[0];
			sscanf(nwtree[begin+1],"%d",&cast2);
			
			cast2++;			

			prof[pc].nA=1;
			prof[pc].nB=1;
			prof[pc].seqprofA = (int*)malloc(sizeof(int));
			prof[pc].seqprofB = (int*)malloc(sizeof(int));
			prof[pc].seqprofA[0] = (cast1 < cast2 ? cast1 : cast2);
			prof[pc].seqprofB[0] = (cast1 > cast2 ? cast1 : cast2);
			prof[pc].typeA=(cast1 < cast2 ? 1 : 0);
			prof[pc].typeB=(cast1 > cast2 ? 1 : 0);

			
			if (rn>1){
				if(array[rn-1][0][0]==0){
					array[rn-1][0][0]=(array[rn][0][0] < cast2 ? array[rn][0][0] : cast2);
					array[rn-1][0][1]=1;
					array[rn][0][0]  =0;
					array[rn][0][1]  =0;					
				} else {
					array[rn-1][1][0]=(array[rn][0][0] < cast2 ? array[rn][0][0] : cast2);
					array[rn-1][1][1]=1;
					array[rn][0][0]  =0;
					array[rn][0][1]  =0;
				}
			}
			
			if (rn==1){
				array[0][cn][0]=prof[pc].seqprofA[0];			
				array[0][cn][1]=1;
				array[rn][0][0]  =0;
				array[rn][0][1]  =0;
				array[rn][1][0]  =0;
				array[rn][1][1]  =0;
				cn++;
			}
			if (rn==1 && ((i+4) < count-1) && (strncmp(nwtree[i+4],"(",1)!=0)){
				sscanf(nwtree[i+4],"%d",&cast1);
				array[0][cn][0]=cast1+1;
				array[0][cn][1]=0;
				cn++;		
			}
			
			--pc;
			return i;

		/* case 3 "((w,z),(x,y))" */

		}	else if (strncmp(nwtree[i],")",1)==0 && (rn != 0)){
			int min=prof[pc+1].seqprofB[0];
			bool myvar2 = false;
			int m=0;

			
			cast1=prof[pc+1].seqprofA[0];
			for (n=pc+2; n <= p; ++n){
				if (min!=prof[n].seqprofA[0] || min!=prof[n].seqprofB[0]){
					myvar2 = true;
				}
				if (myvar2){
					cast2=prof[n].seqprofA[0];
					
					for (m = pc+1; m <= p; ++m){
						if (((prof[m].seqprofA[0]==(cast1 < cast2 ? cast1 : cast2)) && (prof[m].seqprofB[0]== (cast1 > cast2 ? cast1 : cast2))) || (cast1==cast2)){
							myvar2 = false;
							continue;
						}
					}
					if (myvar2){
						n=p;
						continue;
					}
				}
			}
			
			
			
			prof[pc].nA=1;
			prof[pc].nB=1;
			prof[pc].seqprofA = (int*)malloc(sizeof(int));
			prof[pc].seqprofB = (int*)malloc(sizeof(int));
			prof[pc].typeA=1;
			prof[pc].typeB=1;
			prof[pc].seqprofA[0] = (cast1 < cast2 ? cast1 : cast2);
			prof[pc].seqprofB[0] = (cast1 > cast2 ? cast1 : cast2);
			
			if (rn>1){
				if(array[rn-1][0][0]==0){
					array[rn-1][0][0]=(array[rn][0][0] < array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
					array[rn-1][0][1]=1;
					prof[pc].seqprofA[0]=(array[rn][0][0] < array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
					prof[pc].seqprofB[0]=(array[rn][0][0] > array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
					array[rn][0][0]  =0;
					array[rn][0][1]  =0;
					array[rn][1][0]  =0;
					array[rn][1][1]  =0;
				} else {
					array[rn-1][1][0]=(array[rn][0][0] < array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
					array[rn-1][1][1]=1;
					prof[pc].seqprofA[0]=(array[rn][0][0] < array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
					prof[pc].seqprofB[0]=(array[rn][0][0] > array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
					array[rn][0][0]  =0;
					array[rn][0][1]  =0;
					array[rn][1][0]  =0;
					array[rn][1][1]  =0;
	
				}
			}
						
			
			if (rn==1){
				array[0][cn][0]=( array[1][0][0] < array[1][1][0] ? array[1][0][0] : array[1][1][0]); 
				array[0][cn][1]=1;
				prof[pc].seqprofA[0]=(array[rn][0][0] < array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
				prof[pc].seqprofB[0]=(array[rn][0][0] > array[rn][1][0] ? array[rn][0][0] : array[rn][1][0]);
				array[rn][0][0]  =0;
				array[rn][0][1]  =0;
				array[rn][1][0]  =0;
				array[rn][1][1]  =0;
				cn++;
			}
			if (rn==1 && ((i+4) < count-1) && (strncmp(nwtree[i+4],"(",1)!=0)){
				sscanf(nwtree[i+4],"%d",&cast1);
				array[0][cn][0]=cast1+1;
				array[0][cn][1]=0;
				cn++;		
			}
			--pc;
			return i;

		/* case 4 "(A,B,C)" */
		
		} 	else if (strncmp(nwtree[i],")",1)==0 && (rn == 0)){
			
			if ((strncmp(nwtree[i-4],",",1)==0) && array[0][2][0]==0){
				sscanf(nwtree[i-3],"%d",&cast2);
				array[0][2][0]=cast2+1;		
				array[0][2][1]=0;			
			
			}
			
			prof[pc].nA=1;
			prof[pc].nB=1;
			prof[pc].seqprofA = (int*)malloc(sizeof(int));
			prof[pc].seqprofB = (int*)malloc(sizeof(int));
			prof[pc].typeA = (array[0][0][0] < array[0][2][0] ? array[0][0][1] : array[0][2][1]);
			prof[pc].typeB = (array[0][0][0] > array[0][2][0] ? array[0][0][1] : array[0][2][1]);
			prof[pc].seqprofA[0] = (array[0][0][0] < array[0][2][0] ? array[0][0][0] : array[0][2][0]);
			prof[pc].seqprofB[0] = (array[0][0][0] > array[0][2][0] ? array[0][0][0] : array[0][2][0]);
			--pc;
			prof[pc].nA=1;
			prof[pc].nB=1;
			prof[pc].seqprofA = (int*)malloc(sizeof(int));
			prof[pc].seqprofB = (int*)malloc(sizeof(int));
			prof[pc].typeA=( prof[pc+1].seqprofA[0] < array[0][1][0] ? prof[pc+1].typeA : array[0][1][1]);
			prof[pc].typeB=( prof[pc+1].seqprofA[0] > array[0][1][0] ? prof[pc+1].typeA : array[0][1][1]);
			prof[pc].seqprofA[0] = ( prof[pc+1].seqprofA[0] < array[0][1][0] ? prof[pc+1].seqprofA[0] : array[0][1][0]);
			prof[pc].seqprofB[0] = ( prof[pc+1].seqprofA[0] > array[0][1][0] ? prof[pc+1].seqprofA[0] : array[0][1][0]);
			return i;
		}
	}
	return count;
}



/* Function:		
 * 
 * Purpose:		
 *           
 * Args:			
 *
 * Returns:		void
 */

void printpairwise(int seq1,int seq2, char *seqA, char *seqB, int boolean, char *filename){
	int j,fd, len;
	char c='\n';
	char b='>';
	char e=' ';
	char *wname;

	len = strlen(sdata[seq1].ID)+strlen("-vs-");
	
	if (boolean){
		int i;
		for (len=strlen(filename)-1; len >=0 ; --len){
		if(filename[len]=='/')
			break;
		}
		wname = calloc(strlen(filename)-len+1, sizeof(char));
		for (i=0; i < strlen(filename)-len+1; ++i){
			wname[i]=filename[len+i+1];
		}		
	}
	else {
		wname = (char*)calloc((strlen(seqA)+strlen(seqB)+5),sizeof(char));
		for (j=0; j<strlen(sdata[seq1].ID);++j){
			wname[j] = sdata[seq1].ID[j];
		}
		strncat(wname,"-vs-",4);
		for(j=0; j < strlen(sdata[seq2].ID);++j){
			wname[j+len] = sdata[seq2].ID[j];
		}
		wname[j+len]='\0';
	}
	fd = open(wname,O_RDWR|O_CREAT|O_TRUNC, S_IRWXU|S_IRWXG|S_IRWXO);
	free(wname);
	write(fd,&b,sizeof(b));
	write(fd,&e,sizeof(e));
	j = 0;
	while(sdata[seq1].ID[j]!='\0'){ 
		write(fd,&sdata[seq1].ID[j],1);
		++j;
	}
	write(fd,&c,1);

	j=strlen(seqA)-1;
	while(seqA[j]!='\0'){ 
		write(fd,&seqA[j],1);
		--j;
	}
	
	write(fd,&c,1);
	write(fd,&c,1);
	
	write(fd,&b,sizeof(b));
	write(fd,&e,sizeof(e));
	j = 0;
	while(sdata[seq2].ID[j]!='\0'){ 
		write(fd,&sdata[seq2].ID[j],1);
		++j;
	}
	write(fd,&c,1);
	j=strlen(seqB)-1;
	while(seqB[j]!='\0'){ 
		write(fd,&seqB[j],1);
		--j;
	}

	write(fd,&c,1);
	write(fd,&c,1);
	
	close(fd);
	
}




void printstart2end(int begin, int end, struct vector *start_ptr){
	int i;
	for (i = 0; i < begin-1; i++){
		start_ptr=start_ptr->previous_ptr; 
	}
	for (i=begin; i <= end; ++i){
		fprintf(stdout,"%c",start_ptr->base);
		start_ptr=start_ptr->previous_ptr;
	}
	printf("\n");
}

void readsubmatrix(char *filename){
	char buffer[64];
	char *tab, *space, *newline;
	int i,j;
	char *ptr;
	tab="\t";
	space=" ";
	newline="\n";
	FILE *fz;
	fprintf(stdout,"%s\n",filename);
	if ((fz = fopen(filename, "r"))==NULL){
		fprintf(stdout,"unable to open file for some reason\n");
		exit(1); 	
		/*return 0;*/
	} else {
		memset (d, 0, sizeof (d));
		i = 0; 
		j = 0;
		while(fgets(buffer, 64, fz)){
			int len;
			len = strlen(buffer);
			ptr = strtok(buffer, "\n\t ");
			while(ptr != NULL) {
     			if (i % 4 == 0 && i !=0){
					j++;
				}
				sscanf(ptr,"%f",&d[j][i%4]);
				printf("% d. Wort: %s\n",i++,ptr);
				
				
				ptr = strtok(NULL, "\n\t ");
			}
		}
	}
	for (i = 0; i < 4; i++){
		for (j = 0; j < 4; j++){
			fprintf(stdout,"d[%d][%d]:%f\n",i,j,d[i][j]);
		}
	}
	
}	


/* Function:		printtreefile
 * 
 * Purpose:		print the clusters to file
 *           
 * Args:			
 *
 * Returns:		void
 */
void printtreefile(struct profile *prof, int nseqs){
	int fd, j, i;
	char buffer[16];
	char *filename="tree.dat";
		
	fd = open(filename,O_RDWR|O_CREAT|O_TRUNC, S_IRWXU|S_IRWXG|S_IRWXO);
	for (i=nseqs-2; i >=0;--i){
		sprintf(buffer,"Cycle %d %d %d",i,prof[i].pnoA,prof[i].pnoB);
		write(fd, &buffer, strlen(buffer));
		sprintf(buffer,"\n");
		write(fd, &buffer, strlen(buffer));
		sprintf(buffer,"SizeA %d",prof[i].nA);
		write(fd, &buffer, strlen(buffer));
		for (j = 0; j < prof[i].nA; j++){
			sprintf(buffer," %d",prof[i].seqprofA[j]);
			write(fd, &buffer, strlen(buffer));
		}
		sprintf(buffer,"\n");
		write(fd, &buffer, strlen(buffer));
		sprintf(buffer,"SizeB %d",prof[i].nB);
		write(fd, &buffer, strlen(buffer));
		for (j = 0; j < prof[i].nB; j++){
			sprintf(buffer," %d",prof[i].seqprofB[j]);
			write(fd, &buffer, strlen(buffer));
		}
		sprintf(buffer,"\n");
		write(fd, &buffer, strlen(buffer));
	}
	exit(0);
	
	close(fd);
}
