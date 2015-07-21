#define _GNU_SOURCE

#include <dirent.h>
#include "config.h"

#include <stdio.h>
#include <time.h>
#include <ctype.h>

#if HAVE_STRING_H
    #include <string.h>
    #if !HAVE_MEMSET
        #define memset(s, c, n) bzero((s), (n))
    #else
    #endif
    #if !HAVE_STRCHR
        #define strchr index
    #else
    #endif
    #if !HAVE_STRRCHR
    #else
    #endif
    #if !HAVE_STRSTR
        #define strrchr rindex
    #else
    #endif
    #if HAVE_STRDUP
    #else
    #endif
#else
#endif
#if HAVE_STDLIB_H
    #include <stdlib.h>
#else
#endif
#if HAVE_LIMITS_H
    #include <limits.h>
#else
#endif
#if HAVE_UNISTD_H
    #include <unistd.h>
    #if HAVE_GETCWD
    #else
    #endif
#else
#endif
#if HAVE_FCNTL_H
    #include <fcntl.h>
#else
#endif
/* mkdir HAVE_MKDIR */
#if HAVE_SYS_TYPES_H
    #include <sys/types.h>
#else
#endif
#if HAVE_SYS_STAT_H
    #include <sys/stat.h>
#else
#endif
#if HAVE_DIRENT_H
    #include <dirent.h>
#else
#endif


/* STRAL */
#include "stral.h"
#include "align.h"
#include "iofuncs.h"
#include "iterate.h"
#include "create_matrix.h"
#include "config.h"


/* MHASH */
#include "mhash.h"

/* SQUID */
#include "squid.h"
#include "msa.h"


/* RNAfold */
#include "part_func.h"

/* NEIGHBOR */
#include "neighbor.h"

/* BIONJ */
#include "bionj.h"




char *seq;

float alpha=5;
float gapOpenM= 6;
float gapExtM = 1;
float gapOpen = 6;
float gapExt  = 1;
int   njoin   = 0;
int   positive= 1;
int   ribosum = 1;
int   weighted= 0;
int   rooted  = 0;
int   verbose = 0;
int   ppair   = 0;
int   bion    = 0;
int   weigh   = 0;
int   adjust  = 0;
int   upgma   = 0;
int   nofile  = 0;
int   gapmodel= 0;
int   noterm  = 0;
int   palign  = 0;
int 	 fft     = 0;
int   spread  = 0;
int   submat  = 0;
float msascore= 0;







void usagestral(void);
void weighbor(FILE *,FILE*);
float MsaPwIdent(MSA *msa);
void readFile(SQFILE *dbfp, char *seq, SQINFO sqinfo,int type, int nsequences);
/*
static char usage[]  = "\
        Usage: stral [-options] <seqfile>\n\
        Available options:\n\
        -A    : report per-sequence info, not just a summary\n\
        -G    : help; display usage and version\n\
";  

static char experts[] = "\
        --gccomp       : with -a, include GC composition in report (DNA/RNA only)\n\
        --informat <s> : specify sequence file format <s>\n\
        --quiet        : suppress verbose header (used in regression testing)\n\
";
*/
struct opt_s OPTIONS[] = {
        { "-a", TRUE, sqdARG_NONE },    
        { "-h", TRUE, sqdARG_NONE },    
        { "--gccomp",   FALSE, sqdARG_NONE },
        { "--informat", FALSE, sqdARG_STRING },
        { "--quiet",    FALSE, sqdARG_NONE   },
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))




int main(int argc, char **argv){
	char		*seqfile,*seqfile1,*outfile;			/* name of sequence file     */
	SQFILE   *dbfp, *dbfp1;			/* open sequence file        */
	int       fmt;				/* format of seqfile         */
	char     *seq;				/* sequence                  */
	SQINFO    sqinfo;			/* extra info about sequence */
	char     *afile;			/* name of aligned sequence file */
	MSAFILE  *afp;				/* pointer to open alignment file*/
	MSA      *msa;				/* multiple sequence alignment   */  
	int       type,i;			/* kAmino, kDNA, kRNA, or kOtherSeq */
	char     *matrixname;
        
  
  /*
	char		*optname;
	char		*optarg;
	int		optind;
        */
	int		nseqs=0;
	int		nsequences=0;		/*cmp  nseqs*/
	int		nsequences1=0;
	mode_t 	mode = S_IRWXU;
	float	oldmsascore=0;        
	int		iterate = 0;
	int		bfiterate = 0;
	struct	vector *curr_ptr, *old_ptr;
	int		i1,x,nsample=1000;
	float	sum =0.0;
	int		len;
	int 		filebool = 1;
     int		onlytree = 0;
	   /***********************************************
        * Parse command line
        ***********************************************/
	double	temperature = 37.;
	
        
	seqfile1 = NULL;
	outfile = NULL;
	fmt       = SQFILE_UNKNOWN;        /* default: autodetect format  */
	type      = kOtherSeq;                /* just to silence gcc uninit warning */

	if (argc<2)
		usagestral();

	for(i=1;i<argc;i++){
		if(argv[i][0]=='-'){
			if(!(strcmp(&argv[i][1],"-no-freeendgaps"))){
				noterm=0;
			} else if(!(strcmp(&argv[i][1],"-palign"))){
				palign=1;
			} else if(!(strcmp(&argv[i][1],"-groups"))){
				char *g1,*g2;
				i++;
				printf("%d %s\n",i,argv[i]);
				seqfile = argv[i];
				i++;
				printf("%d %s\n",i,argv[i]);
				seqfile1 = argv[i];
				printf("%d %s\n",i,g2);
				filebool = 0;
				/*if(i==argc-1)
					usagestral();
				if(sscanf(argv[++i],"%d",&bfiterate))*/
					;
			} else if(!(strcmp(&argv[i][1],"-fft"))){
				fft=1;
			} else if(!(strcmp(&argv[i][1],"-spread"))){
				spread=1;
			} else if(!(strcmp(&argv[i][1],"-submatrix-file"))){
				matrixname = argv[++i];
				submat=1;
			} else if(!(strcmp(&argv[i][1],"-tree"))){
				onlytree=1;
			} else {
				switch(argv[i][1]){
					case'N':
						njoin=1;
						break;
					case'W':
						weighted=1;
						njoin=1;
						break;
					case'A':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%f",&alpha))
							;		
						else
							usagestral();
						break;
					case'G':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%f",&gapOpen))
							;
						else
							usagestral();
						break;
					case'E':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%f",&gapExt))
							;
						else
							usagestral();
						break;
					case'M':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%f",&gapExtM))
							;
						else
							usagestral();
						break;
					case'O':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%f",&gapOpenM))
							;
						else
							usagestral();
						break;
					case'R':
						ribosum=!ribosum;
						break;	
					case'P':
						positive=!positive;
						break;
					case'V':
						verbose=1;
						break;
					case'r':
						rooted=1;
						njoin=1;
						weighted=1;
						break;
					case'T':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%lf",&temperature))
							;
						else
							usagestral();
						break;
					case'I':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%d",&iterate))
							;
						else
							usagestral();
						break;
					case'B':
						if(argv[i][2]!='\0')
							usagestral();
						if(i==argc-1)
							usagestral();
						if(sscanf(argv[++i],"%d",&bfiterate))
							;
						else
							usagestral();
						break;	
					case'p':
						ppair=!ppair;	
						break;
					case'u':
						njoin=1;
						weighted=1;
						upgma=1;
						break;
					case'b':
						bion=!bion;
						break;
					case'w':
						weigh=!weigh;
						break;
					case'a':
						adjust=!adjust;
						break;
						case'n':
						nofile=!nofile;
						break;
					case'g':
						gapmodel=!gapmodel;
						break;
					case'f':
						if(argv[i][2]!='\0')
							usagestral();
						seqfile = argv[++i];
						filebool = 0;
						break;
					case'h':
						if(argv[i][2]!='\0')
							usagestral();
						seqfile1 = argv[++i];
						filebool = 0;
						break;
					case'i': /* TODO ändern */
						if(argv[i][2]!='\0')
							usagestral();
						outfile = argv[++i];
						break;
					default:
						usagestral();
					}
				}
			}
		}

 
	if (filebool){
		usagestral();	
	}

        if (!njoin && weighted){
                printf("-N/-W resp. -N/-r combination impossible option\n");
                usagestral();
        }
        if (nofile && ppair){
                printf("-p/-n  combination impossible option\n");
                usagestral();
        }        
        if (!iterate && palign){
                printf("please give number of iteration cycles (option -I)\n");
                usagestral();
        }        
        
       
  

        /***********************************************
        * Read the file.
        ***********************************************/
	


	if (seqfile1!=NULL){
		struct profile profile;
		int i;
		
		getcwd(workingdir,200);
		
		dbfp = SeqfileOpen(seqfile, fmt, NULL);
		while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo)){
			nsequences++;
		}
		SeqfileClose(dbfp);
		printf("count %d\n", nsequences);
		
		dbfp = SeqfileOpen(seqfile1, fmt, NULL);
		while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo)){
			nsequences1++;
		}
		SeqfileClose(dbfp);
		printf("count %d\n", nsequences+nsequences1);
		
		sdata =  (struct seqdata *)calloc(nsequences+nsequences1,sizeof(struct seqdata));
		first_ptr = (struct vector **)calloc(nsequences+nsequences1, sizeof(struct vector *));
		last_ptr  = (struct vector **)calloc(nsequences+nsequences1, sizeof(struct vector *));
		
		dbfp = SeqfileOpen(seqfile, fmt, NULL);
		dbfp1 = SeqfileOpen(seqfile1, fmt, NULL);
		
		if (opendir("./probDIR")==NULL){
			mkdir("./probDIR",mode | S_IRGRP | S_IXGRP| S_IXOTH | S_IROTH);
		}
		chdir("./probDIR");

		
		while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo)){
			struct vector *iptr;
			int counter,basecounter;
			char hashfilename[32];
     	     unsigned char	*hash;      
         		char		*structure=NULL;	/*parameter to pf_fold*/      
			
			if (nseqs == 0) 
				type = Seqtype(seq);
			sdata[nseqs].format=dbfp->format;
			sdata[nseqs].ID = strdup(sqinfo.name);
			sdata[nseqs].seqID=nseqs;
			sdata[nseqs].sequence = malloc(sqinfo.len*sizeof(char));
			sdata[nseqs].orisequence = malloc(sqinfo.len*sizeof(char));
			sdata[nseqs].orisequence=strdup(seq);
			seq = remgaps(seq,sqinfo.len);

			ToIUPAC(seq,TRUE);
			ToRNA(seq);
			s2upper(seq);
			ToACGUN(seq);
			sdata[nseqs].sequence =strdup(seq); 
			sdata[nseqs].seqlength=strlen(sdata[nseqs].sequence);
			sdata[nseqs].weight = 1;

			td = mhash_init(MHASH_MD5);
			if (td == MHASH_FAILED) exit(1);
			i = 0;        
			while (i < strlen(sdata[nseqs].sequence)) {
				mhash(td, &sdata[nseqs].sequence[i], 1);
				++i;
			}

			hash = mhash_end(td);

			/*hashfilename = (char *)malloc(mhash_get_block_size(MHASH_MD5)*2*sizeof(char));*/
			for (i = 0; i < mhash_get_block_size(MHASH_MD5); i++) {
				sprintf(&hashfilename[i*2],"%.2x",hash[i]);
			}
                
			if (readfromprobfile(hashfilename, nseqs)){        
				print_data_content2(last_ptr[nseqs]);
#ifdef DEBUG                    
				if (!verbose){
					printf("\n");
					print_data_content1(last_ptr[nseqs]);
					printf("\n");
				}
#endif
			} else {
				pf_fold(sdata[nseqs].sequence, structure,sdata[nseqs].seqID, nsequences);
				print2probfile(last_ptr[nseqs]);
			}
			

			
			iptr = last_ptr[nseqs];
			printf("Sequence: %s\n",sdata[nseqs].orisequence);
			printf("Sequence length: %d\n",sqinfo.len);
			printf("Sequence length: %d\n",sdata[nseqs].seqlength);
			basecounter=0;
			for (counter = 0; counter < sqinfo.len; counter++) {
				printf("Count %d  ",counter);
				printf("baseCount %d\n",basecounter);
				printf("item_ptr %p %c\n",iptr,iptr->base);
				
				/*printf("%c, %d\n",sdata[nseqs].orisequence[counter],counter);*/
				if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter!=0 && basecounter==sdata[nseqs].seqlength){
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->previous_ptr = NULL;
					/*iptr->previous_ptr->next_ptr = new_item_ptr;*/
					iptr->previous_ptr = new_item_ptr;
					new_item_ptr->next_ptr = iptr;
					first_ptr[nseqs]=new_item_ptr;
					iptr=new_item_ptr;
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;
					/*iptr = iptr->previous_ptr;*/
/*					gap einfügen*/
				} else if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter!=0 && counter!=(sqinfo.len-1)){
					/*printf("%c",sdata[nseqs].orisequence[counter]);*/
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->previous_ptr = iptr->previous_ptr;
					iptr->previous_ptr->next_ptr = new_item_ptr;
					iptr->previous_ptr = new_item_ptr;
					new_item_ptr->next_ptr = iptr;
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;
					iptr = iptr->previous_ptr;
/*					gap einfügen*/
				} else if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter==0){
					/*printf("%c",sdata[nseqs].orisequence[counter]);*/
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->previous_ptr = last_ptr[nseqs];
					new_item_ptr->next_ptr = NULL;
					last_ptr[nseqs]->next_ptr= new_item_ptr;
					last_ptr[nseqs] = new_item_ptr;					
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;
					iptr = new_item_ptr;
				} else if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter==sqinfo.len-1){
					/*printf("%c",sdata[nseqs].orisequence[counter]);*/
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->next_ptr = first_ptr[nseqs];
					new_item_ptr->previous_ptr = NULL;
					first_ptr[nseqs]->previous_ptr = new_item_ptr;
					first_ptr[nseqs] = new_item_ptr;
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;		
				} else if (counter==0){
					/*do nothing*/				
				} else {
					/*printf("%c\n",sdata[nseqs].orisequence[counter]);*/
					fflush(stdout);
					if (iptr->previous_ptr!=NULL)
						iptr = iptr->previous_ptr;
					basecounter+=1;	
				}
			}
			print_data_content2(last_ptr[nseqs]);
			printf("\nSeqlength: %d\n",sdata[nseqs].seqlength);
			printf("\n");
			nseqs++;
		}
		SeqfileClose(dbfp);
	

		while (ReadSeq(dbfp1, dbfp1->format, &seq, &sqinfo)){
			struct vector *iptr;
			int counter,basecounter;
			char hashfilename[32];
     	     unsigned char	*hash;      
         		char		*structure=NULL;	/*parameter to pf_fold*/      
			
			if (nseqs == 0) 
				type = Seqtype(seq);
			sdata[nseqs].format=dbfp1->format;
			sdata[nseqs].ID = strdup(sqinfo.name);
			sdata[nseqs].seqID=nseqs;
			sdata[nseqs].sequence = malloc(sqinfo.len*sizeof(char));
			sdata[nseqs].orisequence = malloc(sqinfo.len*sizeof(char));
			sdata[nseqs].orisequence=strdup(seq);
			seq = remgaps(seq,sqinfo.len);

			ToIUPAC(seq,TRUE);
			ToRNA(seq);
			s2upper(seq);
			ToACGUN(seq);
			sdata[nseqs].sequence =strdup(seq); 
			sdata[nseqs].seqlength=strlen(sdata[nseqs].sequence);
			sdata[nseqs].weight = 1;

			td = mhash_init(MHASH_MD5);
			if (td == MHASH_FAILED) exit(1);
			i = 0;        
			while (i < strlen(sdata[nseqs].sequence)) {
				mhash(td, &sdata[nseqs].sequence[i], 1);
				++i;
			}

			hash = mhash_end(td);

			/*hashfilename = (char *)malloc(mhash_get_block_size(MHASH_MD5)*2*sizeof(char));*/
			for (i = 0; i < mhash_get_block_size(MHASH_MD5); i++) {
				sprintf(&hashfilename[i*2],"%.2x",hash[i]);
			}
                
			if (readfromprobfile(hashfilename, nseqs)){        
				print_data_content2(last_ptr[nseqs]);
#ifdef DEBUG                    
			if (!verbose){
				printf("\n");
				print_data_content1(last_ptr[nseqs]);
				printf("\n");
			}
#endif
			} else {
				pf_fold(sdata[nseqs].sequence, structure,sdata[nseqs].seqID, nsequences);
				print2probfile(last_ptr[nseqs]);
			}
			
			iptr = last_ptr[nseqs];
			printf("Sequence: %s\n",sdata[nseqs].orisequence);
			printf("Sequence length: %d\n",sqinfo.len);
			printf("Sequence length: %d\n",sdata[nseqs].seqlength);
			basecounter = 0;
			for (counter = 0; counter < sqinfo.len; counter++) {
				/*printf("%c, %d\n",sdata[nseqs].orisequence[counter],counter);*/
				printf("Count %d  ",counter);
				printf("baseCount %d\n",basecounter);
				fflush(stdout);
				if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter!=0 && basecounter==sdata[nseqs].seqlength){
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->previous_ptr = NULL;
					/*iptr->previous_ptr->next_ptr = new_item_ptr;*/
					iptr->previous_ptr = new_item_ptr;
					new_item_ptr->next_ptr = iptr;
					first_ptr[nseqs]=new_item_ptr;
					iptr=new_item_ptr;
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;
					/*iptr = iptr->previous_ptr;*/
/*					gap einfügen*/
				} else if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter!=0 && counter!=(sqinfo.len-1)){
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->previous_ptr = iptr->previous_ptr;
					iptr->previous_ptr->next_ptr = new_item_ptr;
					iptr->previous_ptr = new_item_ptr;
					new_item_ptr->next_ptr = iptr;
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;
					iptr = iptr->previous_ptr;
/*					gap einfügen*/
				} else if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter==0){
					/*printf("%c",sdata[nseqs].orisequence[counter]);*/
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->previous_ptr = last_ptr[nseqs];
					new_item_ptr->next_ptr = NULL;
					last_ptr[nseqs]->next_ptr= new_item_ptr;
					last_ptr[nseqs] = new_item_ptr;					
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;
					iptr = new_item_ptr;
				} else if (strncmp(&sdata[nseqs].orisequence[counter],"-",1)==0 && counter==sqinfo.len-1){
					/*printf("%c",sdata[nseqs].orisequence[counter]);*/
					struct vector *new_item_ptr;
					new_item_ptr = malloc(sizeof (struct vector));
					new_item_ptr->seqID=iptr->seqID;
					new_item_ptr->next_ptr = first_ptr[nseqs];
					new_item_ptr->previous_ptr = NULL;
					first_ptr[nseqs]->previous_ptr = new_item_ptr;
					first_ptr[nseqs] = new_item_ptr;
					new_item_ptr->base='-';
					new_item_ptr->oribase='-';
					new_item_ptr->replbase=5;
					new_item_ptr->p0=0;
					new_item_ptr->p1=0;
					new_item_ptr->p2=0;		
				} else if (counter==0){
					/*do nothing*/				
				} else {
					/*printf("%c\n",sdata[nseqs].orisequence[counter]);*/
					fflush(stdout);
					if (iptr->previous_ptr!=NULL)
						iptr = iptr->previous_ptr;
					basecounter +=1;
				}
			}
			print_data_content2(last_ptr[nseqs]);
			printf("\nSeqlength: %d\n",sdata[nseqs].seqlength);
			printf("\n");
			nseqs++;
		}
	
		SeqfileClose(dbfp1);

	
		/* initialize the substitution matrix */
		if (submat){
			readsubmatrix(matrixname);
		}
		create_sub_matrix();
		
		/* profile erstellen */
		/* storing clustering information from guide trees*/
		
		profile.seqprofA = (int*)malloc(nsequences*sizeof(int));
		profile.seqprofB = (int*)malloc(nsequences1*sizeof(int));
		profile.nA = nsequences;
		profile.nB = nsequences1;
		profile.typeA = 1;
		profile.typeB = 1;
		for (i = 0; i < nsequences; i++){
			profile.seqprofA[i]=sdata[i].seqID+1;	
		}
		for (i = 0; i < nsequences1; i++){
			profile.seqprofB[i]=sdata[i+nsequences].seqID+1;	
		}
		methodmult(profile, sdata[profile.seqprofA[0]-1].seqlength, sdata[profile.seqprofB[0]-1].seqlength, nseqs, 0);
		
		/* in datei schreiben und nen schicken namen */
		if (!nofile){
			chdir(workingdir);
			if (opendir("./resultDIR")==NULL){
				mkdir("./resultDIR",mode| S_IRGRP | S_IXGRP| S_IXOTH | S_IROTH);
			}
			chdir("./resultDIR");
		}
		if (outfile==NULL){
			outfile="merge1.vie";
		}
		print2file(nseqs, outfile);
		exit(0);
		
	}




        if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
        Die("Failed to open sequence file %s for reading", seqfile);
        
        while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo)){
                nsequences++;
        }
        
        SeqfileClose(dbfp);

	  
        dbfp = SeqfileOpen(seqfile, fmt, NULL);
        
        sdata =  (struct seqdata *)calloc(nsequences,sizeof(struct seqdata));

        first_ptr = (struct vector **)calloc(nsequences, sizeof(struct vector *));
        last_ptr  = (struct vector **)calloc(nsequences, sizeof(struct vector *));

        
        getcwd(workingdir,200);
        
        
        /* get name of input file*/
        for (len=strlen(seqfile)-1; len >=0 ; --len){
                if(seqfile[len]=='/')
                        break;
        }
        afile = calloc(strlen(seqfile)-len+1, sizeof(char));
        for (i=0; i < strlen(seqfile)-len+1; ++i){
                afile[i]=seqfile[len+i+1];
        }






        if (opendir("./probDIR")==NULL){
                mkdir("./probDIR",mode | S_IRGRP | S_IXGRP| S_IXOTH | S_IROTH);
        }
        chdir("./probDIR");


        
     while (ReadSeq(dbfp, dbfp->format, &seq, &sqinfo)){
		int i;
		char hashfilename[32];
          unsigned char	*hash;      
         	char		*structure=NULL;	/*parameter to pf_fold*/      
                
                if (nseqs == 0) 
                        type = Seqtype(seq);
        
                

                sdata[nseqs].format=dbfp->format;
                sdata[nseqs].ID = strdup(sqinfo.name);
                sdata[nseqs].seqID=nseqs;
                sdata[nseqs].sequence = malloc(sqinfo.len*sizeof(char));
                sdata[nseqs].orisequence = malloc(sqinfo.len*sizeof(char));
                seq = remgaps(seq,sqinfo.len);
                sdata[nseqs].orisequence=strdup(seq);
                ToIUPAC(seq,TRUE);
                ToRNA(seq);
                s2upper(seq);
                ToACGUN(seq);
                sdata[nseqs].sequence=strdup(seq); 
                sdata[nseqs].seqlength=strlen(sdata[nseqs].sequence);
                sdata[nseqs].weight = 1;
                
        

        td = mhash_init(MHASH_MD5);
        if (td == MHASH_FAILED) exit(1);
        i = 0;        
                while (i < strlen(sdata[nseqs].sequence)) {
               mhash(td, &sdata[nseqs].sequence[i], 1);
                        ++i;
        }

             hash = mhash_end(td);

		/*hashfilename = (char *)calloc(mhash_get_block_size(MHASH_MD5)*2,sizeof(char));*/
		for (i = 0; i < mhash_get_block_size(MHASH_MD5); i++) {
				sprintf(&hashfilename[i*2],"%.2x",hash[i]);
		}
                

                if (readfromprobfile(hashfilename, nseqs)){        
                        print_data_content2(last_ptr[nseqs]);
#ifdef DEBUG                    
                        if (!verbose){
                                printf("\n");
                                print_data_content1(last_ptr[nseqs]);
                                printf("\n");
                        }
#endif
                } else {
                        pf_fold(sdata[nseqs].sequence, structure,sdata[nseqs].seqID, nsequences);
                        print2probfile(last_ptr[nseqs]);
		}
		nseqs++;
		
	}   
	

	   
/*        fprintf(stdout,"%s\n",matrixname);*/
	if (submat){
		readsubmatrix(matrixname);
	}
	create_sub_matrix();
	
        if (nseqs!=2){
                SeqfileClose(dbfp);
                align(nseqs,seqfile);
                if (!nofile){
                        chdir(workingdir);
                        if (opendir("./resultDIR")==NULL){
                                mkdir("./resultDIR",mode| S_IRGRP | S_IXGRP| S_IXOTH | S_IROTH);
                        }
                        chdir("./resultDIR");
                }
                createDistanceMatrix(array, nseqs);

                
                if (bion){
                        printdistmatrix(nseqs);
                        bionj("dist.phy");
                } else if (weigh){
                        printdistmatrix(nseqs);
                        weighbor(fopen("dist.phy", "r"),fopen("bionjout","w"));
                }

                createTree(nseqs);

                
                
                if (upgma){
                        njoin=0;
                        weighted=0;
                        for (i=0; i < nseqs-1; ++i){
                                free(prof[i].seqprofA);
                                free(prof[i].seqprofB);
                        }
                        free(prof);
                        createTree(nseqs);
                        
                }
                
                
                for (i=0; i < nseqs; ++i){
                        free(distmatrix[i]);
                        free(array[i]);
                }
                free(distmatrix);
                free(array);
#ifdef DEBUG            
                if (!verbose){
                        for (i=nseqs-2; i >=0;--i){
                                printf("\nCycle %d: %s %d(%d) joins %s %d(%d)\n",i,(prof[i].typeA==0 ? "Sequence" : "Cluster" ), prof[i].seqprofA[0],prof[i].pnoA, (prof[i].typeB==0 ? "Sequence" : "Cluster" ),
                                prof[i].seqprofB[0], prof[i].pnoB);
                        }
                }
#endif          
        } else {

                align(nseqs,seqfile);
                prof = (struct profile *)calloc(nseqs-1,sizeof(struct profile));
                prof[0].seqprofA = (int *)calloc(1,sizeof(int));
                prof[0].seqprofA[0] = 1;
                prof[0].typeA = 0;
                prof[0].nA=1;
                prof[0].seqprofB = (int *)calloc(1,sizeof(int));
                prof[0].seqprofB[0] = 2;
                prof[0].typeB = 0;
                prof[0].nB=1;
                if(!nofile){
                        chdir(workingdir);
                        if (opendir("./resultDIR")==NULL){
                                mkdir("./resultDIR",mode| S_IRGRP | S_IXGRP| S_IXOTH | S_IROTH);
                        }
                        chdir("./resultDIR");
                }
        
        }
        
	for (x = 0; x < nsample; x++) {
		int j;
		i1 = CHOOSE(nseqs);
		do {
			j = CHOOSE(nseqs); 
		} while (j == i1); /* make sure j != i1 */
		sum += PairwiseIdentity(sdata[i1].sequence,sdata[j].sequence);
	}

#ifdef DEBUG    
        if (!verbose){
                printf("average pairwise identity=%f\n",sum / (float) nsample);
        }

#endif

	if ((bion || weigh) && nseqs >3 ) {
		parsenewicktree(nseqs);
	}

	if (onlytree){
		int j,k,sizeA,sizeB;
		for (i=nseqs-2; i>=0;--i){
			for(j=i+1; j <= nseqs-2 ;++j){
				if (prof[i].typeA==0 && prof[j].seqprofA[0] == prof[i].seqprofA[0]){
					prof[i].typeA = 1;
				}
			}
		
			if (prof[i].typeA==1){
				for (j=i+1; j <= nseqs-2 ;++j){
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
		
			for(j=i+1; j <= nseqs-2 ;++j){
				if (prof[i].typeB==0 && prof[j].seqprofA[0] == prof[i].seqprofB[0]){
					prof[i].typeB = 1;
				}
			}
			
			if (prof[i].typeB==1){
				for (j=i+1; j <= nseqs-2 ;++j){
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
	
		}
		printtreefile(prof,nseqs);
		exit(0);
	}

	   
	malign(nseqs);

	print2file(nseqs, seqfile);
                
#ifdef DEBUG    
        if(!verbose){
                afp = MSAFileOpen(afile, fmt, NULL);
                msa = MSAFileRead(afp);
                printf("Alignment PWI: %f\n",MsaPwIdent(msa));
                MSAFileClose(afp);
        }
#endif  

        if (adjust){
                afp = MSAFileOpen(afile, fmt, NULL);
                msa = MSAFileRead(afp);
                adjustalpha(MsaPwIdent(msa),nseqs);        
#ifdef DEBUG    
                if(!verbose){
                        afp = MSAFileOpen(afile, fmt, NULL);
                        msa = MSAFileRead(afp);
                        printf("Alignment PWI: %f\n",MsaPwIdent(msa));
                }
                MSAFileClose(afp);
#endif                  
        }
        
        

        
        for (i=0; i < nseqs-1; ++i){
                free(prof[i].seqprofA);
                free(prof[i].seqprofB);
        }
        
        free(prof);


        if (nseqs!=2){
                while(bfiterate > 0){
                        bestfirst(nseqs,1);
                        --bfiterate;
                }
        }
        oldmsascore = msascore;
        srand(time(NULL));
        if (nseqs!=2){
                while (iterate>0){
                        msascore = doiteration(nseqs,5,msascore); 
                        print2file(nseqs, seqfile);                        
                        --iterate;
#ifdef DEBUG    
                        if(!verbose){
                                afp = MSAFileOpen(afile, fmt, NULL);
                                msa = MSAFileRead(afp);
                                printf("Alignment PWI: %f\n",MsaPwIdent(msa));
                                MSAFileClose(afp);
                                printf("\niterations left: %d\n", iterate);
                        }
#endif
                }
        }

        
        
        printf("\nalpha %f\n",alpha);
        if(!nofile){        
                printf("\nprogram exited normally\n");
        } else {
                print2stdout(nseqs);
                remove("./probDIR"); 
                remove(wname);
                remove("outfile");        
                remove("outtree");        
        }
        free(wname);
        for (i = 0; i < nseqs; ++i){
                curr_ptr = last_ptr[i];
                while(curr_ptr!=NULL){
                        old_ptr=curr_ptr->previous_ptr;
                        free(curr_ptr);
                        curr_ptr = old_ptr;
                }
        }
        free(afile);
        free(sdata);
        free(first_ptr);
        free(last_ptr);
        
        return(0);
}



/* Function:            MsaPwIdent
 * 
 * Purpose:             
 *           
 * Args:                        *msa                                        
 *
 * Returns:             pairwise sequence identity of msa file
 */
float MsaPwIdent(MSA *msa)
{
    float  **idmx; /* identity matrix */
    float sum, avgid;
    int i, j;

    MakeIdentityMx(msa->aseq, msa->nseq, &idmx);
    sum=0.0;
    for (i = 0; i < msa->nseq; i++)
        for (j = 0; j < i; j++)
            sum += idmx[i][j];
    avgid = sum / (float) (msa->nseq * (msa->nseq-1)/2.0);
    FMX2Free(idmx);

    return avgid;
}


/* Function:			readFile
 *
 * Purpose:
 *
 * Args:				
 *
 * Returns:			pairwise sequence identity of msa file
 */
void readFile(SQFILE *dbfp, char *seq, SQINFO sqinfo,int type, int nsequences){

}

/* Function:		
 * 
 * Purpose:		
 *           
 * Args:			
 *
 * Returns:		void
 */
void usagestral(void){
    printf("%s version %s\n", PACKAGE, VERSION);
 	 printf("usage: stral [options] -f <sequence-file.vie>\n"
	  		"StrAl"
			"\t[-A <alpha>] structural vs sequence weight <%5.2f>\n"
			"\t[-G <gapOpen>] gap open penalty <%5.2f>\n"
			"\t[-E <gapExt>] gap extension penalty <%5.2f>\n"
			"\t[-M <gapExtMSA>] gap extension penalty in MSA <%5.2f>\n"
			"\t[-O <gapOpenMSA>] gap open penalty in MSA <%5.2f>\n"
			"\t[-N] use neighbor joining <FALSE> upgma default\n"
			"\t[-W] weight sequences(midpoint method) <FALSE> \n"
			"\t[-r] root tree without assigning sequence weights <FALSE>\n"
			"\t[-b] use BIONJ's guide tree <FALSE>\n"
			"\t[-w] use weighbor's guide tree <FALSE>\n"
			"\t[-R] substitution matrix ribosum <TRUE>\n"
			"\t[-P] positive matrix <TRUE>\n"
			"\t[-B <#>] iteration: method best-first\n"
			"\t[-I <#>] iteration: method random split\n"
			"\t[--no-freeendgaps] toggle between full global alignment and alignment with free end gaps\n"
			"\t[--palign] keep consistent parts fixed in iterative step\n"
			"\t[--fft] use fast-fourier transformation\n"
			"\t[-T <temperature>] RNAfold working temperature <37>\n"
			"\t[-p] printing pairwise alignments to file <FALSE>\n"
			"\t[-n] don't print program output to files <FALSE>\n"
			"\t[-V] be verbose <TRUE>\n", alpha,gapOpen,gapExt,gapExtM,gapOpenM);
	exit(1);
}
