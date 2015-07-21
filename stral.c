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




float alpha= 8;
float gapOpenM= 8;
float gapExtM = 0.5;
float gapOpen = 8;
float gapExt  = 0.5;
int   njoin   = 0;
int   positive= 1;
int   ribosum = 1;
int   weighted= 0;
int   rooted  = 0;
int   verbose = 0;
int	ppair   = 0;
int	bion    = 0;
int	weigh   = 0;
int	upgma   = 0;
int 	nofile  = 0;
int palign   = 0;
int noterm   = 0;






void usagestral(void);
void weighbor(FILE *,FILE*);
float MsaPwIdent(MSA *msa);

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

struct opt_s OPTIONS[] = {
	{ "-a", TRUE, sqdARG_NONE },    
	{ "-h", TRUE, sqdARG_NONE },    
	{ "--gccomp",   FALSE, sqdARG_NONE },
	{ "--informat", FALSE, sqdARG_STRING },
	{ "--quiet",    FALSE, sqdARG_NONE   },
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))




int main(int argc, char **argv){
  	char     *seqfile;    	/* name of sequence file     */
  	SQFILE   *dbfp;			/* open sequence file        */
  	int       fmt;				/* format of seqfile         */
  	char     *seq;				/* sequence                  */
  	SQINFO    sqinfo;       /* extra info about sequence */
	char     *afile;              /* name of aligned sequence file */
	MSAFILE  *afp;		/* pointer to open alignment file*/
	MSA      *msa;                /* multiple sequence alignment   */  
  	int       type,i;			/* kAmino, kDNA, kRNA, or kOtherSeq */
  	
	int    allreport;			/* TRUE to do a short table for each sequence */
  	int    be_quiet;			/* TRUE to suppress header */
 	int    do_gccomp;			/* TRUE to include GC composition in per-seq report */
  
  	char  *optname;
  	char  *optarg;
  	int    optind;
	
	int    nseqs=0;
  	int	 nsequences=0;		/*cmp  nseqs*/
  	char  *structure=NULL;	/*parameter to pf_fold*/
	mode_t mode = S_IRWXU;
	struct vector *curr_ptr, *old_ptr;
	int i1,x,nsample=1000;
	float sum =0.0;
	int len;
	
	/***********************************************
	* Parse command line
	***********************************************/
	double temperature = 37.;
	unsigned char *hash;
		
	fmt       = SQFILE_UNKNOWN;	/* default: autodetect format  */
	allreport = FALSE;				/* default: file summary only  */
	be_quiet  = FALSE;				/* show header info by default */
	type      = kOtherSeq;		/* just to silence gcc uninit warning */
	do_gccomp = FALSE;

	if (argc<2)
      usagestral();

	for (i=1; i<argc; i++) {
		if (argv[i][0]=='-') {
      		if  (!(strcmp(&argv[i][1],"-noterminal"))){
				noterm=0;
			}else {
				switch ( argv[i][1] ){
					case 'N':
						njoin=1;
						break;	
					case 'W':
						weighted=1;
						njoin=1;
					break;
					case 'A': 
						if (argv[i][2]!='\0') 
							usagestral();
							if(i==argc-1)
							usagestral();
							if (sscanf(argv[++i], "%f", &alpha))
							;		
						else 
							usagestral();
						break;
					case 'G': 
						if (argv[i][2]!='\0') 
							usagestral();
							if(i==argc-1) 
							usagestral();
							if (sscanf(argv[++i], "%f", &gapOpen))
							;
						else 
							usagestral();
						break;
					case 'E': 
						if (argv[i][2]!='\0') 
							usagestral();
						if(i==argc-1) 
								usagestral();
						if (sscanf(argv[++i], "%f", &gapExt))
							;
							else 
							usagestral();
						break;
					case 'M': 
						if (argv[i][2]!='\0')
							usagestral();
						if(i==argc-1) 
							usagestral();
						if (sscanf(argv[++i], "%f", &gapExtM))
							;
						else 
							usagestral();
						break;
					case 'O': 
						if (argv[i][2]!='\0') 
							usagestral();
						if(i==argc-1) 
							usagestral();
						if (sscanf(argv[++i], "%f", &gapOpenM))
							;
						else 
							usagestral();
						break;
					case 'R':
						ribosum=!ribosum;
						break;	
					case 'P':
						positive=!positive;
						break;
					case 'V':
						verbose=1;
						break;
					case 'r':
						rooted=1;
						njoin=1;
						weighted=1;
						break;
					case 'T': 
						if (argv[i][2]!='\0') 
							usagestral();
						if(i==argc-1) 
							usagestral();
						if (sscanf(argv[++i], "%lf", &temperature))
							;
						else 
							usagestral();
						break;
					case 'p':
						ppair=!ppair;
						break;
					case 'u':
						njoin=1;
						weighted=1;
						upgma=1;
						break;
					case 'b':
						bion=!bion;
						break;
					case 'w':
						weigh=!weigh;
						break;
					case 'n':
						nofile=!nofile;
					break;
					default: 
						usagestral();
				}
			}
  		}
	}

 


	if (!njoin && weighted){
		printf("-N/-W resp. -N/-r combination impossible option\n");
		usagestral();
	}
	if (nofile && ppair){
		printf("-p/-n  combination impossible option\n");
		usagestral();
	}	
	
	argv[1]=argv[argc-1];
	argc = 2;

  
	while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage, &optind, &optname, &optarg)){
		if      (strcmp(optname, "-a")       == 0)  
			allreport = TRUE; 
		else if (strcmp(optname, "--quiet")  == 0)  
			be_quiet  = TRUE; 
		else if (strcmp(optname, "--gccomp") == 0)  
			do_gccomp = TRUE; 
		else if (strcmp(optname, "--informat") == 0) {
			fmt = String2SeqfileFormat(optarg);
			if (fmt == SQFILE_UNKNOWN) 
				Die("unrecognized sequence file format \"%s\"", optarg);
      	} else if (strcmp(optname, "-h") == 0) {
			puts(usage);
			puts(experts);
			exit(EXIT_SUCCESS);
		}
	}
	
	if (argc - optind != 1) Die("%s\n", usage);
		seqfile = argv[argc-1];



	/* Try to work around inability to autodetect from a pipe or .gz:
	* assume FASTA format
	*/
	if (fmt == SQFILE_UNKNOWN && (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
		fmt = SQFILE_FASTA;
	
	/***********************************************
	* Read the file.
	***********************************************/

  
	
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
		char *hashfilename;
		
		
		
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
/*		ToIUPAC(seq);*/
		ToRNA(seq);
		s2upper(seq);
		ToACGUN(seq);
		sdata[nseqs].sequence=strdup(seq); 
		sdata[nseqs].seqlength=strlen(sdata[nseqs].sequence);
		sdata[nseqs].weight = 1;
		
      	free(seq);

       	td = mhash_init(MHASH_MD5);
       	if (td == MHASH_FAILED) exit(1);
       	i = 0;	
		while (i < strlen(sdata[nseqs].sequence)) {
               mhash(td, &sdata[nseqs].sequence[i], 1);
			++i;
       	}

	     hash = mhash_end(td);
		
		/* valgrind ERROR */
       	/*printf("mhash_get_block_size %d\n",mhash_get_block_size(MHASH_MD5));*/
		hashfilename = (char *)malloc(mhash_get_block_size(MHASH_MD5)*2*sizeof(char));
		for (i = 0; i < mhash_get_block_size(MHASH_MD5); i++) {
			 /* valgrind ERROR */
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
		free(hash);
		free(hashfilename);
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


	
	

	
	for (i=0; i < nseqs-1; ++i){
		free(prof[i].seqprofA);
		free(prof[i].seqprofB);
	}
	
	free(prof);



		
	
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
		free(sdata[i].ID);
		free(sdata[i].sequence);
		free(sdata[i].orisequence);
		
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



/* Function:		MsaPwIdent
 * 
 * Purpose:		
 *           
 * Args:			*msa					
 *
 * Returns:		pairwise sequence identity of msa file
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
 	 printf("usage: stral [options] <sequence-file.vie>\n"
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
			"\t[-T <temperature>] RNAfold working temperature <37>\n"
			"\t[-p] printing pairwise alignments to file <FALSE>\n"
			"\t[-n] don't print program output to files <FALSE>\n"
			"\t[-V] be verbose <TRUE>\n", alpha,gapOpen,gapExt,gapExtM,gapOpenM);
	exit(1);
}


