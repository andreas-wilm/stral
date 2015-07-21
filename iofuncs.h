char *wname;
void print2probfile(struct vector *ptr);
void ToACGUN(char *seq);
void print2file(int nseqs, char *filename);
void printdistmatrix(int n);
int readfromprobfile(char *filename, int n);
char *remgaps(char *seq, int len);
void remgapcols1(int n);
int recursion(int begin,char **nwtree, int count, int rn, int ***array, int p);
void parsenewicktree(int nseqs);
void printpairwise(int seq1,int seq2, char *seqA, char *seqB, int boolean, char *filename);
void print2stdout(int nseqs);
void printstart2end(int begin, int end, struct vector *start_ptr);