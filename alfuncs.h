struct blocks {
		int *L;
		int *R;
		int *M;
		int count;
};



float recalculatescorePP(struct vector **ptrarray, int a);
char *symbols(struct vector *i_ptr,struct vector *j_ptr);
char *symbolsPP(struct vector **ptrarrayA,struct vector **ptrarrayB,int a, int b);
char *symbolsPPcpy(struct vector **ptrarrayA,struct vector **ptrarrayB,int a, int b);
float score(struct vector *s1_ptr, struct vector *s2_ptr);
void scorePP(struct vector **ptrarrayA,struct vector **ptrarrayB,struct profile p, float **scorematrix, struct lscore *spA, struct lscore *spB);
float scorePPfft(struct vector **ptrarrayA,struct vector **ptrarrayB,int a, int b);
void scorePPcpy(struct vector **ptrarrayA,struct vector **ptrarrayB,struct profile p, float **scorematrix, struct lscore *spA, struct lscore *spB,int i, int j);
float scoreSP(struct vector **ptrarrayA,int spp);

struct blocks phelix(int spp);
struct blocks findsegment(int spp);
