struct segment {
	int beginA;
	int endA;
	int beginB;
	int endB;
	double score;
	double linkScore;
	int elI;
	int elJ;
	int blocks;
	int boolean;
	struct segment *father;	
	struct segment *child;
};

struct listPath {
	float score;
	int startI;
	int startJ;	
	struct segment *startpoint;
	struct listPath  *child;
	struct listPath *father;
}; 

void fftalign(struct vector *last1_ptr, struct vector *last2_ptr);
void fftmultalign(struct profile p, int lengthSEQ1, int lengthSEQ2, int nseqs);
void realft(float data[], int n, int isign);
void twofft(float data1[], float data2[],float fft1[], float fft2[],int n);
void correl(float data1[],float data2[], int n, float ans[]);
void segmentscorepair(int *,struct vector *,struct vector *);
void segmentscoremult(int *,struct vector **,struct vector **,int , int ,struct profile);
int findmaximizingpaths(int lenA, int lenB);
void lookForChildren(int i, int j,int lenA, int lenB);
int checkIsValidChildren(int i, int j, int k, int l);
int getValidHasChildren( int k, int l);
void retrievePath(int i, int j);
struct segment *getNext(struct segment *treeItem , int count);
int getAllPaths(int counts,int lenA, int lenB);
struct listPath *moveToStart(struct listPath *newItem);
