#include <stdio.h>

float doiteration(int n, int rems, float score);
void adjustalpha(float pwident, int n);
void bestfirst(int n, int rems);
void removegaps(int rn, int rems);
void remgapcols(int *, int );
void removegapsfromcpy(int *removes, int rems);
void remgapcolsfromcpy(int *removes, int n);
void copyvectors(int n);
struct blocks mergeblocks(struct blocks helix, struct blocks segment);
