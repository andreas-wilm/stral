float methodmultcpy(struct profile p, int lenA, int lenB, int nseqs);
float alignpieces(struct profile p, int part, int splits,struct vector ***segptr);
float **retree(int n);
float **data_content_cpy(struct vector *current_ptr, int len, float **arrayX);
