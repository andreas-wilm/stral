/* sre_random.h
 * Header file for sre_random.c
 *
 * SRE, Tue Oct  1 15:24:29 2002
 * CVS $Id: sre_random.h,v 1.1 2006-03-08 22:45:00 wilm Exp $
 */

extern double sre_random(void);
extern void   sre_srandom(int seed);
extern double sre_random_positive(void);
extern double ExponentialRandom(void);
extern double Gaussrandom(double mean, double stddev);
extern int    DChoose(double *p, int N);
extern int    FChoose(float *p, int N);

#define CHOOSE(a)   ((int) (sre_random() * (a)))

  
