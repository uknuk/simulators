#ifndef EF_H
#define EF_H

#include "list.h"

extern double* dalloc(int);
extern int* ialloc(int);


const double AlphaEF = 0.1;


struct EF {
	static double	alpha;
	static int		nA;  // number of actions
	static int		max; // max nummber of nBl for update
	static double	wO;   // outage weight
	static double	wDr;  // drop weight
	int		nBl[2];
	int		nC[2];
	double	r[2];
	double   t[2];
	int		a[2];
	double*	q[2];
	List*		list[2];
	EF();
	~EF();
	void update(int ok, int k, int unit);
	void state(double* );
};

#endif


