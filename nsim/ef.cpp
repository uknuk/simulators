#include "ef.h"
#include <memory.h>
#include <float.h>
#include <math.h>



double	EF::alpha = AlphaEF, EF::wO = 3, EF::wDr = 10;
int		EF::nA = 0, EF::max = 0;


EF::EF()
{
	memset(this, 0, sizeof(EF));
	for (int k = 0; k < 2; k++) { 
		list[k] = new List;
		q[k] = dalloc(EF::nA);
	}
}


EF::~EF()
{
	for (int k = 0; k < 2; k++) {
		delete q[k];
		delete list[k];
	}
}

void EF::update(int ok, int k, int unit)
{
	nC[k] += unit;
	if (!ok)
		nBl[k] += unit;

	if (nC[k] >= max) {
		double rew = nBl[k]/max + pow(wDr,wO*r[k]);
		q[k][a[k]] += alpha*(rew - q[k][a[k]]);
		r[k] = t[k] = 0;
		
		double qMin = DBL_MAX;
		int aOpt;
		for (int n = 0; n < nA; n++) 
			if (q[k][n] < qMin) {
				aOpt = n;
				qMin = q[k][n];
			}
		list[k]->put(new Node(a + k));
		a[k] = aOpt;
		nBl[k] = 0;
		nC[k] = 0;
	}
}

/* void EF::state(mxArray* qA, mxArray* acS, int c)
{
	double* qP = mxGetPr(qA);
	int sub[] = {c,0,0};
	mxArray* acA[2];
	for (int k = 0; k < 2; k++) {
		sub[1] = k;
		for (int a = 0; a < nA; a++) {
			sub[2] = a;
			qP[mxCalcSingleSubscript(qA,3,sub)] = q[k][a];
		}
		acA[k] = list[k]->write();
		mxSetCell(acS,c + 7*k,acA[k]);
	}
} */
