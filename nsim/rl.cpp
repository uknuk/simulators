#include "rl.h"
#include "net.h"
// #include "ef.h"
#include "mob.h"
#include "sim.h"
#include <stdlib.h>
#include "global.h"


double RL::alpha = 0.1, RL::gamma = 0.1; 
double RL::rT = 1e-3, RL::qT = 1e-1, RL::pMax[2] = {10,50}; 
double RL::tau = 0, RL::rTau= 0, RL::lp[2] = {0,0};
int RL::nL = 0, RL::nP = 0;
Sim* RL::sim = 0;
 
																						                
void RL::init(int _c, int sg)
{
	
	c = _c;
	int size[] = {nP, nL};
	int S = nP*nL;
	for (int k = 0; k < 2; k++) {
		create(&av[k], S);
		if (sim->rc) 
			create(&rv[k], nP);
		if (nL > 1 && nP > 1 && sg ) {
			create(&sv[k], S);
			for (int i = 0; i < 2; i++) {
				create(&tet[k][i], size[i]);
			}
		}
	}
}

void RL::reset()
{
	dec = ri = *aps = aps[1] = *rps = rps[1] = ls = 0;
}

RL::~RL(void)
{
	for (int k = 0; k < 2; k++) {
		if (av[k])
			delete av[k]; 
		if (rv[k]) 
			delete rv[k]; 
		if (tet)	{
			delete sv[k];
			for (int i = 0; i < 2; i++)
				delete tet[k][i];
		}
	}
	memset(this, 0, sizeof(RL));
}
			

int RL::admit(struct Mob* m) 
{
	int k, ok = 1;
	int l = 0; 	// load state
	int p[] = {0,0};	// power states
	double g[] = {0,0}; // V(st+1)
	Net* net = sim->net;

	if (nL > 1) {
		if (sim->type > 2) {
			if(nL < 10)
				l = m->sc;
			else // if (!net->rc)
				l = net->bs[c].l - Nmin;
			// else
			//	l = net->bs[c].l + m->b;
		}
		else
			l = net->bs[c].n - Nmin;
	
		if (l >= nL)
			l = nL;	 // if l exceeds nL
	}	
	
	for (k = sim->L[0]; k < sim->L[1]; k++)
		if (ok && (sim->type == FCSIM || 
			(sim->type != FCSIM && net->bs[c].ao[k] < qT))) {
				int& _p = p[k];
				ok = decide(av, k, _p, g + k, l);
		}
		else {
			ok = 0;
	//		net->bl++;
		}
	
	/* if (net->ef)  // admitted if k = 2, otherwise blocked
			net->ef[c].update(k,ok,m->b);	 */
		
	if (dec)	
		for (k = sim->L[0]; k < sim->L[1]; k++)  
			update(av, k, g[k]);
	
	if (ok) {
		dec = 1;
		ls = l;
		for (k = 0; k < 2; k ++) 	 
			aps[k] = p[k];
	}
	
	return ok;
}
		

int RL::rateUp(int k)
{
	int p = 0;
	double g = 0;
	int& _p = p;
	int ok =  decide(rv, k, _p, &g); // l = rc default
	if (ri)	  // rate increase
		update(rv, k, g, 1);
	if (ri = ok)
		rps[k] = p;
	return ok;
}
	
	
int RL::decide(double* v[], int k, int& p,  double* g, int l)
{
	int ok = 1;
	int rc = 0;

	if (l == RC)  {
		rc = 1;
		l	= 0;
	}

	/*if (nP > 1) {
		// if (net->bs[c].p[k] >= net->adp[(k + 1)*nP - 1])	
		//			ok = 0;  // too high power 
		// else // get power state 
			for (p = 0; p < nP; p++) 
				if (sim->net->bs[c].p[k] < psp[k*nP + p]) 
					break;
			if (p == nP)
				p = nP - 1;
	}
	else
		p = 0;  */

	double pow = sim->net->bs[c].p[k];
  	if (!k)
		pow = log10(pow/sim->nu);
  	p = (int) floor(lp[k]*log10(pow));
	if (p > nP -1)
		p = nP - 1;
		
	
	double q = 0;

	if (ok) {
		if (tet[0][0] && !rc && !sv[k][p*nL + l]) {
			int sub[] = {p,l};
			for (int i = 0; i < 2; i++)
				q += tet[k][i][sub[i]];
		}
		else
			q = v[k][p*nL + l];

		/* if (!q && l != RC && net->bs[c].ao[k] > qT)
			ok = 0; */
	
		if (ok && q < 0) {
			double tauV;
			if (!rc) {
				// if (!net->ef)	 /* constant tau */
					tauV = tau; 
				/* else
					tauV = tau[net->ef[c].a[k]];	*/
			}	
			else
				tauV = rTau; /* rate control */
			ok = ran() > 1 - exp(q/tauV);
		}
	}
	*g = ok ? q : 0;
	return ok;
}

void RL::update(double* v[], int k, double g, int rc)
{
	int p, l;
	double rew;
	Net* net = sim->net;
	if (rc) {
		p = rps[k];
		l = 0;
		rew = -net->bs[c].o[k];
	}
	else {
		p = aps[k];
		l = ls;
		rew = rT - net->bs[c].ao[k]; // net->bs[c].out[k];
		/* if (net->ef)
			net->ef[c].r[k] += (net->bs[c].ao[k] - net->ef[c].r[k])/++net->ef[c].t[k];	*/
		net->bs[c].ao[k] = 0;
		if (k)
			net->bs[c].t = 0;
	}
		
	int ind = p*nL + l;

	double delta = alpha*(rew + gamma*g - v[k][ind]);
	v[k][ind] += delta;
	if (tet[0][0] && !rc) {
		int sub[] = {p,l};
		if (!sv[k][p*nL + l])
			sv[k][p*nL + l] = 1;
	
		for (int i = 0; i < 2; i++)
			tet[k][i][sub[i]] += delta;
	}
}

void RL::state(double* v)
{
	/* int subA[4] = {c,0,0,0};
	int subR[3] = {c,0,0};
	int subT[2][3] = {{c,0,0},{c,0,0}};
	mxArray* tetA[2];
	double *tetP[2];
	double* avP = mxGetPr(ar[VA]);
	double *rvP;
	if (ar[RA])
		rvP = mxGetPr(ar[RA]);
	
	if (tet[0][0])
		
		for (int i = 0; i < 2; i++) {
			tetA[i] =  mxGetCell(ar[TET],i);
			tetP[i] = mxGetPr(tetA[i]);
		}
	
	for (int k = 0; k < 2; k++) {
		subA[1] = subR[1] = subT[0][1] = subT[1][1] = k;
		
		for (int p = 0; p < nP; p++) {
				subA[2] = subR[2] = p;
				if (tet[0][0]) {
					subT[0][2] = p;
					tetP[0][mxCalcSingleSubscript(tetA[0],3,subT[0])] = tet[k][0][p];
				}
				if (rv[0])	
					rvP[mxCalcSingleSubscript(ar[RA],3,subR)] = rv[k][p];
				for (int l = 0; l < nL; l++) {
					subA[3] = l;
					avP[c + 7*(k + 2*(p + l*nP))] = av[k][p*nL + l];
					if (tet[0][0]) {
						subT[1][2] = l;
						tetP[1][mxCalcSingleSubscript(tetA[1],3,subT[1])] = 
							tet[k][1][l];
					}
				}
		}
	}
}	 */

	for (int i = 0; i < nP; i++)
		for (int k = 0; k < nL; k++)
			v[i*nL + k] = av[1][i*nL + k];
}




