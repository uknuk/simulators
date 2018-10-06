#include "mob.h"
#include "global.h"
#include "sim.h"
#include "stat.h"
#include "rl.h"
#include "heap.h"
#include <float.h>
#include <stdlib.h>



 

Sim* Mob::sim = 0;
double Mob::spf[3] = {0.1,0.75,0.5};	
// default speed distribution: for voice and audio: 0, 1 & 10 m/s 
// for video and data: 0, 10 m/s 
int Mob::va = 0;
double Mob::nonU = 1.5;  // default non-uniformity
int Mob::cN = 1; 

void Mob::setup(int _sc)
{
	memset(this,0,sizeof(Mob));
	int k;
	sc = _sc;
	genLoc();
	if (!sc && va) {  // voice with voice activity
		st = (int) 0 + 5*(ran() < 0.5);
		for (k = 0; k < 2; k++)
			r[k] = sim->r[k];
	}
	else  {
		cr[DL] = cr[UL] = r[MIN] = sim->r[2*sc];
		if (sim->type == RCSIM  || sc == UDS) /* rate or flow control */
			r[MAX] = sim->r[2*sc + 1];
	}
	b = r[DATA]/sim->r[VOICE];
	
	double pr = rand();

	if (sc < 2)	{
		for (k = 0; k < 2; k++)	
			if (pr > spf[k])
				sp = speed[k];
	}
	else if (pr > spf[2])
		sp = speed[1];

	if (sp)
		t0 = 1/sp;
	shadow();  // to find out the strongest BS

	
}	

void Mob::proceed(double t, double crt)
{
	id = cN++;
	sim->net->heap->put(new Event(crt,REL,this,id));
  	if (va)
		state(t);
	cst = t;
	
	if (!sim->type && sp) {	// init fading generator
		wm = 2*Pi*sp/lambda;
		fi = new double[7][No];
		for (int i = 0; i < 7; i++)
			for (int n = 0; n < No; n++)
				fi[i][n] = 2*Pi*ran();
	}

 }


void Mob::state(double t, Event* ev)
{
	double u;
	if (ev) {
		u = ran();
		st = nTr[st][u > stPr[st]];
	}
	u = ran();
  	double tSt = t + trT[st] - log(u)/trR[st];
	if (!ev)
		ev = new Event(tSt,VATR,this,id);
	else
		ev->t = tSt;
	sim->net->heap->put(ev);
	for (int l = *sim->L; l < sim->L[1]; l++)
			cr[l] = r[bR[st][l]];
}



void Mob::move(double t)
{
	if (!t || (t && fmod(t,t0) < 1e-6)) {
		// fr = (int) t0;
		for (int k = 0; k < 2; k++)	 { // for x and y
			loc[OLD][k] = loc[NEW][k];
			loc[OLD][k] = loc[NEW][k];
			step(k,ran());
		}
	// fr--;
	}

	
		
	
}


void Mob::shadow(double t)	 // decompose ??
{
	int j, k, K; 
	double L[2] = {0,0}, ratio;
	double Lmin = DBL_MAX;

	if (t && sp) {
		ratio	= fmod(t,t0);
		if (ratio < 1e-6)
			move();
	}
							
   for (j = 0; j < 7; j++) { /* calculate loss for j-th bs */
		if (t && sp) // if in call && moving
			K = 2;
		else if (!t) // setup
			K = 1;
	   
		for (k = 0; k < K; k++) 
			L[k] = sim->att[loc[k][0]][loc[k][1]][j];
		
		sl[j] = L[0];
		
		if (t && sp) 
			sl[j] += (L[NEW] - L[OLD])*ratio;
		


		if (sim->type && ((t && sp) || !t)) {	 // no pc
			al[j] = sl[j];
			if (sl[j] < Lmin)	{
				Lmin = sl[j];
				c = j;
			}
		}
	}
	g[0] = g[1] = 0; /* reset */
}

void Mob::step(int k, double rx)  {
	int x, dx;
	
	dx =  - 1 + (int) floor(3*rx);
	x = loc[NEW][k] + dx;
	if (x >= attS[k] || x < 0)
		x = loc[NEW][k] - dx;
	loc[NEW][k] = x;
}

void Mob::genLoc()
{
	double rx;
	int x;
	for (int k = 0; k < 2; k++) { // for x and y
		 rx = ran();
		 if (nonU > 1 & rx < hsA*nonU)
			 x = C0[k] + (int) floor(C[k]*rx);
		 else
			 x = (int) floor(attS[k]*rx);
		if (x == attS[k])
			x = attS[k] - 1;
		loc[NEW][k] = loc[OLD][k] = x;
	}
}

void Mob::receive(Base* bsp, double t)
{
	
	int j, k;
	double I[2] = {0, sim->nu}, Pt;	// sir ratio
	double* loss = new double[7];
	
	Net* net = sim->net;
		
	if (!sim->type && sp) 
		fade(loss,t);
	else
		for (j = 0; j < 7; j++)
			loss[j] = sl[j];
	
	if (sc != 3)
		hm = al[c]*Delta;
	else
		hm = al[c];  // no soft handof for the shared channel

	for (j = 0; j < 7; j++) {
		Pt = 	al[j] <= hm ? of*(net->bs[j].p[DL] - p[DL]) : net->bs[j].p[DL];
		I[DL] += Pt/loss[j];
	}
	// ise of of for the pc simulation ?
	// I[DL] = (I[DL] + net->nu)*loss[c] - of*p[DL];

	if (*sim->L == UL) {	
		if (sim->type) { // no pc	and uplink
			if (!l0)
				l0 = loss[c];
			I[UL] = net->bs[c].p[UL] - p[UL]/l0;
			l0 = 	loss[c];
		}	
		else 
			I[UL] = net->bs[c].p[UL]; 
		// otherwise negative I can occur		
	}				   
	
	for (k = sim->L[0]; k < sim->L[1]; k++) {
		
		if (sim->type || !p[k]) /* power update or open loop */
			p[k] = sirT[k]*I[k]*loss[c];

		if (sim->type) {
			if (p[k] <= pMax && ran() > FER[sc > 0 && sc < 3])
				g[k] = 0;  // for not voice and not UDS FER[1]
			else
				g[k] = 1;
		}
			
		else  {   // power control
			// double sirThr = sirT[k]*sirMargin;  // sir threshold
			double recSir = (p[k]/loss[c])/I[k];  // received sir
			p[k] *= recSir <= sirThr[k] ? up : down;
			double ratio = recSir/sirT[k];
			g[k] += ratio;
		}
		
		if (p[k] > pMax)
			p[k] = pMax;
	}

	if (*sim->L == UL)
		for (j = 0; j < 7; j++)
			bsp[j].p[UL] += p[UL]/loss[j];
	delete loss;
}

void Mob::decode(double t, int ch)
{
	int k, out[2] = {0,0};
	int rel = 0;
	Net* net = sim->net;

	for (k = sim->L[0]; k < sim->L[1]; k++)  {
		if (!sim->type) 
			out[k] = g[k]/sim->nL <= 1;	// mean cir for the frame
		else 
			out[k] = (int) floor(g[k]);
		
		
		if (out[k])
			io[k]++;
		else
			io[k] = 0;

		if (sim->type == 1)
			net->thr[k] += cr[k]; // to check balance of voice activity

		if (sim->rc) {
			if (cr[k] > r[MIN])
					cr[k] -= sim->r[MIN];	/* decrease rate */
			else if (cr[k] < r[MAX]) {
				 if (sim->rc == DRC)
					 cr[k] += sim->r[MIN];
				else if (p[k] < net->bs[c].pMin[k]) {
					net->bs[c].pMin[k] = p[k];
					net->bs[c].chPmin[k] = ch;
				}
			}
				
			th[k] += (cr[k] - r[MIN])/sim->r[VOICE]; 
			// throughput exceeding rMin in rVoice
		}			
	}
	
	if (out[UL] || out[DL])
		co++;

}


void Mob::transmit(int n, Base* bs)
{
	int k;
	
	bs->p[DL] += p[DL];
		
	if (n == sim->nL - 1)
		bs->l += b;
		bs->n++;	 // current load is counted in the end of frame
	  	for (k = sim->L[0]; k < sim->L[1]; k++)  {
			if (io[k] && sim->net->stat)
				bs->fer[k]++;
			if (io[k] > 2) {
				bs->o[k] += b;
				ro[k]++;
			}
		}	
}


void Mob::fade(double* att, double  t)
{	 	
	double w, fd, fn, beta, Xi, Xq;
	double Lmin = DBL_MAX;
	
	for  (int j = 0; j < 7; j++) {
			Xi = Xq = 0;
			fd = cos(wm*t + fi[j][7]);
			for (int n = 0; n < No; n++) {
				w =  wm*cos(2*Pi*(n + 1)/N);
				beta = Pi*(n + 1)/No;
				fn = cos(w*t + fi[j][n]);
				Xi += cos(beta)*fn;
				Xq += sin(beta)*fn;
			}
			Xi = 2*Xi + fd;
			Xq = 2*Xq + fd;
			att[j] = sl[j]*(pow(Xi,2) + pow(Xq,2))/M;
			al[j] = (1 - ff)*al[j] + att[j]*ff;
			if (al[j] < Lmin)	{
				Lmin = al[j];
				c = j;
			}
  	}
}


void Mob::setSIR()
{
	int v;
	if (!sim->type) 
		v = 0;		 // the same sirT for all speeds
	else
		v = sp > 1 ? 1 : 0;
	
	for (int l = sim->L[0]; l < sim->L[1]; l++) 	{
		if (sim->type == RCSIM || (sim->type == FCSIM && sc))	{
			double lw;	
			// 0th polynomial coefficient, link weight
			// if VBR low-delay DL is 2 dB more, for UDS (s = 3) 1 dB 
			if (sc < 3)  // for UDS DL Eb/Io is lower 
				lw = 2;
			else
				lw = 1;
			double x = log(cr[l])/ln2;
			sirT[l] = *Poly + Poly[1]*x + Poly[2]*pow(x,2);
			sirT[l] += (l ? lw : 0) + v;
			sirT[l] = pow(10,sirT[l]/10)*cr[l]*1e3/wUTRA;
		}	
		else
			sirT[l] = sim->cir[4*sc + 2*v + l]*cr[l]*1e3;

		if (!sim->type) { // outer loop
			if (!sirThr[l]) // initialisation
				sirThr[l] = sirT[l]*sirUp;
			/* else
					sirThr[l] *= io[l] ? sirDown : sirUp; */
		}
	}
}
		



