#ifndef MOB_H
#define MOB_H

#include "global.h"
#include "net.h"	  
// for SUN otherwise complains Net type incomplete


struct Mob {
	static	struct	Sim*	sim;
 	static	double	spf[3]; // fraction of fixed and high-speed users
	// for higher data rates (video & UDS) only fixed and cars 
	static   int		va;  // voice activity
	static   double   nonU; // non-uniformity, higher dencity in central area
	static   int cN;   // call number
	int      id;      // number of channel 
	int		c;      /* controlling cell */
	double	p[2];  /* power for both links */
	double	cst;  /* call start time */
	int		co;  /* call outage counter */
	int		io[2]; /* instant outage counter */
	double	sp;     /* speed */
	int		loc[2][2];	 /* coordinates for attenuation */
	double	g[2];  /* received sir counter for pc*/
	double   sirT[2]; // sir target
	double	sirThr[2]; // sir threshold for pc
	double	sl[7]; /* shadowing loss */
	double   al[7]; /* average loss */
	double	l0; // own attenuation from the previous iteration
	double	hm;  /* handoff margin */
	int		st; /* state */
	int		r[2];  /* maximum and minimum rate */ 
	int		cr[2];   /* current rate */
	double	th[2];  /* call throughput */
	int		b;		// rate in voice units
	int		sc;   // service class
	double	(*fi)[No]; // fading generator
	double   wm;  	// fading sinusoid
	double   t0;
	int		ro[2];  // running outage
	int		lro[2]; // last running outage;
	void setup(int sc); 
	void proceed(double t, double crt);
	void state(double t, struct Event* ev = 0);
	void shadow(double t = 0);	 // called from Mob without t
	void step(int, double pr);
	void genLoc();
	void receive(struct Base*, double t);
	void decode(double t, int ch);
	void transmit(int, struct Base*);
	void rateUp(int link);
	void fade(double* att, double t);
	void move(double t = 0);  // called from shadow without t
	void setSIR();
};
#endif

	
