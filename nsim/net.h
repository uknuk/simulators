// header file for net class

#ifndef NET_H
#define NET_H


#include "stat.h"

struct Base {
	double	p[2]; /* power */
	int		l;		// load - sum of minimal rates
	int		n;    // number of users
	double	o[2];  /* cell outage (2 consecutive frame errors) */
	double   fer[2];  // cell frame error rate
	double	g[2];  /* cell cir */
	double	pMin[2]; // minimum DL power for rate control
	int		chPmin[2]; // channel of min power
	double   ao[2];  // average outage for QA
	int		t;	  // counter for average outage
	double   out[2];  // average outage between admissions
	double   lt;  // last admission time;
	void     init(struct Sim* sim); 
};

struct 	Net { 
	struct Sim*		sim;
	Base				bs[7];
	struct Mob*		m;
	struct DMob*   dm;
	struct Stat*	stat;
	double		   thr[2]; /* througput statistics */
	int				nM; // number of mobiles
	double*			cr;  /* call results */
	double 			pMin[2];
	struct PS*		ps;
	struct Heap*	heap;
	int				cN;
	int				fr;
	int            rcN;  // released call Number
	double         t;
	double*        cd;
	int				cds;
			
	Net(Sim* _sim, int n);
	~Net();
	int  run();
	void arrival(double* cd);
	void handle();
	int  admit(Mob* m);
	void radio();
	void transmit(int n, Base* bs);
	void rateUp(int c, int link);
	// void read(char* file);
	// void countOut(int c, double t);
 	void release(Mob*, int value = 0);
	void traceRate(double r) {stat->t[fr] = r; };
 };

#endif
