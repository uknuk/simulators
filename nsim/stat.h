#ifndef STAT_H
#define STAT_H

#include "list.h"


const int lS = 10;  // load statistics size

struct Stat {
	double	(*l)[lS], *r, *t;	// load, radio, throughput
	List<double>*	  list[2];
	int	  tl;
	int     nR;   // number of radio parameters
	int	  cN;	  // number of load entries
	int	  nFr;  // number of frames
	int	  nL;	 // number of links
	int     lN;   // load stat entry number, 
	// needed cause for sType == 4 load stat entry is different from arrivals (only voice are considered)
	Stat(int nC, int nFr, int sType);
	void load(struct Base* base, struct Mob* m);
	void radio(int ptr, struct Base* base, int link);
	void failure(double t, int type);
	void write(struct Sim* sim);

};

#endif


