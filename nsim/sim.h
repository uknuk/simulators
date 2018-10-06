// header file for class Sim

#ifndef SIM_H
#define SIM_H

#include <stdio.h>
#include <memory.h> 


struct Sim {
	struct Net*	net;
	FILE*			cdf;
	FILE*       crf;
	double**		adp;
	double**		rcp;
	int			adN;  // number of admission types
	int			nR;   // number of runs
	int			dls;  // size of data load
	FILE*       dlf;	// data load file
	int			t;		 // sim time
	int			type;	 // sim type
	double		nu;   /* noise */
	double*		cir;
	int			nCh;
	
	const int*	r;  /* array of rates */
 	double**		cd;
	int*			cds; // call data size
	double		**dl;  // data load
	double*     dlv; // data load value
	int         Ncdp;
	double      tFr;
	int         tDr;
	int			Ncrp;  // number of call results parameters
	int			nDu;  // number of data users
	int			trace;  // trace flag
	double		att[260][250][7];
	const int*	L; // number of links
	int			nL; // number of PC loops
	struct RL*	rl;
	struct SLA* sla;
	int         adm;	  // admission
	double		adt[2]; // admission thresholds
	int			rc;
	
	Sim(int argc, char* argv[]);
	void  setPar();
	void	mrun(int m);
	void	run(char*);
	void 	parseAdPar(char* adm);
	int	get(FILE*, double** data);
	double*	parse(char*, const char*);
	void setName(char* );
	void open(char* fname, char arg);
	void quit(char* file = 0);
	double* genLogSp(double* par, int N, int K); // log space
	void setAdm(int n);
	void getState();
	void write(char* fname, double* data, int size);
	void count();
	void genLoad(int load); // int* for 
	~Sim();

};



template <class Type>
	void create(Type** x, int size)
{
	*x = new Type[size];
	memset(*x, 0, size*sizeof(Type));
}

// own random generator to get identical results in PC and UNix
// from Numerical Reciepies in C p. 279
const int		IA = 18807;
const long		IM = 2147483647;
const double 	AM = 1.0/IM;
const long		IQ = 127773;
const int		IR = 2836;

double ran(long seed = 0);

#endif
// Don't go below this line to avoid SUN complains
