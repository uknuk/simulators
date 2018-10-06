#ifndef RL_H
#define RL_H

extern double* dalloc(int);
extern int* ialloc(int);

enum RLtype {RC = -1,AC} ; // rate control flag

struct RL {
	static int	nP; /* power states  first dimension*/
	static int	nL; /* second dimension */
	static double  tau; // tau values
	static double  rTau; // rc tau
	static double	alpha; 
	static double	gamma;
	static double	rT;
	static double	qT;
	static struct  Sim* sim;
	static double  pMax[2];
	static double	lp[2];	// logarithmic power coefficient	

	int		c;
	int		dec; // RL admission decision
	double*	av[2]; // admission values
	int*		sv[2];  // state visited
	double*	rv[2]; // rate control values
	int		ri; /* rate increased */
	int		aps[2]; // admission power state
	int		rps[2]; // rc power state
	int		ls;	  // load state
	double*  tet[2][2];
	
	RL() {}; // for SUN
	~RL();
	void init(int c, int sg);
	int admit(struct Mob* m);
	int rateUp(int link);
	int decide(double* v[], int link, int& p,  double* g, int l = RC);
	void update(double* v[], int link, double g, int rc = 0);
	void state(double* );
	void reset();
};

#endif
