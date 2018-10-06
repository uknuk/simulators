#ifndef DATA_H
#define DATA_H

#include "list.h"
#include "mob.h" 

const int nRates = 32;
const int nPowStates = 10;


const double cSIR = 1.5e3*pow(10,0.6)/3.84e6; // SIR for control channel
const double pMin = 1e-3;  //power boundary for data control channel

enum State {OFF, ON};

struct PS;

struct DMob : public  Mob {
	double (*dl)[2]; // load data
	int nR;	 // 2-nd dimension of DL max number of requests
	double	data; // size of the current request
	double   dataRest;  // pendling data
	double   pkt;   // IP packet
	double   pktRest;
	double   pdu;
	int		req;   // number of processed requests
	double	tSt; // time in state
	double   cp;  // control power
	struct Event*   ev;
	int init(int id, double* data, Heap* heap);
	int decode(PS*, int c, Base* bs); // redefined from Mob
	void transmit(int rate); // redefined from Mob
  	int state(double t);
	void jump(double t);
	void control(Base* bs);
};


struct SLA {  //  RL rate control, originally Stochastic Learning Automata
	  static	double	a;  // a-parameter
	  static double   g;
	  static int		nP;
	  static int      nR;
	  static double	lp;
	  static double	voc; // voice outage cost
	  static double      rMin;
  	  double				(*q)[nRates]; 
	  int					p; // power state
	  int					r;  // last rate (action)
	  
	  SLA();
	  ~SLA() { delete[] q; };
	  int select(int rate, double pow);
	  void update(int rate, double pow, int ack, double out);
	  void getState(int rate, double pow);
	  double* storeState(int*);
	  void reset() {p = r = 0;};
};

struct PS {  // packet scheduler
	struct Sim* sim;
	List<DMob*>	q[7];
	Node<DMob*>*	sm[7];  // scheduled mobile
	Node<DMob*>*   pm[7];  // previous mobile for PM
	int			r[7];   // rate
	int    		a[7];  // acknowledgement
	double      p[7];	 // data power
	int	      req;    // requeests
	
	static double m;
	PS(struct Sim*); 
	void	init(double t, int sc, int id);
	void  setup(DMob* , double t);
	void  setRate(int c);
	void  move(Node<DMob*>*, int to, int from = -1); // default - from off
	void	run(double t, Base*);
};

#endif 
