#include "stat.h"
#include "mob.h"
#include "net.h"
#include "sim.h"

Stat::Stat (int _cN, int _nFr, int sType) : cN(_cN), nFr(_nFr), tl(0), lN(0)
{
	nL = nR = 2;
	if (sType == FCSIM || sType == RCSIM) {
		nL = 1;
		create(&t, nFr);
	}
	else
		t = NULL;
	l = new double[cN][lS];
	r = new double[nFr*nL*nR];
	for (int k = 0; k < 2; k++)
		list[k] = new List<double>;

}

void Stat::radio(int fr, Base* bs, int link)
{
	
	int fptr = fr*nL*nR;
	int lnk, lptr;

	for (int k = 0; k < nL; k++) {
		lptr = fptr + k*nR;
		lnk = link ? link : k;  // fix-up doesn't work for uplink !
		r[lptr] = bs[3].p[lnk];
		r[lptr + 1] = bs[3].fer[lnk];
	}
}
	
void Stat::load(Base* bs, Mob* m)
{
	l[lN][0] = bs[3].l;
	tl += m->b;
	l[lN][1] = tl;
	for (int k = 0; k < 2; k++)
		l[lN][2 + k]  = bs[m->c].out[k];
	lN++;
}



void Stat::write(Sim* sim)
{
	double* data;
	int n, i;
	char* fname[2] = {"stat/block.dat", "stat/drop.dat" };
	// sim->write("stat/load.dat", (double*) l,lS*cN);
	sim->write("stat/radio.dat", r, nL*nFr*nR);
	if (t)
		sim->write("stat/rate.dat", t, nFr);
	for (int k = 0; k < 2; k++)  {
		if (n = list[k]->n) {
			data = new double[n];
			i = 0;
			for (Node<double>* node = list[k]->head; node; node = node->next)
				data[i] = node->item;
		}
		sim->write(fname[k], data, n);
	}
 }



