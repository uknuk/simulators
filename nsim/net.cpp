// function file for Net

#include "global.h"
#include "net.h"
#include "sim.h"
#include "mob.h"
#include "heap.h"
// #include "ef.h"
#include "stat.h"
#include "data.h"
#include "rl.h"




// rates in Kbits



const int EndOfData = -3;





Net::Net(Sim* _sim, int n)
{
	memset(this, 0, sizeof(Net));
  	sim = _sim;
	cds = sim->cds[n];
	int nC = cds/sim->Ncdp;	
	if (sim->trace)
		stat = new Stat(nC, (int) floor(sim->t/sim->tFr), sim->type);

  	create(&cr,nC*sim->Ncrp + 1);
	// cr allocated to maximum number of possible calls

	create(&m, sim->nCh);
	
	heap = new Heap(sim->nCh);
	
	if (sim->type == FCSIM && sim->nDu)	{
		ps = new PS(sim);
		create(&dm, sim->nDu);
		int ptr = 0;
		for (int i = 0; i < sim->nDu; i++) 	
			ptr += 2*dm[i].init(i, sim->dl[n] + ptr, heap) + 1;
	}
																									
	cd = sim->cd[n];
}

int Net::run()
{
	int ptr = 0;		
	while (t < sim->t) {
		if (ptr < cds && cd[ptr] < t)	{
			arrival(cd + ptr);
			ptr += sim->Ncdp;
		}
		if (heap->last > -1)
			handle();
 		radio();
	}

 for (int i = 0; i < nM; i++)
		if (m[i].id)
 			release(m + i);

	int crs = rcN*sim->Ncrp;
  	cr[crs++] = EndOfData;  // for Matlab
	return crs;
}


Net::~Net()
{
	delete m;

   if (stat) {
		stat->write(sim);
		delete stat;
	}
	
	while (*heap->ev) 
		delete heap->get();
	delete heap;
	
	if (ps)	{
		delete ps;
		delete dm;
	}
}

void Net::handle()
{
	Event* ev;
	while (*heap->ev && (*heap->ev)->t < t)	 {
		 ev = heap->get();
		 switch (ev->type) {
			case REL:
				if (ev->id == ev->gen->id)		// if call drop
					release(ev->gen);
				delete ev;
				break;
			case REQ: 
				ps->setup((DMob *) ev->gen, t);
				break;
			case VATR:
				if (ev->id == ev->gen->id)
					ev->gen->state(t,ev);
				else
					delete ev;
				break;
		 }
	}
}

	
void Net::arrival(double* ccd) // current call data 
{

	int sc = 0;
	if (sim->Ncdp > 2)	// multiservice
		sc = (int) *(ccd + 2); // service class

	int ch;
	for (ch = 0; m[ch].id && ch < sim->nCh; ch++);
	if (ch == sim->nCh) {
		grow(m, sim->nCh, INC);
		printf("Channels added \n");
		sim->nCh += INC;
	}
		
	m[ch].setup(sc);
	if (admit(m + ch)) {
			if (ch == nM) // increment nM
				nM = ch + 1;
			m[ch].proceed(t, *(ccd + 1));
			if (stat)
				stat->load(bs, m + ch);
	}
	else 
		release(m + ch,BL);

}

int Net::admit(Mob* nm)
{ 
	int k;
	int ok = 1;
	int c = nm->c;
	
	switch (sim->adm) {
		case NA:
			ok = bs[c].l + nm->b < sim->adt[0];
			break;
		case PA:
			double pow;
			for (k = sim->L[0]; k < sim->L[1]; k++)
				if (ok) {
					if (sim->type == FCSIM && ps)	
						pow = bs[c].p[k] - ps->p[c];
					else
						pow = bs[c].p[k];
					ok = bs[c].p[k] < sim->adt[k];
				}
			break;
		case RLA:
			if ((sim->type < 3 && bs[c].n >= Nmin) || 
					(sim->type >= 3 && bs[c].l >= Nmin))
				ok = sim->rl[c].admit(nm);
			break;
		case QA:
			/* for (k = 0; k < 2; k++) {
			if (ok && bs[c].p[k] >= pMin[k] && bs[c].ao[k] > -*adp*FER[0])	*/
			if ((sim->type < 3 && bs[c].n >= Nmin) || 
					(sim->type >= 3 && bs[c].l >= Nmin))  
			  	 for (k = sim->L[0]; k < sim->L[1]; k++) {
					if (ok &&  bs[c].ao[k] /* bs[c].out[k] */ > -sim->adt[0])
							/* ok = rand()/(double) RAND_MAX > 
								1 - exp((1e-4 - bs[c].ao[k])/-*adp); */
							ok = 0;
					bs[c].ao[k] = 0;
				}
						// bs[c].t = 0;
	}
	return ok;
}

void Net::radio()
{
				  
	int j;
	Base tbs[7];
	double r = 0;
	
	for (int n = 0; n < sim->nL; n++) {
		memset(tbs,0,7*sizeof(Base));
		if (ps)
			ps->run(t,tbs);
			
		for (int i = 0; i < nM; i++) {
			if (m[i].id) {  
				if (!n) { /* frame begins */
					if (m[i].sp)
						m[i].shadow(t);
					m[i].setSIR();  // sirT is set for the whole frame
					/* if (!sim->type)	 
						m[i].g[1] = m[i].g[0] = 0;	*/
				}
				m[i].receive(tbs,t);
								
				if (n == sim->nL - 1)
					m[i].decode(t,i);
							
				for (j = 0; j < 7; j++)	
					if (m[i].al[j] <= m[i].hm) { /* if connected */
						m[i].transmit(n,tbs + j);
					   if (!n && j == 3 && sim->type == RCSIM && stat)	// is it needed somewhere else ?
							r += cr[1];				  
					}

				if ((m[i].io[0] == sim->tDr) || (m[i].io[1] == sim->tDr)) 
					release(m + i,DR);  // DROP = 1
						
			}  
		}  
			
		if (sim->type == RCSIM && stat)
			traceRate(r);
		
		transmit(n,tbs);
	
		if (!sim->type)  // pc
			t = fr*sim->tFr + (n + 1)*tPC;
		// delete(tbs);
	}
	if (stat)
		stat->radio(fr,bs,sim->L[0]); 
	fr++;
	if(sim->type)  // no pc
		t = fr*sim->tFr;
	
}

void Net::transmit(int n, Base* tbs)
{
	int j,k;
	
	for (j = 0; j < 7; j++)	 {
		
	
		if (bs[j].l) {
			if (!sim->type)	 // no pc
				bs[j].p[UL] = (bs[j].l - 1)*tbs[j].p[UL]/(double) bs[j].l + sim->nu;
			else		
				bs[j].p[UL] = tbs[j].p[UL] + sim->nu;
		}
		else
			bs[j].p[UL] = sim->nu; 
		
		bs[j].p[DL] = tbs[j].p[DL] + pilot;
			
		if (n == sim->nL - 1) { // frame ends
			bs[j].l = tbs[j].l;
			bs[j].n = tbs[j].n;
 			// load is recalculated because of handoff

	   	for (k = sim->L[0]; k < sim->L[1]; k++) {
			 	if (bs[j].l) {
					bs[j].o[k] = tbs[j].o[k]/bs[j].l;
					if (stat)
						bs[j].fer[k] = tbs[j].fer[k]/bs[j].l;
					if (bs[j].fer[k] > 1)
						printf("Bullshit");
					if (sim->adm == QA || sim->adm == RLA) {
						if (k == sim->L[0])
							++bs[j].t;	 
						bs[j].ao[k] += (bs[j].o[k] - bs[j].ao[k])/bs[j].t;
						}
					
					if (sim->type == RCSIM && sim->rc > DRC && bs[j].pMin[k] < pMax)
						rateUp(j,k);
		 		}
				else 
					bs[j].o[k] = 0;
		  }
		}
	}
}
		






void Net::rateUp(int c, int k)
{
	int ok = 0;
	switch (sim->rc) {
			case PRC:
				ok = !bs[c].o[k] && bs[c].p[k] < sim->adt[k];
				break;
			case LRC:
				ok = sim->rl[c].rateUp(k);
				break;
			default:
				ok = !bs[c].o[k];	  // CRC - ?
	}

	if (ok) {
		m[bs[c].chPmin[k]].cr[k] += sim->r[MIN];
		bs[c].pMin[k] = pMax;
	}
}

/* void Net::read(char* file)
{
	FILE* stream;
	if (!(stream = fopen(file,"r"))) 
		printf("data.m not found \n");
	else {
		char line[30];
		int i = 0;
		double* par = new double[10];
		while (!feof(stream)) {
			fgets(line,30,stream);
			if (*line != '%')
				par[i++] = atof(line);
		}

		Nmin = (int) par[0];
		RL::pMax[0] = par[1];
		RL::pMax[1] = par[2];
		RL::alpha = par[3];
		RL::gamma = par[4];
		RL::rT = par[5];
		RL::qT = par[6];
		delete par;
	}
}	 



void Net::countOut(int c, double t)
{
	bs[c].out[0] = bs[c].out[1] = 0;
	double dt = (t - bs[c].lt)/tFr;
	double out;
	int l;
	for (int i = 0; i < nM + 1 ; i++) 
		if (m[i] && m[i]->c == c)
			for (l = *L; l < L[1]; l++)	{
				out = (m[i]->ro[l] - m[i]->lro[l])/dt;
				// if (out > 0.05) 
				bs[c].out[l] += out*m[i]->b;
				m[i]->lro[l] = m[i]->ro[l];
			}
	bs[c].lt = t;
	for (l = *L; l < L[1]; l++)
		bs[c].out[l] /= bs[c].l;
}	 */
			




void Net::release(Mob* rm, int reason)
{
	double value, dt;
	double tFr = sim->tFr;
	if (reason < 0) {
		value = reason;
		if (stat)
			stat->list[reason + 2]->put(new Node<double>(t));
	}
	else {
		 dt = t - rm->cst;
		 dt = dt < tFr ? tFr : dt;
		 value = rm->co*tFr/dt;
	}
  	
	if (stat)
		stat->tl -= rm->b;
	
	int ptr = (rcN++)*sim->Ncrp;
	
	cr[ptr] = value;
	if (sim->type > 2)
		cr[ptr + 1] = rm->sc;

	if (sim->type > 3 && value >= 0) 	// only DL
		cr[ptr + 2] = rm->th[1]*tFr/dt;
	
	memset(rm, 0, sizeof(Mob));
}

      


/******************************************************************************/

void Base::init(Sim* sim)
{
	p[UL] = sim->nu;
	p[DL] = pilot;
	for (int k = 0; k < 4; k++)
			pMin[k] = pMax;
}

