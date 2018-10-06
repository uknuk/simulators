#include "data.h"
#include "net.h"
#include "heap.h"
#include <memory.h>
#include <stdlib.h>
#include "sim.h"
#include "stat.h"

const double IPpkt = 10; // IP packet in Kbit
const double df[2] = {0,1.57}; // 0 and Pi for sin and cos

double PS::m = 0;
double SLA::a = 0.1, SLA::g = 0.25, SLA::rMin = 16.0, SLA::lp = 0, SLA::voc = 0;
int SLA::nP = nPowStates, SLA::nR = nRates;

PS::PS(struct Sim* _sim)
{
  memset(this,0,sizeof(PS));
  sim = _sim;
  /* if (sim->rc < LRC) */
    for (int c = 0; c < 7; c++)  {
      r[c] = sim->r[7];  // 512 Kbps
      a[c] = 1;
    }
}

void PS::run(double t, Base* bs)
{
  Node<DMob*> *i, *s;
  DMob *m, *dm;
  int more;
  /* if (fc)
     setRate(); */
  for (int c = 0; c < 7; c++) {
    if (q[c].n) { // send packet to a scheduled mobile

      for (i = q[c].head; i; i = i->next)  {
        dm = i->item;
        dm->move(t);
        dm->shadow(t);
        dm->control(bs); // control channel
      }
      
      setRate(c);
      s = sm[c];
      m = s->item;
      m->transmit(r[c]);
      m->receive(bs,t);
      more = m->decode(this,c,bs);
      if (sim->sla && r[c])
        sim->sla[c].update(r[c], m->cp, a[c], sim->net->bs[c].o[1]);

      if (!more) {  // IP packet sent
        if (s->next) //  schedule a new mobile
          sm[c] = s->next;
        else
          sm[c] = q[c].head;
        
        if (!m->dataRest) {
          q[c].remove(s);
          if (!q[c].n)
            sm[c] = pm[c] = 0;
          
          delete s;
          if (!m->state(t))
            delete m->ev;
        }
        
        else if (m->c != c) {
          move(s,m->c,c);
          pm[c] = 0;
        }
        // handoff for sm if IP was sent and it won't be released and c was changed in shadow()
      }
      else
        p[c] = 0;
    }
  }
}

void PS::setRate(int c)
{

  // if (sim->sla && sim->sla[c].up)

  double pv = sim->net->bs[c].p[1] - p[c];  // voice power in previous frame

  
  int rc = sim->rc;

  if (rc == PM) {
        /* if (pm[c] != sm[c]) {
      pv -= pilot;
      double thr = sim->adt[1] - pilot;
      rate = sim->r[7]*(thr - pv)/thr;
      if (rate < 0)
        r[c] = 0;
      else
        r[c] = (int) ceil(rate/sim->r[6])*sim->r[6];
    } 
     else { */
    
    if (!r[c])
      r[c] = sim->r[1]; // control rate
    pv -= pilot;
    double thr = sim->adt[1] - pilot; // threshold
    double res = (thr - pv)/thr; // residual capacity
    double rate;
    
    if (!p[c])
       rate = sim->r[7]*res;
    else
      rate = r[c]*res*m*pMax/p[c];

    if (rate > sim->r[7])
      r[c] = sim->r[7];
    else if (rate < 0)
      r[c] = 0;
    else
      r[c] = (int) ceil(rate/sim->r[6])*sim->r[6];
  }
  else {
    if (a[c] && /* sim->net->bs[c].o[1]  && */ r[c] < sim->r[7] &&
      ((!rc || ((rc == PL || rc == PLS) && pv < sim->adt[1])) ||
      (rc == LRC && sim->sla[c].select(r[c], sm[c]->item->cp)) ))
         r[c] += sim->r[6];
    
    if (rc == PLS && pv >= sim->adt[1])  // power limit stop
      r[c] = 0;
    
    if (!a[c] && /* sim->net->bs[c].o[1]  &&  */ r[c] >= sim->r[6])
      r[c] -= sim->r[6];
    

  }

}

void PS::setup(DMob* m, double t)
{
  // if (m->tSt)    // not the first request
  //  m->jump(t);
  m->state(t);
  m->shadow(); // to get c
  Node<DMob*>* node = new Node<DMob*>(m);
  move(node,m->c);

}

void PS::move(Node<DMob*>* node, int to, int fr)
{
  if (fr >= 0) {   // -1 on new request
    q[fr].remove(node);
    if (!q[fr].n)
      sm[fr] = 0;
    node->next = node->prev = 0;
  }
  if (!q[to].n)
    sm[to] = node;
  q[to].put(node);
}

/***************************** DMob ********************************/

int DMob::init(int i, double* d, Heap* heap)
{
  setup(UDS);
  id = i;
  nR =  (int) -*d;
  dl = new double[nR][2];
  memcpy(dl, d + 1, 2*nR*sizeof(double));
  ev = new Event(dl[0][0], REQ, this); // ev->id not needed
  heap->put(ev);
  return nR;
}

void DMob::transmit(int rate)
{
  
  if (rate) {
    cr[1] = rate;
    pdu = cr[1]*sim->tFr;
    if (pdu > pktRest) {   // adjust rate to pdu size
      cr[1] = (int) ceil(pktRest/(r[MIN]*sim->tFr))*r[MIN];
      pdu = pktRest;
    }
  }
  else {
    cr[1] = sim->r[1]; // control rate
    pdu = 0;
  }
  setSIR();
}

int DMob::decode(PS* ps, int c0, Base* bs)
{
     Net* net = sim->net;
  // c0 - old cell, m->c can be changed in shadow
  bs[c0].p[1] += p[1];

  ps->a[c0] = (int) floor(1 - g[1]);   // acknowledge

  if (sim->trace && c0 == 3)
    net->traceRate(pdu*g[1]);

  ps->p[c0] = p[1];

  
  if (pdu) {
    if (!g[1]) {
      pktRest -= pdu;
      if (!pktRest) {
        net->thr[0] += pkt;
        dataRest -= pkt;
        if (dataRest)   // take a new packet
          pkt = pktRest = IPpkt < dataRest ? IPpkt : dataRest;
        /* else
          net->thr[1] += (data/(t - tSt) - net->thr[1])/++ps->req;  */
        return 0;
      }
    }
      else
      net->thr[1] += pdu;
  }

  return 1;
}

int DMob::state(double t)
{

  if (req == nR)
    return 0;

  if (st == ON) {
    double tOff = dl[req][0];
    if (tOff && t + tOff < sim->t) {
      ev->t = tOff;
      sim->net->heap->put(ev);
      st = OFF;
    }
    else
      return 0;
  }
  else {
    data = dataRest = dl[req++][1];
    st = ON;
    pkt = pktRest = IPpkt;
    }
  tSt = t;
  return 1;
}

void DMob::jump(double t)
{

  double s = (t - tSt)*sp;
  s = s*ran();
  double angle = ran();
  int x;
  for (int k = 0; k < 2; k++) {
    x = loc[NEW][k] + (int) floor(s*sin(angle + df[k]));
    x = x < attS[k] - 1 ? x : attS[k] - 1;
    x = x >= 0 ? x : 0;
    loc[NEW][k] =  loc[OLD][k] =  x;
  }
}

void DMob::control(Base* bs)
{
  
  double I = 0, Pt; 
  
  Net* net = sim->net;

  int j;
  
  for (j = 0; j < 7; j++) {
    Pt =  j == c ? of*(net->bs[j].p[DL] - cp) : net->bs[j].p[DL];
    I += Pt/sl[j];
  }

  cp = cSIR*I*sl[j];
  cp = cp > pMax ? pMax : cp;
  bs[c].p[1] += cp;
}


/****************************  SLA **************************************/

SLA::SLA() : p(0), r(0)
{
  q = new double[nP][nRates];
  for (int i = 0; i < nP; i++)
    for (int k = 0; k < nR; k++)
       q[i][k] = 0;

}

void SLA::update(int rate, double pow, int ack, double out)
{
  
  getState(rate, pow);
  ack = ack > 0 ? 1 : -1;
  double rew = ack - voc*out;
  q[p][r] += a*(rew - q[p][r]);  

}

int SLA::select(int rate, double pow)
{
  getState(rate, pow);
  int up = 1;
  if (q[p][r + 1] < 0)
    up =  ran() > 1 - exp(q[p][r + 1]/g);
   return up;
}



  // find max for current state
  /* double max = q[p1][0];
  int m = 0; // action with max
  for (int i = 1; i < nR; i++)
    if (q[p1][i] > max) {
      max = q[p1][i];
      m = i;
    }

  // update action value
  if (ps->r[c]) {
    double rew = ps->r[c] - bs->o[1]*bs->l;
    q[p][r] +=  a*(rew + g*max - q[p][r]);
  }
  // select action
  // if (m < nR - 1 && rand()/(double) RAND_MAX < 0.1)
  //  m += 1;

  r = m;
  p = p1;
  return r + 1;

}*/

void SLA::getState(int rate, double pow)
{
  p = (int) floor(lp*log10(pow/pMin));
  if (p > nP - 1)
      p = nP - 1;
  if (p < 0) // when pow = 1
      p = 0;
  r = (int) ceil(rate/rMin) - 1;
}

double*  SLA::storeState(int* s)
{
  *s = nP*nR;
  double* v = new double[*s];

  for (int i = 0; i < nP; i++)
    for (int k = 0; k < nR; k++)
      v[i*nR + k] = q[i][k];
  return v;
}
