// functions of Sim

#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "sim.h"
#include "global.h"
#include "net.h"
#include "mob.h"
#include "data.h"
#include "rl.h"
#include <iostream>

using namespace std;

Sim::Sim(int argc, char* argv[])
{
    int  i, ple = 2, k = 1, out = 0;
  char  arg[80], crfName[25], par[25];
  // call data and data load parameters
  double *spf;
  int disp = 0, load = 0;

  memset(this,0,sizeof(Sim));

  nR = 1;
  strcpy(arg,"");
  strcpy(crfName,"r/");
  while (k < argc) {
    if (*argv[k] == '-')   {
      switch(argv[k][1]) {
        case 's':
          type = atoi(argv[k + 1]);
          break;
        case 't': // simulation time;
          t = atoi(argv[k + 1]);
          break;
        case 'c':
        case 'd':
          strcpy(par,argv[k + 1]);
          open(par, argv[k][1]);
          // t can be assigned from c file name
          break;
        case 'o':
          switch(*argv[k + 1]) {
            case 'g' :
              trace = 1;
            case 'd' :
              disp = 1;
              break;
            default :
              strcat(strcat(crfName,argv[k + 1]),".dat");
          }
          out = 1;
          break;
        case 'a':  // admission
          strcpy(par,argv[k + 1]);
          parseAdPar(par);
          break;
        case 'n': // number of runs
          nR = atoi(argv[k + 1]);
          break;
        case 'm' :
          ple = atoi(argv[k + 1]);
          break;
        case 'v' :
          spf = parse(argv[k + 1],",");
          for (i = 0; i < (int) *spf; i++)
            Mob::spf[i] = spf[i + 1];
          delete spf;
          break;
        case 'u' :
          Mob::nonU = atof(argv[k + 1]);
          break;
        case 'l' : // so far only for voice
          load = atoi(argv[k + 1]);
          break;
        default:
          quit();
        }
     strcat(strcat(arg, argv[k]),argv[k+1]);
     k += 2;
    }
    else
      quit();
  }

  if (!t)
    quit();

  Mob::sim = RL::sim = this;

  create(&cd, nR);
  create(&cds, nR);

  if (dlf)  {
    create(&dl, nR);
    create(&dlv, nR);
  }

  setPar();

  if (load)
    genLoad(load);

  if (!disp) {
    if (!out)
     setName(crfName);
    crf = fopen(crfName,"wb");
    fwrite(arg,sizeof(char),80,crf);
   }

  if (!adN)
    adN = 1;
  mrun(ple);
}

Sim::~Sim()
{
  for (int n = 0; n < nR; n++) {
    delete cd[n];
    if (dl)
      delete[] dl[n];
  }
}

void Sim::mrun(int m)
{
  int k = 0;
  char*  lfiles[4] = {"l2", "l3", "l4","l5"};

    if (!m)
    for (k = 0; k < m; k++)
      run(lfiles[k]);
  else
    run(lfiles[m - 2]);

  if (crf)
    fclose(crf);
  if (dlf)
    fclose(dlf);
  if (rl || sla)  // learning rc
    getState();
}

void  Sim::run(char* m)
{

  int k, n, c, crs;
  char fname[10];

  strcpy(fname,"in/");
  strcat(strcat(fname,m),".dat");
  FILE* af = fopen(fname,"rb");
  fread(att, sizeof(double), AttSize, af);
  fclose(af);

  for (k = 0; k < adN; k++) {
    if (adp)
      setAdm(k);

    for (n = 0; n < nR; n++) {
        // srand((unsigned) 100*(n + 1));  // seed
        ran((n + 1)*100); // any better seed selection ?
        if (!cd[n] && cdf)
          cds[n] = get(cdf, cd + n);
        if (dlf && !dl[n]) {
          fread(dlv + n, sizeof(double), 1, dlf);
            get(dlf, dl + n);
        }

        if (rl)
          for (c = 0; c < 7; c++)
            rl[c].reset();

        if (sla)
          for (c = 0; c < 7; c++)
            sla[c].reset();

        net = new Net(this, n);
        crs = net->run();
        if (crf) {
          if (dlv)
              fwrite(dlv + n,sizeof(double),1,crf);
          fwrite(net->thr,sizeof(double),2,crf);
          fwrite(net->cr,sizeof(double), crs ,crf);
        }
        else
          count();
        delete net;
    }
  }
 }

void Sim::setPar()
{

  double w;  // chip rate and noise bandwidth

  if (type < 2)
     w = wIS95;
  else
    w = wUTRA;

  /* if (type == ACSIM)
    Mob::va = 1;
  else
    Mob::va = 0; */

  if (type != FCSIM && type != RCSIM)
    L = both;
  else
    L = downLink;

  tFr = tFrIS95;
  int nS = 1;
  switch (type) {
    case PCSIM:
    case ACSIM:    // AC without PC
      r = vR;
      break;
    case FCSIM:
        r = dR;
      tFr = tFrUTRA;
      break;
      case MRSIM:   // multi-rate AC
      r = mR;
      nS = 3;
      break;
    case RCSIM:   // multi-rate AC with RC
      r = mRC;
      break;
  }

  tDr = (int) floor(DrT/tFr);
  nu = Boltz*TinK*w;
  if (type != 4) { // all but rate control
    cir = new double[nS*4];
    for (int n = 0; n < nS; n++)
      for (int k = 0; k < 4; k++)
        cir[n*4 + k] = pow(10,Eb2Io[n][k])/w;
  }

  RL::sim = this;
  Mob::sim = this;

  if (type)
    nL = 1;  // number of pc loops within a frame
  else
    nL = (int) floor(tFr/tPC);

  Ncdp = type < 3 ? 2 : 3;
  Ncrp = type < 3 ? 1 : type - 1;
  // for 3 sc is stored as second and for 4 throughput is stored as third
  // read("data.m");
  nCh = maxNch;
}

void Sim::parseAdPar(char* admpar)
{
  char* tok;
  int k = 0;
  char** ad = new char*[maxP];

  tok = strtok(admpar,":");   // parse admission
  while (tok) {
    ad[k++] = strcpy(new char[25],tok);
    tok = strtok(NULL,":");
  }
  adN = k;

  adp = new double*[adN];
  rcp = new double*[adN];
  char* rcpar;
  for (k = 0; k < adN; k++) {
    if (rcpar = strchr(ad[k],'_')) {
      rcp[k] = parse(rcpar + 1,",");
      ad[k] = strtok(ad[k],"_");
    }
    else
      rcp[k] = 0;
    adp[k] = parse(ad[k],",");
  }
}

int Sim::get(FILE* f, double** d)
{
  double val[1];

  fread(val, sizeof(double), 1, f);
  int s = (int) floor(*val);
  *d = new double[s];
  fread(*d, sizeof(double), s, f);
  return s;
}

double* Sim::parse(char* str, const char* term)
{
  char  *tok;
  double* out;
  int n = 1;

  if (str) {
    out = new double[maxP];
    // max expected number of parameters
    tok = strtok(str,term);
    while (tok)
    {
      out[n++] = atof(tok);
      tok = strtok(NULL,term);
    }
    out[0] = n - 1;
  }
  return out;
}

void Sim::setName(char* crfName)
{
   struct tm   *stime;
   long   ltime;
   char date[10];

   time(&ltime);
   stime = gmtime(&ltime);
   strftime(date,8,"%d%H%M",stime);
   strcat(strcat(crfName,date),".dat");
}

void Sim::open(char* fn, char arg)
{
  char fName[25];
  FILE* f;

  strcpy(fName,"in/");
  strcat(strcat(fName,fn),".dat");
  f = fopen(fName,"rb");
  if (!f)
    quit(fName);
  else {
    if (arg == 'c') {
      cdf = f;
      if (!t)
        t = atoi(strtok(fn,"_"));
      // based on -t coming before -c
    }
    else {
      dlf = f;
      /* int dlp[2];
      dls = 1;
      char* tok = strtok(fn,"_");
      for (int i = 0; i < 2; i++) {
        dlp[i] = atoi(tok);
        tok = strtok(NULL,"_");
        dls *= dlp[i];
      } */
      nDu = atoi(strtok(fn,"_"));
      // DMob::nR = dlp[1];
      // dls *= 2;
    }
  }
}

void Sim::quit(char* file)
{
  if (!file) {
    cout << "mandat parameters in mandat order -s (simtype) -t (simtime) -c (call load file) or -l (load)\n";
    cout << "optional -a (admission) -d (data load file) -o (output: g(raphics) d(isplay)) \n";
    cout << " -v (velosity) -u (nonuniformity) -n (runs) -m (ple) \n";
    cout << "Remember spaces ! -s 2 \n";
  }
  else
    cout << "File" << file <<  " not found \n";
  exit(1);
}

void Sim::setAdm(int n)
{

  int c;
  int adpN = (int) adp[n][0];
    double* ad = adp[n] + 1;

  if (adpN < 3) {
    adm = adpN;
    if (adm) {
      if (adm == PA) {
        *adt = pow(10, *ad)*nu;
        adt[1] = pilot*ad[1];
      }
      else {
        *adt = *ad;
        adm = *ad > 0 ? NA : QA;
      }
    }
  }

  else {
    adm = RLA;
    RL::tau = *ad; // first argument - tau
    RL::nP = (int) ad[1];
    RL::nL = (int) ad[2];
    int sg; // state generalisation
    if (adpN > 3)
      sg = (int) ad[3];
    else
      sg = 1;

    if (adN > 4)
      RL::gamma = ad[4];
    if (adN > 5)
      RL::alpha = ad[5];
    /* if (RL::psp)
      delete RL::psp;
    RL::psp = genLogSp(RL::pMax, RL::nP, 2); */

    *RL::lp =  (RL::nP - 1)/log10(*RL::pMax);
    RL::lp[1] = (RL::nP - 1)/log10(RL::pMax[1]/pilot);

    if (rl)
       delete[] rl;

    create(&rl, 7);  // nullset rl

    for (c = 0; c < 7; c++)
      rl[c].init(c, sg);

  }

  if (rcp[n]) {
    if (*rcp[n] > 1)  {
      if (type == RCSIM)
        RL::rTau = rcp[n][1];
      else {
        SLA::g = rcp[n][1];
        SLA::voc = rcp[n][2];
        SLA::lp = (SLA::nP - 1)/log10(pMax/pMin);
        if (sla)
          delete[] sla;
        sla = new SLA[7];
      }
      rc = LRC;
    }
    else {
      if (rcp[n][1] < 1) {
        PS::m = rcp[n][1];
        rc = PM;  // power margin
      }
      else
        rc = (int) rcp[n][1];
    }
  }
}

double* Sim::genLogSp(double* p, int N, int K)
{
  double min[2];
  double max[2];
  if (K > 1) {
    min[UL] = 1; // for PEF is for ffs
    min[DL] = pilot;
    max[UL] = p[UL];
    max[DL] = pilot*p[DL];
  }
  else { // tau or FCSIM
    *min = *p;
    *max = p[1];
  }

  double* v = new double[N*K];
  double dv[2];
  for (int k = 0; k < K; k++) {
    dv[k] = log10(max[k]/min[k])/N;
    for (int n = 0; n < N; n++)  {
      v[k*N + n] = min[k]*pow(10,(n + 1)*dv[k]);
      if (K > 1 && !k)  // fix-up !!!
        v[n] = nu*pow(10,v[n] - 1);
    }
  }
  return v;
}

void Sim::getState()
{

  int arS;
  double *ar;

  if (adm == RLA) {
    arS = nL*RL::nL*RL::nP;
    ar = new double[arS];
  }

  int i;
  /* if (rl->tet[0][0]) {
    ar[TET] = mxCreateCellMatrix(1,2);
    mxArray* tetA[2];
    int size[2] = {P, L};
    int dimT[3] = {7,2,0};
    for (i = 0; i < 2; i++) {
      dimT[2] = size[i];
      tetA[i] = mxCreateNumericArray(3,dimT,mxDOUBLE_CLASS,mxREAL);
      mxSetCell(ar[TET],i,tetA[i]);
    }
  }
  if (ef) {
    int dimA[] = {7,2,EF::nA};
    ar[Q] = mxCreateNumericArray(3,dimA,mxDOUBLE_CLASS,mxREAL);
    ar[A] = mxCreateCellMatrix(7,2);
  } */

  if (adm == RLA) {
    rl[3].state(ar);
    write("stat/rl.dat",ar,arS);
    for (i = 0; i < 7; i++)
      rl[i].~RL();
    }
    /* if (ef) {
      ef[i].state(ar[Q], ar[A],i);
      ef[i].~EF();
    } */
    if (type == FCSIM) {
      int s;
      double* v = sla[3].storeState(&s);
      write("stat/sla.dat", v, s);
      delete[] sla;
    }
}

void Sim::write(char* fname, double* data, int size)
{

  FILE* f = fopen(fname,"wb");
  if (size) {
    fwrite(data, size, sizeof(double), f);
    delete data;
  }
  fclose(f);  // file is flashed if no data
}

void Sim::count()
{
  if (type == FCSIM && dlv)  {
    cout << dlv[0] << "\t" << *net->thr << "\t" << net->thr[1];
    cout << "\t" << *net->thr/dlv[0] << "\n";
  }
  else
    cout <<   *net->thr << "\t" << net->thr[1] << "\n";

  double* cr = net->cr;
  double v, bl, dr, tl, out, pr; // premium traffic
  bl = dr = tl = out = pr = 0;
  int rate = 1;
  int i = 0;
  while ((v = cr[i]) != -3) {
     if (Ncrp > 1)
       rate = r[(int) cr[i + 1]];
     tl += rate;
     if (v == DR)
       dr += *w*rate;
     else if (v == BL)
       bl += rate;
     else if (v > qT)
       out += pow(*w,w[1]*v)*rate;
     if (Ncrp > 2)
       pr += rate*cr[i + 2]*w[2];
     i += Ncrp;
  }

  bl /= tl;
  dr /= tl;
  out /= tl;
  cout <<  bl << "\t";
  cout << dr << "\t";
  cout << out << "\t";
  if (Ncrp == 3)  {
    pr /= tl;
    cout << pr << "\t";
  }
  cout << dr + bl + out + pr << "\n";
}

void Sim::genLoad(int load)
// so far for voice, in the future will be extended to multiservice
{
   ran(1234);  // seed
   double tC = Mu[0]/(7*load); // interval for one call
   int S = (int) ceil(1.1*Ncdp*t/tC); // 1.1 to allow for higher load
   double* d = new double[S];
   double tA = 0; // arrival and release time
   d[0] = 0;
   int i = 0;
   do {
      d[i + 1]  = d[i] - Mu[0]*log(ran()); // tRelease
     // if (Ncdp > 2) generate service class from traffic matrix
    d[i + Ncdp] = d[i] - tC*log(ran()); // tArrival
    i += Ncdp;
   }  while (d[i] < t && i < S);
   cds[0] = i;
   cd[0] = new double[i];
   memcpy(cd[0], d, i*sizeof(double));
   delete d;
}

double ran(long seed)
{
  // seed can't be 0
  static long i;
  long k;
  
  if (seed) 
    i = seed;  // initialisation
  
  k = i/IQ;
  i = IA*(i - k*IQ) - IR*k;
  if (i < 0)
    i += IM;
  return  AM*i;
}
