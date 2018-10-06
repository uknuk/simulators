#ifndef GLOBAL_H
#define GLOBAL_H

#include <math.h>



const double	pilot = 1.0;
const double	pMax = 0.25;

const double	tPC = 1e-3;
const double   DrT = 5;  	// time to drop
const double   FER[2] = {1e-2,1e-4};
const double	of = 0.4;
const double	Delta = pow(10,0.3);	 // handoff margin
const int		maxNch = 1000;
const double	wIS95 =  1.228e6;
const double	wUTRA =  3.84e6;
const double   w[2] = {10,3}; // default cost weights
const double   qT = 0.05; // call quality threshold

const int maxP = 6; // max number of ad parameters
const int AttSize = 455000;
const double defSPF[2] = {0.1, 0.75}; // default spf

const	int vR[] = {8,1};
const int mR[] = {8,64,144};
const int mRC[] = {8,16,32,128,64,256};
const int dR[] = {8,1,0,0,0,0,16,512};
const double Eb2Io[][4] = {{0.5,0.7,0.6,0.8},
		{0.4,0.6,0.5,0.7}, {0.3,0.5,0.4,0.6}};
// 0 or 1m/s: UL, DL 10 m/s: UL, DL

const double	Boltz = 1.23e-23;
const double	TinK = 290;

const double	tFrIS95 =  20e-3;
const double	tFrUTRA =  10e-3;

 
const int both[2] = {0,2};	  // both links
const int downLink[2] = {1,2};	  // downlink


enum SimType {PCSIM,ACSIM,FCSIM,MRSIM,RCSIM};
enum EventType {REL = 1, VATR, REQ}; 
// VATR (Voice Activity State Transition), 0 for cancelled event
enum Adm {NA = 1,PA,QA,PEF,RLA};
enum FailStatType {DR = -2, BL,BLOCK,DROP};
enum Links {UL,DL};
enum CallResStatType {CQ,CT,CRT = 1}; //Call quality, throughput, release time
enum RCtype {DRC = 1, PL = 1, CRC = 2, PLS = 2, PRC = 3, PM = 3, LRC};
// Distributed, Power Limit (flow control), Centralised, Power Limit Stop,
//  Power Margin, Power, Learning: RL or SLA

enum Service {PHONE,AUDIO,VIDEO,UDS};  // Unconstrained Delay Service

const int Nmin = 10;  // minimal mobile numbers for outage statistics

const double trR[6] = {1.1116 ,4.4817, 4.4817, 3.266, 3.266, 1.1116};
// cumulative state transition rates
const double stPr[6] = {0.2506, 0.4813, 0.4813, 0.3196, 0.3196, 0.2506};
// state transition probabilities
const int nTr[6][2] = {{1,3},{5,0},{0,5},{5,0},{0,5},{2,4}};
//	state transition matrix 
const double trT[6] = {0,0,0,0.2,0.2,0};  
// time when state transitions rates are zero 
const int bR[6][2] = {{0,1},{0,0},{0,0},{1,1},{1,1},{1,0}};
// rates for states: UL and DL, 0 - data, 1 control rate 
const double speed[2] = {1,10};
const double Poly[3] = {6.9732,   -0.7077,    0.0208};
// coefficients for uplink 1m/S Eb/Io (dB) = f(log2(r))
// fitted for values 5 4.5	4 3.5 3 2.6 2.3 2 for 2^(3 : 10)
// downlink + 2dB for real-time, +1 dB for best effort 
const double ln2 =  0.6931;
const int C0[2] = {80,75};
const int C[2] = {100,100};	// size of the hot spot
const double hsA = 0.1538; // 100^2/(260*250) area of hot spot



const double Mu[3] = {120,120,120};
// call mean times

const double up = 1.259;  // 1 dB
const double down = 0.794;
const double sirUp = pow(10,0.1); // outer loop
const double sirDown = 1/pow(10,1e-3); 


const int attS[2] = {260, 250};

// Fading constants
const int No = 8;
const int N = 34; // N = 2(2No+1)  
const double M = 17; /* M = <Xi^2> + <Xq^2> = 2(No+0.5) */
const double lambda = 0.16; /* for f = 2 GHz */
const double Pi = 3.14;

const double ff = 0.125; // forgetting factor for averaging

enum Rate {VOICE,DATA = 0,MIN = 0,MAX, CONTROL = 1} ;
enum  Location {OLD,NEW};
 

#endif
