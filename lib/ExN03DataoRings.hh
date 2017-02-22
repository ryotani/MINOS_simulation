#ifndef ExN03DataoRings_h
#define ExN03DataoRings_h 1

#include "TObject.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "ExN03Datai.hh"


using namespace std;

class ExN03DataoRings:public TObject
{
public:
ExN03Datai datai;
vector<double> x_pad,y_pad,t_pad,dt;
vector<double> q_pad;
double Et_pad;
int nb_pads;

int NRings,NSegments;
double sigL,sigT,driftV,Ionis;
double TimeBinSize, ShapingTime;
double Threshold,Gain,NoiseRMS,Theta,FirstRing, LastRing;

public:
ExN03DataoRings();
~ExN03DataoRings();
void ClearEvent();
void ReadConfigurationFile(string);

ClassDef(ExN03DataoRings,1)
};

#endif

