#ifndef ExN03BeamIn_h
#define ExN03BeamIn_h 1

#include <cstdlib>
#include <fstream>
#include <iostream>
#include "TObject.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>


using namespace std;

class ExN03BeamIn:public TObject
{
public:

//initial conditions
double MeanX, SigmaX, MeanY, SigmaY, Z0, MeanMomentumX , SigmaMomentumX, MeanMomentumY, SigmaMomentumY, MomentumZ0, MeanEnergy, SigmaEnergy;
int Z, A;
public:
ExN03BeamIn();
~ExN03BeamIn();
void ReadConfigurationFile(string);

ClassDef(ExN03BeamIn,1)
};

#endif

