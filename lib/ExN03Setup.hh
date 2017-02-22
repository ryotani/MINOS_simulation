#ifndef ExN03Setup_h
#define ExN03Setup_h 1

#include <cstdlib>
#include <fstream>
#include <iostream>
#include "TObject.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>


using namespace std;

class ExN03Setup:public TObject
{
public:

//initial conditions
double TargetRadius, TargetLength, ChamberInnerRadius, ChamberThickness;
double ChamberLength, InnerRohacellThickness, KaptonThickness, OuterRohacellThickness;
double TPCRadiusExt, WindowThickness;

public:
ExN03Setup();
~ExN03Setup();
void ReadConfigurationFile(string);

ClassDef(ExN03Setup,1)
};

#endif

