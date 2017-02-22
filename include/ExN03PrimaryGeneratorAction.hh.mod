#ifndef ExN03PrimaryGeneratorAction_h
#define ExN03PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "/Users/acorsi/codes/MINOS_simulation/lib/ExN03Setup.hh"
#include "/Users/acorsi/codes/MINOS_simulation/lib/ExN03BeamIn.hh"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TRandom.h"

class G4ParticleGun;
class G4Event;
class ExN03DetectorConstruction;
//class ExN03PrimaryGeneratorMessenger;
//class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  ExN03PrimaryGeneratorAction();    
  virtual ~ExN03PrimaryGeneratorAction();

  void SetDefaultPrimaryParticle();
  void GeneratePrimaries(G4Event*);

  void SetMeanX(G4double);
  void SetSigmaX(G4double);
  void SetMeanY(G4double);
  void SetSigmaY(G4double);
  void SetZ0(G4double);
  void SetMeanEnergy(G4double);
  void SetSigmaEnergy(G4double);
  void SetMeanMomentumX(G4double);
  void SetSigmaMomentumX(G4double);
  void SetMeanMomentumY(G4double);
  void SetSigmaMomentumY(G4double);
  void SetMomentumZ0(G4double);
  void SetZ(G4int);  
  void SetA(G4int);
  float lambda(float, float, float);
  TFile *f;
  G4double theta;
  G4double phi;
  TH2F *h_xsec;  
  TFile* xsecFile;
  //TH3F* h_xsec;  
  private:
  G4ParticleGun*                particleGun;	  //pointer a to G4  class
  ExN03DetectorConstruction*    ExN03Detector;    //pointer to the geometry
  //G4ParticleDefinition *particle;
   
  //ExN03PrimaryGeneratorMessenger* gunMessenger;   //messenger of this class

  G4double MeanX;
  G4double SigmaX;
  G4double MeanY;
  G4double SigmaY;
  G4double Z0;
  G4double MomentumZ0;
  G4double MeanEnergy;
  G4double SigmaEnergy;
  G4double MeanMomentumX;
  G4double SigmaMomentumX;
  G4double MeanMomentumY;
  G4double SigmaMomentumY;
  G4double momentum;
  G4int Z,A;
  ExN03BeamIn *BeamIn;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


