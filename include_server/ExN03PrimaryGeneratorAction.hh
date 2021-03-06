#ifndef ExN03PrimaryGeneratorAction_h
#define ExN03PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "../lib/ExN03BeamIn.hh"

class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  ExN03PrimaryGeneratorAction();    
  virtual ~ExN03PrimaryGeneratorAction();

  void SetDefaultPrimaryParticle();
  void GeneratePrimaries(G4Event*);

  
  private:
  G4ParticleGun*                particleGun;	  //pointer a to G4  class

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
  G4int Z,A;
  ExN03BeamIn *BeamIn;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


