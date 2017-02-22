#ifndef ExN03PhysicsList_h
#define ExN03PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03PhysicsList: public G4VUserPhysicsList
{
public:
  ExN03PhysicsList();
  virtual ~ExN03PhysicsList();

  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
 
  void SetCuts();
   
private:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();

  // these methods Construct physics processes and register them
  void ConstructDecay();
  void ConstructEM();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif



