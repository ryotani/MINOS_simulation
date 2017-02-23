#ifndef ExN03SteppingAction_h
#define ExN03SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "RecorderBase.hh"

class ExN03DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03SteppingAction : public G4UserSteppingAction
{
public:
  ExN03SteppingAction(ExN03DetectorConstruction*, RecorderBase*);
  virtual ~ExN03SteppingAction();

  void UserSteppingAction(const G4Step*);
    
private:
  ExN03DetectorConstruction* detector;
  RecorderBase *records;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
