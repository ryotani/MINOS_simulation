#ifndef ExN03RunAction_h
#define ExN03RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "RecorderBase.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

class ExN03RunAction : public G4UserRunAction
{
public:
  ExN03RunAction(RecorderBase *);
  virtual ~ExN03RunAction();

  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
    
private:
  RecorderBase *records;
     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

