#ifndef ExN03EventAction_h
#define ExN03EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "RecorderBase.hh"

class ExN03RunAction;
class ExN03EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03EventAction : public G4UserEventAction
{
public:
  ExN03EventAction(RecorderBase *);
  virtual ~ExN03EventAction();

  void    BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
    
private:
   RecorderBase *records;
   
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
