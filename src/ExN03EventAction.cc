#include "ExN03EventAction.hh"


#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "RecorderBase.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03EventAction::ExN03EventAction(RecorderBase *r)
{
records = r;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03EventAction::~ExN03EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03EventAction::BeginOfEventAction(const G4Event* evt)
{  

  G4int evtNb = evt->GetEventID();
  
    //if(evtNb%1000==0)
    	{G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    CLHEP::HepRandom::showEngineStatus();}
    records->RecordBeginOfEvent();

 
   }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03EventAction::EndOfEventAction(const G4Event* evt)
{
     records->RecordEndOfEvent(); 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
