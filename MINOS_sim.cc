#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "ExN03DetectorConstruction.hh"
#include "ExN03PrimaryGeneratorAction.hh"
#include "ExN03RunAction.hh"
#include "ExN03EventAction.hh"
#include "ExN03SteppingAction.hh"
#include "ExN03PhysicsList.hh"
#include "ExN03TrackingAction.hh"
#include "G4EmCalculator.hh"
#include "QGSP_BERT.hh"
#include "LBE.hh"
#include "QGSP_INCLXX.hh"
#include "MinosModularPhysicsList.hh"
#include "ExN03ROOTuple.hh"
#include "RecorderBase.hh"
  
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#if defined(G4UI_USE_TCSH)
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#elif defined(G4UI_USE_XM)
#include "G4UIXm.hh"
#elif defined(G4UI_USE_WIN32)
#include "G4UIWin32.hh"
#elif defined(G4UI_USE_QT)
#include "G4UIQt.hh"
#include "G4Qt.hh"
#else
#include "G4UIterminal.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
G4cout << "Beginning of programm..." << endl;  
// Choose the Random engine
  //
  
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  //CLHEP::HepRandom::setTheSeed(111);
  
  // Construct the default run manager
  //
  
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  
  ExN03DetectorConstruction* detector = new ExN03DetectorConstruction;
  runManager->SetUserInitialization(detector);

  //runManager->SetUserInitialization(new ExN03PhysicsList);    
  MinosModularPhysicsList *physics = new MinosModularPhysicsList;
  //ExN03PhysicsList *physics = new    ExN03PhysicsList;
  runManager->SetUserInitialization(physics);   

  G4VUserPrimaryGeneratorAction* gen_action = 
                          new ExN03PrimaryGeneratorAction();
  runManager->SetUserAction(gen_action);

  //
  //For data analysis, creating ROOT file, DE, position,...
  RecorderBase * myRecords = new ExN03ROOTuple(detector); 

  ExN03RunAction* run_action = new ExN03RunAction(myRecords);  
  runManager->SetUserAction(run_action);

  ExN03EventAction* event_action = new ExN03EventAction(myRecords);
  runManager->SetUserAction(event_action);
  
  G4UserTrackingAction* tracking_action = new ExN03TrackingAction(myRecords);
  runManager->SetUserAction(tracking_action);
  
  G4UserSteppingAction* stepping_action =
                    new ExN03SteppingAction(detector, myRecords);
  runManager->SetUserAction(stepping_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
G4cout<<" ###SG4 kernel initialized "<<G4endl;
  
#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();      
  
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);    
    }
  else           // interactive mode : define visualization UI terminal
    {
      G4UIsession* session = 0;
#if defined(G4UI_USE_TCSH)
      session = new G4UIterminal(new G4UItcsh);      
#elif defined(G4UI_USE_XM)
      session = new G4UIXm(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_WIN32)
      session = new G4UIWin32();
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_QT)
      session = new G4UIQt(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#else
      session = new G4UIterminal();
#endif
#ifdef G4VIS_USE
      UI->ApplyCommand("/control/execute vis.mac");     
#endif
      session->SessionStart();
      delete session;
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
#ifdef G4VIS_USE
  delete visManager;
#endif                
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
