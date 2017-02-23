#include "ExN03SteppingAction.hh"

#include "ExN03DetectorConstruction.hh"

#include "G4Step.hh"


ExN03SteppingAction::ExN03SteppingAction(ExN03DetectorConstruction* det,
                                         RecorderBase *r)
:detector(det)					 
{ records = r;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03SteppingAction::~ExN03SteppingAction()
{ records  =NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // register track position and energy loss stap by step for protons 
  
  G4VPhysicalVolume* volume 
  = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();      
  G4Track *track = aStep->GetTrack(); 

  if (volume == detector->GetTarget() && ((track->GetDefinition()->GetParticleType() == "nucleus") || (track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron")))
  {
  	if(records!=NULL){records->RecordStepDEtarget(aStep);}
  }
  if (volume == detector->GetTPC())
  {
  	if(records!=NULL){records->RecordStepDEtpc(aStep);}
  }
  if ((volume == detector->GetWindow0() || volume == detector->GetWindow1() || volume == detector->GetWindow2()) && ((track->GetDefinition()->GetParticleType() == "nucleus") || (track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron")))
  {
  	if(records!=NULL){records->RecordStepDEwindow(aStep);}
  }
  if (volume == detector->GetChamber() && ((track->GetDefinition()->GetParticleType() == "nucleus") || (track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron")))
  {
  	if(records!=NULL){records->RecordStepDEchamber(aStep);}
  }
  if (volume == detector->GetInnerRohacell() && ((track->GetDefinition()->GetParticleType() == "nucleus") || (track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron")))
  {
  	if(records!=NULL){records->RecordStepDEInnerRohacell(aStep);}
  }
  if (volume == detector->GetOuterRohacell() && ((track->GetDefinition()->GetParticleType() == "nucleus") || (track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron")))
  {
  	if(records!=NULL){records->RecordStepDEOuterRohacell(aStep);}
  }
  if (volume == detector->GetKapton() && ((track->GetDefinition()->GetParticleType() == "nucleus") || (track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron")))
  {
  	if(records!=NULL){records->RecordStepDEKapton(aStep);}
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
