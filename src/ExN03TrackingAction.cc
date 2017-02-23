#include "ExN03TrackingAction.hh"
#include "ExN03DetectorConstruction.hh"
#include "ExN03EventAction.hh"

ExN03TrackingAction::ExN03TrackingAction(RecorderBase *r)
{ 
   records = r;
}

ExN03TrackingAction::~ExN03TrackingAction()					 
{ 
   records = NULL;
}

void ExN03TrackingAction::PreUserTrackingAction(const G4Track *track)
{

  if(records!=NULL) records->RecordBeginOfTrack(); 

}

void ExN03TrackingAction::PostUserTrackingAction(const G4Track *track)
{
  //if(records!=NULL && (track->GetDefinition()->GetParticleName() == "e-"|| track->GetDefinition()->GetParticleType() == "nucleus"))
  	if(records!=NULL &&(track->GetDefinition()->GetParticleType() == "kaon") ||(track->GetDefinition()->GetParticleType() == "nucleus") || 
  		(track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron"))
  	records->RecordEndOfTrack(track);
  //cout<<"RecordEndOfTrack"<<endl;

 //if(records!=NULL && ((track->GetDefinition()->GetParticleType() == "pi+") || (track->GetDefinition()->GetParticleName() == "pi-"))) records->RecordEndOfTrack(track);
}

