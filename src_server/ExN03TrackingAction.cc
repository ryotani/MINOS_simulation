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
  if(records!=NULL && ((track->GetDefinition()->GetParticleType() == "nucleus") || (track->GetDefinition()->GetParticleName() == "proton") || (track->GetDefinition()->GetParticleName() == "neutron"))) records->RecordEndOfTrack(track);
}

