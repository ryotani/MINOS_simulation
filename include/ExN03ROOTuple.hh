#ifndef ExN03ROOTuple_h
#define ExN03ROOTuple_h 1
#include "RecorderBase.hh"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TMath.h"

#include "G4ios.hh"


#include "globals.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"

//#include "/Users/acorsi/codes/MINOS_simulation/lib/ExN03BeamIn.hh"
//#include "/Users/acorsi/codes/MINOS_simulation/lib/ExN03Datai.hh"
#include "/home/local1/workspace/MINOS_simulation/lib/ExN03BeamIn.hh"
#include "/home/local1/workspace/MINOS_simulation/lib/ExN03Datai.hh"
#include <vector>
#include "G4VProcess.hh"

#define NXY 151


class DetectorConstruction;

class ExN03ROOTuple:public RecorderBase
{

public:
ExN03ROOTuple(ExN03DetectorConstruction*);
virtual ~ExN03ROOTuple();

void RecordBeginOfRun();
void RecordEndOfRun();
void RecordStepDEtpc(const G4Step *);
void RecordStepDEtrigger(const G4Step *);
void RecordStepDEtarget(const G4Step *);
void RecordStepDEchamber(const G4Step *);
void RecordStepDEwindow(const G4Step *);
void RecordBeginOfEvent();
void RecordEndOfEvent();
void RecordBeginOfTrack();
void RecordEndOfTrack(const G4Track *);
void RecordStepDEInnerRohacell(const G4Step *);
void RecordStepDEOuterRohacell(const G4Step *);
void RecordStepDEKapton(const G4Step *);
void WriteDEDXTable(G4int , G4int );
private:
ExN03BeamIn *BeamIn;
ExN03Datai* data;
TFile *rfile;
TTree* tree_tpc;
bool detection;
double Et_tar, Et_win, Et_ch, Et_InnerRohacell, Et_OuterRohacell, Et_Kapton, Et_tpc, Et_trigger;
bool detection_daughter;
bool detection_daughter_PP;
bool detection_daughter_P2P;
bool detection_daughter_P3P_PAlpha;
bool detection_daughter_capture;
bool several_reaction_in_target;
bool detection_neutron;
bool detection_2pi;
int neutron_number;
bool event_P2P;
int det2pi;
int detection_proton, nb_event;
int event, detection_alpha;
int part;
double z[1000], r[1000], Z[1000], A[1000];
bool b_vertex, reaction_in_target;
int nb_events_per_run, nb_runs;
ExN03DetectorConstruction* detector;
double ev_tot, ev_P2P, N0, ev_P2P_bad;

};

#endif

