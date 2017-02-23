#include "ExN03DetectorConstruction.hh"
#include "ExN03ROOTuple.hh"
#include "RecorderBase.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4ios.hh"

#include "G4SDManager.hh"
#include "G4VVisManager.hh"
#include <vector>
using namespace CLHEP;
ExN03ROOTuple::ExN03ROOTuple(ExN03DetectorConstruction* det):detector(det){
BeamIn = new ExN03BeamIn();
data = new ExN03Datai();
static TROOT rootbase("simple","MINOS TPC");
det2pi=0;
}

ExN03ROOTuple::~ExN03ROOTuple(){
//rfile->Close();
//delete rfile;
}

//en debut de run on definit la composition du fichier root
void ExN03ROOTuple::RecordBeginOfRun()
{
  //rfile=new TFile("/Users/acorsi/codes/MINOS_simulation/result/actar_p1bar_cut1keV_r20mm_tpc.root","RECREATE");
  rfile=new TFile("/home/local1/workspace/MINOS_simulation/result/test.root","RECREATE");
  tree_tpc = new TTree("tpc","A ROOT tree for MINOS data");
  tree_tpc->Branch("tpc","ExN03Datai",&data);
  G4cout << "Beginning of run: TTree creation..." << G4endl;
}

void ExN03ROOTuple::RecordEndOfRun()
{
		rfile->Write();
		delete tree_tpc;
		rfile->Close();
		delete rfile;
		G4cout << "End of run: TTree destruction..." << G4endl;
}

void ExN03ROOTuple::RecordBeginOfEvent()
{
//G4cout << "RecordBeginOfEvent..." << G4endl;
data->ClearEvent();
N0 += 1;
data->ClearEvent();
detection_proton = 0;
detection_alpha=0;
detection_daughter = false;
detection_daughter_P2P = false;
detection_daughter_PP = false;
detection_daughter_P3P_PAlpha = false;
detection_daughter_capture = false;
detection_neutron = false;
detection_2pi = false;
neutron_number=0;
event_P2P = false;
several_reaction_in_target=false;
b_vertex=true;
reaction_in_target=true;
event = 0;
part = 0;


}

void ExN03ROOTuple::RecordEndOfEvent()
{
//G4cout << "RecordEndOfEvent..." << G4endl;
if(part==1) event = 0; //pas de reaction
double vertex;
if(part!=1)
{
	for(int ii=0; ii<part; ii++)
	{
		if(z[ii] > -0.4) {vertex = z[ii]; ii=part;}	
	}
	for(int ii=0; ii<part; ii++)
	{
		if(z[ii] > -0.4 && (z[ii] < vertex - 0.01 || z[ii] > vertex + 0.01)) {b_vertex = false;}	
	}
	for(int ii=0; ii<part; ii++)
	{
		if(z[ii] > -0.4 && (z[ii] < vertex - 0.01 || z[ii] > vertex + 0.01) && (z[ii] > 0. && z[ii] < (detector->GetTargetLength()*2.+1.) && r[ii] > 0. && r[ii] < detector->GetTargetRadius())) {several_reaction_in_target = true;}	
	}
	for(int ii=0; ii<part; ii++)
	{
		if(z[ii] > -0.4 && !(z[ii] > 0. && z[ii] < (detector->GetTargetLength()*2.) && r[ii] > 0. && r[ii] < detector->GetTargetRadius())) {reaction_in_target = false;}	
	}
}
if (!b_vertex && part!=1) event = 1; // double reaction

if(b_vertex && part!=1 && reaction_in_target)
{
 if((detection_daughter_P2P) && (detection_proton==2) && (!detection_neutron)) event=5;
 if((detection_daughter_P3P_PAlpha) && (detection_proton==3) && (!detection_neutron) && (detection_alpha==0)) event=7;
 if((detection_daughter_PP) && (detection_proton==1) && (!detection_neutron)) event=3;
 if((detection_daughter_capture)) event=2;
 if((detection_daughter_P2P) && (detection_proton==2) && (detection_neutron)) event=6;
 if((detection_daughter_P3P_PAlpha) && (detection_proton==3) && (detection_neutron) && (detection_alpha==0)) event=8;
 if((detection_daughter_PP) && (detection_proton==1) && (detection_neutron)) event=4;
 if((detection_daughter_P3P_PAlpha) && (detection_proton==1) && (detection_alpha==1) && (!detection_neutron)) event=10;
 if((detection_daughter_P3P_PAlpha) && (detection_proton==1) && (detection_alpha==1) && (detection_neutron)) event=9;
 //if((detection_daughter_PP) && (detection_proton==1) && (neutron_number==2)) event=11;
 if(det2pi) event=11;
}


// add by LA 12/02/13
int nb_reac=0, nb_part[1000], ZZ[1000][1000], AA[1000][1000];//number of reactions, number of particles, charge and mass per reaction
double zz[1000];// Z position of vertex per reaction
for(int i=0; i<1000; i++)
{
	nb_part[i]=0;
}
for(int ii=0; ii<part; ii++)//loop over tracked particles
{
	if(z[ii] > -0.4)//does not consider incoming beam
	{
		if(nb_reac==0) {zz[nb_reac] = z[ii]; ZZ[nb_reac][nb_part[nb_reac]] = Z[ii]; AA[nb_reac][nb_part[nb_reac]] = A[ii]; nb_part[nb_reac]+=1; nb_reac++;}//for the first particle encountered (see the second following comment)
		else
		{
			bool new_reac=true;
			for(int k=0; k<nb_reac; k++)
			{
				if(!(z[ii] < zz[k] - 0.01 || z[ii] > zz[k] + 0.01)) {ZZ[k][nb_part[k]]=Z[ii]; AA[k][nb_part[k]] = A[ii]; nb_part[k]+=1; new_reac=false;} //if not new vertex, select the reaction number (vertex) among vertex already existing for particle ii, increment charge, mass, number of particles for the vertex
			}
			if(new_reac) {zz[nb_reac] = z[ii]; ZZ[nb_reac][nb_part[nb_reac]] = Z[ii]; AA[nb_reac][nb_part[nb_reac]] = A[ii]; nb_part[nb_reac]+=1; nb_reac++;}//if new vertex, new Z poision of vertex, increment charge, mass of the particle, the number of particles for this vertex, increment the number of reaction 
		}
	}
}
bool reac_entrance_window = false;//reaction in Mylar entrance window
bool reac_target = false;//reaction in the target
bool detection_residual = false;
int number_of_protons = 0;
int number_of_neutrons = 0;
for(int i=0; i<nb_reac; i++)//loop over vertex
{
	if(zz[i]<0.) reac_entrance_window = true;//reaction in Mylar entrance window
	if(zz[i]>0. && zz[i] < detector->GetTargetLength()*2.) reac_target = true;//reaction in the target
}
if(reac_entrance_window) N0 -= 1;//event not considered (only reactions in target are interesting)
if(!(reac_entrance_window) && reac_target)//for first reaction in the target
{
	ev_tot += 1;
	int n_min = TMath::LocMin(nb_reac,zz);//select the first reaction which occurs at the minimum Z pos. of vertex
	//if(nb_part[n_min] == 3)//for (p,2p) events: residual + 2p as the first reaction
	if(nb_part[n_min] == 4)
	{
		for(int i=0; i<nb_part[n_min]; i++)
		{
			//if(ZZ[n_min][i] == BeamIn->Z-1 && AA[n_min][i] == BeamIn->A-1) detection_residual = true;//check residual
			if(ZZ[n_min][i] == 1 && AA[n_min][i] == 1) number_of_protons += 1;//check protons
			if(ZZ[n_min][i] == BeamIn->Z && AA[n_min][i] <= BeamIn->A-2) detection_residual = true;
			if(ZZ[n_min][i] == 0 && AA[n_min][i] == 1) number_of_neutrons += 1;//check neutrons
		}
	}
	//if(detection_residual && number_of_protons==2) ev_P2P += 1;//select (p,2p) events
	if(detection_residual && number_of_protons==1 && number_of_neutrons==2) ev_P2P += 1;//select (p,p2n) events
	if(detection_residual && number_of_protons==1 && number_of_neutrons==2 && several_reaction_in_target) {ev_P2P_bad += 1; event=12;}//select (p,p2n) events followed by another reaction in the target
}



data->event = event;








tree_tpc->Fill();


}

void ExN03ROOTuple::RecordBeginOfTrack()
{
//G4cout << "RecordBeginOfTrack..." << G4endl;
	detection = false;
	Et_tar=0.;
	Et_win=0.;
	Et_tpc=0.;
	Et_ch=0.;
	Et_InnerRohacell=0.;
	Et_OuterRohacell=0.;
	Et_Kapton=0.;
	Et_trigger=0.;

}

void ExN03ROOTuple::RecordEndOfTrack(const G4Track* track)
{
//G4cout << "RecordEndOfTrack..." << G4endl;
  	G4int charge = (int) (track->GetDefinition()->GetPDGCharge());
  	G4int mass = track->GetDefinition()->GetBaryonNumber(); 
  	
    G4double z0 = track->GetVertexPosition().z();
  	G4double x0 = track->GetVertexPosition().x();  
  	G4double y0 = track->GetVertexPosition().y();    	

//  	data->trackID.push_back(track->GetTrackID());
//  	data->parentID.push_back(track->GetParentID());  
    	
    	    	  	
    	G4double x0Momentum = track->GetVertexMomentumDirection().x();  
  		G4double y0Momentum = track->GetVertexMomentumDirection().y();
    	G4double z0Momentum = track->GetVertexMomentumDirection().z();
    	G4ThreeVector Momentum = track->GetVertexMomentumDirection();      	
    	
    	G4double theta = acos(z0Momentum)*180./pi;
    	G4double phi=0.;
	if(x0Momentum>0 && y0Momentum>0) phi=180./pi*atan(y0Momentum/x0Momentum);
	if(x0Momentum<0 && y0Momentum>0) phi=180.-180./pi*atan(-1.*y0Momentum/x0Momentum);
	if(x0Momentum<0 && y0Momentum<0) phi=180.+180./pi*atan(y0Momentum/x0Momentum);
	if(x0Momentum>0 && y0Momentum<0) phi=360.-180./pi*atan(-1.*y0Momentum/x0Momentum);
	if(theta<0.000001) phi = -1.;   
    	   	
    	G4double energy = track->GetVertexKineticEnergy();
    	
   	data->theta0.push_back(theta);
    data->phi0.push_back(phi);
    data->energy0.push_back(energy/keV);
    data->x0.push_back(x0);
    data->y0.push_back(y0);
    data->z0.push_back(z0);
    data->Z.push_back(charge);
    data->A.push_back(mass);
    data->detection.push_back(detection);
    data->Et_tar.push_back(Et_tar);
    data->Et_ch.push_back(Et_ch);
    data->Et_win.push_back(Et_win);
    data->Et_tpc.push_back(Et_tpc);
    data->Et_trigger.push_back(Et_trigger);
    data->Et_Kapton.push_back(Et_Kapton);
    data->Et_InnerRohacell.push_back(Et_InnerRohacell);
    data->Et_OuterRohacell.push_back(Et_OuterRohacell);
  	

	bool tst = false;
  	if(charge == (int) BeamIn->Z && track->GetVertexPosition().z()>-0.1) {detection_daughter_PP = true; tst=true;}
  	if(charge == (int) BeamIn->Z - 1) {detection_daughter_P2P = true; tst=true;}
  	if(charge == (int) BeamIn->Z - 2) {detection_daughter_P3P_PAlpha = true; tst=true;}
  	if(charge > (int) BeamIn->Z ) {detection_daughter_capture = true; tst=true;}
  	if(charge == 1 && mass == 1) {detection_proton = detection_proton + 1; tst=true;}
  	if(charge == 0 && mass == 1) {detection_neutron = true; neutron_number+=1; tst=true;} 	
  	if(charge == 2 && mass == 4) {detection_alpha = detection_alpha + 1; tst=true;}
	//cout<<"eot: "<<detection_daughter_P2P<<" "<<detection_proton<<endl;
    	if(tst)
    	{
    		z[part] = z0;
	    	r[part] = sqrt(x0*x0+y0*y0);
    		Z[part] = charge;
    		A[part] = mass;
    		part++;
    	}
    	    	  	
}


void ExN03ROOTuple::RecordStepDEchamber(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
//G4cout << "Recording data step by step in the chamber..." << G4endl; 

/*G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_ch.push_back(x/mm);
data->y_ch.push_back(y/mm);
data->z_ch.push_back(z/mm);
data->e_ch.push_back(edep/eV);

*/
Et_ch += edep/eV;
}
//else G4cout << "Recording data step by step in the chamber: energy too low..." << G4endl; 
}



void ExN03ROOTuple::RecordStepDEtarget(const G4Step* s)
{
G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
//G4cout << "Recording data step by step in the target..." << G4endl; 
G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_tar.push_back(x/mm);
data->y_tar.push_back(y/mm);
data->z_tar.push_back(z/mm);
data->e_tar.push_back(edep/eV);
//data->trackID.push_back(s->GetTrack()->GetTrackID());
//data->parentID.push_back(s->GetTrack()->GetParentID());  
Et_tar += edep/eV;
//data->Et_tpc_tot +=edep/eV;

}
//else G4cout << "Recording data step by step in the target: energy too low..." << G4endl; 
}

void ExN03ROOTuple::RecordStepDEtpc(const G4Step* s)
{
G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){

G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_tpc.push_back(x/mm);
data->y_tpc.push_back(y/mm);
data->z_tpc.push_back(z/mm);
data->e_tpc.push_back(edep/eV);
//data->trackID.push_back(s->GetTrack()->GetTrackID());
//data->parentID.push_back(s->GetTrack()->GetParentID());  
data->Et_tpc_tot +=edep/eV;
Et_tpc += edep/eV;
data->trackID.push_back(s->GetTrack()->GetTrackID());
data->parentID.push_back(s->GetTrack()->GetParentID());  

detection = true;
}
//else G4cout << "Recording data step by step in the TPC: energy too low..." << G4endl; 
}

void ExN03ROOTuple::RecordStepDEtrigger(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
//G4cout << "Recording data step by step in the Mylar window..." << G4endl; 
G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_trigger.push_back(x/mm);
data->y_trigger.push_back(y/mm);
data->z_trigger.push_back(z/mm);
data->e_trigger.push_back(edep/eV);

Et_trigger += edep/eV;
}
//else G4cout << "Recording data step by step in the Mylar window: energy too low..." << G4endl; 
}

void ExN03ROOTuple::RecordStepDEwindow(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
//G4cout << "Recording data step by step in the Mylar window..." << G4endl; 
G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_win.push_back(x/mm);
data->y_win.push_back(y/mm);
data->z_win.push_back(z/mm);
data->e_win.push_back(edep/eV);

Et_win += edep/eV;
}
//else G4cout << "Recording data step by step in the Mylar window: energy too low..." << G4endl; 
}
void ExN03ROOTuple::RecordStepDEInnerRohacell(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
/*G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_InRoh.push_back(x/mm);
data->y_InRoh.push_back(y/mm);
data->z_InRoh.push_back(z/mm);
data->e_InRoh.push_back(edep/eV);
*/
Et_InnerRohacell += edep/eV;
}
}
void ExN03ROOTuple::RecordStepDEOuterRohacell(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
/*G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_OutRoh.push_back(x/mm);
data->y_OutRoh.push_back(y/mm);
data->z_OutRoh.push_back(z/mm);
data->e_OutRoh.push_back(edep/eV);
*/
Et_OuterRohacell += edep/eV;
}
}
void ExN03ROOTuple::RecordStepDEKapton(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
/*G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_Kap.push_back(x/mm);
data->y_Kap.push_back(y/mm);
data->z_Kap.push_back(z/mm);
data->e_Kap.push_back(edep/eV);
*/
Et_Kapton += edep/eV;
}
}


