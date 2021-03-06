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
#include "G4EmCalculator.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

using namespace CLHEP;

ExN03ROOTuple::ExN03ROOTuple(ExN03DetectorConstruction* det)
:detector(det)
{
data = new ExN03Datai();
BeamIn = new ExN03BeamIn();
static TROOT rootbase("simple","MINOS TPC");
nb_runs = 40;//15
nb_events_per_run = 10000;
}

ExN03ROOTuple::~ExN03ROOTuple(){
rfile->Close();
delete rfile;
}

//en debut de run on definit la composition du fichier root
void ExN03ROOTuple::RecordBeginOfRun()
{
//tree_tpc = new TTree("tpc","A ROOT tree for MINOS data");
//tree_tpc->Branch("tpc","ExN03Datai",&data);
//tree_tpc->Branch("tpc","ExN03Datai_eventbyevent_new",&data);
G4cout << "Beginning of run: TTree creation..." << G4endl;
nb_event=0;
ev_P2P=0;
ev_tot=0;
N0=0;
ev_P2P_bad=0;
}

void ExN03ROOTuple::RecordEndOfRun()
{
double density = 6.41*pow(10.,23);
//G4cout << G4endl << "Total cross section: " << ev_tot/N0/density*pow(10.,27) << " mb" << G4endl << "Cross section (p,2p): " << ev_P2P/N0/density*pow(10.,27) << " mb" << G4endl << G4endl;
G4cout << G4endl << "Total cross section: " << ev_tot/N0/density*pow(10.,27) << " mb" << G4endl << "Cross section (p,p2n): " << ev_P2P/N0/density*pow(10.,27) << " mb" << G4endl << G4endl;
G4cout << G4endl << "Number of (p,p2n) events: " << ev_P2P << G4endl << "Number of (p,p2n) events followed by another reaction : " << ev_P2P_bad << G4endl << G4endl;
}

void ExN03ROOTuple::RecordBeginOfEvent()
{
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
neutron_number=0;
event_P2P = false;
several_reaction_in_target=false;
b_vertex=true;
reaction_in_target=true;
event = -1;
vertexcounter = 0;
part = 0;
char s[10];
for(int i=0; i<=nb_runs-1; i++)
{
	sprintf(s, "%d", i);
	if(nb_event == i*nb_events_per_run)
	{
		rfile=new TFile(("/data_tmp/acorsi/result/11Li/simu_modularphys" + string(s) + ".root").c_str(),"RECREATE");
		//rfile=new TFile(("/mnt/hgfs/F/100Sn_2/simu" + string(s) + ".root").c_str(),"RECREATE");
		tree_tpc = new TTree("tpc","A ROOT tree for MINOS data");
		tree_tpc->Branch("tpc","ExN03Datai",&data);
	}
}
}

void ExN03ROOTuple::RecordEndOfEvent()
{
	if(part==0) event = 0; //pas de reaction
	double vertex=0.;
	if(part!=0)
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
			if(z[ii] > -0.4 && (z[ii] < vertex - 0.01 || z[ii] > vertex + 0.01) && (z[ii] > 0. && z[ii] < (detector->GetTargetLength()+1.) && r[ii] > 0. && r[ii] < detector->GetTargetRadius())) {several_reaction_in_target = true;}	
		}
		for(int ii=0; ii<part; ii++)
		{
			if(z[ii] > -0.4 && !(z[ii] > 0. && z[ii] < (detector->GetTargetLength()) && r[ii] > 0. && r[ii] < detector->GetTargetRadius())) {reaction_in_target = false;}	
		}
	}
	if (!b_vertex && part!=0) event = 1; // double reaction

	if(b_vertex && part!=0 && reaction_in_target)
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
 		if((detection_daughter_PP) && (detection_proton==1) && (neutron_number==2)) event=11;
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
		if(zz[i]>0. && zz[i] < detector->GetTargetLength()) reac_target = true;//reaction in the target
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
	data->vertexcounter = vertexcounter;
	tree_tpc->Fill();
	//if(nb_event<50000) tree_tpc->Fill();
	//if (nb_event>=50000 && (event==11 || event==12)) tree_tpc->Fill();

	char s[10];
	for(int i=1; i<=nb_runs; i++)
	{
		sprintf(s, "%d", i);
		if(nb_event == i*nb_events_per_run-1)
		{
			rfile->Write();
			delete tree_tpc;
			rfile->Close();
			delete rfile;
		}
	}
	nb_event++;

	//WriteDEDXTable(BeamIn->A, BeamIn->Z);

}

void ExN03ROOTuple::RecordBeginOfTrack()
{
	detection = false;
	Et_tar=0.;
	Et_win=0.;
	Et_tpc=0.;
	Et_ch=0.;
	Et_InnerRohacell=0.;
	Et_OuterRohacell=0.;
	Et_Kapton=0.;
}

void ExN03ROOTuple::RecordEndOfTrack(const G4Track* track)
{
  	G4int charge = (int) (track->GetDefinition()->GetPDGCharge());
  	G4int mass = track->GetDefinition()->GetBaryonNumber(); 
  	//if(charge>10&&mass>40)
  	if(track->GetParentID()==0)vertexcounter++;
  	bool tst = false;
  	if(charge == (int) BeamIn->Z && track->GetVertexPosition().z()>-0.1) {detection_daughter_PP = true; tst=true;}
  	if(charge == (int) BeamIn->Z - 1) {detection_daughter_P2P = true; tst=true;}
  	if(charge == (int) BeamIn->Z - 2) {detection_daughter_P3P_PAlpha = true; tst=true;}
  	if(charge > (int) BeamIn->Z ) {detection_daughter_capture = true; tst=true;}
  	if(charge == 1 && mass == 1) {detection_proton = detection_proton + 1; tst=true;}
  	if(charge == 0 && mass == 1) {detection_neutron = true; neutron_number+=1; tst=true;} 	
  	if(charge == 2 && mass == 4) {detection_alpha = detection_alpha + 1; tst=true;}

    	G4double z0 = track->GetVertexPosition().z();
  	G4double x0 = track->GetVertexPosition().x();  
  	G4double y0 = track->GetVertexPosition().y();    	
    	
    	if(tst)
    	{
    		z[part] = z0;
	    	r[part] = sqrt(x0*x0+y0*y0);
    		Z[part] = charge;
    		A[part] = mass;
    		part++;
    	}
    	    	  	
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
    	data->energy0.push_back(energy);
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
    	data->Et_Kapton.push_back(Et_Kapton);
    	data->Et_InnerRohacell.push_back(Et_InnerRohacell);
    	data->Et_OuterRohacell.push_back(Et_OuterRohacell);
}

void ExN03ROOTuple::RecordStepDEchamber(const G4Step* s)
{

G4double		 edep = s->GetTotalEnergyDeposit();

if(edep>0.){
//G4cout << "Recording data step by step in the chamber..." << G4endl; 

G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

/*data->x_ch.push_back(x/mm);
data->y_ch.push_back(y/mm);
data->z_ch.push_back(z/mm);
data->e_ch.push_back(edep/eV);*/

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


Et_tar += edep/eV;
}
//else G4cout << "Recording data step by step in the target: energy too low..." << G4endl; 
}

void ExN03ROOTuple::RecordStepDEtpc(const G4Step* s)
{
G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
//G4cout << "Recording data step by step in the TPC..." << G4endl; 
G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_tpc.push_back(x/mm);
data->y_tpc.push_back(y/mm);
data->z_tpc.push_back(z/mm);
data->e_tpc.push_back(edep/eV);

data->Et_tpc_tot +=edep/eV;
Et_tpc += edep/eV;

detection = true;
}
//else G4cout << "Recording data step by step in the TPC: energy too low..." << G4endl; 
}

void ExN03ROOTuple::RecordStepDEwindow(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
//G4cout << "Recording data step by step in the Mylar window..." << G4endl; 
/*G4StepPoint* prepoint = s->GetPreStepPoint();
G4StepPoint* postpoint = s->GetPostStepPoint();

G4double x = (prepoint->GetPosition().x()+postpoint->GetPosition().x())*0.5;
G4double y = (prepoint->GetPosition().y()+postpoint->GetPosition().y())*0.5;
G4double z = (prepoint->GetPosition().z()+postpoint->GetPosition().z())*0.5;

data->x_win.push_back(x/mm);
data->y_win.push_back(y/mm);
data->z_win.push_back(z/mm);
data->e_win.push_back(edep/eV);

*/
Et_win += edep/eV;
}
//else G4cout << "Recording data step by step in the Mylar window: energy too low..." << G4endl; 
}
void ExN03ROOTuple::RecordStepDEInnerRohacell(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
Et_InnerRohacell += edep/eV;
}
}
void ExN03ROOTuple::RecordStepDEOuterRohacell(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
Et_OuterRohacell += edep/eV;
}
}
void ExN03ROOTuple::RecordStepDEKapton(const G4Step* s)
{

G4double edep = s->GetTotalEnergyDeposit();

if(edep>0.){
Et_Kapton += edep/eV;
}
}
void ExN03ROOTuple::WriteDEDXTable(G4int mass, G4int charge)
{
       //   Opening hte output file
       ofstream dedxfile;
       G4String filename="/home/gpfs/manip/mnt23/structurenucleaire/acorsi/MINOS_beam/Inputs/dedx_C_104Sn.txt";
       dedxfile.open(filename)   ;
       G4Material* material = detector->GetTargetMaterial();
       G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
       G4ParticleDefinition* Particle= particleTable->GetIon(charge, mass, 0.);
       G4EmCalculator emCalculator;
       G4double Emin = 0.*mass;
       G4double Emax = 400.*mass;
       for (G4double E=Emin; E < Emax; E+=(Emax-Emin)/10000.)
       {
               G4double dedx = emCalculator.ComputeTotalDEDX(E, Particle, material);
               dedxfile << E/MeV/mass << "\t" << dedx/(MeV/micrometer) << G4endl ;
       }
       dedxfile.close();
}

