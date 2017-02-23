#include "ExN03PrimaryGeneratorAction.hh"

#include "ExN03DetectorConstruction.hh"
//#include "ExN03PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiPropertiesTableAME03.hh"
//#include "G4NucleiPropertiesTable.hh"
#include "G4VIsotopeTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "G4IonTable.hh"
#include <G4KaonPlus.hh>
#include <G4Neutron.hh>
#include <G4Proton.hh>

#define PI 3.141592
using namespace CLHEP;
ExN03PrimaryGeneratorAction::ExN03PrimaryGeneratorAction()
{
  BeamIn = new ExN03BeamIn();
  cout<<"beam in"<<endl;
  //gunMessenger = new ExN03PrimaryGeneratorMessenger(this);

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  //particleGun  = new G4GeneralParticleSource();   
  SetDefaultPrimaryParticle();  
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03PrimaryGeneratorAction::~ExN03PrimaryGeneratorAction()
{
  delete particleGun;
  //delete gunMessenger;
}

void ExN03PrimaryGeneratorAction::SetDefaultPrimaryParticle()
{    
  cout<<"default particle characteristics"<<endl;
  
  MeanX=BeamIn->MeanX;    // position moyenne du faisceau en X
  SigmaX=BeamIn->SigmaX;   // écart-type de la position du faisceau en X
  MeanY=BeamIn->MeanY;    // position moyenne du faisceau en Y
  SigmaY=BeamIn->SigmaY;   // écart-type de la position du faisceau en Y
  Z0=BeamIn->Z0;
  MomentumZ0=BeamIn->MomentumZ0;
  MeanEnergy=BeamIn->MeanEnergy;
  SigmaEnergy=BeamIn->SigmaEnergy;
  MeanMomentumX=BeamIn->MeanMomentumX;    
  SigmaMomentumX=BeamIn->SigmaMomentumX;   
  MeanMomentumY=BeamIn->MeanMomentumY;    
  SigmaMomentumY=BeamIn->SigmaMomentumY;
  Z=BeamIn->Z;
  A=BeamIn->A;
 /* cout<<A<<" "<<Z<<" "<<"reading xsecfile"<<endl;
  xsecFile = new TFile("Inputs/xsec_K.root","READ"); //read TH3F: E1, t1, 
  xsecFile->GetObject("h",h_xsec);*/

  //here you should read xsec from Ogata-san 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // random position of the beam
  G4double x0 = G4RandGauss::shoot(MeanX, SigmaX);
  G4double y0 = G4RandGauss::shoot(MeanY, SigmaY);
  //G4double x0 = 0.;
  //G4double y0 = 0.;
  //G4double z0 = 0; //with INCL
  G4double z0 = Z0; // Z0=-.5 to avoid the vertex of beam inside ot the target.. 
  //G4double z0 = G4RandUniform::shoot()*10*cm; //with Ogata-san
  particleGun->SetParticlePosition(G4ThreeVector(x0*mm,y0*mm,z0*mm));  
  //if xsec from Ogata-san
  //sample 1 event (eg shoot integer random corresponding to index of vector/line)
  //calculate momentumX,momentumY,momentumZ, energy
  //shoot proton1, proton2
  //  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="proton"));//cout<<"primary gen "<<particle[i]<<" "<<*(&particle[i]<<" "<<particle.at(i)<<endl;
  //  particleGun -> SetParticleMomentumDirection(GG4ThreeVector(momentumX,momentumY,momentumZ));
  //  particleGun -> SetParticleEnergy(energy);
  //  particleGun->SetParticlePosition(G4ThreeVector(x0*mm,y0*mm,z0*mm));  
  //  particleGun->GeneratePrimaryVertex(anEvent); 
   //origial part gen
    
  //if xsec from INCL, shoot the beam and let it interact
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //cout<<"G4partab**********************"<<endl;
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  //G4ParticleDefinition *particle = particleTable->GetIon(Z,A,0.);
  G4String str = particle->GetParticleName();
  //G4cerr << str << G4endl;
  particleGun->SetParticleDefinition(particle);  
  G4String str2 = particle->GetParticleType();
  //G4cerr << str2 << G4endl;
  G4int charge = particle->GetPDGCharge();
  G4int mass = particle->GetBaryonNumber();
  G4cout << "Charge = " << charge << " / Mass = " << mass << G4endl;
  
  
  // random momentum of the beam 
  G4double momentumX = G4RandGauss::shoot(MeanMomentumX, SigmaMomentumX);
  G4double momentumY = G4RandGauss::shoot(MeanMomentumY, SigmaMomentumY);
  //G4double momentumX = 0.;
  //G4double momentumY = 0.;
  G4double momentumZ = MomentumZ0;  
  particleGun->SetParticleMomentumDirection(G4ThreeVector(momentumX,momentumY,momentumZ));
  
  // random energy of the beam
  G4double energy = G4RandGauss::shoot(MeanEnergy, SigmaEnergy);
  //G4double energy = MeanEnergy;
  particleGun->SetParticleEnergy(energy*MeV);
  
  //particleGun->GeneratePrimaryVertex(anEvent);
  //this function is called at the begining of event
  //cout<<"G4partab**********************"<<endl;



  //////////////// 2nd Particle ///////////////
  particleGun->SetParticlePosition(G4ThreeVector(x0*mm,y0*mm,z0*mm));  
  //if xsec from INCL, shoot the beam and let it interact
  G4ParticleTable* particleTable2 = G4ParticleTable::GetParticleTable();
  //cout<<"G4partab**********************"<<endl;
  G4String particleName2;
  //G4ParticleDefinition* particle2 = particleTable->FindParticle(particleName="proton");
  G4ParticleDefinition *particle2 = particleTable->GetIon(Z,A,0.);
  //G4String str = particle2->GetParticleName();
  //G4cerr << str << G4endl;
  particleGun->SetParticleDefinition(particle2);
  //G4String str2 = particle2->GetParticleType();
  //G4cerr << str2 << G4endl;
  G4int charge2 = particle2->GetPDGCharge();
  G4int mass2 = particle2->GetBaryonNumber();
  G4cout << "Charge = " << charge2 << " / Mass = " << mass2 << G4endl;
  
  
  // random momentum of the beam 
  //G4double momentumX2 = G4RandGauss::shoot(MeanMomentumX, SigmaMomentumX);
  //G4double momentumY2 = G4RandGauss::shoot(MeanMomentumY, SigmaMomentumY);
  G4double momentumX2 = 0.;
  G4double momentumY2 = 0.;
  G4double momentumZ2 = MomentumZ0;  
  particleGun->SetParticleMomentumDirection(G4ThreeVector(momentumX2,momentumY2,momentumZ2));
  
  // random energy of the beam
  G4double energy2 = G4RandGauss::shoot(MeanEnergy, SigmaEnergy);
  //G4double energy = MeanEnergy;
  particleGun->SetParticleEnergy(energy2*MeV);
  
  particleGun->GeneratePrimaryVertex(anEvent);
  //this function is called at the begining of event
  cout<<"G4partab**********************"<<endl;

  
}

void ExN03PrimaryGeneratorAction::SetMeanX (G4double val )  
{ MeanX = val;}

void ExN03PrimaryGeneratorAction::SetSigmaX (G4double val )  
{ SigmaX = val;}

void ExN03PrimaryGeneratorAction::SetMeanY (G4double val )  
{ MeanY = val;}

void ExN03PrimaryGeneratorAction::SetSigmaY (G4double val )  
{ SigmaY = val;}

void ExN03PrimaryGeneratorAction::SetZ0 (G4double val )  
{ Z0 = val;}

void ExN03PrimaryGeneratorAction::SetMeanEnergy (G4double val )  
{ MeanEnergy = val;}

void ExN03PrimaryGeneratorAction::SetSigmaEnergy (G4double val )  
{ SigmaEnergy = val;}

void ExN03PrimaryGeneratorAction::SetMeanMomentumX (G4double val )  
{ MeanMomentumX = val;}

void ExN03PrimaryGeneratorAction::SetSigmaMomentumX (G4double val )  
{ SigmaMomentumX = val;}

void ExN03PrimaryGeneratorAction::SetMeanMomentumY (G4double val )  
{ MeanMomentumY = val;}

void ExN03PrimaryGeneratorAction::SetSigmaMomentumY (G4double val )  
{ SigmaMomentumY = val;}

void ExN03PrimaryGeneratorAction::SetMomentumZ0 (G4double val )  
{ MomentumZ0 = val;}

void ExN03PrimaryGeneratorAction::SetZ (G4int val )  
{ Z = val;}

void ExN03PrimaryGeneratorAction::SetA (G4int val )  
{ A = val;}

float ExN03PrimaryGeneratorAction::lambda(float x, float y, float z)
{
  return x*x+y*y+z*z-2.*x*y-2.*y*z-2.*z*x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

