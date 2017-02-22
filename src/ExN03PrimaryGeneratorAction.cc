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
    
      //3 body  
  /*  float Afrag1=0.53; float Afrag2=1; float Afrag3=float(A);
    Double_t kinenergy1_LAB=200;Double_t theta1_LAB=40*PI/180.;Double_t theta2_LAB=80*PI/180.;
    Double_t mpion=493; Double_t mom1_LAB=300;
    cout<<"reading K kin"<<endl;
    float Qvalue,gamma_temp ,beta_v,E1_LAB,phi,pv1_LAB[3],p1_LAB,pp0,E0,M0,M23,P2_R23,E2_R23,v2_R23,gamma2_R23,E23,P23;
    float gamma_23,v23,g,costheta12,costheta2prime,AA,BB,DD,p2_LAB_plus,p2_LAB_minus,p2_LAB,pv2_LAB[3];
    float E2_LAB;
    float E3_LAB;
    float pv3_LAB[3],pv3_CM[3];
    float p3_LAB;
    float p3_CM;
    float qpar;
    float qper;
    float q2;
    float p3_LABbis;
    float kinenergy2_LAB;
    float kinenergy3_LAB;
    float Esep;
    float theta3_LAB;
    
  while(0)
  {
    h_xsec->GetRandom2(mom1_LAB,theta1_LAB);
    theta2_LAB=theta1_LAB;
    kinenergy1_LAB=sqrt(mom1_LAB*mom1_LAB*1e6+mpion*mpion)-mpion;
    cout<<"coordinates from xsec file "<<mom1_LAB<<" "<<theta1_LAB<<endl;
    theta1_LAB=theta1_LAB*PI/180.;theta2_LAB=-theta2_LAB*PI/180.;
    // Total beam energy before the target
    //double energy_total = MeanEnergy;        
    //energy_pv=Reconstruct(energy_p,p0[2],"106Zr",mass,"LH2",true)*mass;
    Qvalue=2./931.494;
    gamma_temp = (931.494*A+MeanEnergy)/(A*931.494);            // Gamma before the target
    beta_v = sqrt(1.0 - 1.0/gamma_temp/gamma_temp); // Beta before the target
    E1_LAB=kinenergy1_LAB+Afrag1*931.494;

    phi=0; //to be fixed later
    //compute p1
    pv1_LAB[3]; 
    p1_LAB=sqrt(E1_LAB*E1_LAB-Afrag1*Afrag1*931.494*931.494);
    pv1_LAB[0]=p1_LAB*sin(theta1_LAB)*cos(phi);
    pv1_LAB[1]=p1_LAB*sin(theta1_LAB)*sin(phi);
    pv1_LAB[2]=p1_LAB*cos(theta1_LAB);
    //compute p2
    pp0=sqrt(MeanEnergy*MeanEnergy+2*931*A*MeanEnergy);
    //cout<<"pp0 "<<sqrt(gamma_temp*gamma_temp-1.)*A*931.494<<" "<<sqrt(MeanEnergy*MeanEnergy+2*931*MeanEnergy)<<endl;
    E0=MeanEnergy+(A+1)*931.494; //1 is the mass of the target
    M0=sqrt(E0*E0-pp0*pp0); 
//cout<<E1_LAB<<" "<<p1_LAB<<" "<<E0<<" "<<pp0<<endl;
    M23=sqrt(M0*M0+Afrag1*Afrag1*931.494*931.494-2*E1_LAB*E0+2*pp0*p1_LAB*cos(theta1_LAB))/931.494;
    P2_R23=sqrt(lambda(M23*M23,Afrag2*Afrag2,Afrag3*Afrag3))/(2*M23)*931.494;
//cout<< "lambda "<<lambda(M23*M23,Afrag2*Afrag2,Afrag3*Afrag3)<<" "<< M23<<" "<<  P2_R23<<endl;
    E2_R23=(M23*M23+Afrag2*Afrag2-Afrag3*Afrag3)/(2*M23)*931.494;
    v2_R23=P2_R23/E2_R23;
    gamma2_R23=E2_R23/(Afrag2*931.494);    
//    cout<<sqrt(1.0 - 1.0/gamma2_R23/gamma2_R23)<<endl; 
    E23=E0-E1_LAB;
    P23=sqrt(pp0*pp0+p1_LAB*p1_LAB-2*pp0*p1_LAB*cos(theta1_LAB));
    gamma_23=E23/(M23*931.494);
    v23=P23/(M23*gamma_23*931.494);
    g=v23/v2_R23;
   // float sintheta23=p1_LAB*sin(theta1_LAB)/P23;
   // float costheta23=(pp0-p1_LAB*cos(theta1_LAB))/P23;
    costheta12=cos(theta1_LAB)*cos(theta2_LAB)+sin(theta1_LAB)*sin(theta2_LAB)*cos(0);// 0=phi1_phi2
    costheta2prime=(pp0*cos(theta2_LAB)-p1_LAB*costheta12)/P23;
//cout<<"angles "<<costheta12<<" "<<costheta2prime<<endl;
    AA=gamma_23*(1-v23*v23*costheta2prime*costheta2prime);
    BB=g*costheta2prime;
    DD=gamma_23*gamma_23*(1-g*g)+g*g*pow(gamma_23/gamma2_R23,2.)*costheta2prime*costheta2prime;
        //if(D>0)
  /*{cout<<"D"<<endl;
  cout<<gamma2_R23<<" "<<gamma_23<<" "<<g<<" "<<P23<<" "<<v23<<" "<<v2_R23<<endl;
        cout<<gamma_23*gamma_23*(1-g*g)<<" "<<g*g*pow(gamma_23/gamma2_R23,2.)*costheta2prime*costheta2prime<<endl;
  cout<<AA<<" "<<BB<<" "<<DD<<endl;}
    p2_LAB_plus=P2_R23*(BB+sqrt(DD))/AA;
    p2_LAB_minus=P2_R23*(BB-sqrt(DD))/AA;
    
    
    if(p2_LAB_plus>0&&p2_LAB_minus>0)p2_LAB=min(p2_LAB_plus,p2_LAB_minus);
    if(p2_LAB_plus<0&&p2_LAB_minus>0)p2_LAB=p2_LAB_minus;
    if(p2_LAB_plus>0&&p2_LAB_minus<0)p2_LAB=p2_LAB_plus;
    
    pv2_LAB[0]=p2_LAB*sin(theta2_LAB)*cos(phi);
    pv2_LAB[1]=p2_LAB*sin(theta2_LAB)*sin(phi);
    pv2_LAB[2]=p2_LAB*cos(theta2_LAB);
    //cout<<"p1_LAB: "<<p1_LAB<<", p2_LAB: "<<p2_LAB_plus<<","<<p2_LAB_minus<<endl;
    E2_LAB=sqrt(p2_LAB*p2_LAB+Afrag2*Afrag2*931.494*931.494);
    E3_LAB=E0-E1_LAB-E2_LAB;
 
    pv3_LAB[0]=-pv1_LAB[0]-pv2_LAB[0];
    pv3_LAB[1]=-pv1_LAB[1]-pv2_LAB[1];
    pv3_LAB[2]=pp0-pv1_LAB[2]-pv2_LAB[2];
    p3_LAB=sqrt(pv3_LAB[0]*pv3_LAB[0]+pv3_LAB[1]*pv3_LAB[1]+pv3_LAB[2]*pv3_LAB[2]);
    pv3_CM[0]=pv3_LAB[0];
    pv3_CM[1]=pv3_LAB[1];
    pv3_CM[2]=-E3_LAB*gamma_temp*beta_v+pv3_LAB[2]*gamma_temp;
    p3_CM=sqrt(pv3_CM[0]*pv3_CM[0]+pv3_CM[1]*pv3_CM[1]+pv3_CM[2]*pv3_CM[2]);
    qpar=1/gamma_temp*(pv2_LAB[2]+pv1_LAB[2]-beta_v*gamma_temp*931.494);
    qper=pv2_LAB[0]+pv1_LAB[0];
    q2=qpar*qpar+qper*qper;
    //cout<<"q2 "<<q2/(2*mass_f3*931.494)<<"="<<p3_CM*p3_CM/(2*mass_f3*931.494)<<endl;
    p3_LABbis=sqrt(pp0*pp0+p1_LAB*p1_LAB+p2_LAB*p2_LAB-2*pp0*p1_LAB*cos(theta1_LAB)-2*pp0*p2_LAB*cos(theta2_LAB)+2*p1_LAB*p2_LAB*costheta12);
//cout<< p3_LAB<<"="<<p3_CM<<endl;   
//cout<<"check conservation"<<endl;
    cout<<E0<<"="<<E1_LAB+E2_LAB+E3_LAB<<endl;
    kinenergy2_LAB=E2_LAB-Afrag2*931.494;
    kinenergy3_LAB=E3_LAB-Afrag3*931.494;
    Esep=MeanEnergy/(AA*931)-gamma_temp*(kinenergy1_LAB+kinenergy2_LAB)-2*(gamma_temp-1)*931.494+beta_v*gamma_temp*(pv1_LAB[2]+pv2_LAB[2])-p3_CM*p3_CM/(2*Afrag2*931.494);
    //cout<<" Esep: "<<Esep<<"="<<MeanEnergy/(A*931)<<"+"<<-gamma_temp*(kinenergy1_LAB+kinenergy2_LAB)<<"+"<<-2*(gamma_temp-1)*931.494<<"+"<<beta_v*gamma_temp*(pv1_LAB[2]+pv2_LAB[2])<<"+"<<-p3_CM*p3_CM/(2*Afrag3*931.494)<<endl;
//cout<<"Details  pv1_LAB[2] "<<pv1_LAB[2]<<"  pv2_LAB[2] "<<pv2_LAB[2]<<" qpar "<<qpar<<" qper "<<qper<<endl;
    theta1_LAB=theta1_LAB*180./PI;
    theta2_LAB=theta2_LAB*180./PI;
    theta3_LAB=acos(pv3_LAB[2]/p3_LAB)*180./PI;
    if(kinenergy1_LAB>10 && kinenergy2_LAB>10 && kinenergy3_LAB>10) break;
  }
  */

  // random position of the beam
  G4double x0 = G4RandGauss::shoot(MeanX, SigmaX);
  G4double y0 = G4RandGauss::shoot(MeanY, SigmaY);
  //G4double x0 = 0.;
  //G4double y0 = 0.;
  G4double z0 = 0; //with INCL
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
  //G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  G4ParticleDefinition *particle = particleTable->GetIon(Z,A,0.);
  G4String str = particle->GetParticleName();
  //G4cerr << str << G4endl;
  particleGun->SetParticleDefinition(particle);  
  G4String str2 = particle->GetParticleType();
  //G4cerr << str2 << G4endl;
  	G4int charge = particle->GetPDGCharge();
  	G4int mass = particle->GetBaryonNumber();
  	//G4cout << "Charge = " << charge << " / Mass = " << mass << G4endl;
  	
    

  
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
  
  particleGun->GeneratePrimaryVertex(anEvent);
  //this function is called at the begining of event
  //cout<<"G4partab**********************"<<endl;


//h_atomicBackground->GetRandom2(atomicBGTheta,atomicBGEnery);
//to shoot protons at 90deg
/*

  //definition of beam
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //cout<<"G4partab**********************"<<endl;
  G4String particleName;
  G4ParticleDefinition *particle = particleTable->GetIon(3,9,0.);
  //G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  //G4ParticleDefinition *particle = particleTable->GetIon(1,1,0.);
  G4String str = particle->GetParticleName();
  //G4cerr << str << G4endl;
  particleGun->SetParticleDefinition(particle);  
  G4String str2 = particle->GetParticleType();
  //G4cerr << str2 << G4endl;
  	G4int charge = particle->GetPDGCharge();
  	G4int mass = particle->GetBaryonNumber();
  	//G4cout << "Charge = " << charge << " / Mass = " << mass << G4endl;
  G4double x0 = G4RandGauss::shoot(MeanX, SigmaX);
  G4double y0 = G4RandGauss::shoot(MeanY, SigmaY);
  //G4double x0 = 0.;
  //G4double y0 = 0.;
  //G4double z0 = CLHEP::RandFlat::shoot(1)*15 *cm;
  G4double z0 = 0;
  particleGun->SetParticlePosition(G4ThreeVector(x0*mm,y0*mm,z0*mm));  
  	
    
   // float cos_theta = 1-CLHEP::RandFlat::shoot(0.,2.);
  double costhetamin=0.99;
  double costhetamax=1;
  double costheta=costhetamax-2*(costhetamax-costhetamin)*CLHEP::RandFlat::shoot(1);
  //double theta=acos(costheta);
  double theta=0;
  double phitemp=CLHEP::RandFlat::shoot(6.28);//cout<<"***phi*** "<<phitemp<<endl;
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*sin(phitemp),sin(theta)*cos(phitemp),cos(theta)));
  
  // random energy of the beam
  G4double energy = G4RandGauss::shoot(MeanEnergy, SigmaEnergy);
  //G4double energy = MeanEnergy;
  particleGun->SetParticleEnergy(250*MeV);
  
  particleGun->GeneratePrimaryVertex(anEvent);
  //this function is called at the begining of event
  //cout<<"G4partab**********************"<<endl;*/
  
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

