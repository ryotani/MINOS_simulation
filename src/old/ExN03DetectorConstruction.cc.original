#include "ExN03DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4UserLimits.hh"
#include "MyMagneticField.hh"

#include "G4MagIntegratorStepper.hh"

#include "G4UserLimits.hh"
using namespace CLHEP;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorConstruction::ExN03DetectorConstruction()
{
  
  setup = new ExN03Setup();
  TargetLength = setup->TargetLength/2.*mm;
  TargetRadius=setup->TargetRadius*mm;
  WindowThickness=setup->WindowThickness/2.*mm;
  ChamberLength=setup->ChamberLength/2.*mm;
  ChamberInnerRadius=setup->ChamberInnerRadius*mm;
  TPCRadiusExt=setup->TPCRadiusExt*mm; 
  ChamberThickness =setup->ChamberThickness*mm;
  InnerRohacellThickness =setup->InnerRohacellThickness*mm;
  OuterRohacellThickness =setup->OuterRohacellThickness*mm;
  KaptonThickness =setup->KaptonThickness*mm;
  
  // materials

  DefineMaterials();
  SetTargetMaterial("LH2");
  SetChamberMaterial("Inox");
  SetTPCMaterial("LH2"); //SetTPCMaterial("mix"); 
  SetWindowMaterial("Mylar");  
  SetKaptonMaterial("Kapton");  
  SetInnerRohacellMaterial("Rohacell");
  SetOuterRohacellMaterial("Rohacell");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN03DetectorConstruction::~ExN03DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN03DetectorConstruction::ConstructMINOS()
{
  return Construct();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN03DetectorConstruction::DefineMaterials()
{ 
//This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  

G4int ncomponents, natoms;
G4double fractionmass;

//
// define Elements
//

G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
G4Element* F  = new G4Element("Fluorin"  ,symbol="F" , z= 9., a= 19.0*g/mole);
G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
G4Element* Ar = new G4Element("Argon", symbol="Ar", z=18, a=39.948*g/mole);
G4Element* Fe = new G4Element("Fer",symbol="Fe" , z= 26., a= 55.9*g/mole);
G4Element* Cr = new G4Element("Chrome",symbol="Cr" , z= 24., a= 52.*g/mole);


new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);
new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
new G4Material("Lead"     , z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
new G4Material("Silicium", z=14., a=28.09*g/mole, density=2.330*g/cm3);
new G4Material("Titanium", z=22., a=47.87*g/mole, density=4.510*g/cm3);

G4Material* iso = new G4Material("isobutane", density=0.002506*g/cm3, ncomponents=2);
iso->AddElement(C, natoms=4);
iso->AddElement(H, natoms=10);
G4Material* CF4 = 
new G4Material("CF4", density= 0.0036586*g/cm3, ncomponents=2);
CF4->AddElement(C, natoms=1);
CF4->AddElement(F, natoms=4);
// overwrite computed meanExcitationEnergy with ICRU recommended value 
CF4->GetIonisation()->SetMeanExcitationEnergy(20.0*eV);
G4Material* mix = 
new G4Material("mix", density= 0.0019836*g/cm3, ncomponents=3);
mix->AddMaterial(CF4, fractionmass=15.*perCent);
mix->AddMaterial(iso, fractionmass=3.*perCent);
mix->AddElement(Ar, fractionmass=82.*perCent);
// overwrite computed meanExcitationEnergy with ICRU recommended value 
mix->GetIonisation()->SetMeanExcitationEnergy(25.0*eV);

G4Material* Ar_CF4_95_5 = 
new G4Material("Ar_CF4_95_5", density= 0.0017611*g/cm3, ncomponents=2);
Ar_CF4_95_5->AddMaterial(CF4, fractionmass=5.*perCent);
Ar_CF4_95_5->AddElement(Ar, fractionmass=95.*perCent);
Ar_CF4_95_5->GetIonisation()->SetMeanExcitationEnergy(30.0*eV);

G4Material* Ar_CF4_90_10 = 
new G4Material("Ar_CF4_90_10", density= 0.0018610*g/cm3, ncomponents=2);
Ar_CF4_90_10->AddMaterial(CF4, fractionmass=10.*perCent);
Ar_CF4_90_10->AddElement(Ar, fractionmass=90.*perCent);
Ar_CF4_90_10->GetIonisation()->SetMeanExcitationEnergy(25.0*eV);

G4Material* Ar_iso_97_3 = 
new G4Material("Ar_iso_97_3", density= 0.0016838*g/cm3, ncomponents=2);
Ar_iso_97_3->AddMaterial(iso, fractionmass=3.*perCent);
Ar_iso_97_3->AddElement(Ar, fractionmass=97.*perCent);
Ar_iso_97_3->GetIonisation()->SetMeanExcitationEnergy(25.0*eV);

G4Material* Ar_iso_95_5 = 
new G4Material("Ar_iso_95_5", density= 0.0016990*g/cm3, ncomponents=2);
Ar_iso_95_5->AddMaterial(iso, fractionmass=5.*perCent);
Ar_iso_95_5->AddElement(Ar, fractionmass=95.*perCent);
Ar_iso_95_5->GetIonisation()->SetMeanExcitationEnergy(25.0*eV);

G4Material* LH2 = 
//new G4Material("LH2", density= 0.0715*g/cm3, ncomponents=1);
new G4Material("LH2", density= 8.3e-5*g/cm3, ncomponents=1);
//new G4Material("LH2", density= 4.1906e-5*g/cm3, ncomponents=1);
LH2->AddElement(H, natoms=2);

G4Material* Myl = 
new G4Material("Mylar", density= 1.397*g/cm3, ncomponents=3);
Myl->AddElement(C, natoms=10);
Myl->AddElement(H, natoms= 8);
Myl->AddElement(O, natoms= 4);

G4Material* Epo = 
new G4Material("Epoxy", density= 1.85*g/cm3, ncomponents=3); //density of FR4 (Wikipedia)
Epo->AddElement(C, natoms=18);
Epo->AddElement(H, natoms= 20);
Epo->AddElement(O, natoms= 3);

G4Material* Air = 
new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=3);
Air->AddElement(N, fractionmass=0.781);
Air->AddElement(O, fractionmass=0.21);
Air->AddElement(Ar, fractionmass=0.009);


G4Material* Vacuum =
new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 2.73*kelvin, 3.e-18*pascal);

G4Material* beam = 
new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,
                       kStateGas, STP_Temperature, 2.e-2*bar);
beam->AddMaterial(Air, fractionmass=1.);

G4Material* Inox = 
new G4Material("Inox", density= 8.02*g/cm3, ncomponents=3);
Inox->AddElement(C, fractionmass=0.001);
Inox->AddElement(Fe, fractionmass=0.829);
Inox->AddElement(Cr, fractionmass=0.17);

G4Material* Kapton = 
new G4Material("Kapton", density= 1.42*g/cm3, ncomponents=4);
Kapton->AddElement(C, fractionmass=0.691133);
Kapton->AddElement(H, fractionmass=0.026362);
Kapton->AddElement(O, fractionmass=0.209235);
Kapton->AddElement(N, fractionmass=0.073270);
   
G4Material* Rohacell = 
new G4Material("Rohacell", density= 0.075*g/cm3, ncomponents=4);
Rohacell->AddElement(C, fractionmass=0.6014);
Rohacell->AddElement(H, fractionmass=0.0805);
Rohacell->AddElement(O, fractionmass=0.3154);
Rohacell->AddElement(N, fractionmass=0.00276);

   G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//default materials of the World
defaultMaterial  = Vacuum;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN03DetectorConstruction::Construct()
{

  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
 // G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeMINOSParameters();
   
 
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
                   WorldSizeXY,WorldSizeXY,WorldSizeZ);	//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 logicWorld,		//its logical volume				 
                                 "World",		//its name
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                               
  // Target
  //  
  solidTarget=0; logicTarget=0; physiTarget=0;
  //solidChamber=0; logicChamber=0; physiChamber=0;
  solidTPC=0; logicTPC=0; physiTPC=0;
  
     solidTarget = new G4Tubs("Target",		//its name
    		       0.,TargetRadius,TargetLength,0,360.);//size
    			     
      logicTarget = new G4LogicalVolume(solidTarget,	//its solid
      				       TargetMaterial,	//its material
      				       "Target");	//its name
    				       
      physiTarget = new G4PVPlacement(0,			//no rotation
                                     G4ThreeVector(0,0,TargetLength),	//at (0,0,0)
                                     logicTarget,	//its logical volume
                                     "Target",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
   //                                 
  // Chamber
  //
     /* solidChamber = new G4Tubs("Chamber",			//its name
                       ChamberInnerRadius,ChamberInnerRadius+ChamberThickness,ChamberLength,0,360.); //size
                       
      logicChamber = new G4LogicalVolume(solidChamber,	//its solid
                                       ChamberMaterial,	//its material
                                       "Chamber");	//its name
        physiChamber = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,ChamberLength),	//at (0,0,0)
                                     logicChamber,	//its logical volume
                                     "Chamber",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  
  //                                 
  // Inner Rohacell
  //
      solidInnerRohacell = new G4Tubs("InnerRohacell",			//its name
                       ChamberInnerRadius + ChamberThickness,ChamberInnerRadius + ChamberThickness+InnerRohacellThickness,ChamberLength,0,360.); //size
                       
      logicInnerRohacell = new G4LogicalVolume(solidInnerRohacell,	//its solid
                                       InnerRohacellMaterial,	//its material
                                       "InnerRohacell");	//its name
        physiInnerRohacell = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,ChamberLength),	//at (0,0,0)
                                     logicInnerRohacell,	//its logical volume
                                     "InnerRohacell",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  //                                 
  // Outer Rohacell
  //
      solidOuterRohacell = new G4Tubs("OuterRohacell",			//its name
                       ChamberInnerRadius + ChamberThickness + InnerRohacellThickness + KaptonThickness,ChamberInnerRadius + ChamberThickness + InnerRohacellThickness + KaptonThickness+OuterRohacellThickness,ChamberLength,0,360.); //size
                       
      logicOuterRohacell = new G4LogicalVolume(solidOuterRohacell,	//its solid
                                       OuterRohacellMaterial,	//its material
                                       "OuterRohacell");	//its name
        physiOuterRohacell = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,ChamberLength),	//at (0,0,0)
                                     logicOuterRohacell,	//its logical volume
                                     "OuterRohacell",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  //                                 
  // Kapton
  //
      solidKapton = new G4Tubs("Kapton",			//its name
                       ChamberInnerRadius + ChamberThickness+InnerRohacellThickness,ChamberInnerRadius + ChamberThickness+InnerRohacellThickness+KaptonThickness,ChamberLength,0,360.); //size
                       
      logicKapton = new G4LogicalVolume(solidKapton,	//its solid
                                       KaptonMaterial,	//its material
                                       "Kapton");	//its name
        physiKapton = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,ChamberLength),	//at (0,0,0)
                                     logicKapton,	//its logical volume
                                     "Kapton",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number*/
  //                               
  // TPC
  
      solidTPC = new G4Tubs("TPC",		//its name
                          //ChamberInnerRadius + ChamberThickness + InnerRohacellThickness + KaptonThickness + OuterRohacellThickness,TPCRadiusExt,ChamberLength,0,360.); 
                         TargetRadius,TPCRadiusExt,ChamberLength,0,360.); 
                           
      logicTPC = new G4LogicalVolume(solidTPC,    //its solid
      			                  TPCMaterial, //its material
      			                  "TPC"); //name
      			                  
      physiTPC = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,ChamberLength),	//at (0,0,0)
                                     logicTPC,	//its logical volume
                                     "TPC",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  //                               
  // window
  //
      /*solidWindow0 = new G4Tubs("WindowTube",		//its name
                          TargetRadius,TargetRadius+WindowThickness*2.,TargetLength,0,360.);
                          
      logicWindow0 = new G4LogicalVolume(solidWindow0,    //its solid
      			                  WindowMaterial, //its material
      			                  "WindowTube"); //name
      			                  
      physiWindow0 = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,TargetLength),	//at (0,0,0)
                                     logicWindow0,	//its logical volume
                                     "WindowTube",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  
      solidWindow1 = new G4Tubs("WindowEntrance",		//its name
                          0.,TargetRadius+2.*WindowThickness,WindowThickness,0,360.); 
                          
      logicWindow1 = new G4LogicalVolume(solidWindow1,    //its solid
      			                  WindowMaterial, //its material
      			                  "WindowEntrance"); //name
      			                  
      physiWindow1 = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,-1.*WindowThickness),	//at (0,0,0)
                                     logicWindow1,	//its logical volume
                                     "WindowEntrance",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number
  
      solidWindow2 = new G4Tubs("WindowOutcoming",		//its name
                          0.,TargetRadius+2.*WindowThickness,WindowThickness,0,360.); 
                          
      logicWindow2 = new G4LogicalVolume(solidWindow2,    //its solid
      			                  WindowMaterial, //its material
      			                  "WindowOutcoming"); //name
      			                  
      physiWindow2 = new G4PVPlacement(0,		//its name
                                    G4ThreeVector(0,0,2.*TargetLength+WindowThickness),	//at (0,0,0)
                                     logicWindow2,	//its logical volume
                                     "WindowOutcoming",	//its name
                                     logicWorld,	//its mother  volume
                                     false,		//no boolean operation
                                     0);		//copy number*/
     //-------------------------------------------------------------------------
   // add also My Magnetic field
   //-------------------------------------------------------------------------


   static G4bool fieldIsInitialized = false;

   if(!fieldIsInitialized)
     {


       MyMagneticField* myField = new MyMagneticField(G4ThreeVector(0.,0.,0.));
      
       G4FieldManager* fieldMgr
         = G4TransportationManager::GetTransportationManager()
         ->GetFieldManager();
       fieldMgr->SetDetectorField(myField);

  
       
         G4MagIntegratorStepper *pItsStepper;
         G4ChordFinder* pChordFinder= new G4ChordFinder(myField,
                                                        1.0e-2*mm,  // stepper size
                    pItsStepper=0);
                    fieldMgr->SetChordFinder(pChordFinder);
          

       fieldMgr->CreateChordFinder(myField);

       fieldIsInitialized = false;
     }

 
  G4Region* aRegion1 = new G4Region("TPCLog");
  logicTPC -> SetRegion(aRegion1);
  aRegion1 -> AddRootLogicalVolume(logicTPC);
  G4Region* aRegion2 = new G4Region("TargetLog");
  logicTarget -> SetRegion(aRegion2);
  aRegion2 -> AddRootLogicalVolume(logicTarget);  
  //                                        
  // Visualization attributes
  //
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);



  // Below are vis attributes that permits someone to test / play 
  // with the interactive expansion / contraction geometry system of the
  // vis/OpenInventor driver :
 /*{G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0,1,0));
  simpleBoxVisAtt->SetVisibility(true);
  logicChamber->SetVisAttributes(simpleBoxVisAtt);}*/

 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(1.,1.,0.6));
  logicTPC->SetVisAttributes(atb);}
  
 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.6,1.,1.));
  atb->SetForceSolid(true);
  logicTarget->SetVisAttributes(atb);}
  
  /*{G4VisAttributes* atb= new G4VisAttributes(G4Colour(0,0,1));
  atb->SetForceSolid(true);
  logicWindow0->SetVisAttributes(atb);}
  
  {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0,0,1));
  atb->SetForceSolid(true);
  logicWindow1->SetVisAttributes(atb);}  
  
  {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0,0,1));
  atb->SetForceSolid(true);
  logicWindow2->SetVisAttributes(atb);}  */
  
  //
  //always return the physical World
  //
  
  return physiWorld;
}

void ExN03DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) TargetMaterial = pttoMaterial;
}
void ExN03DetectorConstruction::SetChamberMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) ChamberMaterial = pttoMaterial;
}

void ExN03DetectorConstruction::SetTPCMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) TPCMaterial = pttoMaterial;
}

void ExN03DetectorConstruction::SetWindowMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) WindowMaterial = pttoMaterial;
}
void ExN03DetectorConstruction::SetInnerRohacellMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) InnerRohacellMaterial = pttoMaterial;
}
void ExN03DetectorConstruction::SetOuterRohacellMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) OuterRohacellMaterial = pttoMaterial;
}
void ExN03DetectorConstruction::SetKaptonMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) KaptonMaterial = pttoMaterial;
}

