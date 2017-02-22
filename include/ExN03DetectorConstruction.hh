#ifndef ExN03DetectorConstruction_h
#define ExN03DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "/Users/acorsi/codes/MINOS_simulation/lib/ExN03Setup.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    ExN03DetectorConstruction();
   ~ExN03DetectorConstruction();

  public:
     
     void SetTPCMaterial (G4String);     

     void SetChamberMaterial (G4String);     
     
     void SetTargetMaterial (G4String);     
           
     void SetInnerRohacellMaterial (G4String);     

     void SetKaptonMaterial (G4String);     

     void SetWindowMaterial (G4String);     
    
     void SetOuterRohacellMaterial (G4String);     

     G4VPhysicalVolume* Construct();

  public:
  
     
     G4double    GetTargetLength()      {return TargetLength*2.;};
     G4Material* GetTargetMaterial()    {return TargetMaterial;};
     
     G4double    GetTargetRadius()      {return TargetRadius;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetTarget()     {return physiTarget;};
     const G4VPhysicalVolume* GetChamber()    {return physiChamber;};
     const G4VPhysicalVolume* GetTPC()        {return physiTPC;};
     const G4VPhysicalVolume* GetWindow0()     {return physiWindow0;};
     const G4VPhysicalVolume* GetWindow1()     {return physiWindow1;};                 
     const G4VPhysicalVolume* GetWindow2()     {return physiWindow2;}; 
     const G4VPhysicalVolume* GetInnerRohacell()        {return physiInnerRohacell;};
     const G4VPhysicalVolume* GetOuterRohacell()        {return physiOuterRohacell;};
     const G4VPhysicalVolume* GetKapton()        {return physiKapton;};
     
  private:
     
     G4Material*        TargetMaterial;
     G4double           TargetRadius;
     G4double           TargetLength;
     
     G4Material*        WindowMaterial;
     G4double           WindowThickness;
     
     G4Material*        ChamberMaterial;
     G4double           ChamberInnerRadius;
     G4double           ChamberLength;
     G4double           ChamberThickness;
    
     G4Material*        InnerRohacellMaterial;
     G4double           InnerRohacellThickness;

     G4Material*        OuterRohacellMaterial;
     G4double           OuterRohacellThickness;
     
     G4Material*        KaptonMaterial;
     G4double           KaptonThickness;
     
     G4Material*        TPCMaterial;
     G4double           TPCRadiusExt;
         
     G4Material*        defaultMaterial;
     G4double           WorldSizeXY;
     G4double           WorldSizeZ;
            
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Tubs*             solidTarget;   
     G4LogicalVolume*   logicTarget;   
     G4VPhysicalVolume* physiTarget;   
     
     G4Tubs*             solidChamber;  
     G4LogicalVolume*   logicChamber;  
     G4VPhysicalVolume* physiChamber;  
         
     G4Tubs*             solidTPC; 
     G4LogicalVolume*   logicTPC; 
     G4VPhysicalVolume* physiTPC; 
     
     G4Tubs*             solidWindow0; 
     G4LogicalVolume*   logicWindow0; 
     G4VPhysicalVolume* physiWindow0; 
     
     G4Tubs*             solidWindow1; 
     G4LogicalVolume*   logicWindow1; 
     G4VPhysicalVolume* physiWindow1; 

     G4Tubs*             solidWindow2; 
     G4LogicalVolume*   logicWindow2; 
     G4VPhysicalVolume* physiWindow2; 
    
     G4Tubs*             solidInnerRohacell;   
     G4LogicalVolume*   logicInnerRohacell;   
     G4VPhysicalVolume* physiInnerRohacell;   
     
     G4Tubs*             solidOuterRohacell;   
     G4LogicalVolume*   logicOuterRohacell;   
     G4VPhysicalVolume* physiOuterRohacell;   
     
     G4Tubs*             solidKapton;   
     G4LogicalVolume*   logicKapton;   
     G4VPhysicalVolume* physiKapton;   
     
      
  private:
    
     void DefineMaterials();
     void ComputeMINOSParameters();
     G4VPhysicalVolume* ConstructMINOS();     
     ExN03Setup *setup;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void ExN03DetectorConstruction::ComputeMINOSParameters()
{
     WorldSizeZ = 1.2*2.*ChamberLength; WorldSizeXY = 4*TPCRadiusExt;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

