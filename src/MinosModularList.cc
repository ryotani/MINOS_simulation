#include "MinosModularPhysicsList.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4HadronPhysicsINCLXX.hh"
//#include "HadronPhysicsQGSP_INCLXX.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
//#include "G4HadronPhysicsQGSP_BERT_CHIPS.hh"


MinosModularPhysicsList::MinosModularPhysicsList():G4VModularPhysicsList()
{
  defaultCutValue = 1000*mm;
  SetVerboseLevel(1);
  this->RegisterPhysics(new G4EmStandardPhysics);
  this->RegisterPhysics(new G4DecayPhysics);
  this->RegisterPhysics(new G4HadronPhysicsINCLXX);
  this->RegisterPhysics(new G4IonINCLXXPhysics);
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MinosModularPhysicsList::~MinosModularPhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MinosModularPhysicsList::SetCuts()
{
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "alpha");//Ajout AO
  SetCutValue(defaultCutValue, "proton");//Ajout AO

  //G4double lowlimit=25.*eV; //mix
  G4double lowlimit=1.*keV; //H2
  G4ProductionCutsTable::GetProductionCutsTable() ->SetEnergyRange(lowlimit, 100.*GeV);

  G4String regionName1 = "TPCLog";
  G4Region* region1 = G4RegionStore::GetInstance()->GetRegion(regionName1);
  G4String regionName2 = "TargetLog";
  G4Region* region2 = G4RegionStore::GetInstance()->GetRegion(regionName2);
  G4ProductionCuts* cuts = new G4ProductionCuts ;
  G4double regionCut = 0.05*mm;
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("gamma"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e-"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("e+"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("alpha"));
  cuts -> SetProductionCut(regionCut,G4ProductionCuts::GetIndex("proton"));
  region1 -> SetProductionCuts(cuts);
  region2 -> SetProductionCuts(cuts);

  if (verboseLevel>0) DumpCutValuesTable();

}
