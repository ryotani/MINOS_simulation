#ifndef MinosModularPhysicsList_h
#define MinosModularPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class MinosModularPhysicsList: public G4VModularPhysicsList
{
public:
  MinosModularPhysicsList();
  virtual ~MinosModularPhysicsList();
 
  void SetCuts();
  };

#endif



