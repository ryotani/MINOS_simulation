#include "TMinosHough.h"

ClassImp(TMinosHough)

TMinosHough::TMinosHough()
: TObject(), x_mm(0), y_mm(0), z_mm(0), Chargemax(0), n_Cluster(0){
// Constructor
}
TMinosHough::TMinosHough(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t ncluster)
  : TObject()
{
  // Constructor
  x_mm = xmm;
  y_mm= ymm;
  z_mm = zmm;
  Chargemax = chargemax;
  n_Cluster=ncluster;
}
TMinosHough::~TMinosHough()
{
// Destructor
}

void TMinosHough::Set(Double_t xmm, Double_t ymm, Double_t zmm, Double_t chargemax, Double_t ncluster)
{
  x_mm = xmm;
  y_mm= ymm;
  z_mm = zmm;
  Chargemax = chargemax;
  n_Cluster=ncluster;
}

void TMinosHough::SetZ(Double_t zmm)
{
  z_mm = zmm;
}
