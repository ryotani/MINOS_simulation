//---------- Written by C. Santamaria / CEA Saclay !----------
#include "TMinosClust.h"
#include <iostream>

ClassImp(TMinosClust)

TMinosClust::TMinosClust()
: TObject(), x_mm(0), y_mm(0), z_mm(0), Chargemax(0), n_Cluster(0), n_Pads(0), RingBool(0){
// Constructor
}
TMinosClust::TMinosClust(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t ncluster, Double_t npads, Double_t ringbool)
  : TObject()
{
  // Constructor
  x_mm = xmm;
  y_mm= ymm;
  //t_ns = tns;
  z_mm = zmm;
  Chargemax = chargemax;
  n_Cluster = ncluster;
  n_Pads = npads;
  RingBool = ringbool;
}
TMinosClust::~TMinosClust()
{
// Destructor
}

void TMinosClust::Set(Double_t xmm, Double_t ymm, Double_t zmm, Double_t chargemax, Double_t ncluster, Double_t npads, Double_t ringbool)
{
  x_mm = xmm;
  y_mm= ymm;
  //t_ns = tns;
  z_mm = zmm;
  Chargemax = chargemax;
  n_Cluster = ncluster;
  n_Pads = npads;
  RingBool = ringbool;
}

void TMinosClust::SetZ(Double_t zmm)
{
  z_mm = zmm;
}
