#include "TMinosTrack.h"

ClassImp(TMinosTrack)

TMinosTrack::TMinosTrack()
: TObject(), x_mm(0), y_mm(0),  z_mm(0), Chargemax(0), Phi(0), Begin(0) , n_Cluster(0), n_Pads(0), RingBool(0){
// Constructor
}
TMinosTrack::TMinosTrack(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t phi, Double_t begin, Double_t ncluster, Double_t npads, Double_t ringbool)
  : TObject()
{
  // Constructor
  x_mm = xmm;
  y_mm= ymm;
  z_mm = zmm;
  Chargemax = chargemax;
  Phi=phi;
  Begin = begin;
  n_Cluster=ncluster;
  n_Pads = npads;
  RingBool = ringbool;
}
TMinosTrack::~TMinosTrack()
{
// Destructor
}

void TMinosTrack::Set(Double_t xmm, Double_t ymm, Double_t zmm, Double_t chargemax, Double_t phi, Double_t begin, Double_t ncluster, Double_t npads, Double_t ringbool)
{
  x_mm = xmm;
  y_mm= ymm;
  z_mm = zmm;
  Chargemax = chargemax;
  Phi=phi;
  Begin = begin;
  n_Cluster=ncluster;
  n_Pads = npads;
  RingBool = ringbool;
}
