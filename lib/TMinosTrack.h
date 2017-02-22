#ifndef TMINOSTRACK_H
#define TMINOSTRACK_H

//#include "TClonesArray.h"
#include "TObject.h"

class TMinosTrack : public TObject {
 public:
  TMinosTrack();
  TMinosTrack(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t phi, Double_t begin, Double_t ncluster, Double_t npads, Double_t ringbool);
  void Set(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t phi, Double_t begin, Double_t ncluster, Double_t npads, Double_t ringbool);
  virtual ~TMinosTrack();
 public:
  Double_t x_mm;
  Double_t y_mm;
  Double_t z_mm;
  Double_t Chargemax;
  Double_t Phi;
  Double_t Begin;
  Double_t n_Cluster;
  Double_t n_Pads;
  Double_t RingBool;

  ClassDef(TMinosTrack,1);
};
#endif // end of #ifndef TMINOSTRACK_H
