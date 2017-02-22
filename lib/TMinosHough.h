#ifndef TMINOSHOUGH_H
#define TMINOSHOUGH_H

//#include "TClonesArray.h"
#include "TObject.h"

class TMinosHough : public TObject {
 public:
  TMinosHough();
  TMinosHough(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t ncluster);
  void Set(Double_t xmm, Double_t ymm, Double_t zmm , Double_t chargemax, Double_t ncluster);
  void SetZ(Double_t zmm);
  virtual ~TMinosHough();
 public:
  Double_t x_mm;
  Double_t y_mm;
  Double_t z_mm;
  Double_t Chargemax;
  Double_t n_Cluster;

  ClassDef(TMinosHough,1);
};
#endif // end of #ifndef TMINOSHOUGH_H
