#ifndef TMINOSTPCSIMDATA_H
#define TMINOSTPCSIMDATA_H

//#include "TClonesArray.h"
#include "TObject.h"

class TMinosTPCSimData : public TObject {
 public:
  TMinosTPCSimData();
  TMinosTPCSimData(Double_t xmm , Double_t ymm, Double_t tpad, Double_t dtpad, Double_t qpad);
  void Set(Double_t xmm , Double_t ymm, Double_t tpad, Double_t dtpad, Double_t qpad);
  virtual ~TMinosTPCSimData();
 public:
  Double_t x_mm;
  Double_t y_mm; 
  Double_t t_pad;
  Double_t dt_pad; //difference between fitted time and time barycenter of the primary electrons
  Double_t q_pad;

  ClassDef(TMinosTPCSimData,1);
};
#endif // end of #ifndef TMINOSTPCSIMDATA_H
