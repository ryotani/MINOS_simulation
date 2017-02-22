#include "TMinosTPCSimData.h"

ClassImp(TMinosTPCSimData)

TMinosTPCSimData::TMinosTPCSimData()
  : TObject(), x_mm(0), y_mm(0), t_pad(0), dt_pad(0), q_pad(0) {
// Constructor
}
TMinosTPCSimData::TMinosTPCSimData(Double_t xmm , Double_t ymm, Double_t tpad, Double_t dtpad, Double_t qpad)
  : TObject()
{
// Constructor
  x_mm = xmm;
  y_mm = ymm;
  t_pad = tpad;
  dt_pad = dtpad;
  q_pad = qpad;
}
TMinosTPCSimData::~TMinosTPCSimData()
{
// Destructor
}

void TMinosTPCSimData::Set(Double_t xmm , Double_t ymm, Double_t tpad, Double_t dtpad, Double_t qpad)
{
  x_mm = xmm;
  y_mm = ymm;
  t_pad = tpad;
  dt_pad = dtpad;
  q_pad = qpad;
}
