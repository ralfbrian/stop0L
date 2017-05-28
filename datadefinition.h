#ifndef _DATADEFINITION_H_
#define _DATADEFINITION_H_
#include "TObject.h"

class objPro : public TObject{
  public :
  float yylcharge = -1e4;
  double yylpt = -1e4;
  double yyleta = -1e4;
  double yylphi = -1e4;
  double yylm = -1e4;
  double yyle = -1e4;
  double yylrapidity = -1e4;
  float yylbtagweight = -1e4;
  // value 1 or 0 is initialzied. -100 means non-initialized
  uint8_t signal = 0;
  uint8_t baseline = 0;
  uint8_t bad = 0;
  uint8_t cosmic = 0;
  uint8_t btagged = 0;
  uint8_t passOR = 0;
  uint8_t trgmatch_2015 = 0;
  uint8_t trgmatch_2016 = 0;
  uint32_t yylnTrk = 0;
};

class metPro : public TObject{
  public :
  double yylmpx = -1e4;
  double yylmpy = -1e4;
  double yylmet = -1e4;
  double yylphi = -1e4;
  double yylsumet = -1e4;
};

#endif
