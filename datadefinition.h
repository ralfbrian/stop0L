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
  bool signal = false;
  bool baseline = false;
  bool bad = false;
  bool cosmic = false;
  bool btagged = false;
  bool passOR = false;
  bool trgmatch_2015 = false;
  bool trgmatch_2016 = false;
  // this is needed to distribute the algorithm to the workers
};

class metPro : public TObject{
  public :
  double yylmpx = -1e4;
  double yylmpy = -1e4;
  double yylmet = -1e4;
  double yylphi = -1e4;
  double yylsumet = -1e4;
  // this is needed to distribute the algorithm to the workers
};

