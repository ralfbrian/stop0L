#ifdef __CINT__
#pragma link C++ class vector<objPro>+;
#pragma link C++ class map<string,metPro>+;
#pragma link C++ class map<string,string>+;
#pragma link C++ class map<string,double>+;
#endif
#ifndef _PREANALYSIS_H_
#define _PREANALYSIS_H_
using namespace std;
class PreAnalysis : public HistInit{
  public:
  PreAnalysis();

  enum sel {
    signal,
    baseline,
    bad,
    cosmic,
    goodpt,
    btagged,
    trgmatch,
    passOR
  };

  void SetTree(TTree *tt);

  void SetName(string fname);

  void init();

  void FillLz(TLorentzVector& vec, objPro *obj);

  bool chi2Top();

  void CalMT2();

  public:
  // event info
  string filename;
  int32_t runnumber;
  int64_t eventnumber;
  bool isMC;
  bool isSignal;
  bool d16;
  bool d15;
  // weight
  double w0;//mceventweight
  double xsec;
  double w1;//pileup
  double w2;// jet SF
  double w3;// jvt SF
  double w4;// b tag SF
  double w5;
  double w6;
  double w7;
  double w8;
  double w9;
  double w10;

  TLorentzVector tst;
  TLorentzVector trk;
  TLorentzVector firstjet;
  TLorentzVector secondjet;
  TLorentzVector topCand0;
  TLorentzVector topCand1;
  TLorentzVector WCand0;
  TLorentzVector WCand1;
  double Chi2min;
  double MT2;
  int32_t jet20;
  int32_t jet40;
  int32_t jet80;
  std::vector<float> jetptsort; // GeV
  std::vector<double> j12sort; // MeV
  std::vector<double> j08sort; // MeV
  float Mtbmin;//MeV
  float fabsdphi;
  float dRbb;
  std::vector<int> jet_idx;
  std::vector<int> ele_idx;
  std::vector<int> mu_idx;
  std::vector<int> tau_idx;
  std::vector<objPro*> getbasejet;
  std::vector<objPro*> getbaseele;
  std::vector<objPro*> getbasemu;
  std::vector<objPro*> getbjet;

  // pre selection
  bool passGRL;
  bool PV;
  bool Dstatus;
  bool calmt2;

  private:
  
  void jetcal();
  
  void ecal();
  
  void mcal();

  void taucal();
  
  void j12cal();

  void j08cal();

  void RClustering(std::vector<objPro*>& obj,std::vector<objPro>& rcjet,float dr, float ptmin);

  void Bsort(vector<objPro*>& v_bjets);
  
  public:
    std::vector<objPro> *econtainer;
    std::vector<objPro> *mcontainer;
    std::vector<objPro> *pcontainer;
    std::vector<objPro> *tcontainer;
    std::vector<objPro> *jcontainer;
    std::vector<objPro> *j12container;
    std::vector<objPro> *j08container;
    std::map<std::string,double> *m_weight;
    std::map<std::string,metPro> *metcontainer;
    std::map<std::string,std::string> *fileinfo;
};

#endif
