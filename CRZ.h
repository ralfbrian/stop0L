
using namespace std;
class CRZ 
{
  public:
    CRZ(PreAnalysis *A,string name = "");
    ~CRZ();
    vector<bool> pass();
    vector<bool> cutflow;
    vector<double> w;
    void reset();
    double getmetp();
    double getZmass();
    double getZpt();
    double get1stmpt();
    double get2ndmpt();
    double get1stept();
    double get2ndept();
    bool findkey(vector<string>& a, string key);
    // for developing new method
    void checkCF(uint32_t pos);
    void FillCF(size_t ttype);
    bool firstrun;
    uint32_t order;
    vector<vector<pair<string,int64_t>>> nCFs;
    vector<vector<pair<string,double>>> dCFs;
  public:
  // cut flow 
    CRZ& cut_grl();
    CRZ& cut_dstatus();
    CRZ& cut_trigger();
    CRZ& cut_trigger_benchmark();
    CRZ& cut_PV();
    CRZ& cut_cleanjet();
    CRZ& cut_cosmics();
    CRZ& cut_cleanmuon();
    CRZ& cut_nlepton();
    CRZ& cut_lpt();
    CRZ& cut_zmass();
    CRZ& cut_njet();
    CRZ& cut_jpt();
    CRZ& cut_met();
    CRZ& cut_metp();
    CRZ& cut_btag();
    CRZ& cut_onebtag();
    CRZ& cut_makt120();
    CRZ& cut_makt121();
  private:
    void NameInit();
  private:
    bool ismuon;
    bool isele;
    bool runmetp;
    bool runzmass;
    TLorentzVector e1;
    TLorentzVector e2;
    TLorentzVector m1;
    TLorentzVector m2;
    TLorentzVector Zmass;
    TLorentzVector METp;
    PreAnalysis *A;
    string Name;
    vector<uint32_t> orders;
    vector<string> CRZnames;
};

CRZ::CRZ(PreAnalysis *Ana, string name) : A(Ana),Name(name),cutflow(Ana->base.size(),true),w(Ana->base.size(),1.),nCFs(Ana->base.size()),dCFs(Ana->base.size()),ismuon(false),isele(false),runmetp(false),runzmass(false),firstrun(true),order(0){
  m1.SetPtEtaPhiM(0,0,0,0);
  m2.SetPtEtaPhiM(0,0,0,0);
  e1.SetPtEtaPhiM(0,0,0,0);
  e2.SetPtEtaPhiM(0,0,0,0);
  Zmass.SetPtEtaPhiM(0,0,0,0);
  METp.SetPtEtaPhiM(0,0,0,0);
  NameInit();
}

CRZ::~CRZ(){
  stringmanage sg;
  string treename = "CutflowCRZ";
  if (Name.size()!=0){
    treename += "_" + Name;
  }
  string checkname = treename;
  int32_t counter = 1;
  while (A->fout->FindObjectAny(checkname.c_str())){
    checkname = treename + "_" + sg.ntos(counter);
    counter++;
  }
  if (counter!=1){
    treename += "_" + sg.ntos(counter);
  }
  A->fout->cd();
  TTree *tout = new TTree(treename.c_str(),treename.c_str());
  tout->SetDirectory(A->fout);
  for (size_t i = 0; i < A->base.size(); i++){
    tout->Branch(("n" + A->base[i]).c_str(),&nCFs[i]);
    tout->Branch(("d" + A->base[i]).c_str(),&dCFs[i]);
  }
  tout->Fill();
}

bool CRZ::findkey(vector<string>& a, string key){
  for (auto& aa : a){
    if (strcmp(aa.c_str(),key.c_str())==0)
      return true;
  }
  return false;
}

void CRZ::reset(){
  cutflow.clear();
  cutflow.resize(A->base.size(),true);
  w.clear();
  w.resize(A->base.size(),1.);
  ismuon = false;
  isele= false;
  runmetp = false;
  runzmass = false;
  firstrun = false;
  order = 0;
  m1.SetPtEtaPhiM(0,0,0,0);
  m2.SetPtEtaPhiM(0,0,0,0);
  e1.SetPtEtaPhiM(0,0,0,0);
  e2.SetPtEtaPhiM(0,0,0,0);
  Zmass.SetPtEtaPhiM(0,0,0,0);
  METp.SetPtEtaPhiM(0,0,0,0);
}


vector<bool> CRZ::pass(){
  return cutflow;
}


double CRZ::getmetp(){
  if (!runmetp){
    if (ismuon){
      METp = m1 + m2 + A->tst;
    }
    if (isele){
      METp = e1 + e2 + A->tst;
    }
  }
  runmetp = true;
  return METp.Pt();
}

double CRZ::getZmass(){
  if (!runzmass){
    if (ismuon){
      Zmass = m1 + m2;
    }
    if (isele){
      Zmass = e1 + e2;
    }
  }
  runzmass = true;
  return Zmass.M();
}

double CRZ::getZpt(){
  if (!runzmass){
    if (ismuon){
      Zmass = m1 + m2;
    }
    if (isele){
      Zmass = e1 + e2;
    }
  }
  runzmass = true;
  return Zmass.Pt();
}

double CRZ::get1stmpt(){
  return m1.Pt();
}

double CRZ::get2ndmpt(){
  return m2.Pt();
}

double CRZ::get1stept(){
  return e1.Pt();
}

double CRZ::get2ndept(){
  return e2.Pt();
}

CRZ& CRZ::cut_grl(){
  checkCF(0);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] = A->w0*A->xsec*A->w1;
    }
    if (A->passGRL){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_dstatus(){
  checkCF(1);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->Dstatus){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_trigger(){
  checkCF(2);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (A->base[i].find("2015")!=std::string::npos){
      if((*(A->fileinfo))["triggermatch_el_2015"] == "true" || (*(A->fileinfo))["triggermatch_mu_2015"] == "true"){
        FillCF(i);
      }
      else{
        cutflow[i] = false;
      }
    }
    if (A->base[i].find("2016")!=std::string::npos){
      if((*(A->fileinfo))["triggermatch_el_2016"] == "true" || (*(A->fileinfo))["triggermatch_mu_2016"] == "true"){
        FillCF(i);
      }
      else{
        cutflow[i] = false;
      }
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_trigger_benchmark(){
  checkCF(2);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if((*(A->fileinfo))["triggermatch_el_2015"] == "true" || (*(A->fileinfo))["triggermatch_mu_2015"] == "true" || (*(A->fileinfo))["triggermatch_el_2016"] == "true" || (*(A->fileinfo))["triggermatch_mu_2016"] == "true"){
      FillCF(i);
    }
    else{
      cutflow[i] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_PV(){
  checkCF(3);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->PV){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_cleanjet(){
  checkCF(4);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->jet_idx[PreAnalysis::bad] == 0){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_cosmics(){
  checkCF(5);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->mu_idx[PreAnalysis::cosmic] == 0){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_cleanmuon(){
  checkCF(6);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->mu_idx[PreAnalysis::bad] == 0){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_nlepton(){
  checkCF(7);
  int32_t pele(0);
  int32_t mele(0);
  int32_t pmuon(0);
  int32_t mmuon(0);
  bool isreturn = false;
  for (size_t i = 0; i < cutflow.size();i++){
    isreturn |= cutflow[i];
    if (A->isMC){
      if (A->base[i].find("2015")!=std::string::npos){
        w[i] *= A->w5*A->w7;
      }
      if (A->base[i].find("2016")!=std::string::npos){
        w[i] *= A->w6*A->w8;
      }
    }
  }
  if (!isreturn){
    order++;
    return *this;
  }
  // muon selection
  for (auto& muon : A->getbasemu){
      if (muon->passOR == 1){
        if (muon->yylcharge > 0.9 && muon->yylcharge < 1.1){
          m1.SetPtEtaPhiM(muon->yylpt,muon->yyleta,muon->yylphi,muon->yylm);
          pmuon++;
        }
        if (muon->yylcharge > -1.1 && muon->yylcharge < -0.9){
          m2.SetPtEtaPhiM(muon->yylpt,muon->yyleta,muon->yylphi,muon->yylm);
          mmuon++;
        }
      }
  }
  // electron selection
  for (auto& electron : A->getbaseele){
      if (electron->passOR == 1 ){
        if (electron->yylcharge > 0.9 && electron->yylcharge < 1.1){
          e1.SetPtEtaPhiM(electron->yylpt,electron->yyleta,electron->yylphi,electron->yylm);
          pele++;
        }
        if (electron->yylcharge > -1.1 && electron->yylcharge < -0.9){
          e2.SetPtEtaPhiM(electron->yylpt,electron->yyleta,electron->yylphi,electron->yylm);
          mele++;
        }
      }
  }
  // 2 oppite charge lepton
  if((pmuon == 1 && mmuon == 1 && pele == 0 && mele == 0) || (pmuon == 0 && mmuon == 0 && pele == 1 && mele == 1)){
      if(pmuon == 1 && mmuon ==1){
        ismuon = true;
      }
      if(pele == 1 && mele ==1){
        isele = true;
      }
  }
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if(ismuon || isele){
      FillCF(i);
    }else{cutflow[i] = false;}
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_lpt(){
  checkCF(8);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    // muon selection
    if (ismuon){
      int32_t mcheckpt = 0;
      for (auto& muon : A->getbasemu){
        if (A->base[i].find("2015")!=std::string::npos){
          if ((muon->signal == 1 && muon->yylpt > 20000.) || (muon->trgmatch_2015 == 1 && muon->yylpt > 25000.)){
            mcheckpt++;
          }
        }
        if (A->base[i].find("2016")!=std::string::npos){
          if ((muon->signal == 1 && muon->yylpt > 20000.) || (muon->trgmatch_2016 == 1 && muon->yylpt > 25000.)){
            mcheckpt++;
          }
        }
      }
      if (mcheckpt!=2) cutflow[i] = false;
    }
    // electron selection
    if (isele){
      int32_t echeckpt = 0;
      for (auto& electron : A->getbaseele){
        if (A->base[i].find("2015")!=std::string::npos){
          if ((electron->signal == 1 && electron->yylpt > 20000.) || (electron->trgmatch_2015 == 1 && electron->yylpt > 30000.) ){
            echeckpt++;
          }
        }
        if (A->base[i].find("2016")!=std::string::npos){
          if ((electron->signal == 1 && electron->yylpt > 20000.) || (electron->trgmatch_2016 == 1 && electron->yylpt > 30000.)){
            echeckpt++;
          }
        }
      }
      if (echeckpt!=2) cutflow[i] = false;
    }
    if (cutflow[i])
      FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_zmass(){
  checkCF(9);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (!runzmass){
      if (ismuon){
        Zmass = m1 + m2;
      }
      if (isele){
        Zmass = e1 + e2;
      }
    }
    runzmass = true;
    // lepton invariant mass
    if(!(Zmass.M() > 86000. && Zmass.M() < 96000.)) cutflow[i] = false;
    else FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_njet(){
  checkCF(10);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (A->isMC){
      w[i] *= A->w3;
    }
    // jet multiplicity
    if (!(A->getbasejet.size()>= 4)) cutflow[i] = false;
    else FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_jpt(){
  checkCF(11);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    // jet pt requirement
    if(!(A->jet40 >= 4 && A->jet80 >= 2)) cutflow[i] = false;
    else FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_met(){
  checkCF(12);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    // MET
    if (!(A->tst.Pt() < 50000.)) cutflow[i] = false;
    else FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_metp(){
  checkCF(13);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (!runmetp){
      if (ismuon){
        METp = m1 + m2 + A->tst;
      }
      if (isele){
        METp = e1 + e2 + A->tst;
      }
    }
    runmetp = true;
    // MET with invisible leptons
    if (!(METp.Pt() > 70000. )) cutflow[i] = false;
    else FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_btag(){
  checkCF(14);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (A->isMC){
      w[i] *= A->w4;
    }
    // bjet
    if (!(A->jet_idx[PreAnalysis::btagged]>= 2 )) cutflow[i] = false;
    else FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_onebtag(){
  if (firstrun){
    for (size_t i = 0; i < A->base.size(); i++){
      nCFs[i].push_back(make_pair("one b jet",0));
      dCFs[i].push_back(make_pair("one b jet",0.));
    }
  }
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (A->isMC){
      w[i] *= A->w4;
    }
    // bjet
    if (!(A->jet_idx[PreAnalysis::btagged]== 1 )) cutflow[i] = false;
    else FillCF(i);
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_makt120(){
  checkCF(15);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j12sort.size() >= 1 && A->j12sort.at(0) > 120000.){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

CRZ& CRZ::cut_makt121(){
  checkCF(16);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j12sort.size() >= 2 && A->j12sort.at(1) > 60000.){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

void CRZ::NameInit(){
  CRZnames = {"GRL",
              "LAr and Tile",
              "trigger",
              "PV",
              "clean jet",
              "cosmic",
              "clean muon",
              "2 leptons",
              "lepton pt",
              "86GeV < Z mass < 96GeV",
              "at least 4 jets",
              "jet pt(80,80,40,40)",
              "MET < 50GeV",
              "METp > 70GeV",
              "at least 2 b jets",
              "mass of first R = 1.2 jet > 120 GeV",
              "mass of second R = 1.2 jet > 60 GeV"};
  return;
}

void CRZ::checkCF(uint32_t pos){
  if (firstrun){
    for (size_t i = 0; i < A->base.size(); i++){
      nCFs[i].push_back(make_pair(CRZnames[pos],0));
      dCFs[i].push_back(make_pair(CRZnames[pos],0.));
    }
    orders.push_back(pos);
  }
  return;
}

void CRZ::FillCF(size_t ttype){
  nCFs[ttype][order].second++;
  dCFs[ttype][order].second +=w[ttype];
  return;
}
