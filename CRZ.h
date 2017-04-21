
using namespace std;
class CRZ 
{
  public:
    CRZ(PreAnalysis *A);
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
  // cut flow 
    CRZ& cut_grl();
    CRZ& cut_dstatus();
    CRZ& cut_trigger();
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
  private:
    bool ismuon;
    bool isele;
    bool runmetp;
    bool runzmass;
    int pele;
    int mele;
    int pmuon;
    int mmuon;
    TLorentzVector e1;
    TLorentzVector e2;
    TLorentzVector m1;
    TLorentzVector m2;
    TLorentzVector Zmass;
    TLorentzVector METp;
    PreAnalysis *A;
};

CRZ::CRZ(PreAnalysis *Ana) : A(Ana),pele(0),mele(0),pmuon(0),mmuon(0),cutflow(Ana->base.size(),true),w(Ana->base.size(),1.),ismuon(false),isele(false),runmetp(false),runzmass(false){
  m1.SetPtEtaPhiM(0,0,0,0);
  m2.SetPtEtaPhiM(0,0,0,0);
  e1.SetPtEtaPhiM(0,0,0,0);
  e2.SetPtEtaPhiM(0,0,0,0);
  Zmass.SetPtEtaPhiM(0,0,0,0);
  METp.SetPtEtaPhiM(0,0,0,0);
}

CRZ::~CRZ(){
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
  pele = 0;
  mele = 0;
  pmuon = 0;
  mmuon = 0;
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
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] = A->w0*A->xsec*A->w1;
    }
    if (A->passGRL){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

CRZ& CRZ::cut_dstatus(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->Dstatus){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

CRZ& CRZ::cut_trigger(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
    if (A->base[i].find("2015")!=std::string::npos){
      if((*(A->fileinfo))["triggermatch_el_2015"] == "true" || (*(A->fileinfo))["triggermatch_mu_2015"] == "true"){
      }
      else{
        cutflow[i] = false;
      }
    }
    if (A->base[i].find("2016_1")!=std::string::npos){
      if((*(A->fileinfo))["triggermatch_el_2016"] == "true" || (*(A->fileinfo))["triggermatch_mu_2016"] == "true"){
      }
      else{
        cutflow[i] = false;
      }
    }
  }
  return *this;
}

CRZ& CRZ::cut_PV(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->PV){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

CRZ& CRZ::cut_cleanjet(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->jet_idx[PreAnalysis::bad] == 0){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

CRZ& CRZ::cut_cosmics(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->mu_idx[PreAnalysis::cosmic] == 0){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

CRZ& CRZ::cut_cleanmuon(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->mu_idx[PreAnalysis::bad] == 0){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

CRZ& CRZ::cut_nlepton(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
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
    }else{cutflow[i] = false;}
  }
  return *this;
}

CRZ& CRZ::cut_lpt(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
    // muon selection
    if (ismuon){
      int mcheckpt = 0;
      for (auto& muon : A->getbasemu){
        if ((muon->signal == 1 && muon->yylpt > 20000.) || (muon->trgmatch_2015 == 1 && muon->yylpt > 25000.) || (muon->trgmatch_2016 == 1 && muon->yylpt > 25000.)){
          mcheckpt++;
        }
      }
      if (mcheckpt!=2) cutflow[i] = false;
    }
    // electron selection
    if (isele){
      int echeckpt = 0;
      for (auto& electron : A->getbaseele){
        if ((electron->signal == 1 && electron->yylpt > 20000.) || (electron->trgmatch_2015 == 1 && electron->yylpt > 30000.) || (electron->trgmatch_2016 == 1 && electron->yylpt > 30000.)){
          echeckpt++;
        }
      }
      if (echeckpt!=2) cutflow[i] = false;
    }
  }
  return *this;
}

CRZ& CRZ::cut_zmass(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
    if (!runzmass){
      if (ismuon){
        Zmass = m1 + m2;
      }
      if (isele){
        Zmass = e1 + e2;
      }
    }
    runzmass = false;
    // lepton invariant mass
    if(!(Zmass.M() > 86000. && Zmass.M() < 96000.)) cutflow[i] = false;
  }
  return *this;
}

CRZ& CRZ::cut_njet(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
    // jet multiplicity
    if (!(A->getbasejet.size()>= 4)) cutflow[i] = false;
  }
  return *this;
}

CRZ& CRZ::cut_jpt(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
    // jet pt requirement
    if(!(A->jet20 >= 4 && A->jet40 >= 2)) cutflow[i] = false;
  }
  return *this;
}

CRZ& CRZ::cut_met(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
    // MET
    if (!(A->tst.Pt() < 50000.)) cutflow[i] = false;
  }
  return *this;
}

CRZ& CRZ::cut_metp(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
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
  }
  return *this;
}

CRZ& CRZ::cut_btag(){
  for (int i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) return *this;
    // bjet
    if (!(A->jet_idx[PreAnalysis::btagged]>= 2 )) cutflow[i] = false;
  }
  return *this;
}

