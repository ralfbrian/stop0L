#ifndef _SRA_H_
#define _SRA_H_
#ifndef _PREANALYSIS_H_
#include "PreAnalysis.h"
#endif
using namespace std;
class SRA{
  public: 
    SRA(PreAnalysis *Ana,string name = "TT");
    ~SRA();
    vector<bool> cutflow;
    vector<double> w;
    vector<bool> pass();
    void reset();
    vector<string> SRAnames;
    double dpmm;
    double dp2m1j;
    double dp2m2j;
    double dp2m3j;
    bool findkey(vector<string>& a, string key);
    // for developing new method
    void checkCF(uint32_t pos);
    void FillCF(size_t ttype);
    uint32_t order;
    bool firstrun;
    vector<vector<pair<string,int64_t>>> nCFs;
    vector<vector<pair<string,double>>> dCFs;
  public:
    SRA& cut_grl();
    SRA& cut_dstatus();
    SRA& cut_trigger();
    SRA& cut_PV();
    SRA& cut_cleanjet();
    SRA& cut_cosmics();
    SRA& cut_cleanmuon();
    SRA& cut_leptonveto();
    SRA& cut_jetpt();
    SRA& cut_oneb();
    SRA& cut_met();
    SRA& cut_dpmj3();
    SRA& cut_trk();
    SRA& cut_dptrttst();
    SRA& cut_makt120();
    SRA& cut_makt121();
    SRA& cut_makt080();
    SRA& cut_mtbmin();
    SRA& cut_twob();
    SRA& cut_tauveto();
    SRA& cut_met400();
    SRA& cut_drbb();
    SRA& cut_MT2();
  public:
    double makt120;
    double makt121;
    double makt120u;
    double makt121u;
    double makt080;
    double mtbmin;
    double met400;

  private:
    PreAnalysis *A;
    string Name;
    vector<uint32_t> orders;
};

SRA::SRA(PreAnalysis *Ana, string name):A(Ana),Name(name),cutflow(Ana->base.size(),true),w(Ana->base.size(),1.),nCFs(Ana->base.size()),dCFs(Ana->base.size()),dpmm(10.),dp2m1j(0.),dp2m2j(0.),dp2m3j(0.),order(0),firstrun(true){
  SRAnames = {"GRL",
              "LAr and Tile",
              "trigger",
              "PV",
              "clean jet",
              "cosmic",
              "clean muon",
              "lepton veto",
              "jet pt",
              "least one b jet",
              "MET > 250",
              "dphiMetJetMin2 > 0.4",
              "track met > 30 GeV",
              "delta phi between track met and tst met",
              "least two b jets",
              "mass of first R = 1.2 jet > 120 GeV",
              "mass of second R = 1.2 jet > 120 GeV",
              "mass of R = 0.8 jet > 60 GeV",
              "mass between b jet and MET",
              "tau veto",
              "MET > 400",
              "DRBB > 1",
              "MT2 > 400"
              };
  makt120 = 120000.;
  makt121 = 120000.;
  makt080 = 60000.;
  mtbmin = 200000.;
  met400 = 400000.;
  makt120u = -1.;
  makt121u = -1.;
  for(map<string,double>::iterator itr = Ana->mapd.begin();itr!=Ana->mapd.end();itr++){
      if ((*itr).first.find("SRA_makt120")==0) makt120 = (*itr).second;
      if ((*itr).first.find("SRA_makt121")==0) makt121 = (*itr).second;
      if ((*itr).first.find("SRA_makt120u")==0) makt120u = (*itr).second;
      if ((*itr).first.find("SRA_makt121u")==0) makt121u = (*itr).second;
      if ((*itr).first.find("SRA_makt080")==0) makt080 = (*itr).second;
      if ((*itr).first.find("SRA_mtbmin")==0) mtbmin = (*itr).second;
      if ((*itr).first.find("SRA_met400")==0) met400 = (*itr).second;
  }
}

SRA::~SRA(){
  stringmanage sg;
  string treename = "CutflowSRA";
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

bool SRA::findkey(vector<string>& a, string key){
  for (auto& aa : a){
    if (strcmp(aa.c_str(),key.c_str())==0)
      return true;
  }
  return false;
}

void SRA::reset(){
  cutflow.clear();
  cutflow.resize(A->base.size(),true);
  w.clear();
  w.resize(A->base.size(),1.);
  order = 0;
  firstrun = false;
  return;
}

vector<bool> SRA::pass(){
  return cutflow;
}

SRA& SRA::cut_grl(){
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

SRA& SRA::cut_dstatus(){
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
SRA& SRA::cut_trigger(){
  checkCF(2);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->base[t].find("2015")!=std::string::npos ){
      if ((*(A->fileinfo))["triggermatch_xe70"] == "HLT_xe70_tc_lcw"){
        FillCF(t);
      }
      else{
        cutflow[t] = false;
      }
    }
    if (A->base[t].find("2016_1")!=std::string::npos ){
      if ((*(A->fileinfo))["triggermatch_xe90"] == "HLT_xe90_mht_L1XE50"){
      //if ((*(A->fileinfo))["triggermatch_xe100"] == "HLT_xe100_mht_L1XE50"){
        FillCF(t);
      }
      else{
        cutflow[t] = false;
      }
    }
    if (A->base[t].find("2016_2")!=std::string::npos ){
      if ((*(A->fileinfo))["triggermatch_xe100"] == "HLT_xe100_mht_L1XE50"){
        FillCF(t);
      }
      else{
        cutflow[t] = false;
      }
    }
    if (A->base[t].find("2016_3")!=std::string::npos ){
      if ((*(A->fileinfo))["triggermatch_xe110"] == "HLT_xe110_mht_L1XE50"){
        FillCF(t);
      }
      else{
        cutflow[t] = false;
      }
    }
  }
  order++;
  return *this;
}
SRA& SRA::cut_PV(){
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

SRA& SRA::cut_cleanjet(){
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

SRA& SRA::cut_cosmics(){
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

SRA& SRA::cut_cleanmuon(){
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

SRA& SRA::cut_leptonveto(){
  checkCF(7);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->mu_idx[PreAnalysis::baseline] == 0 && A->ele_idx[PreAnalysis::baseline] == 0){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_jetpt(){
  checkCF(8);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] *= A->w3;
    }
    if (A->jet40 >= 4 && A->jet80 >= 2){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_oneb(){
  checkCF(9);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] *= A->w4;
    }
    if (A->jet_idx[PreAnalysis::btagged] >= 1){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_met(){
  checkCF(10);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->tst.Pt() > 250000.){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_dpmj3(){
  checkCF(11);
  bool isreturn = false;
  for (size_t i = 0; i < cutflow.size();i++){
    isreturn |= cutflow[i];
  }
  if (!isreturn){
    order++;
    return *this;
  }
  dp2m1j = fabs(TVector2::Phi_mpi_pi(A->firstjet.Phi()-A->tst.Phi()));
  dp2m2j = fabs(TVector2::Phi_mpi_pi(A->secondjet.Phi()-A->tst.Phi()));
  dp2m3j = fabs(TVector2::Phi_mpi_pi(A->getbasejet[2]->yylphi-A->tst.Phi()));
  bool mcutflow = false;
  if (dp2m1j > 0.4 && dp2m2j > 0.4 && dp2m3j > 0.4)mcutflow = true;
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    //double dj1 = fabs(A->firstjet.Phi()-A->tst.Phi());
    //double dj2 = fabs(A->secondjet.Phi()-A->tst.Phi());
    //dp2mfj = dj1 > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dj1) : dj1;
    //dp2msj = dj2 > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dj2) : dj2;
    if (mcutflow){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_trk(){
  checkCF(12);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->trk.Pt() > 30000.){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_dptrttst(){
  checkCF(13);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    double dp = fabs(A->tst.Phi()-A->trk.Phi());
    dpmm = dp > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dp) : dp;
    if (dpmm < std::atan(1.0)*4/3){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_makt120(){
  checkCF(15);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j12sort.size() >= 1 && A->j12sort.at(0) > makt120 && (makt120u < 0. || A->j12sort.at(0) < makt120u)){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_makt121(){
  checkCF(16);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j12sort.size() >= 2 && A->j12sort.at(1) > makt121 && (makt121u < 0. ||A->j12sort.at(1) < makt121u)){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_makt080(){
  checkCF(17);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j08sort.size() >= 1 && A->j08sort.at(0) > makt080){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_mtbmin(){
  checkCF(18);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->Mtbmin > mtbmin){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_twob(){
  checkCF(14);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->jet_idx[PreAnalysis::btagged] >= 2){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_tauveto(){
  checkCF(19);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->tau_idx[PreAnalysis::baseline] == 0){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_met400(){
  checkCF(20);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->tst.Pt() > met400){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_drbb(){
  checkCF(21);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->dRbb > 1.){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

SRA& SRA::cut_MT2(){
  checkCF(22);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    A->CalMT2();
    if (A->MT2> 400000.){
      FillCF(t);
    }
    else{
      cutflow[t] = false;
    }
  }
  order++;
  return *this;
}

void SRA::checkCF(uint32_t pos){
  if (firstrun){
    for (size_t i = 0; i < A->base.size(); i++){
      nCFs[i].push_back(make_pair(SRAnames[pos],0));
      dCFs[i].push_back(make_pair(SRAnames[pos],0.));
    }
    orders.push_back(pos);
  }
  return;
}

void SRA::FillCF(size_t ttype){
  nCFs[ttype][order].second++;
  dCFs[ttype][order].second +=w[ttype];
  return;
}

#endif
