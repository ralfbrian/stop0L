using namespace std;
class SRA{
  public: 
    SRA(PreAnalysis *Ana);
    vector<bool> cutflow;
    vector<double> w;
    vector<bool> pass();
    void reset();
    vector<string> sranames;
    double dpmm;
    double dp2mfj;
    double dp2msj;
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
    SRA& cut_dpmj2();
    SRA& cut_trk();
    SRA& cut_dptrttst();
    SRA& cut_makt120();
    SRA& cut_makt121();
    SRA& cut_makt080();
    SRA& cut_mtbmin();
    SRA& cut_twob();
    SRA& cut_tauveto();
    SRA& cut_met400();

    private:
      PreAnalysis *A;
};

SRA::SRA(PreAnalysis *Ana):A(Ana),cutflow(2,true),w(4,1.),dpmm(10.),dp2mfj(0.),dp2msj(0.){
  sranames = {"GRL",
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
              "MET > 400"
              };
}

void SRA::reset(){
  cutflow.clear();
  cutflow.resize(2,true);
  w.clear();
  w.resize(4,1.);
  return;
}

vector<bool> SRA::pass(){
  return cutflow;
}

SRA& SRA::cut_grl(){
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

SRA& SRA::cut_dstatus(){
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
SRA& SRA::cut_trigger(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (t == 0  ){
      if ((*(A->fileinfo))["triggermatch_xe70"] == "HLT_xe70_tc_lcw"){
      }
      else{
        cutflow[t] = false;
      }
    }
    if (t == 1 ){
      if ((*(A->fileinfo))["triggermatch_xe100"] == "HLT_xe100_mht_L1XE50"){
      }
      else{
        cutflow[t] = false;
      }
    }
  }
  return *this;
}
SRA& SRA::cut_PV(){
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

SRA& SRA::cut_cleanjet(){
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

SRA& SRA::cut_cosmics(){
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

SRA& SRA::cut_cleanmuon(){
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

SRA& SRA::cut_leptonveto(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->mu_idx[PreAnalysis::baseline] == 0 && A->ele_idx[PreAnalysis::baseline] == 0){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_jetpt(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] *= A->w3;
    }
    if (A->jet40 >= 4 && A->jet80 >= 2){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_oneb(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] *= A->w4;
    }
    if (A->jet_idx[PreAnalysis::btagged] >= 1){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_met(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->tst.Pt() > 250000.){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_dpmj2(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    double dj1 = fabs(A->firstjet.Phi()-A->tst.Phi());
    double dj2 = fabs(A->secondjet.Phi()-A->tst.Phi());
    dp2mfj = dj1 > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dj1) : dj1;
    dp2msj = dj2 > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dj2) : dj2;
    if (dp2mfj > 0.4 && dp2msj > 0.4){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_trk(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->trk.Pt() > 30000.){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_dptrttst(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    double dp = fabs(A->tst.Phi()-A->trk.Phi());
    dpmm = dp > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dp) : dp;
    if (dpmm < std::atan(1.0)*4/3){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_makt120(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j12sort.size() >= 1 && A->j12sort.at(0) > 120000.){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_makt121(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j12sort.size() >= 2 && A->j12sort.at(1) > 120000.){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_makt080(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->j08sort.size() >= 1 && A->j08sort.at(0) > 60000.){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_mtbmin(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->Mtbmin > 200000.){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_twob(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->jet_idx[PreAnalysis::btagged] >= 2){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_tauveto(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->tau_idx[PreAnalysis::baseline] == 0){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}

SRA& SRA::cut_met400(){
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->tst.Pt() > 400000.){
    }
    else{
      cutflow[t] = false;
    }
  }
  return *this;
}
