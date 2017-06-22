using namespace std;

class OneLepton{
  public:
    OneLepton(PreAnalysis *Ana, string name = "");
    ~OneLepton();
    vector<bool> cutflow;
    vector<bool> pass();
    vector<double> w;
    void reset();
    float getmtlep();
    float getdrbl();
    TLorentzVector lep1;
    vector<vector<pair<string,int64_t>>> nCFs;
    vector<vector<pair<string,double>>> dCFs;
    void SetCR(string CR);
    enum Region{
      CRT,
      CRST,
      CRW
    };
    bool firstrun;
    uint32_t order;
    void checkCF(uint32_t pos);
    void FillCF(size_t ttype);
  public:
    OneLepton& cut_GRL();
    OneLepton& cut_LTerror();
    OneLepton& cut_trigger();
    OneLepton& cut_PV();
    OneLepton& cut_cleanjets();
    OneLepton& cut_cosmics();
    OneLepton& cut_cleanmuons();
    OneLepton& cut_nlepton();
    OneLepton& cut_lpt();
    OneLepton& cut_njet();
    OneLepton& cut_jpt();
    OneLepton& cut_MET250();
    OneLepton& cut_twobtags();
    OneLepton& cut_onebtags();
    OneLepton& cut_dp2mj2();// dphiMETJetMin2 > 0.4
    OneLepton& cut_mtlep();
    OneLepton& cut_makt120();
    OneLepton& cut_makt120_ben();
    OneLepton& cut_makt121();
    OneLepton& cut_mtbmin();
    OneLepton& cut_drbl();
    OneLepton& cut_makt080();
    OneLepton& cut_drbb();
    OneLepton& cut_MET();
  private:
    void bjetcal();
    void NameInit();
    void mtlepcal();
  private:
    bool isbjetcal;
    bool ismtlepcal;
    vector<string> CRTnames;
    vector<string> *Currnames;
    PreAnalysis *A;
    TLorentzVector mtlep;
    string Name;
    vector<string> Tnames;
    int32_t region;
    float m_mtlep;
    float dRminbl;
};

OneLepton::OneLepton(PreAnalysis *Ana, string name):region(Region::CRT),A(Ana),Name(name),cutflow(Ana->base.size(),true),w(Ana->base.size(),1.),nCFs(Ana->base.size()),dCFs(Ana->base.size()),isbjetcal(false),firstrun(true),m_mtlep(0.),ismtlepcal(false),dRminbl(-1.),order(0){
  lep1.SetPtEtaPhiM(0,0,0,0);
  mtlep.SetPtEtaPhiM(0,0,0,0);
  NameInit();
  Currnames = &CRTnames;
}

OneLepton::~OneLepton(){
  stringmanage sg;
  string treename = "Cutflow" + Tnames[region];
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

void OneLepton::SetCR(string CR){
  if (strcmp("CRT",CR.c_str())==0){
    region = Region::CRT;
    Currnames = &CRTnames;
  }else if (strcmp("CRST",CR.c_str())==0){
    region = Region::CRST;
  }else if (strcmp("CRW",CR.c_str())==0){
    region = Region::CRW;
  }else cout << "OneLepton : Error! no matched control region!" << endl;
  return;
}

void OneLepton::reset(){
  cutflow.clear();
  cutflow.resize(A->base.size(),true);
  w.clear();
  w.resize(A->base.size(),1.);
  lep1.SetPtEtaPhiM(0,0,0,0);
  mtlep.SetPtEtaPhiM(0,0,0,0);
  isbjetcal = false;
  ismtlepcal = false;
  firstrun = false;
  order = 0;
}

vector<bool> OneLepton::pass(){
  return cutflow;
}

bool compweight(objPro* aa, objPro* bb){
  return aa->yylbtagweight < bb->yylbtagweight;
}

void OneLepton::bjetcal(){
  if (A->getbjet.size()>=2){
    std::sort(A->getbjet.rbegin(),A->getbjet.rend(),compweight);
    float dR1 = sqrt(pow(A->getbjet.at(0)->yyleta-lep1.Eta(),2)+pow(A->getbjet.at(0)->yylphi-lep1.Phi(),2));
    float dR2 = sqrt(pow(A->getbjet.at(1)->yyleta-lep1.Eta(),2)+pow(A->getbjet.at(1)->yylphi-lep1.Phi(),2));
    if (dR1 > dR2){dRminbl=dR2;}
    else{dRminbl=dR1;}
  }
  isbjetcal = true;
  return;
}
float OneLepton::getdrbl(){
  if (!isbjetcal) bjetcal();
  return dRminbl;
}

void OneLepton::mtlepcal(){
  m_mtlep = sqrt(2*(lep1.Pt())*(A->tst.Pt())*(1-cos(lep1.Phi()-A->tst.Phi())));
  ismtlepcal = true;
  return;
}

float OneLepton::getmtlep(){
  if (!ismtlepcal)mtlepcal();
  return m_mtlep;
}

OneLepton& OneLepton::cut_GRL(){
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

OneLepton& OneLepton::cut_LTerror(){
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

OneLepton& OneLepton::cut_trigger(){
  checkCF(2);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (A->base[i].find("2015")!=std::string::npos){
      if ((*(A->fileinfo))["triggermatch_xe70"] == "HLT_xe70_tc_lcw"){ 
        FillCF(i);
      }else cutflow[i] = false;
    }
    if (A->base[i].find("2016_1")!=std::string::npos ){
      if ((*(A->fileinfo))["triggermatch_xe90"] == "HLT_xe90_mht_L1XE50"){
      //if ((*(A->fileinfo))["triggermatch_xe100"] == "HLT_xe100_mht_L1XE50"){
        FillCF(i);
      }else cutflow[i] = false;
    }
    if (A->base[i].find("2016_2")!=std::string::npos){
      if ((*(A->fileinfo))["triggermatch_xe100"] == "HLT_xe100_mht_L1XE50"){
        FillCF(i);
      }else cutflow[i] = false;
    }
    if (A->base[i].find("2016_3")!=std::string::npos ){
      if ((*(A->fileinfo))["triggermatch_xe110"] == "HLT_xe110_mht_L1XE50"){
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

OneLepton& OneLepton::cut_PV(){
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

OneLepton& OneLepton::cut_cleanjets(){
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

OneLepton& OneLepton::cut_cosmics(){
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

OneLepton& OneLepton::cut_cleanmuons(){
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

OneLepton& OneLepton::cut_nlepton(){
  checkCF(7);
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if (A->isMC){
      if (A->base[i].find("2015")!=std::string::npos){
        w[i] *= A->w5*A->w7;
      }
      if (A->base[i].find("2016")!=std::string::npos){
        w[i] *= A->w6*A->w8;
      }
    }
    if ((A->mu_idx[PreAnalysis::signal] == 1 && A->mu_idx[PreAnalysis::baseline] == 1 && A->ele_idx[PreAnalysis::baseline] == 0) || (A->ele_idx[PreAnalysis::signal] == 1 && A->mu_idx[PreAnalysis::baseline] == 0 && A->ele_idx[PreAnalysis::baseline] == 1)){
      if (A->mu_idx[PreAnalysis::signal] == 1){ 
        auto mu = A->getbasemu[0];
        lep1.SetPtEtaPhiM(mu->yylpt,mu->yyleta,mu->yylphi,mu->yylm);
      }
      if (A->ele_idx[PreAnalysis::signal] == 1){ 
        auto el = A->getbaseele[0];
        lep1.SetPtEtaPhiM(el->yylpt,el->yyleta,el->yylphi,el->yylm);
      }
      FillCF(i);
    }else{cutflow[i] = false;}
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_lpt(){
  checkCF(8);
  bool isreturn = false;
  for (size_t i = 0; i < cutflow.size();i++){
    isreturn |= cutflow[i];
  }
  if (!isreturn){
    order++;
    return *this;
  }
  bool mcutflow = false;
  for (auto& mu : A->getbasemu){
    lep1.SetPtEtaPhiM(mu->yylpt,mu->yyleta,mu->yylphi,mu->yylm);
    if (mu->yylpt > 20000.)
      mcutflow = true;
  }
  for (auto& el : A->getbaseele){
    lep1.SetPtEtaPhiM(el->yylpt,el->yyleta,el->yylphi,el->yylm);
    if (el->yylpt > 20000.)
      mcutflow = true;
  }
  for (size_t i = 0; i < cutflow.size();i++){
    if (!cutflow[i]) continue;
    if(mcutflow){
      FillCF(i);
    }else{cutflow[i] = false;}
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_njet(){
  checkCF(9);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->isMC){
      w[i] *= A->w3;
    }
    if (A->jet_idx[PreAnalysis::signal] + 1 >=4){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_jpt(){
  checkCF(10);
  bool con = false;
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    con = true;
  }
  if (!con){
    order++;
    return *this;
  }
  int32_t j40 = A->jet40;
  int32_t j80 = A->jet80;
  if (lep1.Pt() > 40000 ) j40++;
  if (lep1.Pt() > 80000 ) j80++;
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (j40 >= 4 && j80 >= 2){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_MET250(){
  checkCF(11);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->tst.Pt() > 250000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_twobtags(){
  checkCF(12);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->isMC){
      w[i] *= A->w4;
    }
    if (A->jet_idx[PreAnalysis::btagged] >= 2){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_onebtags(){
  if (firstrun){
    for (size_t i = 0; i < A->base.size(); i++){
      nCFs[i].push_back(make_pair("at least one b jet",0));
      dCFs[i].push_back(make_pair("at least one b jet",0.));
    }
  }
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->isMC){
      w[i] *= A->w4;
    }
    if (A->jet_idx[PreAnalysis::btagged] >= 1){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_dp2mj2(){
  checkCF(13);
  bool isreturn = false;
  for (size_t i = 0; i < cutflow.size();i++){
    isreturn |= cutflow[i];
  }
  if (!isreturn){
    order++;
    return *this;
  }
  //float dp2mfj = fabs(A->firstjet.Phi()-A->tst.Phi());
  //float dp2msj = fabs(A->secondjet.Phi()-A->tst.Phi());
  //dp2mfj = dp2mfj > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dp2mfj) : dp2mfj;
  //dp2msj = dp2msj > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dp2msj) : dp2msj;
  float dp2mfj = fabs(TVector2::Phi_mpi_pi(A->firstjet.Phi()-A->tst.Phi()));
  float dp2msj = fabs(TVector2::Phi_mpi_pi(A->secondjet.Phi()-A->tst.Phi()));
  bool mcutflow = false;
  if (dp2mfj > 0.4 && dp2msj > 0.4)mcutflow = true;
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (mcutflow){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_mtlep(){
  checkCF(14);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (!ismtlepcal)mtlepcal();
    if (m_mtlep > 30000. && m_mtlep < 120000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_makt120(){
  checkCF(15);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->j12sort.size() >= 1 && A->j12sort.at(0) > 120000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_makt120_ben(){
  checkCF(15);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if ((A->j12sort.size() >= 1 && A->j12sort.at(0) > 70000.)|| lep1.Pt()> 70000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_makt121(){
  checkCF(16);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->j12sort.size() >= 2 && A->j12sort.at(1) > 120000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_mtbmin(){
  checkCF(17);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->Mtbmin > 100000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}


OneLepton& OneLepton::cut_drbl(){
  checkCF(18);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (!isbjetcal) bjetcal();
    if (dRminbl < 1.5){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_makt080(){
  checkCF(19);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->j08sort.size()>=1 && A->j08sort.at(0) > 60000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_drbb(){
  checkCF(20);
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i]) continue;
    if (A->dRbb > 1.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_MET(){
  order++;
  return *this;
}

void OneLepton::NameInit(){
  CRTnames = {"GRL passed",
             "LAr and Tile passed",
             "trigger passed",
             "PV passed",
             "clean jet",
             "cosmic",
             "clean muon",
             "one lepton",
             "lepton pt > 20GeV",
             "4 jets at least",
             "jet pt (80, 80 , 20, 20)GeV",
             "MET > 250GeV",
             "btag >= 2",
             "dphiMetJetMin2 > 0.4",
             "30GeV < mT of the lepton and MET < 120GeV",
             "mass of first R = 1.2 jet > 120 GeV",
             "mass of second R = 1.2 jet > 120 GeV",
             "mass between b jet and MET > 100GeV",
             "min dR between lepton and b jet < 1.5",
             "mass of first R = 0.8 jet > 60 GeV",
             "dR(b,b) > 1",
             };
  Tnames = {"CRT","CRST","CRW"};
  return;
}

void OneLepton::checkCF(uint32_t pos){
  if (firstrun){
    for (size_t i = 0; i < A->base.size(); i++){
      nCFs[i].push_back(make_pair((*Currnames)[pos],0));
      dCFs[i].push_back(make_pair((*Currnames)[pos],0.));
    }
  }
  return;
}

void OneLepton::FillCF(size_t ttype){
  nCFs[ttype][order].second++;
  dCFs[ttype][order].second +=w[ttype];
  return;
}
