using namespace std;

class OneLepton{
  public:
    OneLepton(PreAnalysis *Ana, string name = "");
    ~OneLepton();
    vector<bool> cutflow;
    vector<bool> pass();
    vector<double> w;
    void reset();
    double getmtlep();
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
    OneLepton& cut_dp2mj2();// dphiMETJetMin2 > 0.4
    OneLepton& cut_mtlep();
    OneLepton& cut_makt120();
    OneLepton& cut_mtbmin();
    OneLepton& cut_drbl();
  private:
    void bjetcal();
    void NameInit();
    void mtlepcal();
  private:
    bool isbjetcal;
    vector<string> CRTnames;
    vector<string> *Currnames;
    PreAnalysis *A;
    TLorentzVector lep1;
    TLorentzVector mtlep;
    string Name;
    string Tnames;
    int32_t region;
};

OneLepton::OneLepton(PreAnalysis *Ana, string name = ""):region(Region::CRT),A(Ana),Name(name),cutflow(Ana->base.size(),true),w(Ana->base.size(),1.),nCFs(Ana->base.size()),dCFs(Ana->base.size()),isbjetcal(false),firstrun(true){
  lep1.SetPtEtaPhiM(0,0,0,0);
  mtlep.SetPtEtaPhiM(0,0,0,0);
  NameInit();
  Currnames = &CRTnames;
}

OneLepton::~OneLepton(){
  stringmanage sg;
  string treename = "Cutflow" + Tnames(region);
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
  firstrun = false;
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
    float dR1 = sqrt(pow(getbjet.at(0)->yyleta-lep1.Eta(),2)+pow(getbjet.at(0)->yylphi-lep1.Phi(),2));
    float dR2 = sqrt(pow(getbjet.at(1)->yyleta-lep1.Eta(),2)+pow(getbjet.at(1)->yylphi-lep1.Phi(),2));
    if (dR1 > dR2){dRminbl=dR2;}
    else{dRminbl=dR1;}
  }
  isbjetcal = true;
  return;
}
void OneLepton::mtlepcal(){

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
    if (A->base[i].find("2016_1")!=std::string::npos){
      if ((*(A->fileinfo))["triggermatch_xe100"] == "HLT_xe100_mht_L1XE50"){
        FillCF(i);
      }else cutflow[i] = false;
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
    if ((mu_idx[OneLepton::baseline] == 1 && el_idx[OneLepton::baseline] == 0) || (mu_idx[OneLepton::baseline] == 0 && el_idx[OneLepton::baseline] == 1)){
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
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] *= A->w3;
    }
    if (jet_idx[OneLepton::baseline]>=4){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_jpt(){
  checkCF(9);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (jet40 >= 4 && jet80 >= 2){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_MET250(){
  checkCF(10);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->tst.Pt() > 250000.){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_twobtags(){
  checkCF(11);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (A->isMC){
      w[t] *= A->w4;
    }
    if (jet_idx[OneLepton::btagged] >= 2){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_dp2mj2(){
  checkCF(12);
  bool isreturn = false;
  for (size_t i = 0; i < cutflow.size();i++){
    isreturn |= cutflow[i];
  }
  if (!isreturn){
    order++;
    return *this;
  }
  float dp2mfj = fabs(A->firstjet.Phi()-A->tst.Phi());
  float dp2msj = fabs(A->secondjet.Phi()-A->tst.Phi());
  dp2mfj = dp2mfj > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dp2mfj) : dp2mfj;
  dp2msj = dp2msj > (std::atan(1.0)*4) ? (std::atan(1.0)*8-dp2msj) : dp2msj;
  bool mcutflow = false;
  if (dp2mfj > 0.4 && dp2msj > 0.4)mcutflow = true;
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
    if (mcutflow){
      FillCF(i);
    }else cutflow[i] = false;
  }
  order++;
  return *this;
}

OneLepton& OneLepton::cut_mtlep(){
  checkCF(13);
  for (size_t t = 0; t < cutflow.size(); t++){
    if (!cutflow[t]) continue;
  mtlep = lep1 + tst;
  float m_mtlep = sqrt(2*(lep1.Pt())*(tst.Pt())*(1-cos(lep1.Phi()-tst.Phi())));
  if (m_mtlep > 30000. && m_mtlep < 120000.){
    if(isMC){
      nOneLepton16[14].second++;
      nOneLepton15[14].second++;
      dOneLepton16[14].second+=weight16;
      dOneLepton15[14].second+=weight15;
    }
    else{
      if (isdata15){
        nOneLepton15[14].second++;
        dOneLepton15[14].second+=weight15;
      }
      if (isdata16){
        nOneLepton16[14].second++;
        dOneLepton16[14].second+=weight16;
      }
    }
  }else{cutflow = false;}
  return *this;
}

OneLepton& OneLepton::cut_makt120(){
  if (!cutflow) return *this;
  if (!isjet12cal) jet12cal();
  if (j12sort.size() >= 1 && j12sort.at(0) > 70000.){
    if(isMC){
      nOneLepton16[15].second++;
      nOneLepton15[15].second++;
      dOneLepton16[15].second+=weight16;
      dOneLepton15[15].second+=weight15;
    }
    else{
      if (isdata15){
        nOneLepton15[15].second++;
        dOneLepton15[15].second+=weight15;
      }
      if (isdata16){
        nOneLepton16[15].second++;
        dOneLepton16[15].second+=weight16;
      }
    }
  }else{cutflow = false;}
  return *this;
}

OneLepton& OneLepton::cut_mtbmin(){
  if (!cutflow) return *this;
  if (!isjetcal) jetcal();
  if (Mtbmin > 100000.){
    if(isMC){
      nOneLepton16[16].second++;
      nOneLepton15[16].second++;
      dOneLepton16[16].second+=weight16;
      dOneLepton15[16].second+=weight15;
    }
    else{
      if (isdata15){
        nOneLepton15[16].second++;
        dOneLepton15[16].second+=weight15;
      }
      if (isdata16){
        nOneLepton16[16].second++;
        dOneLepton16[16].second+=weight16;
      }
    }
  }else{cutflow = false;}
  return *this;
}


OneLepton& OneLepton::cut_drbl(){
  if (!cutflow) return *this;
  if (!isbjetcal) bjetcal();
  if (dRminbl < 1.5){
    if(isMC){
      nOneLepton16[17].second++;
      nOneLepton15[17].second++;
      dOneLepton16[17].second+=weight16;
      dOneLepton15[17].second+=weight15;
    }
    else{
      if (isdata15){
        nOneLepton15[17].second++;
        dOneLepton15[17].second+=weight15;
      }
      if (isdata16){
        nOneLepton16[17].second++;
        dOneLepton16[17].second+=weight16;
      }
    }
  }else{cutflow = false;}
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
             "dphiMetJetMin2 <= 0.4",
             "30GeV < transverse mass of the lepton and MET < 120GeV",
             "mass of first R = 1.2 jet > 70 GeV",
             "mass between b jet and MET > 100GeV",
             "min dR between lepton and b jet < 1.5"};
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
