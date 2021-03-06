using namespace std;
class HistInit{
  public:
    HistInit();
    struct HistInfo{
      int type;
      string name;
      string title;
      int32_t nbins;
      double x1;
      double x2;
    };
    TFile *fout;
    //cutflow
    TH1 *h_SRATT[9];
    TH1 *h_SRATT_NEW[9];
    TH1 *h_CRZ[9];
    TH1 *h_CRZ_NEW[9];
    TH1 *h_CRT[9];
    TH1 *h_CRT_NEW[9];
    // CRZ
    TH1 *h_CRZ_njetr[9];
    TH1 *h_CRZ_eptr[9];
    TH1 *h_CRZ_mptr[9];
    TH1 *h_CRZ_jet1r[9];
    TH1 *h_CRZ_jet2r[9];
    TH1 *h_CRZ_jet3r[9];
    TH1 *h_CRZ_jet4r[9];
    TH1 *h_CRZ_jet_121r[9];
    TH1 *h_CRZ_jet_122r[9];
    TH1 *h_CRZ_jet_08r[9];
    TH1 *h_CRZ_Zr[9];
    TH1 *h_CRZ_METr[9];
    TH1 *h_CRZ_METpr[9];
    TH1 *h_CRZ_nbjetr[9];
    TH1 *h_CRZ_nbjetr_NEW[9];
    TH1 *h_CRZ_njet[9];
    TH1 *h_CRZ_ept[9];
    TH1 *h_CRZ_mpt[9];
    TH1 *h_CRZ_jet1[9];
    TH1 *h_CRZ_jet2[9];
    TH1 *h_CRZ_jet3[9];
    TH1 *h_CRZ_jet4[9];
    TH1 *h_CRZ_jet_121[9];
    TH1 *h_CRZ_jet_122[9];
    TH1 *h_CRZ_jet_08[9];
    TH1 *h_CRZ_Z[9];
    TH1 *h_CRZ_MET[9];
    TH1 *h_CRZ_METp[9];
    TH1 *h_CRZ_nbjet[9];
    TH1 *h_CRZ_nbjet_NEW[9];
    TH1 *h_CRZ_F_Z[9];
    TH1 *h_CRZ_F_MET[9];
    TH1 *h_CRZ_F_METp[9];
    TH1 *h_CRZ_F_nbjet[9];
    TH1 *h_CRZ_Zpt[9];
    TH1 *h_CRZ_F_Zpt[9]; 
    TH1 *h_CRZ_ne[9];
    TH1 *h_CRZ_nm[9];
    TH1 *h_CRZ_makt120[9];
    TH1 *h_CRZ_makt121[9];
    TH1 *h_CRZ_makt120_1b[9];
    TH1 *h_CRZ_makt121_1b[9];
    // CRT
    TH1 *h_CRT_08e[9];
    TH1 *h_CRT_08m[9];
    TH1 *h_CRT_09[9]; 
    TH1 *h_CRT_101[9];
    TH1 *h_CRT_102[9];
    TH1 *h_CRT_103[9];
    TH1 *h_CRT_104[9];
    TH1 *h_CRT_11[9]; 
    TH1 *h_CRT_12[9]; 
    TH1 *h_CRT_14[9];
    TH1 *h_CRT_14_p[9];
    TH1 *h_CRT_15[9];
    TH1 *h_CRT_16[9]; 
    TH1 *h_CRT_17[9]; 
    TH1 *h_CRT_18[9]; 
    TH1 *h_CRT_MET_F[9];
    TH1 *h_CRT_Makt08[9];
    // SRA
    TH1 *h_cut11[9]; 
    TH1 *h_cut12[9]; 
    TH1 *h_cut13[9]; 
    TH1 *h_cut151[9];
    TH1 *h_cut152[9];
    TH1 *h_cut153[9];
    TH1 *h_cut154[9];
    TH1 *h_cut16[9]; 
    TH1 *h_cut17[9]; 
    TH1 *h_cut18[9];
    TH1 *h_cut19[9];
    TH1 *h_cut20[9];
    TH1 *h_cut21[9];
    TH1 *h_cut22[9];
    TH1 *h_cut23[9];
    TH1 *h_cut24[9];
    TH1 *h_cut25[9];
    TH1 *h_cut27[9]; 
    TH1 *h_cut28_p[9]; 
    TH1 *h_cut28[9]; 
    TH1 *h_cut29_p[9]; 
    TH1 *h_cut29[9]; 
    TH1 *h_cut22_F[9];
    TH1 *h_cut23_F[9];
    TH1 *h_cut22_new[9];
    TH1 *h_cut23_new[9];
    TH1 *h_cut24_new[9];
    TH1 *h_cut25_new[9];
    TH1 *h_cut27_new[9]; 
    TH1 *h_cut22_new_F[9];
    TH1 *h_cut23_new_F[9];
    TH1 *h_mass[9];
    TH1 *h_mass_e[9];
    TH1 *h_mass_mu[9];
    TH1 *h_mass_jet[9];
    TH1 *h_check_bjet[9];

    TH1* NEWTH1F(string name, string title, int32_t bins, double fbin, double lbin); // type = 1
    TH1* NEWTH1D(string name, string title, int32_t bins, double fbin, double lbin);// type = 2
    TH1* NEWTH1I(string name, string title, int32_t bins, double fbin, double lbin);// type = 3
    void Init(TH1 *(*h)[9],struct HistInfo HI);
    TH1** GH(TH1 *(*h)[9]);
    vector<string> base;
    map<string,double> mapd;
    vector<pair<TH1*(*)[9],struct HistInfo>> histinfo;
    struct HistInfo SetInfo(string name, string title, int32_t bins, double fbin, double lbin, int32_t type);
    vector<TH1*> hall;
    void SetTFile(TFile *file);
    void Readconfig();

};

HistInit::HistInit():fout(0){
  Readconfig();
  base = {"_2015","_2016_1","_2016_2","_2016_3"};
  auto itr = mapd.find("BASE");
  if (itr!=mapd.end() && itr->second > 0.6) base.resize((uint32_t)(itr->second+0.5));
  //base = {"_2015","_2016_1"};
  //base = {"_2015"};
  Init(&h_CRZ, SetInfo("CRZ","Cut Flow of CRZ",17,0.5,17.5,1));
  Init(&h_CRZ_NEW, SetInfo("CRZ_NEW","Cut Flow of CRZ",17,0.5,17.5,1));
  Init(&h_CRZ_njetr, SetInfo("CRZnjetr", "number of jet reference;Njet", 16,-0.5,15.5,1));
  Init(&h_CRZ_eptr, SetInfo("CRZeptr","electron pt reference;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_mptr, SetInfo("CRZmptr","muon pt reference;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_jet1r, SetInfo("CRZjet1r","Jet 1 pt reference;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRZ_jet2r, SetInfo("CRZjet2r","Jet 2 pt reference;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRZ_jet3r, SetInfo("CRZjet3r","Jet 3 pt reference;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRZ_jet4r, SetInfo("CRZjet4r","Jet 4 pt reference;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRZ_jet_121r, SetInfo("CRZjet_121r","R1.2 jet 1 reference;m[GeV];Events/25GeV",18,0,450,1));
  Init(&h_CRZ_jet_122r, SetInfo("CRZjet_122r","R1.2 jet 2 reference;m[GeV];Events/25GeV",10,0,250,1));
  Init(&h_CRZ_jet_08r, SetInfo("CRZjet_08r","R0.8 jet 1 reference;m[GeV];Events/30GeV",15,0,450,1));
  Init(&h_CRZ_Zr, SetInfo("CRZZr", "Zmass reference;m [GeV];Events/0.2GeV", 50,86,96,1));
  Init(&h_CRZ_METr, SetInfo("CRZMETr", "MET reference;ETmiss [GeV];Events/20GeV", 20,0,400,1));
  Init(&h_CRZ_METpr, SetInfo("CRZMETpr", "MET(p) reference;ETmiss [GeV];Events/20GeV", 40,0,800,1));
  Init(&h_CRZ_nbjetr, SetInfo("CRZnbjetr", "number of bjet reference;Nbjet", 11,-0.5,10.5,1));
  Init(&h_CRZ_nbjetr_NEW, SetInfo("CRZnbjetr_NEW", "number of bjet reference;Nbjet", 11,-0.5,10.5,1));
  Init(&h_CRZ_njet, SetInfo("CRZnjet", "number of jet;Njet", 16,-0.5,15.5,1));
  Init(&h_CRZ_ept, SetInfo("CRZept","electron pt;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_mpt, SetInfo("CRZmpt","muon pt;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_jet1, SetInfo("CRZjet1","Jet 1 pt;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRZ_jet2, SetInfo("CRZjet2","Jet 2 pt;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRZ_jet3, SetInfo("CRZjet3","Jet 3 pt;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRZ_jet4, SetInfo("CRZjet4","Jet 4 pt;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRZ_jet_121, SetInfo("CRZjet_121","R1.2 jet 1;m[GeV];Events/25GeV",18,0,450,1));
  Init(&h_CRZ_jet_122, SetInfo("CRZjet_122","R1.2 jet 2;m[GeV];Events/25GeV",10,0,250,1));
  Init(&h_CRZ_jet_08, SetInfo("CRZjet_08","R0.8 jet 1;m[GeV];Events/30GeV",15,0,450,1));
  Init(&h_CRZ_Z, SetInfo("CRZZ", "Zmass ;m [GeV];Events/0.2GeV", 50,86,96,1));
  Init(&h_CRZ_MET, SetInfo("CRZMET", "MET;ETmiss [GeV];Events/20GeV", 20,0,400,1));
  Init(&h_CRZ_METp, SetInfo("CRZMETp", "MET(p) reference;ETmiss [GeV];Events/20GeV", 40,0,800,1));
  Init(&h_CRZ_nbjet, SetInfo("CRZnbjet", "number of bjet;Nbjet", 11,-0.5,10.5,1));
  Init(&h_CRZ_nbjet_NEW, SetInfo("CRZnbjet_NEW", "number of bjet;Nbjet", 11,-0.5,10.5,1));
  Init(&h_CRZ_F_Z, SetInfo("CRZZ_F", "Zmass ;m [GeV];Events/0.2GeV", 50,86,96,1));
  Init(&h_CRZ_F_MET, SetInfo("CRZMET_F", "MET;ETmiss [GeV];Events/20GeV", 20,0,400,1));
  Init(&h_CRZ_F_METp, SetInfo("CRZMETp_F", "MET(p) reference;ETmiss [GeV];Events/20GeV", 40,0,800,1));
  Init(&h_CRZ_F_nbjet, SetInfo("CRZnbjet_F", "number of bjet;Nbjet", 11,-0.5,10.5,1));
  Init(&h_CRZ_Zpt, SetInfo("CRZZpt","Z pt;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_F_Zpt, SetInfo("CRZZpt_F","Z pt;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_ne, SetInfo("CRZne","electron size;N. electron;Events/1",11,-0.5,10.5,1));
  Init(&h_CRZ_nm, SetInfo("CRZnm","muon size;N. muon;Events/1",11,-0.5,10.5,1));
  Init(&h_CRZ_makt120, SetInfo("CRZmakt120", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_CRZ_makt121, SetInfo("CRZmakt121", "Makt12,1 > 120GeV;Makt12,1 [GeV]", 52,0,780,1));
  Init(&h_CRZ_makt120_1b, SetInfo("CRZmakt120_1b", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_CRZ_makt121_1b, SetInfo("CRZmakt121_1b", "Makt12,1 > 120GeV;Makt12,1 [GeV]", 52,0,780,1));
  // CRT
  Init(&h_CRT, SetInfo("CRT","Cut Flow of CRT",21,0.5,21.5,1));
  Init(&h_CRT_NEW, SetInfo("CRT_NEW","Cut Flow of CRT",18,0.5,18.5,1));
  Init(&h_CRT_08e, SetInfo("CRT_08e","electron pt;pt[GeV];Events/20GeV",15,0,300,1));
  Init(&h_CRT_08m, SetInfo("CRT_08m","muon pt;pt[GeV];Events/20GeV",15,0,300,1));
  Init(&h_CRT_09, SetInfo("CRT_09", "number of jet;Njet", 15,0.5,15.5,1));
  Init(&h_CRT_101, SetInfo("CRT_101","Jet 1 pt;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRT_102, SetInfo("CRT_102","Jet 2 pt;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRT_103, SetInfo("CRT_103","Jet 3 pt;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRT_104, SetInfo("CRT_104","Jet 4 pt;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRT_11, SetInfo("CRT_11", "MET;ETmiss [GeV];Events/20GeV", 50,0,1000,1));
  Init(&h_CRT_12, SetInfo("CRT_12", "number of bjet;Nbjet", 10,0.5,10.5,1));
  Init(&h_CRT_14_p, SetInfo("CRT_14_p", "Mt lep ;[GeV]", 26,0,390,1));
  Init(&h_CRT_14, SetInfo("CRT_14", "Mt lep ;[GeV]", 26,0,390,1));
  Init(&h_CRT_15, SetInfo("CRT_15", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_CRT_16, SetInfo("CRT_16", "Makt12,1 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_CRT_17, SetInfo("CRT_17", "MtBMin > 100GeV;MtBMin [GeV]", 20,0,1000,1));
  Init(&h_CRT_18, SetInfo("CRT_18", "dR(b,l) < 1.5;MtBMin [GeV]", 100,0,2,1));
  Init(&h_CRT_Makt08, SetInfo("CRT_Makt08", "Makt08 > 60GeV;Makt08,0 [GeV]", 20,0,600,1));
  Init(&h_CRT_MET_F, SetInfo("CRT_MET_F", "MET;ETmiss [GeV];Events/20GeV", 50,0,1000,1));
  // SRA
  Init(&h_cut11, SetInfo("cut11", "clean jet;Njet", 15,0.5,15.5,1));
  Init(&h_cut12, SetInfo("cut12", "cosmic;Nmuon", 15,0.5,15.5,1));
  Init(&h_cut13, SetInfo("cut13", "clean muon;Nmuon", 15,0.5,15.5,1));
  Init(&h_cut151, SetInfo("cut151", "signal jet (jet1);pT [GeV]", 10,0,400,1));
  Init(&h_cut152, SetInfo("cut152", "signal jet (jet2);pT [GeV]", 10,0,400,1));
  Init(&h_cut153, SetInfo("cut153", "signal jet (jet3);pT [GeV]", 10,0,400,1));
  Init(&h_cut154, SetInfo("cut154", "signal jet (jet4);pT [GeV]", 10,0,400,1));
  Init(&h_cut16, SetInfo("cut16", "one b jet;Nbtag", 6,0.5,6.5,1));
  Init(&h_cut17, SetInfo("cut17", "MET > 250GeV ;ETmiss [GeV]", 20,0,1000,1));
  Init(&h_cut18, SetInfo("cut18", "d phi between MET and 2 jets > 0.4 ;d phi [rad]", 40,0,4,1));
  Init(&h_cut19, SetInfo("cut19", "MET track  > 30GeV;ETmiss,track [GeV]", 40,0,400,1));
  Init(&h_cut20, SetInfo("cut20", "d phi between MET track and MET < pi/3;d phi [rad]", 32,0,3.2,1));
  Init(&h_cut21, SetInfo("cut21", "two b tag;Nbtag", 6,0.5,6.5,1));
  Init(&h_cut22, SetInfo("cut22", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_cut23, SetInfo("cut23", "Makt12,1 > 120GeV;Makt12,1 [GeV]", 52,0,780,1));
  Init(&h_cut24, SetInfo("cut24", "Makt08 > 60GeV;Makt08,0 [GeV]", 20,0,600,1));
  Init(&h_cut25, SetInfo("cut25", "MtBMin > 200GeV;MtBMin [GeV]", 50,0,1000,1));
  Init(&h_cut27, SetInfo("cut27", "MET > 400GeV;ETmiss [GeV]", 20,0,1000,1)); 
  Init(&h_cut28_p, SetInfo("cut28_p", "dRbb;dR", 50,0,10,1)); 
  Init(&h_cut28, SetInfo("cut28", "dRbb;dR", 50,0,10,1)); 
  Init(&h_cut29_p, SetInfo("cut29_p", "MT2; [GeV]", 40,0,1000,1)); 
  Init(&h_cut29, SetInfo("cut29", "MT2;[GeV]", 40,0,1000,1)); 
  Init(&h_cut22_F, SetInfo("cut22_F", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_cut23_F, SetInfo("cut23_F", "Makt12,1 > 120GeV;Makt12,1 [GeV]", 52,0,780,1));
  Init(&h_SRATT, SetInfo("SRATT_O","Cut Flow of SRA-TT",23,0,23,1));
  Init(&h_cut22_new, SetInfo("cut22_new", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_cut23_new, SetInfo("cut23_new", "Makt12,1 > 120GeV;Makt12,1 [GeV]", 52,0,780,1));
  Init(&h_cut24_new, SetInfo("cut24_new", "Makt08 > 60GeV;Makt08,0 [GeV]", 20,0,600,1));
  Init(&h_cut25_new, SetInfo("cut25_new", "MtBMin > 200GeV;MtBMin [GeV]", 50,0,1000,1));
  Init(&h_cut27_new, SetInfo("cut27_new", "MET > 400GeV;ETmiss [GeV]", 20,0,1000,1)); 
  Init(&h_cut22_new_F, SetInfo("cut22_new_F", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 52,0,780,1));
  Init(&h_cut23_new_F, SetInfo("cut23_new_F", "Makt12,1 > 120GeV;Makt12,1 [GeV]", 52,0,780,1));
  Init(&h_SRATT_NEW, SetInfo("SRATT_NEW","Cut Flow of New SRA-TT",20,0,20,1));
  // interest
  Init(&h_mass, SetInfo("MASS","MASS SPECTRUM;mass [GeV]",10000,0,1000,1));
  Init(&h_mass_e, SetInfo("MASS_E","MASS SPECTRUM;e mass [GeV]",1000,0,10,1));
  Init(&h_mass_mu, SetInfo("MASS_MU","MASS SPECTRUM;mu mass [GeV]",1000,0,10,1));
  Init(&h_mass_jet, SetInfo("MASS_JET","MASS SPECTRUM;jet mass [GeV]",10000,0,1000,1));
  Init(&h_check_bjet, SetInfo("Check_Bjet", "b tag jets;Nbtag", 11,-0.5,10.5,1));

}

TH1* HistInit::NEWTH1F(string name, string title, int32_t bins, double fbin, double lbin){
  fout->cd();
  TH1 *hh = new TH1F(name.c_str(),title.c_str(),bins,fbin,lbin);
  hh->SetDirectory(fout);
  hh->Sumw2();
  hall.push_back(hh);
  return hh;
}

TH1* HistInit::NEWTH1D(string name, string title, int32_t bins, double fbin, double lbin){
  fout->cd();
  TH1 *hh = new TH1D(name.c_str(),title.c_str(),bins,fbin,lbin);
  hh->SetDirectory(fout);
  hh->Sumw2();
  hall.push_back(hh);
  return hh;
}

TH1* HistInit::NEWTH1I(string name, string title, int32_t bins, double fbin, double lbin){
  fout->cd();
  TH1 *hh = new TH1I(name.c_str(),title.c_str(),bins,fbin,lbin);
  hh->SetDirectory(fout);
  hh->Sumw2();
  hall.push_back(hh);
  return hh;
}

void HistInit::Init(TH1 *(*h)[9], HistInit::HistInfo HI){
  for(size_t i = 0; i < base.size(); i++){
    (*h)[i] = 0;
  }
  histinfo.push_back(make_pair(h,HI));
  return;
}

TH1** HistInit::GH(TH1 *(*h)[9]){
  if (!(*h)[0]){
    for(size_t i = 0; i < histinfo.size(); i++){
      if (histinfo[i].first == h){
        string title = histinfo[i].second.title;
        int32_t nbins = histinfo[i].second.nbins;
        double x1 = histinfo[i].second.x1;
        double x2 = histinfo[i].second.x2;
        for(size_t j = 0; j< base.size();j++){
          string name = histinfo[i].second.name + base[j];
          switch(histinfo[i].second.type){
            case 1:
              (*h)[j] = NEWTH1F(name,title,nbins,x1,x2);
              break;
            case 2:
              (*h)[j] = NEWTH1D(name,title,nbins,x1,x2);
              break;
            case 3:
              (*h)[j] = NEWTH1D(name,title,nbins,x1,x2);
              break;
          }
        }
      }
    }
  }
  return *h;
}

HistInit::HistInfo HistInit::SetInfo(string name, string title, int32_t bins, double fbin, double lbin, int32_t type){
  struct HistInfo hs;
  hs.name = name;
  hs.title = title;
  hs.nbins = bins;
  hs.x1 = fbin;
  hs.x2 = lbin;
  hs.type = type;
  return hs;
}

void HistInit::SetTFile(TFile *file){
  fout = file;
  return;
}

void HistInit::Readconfig(){
  string currDir = gSystem->pwd();
  fstream fp;
  fp.open("YYL.conf",ios::in);
  if (!fp){
    cout << "======NO config File. Switch to default setting!=======" <<endl;
    return;
  }
  string line;
  vector<string> tokens;
  cout << "Reading config file..." <<endl;
  while(getline(fp,line)){
    stringstream ss(line);
    string token;
    while (getline(ss,token,' ')){
      if (token.size()!=0)
        tokens.push_back(token);
    }
    if (tokens.size()>=1 && tokens[0].find("#")==0){
      tokens.clear();
    }
    if (tokens.size()>=3){
      mapd[tokens[0]] = stod(tokens[2]);
    }
    tokens.clear();
  }
  return;
}
