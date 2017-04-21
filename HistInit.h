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
    TH1 *h_SRATT[4];
    TH1 *h_SRATT_NEW[4];
    TH1 *h_CRZ[4];
    // CRZ
    TH1 *h_CRZ_njetr[4];
    TH1 *h_CRZ_eptr[4];
    TH1 *h_CRZ_mptr[4];
    TH1 *h_CRZ_jet1r[4];
    TH1 *h_CRZ_jet2r[4];
    TH1 *h_CRZ_jet3r[4];
    TH1 *h_CRZ_jet4r[4];
    TH1 *h_CRZ_jet_121r[4];
    TH1 *h_CRZ_jet_122r[4];
    TH1 *h_CRZ_jet_08r[4];
    TH1 *h_CRZ_Zr[4];
    TH1 *h_CRZ_METr[4];
    TH1 *h_CRZ_METpr[4];
    TH1 *h_CRZ_nbjetr[4];
    TH1 *h_CRZ_njet[4];
    TH1 *h_CRZ_ept[4];
    TH1 *h_CRZ_mpt[4];
    TH1 *h_CRZ_jet1[4];
    TH1 *h_CRZ_jet2[4];
    TH1 *h_CRZ_jet3[4];
    TH1 *h_CRZ_jet4[4];
    TH1 *h_CRZ_jet_121[4];
    TH1 *h_CRZ_jet_122[4];
    TH1 *h_CRZ_jet_08[4];
    TH1 *h_CRZ_Z[4];
    TH1 *h_CRZ_MET[4];
    TH1 *h_CRZ_METp[4];
    TH1 *h_CRZ_nbjet[4];
    TH1 *h_CRZ_F_Z[4];
    TH1 *h_CRZ_F_MET[4];
    TH1 *h_CRZ_F_METp[4];
    TH1 *h_CRZ_F_nbjet[4];
    TH1 *h_CRZ_Zpt[4];
    TH1 *h_CRZ_F_Zpt[4]; 
    TH1 *h_CRZ_ne[4];
    TH1 *h_CRZ_nm[4];
    // CRT
    TH1 *h_CRT_08e[4];
    TH1 *h_CRT_08m[4];
    TH1 *h_CRT_09[4]; 
    TH1 *h_CRT_101[4];
    TH1 *h_CRT_102[4];
    TH1 *h_CRT_103[4];
    TH1 *h_CRT_104[4];
    TH1 *h_CRT_11[4]; 
    TH1 *h_CRT_12[4]; 
    TH1 *h_CRT_15[4];
    TH1 *h_CRT_16[4]; 
    // SRA
    TH1 *h_cut11[4]; 
    TH1 *h_cut12[4]; 
    TH1 *h_cut13[4]; 
    TH1 *h_cut151[4];
    TH1 *h_cut152[4];
    TH1 *h_cut153[4];
    TH1 *h_cut154[4];
    TH1 *h_cut16[4]; 
    TH1 *h_cut17[4]; 
    TH1 *h_cut18[4];
    TH1 *h_cut19[4];
    TH1 *h_cut20[4];
    TH1 *h_cut21[4];
    TH1 *h_cut22[4];
    TH1 *h_cut23[4];
    TH1 *h_cut24[4];
    TH1 *h_cut25[4];
    TH1 *h_cut27[4]; 
    TH1 *h_mass[4];
    TH1 *h_mass_e[4];
    TH1 *h_mass_mu[4];
    TH1 *h_mass_jet[4];

    TH1* NEWTH1F(string name, string title, int32_t bins, double fbin, double lbin); // type = 1
    TH1* NEWTH1D(string name, string title, int32_t bins, double fbin, double lbin);// type = 2
    TH1* NEWTH1I(string name, string title, int32_t bins, double fbin, double lbin);// type = 3
    void Init(TH1 *(*h)[4],struct HistInfo HI);
    TH1** GH(TH1 *(*h)[4]);
    vector<string> base;
    vector<pair<TH1*(*)[4],struct HistInfo>> histinfo;
    struct HistInfo SetInfo(string name, string title, int32_t bins, double fbin, double lbin, int32_t type);
    vector<TH1*> hall;
    void SetTFile(TFile *file);

};

HistInit::HistInit():fout(0){
  //base = {"_2015","_2016_1","_2016_2","_2016_3"};
  base = {"_2015","_2016_1"};
  Init(&h_CRZ, SetInfo("CRZ","Cut Flow of CRZ",14,0.5,14.5,1));
  Init(&h_CRZ_njetr, SetInfo("CRZnjetr", "number of jet reference;Njet", 15,0.5,15.5,1));
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
  Init(&h_CRZ_nbjetr, SetInfo("CRZnbjetr", "number of bjet reference;Nbjet", 10,0.5,10.5,1));
  Init(&h_CRZ_njet, SetInfo("CRZnjet", "number of jet;Njet", 15,0.5,15.5,1));
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
  Init(&h_CRZ_nbjet, SetInfo("CRZnbjet", "number of bjet;Nbjet", 10,0.5,10.5,1));
  Init(&h_CRZ_F_Z, SetInfo("CRZZ_F", "Zmass ;m [GeV];Events/0.2GeV", 50,86,96,1));
  Init(&h_CRZ_F_MET, SetInfo("CRZMET_F", "MET;ETmiss [GeV];Events/20GeV", 20,0,400,1));
  Init(&h_CRZ_F_METp, SetInfo("CRZMETp_F", "MET(p) reference;ETmiss [GeV];Events/20GeV", 40,0,800,1));
  Init(&h_CRZ_F_nbjet, SetInfo("CRZnbjet_F", "number of bjet;Nbjet", 10,0.5,10.5,1));
  Init(&h_CRZ_Zpt, SetInfo("CRZZpt","Z pt;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_F_Zpt, SetInfo("CRZZpt_F","Z pt;pt[GeV];Events/20GeV",30,0,600,1));
  Init(&h_CRZ_ne, SetInfo("CRZne","electron size;N. electron;Events/1",11,-0.5,10.5,1));
  Init(&h_CRZ_nm, SetInfo("CRZnm","muon size;N. muon;Events/1",11,-0.5,10.5,1));
  // CRT
  Init(&h_CRT_08e, SetInfo("CRT_08e","electron pt;pt[GeV];Events/20GeV",15,0,300,1));
  Init(&h_CRT_08m, SetInfo("CRT_08m","muon pt;pt[GeV];Events/20GeV",15,0,300,1));
  Init(&h_CRT_09, SetInfo("CRT_09", "number of jet;Njet", 15,0.5,15.5,1));
  Init(&h_CRT_101, SetInfo("CRT_101","Jet 1 pt;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRT_102, SetInfo("CRT_102","Jet 2 pt;pt[GeV];Events/40GeV",15,0,600,1));
  Init(&h_CRT_103, SetInfo("CRT_103","Jet 3 pt;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRT_104, SetInfo("CRT_104","Jet 4 pt;pt[GeV];Events/40GeV",10,0,400,1));
  Init(&h_CRT_11, SetInfo("CRT_11", "MET;ETmiss [GeV];Events/20GeV", 50,0,1000,1));
  Init(&h_CRT_12, SetInfo("CRT_12", "number of bjet;Nbjet", 10,0.5,10.5,1));
  Init(&h_CRT_15, SetInfo("CRT_15", "Makt12,0 > 70GeV;Makt12,0 [GeV]", 26,0,390,1));
  Init(&h_CRT_16, SetInfo("CRT_16", "MtBMin > 100GeV;MtBMin [GeV]", 20,0,1000,1));
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
  Init(&h_cut22, SetInfo("cut22", "Makt12,0 > 120GeV;Makt12,0 [GeV]", 26,0,390,1));
  Init(&h_cut23, SetInfo("cut23", "Makt12,1 > 120GeV;Makt12,1 [GeV]", 26,0,390,1));
  Init(&h_cut24, SetInfo("cut24", "Makt08 > 60GeV;Makt08,0 [GeV]", 15,0,450,1));
  Init(&h_cut25, SetInfo("cut25", "MtBMin > 200GeV;MtBMin [GeV]", 20,0,1000,1));
  Init(&h_cut27, SetInfo("cut27", "MET > 400GeV;ETmiss [GeV]", 20,0,1000,1)); 
  Init(&h_SRATT, SetInfo("SRATT_O","Cut Flow of SRA-TT",21,0,21,1));
  Init(&h_SRATT_NEW, SetInfo("SRATT_NEW","Cut Flow of New SRA-TT",20,0,20,1));
  // interest
  Init(&h_mass, SetInfo("MASS","MASS SPECTRUM;mass [GeV]",10000,0,1000,1));
  Init(&h_mass_e, SetInfo("MASS_E","MASS SPECTRUM;e mass [GeV]",1000,0,10,1));
  Init(&h_mass_mu, SetInfo("MASS_MU","MASS SPECTRUM;mu mass [GeV]",1000,0,10,1));
  Init(&h_mass_jet, SetInfo("MASS_JET","MASS SPECTRUM;jet mass [GeV]",10000,0,1000,1));

}

TH1* HistInit::NEWTH1F(string name, string title, int32_t bins, double fbin, double lbin){
  fout->cd();
  TH1 *hh = new TH1F(name.c_str(),title.c_str(),bins,fbin,lbin);
  //hh->SetDirectory(0);
  hh->Sumw2();
  hall.push_back(hh);
  return hh;
}

TH1* HistInit::NEWTH1D(string name, string title, int32_t bins, double fbin, double lbin){
  fout->cd();
  TH1 *hh = new TH1D(name.c_str(),title.c_str(),bins,fbin,lbin);
  //hh->SetDirectory(fout);
  hh->Sumw2();
  hall.push_back(hh);
  return hh;
}

TH1* HistInit::NEWTH1I(string name, string title, int32_t bins, double fbin, double lbin){
  fout->cd();
  TH1 *hh = new TH1I(name.c_str(),title.c_str(),bins,fbin,lbin);
  //hh->SetDirectory(fout);
  hh->Sumw2();
  hall.push_back(hh);
  return hh;
}

void HistInit::Init(TH1 *(*h)[4], HistInit::HistInfo HI){
  for(size_t i = 0; i < base.size(); i++){
    (*h)[i] = 0;
  }
  histinfo.push_back(make_pair(h,HI));
  return;
}

TH1** HistInit::GH(TH1 *(*h)[4]){
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
