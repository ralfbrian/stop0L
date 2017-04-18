using namespace std;
class Histogram{
  public:
    Histogram();
    //cutflow
    TH1 *h_SRATT[4];
    TH1 *h_NEW_SRATT[4];
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

    TH1* NEWTH1F(const char* name, const char* title, int32_t bins, double fbin, double lbin);
    TH1* NEWTH1D(const char* name, const char* title, int32_t bins, double fbin, double lbin);
    TH1* NEWTH1I(const char* name, const char* title, int32_t bins, double fbin, double lbin);
    vector<string> base;
    vector<TH1*> h;

};

Histogram::Histogram(){
    base = {"_2015","2016_1","_2016_2","_2016_3"};
  for(size_t i = 0; i < base.size(); i++){
    h_CRZ_njetr[i] = NEWTH1I(("CRZnjetr" + base[i]).c_str(), "number of jet reference;Njet", 15,0.5,15.5);
    h_CRZ_eptr[i] = NEWTH1F(("CRZeptr" + base[i]).c_str(),"electron pt reference;pt[GeV];Events/20GeV",30,0,600);
    h_CRZ_mptr[i] = NEWTH1F(("CRZmptr" + base[i]).c_str(),"muon pt reference;pt[GeV];Events/20GeV",30,0,600);
    h_CRZ_jet1r[i] = NEWTH1F(("CRZjet1r" + base[i]).c_str(),"Jet 1 pt reference;pt[GeV];Events/40GeV",15,0,600);
    h_CRZ_jet2r[i] = NEWTH1F(("CRZjet2r" + base[i]).c_str(),"Jet 2 pt reference;pt[GeV];Events/40GeV",15,0,600);
    h_CRZ_jet3r[i] = NEWTH1F(("CRZjet3r" + base[i]).c_str(),"Jet 3 pt reference;pt[GeV];Events/40GeV",10,0,400);
    h_CRZ_jet4r[i] = NEWTH1F(("CRZjet4r" + base[i]).c_str(),"Jet 4 pt reference;pt[GeV];Events/40GeV",10,0,400);
    h_CRZ_jet_121r[i] = NEWTH1F(("CRZjet_121r" + base[i]).c_str(),"R1.2 jet 1 reference;m[GeV];Events/25GeV",18,0,450);
    h_CRZ_jet_122r[i] = NEWTH1F(("CRZjet_122r" + base[i]).c_str(),"R1.2 jet 2 reference;m[GeV];Events/25GeV",10,0,250);
    h_CRZ_jet_08r[i] = NEWTH1F(("CRZjet_08r" + base[i]).c_str(),"R0.8 jet 1 reference;m[GeV];Events/30GeV",15,0,450);
    h_CRZ_Zr[i] = NEWTH1F(("CRZZr" + base[i]).c_str(), "Zmass reference;m [GeV];Events/0.2GeV", 50,86,96);
    h_CRZ_METr[i] = NEWTH1F(("CRZMETr" + base[i]).c_str(), "MET reference;ETmiss [GeV];Events/20GeV", 20,0,400);
    h_CRZ_METpr[i] = NEWTH1F(("CRZMETpr" + base[i]).c_str(), "MET(p) reference;ETmiss [GeV];Events/20GeV", 40,0,800);
    h_CRZ_nbjetr[i] = NEWTH1I(("CRZnbjetr" + base[i]).c_str(), "number of bjet reference;Nbjet", 10,0.5,10.5);
    h_CRZ_njet[i] = NEWTH1I(("CRZnjet" + base[i]).c_str(), "number of jet;Njet", 15,0.5,15.5);
    h_CRZ_ept[i] = NEWTH1F(("CRZept" + base[i]).c_str(),"electron pt;pt[GeV];Events/20GeV",30,0,600);
    h_CRZ_mpt[i] = NEWTH1F(("CRZmpt" + base[i]).c_str(),"muon pt;pt[GeV];Events/20GeV",30,0,600);
    h_CRZ_jet1[i] = NEWTH1F(("CRZjet1" + base[i]).c_str(),"Jet 1 pt;pt[GeV];Events/40GeV",15,0,600);
    h_CRZ_jet2[i] = NEWTH1F(("CRZjet2" + base[i]).c_str(),"Jet 2 pt;pt[GeV];Events/40GeV",15,0,600);
    h_CRZ_jet3[i] = NEWTH1F(("CRZjet3" + base[i]).c_str(),"Jet 3 pt;pt[GeV];Events/40GeV",10,0,400);
    h_CRZ_jet4[i] = NEWTH1F(("CRZjet4" + base[i]).c_str(),"Jet 4 pt;pt[GeV];Events/40GeV",10,0,400);
    h_CRZ_jet_121[i] = NEWTH1F(("CRZjet_121" + base[i]).c_str(),"R1.2 jet 1;m[GeV];Events/25GeV",18,0,450);
    h_CRZ_jet_122[i] = NEWTH1F(("CRZjet_122" + base[i]).c_str(),"R1.2 jet 2;m[GeV];Events/25GeV",10,0,250);
    h_CRZ_jet_08[i] = NEWTH1F(("CRZjet_08" + base[i]).c_str(),"R0.8 jet 1;m[GeV];Events/30GeV",15,0,450);
    h_CRZ_Z[i] = NEWTH1F(("CRZZ" + base[i]).c_str(), "Zmass ;m [GeV];Events/0.2GeV", 50,86,96);
    h_CRZ_MET[i] = NEWTH1F(("CRZMET" + base[i]).c_str(), "MET;ETmiss [GeV];Events/20GeV", 20,0,400);
    h_CRZ_METp[i] = NEWTH1F(("CRZMETp" + base[i]).c_str(), "MET(p) reference;ETmiss [GeV];Events/20GeV", 40,0,800);
    h_CRZ_nbjet[i] = NEWTH1I(("CRZnbjet" + base[i]).c_str(), "number of bjet;Nbjet", 10,0.5,10.5);
    h_CRZ_F_Z[i] = NEWTH1F(("CRZZ_F" + base[i]).c_str(), "Zmass ;m [GeV];Events/0.2GeV", 50,86,96);
    h_CRZ_F_MET[i] = NEWTH1F(("CRZMET_F" + base[i]).c_str(), "MET;ETmiss [GeV];Events/20GeV", 20,0,400);
    h_CRZ_F_METp[i] = NEWTH1F(("CRZMETp_F" + base[i]).c_str(), "MET(p) reference;ETmiss [GeV];Events/20GeV", 40,0,800);
    h_CRZ_F_nbjet[i] = NEWTH1I(("CRZnbjet_F" + base[i]).c_str(), "number of bjet;Nbjet", 10,0.5,10.5);
    h_CRZ_Zpt[i] = NEWTH1F(("CRZZpt" + base[i]).c_str(),"Z pt;pt[GeV];Events/20GeV",30,0,600);
    h_CRZ_F_Zpt[i] = NEWTH1F(("CRZZpt_F" + base[i]).c_str(),"Z pt;pt[GeV];Events/20GeV",30,0,600);
    h_CRZ_ne[i] = NEWTH1F(("CRZne" + base[i]).c_str(),"electron size;N. electron;Events/1",11,-0.5,10.5);
    h_CRZ_nm[i] = NEWTH1F(("CRZnm" + base[i]).c_str(),"muon size;N. muon;Events/1",11,-0.5,10.5);
    // CRT
    h_CRT_08e[i] = NEWTH1F(("CRT_08e" + base[i]).c_str(),"electron pt;pt[GeV];Events/20GeV",15,0,300);
    h_CRT_08m[i] = NEWTH1F(("CRT_08m" + base[i]).c_str(),"muon pt;pt[GeV];Events/20GeV",15,0,300);
    h_CRT_09[i] = NEWTH1I(("CRT_09" + base[i]).c_str(), "number of jet;Njet", 15,0.5,15.5);
    h_CRT_101[i] = NEWTH1F(("CRT_101" + base[i]).c_str(),"Jet 1 pt;pt[GeV];Events/40GeV",15,0,600);
    h_CRT_102[i] = NEWTH1F(("CRT_102" + base[i]).c_str(),"Jet 2 pt;pt[GeV];Events/40GeV",15,0,600);
    h_CRT_103[i] = NEWTH1F(("CRT_103" + base[i]).c_str(),"Jet 3 pt;pt[GeV];Events/40GeV",10,0,400);
    h_CRT_104[i] = NEWTH1F(("CRT_104" + base[i]).c_str(),"Jet 4 pt;pt[GeV];Events/40GeV",10,0,400);
    h_CRT_11[i] = NEWTH1F(("CRT_11" + base[i]).c_str(), "MET;ETmiss [GeV];Events/20GeV", 50,0,1000);
    h_CRT_12[i] = NEWTH1I(("CRT_12" + base[i]).c_str(), "number of bjet;Nbjet", 10,0.5,10.5);
    h_CRT_15[i] = NEWTH1F(("CRT_15" + base[i]).c_str(), "Makt12,0 > 70GeV;Makt12,0 [GeV]", 26,0,390);
    h_CRT_16[i] = NEWTH1F(("CRT_16" + base[i]).c_str(), "MtBMin > 100GeV;MtBMin [GeV]", 20,0,1000);
    // SRA
    //h_cut11[i] = NEWTH1I(("cut11" + base[i]).c_str(), "clean jet;Njet", 15,0.5,15.5);
    //h_cut12[i] = NEWTH1I(("cut12" + base[i]).c_str(), "cosmic;Nmuon", 15,0.5,15.5);
    //h_cut13[i] = NEWTH1I(("cut13" + base[i]).c_str(), "clean muon;Nmuon", 15,0.5,15.5);
    h_cut151[i] = NEWTH1F(("cut151" + base[i]).c_str(), "signal jet (jet1);pT [GeV]", 10,0,400);
    h_cut152[i] = NEWTH1F(("cut152" + base[i]).c_str(), "signal jet (jet2);pT [GeV]", 10,0,400);
    h_cut153[i] = NEWTH1F(("cut153" + base[i]).c_str(), "signal jet (jet3);pT [GeV]", 10,0,400);
    h_cut154[i] = NEWTH1F(("cut154" + base[i]).c_str(), "signal jet (jet4);pT [GeV]", 10,0,400);
    h_cut16[i] = NEWTH1I(("cut16" + base[i]).c_str(), "one b jet;Nbtag", 6,0.5,6.5);
    h_cut17[i] = NEWTH1F(("cut17" + base[i]).c_str(), "MET > 250GeV ;ETmiss [GeV]", 20,0,1000);
    h_cut18[i] = NEWTH1F(("cut18" + base[i]).c_str(), "d phi between MET and 2 jets > 0.4 ;d phi [rad]", 40,0,4);
    h_cut19[i] = NEWTH1F(("cut19" + base[i]).c_str(), "MET track  > 30GeV;ETmiss,track [GeV]", 40,0,400);
    h_cut20[i] = NEWTH1F(("cut20" + base[i]).c_str(), "d phi between MET track and MET < pi/3;d phi [rad]", 32,0,3.2);
    h_cut21[i] = NEWTH1F(("cut21" + base[i]).c_str(), "two b tag;Nbtag", 6,0.5,6.5);
    h_cut22[i] = NEWTH1F(("cut22" + base[i]).c_str(), "Makt12,0 > 120GeV;Makt12,0 [GeV]", 26,0,390);
    h_cut23[i] = NEWTH1F(("cut23" + base[i]).c_str(), "Makt12,1 > 120GeV;Makt12,1 [GeV]", 26,0,390);
    h_cut24[i] = NEWTH1F(("cut24" + base[i]).c_str(), "Makt08 > 60GeV;Makt08,0 [GeV]", 15,0,450);
    h_cut25[i] = NEWTH1F(("cut25" + base[i]).c_str(), "MtBMin > 200GeV;MtBMin [GeV]", 20,0,1000);
    h_cut27[i] = NEWTH1F(("cut27" + base[i]).c_str(), "MET > 400GeV;ETmiss [GeV]", 20,0,1000); 
    h_SRATT[i] = NEWTH1F(("SRATT_O" + base[i]).c_str(),"Cut Flow of SRA-TT",21,0,21);
    h_NEW_SRATT[i] = NEWTH1F(("SRATT_NEW" + base[i]).c_str(),"Cut Flow of New SRA-TT",20,0,20);
  }
}

TH1* Histogram::NEWTH1F(const char* name, const char* title, int32_t bins, double fbin, double lbin){
  TH1 *hh = new TH1F(name,title,bins,fbin,lbin);
  //hh->SetDirectory(0);
  hh->Sumw2();
  h.push_back(hh);
  return hh;
}

TH1* Histogram::NEWTH1D(const char* name, const char* title, int32_t bins, double fbin, double lbin){
  TH1 *hh = new TH1D(name,title,bins,fbin,lbin);
  //hh->SetDirectory(0);
  hh->Sumw2();
  h.push_back(hh);
  return hh;
}

TH1* Histogram::NEWTH1I(const char* name, const char* title, int32_t bins, double fbin, double lbin){
  TH1 *hh = new TH1I(name,title,bins,fbin,lbin);
  //hh->SetDirectory(0);
  hh->Sumw2();
  h.push_back(hh);
  return hh;
}
