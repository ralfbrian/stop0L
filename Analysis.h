#ifndef _PREANALYSIS_H_
#include "PreAnalysis.h"
#endif
#ifndef _SRA_H__
#include "SRA.h"
#endif

class Analysis{
  public:
    Analysis();
    //void Readconfig();
    //template<class X>
    //void Assign(X *region);
    void SetPA(PreAnalysis *Ana);
    void SetAC(ApplyCut *Ac);
    void SetupCutflow(std::vector<bool>& cf);
    void runSRA(SRA *sra);
    void runNEWSRA(SRA *sra);
    void runCRZ(CRZ *crz);
    void runNEWCRZ(CRZ *crz);
    void runCRT(OneLepton *crt);
    void runNEWCRT(OneLepton *crt);
    void runBjets(SRA *sra);
    int32_t n0;
    int32_t n1;
    int32_t n2;
  private:
    PreAnalysis *A;
    ApplyCut *ac;
  private:
    //map<string,double> mapd;
    //map<string,uint32_t> mapi;
};

Analysis::Analysis(){
  n0 = 0;
  n1 = 0;
  n2 = 0;
  //Readconfig();
}

void Analysis::SetPA(PreAnalysis *Ana){
  A = Ana;
  return;
}
/*
void Analysis::Readconfig(){
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
*/
void Analysis::SetAC(ApplyCut *Ac){
  ac = Ac;
  return;
}
/*
template<class X>
void Analysis::Assign(X *region){
  for(map<string,double>::iterator itr = mapd.begin();itr!=mapd.end();itr++){
    if (typeid(*region)==typeid(SRA)){
      if ((*itr).first.find("SRA_makt120")==0) region->makt120 = (*itr).second;
      if ((*itr).first.find("SRA_makt121")==0) region->makt121 = (*itr).second;
      if ((*itr).first.find("SRA_makt080")==0) region->makt080 = (*itr).second;
      if ((*itr).first.find("SRA_mtbmin")==0) region->mtbmin = (*itr).second;
      if ((*itr).first.find("SRA_met400")==0) region->met400 = (*itr).second;
    } 
    if (typeid(*region)==typeid(CRZ)){
    } 
    if (typeid(*region)==typeid(OneLepton)){
    } 
  }
  return;
}
*/
void Analysis::SetupCutflow(std::vector<bool>& cf){
  if (!A->isMC){
    for (size_t i = 0; i < A->base.size(); i++){
      cf[i] = false;
      if ((A->base[i].find("2015")!=std::string::npos && A->d15) ||
          (A->base[i].find("2016_1")!=std::string::npos && A->d16 && A->runnumber <= 302872) ||
          (A->base[i].find("2016_2")!=std::string::npos && A->d16 && A->runnumber <= 304008 && A->runnumber >= 302919) ||
          (A->base[i].find("2016_3")!=std::string::npos && A->d16 && A->runnumber >= 304128) 
         ){
        cf[i] = true;
      }
    }
  }
  return;
}

void Analysis::runBjets(SRA *sra){
  SetupCutflow(sra->cutflow);
  ac->AddCut(sra->cut_grl().
              cut_dstatus().
              cut_trigger().
              cut_PV().
              cut_cleanjet().
              cut_cosmics().
              cut_cleanmuon().
              cut_leptonveto().
              cut_jetpt().cut_met().pass());
  ac->FillHist(A->GH(&A->h_check_bjet),A->jet_idx[PreAnalysis::btagged],&(sra->w));
  sra->reset();
  return;
}

void Analysis::runSRA(SRA *sra){
  //Assign(sra);
  SetupCutflow(sra->cutflow);

  ac->AddCut(sra->cut_grl().pass());
  ac->FillHist(A->GH(&A->h_SRATT),0.5,&(sra->w));

  ac->AddCut(sra->cut_dstatus().pass());
  ac->FillHist(A->GH(&A->h_SRATT),1.5,&(sra->w));

  ac->AddCut(sra->cut_trigger().pass());
  ac->FillHist(A->GH(&A->h_SRATT),2.5,&(sra->w));

  ac->AddCut(sra->cut_PV().pass());
  ac->FillHist(A->GH(&A->h_SRATT),3.5,&(sra->w));

  ac->AddCut(sra->cut_cleanjet().pass());
  ac->FillHist(A->GH(&A->h_SRATT),4.5,&(sra->w));

  ac->AddCut(sra->cut_cosmics().pass());
  ac->FillHist(A->GH(&A->h_SRATT),5.5,&(sra->w));

  ac->AddCut(sra->cut_cleanmuon().pass());
  ac->FillHist(A->GH(&A->h_SRATT),6.5,&(sra->w));
  //for (auto& mu : A->getbasemu) ac->FillHist(A->GH(&A->h_mass_mu),mu->yylm*0.001,&(sra->w));
  //for (auto& el : A->getbaseele) ac->FillHist(A->GH(&A->h_mass_e),el->yylm*0.001,&(sra->w));

  ac->AddCut(sra->cut_leptonveto().pass());
  ac->FillHist(A->GH(&A->h_SRATT),7.5,&(sra->w));

  ac->AddCut(sra->cut_jetpt().pass());
  ac->FillHist(A->GH(&A->h_SRATT),8.5,&(sra->w));
  if (A->jetptsort.size()>=4){
    ac->FillHist(A->GH(&A->h_cut151),A->jetptsort[0],&(sra->w));
    ac->FillHist(A->GH(&A->h_cut152),A->jetptsort[1],&(sra->w));
    ac->FillHist(A->GH(&A->h_cut153),A->jetptsort[2],&(sra->w));
    ac->FillHist(A->GH(&A->h_cut154),A->jetptsort[3],&(sra->w));
  }

  ac->AddCut(sra->cut_oneb().pass());
  ac->FillHist(A->GH(&A->h_SRATT),9.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut16),A->jet_idx[PreAnalysis::btagged],&(sra->w));

  ac->AddCut(sra->cut_met().pass());
  ac->FillHist(A->GH(&A->h_SRATT),10.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut17),A->tst.Pt()*0.001,&(sra->w));
  //for (auto& jet : A->getbasejet) ac->FillHist(A->GH(&A->h_mass_jet),jet->yylm*0.001,&(sra->w));

  ac->AddCut(sra->cut_dpmj3().pass());
  ac->FillHist(A->GH(&A->h_SRATT),11.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut18),sra->dp2m1j,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut18),sra->dp2m2j,&(sra->w));

  ac->AddCut(sra->cut_trk().pass());
  ac->FillHist(A->GH(&A->h_SRATT),12.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut19),A->trk.Pt()*0.001,&(sra->w));

  ac->AddCut(sra->cut_dptrttst().pass());
  ac->FillHist(A->GH(&A->h_SRATT),13.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut20),sra->dpmm,&(sra->w));

  ac->AddCut(sra->cut_twob().pass());
  ac->FillHist(A->GH(&A->h_SRATT),14.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut21),A->jet_idx[PreAnalysis::btagged],&(sra->w));

  ac->AddCut(sra->cut_makt120().pass());
  ac->FillHist(A->GH(&A->h_SRATT),15.5,&(sra->w));
  if (A->j12sort.size()>0){
    ac->FillHist(A->GH(&A->h_cut22),A->j12sort[0]*0.001,&(sra->w));
  }

  ac->AddCut(sra->cut_makt121().pass());
  ac->FillHist(A->GH(&A->h_SRATT),16.5,&(sra->w));
  if (A->j12sort.size()>1){
    ac->FillHist(A->GH(&A->h_cut23),A->j12sort[1]*0.001,&(sra->w));
  }

  ac->AddCut(sra->cut_makt080().pass());
  ac->FillHist(A->GH(&A->h_SRATT),17.5,&(sra->w));
  if (A->j08sort.size()>0){
    ac->FillHist(A->GH(&A->h_cut24),A->j08sort[0]*0.001,&(sra->w));
  }

  ac->AddCut(sra->cut_mtbmin().pass());
  ac->FillHist(A->GH(&A->h_SRATT),18.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut25),A->Mtbmin*0.001,&(sra->w));

  ac->AddCut(sra->cut_tauveto().pass());
  ac->FillHist(A->GH(&A->h_SRATT),19.5,&(sra->w));

  ac->AddCut(sra->cut_met400().pass());
  ac->FillHist(A->GH(&A->h_SRATT),20.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut27),A->tst.Pt()*0.001,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut28),A->dRbb,&(sra->w));
  if (A->j12sort.size()>1){
    ac->FillHist(A->GH(&A->h_cut22_F),A->j12sort[0]*0.001,&(sra->w));
    ac->FillHist(A->GH(&A->h_cut23_F),A->j12sort[1]*0.001,&(sra->w));
  }

  ac->AddCut(sra->cut_drbb().pass());
  ac->FillHist(A->GH(&A->h_SRATT),21.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut28),A->dRbb,&(sra->w));
  A->CalMT2();
  ac->FillHist(A->GH(&A->h_cut29_p),A->MT2,&(sra->w));

  ac->AddCut(sra->cut_MT2().pass());
  ac->FillHist(A->GH(&A->h_SRATT),22.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut29),A->MT2,&(sra->w));
  
  sra->reset();
  return;
}

void Analysis::runNEWSRA(SRA *sra){
  SetupCutflow(sra->cutflow);

  ac->AddCut(sra->cut_grl().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),0.5,&(sra->w));

  ac->AddCut(sra->cut_dstatus().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),1.5,&(sra->w));

  ac->AddCut(sra->cut_trigger().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),2.5,&(sra->w));

  ac->AddCut(sra->cut_PV().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),3.5,&(sra->w));

  ac->AddCut(sra->cut_cleanjet().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),4.5,&(sra->w));

  ac->AddCut(sra->cut_cosmics().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),5.5,&(sra->w));

  ac->AddCut(sra->cut_cleanmuon().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),6.5,&(sra->w));

  ac->AddCut(sra->cut_leptonveto().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),7.5,&(sra->w));

  ac->AddCut(sra->cut_jetpt().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),8.5,&(sra->w));

  for (size_t newcut = 0; newcut < sra->cutflow.size(); newcut++){
    if (A->jet_idx[PreAnalysis::btagged]!=1){
      sra->cutflow[newcut] = false;
    }
  }
  ac->AddCut(sra->pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),9.5,&(sra->w));

  ac->AddCut(sra->cut_met().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),10.5,&(sra->w));

  ac->AddCut(sra->cut_dpmj3().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),11.5,&(sra->w));

  ac->AddCut(sra->cut_trk().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),12.5,&(sra->w));

  ac->AddCut(sra->cut_dptrttst().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),13.5,&(sra->w));

  ac->AddCut(sra->cut_makt120().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),14.5,&(sra->w));
  if (A->j12sort.size()>0){
    ac->FillHist(A->GH(&A->h_cut22_new),A->j12sort[0]*0.001,&(sra->w));
  }

  ac->AddCut(sra->cut_makt121().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),15.5,&(sra->w));
  if (A->j12sort.size()>1){
    ac->FillHist(A->GH(&A->h_cut23_new),A->j12sort[1]*0.001,&(sra->w));
  }

  ac->AddCut(sra->cut_makt080().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),16.5,&(sra->w));
  if (A->j08sort.size()>0){
    ac->FillHist(A->GH(&A->h_cut24_new),A->j08sort[0]*0.001,&(sra->w));
  }

  ac->AddCut(sra->cut_mtbmin().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),17.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut25_new),A->Mtbmin*0.001,&(sra->w));

  ac->AddCut(sra->cut_tauveto().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),18.5,&(sra->w));

  ac->AddCut(sra->cut_met400().pass());
  ac->FillHist(A->GH(&A->h_SRATT_NEW),19.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut27_new),A->tst.Pt()*0.001,&(sra->w));
  if (A->j12sort.size()>1){
    ac->FillHist(A->GH(&A->h_cut22_new_F),A->j12sort[0]*0.001,&(sra->w));
    ac->FillHist(A->GH(&A->h_cut23_new_F),A->j12sort[1]*0.001,&(sra->w));
  }
  
  sra->reset();
  return;
}

void Analysis::runCRZ(CRZ *crz){
  SetupCutflow(crz->cutflow);

  ac->AddCut(crz->cut_grl().pass());
  ac->FillHist(A->GH(&A->h_CRZ),1,&(crz->w));

  ac->AddCut(crz->cut_dstatus().pass());
  ac->FillHist(A->GH(&A->h_CRZ),2,&(crz->w));

  ac->AddCut(crz->cut_trigger().pass());
  ac->FillHist(A->GH(&A->h_CRZ),3,&(crz->w));

  //ac->AddCut(crz->cut_trigger_benchmark().pass());
  //ac->FillHist(A->GH(&A->h_CRZ),3,&(crz->w));

  ac->AddCut(crz->cut_PV().pass());
  ac->FillHist(A->GH(&A->h_CRZ),4,&(crz->w));

  ac->AddCut(crz->cut_cleanjet().pass());
  ac->FillHist(A->GH(&A->h_CRZ),5,&(crz->w));

  ac->AddCut(crz->cut_cosmics().pass());
  ac->FillHist(A->GH(&A->h_CRZ),6,&(crz->w));

  ac->AddCut(crz->cut_cleanmuon().pass());
  ac->FillHist(A->GH(&A->h_CRZ),7,&(crz->w));
    
  ac->AddCut(crz->cut_nlepton().pass());
  ac->FillHist(A->GH(&A->h_CRZ),8,&(crz->w));

  ac->AddCut(crz->cut_lpt().pass());
  ac->FillHist(A->GH(&A->h_CRZ),9,&(crz->w));
  for(auto& el : A->getbaseele) ac->FillHist(A->GH(&A->h_CRZ_ept),el->yylpt*0.001,&(crz->w));
  for(auto& mu : A->getbasemu) ac->FillHist(A->GH(&A->h_CRZ_mpt),mu->yylpt*0.001,&(crz->w));

  ac->AddCut(crz->cut_met().pass());
  ac->FillHist(A->GH(&A->h_CRZ),10,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ_MET),A->tst.Pt()*0.001,&(crz->w));

  ac->AddCut(crz->cut_njet().pass());
  ac->FillHist(A->GH(&A->h_CRZ),11,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ_njet),A->getbasejet.size(),&(crz->w));

  ac->AddCut(crz->cut_jpt().pass());
  ac->FillHist(A->GH(&A->h_CRZ),12,&(crz->w));

  ac->AddCut(crz->cut_zmass().pass());
  ac->FillHist(A->GH(&A->h_CRZ),13,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ_Z),crz->getZmass()*0.001,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ_Zpt),crz->getZpt()*0.001,&(crz->w));

  ac->AddCut(crz->cut_metp().pass());
  ac->FillHist(A->GH(&A->h_CRZ),14,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ_METp),crz->getmetp()*0.001,&(crz->w));

  ac->AddCut(crz->cut_btag().pass());
  ac->FillHist(A->GH(&A->h_CRZ),15,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ),A->jet_idx[PreAnalysis::btagged],&(crz->w));

  ac->AddCut(crz->cut_makt120().pass());
  ac->FillHist(A->GH(&A->h_CRZ),16,&(crz->w));
  if (A->j12sort.size()>0){
    ac->FillHist(A->GH(&A->h_CRZ_makt120),A->j12sort[0]*0.001,&(crz->w));
  }

  ac->AddCut(crz->cut_makt121().pass());
  ac->FillHist(A->GH(&A->h_CRZ),17,&(crz->w));
  if (A->j12sort.size()>1){
    ac->FillHist(A->GH(&A->h_CRZ_makt121),A->j12sort[1]*0.001,&(crz->w));
  }

  crz->reset();
  return;
}

void Analysis::runNEWCRZ(CRZ *crz){
  SetupCutflow(crz->cutflow);

  ac->AddCut(crz->cut_grl().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),1,&(crz->w));

  ac->AddCut(crz->cut_dstatus().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),2,&(crz->w));

  ac->AddCut(crz->cut_trigger().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),3,&(crz->w));

  ac->AddCut(crz->cut_PV().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),4,&(crz->w));

  ac->AddCut(crz->cut_cleanjet().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),5,&(crz->w));

  ac->AddCut(crz->cut_cosmics().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),6,&(crz->w));

  ac->AddCut(crz->cut_cleanmuon().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),7,&(crz->w));
    
  ac->AddCut(crz->cut_nlepton().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),8,&(crz->w));

  ac->AddCut(crz->cut_lpt().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),9,&(crz->w));

  ac->AddCut(crz->cut_met().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),10,&(crz->w));

  ac->AddCut(crz->cut_njet().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),11,&(crz->w));

  ac->AddCut(crz->cut_jpt().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),12,&(crz->w));

  ac->AddCut(crz->cut_zmass().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),13,&(crz->w));

  ac->AddCut(crz->cut_metp().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),14,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ_nbjetr_NEW),A->jet_idx[PreAnalysis::btagged],&(crz->w));

  ac->AddCut(crz->cut_onebtag().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),15,&(crz->w));
  ac->FillHist(A->GH(&A->h_CRZ_nbjet_NEW),A->jet_idx[PreAnalysis::btagged],&(crz->w));

  ac->AddCut(crz->cut_makt120().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),16,&(crz->w));
  if (A->j12sort.size()>0){
    ac->FillHist(A->GH(&A->h_CRZ_makt120_1b),A->j12sort[0]*0.001,&(crz->w));
  }

  ac->AddCut(crz->cut_makt121().pass());
  ac->FillHist(A->GH(&A->h_CRZ_NEW),17,&(crz->w));
  if (A->j12sort.size()>1){
    ac->FillHist(A->GH(&A->h_CRZ_makt121_1b),A->j12sort[1]*0.001,&(crz->w));
  }

  crz->reset();
  return;
}

void Analysis::runCRT(OneLepton *crt){
  SetupCutflow(crt->cutflow);

  ac->AddCut(crt->cut_GRL().pass());
  ac->FillHist(A->GH(&A->h_CRT),1,&(crt->w));

  ac->AddCut(crt->cut_LTerror().pass());
  ac->FillHist(A->GH(&A->h_CRT),2,&(crt->w));

  ac->AddCut(crt->cut_trigger().pass());
  ac->FillHist(A->GH(&A->h_CRT),3,&(crt->w));

  ac->AddCut(crt->cut_PV().pass());
  ac->FillHist(A->GH(&A->h_CRT),4,&(crt->w));

  ac->AddCut(crt->cut_cleanjets().pass());
  ac->FillHist(A->GH(&A->h_CRT),5,&(crt->w));

  ac->AddCut(crt->cut_cosmics().pass());
  ac->FillHist(A->GH(&A->h_CRT),6,&(crt->w));

  ac->AddCut(crt->cut_cleanmuons().pass());
  ac->FillHist(A->GH(&A->h_CRT),7,&(crt->w));

  ac->AddCut(crt->cut_nlepton().pass());
  ac->FillHist(A->GH(&A->h_CRT),8,&(crt->w));

  //ac->AddCut(crt->cut_lpt().pass());
  //ac->FillHist(A->GH(&A->h_CRT),9,&(crt->w));
  if (A->mu_idx[PreAnalysis::baseline]==1) ac->FillHist(A->GH(&A->h_CRT_08m),crt->lep1.Pt()*0.001,&(crt->w));
  if (A->ele_idx[PreAnalysis::baseline]==1) ac->FillHist(A->GH(&A->h_CRT_08e),crt->lep1.Pt()*0.001,&(crt->w));

  ac->AddCut(crt->cut_njet().pass());
  ac->FillHist(A->GH(&A->h_CRT),9,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_09),A->jet_idx[PreAnalysis::baseline],&(crt->w));

  ac->AddCut(crt->cut_jpt().pass());
  ac->FillHist(A->GH(&A->h_CRT),10,&(crt->w));
  if (A->jetptsort.size()>=4){
    ac->FillHist(A->GH(&A->h_CRT_101),A->jetptsort.at(0),&(crt->w));
    ac->FillHist(A->GH(&A->h_CRT_102),A->jetptsort.at(1),&(crt->w));
    ac->FillHist(A->GH(&A->h_CRT_103),A->jetptsort.at(2),&(crt->w));
    ac->FillHist(A->GH(&A->h_CRT_104),A->jetptsort.at(3),&(crt->w));
  }

 // ac->AddCut(crt->cut_onebtags().pass());

  ac->AddCut(crt->cut_dp2mj2().pass());
  ac->FillHist(A->GH(&A->h_CRT),11,&(crt->w));

  ac->AddCut(crt->cut_MET250().pass());
  ac->FillHist(A->GH(&A->h_CRT),12,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_11),A->tst.Pt()*0.001,&(crt->w));

  //ac->AddCut(crt->cut_mtbmin().pass());
  //ac->AddCut(crt->cut_makt120_ben().pass());
  
  ac->AddCut(crt->cut_makt120().pass());
  ac->FillHist(A->GH(&A->h_CRT),13,&(crt->w));
  if(A->j12sort.size() >= 1 )ac->FillHist(A->GH(&A->h_CRT_15),A->j12sort.at(0)*0.001,&(crt->w));

  ac->AddCut(crt->cut_makt121().pass());
  ac->FillHist(A->GH(&A->h_CRT),14,&(crt->w));
  if(A->j12sort.size() >= 2 )ac->FillHist(A->GH(&A->h_CRT_15),A->j12sort.at(0)*0.001,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_14_p),crt->getmtlep(),&(crt->w));

  ac->AddCut(crt->cut_mtlep().pass());
  ac->FillHist(A->GH(&A->h_CRT),15,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_14),crt->getmtlep(),&(crt->w));

  ac->AddCut(crt->cut_twobtags().pass());
  ac->FillHist(A->GH(&A->h_CRT),16,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_12),A->jet_idx[PreAnalysis::btagged],&(crt->w));

  ac->AddCut(crt->cut_mtbmin().pass());
  ac->FillHist(A->GH(&A->h_CRT),17,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_16),A->Mtbmin,&(crt->w));

  ac->AddCut(crt->cut_drbl().pass());
  ac->FillHist(A->GH(&A->h_CRT),18,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_17),crt->getdrbl(),&(crt->w));

  ac->AddCut(crt->cut_makt080().pass());
  ac->FillHist(A->GH(&A->h_CRT),19,&(crt->w));
  if (A->j08sort.size() >= 1)
    ac->FillHist(A->GH(&A->h_CRT_Makt08),A->j08sort[0]*0.001,&(crt->w));

  ac->AddCut(crt->cut_drbb().pass());
  ac->FillHist(A->GH(&A->h_CRT),20,&(crt->w));
  ac->FillHist(A->GH(&A->h_CRT_MET_F),A->tst.Pt()*0.001,&(crt->w));

  crt->reset();
  return;
}
