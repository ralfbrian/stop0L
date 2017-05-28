#ifndef _PREANALYSIS_H_
#include "PreAnalysis.h"
#endif
#ifndef _SRA_H__
#include "SRA.h"
#endif

class Analysis{
  public:
    Analysis();
    void SetPA(PreAnalysis *Ana);
    void SetAC(ApplyCut *Ac);
    void runSRA(SRA *sra);
    void runNEWSRA(SRA *sra);
    void runCRZ(CRZ *crz);
    void runNEWCRZ(CRZ *crz);
    int32_t n0;
    int32_t n1;
    int32_t n2;
  private:
    PreAnalysis *A;
    ApplyCut *ac;
};

Analysis::Analysis(){
  n0 = 0;
  n1 = 0;
  n2 = 0;
}

void Analysis::SetPA(PreAnalysis *Ana){
  A = Ana;
  return;
}

void Analysis::SetAC(ApplyCut *Ac){
  ac = Ac;
  return;
}

void Analysis::runSRA(SRA *sra){
  if (!A->isMC){
    if (A->d15){
      sra->cutflow[1] = false;
    }
    else{
      sra->cutflow[0] = false;
    }
  }
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

  //ac->AddCut(sra->cut_oneb().pass());
  //ac->FillHist(A->GH(&A->h_SRATT),9.5,&(sra->w));
  //ac->FillHist(A->GH(&A->h_cut16),A->jet_idx[PreAnalysis::btagged],&(sra->w));

  ac->AddCut(sra->cut_met().pass());
  ac->FillHist(A->GH(&A->h_SRATT),10.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut17),A->tst.Pt()*0.001,&(sra->w));
  //for (auto& jet : A->getbasejet) ac->FillHist(A->GH(&A->h_mass_jet),jet->yylm*0.001,&(sra->w));
  if (sra->pass()[0]){
    if (A->jet_idx[PreAnalysis::btagged]==0){
      n0++;
    }
    else if (A->jet_idx[PreAnalysis::btagged]==1){
      n1++;
    }
    else{
      n2++;
    }
  }

  ac->AddCut(sra->cut_dpmj2().pass());
  ac->FillHist(A->GH(&A->h_SRATT),11.5,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut18),sra->dp2mfj,&(sra->w));
  ac->FillHist(A->GH(&A->h_cut18),sra->dp2msj,&(sra->w));

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
  if (A->j12sort.size()>1){
    ac->FillHist(A->GH(&A->h_cut22_F),A->j12sort[0]*0.001,&(sra->w));
    ac->FillHist(A->GH(&A->h_cut23_F),A->j12sort[1]*0.001,&(sra->w));
  }
  
  sra->reset();
  return;
}

void Analysis::runNEWSRA(SRA *sra){
  if (!A->isMC){
    if (A->d15){
      sra->cutflow[1] = false;
    }
    else{
      sra->cutflow[0] = false;
    }
  }
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

  ac->AddCut(sra->cut_dpmj2().pass());
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
  if (!A->isMC){
    if (A->d15){
      crz->cutflow[1] = false;
    }
    else{
      crz->cutflow[0] = false;
    }
  }

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

  crz->reset();
  return;
}

void Analysis::runNEWCRZ(CRZ *crz){
  if (!A->isMC){
    if (A->d15){
      crz->cutflow[1] = false;
    }
    else{
      crz->cutflow[0] = false;
    }
  }

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

  crz->reset();
  return;
}
