#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include "HistInit.h"
#include "datadefinition.h"
#include "PreAnalysis.h"


PreAnalysis::PreAnalysis():HistInit(),
                           jet_idx(8,0),
                           ele_idx(8,0),
                           mu_idx(8,0),
                           econtainer(0),
                           mcontainer(0),
                           pcontainer(0),
                           tcontainer(0),
                           jcontainer(0),
                           j12container(0),
                           j08container(0),
                           m_weight(0),
                           metcontainer(0),
                           fileinfo(0),
                           filename(""),
                           isSignal(false)
{
  firstjet.SetPxPyPzE(0.,0.,0.,0.);
  secondjet.SetPxPyPzE(0.,0.,0.,0.);
}

void PreAnalysis::init(){
  w0 = 1.;
  if (!isSignal)xsec = 1.;
  w1 = 1.;
  w2 = 1.;
  w3 = 1.;
  w4 = 1.;
  w5 = 1.;
  w6 = 1.;
  w7 = 1.;
  w8 = 1.;
  w9 = 1.;
  w10 = 1.;
  if(isMC){
    w0 = (*m_weight)["mcEventWeight"];
    if (!isSignal) xsec = (*m_weight)["Xsec"];
    w1 = (*m_weight)["PileupWeight"];
    w2 = (*m_weight)["JetSF"];
    w3 = (*m_weight)["JvtSF"];
    w4 = (*m_weight)["BtagSF"];
    w5 = (*m_weight)["ElectronSF2015"];
    w6 = (*m_weight)["ElectronSF2016"];
    w7 = (*m_weight)["MuonSF2015"];
    w8 = (*m_weight)["MuonSF2016"];
    w9 = (*m_weight)["SumofWeight"];
    w10 = (*m_weight)["SumofSquaredWeight"];
  }
  else{
    if((*fileinfo)["RunNumber"].size()>5){
      runnumber = stoi((*fileinfo)["RunNumber"]);
      eventnumber = stol((*fileinfo)["EventNumber"]);
    }
  }
  // MET
  metPro* met = &(*metcontainer)["TST"];
  tst.SetPxPyPzE(met->yylmpx,met->yylmpy,0.,met->yylsumet);
  met = &(*metcontainer)["TRK"];
  trk.SetPxPyPzE(met->yylmpx,met->yylmpy,0.,met->yylsumet);
  
  jetcal();
  
  ecal();
 
  mcal();

  taucal();
  
  j12cal();

  j08cal();
  
  passGRL = false;
  PV = false;
  Dstatus = false;
  if ((*fileinfo)["GRL"] == "isMC" || (*fileinfo)["GRL"] == "2015" || (*fileinfo)["GRL"] == "2016"){
    passGRL = true;
  }
  if ((*fileinfo)["Dstatus"] == "isMC" || (*fileinfo)["Dstatus"] != "false"){
    Dstatus = true;
  }
  if ((*fileinfo)["primaryvertex"] == "true"){
    PV = true;
  }
  return;
}

void PreAnalysis::jetcal(){
  bool firstrun = true;
  bool firstrun2= true;
  jet80 = 0;
  jet40 = 0;
  jet20 = 0;
  jetptsort.clear();
  getbasejet.clear();
  getbjet.clear();
  jet_idx.clear();
  jet_idx.resize(8,0);
  objPro *fjet = 0;
  objPro *sjet = 0;
  for (auto& jet : *jcontainer){
    if (jet.passOR == 1   && jet.baseline == 1){
      jet_idx[PreAnalysis::baseline]++;
      if (jet.bad == 1 ){
        jet_idx[PreAnalysis::bad]++;
      }
      else{
        if (jet.signal == 0)continue;
        jetptsort.push_back(jet.yylpt*0.001);
        getbasejet.push_back(&jet);
        if (jet.yylpt > 20000.) jet20++;
        if (jet.yylpt > 40000.){
          jet40++;
          if (jet.yylpt > 80000.){
            jet80++;
            if (firstrun){
              fjet = &jet;
              firstrun = false;
            }
            else{
              if (jet.yylpt > fjet->yylpt){
                    sjet = fjet;
                    fjet = &jet;
              }
              else{
                if (sjet != NULL){
                  if (jet.yylpt > sjet->yylpt) sjet = &jet;
                }
                else
                  sjet = &jet;
              }
            }
          }
        }
        if (jet.btagged == 1){
          getbjet.push_back(&jet);
          jet_idx[PreAnalysis::btagged]++;
          if (firstrun2){
            Mtbmin = sqrt(2*(jet.yylpt)*(tst.Pt())*(1-cos(jet.yylphi-tst.Phi())));
            fabsdphi = fabs(jet.yylphi-tst.Phi());
            firstrun2 = false;
          }
          else {
            if (fabs(jet.yylphi-tst.Phi()) < fabsdphi){
              Mtbmin = sqrt(2*(jet.yylpt)*(tst.Pt())*(1-cos(jet.yylphi-tst.Phi())));
              fabsdphi = fabs(jet.yylphi-tst.Phi());
            }
          }
        }
      }
    }
  }
  std::sort(jetptsort.rbegin(),jetptsort.rend());
  if (fjet)firstjet.SetPtEtaPhiM(fjet->yylpt,fjet->yyleta,fjet->yylphi,fjet->yylm);
  if (sjet)secondjet.SetPtEtaPhiM(sjet->yylpt,sjet->yyleta,sjet->yylphi,sjet->yylm);
  return;
}

void PreAnalysis::ecal(){
  getbaseele.clear();
  ele_idx.clear();
  ele_idx.resize(8,0);
  for (auto& electron : *econtainer){
    if (electron.passOR == 1) {
      if (electron.baseline == 1){
        getbaseele.push_back(&electron);
        ele_idx[PreAnalysis::baseline]++;
      }
      if (electron.signal == 1 ){
        ele_idx[PreAnalysis::signallep]++;
      } 
    }
  }
  return;
}

void PreAnalysis::mcal(){
  mu_idx.clear();
  mu_idx.resize(8,0);
  getbasemu.clear();
  for (auto& muon : *mcontainer){
    if (muon.bad == 1 && muon.baseline == 1)
      mu_idx[PreAnalysis::bad]++;
    if (muon.passOR == 1 && muon.baseline == 1){
      mu_idx[PreAnalysis::baseline]++;
      getbasemu.push_back(&muon);
      if (muon.cosmic == 1 )
        mu_idx[PreAnalysis::cosmic]++;
    }
  }
  return;
}

void PreAnalysis::j12cal(){
  j12sort.clear();
  for (auto& jet12 : *j12container){
    j12sort.push_back(jet12.yylm);
  }
  std::sort(j12sort.rbegin(),j12sort.rend());
  return;
}

void PreAnalysis::j08cal(){
  j08sort.clear();
  for (auto& jet08 : *j08container){
    j08sort.push_back(jet08.yylm);
  }
  std::sort(j08sort.rbegin(),j08sort.rend());
  return;
}

void PreAnalysis::taucal(){
  tau_idx.clear();
  tau_idx.resize(8,0);
/*
  for (auto& tau : *tcontainer){
    if (tau.passOR == 1 && tau.baseline ==1)
      tau_idx[PreAnalysis::baseline]++;
  } 
*/
  objPro *taujet = 0;
  float minphi = std::numeric_limits<float>::max();
  for (auto& jet : getbasejet){
    if (fabs(jet->yylphi-tst.Phi()) < minphi &&jet->btagged == 0 && jet->yylnTrk <= 4 && jet->yyleta < 2.5){
      taujet = jet;
      minphi = fabs(jet->yylphi-tst.Phi());
    }
  }
  if (taujet && minphi < std::atan(1.0)*4/5){
    tau_idx[PreAnalysis::baseline]++;
  }
  return;
}

void PreAnalysis::SetTree(TTree *Yggdrasil){
  Yggdrasil->SetBranchAddress("fileinfo",&fileinfo);
  Yggdrasil->GetEntry();
  if ((*fileinfo)["GRL"] == "isMC"){
    isMC = true;
    Yggdrasil->SetBranchAddress("weight",&m_weight);
  }
  else{
    isMC = false;
    if ((*fileinfo)["GRL"] == "2015"){
      d15 = true;
      d16 = false;
    }
    if ((*fileinfo)["GRL"] == "2016"){
      d15 = false;
      d16 = true;
    }
  }
  Yggdrasil->SetBranchAddress("econtainer",&econtainer);
  Yggdrasil->SetBranchAddress("mcontainer",&mcontainer);
  Yggdrasil->SetBranchAddress("tcontainer",&tcontainer);
  Yggdrasil->SetBranchAddress("jcontainer",&jcontainer);
  Yggdrasil->SetBranchAddress("j12container",&j12container);
  Yggdrasil->SetBranchAddress("j08container",&j08container);
  Yggdrasil->SetBranchAddress("metcontainer",&metcontainer);
  return;
}

void PreAnalysis::SetName(string fname){
  size_t pos = fname.rfind("/");
  if (pos!=std::string::npos){
    filename=fname.substr(pos+1);
  }
  else{
    filename = fname;
  }
  size_t pos1 = filename.find("directTT");
  if (pos1!= std::string::npos){ 
    isSignal = true;
    size_t pos2 = filename.find("_",pos1+1);
    size_t pos3 = filename.find("_",pos2+1);
    int32_t mass = 0;
    double mxsec = 0.;
    double muncert = 0.;
    fstream fp;
    fp.open("/home/yilin/CrossSections/xsecs.txt",ios::in);
    if (!fp){
      cout << "PreAnalysis : no file of xsec is found" << endl;
      return;
    }
    string line = "";
    while(getline(fp,line)){
      stringstream ss(line);
      ss >> mass >> mxsec >> muncert;
      if (mass == stoi(filename.substr(pos2+1,pos3-pos2-1))){
        xsec = mxsec;
      }
    }
  }
  return;
}
