#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <stdlib.h>
#include <stdio.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include "HistInit.h"
#include "datadefinition.h"
#include "stringmanage_root.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "MT2.h"
#include "MT2_ROOT.h"
#include "PreAnalysis.h"
#include "SRA.h"
#include "CRZ.h"
#include "OneLepton.h"
#include "ApplyCut.h"
#include "Analysis.h"


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
  dRbb = -100.;
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
  calmt2 = false;
  MT2 = 0.;
  return;
}

bool compjet(objPro *a, objPro *b){return a->yylpt > b->yylpt;}
bool compbjet(objPro *a, objPro *b){return a->yylbtagweight >= b->yylbtagweight;}

void PreAnalysis::Bsort(vector<objPro*>& v_bjets){
  if (v_bjets.size()<=1) return;
  std::vector<objPro*>::iterator itr = v_bjets.begin();
  std::vector<objPro*>::iterator p_itr;
  while(itr !=v_bjets.end()){
    if (itr == v_bjets.begin()) {
      p_itr=itr;
      itr++;
    }
    else{
      p_itr = itr-1;
    }
    if (!compbjet((*p_itr),(*itr))){
      objPro *p_ptr = *p_itr;
      objPro *ptr = *itr;
      (*itr) = p_ptr;
      (*p_itr) = ptr;
      itr--;
    }
    else{
      itr++;
    }
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
  fabsdphi = std::numeric_limits<float>::max();
  for (auto& jet : *jcontainer){
    if (jet.passOR != 0   && jet.baseline == 1){
      jet_idx[PreAnalysis::baseline]++;
      if (jet.bad == 1 ){
        jet_idx[PreAnalysis::bad]++;
      }
      else{
        if (jet.signal == 0)continue;
        jet_idx[PreAnalysis::signal]++;
        jetptsort.push_back(jet.yylpt*0.001);
        getbasejet.push_back(&jet);
        if (jet.yylpt > 20000.) jet20++;
        if (jet.yylpt > 40000.) jet40++;
        if (jet.yylpt > 80000.) jet80++;
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
        if (jet.btagged == 1 ){
          jet_idx[PreAnalysis::btagged]++;
          if (jet.yyleta < 2.5){
            getbjet.push_back(&jet);
            float phitemp = fabs(TVector2::Phi_mpi_pi(jet.yylphi-tst.Phi()));
            if (phitemp < fabsdphi){
              Mtbmin = sqrt(2*(jet.yylpt)*(tst.Pt())*(1-cos(jet.yylphi-tst.Phi())));
              fabsdphi = phitemp ;
            }
          }
        }
      }
    }
  }
  std::sort(jetptsort.rbegin(),jetptsort.rend());
  std::sort(getbasejet.begin(),getbasejet.end(),compjet);
  //std::sort(getbjet.begin(),getbjet.end(),compjet);
  Bsort(getbjet);
  if (getbjet.size() >= 2){
    if (getbjet[0]->yylbtagweight < getbjet[1]->yylbtagweight)
      cout << getbjet[0]->yylbtagweight << " " << getbjet[1]->yylbtagweight<<endl;
    dRbb = sqrt(pow(getbjet[0]->yyleta-getbjet[1]->yyleta,2)+pow(TVector2::Phi_mpi_pi(getbjet[0]->yylphi-getbjet[1]->yylphi),2));
  }
  if (getbasejet.size()>=1){
    fjet = getbasejet[0];
    firstjet.SetPtEtaPhiM(fjet->yylpt,fjet->yyleta,fjet->yylphi,fjet->yylm);
  }
  if (getbasejet.size()>=2){
    sjet = getbasejet[1];
    secondjet.SetPtEtaPhiM(sjet->yylpt,sjet->yyleta,sjet->yylphi,sjet->yylm);
  }
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
        ele_idx[PreAnalysis::signal]++;
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
      if (muon.signal == 1 )
        mu_idx[PreAnalysis::signal]++;
    }
  }
  return;
}

void PreAnalysis::j12cal(){
  j12sort.clear();
  std::vector<objPro> temp;
  RClustering(getbasejet,temp,1.2,5000.);
  for (auto& jet12 : temp){
    j12sort.push_back(jet12.yylm);
  }
  //for (auto& jet12 : *j12container){
  //  j12sort.push_back(jet12.yylm);
  //}
  //std::sort(j12sort.rbegin(),j12sort.rend());
  return;
}

void PreAnalysis::j08cal(){
  j08sort.clear();
  std::vector<objPro> temp;
  RClustering(getbasejet,temp,0.8,5000.);
  for (auto& jet08 : temp){
    j08sort.push_back(jet08.yylm);
  }
  //std::sort(j08sort.rbegin(),j08sort.rend());
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
    float phitemp = fabs(TVector2::Phi_mpi_pi(jet->yylphi-tst.Phi()));
    if (phitemp < minphi &&jet->btagged == 0 && jet->yylnTrk <= 4 && fabs(jet->yyleta) < 2.5){
      taujet = jet;
      minphi = phitemp;
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
  size_t pos4 = filename.find(".");
  size_t pos5 = filename.find(".",pos4+1);
  runnumber = stoi(filename.substr(pos4+1,pos5-pos4-1));
  return;
}

void PreAnalysis::FillLz(TLorentzVector& vec, objPro *obj){
  vec.SetPtEtaPhiM(obj->yylpt,obj->yyleta,obj->yylphi,obj->yylm);
  return;
}

void PreAnalysis::RClustering(std::vector<objPro*>& obj,std::vector<objPro>& rcjet,float dr, float ptmin){
		std::vector<fastjet::PseudoJet> inputPseudoJets;
		std::vector<fastjet::PseudoJet> rcPseudoJets;
		std::vector<fastjet::PseudoJet> constPseudoJets;
		for(size_t j=0;j<obj.size();j++){
      TLorentzVector tmpjet;
      tmpjet.SetPtEtaPhiM(obj[j]->yylpt,obj[j]->yyleta,obj[j]->yylphi,obj[j]->yylm);
			fastjet::PseudoJet pseudoJet(tmpjet.Px(), tmpjet.Py(), tmpjet.Pz(), tmpjet.E());
			pseudoJet.set_user_index(j);
			inputPseudoJets.push_back(pseudoJet);
		}
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, dr);

		fastjet::ClusterSequence* currentClusSeq = new fastjet::ClusterSequence(inputPseudoJets, jetDef);

		rcPseudoJets = fastjet::sorted_by_pt(currentClusSeq->inclusive_jets(ptmin));
		for (unsigned int rcPJetsIndex = 0; rcPJetsIndex < rcPseudoJets.size(); rcPJetsIndex++){
      TLorentzVector tmpjet;
      double px = rcPseudoJets[rcPJetsIndex].px();
      double py = rcPseudoJets[rcPJetsIndex].py();
      double pz = rcPseudoJets[rcPJetsIndex].pz();
      double e = rcPseudoJets[rcPJetsIndex].e();
      tmpjet.SetPxPyPzE(px,py,pz,e);
      objPro opro;
      opro.yylpt = tmpjet.Pt(); 
      opro.yyleta = tmpjet.Eta();
      opro.yylphi = tmpjet.Phi();
      opro.yyle = tmpjet.E();
      opro.yylm = tmpjet.M();
      opro.yylrapidity= tmpjet.Rapidity();
      rcjet.push_back(opro);
    }
		delete currentClusSeq;
 
  return;
}

bool PreAnalysis::chi2Top(){
		if(getbasejet.size()< 4) return false;
		if(getbjet.size()< 2) return false;

		float realWMass = 80385;
		float realTopMass = 173210;
		
		//Chi2 method
		Chi2min = DBL_MAX;
		int W1j1_low = -1,W1j2_low = -1,W2j1_low = -1,W2j2_low = -1,b1_low = -1,b2_low = -1;

		// First separate light from b-jets
		std::vector<TLorentzVector> bjets(2,TLorentzVector());  //fill with TLorentzVectors of the two signal jets with the highest b-tag
		std::vector<TLorentzVector> ljets;  //fill with TLorentzVectors of ALL OTHER signal jets (even if they are b-tagged)

		for(size_t j=0;j<getbasejet.size();j++){
			TLorentzVector tempTLV;
      FillLz(tempTLV,getbasejet[j]);
			if(getbasejet[j]==getbjet[0]){
				bjets[0] = tempTLV;
			}
			else if(getbasejet[j]==getbjet[1]){
				bjets[1] = tempTLV;
			}
			else {
				ljets.push_back(tempTLV);
			}
		}
		for(size_t W1j1=0; W1j1<ljets.size(); W1j1++) {// <------------------This lines has to be replaced
			for(size_t W2j1=0;W2j1<ljets.size(); W2j1++) {// <------------------This lines has to be replaced
				if (W2j1==W1j1) continue;// <------------------This lines has to be added
				//		for(int W1j1=0; W1j1<(int)ljets.size()-1; W1j1++) {
				//		for(int W2j1=W1j1+1;W2j1<(int)ljets.size(); W2j1++) {
				for(size_t b1=0;b1<bjets.size();b1++){
					for(size_t b2=0;b2<bjets.size();b2++){
						if(b2==b1) continue;
						double chi21, chi22, mW1, mW2, mt1, mt2;
						//try bl,bl top candidates
						if(W2j1>W1j1){
							mW1 = ljets.at(W1j1).M();
							mW2 = ljets.at(W2j1).M();
							mt1 = (ljets.at(W1j1) + bjets.at(b1)).M();
							mt2 = (ljets.at(W2j1) + bjets.at(b2)).M();
							chi21 = (mW1-realWMass)*(mW1-realWMass)/realWMass + (mt1-realTopMass)*(mt1-realTopMass)/realTopMass;
							chi22 = (mW2-realWMass)*(mW2-realWMass)/realWMass + (mt2-realTopMass)*(mt2-realTopMass)/realTopMass;
							if(Chi2min > (chi21 + chi22)){
								Chi2min = chi21 + chi22;
								if(chi21 < chi22){
									W1j1_low = W1j1;
									W1j2_low = -1;
									W2j1_low = W2j1;
									W2j2_low = -1;
									b1_low = b1;
									b2_low = b2;
								}
								else{
									W2j1_low = W1j1;
									W2j2_low = -1;
									W1j1_low = W2j1;
									W1j2_low = -1;
									b2_low = b1;
									b1_low = b2;
								}
							}
						}
						if(ljets.size() < 3) continue;
						for(size_t W1j2=W1j1+1;W1j2<ljets.size(); W1j2++) {
							if(W1j2==W2j1) continue;
							
							//try bll,bl top candidates
							mW1 = (ljets.at(W1j1) + ljets.at(W1j2)).M();
							mW2 = ljets.at(W2j1).M();
							mt1 = (ljets.at(W1j1) + ljets.at(W1j2) + bjets.at(b1)).M();
							mt2 = (ljets.at(W2j1) + bjets.at(b2)).M();
							chi21 = (mW1-realWMass)*(mW1-realWMass)/realWMass + (mt1-realTopMass)*(mt1-realTopMass)/realTopMass;
							chi22 = (mW2-realWMass)*(mW2-realWMass)/realWMass + (mt2-realTopMass)*(mt2-realTopMass)/realTopMass;
							if(Chi2min > (chi21 + chi22)){
								Chi2min = chi21 + chi22;
								if(chi21 < chi22){
									W1j1_low = W1j1;
									W1j2_low = W1j2;
									W2j1_low = W2j1;
									W2j2_low = -1;
									b1_low = b1;
									b2_low = b2;
								}
								else{
									W2j1_low = W1j1;
									W2j2_low = W1j2;
									W1j1_low = W2j1;
									W1j2_low = -1;
									b2_low = b1;
									b1_low = b2;
								}
							}
							if(ljets.size() < 4)continue;
							//try bll, bll top candidates
							for(size_t W2j2=W2j1+1;W2j2<ljets.size(); W2j2++){
								if((W2j2==W1j1) || (W2j2==W1j2)) continue;
								if(W2j1<W1j1) continue;  //runtime reasons, we don't want combinations checked twice <--------------------This line should be added
								mW1 = (ljets.at(W1j1) + ljets.at(W1j2)).M();
								mW2 = (ljets.at(W2j1) + ljets.at(W2j2)).M();
								mt1 = (ljets.at(W1j1) + ljets.at(W1j2) + bjets.at(b1)).M();
								mt2 = (ljets.at(W2j1) + ljets.at(W2j2) + bjets.at(b2)).M();
								chi21 = (mW1-realWMass)*(mW1-realWMass)/realWMass + (mt1-realTopMass)*(mt1-realTopMass)/realTopMass;
								chi22 = (mW2-realWMass)*(mW2-realWMass)/realWMass + (mt2-realTopMass)*(mt2-realTopMass)/realTopMass;
								if(Chi2min > (chi21 + chi22)){
									Chi2min = chi21 + chi22;
									if(chi21 < chi22){
										W1j1_low = W1j1;
										W1j2_low = W1j2;
										W2j1_low = W2j1;
										W2j2_low = W2j2;
										b1_low = b1;
										b2_low = b2;
									}
									else{
										W2j1_low = W1j1;
										W2j2_low = W1j2;
										W1j1_low = W2j1;
										W1j2_low = W2j2;
										b2_low = b1;
										b1_low = b2;
									}
								}
							}
						}
					}
				}
			}
		}
		if(W1j2_low == -1) WCand0 = ljets.at(W1j1_low);
		else WCand0 = ljets.at(W1j1_low) + ljets.at(W1j2_low);
		topCand0 = WCand0 + bjets.at(b1_low);
		if(W2j2_low == -1) WCand1 = ljets.at(W2j1_low);
		else WCand1 = ljets.at(W2j1_low) + ljets.at(W2j2_low);
		topCand1 = WCand1 + bjets.at(b2_low);
		return true;
}

void PreAnalysis::CalMT2(){
  if (calmt2) return;
  if (!chi2Top()){
    calmt2 = true;
    return;
  }
  TLorentzVector mtop0;
  TLorentzVector mtop1;
  mtop0.SetPtEtaPhiM(topCand0.Pt(),0,topCand0.Phi(),173210);
  mtop1.SetPtEtaPhiM(topCand1.Pt(),0,topCand1.Phi(),173210);
  MT2 = ComputeMT2(mtop0,mtop1,tst).Compute();
  calmt2 = true;
  return;
}
