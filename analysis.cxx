//#include "CRZ.h"
//#include "CRT.h"

using namespace std;

void analysis(const std::string& filetxt = "", const std::string fname = "", const std::string path = "",int64_t MaxEvent = -1, int64_t SkipEvent = 0){

  gROOT->ProcessLine(".L ~/mylib/m_loader.c+");
  //gSystem->Load("~/mylib/m_loader_c.so");
  std::cout << "starting.."<<std::endl;
  // basic declaration
  std::vector<std::string> wholefiles = {"Your time is up.You cant see me.My time is now."};
  wholefiles.clear();
  std::string dataname="";
  //std::string datapath = "/home/yilin/data/data15/";
  std::vector<std::string> filename = {"file list"};
  std::vector<std::string> fullpath={"nothing"};
  fullpath.clear();
  std::size_t foundW;
  const char *checkname = "forinitializing";
  checkname = "forinitializing";
  void *wdir = gSystem->OpenDirectory((path + fname+ "/").c_str());
  if (wdir){
    while (checkname != NULL){
      checkname = gSystem->GetDirEntry (wdir);
      if (checkname != NULL){
        std::string checknames = checkname;
        foundW = checknames.find(".root");
        if (foundW!=std::string::npos){
          filename.push_back(checknames);
          fullpath.push_back(path+fname+"/"+checknames);
        }
      }
    }
  }
  else{
      Error("No file","%s cannot be found.",(path + fname + "/").c_str());
      return;
  }
  Info("Summary!", "======================================" );
  for (const auto& file: filename){
    Info("Summary!", "%s",file.c_str() );
  }
  Info("Summary!", "======================================" );
  // create fold
  
  if (gSystem->MakeDirectory (fname.c_str()) != 0)
    std::cout<< ("could not create output directory " + fname) << std::endl;
  //if (gSystem->OpenDirectory((workinglocation + "/" + submitDir + "/").c_str())){
  //  std::cout<< ("could not open output directory " + submitDir) << std::endl;
  //  return;
  //}
  // create outputfile and tree
  TFile *output = new TFile(("./" + fname + "/" + fname +".root").c_str(),"RECREATE");
  TTree *tout = new TTree("ProcessedInfo","ProcessedInfo");
  ULong64_t proc_events = 0;
  tout->Branch("P_Events",&proc_events);
  PreAnalysis *A = new PreAnalysis;
  A->SetTFile(output);
  SRA *sra = new SRA(A);
  SRA *sra2 = new SRA(A,"TT_NEW");
  SRA *checkb = new SRA(A,"For_B");
  CRZ *crz = new CRZ(A);
  CRZ *crz2 = new CRZ(A,"OneB");
  OneLepton *crt = new OneLepton(A);
  crt->SetCR("CRT");
  ApplyCut *ac = new ApplyCut(A);
  Analysis *m_analysis = new Analysis;
  m_analysis->SetPA(A);
  m_analysis->SetAC(ac);
  bool limitevent = true;
  if(MaxEvent == -1){
    limitevent = false;
  }
  uint64_t nevents = 0;
  uint64_t m_eventcounter = 0;
  bool activate = false;
  for (auto& fp : fullpath){
    if (limitevent && m_eventcounter > MaxEvent){
       break;
    }
    TFile *file = new TFile(fp.c_str(),"READ");
    TTree *Yggdrasil = 0;
    Long64_t runevent = 0;
    if(file->FindKey("Yggdrasil")){ 
      Yggdrasil = dynamic_cast<TTree*>((file)->Get("Yggdrasil"));
      A->SetTree(Yggdrasil);
      A->SetName(fp);
      nevents = Yggdrasil->GetEntries();
    }
    else{
      Error("Yggdrasil","not found in %s",fp.c_str());
      return;
    }
    if (!activate){
      if (SkipEvent >= nevents + m_eventcounter){
        delete file;
        m_eventcounter += nevents;
        continue;
      }
      else{
        runevent = SkipEvent-m_eventcounter;
        m_eventcounter = SkipEvent;
        activate = true;
      }
    }
    for (Long64_t event = runevent; event < Yggdrasil->GetEntries() ; event++){
      if (limitevent && m_eventcounter > MaxEvent-1){
         break;
      }
      if( ( m_eventcounter % 1000) ==0 ) Info("\033[32mexecute()", "\033[0mEvent number = %lu", m_eventcounter);
      m_eventcounter++;
      Yggdrasil->GetEntry(event);
      A->init();
      m_analysis->runSRA(sra);
      //m_analysis->runBjets(checkb);
      m_analysis->runNEWSRA(sra2);
      m_analysis->runCRZ(crz);
      m_analysis->runNEWCRZ(crz2);
      m_analysis->runCRT(crt);
      proc_events++;
    }// end of a tree
    file->Close();
    delete file;
  }// end of whole files
  delete m_analysis;
  delete ac;
  delete crz;
  delete crz2;
  delete sra;
  delete sra2;
  delete crt;
  delete checkb;
  delete A;
  tout->Fill();
  output->Write();
  output->Close();
  delete output;
  std::cout << "The fold created automatically is \"\033[31m" << fname<< "\033[0m\"."<<std::endl;
  return;
}

