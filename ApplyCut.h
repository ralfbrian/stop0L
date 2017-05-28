#ifndef _PREANALYSIS_H_
#include "PreAnalysis.h"
#endif
using namespace std;
class ApplyCut{
  public:
    ApplyCut(PreAnalysis *Ana);
    void AddCut(vector<bool> cut);
    template <class X>
    void FillHist(TH1 *h[4],X data, vector<double> *weight);

  private:
    vector<bool> cutflow;
    PreAnalysis *A;
};

ApplyCut::ApplyCut(PreAnalysis *Ana):A(Ana){
}

void ApplyCut::AddCut(vector<bool> cut){
  cutflow.clear();
  cutflow.swap(cut);
  return;
}

template <class X>
void ApplyCut::FillHist(TH1 *h[4],X data, vector<double> *weight){
  for (size_t i = 0; i < cutflow.size(); i++){
    if (!cutflow[i])continue;
    if (A->isMC){
      h[i]->Fill(data,(*weight)[i]);
    }
    else{
      h[i]->Fill(data);
    }
  }
  return;
}
