#include "stdlib.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#include <map>
#include <set>
#include <vector>

#include "WCPLEEANA/TOsc.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"
#include "TParameter.h"
#include <algorithm>
#include <chrono> // timer

/// Obtain the TGraph(s) corresponding to a particular contour level for a TH2
///
/// \param h2     The TH2 to examine
/// \param level  Fractional level in [0, 1]
/// \return       Vector of TGraphs that represent the whole contour (multiple
/// if contour has disjoint pieces)
std::vector<TGraph *> GetContourGraphs(TH2 *h2, double level) {
  std::vector<TGraph *> ret;

  std::unique_ptr<TH2> surf(
      dynamic_cast<TH2 *>(h2->Clone("tmp_h2_for_drawing_graphs")));

  TVirtualPad *bak = gPad;

  const bool wasbatch = gROOT->IsBatch();
  gROOT->SetBatch(); // User doesn't want to see our temporary canvas
  TCanvas tmp;

  gStyle->SetOptStat(0);

  surf->SetContour(1, &level);
  surf->Draw("cont list");

  tmp.Update();
  tmp.Paint();

  gROOT->SetBatch(wasbatch);
  gPad = bak;

  // The graphs we need (contained inside TLists, contained inside
  // TObjArrays) are in the list of specials. But we need to be careful about
  // types, because other stuff can get in here too (TDatabasePDG for
  // example).
  TCollection *specs = gROOT->GetListOfSpecials();

  TIter nextSpec(specs);
  while (TObject *spec = nextSpec()) {
    if (!spec->InheritsFrom(TObjArray::Class()))
      continue;
    auto conts = dynamic_cast<TObjArray *>(spec);

    if (conts->IsEmpty())
      continue;

    if (!conts->At(0)->InheritsFrom(TList::Class()))
      continue;
    auto cont = dynamic_cast<TList *>(conts->At(0));

    TIter nextObj(cont);
    // Contour could be split into multiple pieces
    std::size_t piece = 0;
    while (TObject *obj = nextObj()) {
      if (!obj->InheritsFrom(TGraph::Class()))
        continue;

      ret.push_back(dynamic_cast<TGraph *>(obj->Clone(
          Form("%s_contour%f_piece%zu", obj->GetName(), level, piece))));
      piece++;
    } // end for obj
  } // end for spec


  return ret;
}
void get_CL_curve(TH2 *h2_CL_input, TGraph *gh_CL_curve, int flag_index)
{
  gROOT->SetBatch( 1 );

  TString roostr = "";

  roostr = TString::Format("canv_list_95_%04d", flag_index);
  TCanvas *canv_list_95 = new TCanvas(roostr, roostr, 800, 600);
  cout<<endl<<" ---> "<<roostr<<endl;

  double user_contours_95[1] = {0.95};
  //double user_contours_95[1] = {0.683};
  h2_CL_input->SetContour(1, user_contours_95);

  // Draw contours as filled regions, and Save points
  h2_CL_input->Draw("CONT Z LIST");
  canv_list_95->Update();

  // Get Contours
  TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");

  if (!conts){
    printf("*** No Contours Were Extracted!\n");
  }

  TList* contLevel = nullptr;
  TGraph* curv     = nullptr;

  Int_t TotalConts = conts->GetSize();

  printf("TotalConts = %d\n", TotalConts);

  for(int i = 0; i < TotalConts; i++){
    contLevel = (TList*)conts->At(i);
    printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
  }

  ///////

  int j_maxsize_index = 0;
  int j_maxsize = -1;

  int num_graphs = contLevel->GetSize();

  for(int i = 0; i < TotalConts; i++) {
    contLevel = (TList*)conts->At(i);

    if(i!=0) continue;// only use the one having the most points

    for(int j = 0; j < contLevel->GetSize(); j++) {
      curv = (TGraph*)contLevel->At(j);
      int Npoints = curv->GetN();
      printf(" ---> graph %2d, size %3d\n", j+1, Npoints);

      if( Npoints >=j_maxsize  ) {
	j_maxsize = Npoints;
	j_maxsize_index = j;
      }

    }// for(int j = 0; j < contLevel->GetSize(); j++)
  }// for(int i = 0; i < TotalConts; i++)

  ///////

  for(int i = 0; i < TotalConts; i++) {
    contLevel = (TList*)conts->At(i);

    if(i!=0) continue;// only use the one having the most points

    for(int j = 0; j < contLevel->GetSize(); j++) {
      curv = (TGraph*)contLevel->At(j);
      int Npoints = curv->GetN();
      printf(" ---> graph %2d, size %3d\n", j+1, Npoints);


      if( num_graphs==2 ) {
	for(int ip=0; ip<Npoints; ip++) {
	  double dm2_val, ttt_val;
	  curv->GetPoint(ip, ttt_val, dm2_val);

	  //gh_CL_curve->SetPoint( gh_CL_curve->GetN(), pow(10, ttt_val), pow(10, dm2_val) );

	  gh_CL_curve->SetPoint( gh_CL_curve->GetN(), ttt_val, dm2_val );
	}// for(int ip=0; ip<Npoints; ip++)
      }
      else {

	if(j==j_maxsize_index) {
	  for(int ip=0; ip<Npoints; ip++) {
	    double dm2_val, ttt_val;
	    curv->GetPoint(ip, ttt_val, dm2_val);

	    //gh_CL_curve->SetPoint( gh_CL_curve->GetN(), pow(10, ttt_val), pow(10, dm2_val) );

	    gh_CL_curve->SetPoint( gh_CL_curve->GetN(), ttt_val, dm2_val );
	  }// for(int ip=0; ip<Npoints; ip++)
	}// if(j==j_maxsize_index)

      }


    }// for(int j = 0; j < contLevel->GetSize(); j++)
  }// for(int i = 0; i < TotalConts; i++)

  ///////

  // for(int i = 0; i < TotalConts; i++) {
  //   contLevel = (TList*)conts->At(i);

  //   if(i!=0) continue;// only use the one having the most points

  //   for(int j = 0; j < contLevel->GetSize(); j++) {
  //     curv = (TGraph*)contLevel->At(j);
  //     int Npoints = curv->GetN();
  //     printf(" ---> graph %2d, size %3d\n", j+1, Npoints);

  //     if(i==0) {
  // 	for(int ip=0; ip<Npoints; ip++) {
  // 	  double dm2_val, ttt_val;
  // 	  curv->GetPoint(ip, ttt_val, dm2_val);

  // 	  //gh_CL_curve->SetPoint( gh_CL_curve->GetN(), pow(10, ttt_val), pow(10, dm2_val) );

  // 	  gh_CL_curve->SetPoint( gh_CL_curve->GetN(), ttt_val, dm2_val );

  // 	}// for(int ip=0; ip<Npoints; ip++)
  //     }

  //   }// for(int j = 0; j < contLevel->GetSize(); j++)
  // }// for(int i = 0; i < TotalConts; i++)

  // cout<<endl;

  gROOT->SetBatch( 0 );
}
int main(void) {
  int N_idm2 = 60;
  int N_it14 = 60;
  double contour_level = 0.95;

  std::vector<std::vector<double>> cls_cache;
  // Number of files
  cls_cache.resize(3600);
  for (int idm2 = 0; idm2 < 60; ++idm2) {
    for (int it14 = 0; it14 < 60; ++it14) {

      int idx = idm2 * N_it14 + it14;

      TString roostr = TString::Format("output/size_1k_cls_out_dm2_ttt_%03d_%03d.root", idm2+1, it14+1);
      TFile file(roostr, "READ");

      TTree* tree = (TTree*)file.Get("tree");
      std::vector<double>* invecconf = nullptr;
      tree->SetBranchAddress("vec_confidence", &invecconf);
      tree->GetEntry(0);

      cls_cache[idx] = *invecconf;
      // file closes automatically here (scope end)
    }
  }
  // for (auto& v : vec4v_cache) std::sort(v.begin(), v.end());
  // for (auto& v : vec3v_cache) std::sort(v.begin(), v.end());
  cout << "Finish Presave all out*.root files\n";
  const int N_universe = 2002;
  for (int universe = 0; universe < N_universe; universe++) {
    cout <<"START UNIVERSE " << universe << endl;
    // Create a Histogram
    TH2D histConfidence("histConfidence", "Confidence Grid", 60,1,60,60,1,60);
    // Loop through and fill out the histogram
    for (int idm2 = 0; idm2 < 60; ++idm2) {
      for (int it14 = 0; it14 < 60; ++it14) {
        int idx = idm2 * N_it14 + it14;
        double confidence = cls_cache[idx][universe];
        // cout << cls_cache[idx].size() << endl;
        histConfidence.SetBinContent(it14+1, idm2+1, confidence);
      }
    }
    // Run the Contour Finder and report
    TGraph *gh = new TGraph();
    get_CL_curve(&histConfidence, gh, universe);
  }


}
