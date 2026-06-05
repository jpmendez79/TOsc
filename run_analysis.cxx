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
#include <chrono> // timer

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
  const int NUM_dm2 = 60;
  const int NUM_ttt = 60;

  // TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -2, 1);
  // TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -2, 0);
  // // First load 3600 files

  // Pre-create all Vectors
  TH1::AddDirectory(false);
  std::vector<TH2D *> vec_cls_universe;
  for (int k = 0; k < 2002; ++k) {
    vec_cls_universe.push_back(
        new TH2D(Form("hcls_%04d", k), "", 60, 1, 60, 60, 1, 60));
  }

  for (int idm2 = 1; idm2 <= 60; idm2++) {
    for (int ittt = 1; ittt <= 60; ittt++) {

      // Open the file
      TString roostr = TString::Format("output/cls_out_dm2_ttt_%03d_%03d.root", idm2, ittt);
      TFile f(roostr, "READ");

      // Get the tree
      TTree *tree = nullptr;
      f.GetObject("tree", tree);

      // Disable all unused branches
      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus("vec_confidence", 1);

      // Connect the branch
      std::vector<double> *vec_confidence = nullptr;
      tree->SetBranchAddress("vec_confidence", &vec_confidence);
      tree->GetEntry(0);

      // Fill histograms
      int num_universe = vec_confidence->size();
      for (int universe = 0; universe < num_universe; universe++) {
        vec_cls_universe[universe]->SetBinContent(ittt, idm2, (*vec_confidence)[universe]);
      }
      f.Close();
    }
  }
  // Now Do the graph interpretation
  std::vector<TGraph*> cl_curves(vec_cls_universe.size());

  for (int u = 0; u < vec_cls_universe.size(); u++) {

  cl_curves[u] = new TGraph();

  get_CL_curve(vec_cls_universe[u], cl_curves[u],u);
  }

  // Now save everything
  TFile out("cls_95_60x60_vanilla_numu.root", "RECREATE");

  TDirectory *dh = out.mkdir("histograms");
  TDirectory *dg = out.mkdir("cl_curves");
  // histograms
dh->cd();
for (int u = 0; u < vec_cls_universe.size(); u++) {
  vec_cls_universe[u]->Write(Form("h_%04d", u));
}

// graphs
dg->cd();
for (int u = 0; u < vec_cls_universe.size(); u++) {
  cl_curves[u]->Write(Form("g_%04d", u));
}
out.Close();
}
