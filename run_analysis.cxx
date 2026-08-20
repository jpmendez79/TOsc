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
  const double DM2_LO = -1, DM2_HI = 2;
  const double TTT_LO = -3, TTT_HI = 0;
  double xbins[61], ybins[61];
  // Construct log10 bin boundaries
  for (int i = 0; i <= 60; i++) {
    for (int j = 0; j <= 60; j++) {
      xbins[i] = pow(10, TTT_LO + i*(TTT_HI-TTT_LO)/NUM_ttt);  // theta bins
      ybins[j] = pow(10, DM2_LO + j*(DM2_HI-DM2_LO)/NUM_dm2);  // dm2 bins
    }
  }
  //
  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, &ybins[0]);
  TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, &xbins[0]);

   // TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -1, 2);
   // TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -3, 0);

  // First load 3600 files

  // Pre-create all Vectors
  TH1::AddDirectory(false);
  std::vector<TH2D *> vec_cls_universe;

  for (int k = 0; k < 2002; ++k) {
    // vec_cls_universe.push_back(
        // new TH2D(Form("hcls_%04d", k), "", 60, TTT_LO, TTT_HI, 60, DM2_LO, DM2_HI));
    vec_cls_universe.push_back(   new TH2D(Form("hcls_%04d", k), "", NUM_ttt, &xbins[0], NUM_dm2, &ybins[0]));//
  }
  std::vector<TGraph *> cl_curves(vec_cls_universe.size());
  std::vector<TGraph*> cl_curves_invert(vec_cls_universe.size());
  for (int idm2 = 1; idm2 <= 60; idm2++) {
    for (int ittt = 1; ittt <= 60; ittt++) {

      // TString roostr = TString::Format("output/out_dm2_ttt_%03d_%03d.root", idm2, ittt);
      TString roostr = TString::Format("output/BNBvanilla_numu_disp_grid_60x60_dm2_ttt_%03d_%03d.root", idm2, ittt);
      TFile f(roostr, "READ");

      // Get the tree
      TTree *tree = nullptr;
      f.GetObject("tree", tree);

      // Disable all unused branches
      tree->SetBranchStatus("*", 0);
      tree->SetBranchStatus("vec_confidence", 1);
    // tree->SetBranchStatus("vec_dchi2_with_data", 1);
      // Connect the branch
      std::vector<double> *vec_confidence = nullptr;
      tree->SetBranchAddress("vec_confidence", &vec_confidence);
      // tree->SetBranchAddress("vec_dchi2_with_data", &vec_confidence);
tree->GetEntry(0);

      // Fill histograms
      // int num_universe = vec_confidence->size();
      int num_universe = 2002;
      for (int universe = 0; universe < num_universe; universe++) {
        vec_cls_universe[universe]->SetBinContent(ittt, idm2,
                                                  (*vec_confidence)[universe]);
      }
      f.Close();
    }
  }

  // Now save everything
  TFile out("cls_map_BNB_vanilla_disp_60x60-notfancy.root", "RECREATE");

  // histograms

for (int u = 0; u < vec_cls_universe.size(); u++) {
  vec_cls_universe[u]->Write(Form("h2_%04d", u));
  cl_curves[u] = new TGraph();
  cl_curves_invert[u] = new TGraph();
  get_CL_curve(vec_cls_universe[u], cl_curves[u], u);
  int size = cl_curves[u]->GetN();
  for(int idx=0; idx<size; idx++) {
    double xx, yy;
    cl_curves[u]->GetPoint(idx, xx, yy);
    cl_curves_invert[u]->SetPoint(idx, yy, xx);
    }

  cl_curves[u]->Write(Form("gr_%04d", u));
  cl_curves_invert[u]->Write(Form("iv_%04d", u));
}

// l2sigma, l1sigma, m, u1sigma, u2sigma
double p[5] = {0.023, 0.159, 0.5, .841, .977};
double l2sigma[NUM_dm2];
double l1sigma[NUM_dm2];
double median[NUM_dm2];
double u1sigma[NUM_dm2];
double u2sigma[NUM_dm2];
for (int idm2 = 0; idm2 < NUM_dm2; idm2++) {
  std::vector<double> xvals;
  for (int universe = 0; universe < 2002; universe++) {
    xvals.push_back(cl_curves_invert[universe]->Eval(ybins[idm2]));
  }
  std::sort(xvals.begin(), xvals.end());

  int n = xvals.size();
  auto percentile_val = [&](double frac) {
    int idx = std::min(std::max((int)(frac * n), 0), n - 1);
    return xvals[idx];
  };

  double v_l2 = percentile_val(p[0]);
  double v_l1 = percentile_val(p[1]);
  double v_m  = percentile_val(p[2]);
  double v_u1 = percentile_val(p[3]);
  double v_u2 = percentile_val(p[4]);

  median[idm2]  = v_m;
  l1sigma[idm2] = v_m  - v_l1;
  u1sigma[idm2] = v_u1 - v_m;
  l2sigma[idm2] = v_m  - v_l2;
  u2sigma[idm2] = v_u2 - v_m;
}
TGraphAsymmErrors *sigma1 = new TGraphAsymmErrors(NUM_dm2, &median[0], &ybins[0], &l1sigma[0], &u1sigma[0]);
TGraphAsymmErrors *sigma2 = new TGraphAsymmErrors(NUM_dm2, &median[0], &ybins[0], &l2sigma[0], &u2sigma[0]);
sigma1->Write("sigma1");
sigma2->Write("sigma2");
// graphs
// dg->cd();
// for (int u = 0; u < vec_cls_universe.size(); u++) {
//
// }
out.Close();
}
