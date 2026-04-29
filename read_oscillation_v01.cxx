#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<vector>
#include<map>
#include<set>

#include "WCPLEEANA/TOsc.h"

#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"
#include "TParameter.h"
#include <chrono> // timer

////////////////////////////////// SUB /////////////////////////////////////////////

std::vector<TGraph*> GetContourGraphs(TH2* h2, double level);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A Hello story ..."<<endl<<endl;
  int analysis = 0;
  int ifile = -1;
  int obs_throw = 0;
  int fcls = 0;
  int xthrow = 0;
  int draw_confidence_map = 1;
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;


  int sensitivity_study = 0;
  int it14 = 0;
  int idm2 = 0;
  int it24 = 0;

  const int NUM_dm2 = 60;
  const int NUM_ttt = 60;
  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -2, 1);
  TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -2, 0);



  // double it14 = 0;
  // double idm2 = 0;
  // double it24 = 0;

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-ifile")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-fobsthrow")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>obs_throw ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-fcls")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>fcls ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
    if (strcmp(argv[i], "-option") == 0) {
      if (strcmp(argv[i+1], "analysis") == 0) {
        analysis = 1;
        // cout << "ANALYSIS MODE \n";
                }
    }
    if (strcmp(argv[i], "-flag") == 0) {
      if (strcmp(argv[i+1], "sens") == 0) {
        sensitivity_study = 1;
        // cout << "ANALYSIS MODE \n";
      }
      if (strcmp(argv[i+1], "xfile") == 0) {
        xthrow = 1;
        // cout << "ANALYSIS MODE \n";
      }
    }
    if( strcmp(argv[i],"-pbnb")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_BNB ) ) { cerr<<" ---> Error scaleF_POT_BNB !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-pnumi")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_NuMI ) ) { cerr<<" ---> Error scaleF_POT_NuMI !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-d")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>display ) ) { cerr<<" ---> Error display !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it14")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it14 ) ) { cerr<<" ---> Error it14 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-idm2")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>idm2 ) ) { cerr<<" ---> Error idm2 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it24")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it24 ) ) { cerr<<" ---> Error it24 !"<<endl; exit(1); }
    }
  }

  ///////////////////////////////////////////////////////////

  if( !display ) {
    gROOT->SetBatch( 1 );
  }

  TApplication theApp("theApp",&argc,argv);

  /////////////////////////////////////////////////////////// Draw style

  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ///////////////////////////////////////////////////////////

  TOsc *osc_test = new TOsc();

  ///////

  osc_test->tosc_scaleF_POT_BNB  = scaleF_POT_BNB;
  osc_test->tosc_scaleF_POT_NuMI = scaleF_POT_NuMI;

  ///////

  osc_test->flag_apply_oscillation_BNB  = Configure_Osc::flag_apply_oscillation_BNB;
  osc_test->flag_apply_oscillation_NuMI = Configure_Osc::flag_apply_oscillation_NuMI;

  osc_test->flag_goodness_of_fit_CNP    = Configure_Osc::flag_goodness_of_fit_CNP;

  ///////

  osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;

  ///////

  osc_test->flag_NuMI_nueCC_from_intnue         = Configure_Osc::flag_NuMI_nueCC_from_intnue;
  osc_test->flag_NuMI_nueCC_from_overlaynumu    = Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
  osc_test->flag_NuMI_nueCC_from_appnue         = Configure_Osc::flag_NuMI_nueCC_from_appnue;
  osc_test->flag_NuMI_nueCC_from_appnumu        = Configure_Osc::flag_NuMI_nueCC_from_appnumu;
  osc_test->flag_NuMI_nueCC_from_overlaynueNC   = Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
  osc_test->flag_NuMI_nueCC_from_overlaynumuNC  = Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;

  osc_test->flag_NuMI_numuCC_from_overlaynumu   = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue    = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue        = Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu       = Configure_Osc::flag_NuMI_numuCC_from_appnumu;
  osc_test->flag_NuMI_numuCC_from_overlaynumuNC = Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
  osc_test->flag_NuMI_numuCC_from_overlaynueNC  = Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;

  osc_test->flag_NuMI_CCpi0_from_overlaynumu    = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue         = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  osc_test->flag_NuMI_CCpi0_from_overlaynumuNC  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_CCpi0_from_overlaynueNC   = Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;

  osc_test->flag_NuMI_NCpi0_from_overlaynumu    = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue         = Configure_Osc::flag_NuMI_NCpi0_from_appnue;
  osc_test->flag_NuMI_NCpi0_from_overlaynumuNC  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_NCpi0_from_overlaynueNC   = Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;


  ///////

  osc_test->flag_BNB_nueCC_from_intnue         = Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu    = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue         = Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu        = Configure_Osc::flag_BNB_nueCC_from_appnumu;
  osc_test->flag_BNB_nueCC_from_overlaynueNC   = Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
  osc_test->flag_BNB_nueCC_from_overlaynumuNC  = Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;

  osc_test->flag_BNB_numuCC_from_overlaynumu   = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue    = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue        = Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu       = Configure_Osc::flag_BNB_numuCC_from_appnumu;
  osc_test->flag_BNB_numuCC_from_overlaynumuNC = Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
  osc_test->flag_BNB_numuCC_from_overlaynueNC  = Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;

  osc_test->flag_BNB_CCpi0_from_overlaynumu    = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue         = Configure_Osc::flag_BNB_CCpi0_from_appnue;
  osc_test->flag_BNB_CCpi0_from_overlaynumuNC  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_CCpi0_from_overlaynueNC   = Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;

  osc_test->flag_BNB_NCpi0_from_overlaynumu    = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue         = Configure_Osc::flag_BNB_NCpi0_from_appnue;
  osc_test->flag_BNB_NCpi0_from_overlaynumuNC  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_NCpi0_from_overlaynueNC   = Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;

  /////// set only one time
  /////// set only one time

  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
                               Configure_Osc::default_dirtadd_file,
                               Configure_Osc::default_mcstat_file,
                               Configure_Osc::default_fluxXs_dir,
                               Configure_Osc::default_detector_dir);

  osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);

  /////////////////////////////////////////////////////////// choose to use the dataset of BNB+NuMI, BNB, or NuMI??
  //
  // A trick way: do it in the function of  "TOsc::FCN(const double *par)".
  // The lines below "modify".
  //
  /////////////////////////////////////////////////////////// change the oscillation prob equations??
  //
  // Do it in the function of "TOsc::Prob_oscillaion(double Etrue, double baseline, int strflag_osc)"
  // There are examples of "dm2 vs. sin2_2Tue vs. t24", "dm2 vs. sin2_2Tee vs. t24", and "dm2 vs. t14 vs. t24".
  //
  //
  // The class-data-members used in the TOsc for the oscillation analysis are
  //
  // double tosc_dm2_41
  // double tosc_si2n_2theta_14
  // double tosc_sin2_theta_24
  // double tosc_sin2_theta_34
  //
  // which are defined in "TOsc.h"
  //
  ///////////////////////////////////////////////////////////

  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)

  // double val_dm2_41         = 0;
  // double val_sin2_2theta_14 = it14;
  // double val_sin2_theta_24  = 0;
  // double val_sin2_theta_34  = 0;

  /// standard order
  /// standard order

  // val_dm2_41         = 0;
  // // val_sin2_2theta_14 = 0.236;
  // val_sin2_theta_24  = 0.3;
  // osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);
  // osc_test->Apply_oscillation();
  // osc_test->Set_apply_POT();// meas, CV, COV: all ready
  // osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"



  if (fcls == 1) { // Creating the pdf distribution for each point individually
    // Create the data structures
    vector<double> vec_chi2_4v_with_3vToy;
    vector<double> vec_chi2_3v_with_3vToy;
    vector<double> vec_dchi2_with_3vToy;
    vector<double> vec_chi2_4v_with_4vToy;
    vector<double> vec_chi2_3v_with_4vToy;
    vector<double> vec_dchi2_with_4vToy;


    const int NUM_TOYS = 10000;

    TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -2, 1);
    TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -2, 0);

    // convert command line flags to osc parameters
    double val_obj_dm2 = h1d_dm2->GetBinCenter(idm2);
    val_obj_dm2 = pow(10, val_obj_dm2);
    double val_obj_ttt = h1d_ttt->GetBinCenter(it14);
    val_obj_ttt = pow(10, val_obj_ttt);
    double grid_4v[4] = {val_obj_dm2, val_obj_ttt, 0, 0};
    double pars_3v_small[4] = {0, 0.10, 0.11, 0};
    cout << "Generating 3v \n";
      osc_test->Set_oscillation_pars(0, 0.1, 0.11, 0);
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT(); // meas, CV, COV: all ready
      osc_test->Set_toy_variations(NUM_TOYS);
      for (int i = 0; i < NUM_TOYS; i++) {
      //
      osc_test->Set_toy2fitdata(i + 1);
      double chi2_4v = osc_test->FCN(grid_4v);
      double chi2_3v = osc_test->FCN(pars_3v_small);
      double deltachi = chi2_4v - chi2_3v;
      vec_chi2_4v_with_3vToy.push_back(chi2_4v);
      vec_chi2_3v_with_3vToy.push_back(chi2_3v);
      vec_dchi2_with_3vToy.push_back(deltachi);
    }
      cout << "Generating 4v \n";
      osc_test->Set_oscillation_pars(grid_4v[0], grid_4v[1], grid_4v[2], grid_4v[3]);
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT(); // meas, CV, COV: all ready
      osc_test->Set_toy_variations(NUM_TOYS);
      for (int i = 0; i < NUM_TOYS; i++) {
      //
      osc_test->Set_toy2fitdata(i + 1);
      double chi2_4v = osc_test->FCN(grid_4v);
      double chi2_3v = osc_test->FCN(pars_3v_small);
      double deltachi = chi2_4v - chi2_3v;
      vec_chi2_4v_with_4vToy.push_back(chi2_4v);
      vec_chi2_3v_with_4vToy.push_back(chi2_3v);
      vec_dchi2_with_4vToy.push_back(deltachi);
      }

      // Create the tree Structure
      roostr = TString::Format("60x60_dm2_ttt_%03d_%03d.root", idm2, it14);
	  TFile *output_file = new TFile(roostr, "recreate");
	  TTree *tree = new TTree("tree", "tree");
      int grid_idx_dm2 = idm2;
      int grid_idx_ttt = it14;

	  tree->Branch("grid_idx_dm2",  &grid_idx_dm2,  "grid_idx_dm2/I" );
	  tree->Branch("grid_idx_ttt",  &grid_idx_ttt,  "grid_idx_ttt/I" );

	  tree->Branch("vec_chi2_4v_with_4vToy",   &vec_chi2_4v_with_4vToy);
	  tree->Branch("vec_chi2_3v_with_4vToy",   &vec_chi2_3v_with_4vToy);
	  tree->Branch("vec_dchi2_with_4vToy",     &vec_dchi2_with_4vToy);

	  tree->Branch("vec_chi2_4v_with_3vToy",   &vec_chi2_4v_with_3vToy);
	  tree->Branch("vec_chi2_3v_with_3vToy",   &vec_chi2_3v_with_3vToy);
	  tree->Branch("vec_dchi2_with_3vToy",     &vec_dchi2_with_3vToy);
      tree->Fill();

	  output_file->cd();
	  tree->Write();
	  output_file->Close();


  } ///
  // int obsdebug = 1;
  ///
  // if (obsdebug) {
  //   cout << "XTHROW DEBUG\n";
  //   double pars_3v_small[4] = {0, 0.10, 0.11, 0};
  //   osc_test->Set_oscillation_pars(0, 0.1, 0.11, 0);
  //   osc_test->Apply_oscillation();
  //   osc_test->Set_apply_POT(); // meas, CV, COV: all ready
  //   osc_test->Set_toy_variations(1);
  //   osc_test->Set_toy2fitdata(1);
  //   for (int universe = 1; universe <= 200; universe++) {
  //     TString fname = TString::Format("output/60x60-deltachi2_obs-xtoyuniverse-%i.root", universe);
  //     // Load 3v toy as fit data
  //     TString xpath = "input/presave_3v_hypothesis_toydata_01_cv.root";
  //     TFile *inputfile_toydata_cv = new TFile(xpath, "read");
  //     TTree *tree_toydata = (TTree*)inputfile_toydata_cv->Get("tree_toydata");
  //     // Declaration of leaf types
  //     Int_t i_toydata;
  //     vector<double> *vec_toydata_spectrum = nullptr;
  //     // List of branches
  //     TBranch        *b_i_toydata;   //!
  //     TBranch        *b_vec_toydata_spectrum;   //!
  //     // Set branch addresses and branch pointers
  //     tree_toydata->SetBranchAddress("i_toydata", &i_toydata, &b_i_toydata);
  //     tree_toydata->SetBranchAddress("vec_toydata_spectrum", &vec_toydata_spectrum, &b_vec_toydata_spectrum);
  //     tree_toydata->GetEntry(ifile+1);
  //     cout << "Loading " << i_toydata << " from Xiangpan (first pseudo=3)"  << "\n";
  //     for (int i = 0; i <= 181; i++) {
  //       double content = vec_toydata_spectrum->at(i);
  //       osc_test->matrix_tosc_fitdata_newworld(0, i) = content;
  //     }
  //     double chi2_3v = osc_test->FCN(pars_3v_small); //
  //     cout << universe << " " chi2_3v << endl;
  //   }
  // }


  if (obs_throw) { // Creating a series of reference deltachi2 parallel
    // Create 3v measurement
    double pars_3v_small[4] = {0, 0.10, 0.11, 0};
    osc_test->Set_oscillation_pars(0, 0.1, 0.11, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT(); // meas, CV, COV: all ready
    osc_test->Set_toy_variations(1);
    osc_test->Set_toy2fitdata(1);
    TString fname = TString::Format("output/60x60-deltachi2_obs-toyuniverse-%i.root", ifile);
    if (xthrow) {
      cout << "XTHROW\n";
      fname = TString::Format("output/60x60-deltachi2_obs-xtoyuniverse-%i.root", ifile);
      // Load 3v toy as fit data
      TString xpath = "input/presave_3v_hypothesis_toydata_01_cv.root";
      TFile *inputfile_toydata_cv = new TFile(xpath, "read");
      TTree *tree_toydata = (TTree*)inputfile_toydata_cv->Get("tree_toydata");
      // Declaration of leaf types
      Int_t i_toydata;
      vector<double> *vec_toydata_spectrum = nullptr;
      // List of branches
      TBranch        *b_i_toydata;   //!
      TBranch        *b_vec_toydata_spectrum;   //!
      // Set branch addresses and branch pointers
      tree_toydata->SetBranchAddress("i_toydata", &i_toydata, &b_i_toydata);
      tree_toydata->SetBranchAddress("vec_toydata_spectrum", &vec_toydata_spectrum, &b_vec_toydata_spectrum);
      tree_toydata->GetEntry(ifile+1);
      cout << "Loading " << i_toydata << " from Xiangpan (first pseudo=3)"  << "\n";
      for (int i = 0; i <= 181; i++) {
        double content = vec_toydata_spectrum->at(i);
        osc_test->matrix_tosc_fitdata_newworld(0, i) = content;
      }
      cout << "Finished loading fit data " << endl;
    }

    // Caclulate chi2_3v
    double chi2_3v = osc_test->FCN(pars_3v_small); //

    // Create the data structures needed
    TMatrixD obs_map(NUM_ttt,NUM_dm2);
    // Convert the grid into actual values

    // TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -2, 1);
    // TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -2, 0);
    cout << "Creating universe " << ifile << " \n";
    for (int indexdm2 = 0; indexdm2 < NUM_dm2; indexdm2++) {
      double val_obj_dm2 = h1d_dm2->GetBinCenter(indexdm2);
      val_obj_dm2 = pow(10, val_obj_dm2);
      for (int indextheta = 0; indextheta < NUM_ttt; indextheta++) {
        // Convert the index into real values
        double val_obj_ttt = h1d_ttt->GetBinCenter(indextheta);
        val_obj_ttt = pow(10, val_obj_ttt);
        double grid_4v[4] = {val_obj_dm2, val_obj_ttt, 0, 0};
        // Caluclate chi2_4v
        double chi2_4v = osc_test->FCN(grid_4v);
        // Calculate deltachi2
        double deltachi2_obs = chi2_4v - chi2_3v;
        obs_map(indexdm2,indextheta) = deltachi2_obs;
      }
    }

    TFile *rootfile = new TFile(fname, "recreate");
    obs_map.Write("obs_map");
    rootfile->Close();
  }

  if (analysis == 1) {
    if (ifile < 0) {
      cerr << "Specify number of universe input files with -file\n";
      exit(1);
    }
    int num_universe = ifile;
    int probdist_size = 10000;

    std::vector<TMatrixD> universe_mat_obs;
    std::vector<TMatrixD> confidence_map;
    for (int i = 1; i <= num_universe; i++) {
      // TString universestring = TString::Format("xj-recreate_confidence-map-universe-%i.root", i);
      TString universestring = TString::Format("input/60x60-deltachi2_obs-toyuniverse-%i.root", i);
      if (xthrow) {
        universestring = TString::Format("input/60x60-deltachi2_obs-xtoyuniverse-%i.root", i);
      }
      TFile deltafile(universestring, "READ");
      TMatrixD* matrix = nullptr;
      // Load universe reference for map
      deltafile.GetObject("obs_map", matrix);
      universe_mat_obs.push_back(*matrix);
      deltafile.Close();
      TMatrixD confidence(NUM_dm2, NUM_ttt);
      confidence_map.push_back(confidence);
    }
    cout << "Finished loading universe\n";//

    // Create pointers once then change the pointer for each file
    TTree *tree = nullptr;
    vector<double> *vec_dchi2_with_3vToy = nullptr;
      for (int xindex = 0; xindex < 60; xindex++) {
        for (int yindex = 0; yindex < 60; yindex++) {
          // Load deltachi2_3v, deltachi2_4v
          roostr = TString::Format("input/60x60_dm2_ttt_%03d_%03d.root", yindex + 1, xindex + 1);
          TFile *rootfile = new TFile(roostr, "read");
          rootfile->GetObject("tree", tree);
          tree->SetBranchStatus("*", 0);
          tree->SetBranchStatus("vec_dchi2_with_3vToy", 1);
          tree->SetBranchAddress("vec_dchi2_with_3vToy",&vec_dchi2_with_3vToy);
          tree->GetEntry(0);
          rootfile->Close();
          for (int obs = 0; obs < num_universe; obs++) {
              double count3v = 0;
              double count4v = 0;
              for (int j = 0; j < probdist_size; j++) {
                    double dchi3v = vec_dchi2_with_3vToy->at(j);
                    double dchi4v = vec_dchi2_with_3vToy->at(j+probdist_size);

                    if(dchi3v >= universe_mat_obs[obs](yindex,xindex)) {
                      count3v++;
                    }
                    if(dchi4v >= universe_mat_obs[obs](yindex,xindex)) {
                      count4v++;
                    }
              }
              double cls = 0;
              if( count3v == 0 ) {
                if( count4v == 0 ) cls = 0;
                else cls = 1;
              }
              else cls = count4v / count3v;
              if( count4v>=count3v ) cls = 1;
              double confidence = 1 - cls;
              confidence_map[obs](yindex,xindex) = confidence;
          }
          vec_dchi2_with_3vToy->clear();

        } // yindex
        cout << "Finished x index: " << xindex+1 << "\n";
      } // xindex
      TString outname = "output/confidence_maps_grid.root";
      TFile outfile(outname, "RECREATE");
      int i = 1;
      for (const auto& m : confidence_map) {
        outfile.WriteObject(&m, Form("universe_%d", i++));
      }
      // i = 1;
      // for (const auto& m : universe_mat_obs) {
      //   outfile.WriteObject(&m, Form("uobs_%d", i++));
      // }
      outfile.Close(); //
      cout << "Wrote " << outname << "\n";
  }

  if (sensitivity_study == 0 && idm2 != 0) {

    // #include "MendezStyle.h"
    // void test() {
    //   int dm2_val = 40;
    //   vector<double> xvals;
    //   cout << "Sensitivity across constant dm2\n";
    //   double cls_target = 0.95;
    //   TString confidencein = "input/confidence_maps_grid.root";
    //   TFile confidence_maps(confidencein, "READ");
    //   int num_universe = confidence_maps.GetNkeys();
    //   TMatrixD* matrix = nullptr;
    //   confidence_maps.GetObject("universe_1", matrix);
    //   int ncols = (*matrix).GetNcols();
    //   int nrows = (*matrix).GetNrows();
    //   TH1D hsens("h1d_ttt", "h1d_ttt", 60, -2, 0);
    //   for (int i=1; i<=num_universe; i++) {
    //     TString matname = TString::Format("universe_%i", i);
    //     // Load universe reference for map
    //     confidence_maps.GetObject(matname, matrix);
    //     // Create a histogram from the matrix
    //     TH2D htemp("htemp", "Parameter Map;X;Y", ncols, -2, 0, nrows, -2, 1);
    //     for (int xindex=0; xindex<ncols; xindex++) {
    //       for (int yindex=0; yindex<nrows; yindex++) {
    //         htemp.SetBinContent(xindex+1,yindex+1,(*matrix)(yindex,xindex));
    //       }
    //     }
    //     std::vector<TGraph*> contour = mendezstyle::GetContourGraphs(&htemp,
    //     cls_target); int graph_index = 0; int index_biggest = 0; int
    //     largest_size = 0; for (const auto& graph : contour ) {
    //       if (graph->GetN()  > largest_size) {
    //         index_biggest = graph_index;
    //         largest_size = graph->GetN();
    //       }
    //       graph_index++;
    //     }
    //     // Pull Out graph and invert
    //     TGraph temp(contour[index_biggest]->GetN(),
    //     contour[index_biggest]->GetY(), contour[index_biggest]->GetX());
    //     double x = temp.Eval(dm2_val);
    //     hsens.Fill(x);
    //     xvals.push_back(x);
    //   }
    //   std::sort(xvals.begin(),xvals.end());//Sorting the vector
    //   double p[5] = {0.023, 0.159, 0.5, .841, .977};
    //   double qv[5];
    //   for (int i=0; i<5; i++) {
    //     qv[i] = p[i]*xvals.size();
    //     cout << qv[i] << "\n";
    //   }
  }
//   if (0) { // Debug Output
//
//     const int xdims = 61;
//     const int ydims = 81;
//
//     TH1D* h1_dm2 = new TH1D("dm2", "dm2", ydims, -2,2);
//     TH1D *h1_sin2_2tuu = new TH1D("sin2_2tuu", "sin2_2tuu", xdims, -2, 0);
//     double dm2_41_grid = h1_dm2->GetBinCenter(idm2);
//     dm2_41_grid = pow(10.0, dm2_41_grid);
//     double theta_grid = h1_sin2_2tuu->GetBinCenter(it14);
//     theta_grid = pow(10.0, theta_grid);
//     double pars_4v_grid[4] ={dm2_41_grid, theta_grid, 0.0045, 0};
//     double pars_3v_small[4] = {0, 0.1, 0.1, 0};
//     const int NUM_TOYS = 5;
//     double deltachi2_obs[NUM_TOYS];
//     // Generate 3v psuedo
//     osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1], pars_3v_small[2], pars_3v_small[3]);
//     osc_test->Apply_oscillation();
//     osc_test->Set_apply_POT();// meas, CV, COV: all ready
//     osc_test->Set_asimov2fitdata(); // set the "asimov toy sample" as the
//     TString fname = "asimov_prediction.root";
//
//     TFile outfile(fname, "RECREATE");
//     // matrix_tosc_eff_newworld_noosc
//     outfile.WriteObject(&osc_test->matrix_tosc_fitdata_newworld, "asimov_mat");//
//     osc_test->Set_toy_variations(NUM_TOYS); // produce NUM_TOYS pseudo
//   }
//   cout<<endl;
  // cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
  // cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

  // cout<<endl;
  // cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f, theta14(-it14) %d, dm2(-idm2) %d, theta24(-it24) %d",
  // 			display, ifile, osc_test->tosc_scaleF_POT_BNB, osc_test->tosc_scaleF_POT_NuMI, it14, idm2, it24)<<endl;
  // cout<<endl;

  // cout<<" ---> flag_apply_oscillation_BNB  "<< osc_test->flag_apply_oscillation_BNB<<endl;
  // cout<<" ---> flag_apply_oscillation_NuMI "<< osc_test->flag_apply_oscillation_NuMI<<endl;
  // cout<<endl;

  // cout<<" ---> flag_goodness_of_fit_CNP    "<< osc_test->flag_goodness_of_fit_CNP<<endl;
  // cout<<endl;

  // cout<<TString::Format(" ---> flag_syst_dirt    %d", osc_test->flag_syst_dirt)<<endl;
  // cout<<TString::Format(" ---> flag_syst_mcstat  %d", osc_test->flag_syst_mcstat)<<endl;
  // cout<<TString::Format(" ---> flag_syst_flux    %d", osc_test->flag_syst_flux)<<endl;
  // cout<<TString::Format(" ---> flag_syst_geant   %d", osc_test->flag_syst_geant)<<endl;
  // cout<<TString::Format(" ---> flag_syst_Xs      %d", osc_test->flag_syst_Xs)<<endl;
  // cout<<TString::Format(" ---> flag_syst_det     %d", osc_test->flag_syst_det)<<endl;

  // cout<<endl;
  cout<<" ---> Finished successfully"<<endl;


  cout<<endl;
  if( display ) {
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<endl;
    theApp.Run();
  }

  return 0;
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

  for(int i = 0; i < TotalConts; i++) {
    contLevel = (TList*)conts->At(i);

    if(i!=0) continue;// only use the one having the most points

    for(int j = 0; j < contLevel->GetSize(); j++) {
      curv = (TGraph*)contLevel->At(j);
      int Npoints = curv->GetN();
      printf(" ---> graph %2d, size %3d\n", j+1, Npoints);

      if(i==0) {
	for(int ip=0; ip<Npoints; ip++) {
	  double dm2_val, ttt_val;
	  curv->GetPoint(ip, ttt_val, dm2_val);

	  //gh_CL_curve->SetPoint( gh_CL_curve->GetN(), pow(10, ttt_val), pow(10, dm2_val) );

	  gh_CL_curve->SetPoint( gh_CL_curve->GetN(), ttt_val, dm2_val );

	}// for(int ip=0; ip<Npoints; ip++)
      }

    }// for(int j = 0; j < contLevel->GetSize(); j++)
  }// for(int i = 0; i < TotalConts; i++)

  cout<<endl;

  gROOT->SetBatch( 0 );
}





  /// Obtain the TGraph(s) corresponding to a particular contour level for a TH2
  ///
  /// \param h2     The TH2 to examine
  /// \param level  Fractional level in [0, 1]
  /// \return       Vector of TGraphs that represent the whole contour (multiple if contour has disjoint pieces)
  std::vector<TGraph*> GetContourGraphs(TH2* h2, double level)
  {
    std::vector<TGraph*> ret;

    std::unique_ptr<TH2> surf(dynamic_cast<TH2*>(h2->Clone("tmp_h2_for_drawing_graphs")));

    TVirtualPad* bak = gPad;

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
    TCollection* specs = gROOT->GetListOfSpecials();

    TIter nextSpec(specs);
    while(TObject* spec = nextSpec()){
      if(!spec->InheritsFrom(TObjArray::Class())) continue;
      auto conts = dynamic_cast<TObjArray*>(spec);

      if(conts->IsEmpty()) continue;

      if(!conts->At(0)->InheritsFrom(TList::Class())) continue;
      auto cont = dynamic_cast<TList*>(conts->At(0));

      TIter nextObj(cont);
      // Contour could be split into multiple pieces
      std::size_t piece = 0;
      while(TObject* obj = nextObj()){
        if(!obj->InheritsFrom(TGraph::Class())) continue;

        ret.push_back(dynamic_cast<TGraph*>(obj->Clone(Form("%s_contour%f_piece%zu", obj->GetName(), level, piece))));
        piece++;
      } // end for obj
    } // end for spec

    return ret;
  }
