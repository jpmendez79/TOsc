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


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A story ..."<<endl<<endl;

  int ifile = 1;
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;

  int it14 = 0;
  int idm2 = 0;
  int it24 = 0;
  // double it14 = 0;
  // double idm2 = 0;
  // double it24 = 0;  
  
  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
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

  double val_dm2_41         = 0;
  double val_sin2_2theta_14 = it14;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;

  /// standard order
  /// standard order
  
  val_dm2_41         = 0;
  // val_sin2_2theta_14 = 0.236;
  val_sin2_theta_24  = 0.3;  
  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
  //osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"

  /// one example: calculate a chi2
  /// one example: calcualte a chi2


  if ( 0 ) { // Heat map generation Jesse Mendez Current Study
    int dimensions = 60;
    TH1D* h1_dm2 = new TH1D("dm2", "dm2", 80, -2,2);
    TH1D* h1_sin2_2tuu = new TH1D("sin2_2tuu", "sin2_2tuu", dimensions, -2,0);

    cout << "Heat Map Generation\n";
    double dm2_41_grid = h1_dm2->GetBinCenter(idm2);
    dm2_41_grid = pow(10.0, dm2_41_grid);
    double theta_grid = h1_sin2_2tuu->GetBinCenter(it14);
    theta_grid = pow(10.0, theta_grid);
    double pars_4v_grid[4] ={dm2_41_grid, theta_grid, 0.0045, 0};
    double pars_3v_small[4] = {0, 0.1, 0.1, 0};
    const int NUM_TOYS = 2500;
  
    double step_size = 0.0005; // profiling step size
    int steps = 2000; // profiling steps

    // Profiling data structures
    double profiled_sin2_theta24 = 0;
    double test_sin2_theta24 = 0;
    double chi2_min = pow(10,6);	// set very large to start

    // Data Structures for Chi^2 arrays
    // 4v stuctures
    double array_4vhyp_deltachisquare_grid_toy[NUM_TOYS];
    double array_4vhyp_4nu_chisquare_grid[NUM_TOYS];
    double array_4vhyp_4nu_chisquare_3nu[NUM_TOYS];

    double array_3vhyp_4nu_chisquare_grid[NUM_TOYS];
    double array_3vhyp_4nu_chisquare_3nu[NUM_TOYS];
    double array_3vhyp_deltachisquare_grid_toy[NUM_TOYS];

    // Reference deltachisquare for CLs
    double ref = 0;

    // cout << "Profile \n";
    // // Grid scan
    // for(int i = 0; i < steps; i++) {
    //   test_sin2_theta24 = i * step_size;
    //   if( theta_grid > test_sin2_theta24 ) {
    // 	continue; 		// Need to because I am inputing mixing ange directly
    //   }	
    //   double pars_4nu[4] = {dm2_41_grid, theta_grid, test_sin2_theta24, 0};// dm2, t14, t24, t34
    //   double chi2 = osc_test->FCN( pars_4nu );
    //   if(chi2 < chi2_min) {
    // 	chi2_min = chi2;
    // 	profiled_sin2_theta24 = test_sin2_theta24;
    //   }
    // }
    // pars_4v_grid[2] = profiled_sin2_theta24;
    cout << "4nu \n";
    // Generating psuedo experiments
    osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    // osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
    osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments
    for (int i = 0; i < NUM_TOYS; i++) {
      osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
      array_4vhyp_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
      array_4vhyp_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
      array_4vhyp_deltachisquare_grid_toy[i] = array_4vhyp_4nu_chisquare_grid[i] - array_4vhyp_4nu_chisquare_3nu[i];
    }
    cout << "3nu \n";
    // Generating psuedo experiments
    osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1], pars_3v_small[2], pars_3v_small[3]);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
    osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments
    for (int i = 0; i < NUM_TOYS; i++) {
      osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
      array_3vhyp_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
      array_3vhyp_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
      array_3vhyp_deltachisquare_grid_toy[i] = array_3vhyp_4nu_chisquare_grid[i] - array_3vhyp_4nu_chisquare_3nu[i];
      if (isnan(array_3vhyp_4nu_chisquare_grid[i]) || isnan(array_3vhyp_4nu_chisquare_3nu[i] ) || isnan(array_3vhyp_deltachisquare_grid_toy[i])) {
	cout <<"NaaaaaaNNNNNNNNN\n";
      } 
    }
    cout << "Delta Chisquare Obs\n";
    osc_test->Set_oscillation_pars(0, 0.1, 0.1, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    // osc_test->Set_asimov2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    double chi2_testA = osc_test->FCN( pars_3v_small );
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    // osc_test->Set_asimov2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    double chi2_testC = osc_test->FCN( pars_4v_grid );
    ref = chi2_testC - chi2_testA;
    // Create histograms for 3v/4v delta_chi2_grid_toy
    // double min3 = 1e30, max3 = -1e30;
    // double min4 = 1e30, max4 = -1e30;
    // for (int j = 0; j < NUM_TOYS; ++j) {
    //   min4 = std::min(min4, array_4vhyp_deltachisquare_grid_toy[j]); max4 = std::max(max4, array_4vhyp_deltachisquare_grid_toy[j]);
    //   min3 = std::min(min3, array_3vhyp_deltachisquare_grid_toy[j]); max3 = std::max(max3, array_3vhyp_deltachisquare_grid_toy[j]);
    // }
    // TH1D *h1 = new TH1D("h1", "p3v;Delta chi2;Entries", 100, min3, max3);
    // TH1D *h2 = new TH1D("h2", "p4v;Delta chi2;Entries", 100, min4, max4);
    // for (int j = 0; j < NUM_TOYS; ++j) {
    //   // cout << array_3vhyp_deltachisquare_grid_toy[j] << ":" << array_4vhyp_deltachisquare_grid_toy[j] <<"\n";
    //   h1->Fill(array_3vhyp_deltachisquare_grid_toy[j]);
    //   h2->Fill(array_4vhyp_deltachisquare_grid_toy[j]);
    // }
      
    // h1->Scale(1.0 / h1->Integral());
    // h2->Scale(1.0 / h2->Integral());
    // // Calculate CLs from histogram bin
    // int bin4v = h2->FindBin(ref);
    // int bin3v = h1->FindBin(ref);
    // double p3v = h1->Integral(bin3v, h1->GetNbinsX());
    // double p4v = h2->Integral(bin4v, h2->GetNbinsX());
    // double val_CLs = p4v / p3v;
    // Save Everything I need to files
    TString fname = TString::Format("./cls_inv_decay_bnb_only_g2_25pi_grid_60x80_4Uu4square_%d_dm2_%d.root", it14, idm2);
    TFile *rootfile = new TFile(fname, "recreate");
    // Create objects for cpp variables
    TParameter<double> deltachi2_ref("deltachi2_ref", ref);
    TParameter<int> psuedos("num_toys", NUM_TOYS);
    TParameter<double> theta("theta_grid", it14);
    TParameter<double> dm2("dm2_grid", idm2);
    // TParameter<double> cls("CLs", val_CLs);
    // Create TArray objects for all cpp arrays
    // 4v stuctures
    TVectorD arr4d(NUM_TOYS, array_4vhyp_deltachisquare_grid_toy);
    TArrayD arr4g(NUM_TOYS, array_4vhyp_4nu_chisquare_grid);
    TArrayD arr43(NUM_TOYS, array_4vhyp_4nu_chisquare_3nu);
    TVectorD arr3d(NUM_TOYS, array_3vhyp_deltachisquare_grid_toy);
    TArrayD arr3g(NUM_TOYS, array_3vhyp_4nu_chisquare_grid);
    TArrayD arr33(NUM_TOYS, array_3vhyp_4nu_chisquare_3nu);
    deltachi2_ref.Write();
    psuedos.Write();
    // cls.Write();
    theta.Write();
    dm2.Write();
    arr4d.Write("4v_deltachi2");
    arr3d.Write("3v_deltachi2");
    // h1->Write();
    // h2->Write();
    rootfile->Write();
    rootfile->Close();
    delete rootfile;
  }// if

  ///////////////////////////////////////////////////////////
  /// Brazil Band Code                               /// 
  ///////////////////////////////////////////////////////////

  if (1) {
    // Read in the parameters
    int xdims = 60;
    int ydims = 80;
    TH1D* h1_dm2 = new TH1D("dm2", "dm2", ydims, -2,2);
    TH1D* h1_sin2_2tuu = new TH1D("sin2_2tuu", "sin2_2tuu", xdims, -2,0);

    cout << "Deltachi2 Vector\n";
    double dm2_41_grid = h1_dm2->GetBinCenter(idm2);
    dm2_41_grid = pow(10.0, dm2_41_grid);
    double theta_grid = h1_sin2_2tuu->GetBinCenter(it14);
    theta_grid = pow(10.0, theta_grid);
    double pars_4v_grid[4] ={dm2_41_grid, theta_grid, 0.0045, 0};
    double pars_3v_small[4] = {0, 0.1, 0.1, 0};
    const int NUM_TOYS = 1000;
    double deltachi2_obs[NUM_TOYS];
    // Generate 3v psuedo
    osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_asimov2fitdata(); // set the "asimov toy sample" as the
    osc_test->Set_toy_variations( NUM_TOYS ); // produce NUM_TOYS pseudo experiments
    // Calculate the deltachisquare and save
    for (int i = 0; i < NUM_TOYS; i++) {
      osc_test->Set_toy2fitdata(i + 1);
      double chi2_3v = osc_test->FCN(pars_3v_small);
      double chi2_4v = osc_test->FCN( pars_4v_grid );
      deltachi2_obs[i] = chi2_4v - chi2_3v;
    }
    // Write out everything to file
    TString fname = TString::Format("./vanilla-numu_grid_60x80_sinsquare_theta_uu_%d_dm2_%d-obs_deltachi2.root", it14, idm2);
    TFile *rootfile = new TFile(fname, "recreate");
    // Create objects for cpp variables
    TParameter<int> psuedos("num_toys", NUM_TOYS);
    TParameter<double> theta("theta_grid", it14);
    TParameter<double> dm2("dm2_grid", idm2);
    TVectorD arrdeltachi2(NUM_TOYS, deltachi2_obs);
    psuedos.Write();
    theta.Write();
    dm2.Write();
    arrdeltachi2.Write("deltachi2_obs");
    // h1->Write();
    // h2->Write();
    rootfile->Write();
    rootfile->Close();
    delete rootfile;

  }// End Brazil band code
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////

  cout<<endl;
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
