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

  cout<<endl<<" ---> A Hello story ..."<<endl<<endl;

  int ifile = -1;
  int obs_throw = 0;
  int fcls = 0;
  int xthrow = 1;
  int draw_confidence_map = 0;
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;

  int it14 = 0;
  int idm2 = 0;
  int it24 = 0;
  int option = 0;
  // double it14 = 0;
  // double idm2 = 0;
  // double it24 = 0;

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-file")==0 ) {
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
    if( strcmp(argv[i],"-o")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>option ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
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


//   if (option == 0) {
//     cout << "Enter an option flag with -o N where N is a flag from 1->N\n";
//     cout << "Creating a set of universes using xiangpan\n";
//     map<int, TMatrixD>map_toydata_spectrum;
//     map<int, double>map_prof_t24_val;
//     int rows = -1;
//
//     TString xpath = "presave_3v_hypothesis_toydata_01_cv.root";
//     TFile *inputfile_toydata_cv = new TFile(xpath, "read");
//     TTree *tree_toydata = (TTree*)inputfile_toydata_cv->Get("tree_toydata");
//
//     // Declaration of leaf types
//     Int_t           i_toydata;
//     vector<double>  *vec_toydata_spectrum;
//
//     // List of branches
//     TBranch        *b_i_toydata;   //!
//     TBranch        *b_vec_toydata_spectrum;   //!
//
//     // Set object pointer
//     vec_toydata_spectrum = 0;
//
//     // Set branch addresses and branch pointers
//     tree_toydata->SetBranchAddress("i_toydata", &i_toydata, &b_i_toydata);
//     tree_toydata->SetBranchAddress("vec_toydata_spectrum", &vec_toydata_spectrum, &b_vec_toydata_spectrum);
//
//
//     int entries = tree_toydata->GetEntries();
//     cout<<endl<<" ---> reading tree_toydata, entries "<<entries<<endl<<endl;
//     tree_toydata->GetEntry(ifile);
//     rows = vec_toydata_spectrum->size();
//
//     cout << "Creating universe " << ifile<< " \n";
// //     const int xdims = 60;
// //     const int ydims = 60;
// //     TH1D* h1_dm2 = new TH1D("dm2", "dm2", ydims, -2,2);
// //     TH1D *h1_sin2_2tuu = new TH1D("sin2_2tuu", "sin2_2tuu", xdims, -3, 0);
// //     double obs_map[xdims+1][ydims+1];
// //     double pars_3v_small[4] = {0, 0.1, 0.1, 0};
// //     osc_test->Set_oscillation_pars(0, 0.1, 0.1, 0);
// //     osc_test->Apply_oscillation();
// //     osc_test->Set_apply_POT(); // meas, CV, COV: all ready
// //     // osc_test->Set_toy_variations(10);
// //     // cout << "DEBUG 3v" << endl;
// //     // double chi2_3v = osc_test->FCN(pars_3v_small); //
// //     for (int indexdm2 = 0; indexdm2 <= ydims; indexdm2++) {
// //       for (int indextheta = 0; indextheta <= xdims; indextheta++) {
// //     	double dm2_41_grid = h1_dm2->GetBinCenter(indexdm2);
// //     	dm2_41_grid = pow(10.0, dm2_41_grid);
// //     	double theta_grid = h1_sin2_2tuu->GetBinCenter(indextheta);
// //     	theta_grid = pow(10.0, theta_grid);
// //         double pars_4v_grid[4] = {dm2_41_grid, theta_grid, 0.0045, 0};
// //         tree_toydata->GetEntry( ifile );
// //         for(int idx=0; idx<rows; idx++) {
// //           double content = vec_toydata_spectrum->at(idx);
// //           osc_test->matrix_tosc_fitdata_newworld(0,idx) = content;
// // }
//         // cout << "Loading " << ifile << " from Xiangpan" << "/n";
//     	// osc_test->Set_toy2fitdata(1);
//     	// osc_test->Set_asimov2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//         double chi2_3v = osc_test->FCN(pars_3v_small); //
//         // tree_toydata->GetEntry( ifile );
//         // for(int idx=0; idx<rows; idx++)
//         // cout << "Loading " << ifile << " from Xiangpan" << "/n";
//         // osc_test->Set_toy2fitdata(1);
//     	// osc_test->Set_asimov2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//     	double chi2_4v = osc_test->FCN( pars_4v_grid );
//     	double ref = chi2_4v - chi2_3v;
//         obs_map[indextheta][indexdm2] = ref;
//       }
//     }
// //     TMatrixD matrix(xdims+1, ydims+1);
// //     for (int j = 0; j <= ydims; j++)
// //       for (int i = 0; i <= xdims; i++)
// //         matrix(i,j) = obs_map[i][j];
// //     // obs_tree.Branch("obs_map", &obs_map[xdims][ydims]);
// //     // obs_tree.Fill();
// //     TString fname = TString::Format("xiangpan-test-60x60/vanilla-numu_grid_60x60_xiangpan-toyuniverse_%i.root", ifile);
// //     TFile *rootfile = new TFile(fname, "recreate");
// //     matrix.Write("obs_map");
// //     rootfile->Close();
//   }
//
//   if (option == 1) {
//     cout << "Option " << option << "\n";
//
//
//   }

  //osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"

  /// one example: calculate a chi2
  /// one example: calcualte a chi2


//   if (0) { // Heat map generation Jesse Mendez Current Study
//     int dimensions = 60;
//     TH1D* h1_dm2 = new TH1D("dm2", "dm2", 60, -2,2);
//     TH1D* h1_sin2_2tuu = new TH1D("sin2_2tuu", "sin2_2tuu", dimensions, -3,0);
//
//     cout << "Heat Map Generation\n";
//     double dm2_41_grid = h1_dm2->GetBinCenter(idm2);
//     dm2_41_grid = pow(10.0, dm2_41_grid);
//     double theta_grid = h1_sin2_2tuu->GetBinCenter(it14);
//     theta_grid = pow(10.0, theta_grid);
//     double pars_4v_grid[4] ={dm2_41_grid, theta_grid, 0.0045, 0};
//     double pars_3v_small[4] = {0, 0.1, 0.1, 0};
//     const int NUM_TOYS = 10000;
//
//     double step_size = 0.0005; // profiling step size
//     int steps = 2000; // profiling steps
//
//     // Profiling data structures
//     double profiled_sin2_theta24 = 0;
//     double test_sin2_theta24 = 0;
//     double chi2_min = pow(10,6);	// set very large to start
//
//     // Data Structures for Chi^2 arrays
//     // 4v stuctures
//     double array_4vhyp_deltachisquare_grid_toy[NUM_TOYS];
//     double array_4vhyp_4nu_chisquare_grid[NUM_TOYS];
//     double array_4vhyp_4nu_chisquare_3nu[NUM_TOYS];
//
//     double array_3vhyp_4nu_chisquare_grid[NUM_TOYS];
//     double array_3vhyp_4nu_chisquare_3nu[NUM_TOYS];
//     double array_3vhyp_deltachisquare_grid_toy[NUM_TOYS];
//
//     // Reference deltachisquare for CLs
//     // double ref = 0;
//
//     // cout << "Profile \n";
//     // // Grid scan
//     // for(int i = 0; i < steps; i++) {
//     //   test_sin2_theta24 = i * step_size;
//     //   if( theta_grid > test_sin2_theta24 ) {
//     // 	continue; 		// Need to because I am inputing mixing ange directly
//     //   }
//     //   double pars_4nu[4] = {dm2_41_grid, theta_grid, test_sin2_theta24, 0};// dm2, t14, t24, t34
//     //   double chi2 = osc_test->FCN( pars_4nu );
//     //   if(chi2 < chi2_min) {
//     // 	chi2_min = chi2;
//     // 	profiled_sin2_theta24 = test_sin2_theta24;
//     //   }
//     // }
//     // pars_4v_grid[2] = profiled_sin2_theta24;
//     cout << "4nu \n";
//     // Generating psuedo experiments
//     osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);
//     osc_test->Apply_oscillation();
//     osc_test->Set_apply_POT();// meas, CV, COV: all ready
//     // osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//     osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
//     osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments
//     for (int i = 0; i < NUM_TOYS; i++) {
//       osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
//       array_4vhyp_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
//       array_4vhyp_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
//       array_4vhyp_deltachisquare_grid_toy[i] = array_4vhyp_4nu_chisquare_grid[i] - array_4vhyp_4nu_chisquare_3nu[i];
//     }
//     cout << "3nu \n";
//     // Generating psuedo experiments
//     osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1], pars_3v_small[2], pars_3v_small[3]);
//     osc_test->Apply_oscillation();
//     osc_test->Set_apply_POT();// meas, CV, COV: all ready
//     //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//     osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
//     osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments
//     for (int i = 0; i < NUM_TOYS; i++) {
//       osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
//       array_3vhyp_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
//       array_3vhyp_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
//       array_3vhyp_deltachisquare_grid_toy[i] = array_3vhyp_4nu_chisquare_grid[i] - array_3vhyp_4nu_chisquare_3nu[i];
//       if (isnan(array_3vhyp_4nu_chisquare_grid[i]) || isnan(array_3vhyp_4nu_chisquare_3nu[i] ) || isnan(array_3vhyp_deltachisquare_grid_toy[i])) {
// 	    cout <<"NaaaaaaNNNNNNNNN\n";
//       }
//     }
//     // cout << "Delta Chisquare Obs\n";
//     // osc_test->Set_oscillation_pars(0, 0.1, 0.1, 0);
//     // osc_test->Apply_oscillation();
//     // osc_test->Set_apply_POT();// meas, CV, COV: all ready
//     // osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//     // // osc_test->Set_asimov2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//     // double chi2_testA = osc_test->FCN( pars_3v_small );
//     // osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//     // // osc_test->Set_asimov2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//     // double chi2_testC = osc_test->FCN( pars_4v_grid );
//     // ref = chi2_testC - chi2_testA;
//     // Save Everything I need to files
//     TString fname = TString::Format("./202604-vanilla-numu-disappear-60x60-grid-sinsquare-twothetauu_%i-deltamsquare41_%i.root.probdist", it14, idm2);
//     TFile *rootfile = new TFile(fname, "recreate");
//     // Create objects for cpp variables
//     // TParameter<double> deltachi2_ref("deltachi2_ref", ref);
//     TParameter<int> psuedos("pdf_size", NUM_TOYS);
//     TParameter<double> theta("theta_grid", it14);
//     TParameter<double> dm2("dm2_grid", idm2);
//     // TParameter<double> cls("CLs", val_CLs);
//     // Create TArray objects for all cpp arrays
//     // 4v stuctures
//     TVectorD arr4d(NUM_TOYS, array_4vhyp_deltachisquare_grid_toy);
//     TVectorD arr3d(NUM_TOYS, array_3vhyp_deltachisquare_grid_toy);
//     // deltachi2_ref.Write();
//     psuedos.Write();
//     // cls.Write();
//     theta.Write();
//     dm2.Write();
//     arr4d.Write("4v_deltachi2");
//     arr3d.Write("3v_deltachi2");
//     rootfile->Write();
//     rootfile->Close();
//     delete rootfile;
//   }// if
//
//   ///////////////////////////////////////////////////////////
//   /// Brazil Band Code                               ///
//   ///////////////////////////////////////////////////////////
//
//   if (0) {
//     // Read in the parameters
//     int xdims = 60;
//     int ydims = 80;
//     TH1D* h1_dm2 = new TH1D("dm2", "dm2", ydims, -2,2);
//     TH1D* h1_sin2_2tuu = new TH1D("sin2_2tuu", "sin2_2tuu", xdims, -2,0);
//
//     cout << "Deltachi2 Vector\n";
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
//     osc_test->Set_toy_variations( NUM_TOYS ); // produce NUM_TOYS pseudo experiments
//     // Calculate the deltachisquare and save
//     for (int i = 0; i < NUM_TOYS; i++) {
//       osc_test->Set_toy2fitdata(i + 1);
//       double chi2_3v = osc_test->FCN(pars_3v_small);
//       double chi2_4v = osc_test->FCN( pars_4v_grid );
//       deltachi2_obs[i] = chi2_4v - chi2_3v;
//     }
//     // Write out everything to file
//     TString fname = TString::Format("./vanilla-numu_grid_60x80_sinsquare_theta_uu_%d_dm2_%d.root.deltachi2", it14, idm2);
//     TFile *rootfile = new TFile(fname, "recreate");
//     // Create objects for cpp variables
//     TParameter<int> psuedos("num_toys", NUM_TOYS);
//     TParameter<double> theta("theta_grid", it14);
//     TParameter<double> dm2("dm2_grid", idm2);
//     TVectorD arrdeltachi2(NUM_TOYS, deltachi2_obs);
//     psuedos.Write();
//     theta.Write();
//     dm2.Write();
//     arrdeltachi2.Write("deltachi2_obs");
//     // h1->Write();
//     // h2->Write();
//     rootfile->Write();
//     rootfile->Close();
//     delete rootfile;
//
//   } // End Brazil band code
//
  // if (0) { // Creating a series of reference deltachi2


    // cout << "Creating tree of deltachi2_obs \n";
  //   TFile outfile("deltachi2_universe.root", "RECREATE");
  //   TTree obs_tree("obs_tree", "Tree of deltachi2 obs values");
  //     const int xdims = 60;
  //     const int ydims = 80;
  //     TH1D* h1_dm2 = new TH1D("dm2", "dm2", ydims, -2,2);
  //     TH1D *h1_sin2_2tuu = new TH1D("sin2_2tuu", "sin2_2tuu", xdims, -3, 0);
  //
  //     double pars_3v_small[4] = {0, 0.1, 0.1, 0};
  //     const int NUM_TOYS = 1000;
  //     double deltachi2_obs[NUM_TOYS];
  //     // Generate 3v psuedo
  //     osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1],
  //     pars_3v_small[2], pars_3v_small[3]); osc_test->Apply_oscillation();
  //     osc_test->Set_apply_POT();// meas, CV, COV: all ready
  //     osc_test->Set_asimov2fitdata(); // set the "asimov toy sample" as the
  //     osc_test->Set_toy_variations(NUM_TOYS); // produce NUM_TOYS pseudo
  //
  //     for (int universe = 0; universe < NUM_TOYS; universe++) {
  //       cout << "Starting toy " << universe << "\n";
  //       double obs_map[xdims][ydims];
  //       for (int indexdm2 = 0; indexdm2 < ydims; indexdm2++) {
  // 	    for (int indextheta = 0; indextheta < xdims; indextheta++) {
  // 	      double dm2_41_grid = h1_dm2->GetBinCenter(indexdm2);
  // 	      dm2_41_grid = pow(10.0, dm2_41_grid);
  // 	      double theta_grid = h1_sin2_2tuu->GetBinCenter(indextheta);
  // 	      theta_grid = pow(10.0, theta_grid);
  //           double pars_4v_grid[4] = {dm2_41_grid, theta_grid, 0.0045, 0};
  //
  //           osc_test->Set_toy2fitdata(universe + 1);
  // 	      double chi2_3v = osc_test->FCN(pars_3v_small);
  // 	      double chi2_4v = osc_test->FCN( pars_4v_grid );
  //           obs_map[indextheta][indexdm2] = chi2_4v - chi2_3v;
  //         }
  //       }
  //       obs_tree.Branch("obs_map", &obs_map[xdims][ydims]);
  //       obs_tree.Fill();
  //       cout << "Finishing toy " << universe << "\n";
  //     }
  //     obs_tree.Write();
  //     outfile.Close();
  //   }    //
  //   ///////////////////////////////////////////////////////////
  //   ///////////////////////////////////////////////////////////
  if (fcls == 1) { // Creating the pdf distribution for each point individually
    // Create the data structures
    vector<double> vec_chi2_4v_with_3vToy;
    vector<double> vec_chi2_3v_with_3vToy;
    vector<double> vec_dchi2_with_3vToy;
    vector<double> vec_chi2_4v_with_4vToy;
    vector<double> vec_chi2_3v_with_4vToy;
    vector<double> vec_dchi2_with_4vToy;

    const int NUM_dm2 = 60;
    const int NUM_ttt = 60;
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


  }///
  if (obs_throw) { // Creating a series of reference deltachi2 parallel
    // Create 3v measurement
    double pars_3v_small[4] = {0, 0.10, 0.11, 0};
    osc_test->Set_oscillation_pars(0, 0.1, 0.11, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT(); // meas, CV, COV: all ready
    osc_test->Set_toy_variations(1);
    osc_test->Set_toy2fitdata(1);
    // Caclulate chi2_3v
    double chi2_3v = osc_test->FCN(pars_3v_small); //

    // Create the data structures needed
    // TTree obs_tree("obs_tree", "Tree of deltachi2 obs values");
    const int NUM_dm2 = 60;
    const int NUM_ttt = 60;
    // double obs_map[NUM_ttt][NUM_dm2];
    TMatrixD obs_map(NUM_ttt,NUM_dm2);
    // Convert the grid into actual values

    TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -2, 1);
    TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -2, 0);
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
    TString fname = TString::Format("60x60-deltachi2_obs-toyuniverse-%i.root", ifile);
    TFile *rootfile = new TFile(fname, "recreate");
    obs_map.Write("obs_map");
    rootfile->Close();
  }

  if (xthrow) {
    // Create 3v measurement
    double pars_3v_small[4] = {0, 0.10, 0.11, 0};
    osc_test->Set_oscillation_pars(0, 0.1, 0.11, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT(); // meas, CV, COV: all ready
    osc_test->Set_toy_variations(1);
    osc_test->Set_toy2fitdata(1);
      // Load 3v toy as fit data
  TString xpath = "presave_3v_hypothesis_toydata_01_cv.root";
  TFile *inputfile_toydata_cv = new TFile(xpath, "read");
  TTree *tree_toydata = (TTree*)inputfile_toydata_cv->Get("tree_toydata");
  // Declaration of leaf types
  Int_t i_toydata;
  vector<double> *vec_toydata_spectrum = nullptr;
  // vector<double>  *vec_toydata_spectrum;
  // List of branches
  TBranch        *b_i_toydata;   //!
  TBranch        *b_vec_toydata_spectrum;   //!
  // Set branch addresses and branch pointers
  tree_toydata->SetBranchAddress("i_toydata", &i_toydata, &b_i_toydata);
  tree_toydata->SetBranchAddress("vec_toydata_spectrum", &vec_toydata_spectrum, &b_vec_toydata_spectrum);
  // int entries = tree_toydata->GetEntries();
  // cout << endl<< " ---> reading tree_toydata, entries " << entries << endl << endl;c
  cout << "Loading " << ifile << " from Xiangpan"  << "\n";
  tree_toydata->GetEntry(ifile);

    // int rows = vec_toydata_spectrum->size();
    for (int i = 0; i <= 181; i++) {
      double content = vec_toydata_spectrum->at(i);
      osc_test->matrix_tosc_fitdata_newworld(0, i) = content;
      // cout << osc_test->matrix_tosc_fitdata_newworld(0, i) << " " << content <<"\n";
    }
    cout << "Finished loading fit data " << endl;
    // Caclulate chi2_3v
    double chi2_3v = osc_test->FCN(pars_3v_small); //

    // Create the data structures needed
    // TTree obs_tree("obs_tree", "Tree of deltachi2 obs values");
    const int NUM_dm2 = 60;
    const int NUM_ttt = 60;
    // double obs_map[NUM_ttt][NUM_dm2];
    TMatrixD obs_map(NUM_ttt,NUM_dm2);
    // Convert the grid into actual values

    TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -2, 1);
    TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -2, 0);
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
    TString fname = TString::Format("60x60-deltachi2_obs-xtoyuniverse-%i.root", ifile);
    TFile *rootfile = new TFile(fname, "recreate");
    obs_map.Write("obs_map");
    rootfile->Close();

  }// xthrow

  if (draw_confidence_map) {
    const int NUM_dm2 = 60;
    const int NUM_ttt = 60;
    int start_universe_index = 1;
    int end_universe_index = 200;

    for (int i = start_universe_index; i <= end_universe_index; i++) {
      // Load deltachi2_3v, deltachi2_4v, deltachi2_obs

      // Calculate pvalue4v

      // Calculate pvalue3v

      // Calculate cls w/ error protection for 0/0, inf, NaN values

      // Caclulate confidence value

      // Write to location in confidence map
    }
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
