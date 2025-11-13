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
// double profile_sin2_theta24(double pars, double profile_steps) {
//     cout << "Profile \n";
//     double test_sin2_theta24 = 0.001;
//     double profiled_sin2_theta24 = -1;
//     double chi2_min = 1e6;
//     // auto profile_forloop_start = chrono::high_resolution_clock::now();
//     // Grid scan
//     double step_size = sin24_range / profile_steps;
//     for(int i = 0; i < profile_steps; i++) {
//       test_sin2_theta24 = test_sin2_theta24 + (i * step_size);
//       if( sin2_2theta_mue_grid > test_sin2_theta24 ) {
//         continue; 		// Need to because I am inputing mixing ange directly
//       }	
//       pars[2] = test_sin2_theta24;// dm2, t14, t24, t34
//       double chi2 = osc_test->FCN( pars );
//       if( chi2 < chi2_min) {
//         chi2_min = chi2;
//         profiled_sin2_theta24 = test_sin2_theta24;
//       }
//     }
//     auto profile_forloop_end = chrono::high_resolution_clock::now();
//     cout << "Profiled sin^2 theta24: " << profiled_sin2_theta24 << endl;
//     cout << "Chi Square: " << chi2_min << endl;
//     auto total_profile_duration = chrono::duration_cast<chrono::seconds>(profile_forloop_end - profile_forloop_start);
//     // cout << "Total Profile Loop duration: " << total_profile_duration.count() << " microseconds" << endl;
//     return profiled_sin2_theta24;
//   }    
//     double calculate_ref_chi2(int data_or_asimov, double pars) {
//       cout << "Delta Chisquare Obs\n";
//       if (data_or_asimov == 0)
// 	{ // data
// 	  osc_test->Set_oscillation_pars(pars[0], pars[1], pars[2], pars[3]);
// 	  osc_test->Apply_oscillation();
// 	  osc_test->Set_apply_POT();// meas, CV, COV: all ready
// 	  osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
// 	}
//       // osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
    
//       osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
//       double chi2_testA = osc_test->FCN( pars_3v_small );
//       cout<<endl<<" ---> chi2_test3v "<<chi2_testA<<endl;

//       }
double calculate_prob_area(double pars, double ref_chi2);
double calculate_prob_delta_chi2(double pars_h0, double pars_h1, int num_toy);

double find_midpoint(double param_start, double param_end) {
  double midpoint = (param_end - param_start) / 2;
  return midpoint;
}

int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A story ..."<<endl<<endl;

  int ifile = 1;
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;

  // int it14 = 0;
  // int idm2 = 0;
  // int it24 = 0;
  double it14 = 0;
  double idm2 = 0;
  double it24 = 0;  
  
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

  if( 0 ) {
    
    val_dm2_41         = 1.19;
    val_sin2_2theta_14 = 0.8;
    val_sin2_theta_24  = 0.015;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
    
    double pars_4v_testA[4] = {2, 0.3, 0.1, 0};// dm2, sin_theta14, sin_theta24, other
    double chi2_testA = osc_test->FCN( pars_4v_testA );
    cout<<endl<<" ---> chi2_testA "<<chi2_testA<<endl<<endl;
    
    double pars_4v_testB[4] = {1.19, 0.8, 0.015};// dm2, sin_theta14, sin_theta24, other
    double chi2_testB = osc_test->FCN( pars_4v_testB );
    cout<<endl<<" ---> chi2_testB "<<chi2_testB<<endl<<endl;
    
    ///////
    
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    double pars_4v_testC[4] = {1.295, 0.936, 0};// dm2, sin_theta14, sin_theta24, other
    double chi2_testC = osc_test->FCN( pars_4v_testC );
    cout<<endl<<" ---> chi2_testC "<<chi2_testC<<endl<<endl;
        
  }
  
  /// one example: calculate a chi2 with pseudo experiment(s)
  /// one example: calcualte a chi2 with pseudo experiment(s)

  if( 0 ) {
  
    val_dm2_41         = 1.19;
    val_sin2_2theta_14 = 0.8;
    val_sin2_theta_24  = 0.015;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
    
    double pars_4v_testB[4] = {1.19, 0.8, 0.015};// dm2, sin_theta14, sin_theta24, other
    double chi2_testB = osc_test->FCN( pars_4v_testB );
    cout<<endl<<" ---> chi2_testB "<<chi2_testB<<endl<<endl;

    /////// 

    val_dm2_41         = 1.19;
    val_sin2_2theta_14 = 0.8;
    val_sin2_theta_24  = 0.015;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
        
    int NUM_TOYS = 10;
    osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments

    osc_test->Set_toy2fitdata( 1 );// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
    double chi2_testB1 = osc_test->FCN( pars_4v_testB );
    cout<<endl<<" ---> chi2_testB1 "<<chi2_testB1<<endl<<endl;
    
    osc_test->Set_toy2fitdata( 2 );// use the 2nd pseudo experiment as the "data", which will be compared with the "pred"
    double chi2_testB2 = osc_test->FCN( pars_4v_testB );
    cout<<endl<<" ---> chi2_testB2 "<<chi2_testB2<<endl<<endl;
    
    osc_test->Set_toy2fitdata( 9 );// use the 9th pseudo experiment as the "data", which will be compared with the "pred"
    double chi2_testB9 = osc_test->FCN( pars_4v_testB );
    cout<<endl<<" ---> chi2_testB9 "<<chi2_testB9<<endl<<endl;
    
  }

  /////// a test of the minimization
  /////// a test of the minimization

  if( 0 ) {
    cout << "a test of the minimization \n";
    //osc_test->Plot_user();  
    val_dm2_41         = 0;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.1;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    // osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"

    osc_test->Minimization_OscPars_FullCov(1.1, 0.7, 0.02, 0, "str_flag_fixpar");

  }

  if( 0 ) {
    cout << "Profile sinsquaretheta24";
    //osc_test->Plot_user();  
    val_dm2_41         = 0;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.1;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    // osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"

    osc_test->Minimization_OscPars_FullCov(55.38, 0.002330, 0.02, 0, "dm2_41_sin2_theta_14");

  }

  /////////////////////////////////////////////////////////// cross-check GoF
  /////////////////////////////////////////////////////////// cross-check GoF

  if( 0 ) {
      
    /// standard order
    val_dm2_41         = 0;
    val_sin2_2theta_14 = 0.2;
    val_sin2_theta_24  = 0.3;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_meas2fitdata();

    if( 1 ) {
      int rows = 25;
      int num_Y = 25;

      TMatrixD matrix_gof_trans( osc_test->Get_newworld_rows(), rows );// oldworld, newworld
      for( int ibin=1; ibin<=25; ibin++) matrix_gof_trans(26*7 + ibin-1, ibin-1) = 1; 

      TMatrixD matrix_gof_trans_T = matrix_gof_trans.T(); matrix_gof_trans.T(); 
      TMatrixD matrix_gof_pred = osc_test->matrix_tosc_eff_newworld_pred * matrix_gof_trans;
      TMatrixD matrix_gof_data = osc_test->matrix_tosc_fitdata_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_syst = matrix_gof_trans_T * (osc_test->matrix_tosc_eff_newworld_abs_syst_total) * matrix_gof_trans;     
      osc_test->Exe_Goodness_of_fit( num_Y, rows-num_Y, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 501 );
    }    

    if( 0 ) {
      int rows = 25;
      int num_Y = 25;

      TMatrixD matrix_gof_trans( osc_test->Get_newworld_rows(), rows );// oldworld, newworld
      for( int ibin=1; ibin<=25; ibin++) matrix_gof_trans(26*8 + ibin-1, ibin-1) = 1; 

      TMatrixD matrix_gof_trans_T = matrix_gof_trans.T(); matrix_gof_trans.T();
      TMatrixD matrix_gof_pred = osc_test->matrix_tosc_eff_newworld_pred * matrix_gof_trans;
      TMatrixD matrix_gof_data = osc_test->matrix_tosc_fitdata_newworld * matrix_gof_trans;
      TMatrixD matrix_gof_syst = matrix_gof_trans_T * (osc_test->matrix_tosc_eff_newworld_abs_syst_total) * matrix_gof_trans;
      osc_test->Exe_Goodness_of_fit( num_Y, rows-num_Y, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 502 );
    }
   
  }

  if ( 0 ) {
    cout << "Profile \n";
    double dm2_41_grid = 55.38;
    double sin2_2theta_mue_grid = 0.002330;
    double step_size = 0.0005; // profiling step size
    int steps = 2000; // profiling steps
   	  
    // Profiling data structures
    double profiled_sin2_theta24 = 0;
    double test_sin2_theta24 = 0;
    double chi2_min = pow(10,6);	// set very large to start

    
    // Grid scan

    for(int i = 0; i < steps; i++) {
      test_sin2_theta24 = i * step_size;
      if( sin2_2theta_mue_grid > test_sin2_theta24 ) {
	continue; 		// Need to because I am inputing mixing ange directly
      }	
      double pars_4nu[4] = {dm2_41_grid, sin2_2theta_mue_grid, test_sin2_theta24, 0};// dm2, t14, t24, t34
      double chi2 = osc_test->FCN( pars_4nu );
      if(chi2 < chi2_min) {
	chi2_min = chi2;
	profiled_sin2_theta24 = test_sin2_theta24;
      }
    }

    cout << "Profiled sin^2 theta24: " << profiled_sin2_theta24;
    cout << "Chi Square: " << chi2_min;

  }

  if ( 0 ) {
    cout << "4nu \n";
    // 4v Hypothesis

    //  3. For each 4\nu psuedoexperiment
    // - Calculate the \chi^{2}_{grid} by comparing the psuedoexperiment to the profiled grid point oscillation probabilbility equation

    // - Calculate the \chi^{2}_{3\nu{}} by comparing the psuedoexperiment to the 3\nu{} oscillation probability equation

    // - Calculate \Delta\chi^{2}_{grid toy} by taking the difference between \chi^{2}_{grid} and \chi^{2}_{3\nu{}}

    double dm2_41_grid = 55.38;
    double sin2_2theta_mue_grid = 0.002330;

    double pars_4v_grid[4] ={dm2_41_grid, sin2_2theta_mue_grid, 0.0045, 0};
    double pars_3v_small[4] = {0, 0.1, 0.1, 0};

    // Generating psuedo experiments
    osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
        
    const int NUM_TOYS = 500;
    osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments
    double array_4nu_chisquare_grid[NUM_TOYS];
    double array_4nu_chisquare_3nu[NUM_TOYS];
    double array_deltachisquare_grid_toy[NUM_TOYS];

    for (int i = 0; i < NUM_TOYS; i++) {
      osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
      array_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
      array_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
      array_deltachisquare_grid_toy[i] = array_4nu_chisquare_grid[i] - array_4nu_chisquare_3nu[i];
    }

    
    TFile *rootfile = new TFile("mendez_4nu_psuedo.root", "recreate");
    TTree *ttree = new TTree("ttree", "ttree");
    ttree->Branch( "chi2_grid",  &array_4nu_chisquare_grid, Form("array_4nu_chisquare_grid[%d]/D", NUM_TOYS));
    ttree->Branch( "chi2_3v",  &array_4nu_chisquare_3nu, Form("array_4nu_chisquare_3nu[%d]/D", NUM_TOYS));
    ttree->Branch( "delta_chi2_grid_toy",  &array_deltachisquare_grid_toy, Form("array_deltachisquare_grid_toy[%d]/D", NUM_TOYS));
    ttree->Fill();
    rootfile->Write();
    rootfile->Close();
    delete rootfile;
  }    

  if ( 0 ) {
    cout << "3nu \n";
    // 4v Hypothesis

    //  3. For each 4\nu psuedoexperiment
    // - Calculate the \chi^{2}_{grid} by comparing the psuedoexperiment to the profiled grid point oscillation probabilbility equation

    // - Calculate the \chi^{2}_{3\nu{}} by comparing the psuedoexperiment to the 3\nu{} oscillation probability equation

    // - Calculate \Delta\chi^{2}_{grid toy} by taking the difference between \chi^{2}_{grid} and \chi^{2}_{3\nu{}}

    double dm2_41_grid = 55.38;
    double sin2_2theta_mue_grid = 0.002330;

    double pars_4v_grid[4] ={dm2_41_grid, sin2_2theta_mue_grid, 0.0045, 0};
    double pars_3v_small[4] = {0, 0.1, 0.1, 0};

    // Generating psuedo experiments
    osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1], pars_3v_small[2], pars_3v_small[3]);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
        
    const int NUM_TOYS = 500;
    osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments
    double array_4nu_chisquare_grid[NUM_TOYS];
    double array_4nu_chisquare_3nu[NUM_TOYS];
    double array_deltachisquare_grid_toy[NUM_TOYS];
    
    for (int i = 0; i < NUM_TOYS; i++) {
      osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
      array_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
      array_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
      array_deltachisquare_grid_toy[i] = array_4nu_chisquare_grid[i] - array_4nu_chisquare_3nu[i];
    } 
    
    TFile *rootfile = new TFile("mendez_3nu_psuedo.root", "recreate");
    TTree *ttree = new TTree("ttree", "ttree");
    ttree->Branch( "chi2_grid",  &array_4nu_chisquare_grid, Form("array_4nu_chisquare_grid[%d]/D", NUM_TOYS));
    ttree->Branch( "chi2_3v",  &array_4nu_chisquare_3nu, Form("array_4nu_chisquare_3nu[%d]/D", NUM_TOYS));
    ttree->Branch( "delta_chi2_grid_toy",  &array_deltachisquare_grid_toy, Form("array_deltachisquare_grid_toy[%d]/D", NUM_TOYS));
    ttree->Fill();
    rootfile->Write();
    rootfile->Close();
    delete rootfile;
  }    
    
  if( 0 ) {
    double dm2_41_grid = 55.38;
    double sin2_2theta_mue_grid = 0.002330;

    double pars_4v_grid[4] ={dm2_41_grid, sin2_2theta_mue_grid, 0.0045, 0};
    double pars_3v_small[4] = {0, 0.1, 0.1, 0};
   
    val_dm2_41         = 1.19;
    val_sin2_2theta_14 = 0.8;
    val_sin2_theta_24  = 0.015;  
    osc_test->Set_oscillation_pars(0, 0.1, 0.1, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    // osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
    
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    double chi2_testA = osc_test->FCN( pars_3v_small );
    cout<<endl<<" ---> chi2_test3v "<<chi2_testA<<endl<<endl;
    
  
    ///////
    
    osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    double chi2_testC = osc_test->FCN( pars_4v_grid );
    cout<<endl<<" ---> chi2_test4v "<<chi2_testC<<endl<<endl;
        
  }

  if ( 1 ) { 			// Complete CLs in one run using command line functions as input
    auto total_cls_start = chrono::high_resolution_clock::now();
    //auto time_stop = chrono::high_resolution_clock::now();
    //auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);

    double dm2_41_grid = idm2;
    double sin2_2theta_mue_grid = it14;
    cout << "CLs for parameters (" << it14 << "," << idm2 << ")\n";

    double pars_4v_grid[4] ={dm2_41_grid, sin2_2theta_mue_grid, 0.0045, 0};
    double pars_3v_small[4] = {0, 0.1, 0.1, 0};
    const int NUM_TOYS = 500;
  
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
    if ( 1 ) {
      cout << "Profile \n";

      auto profile_forloop_start = chrono::high_resolution_clock::now();
      // Grid scan
      for(int i = 0; i < steps; i++) {
	test_sin2_theta24 = i * step_size;
	if( sin2_2theta_mue_grid > test_sin2_theta24 ) {
	  continue; 		// Need to because I am inputing mixing ange directly
	}	
	double pars_4nu[4] = {dm2_41_grid, sin2_2theta_mue_grid, test_sin2_theta24, 0};// dm2, t14, t24, t34
	double chi2 = osc_test->FCN( pars_4nu );
	if(chi2 < chi2_min) {
	  chi2_min = chi2;
	  profiled_sin2_theta24 = test_sin2_theta24;
	}
      }
      auto profile_forloop_end = chrono::high_resolution_clock::now();
      cout << "Profiled sin^2 theta24: " << profiled_sin2_theta24 << endl;
      cout << "Chi Square: " << chi2_min << endl;
      auto total_profile_duration = chrono::duration_cast<chrono::seconds>(profile_forloop_end - profile_forloop_start);
      pars_4v_grid[2] = profiled_sin2_theta24;
      cout << "Total Profile Loop duration: " << total_profile_duration.count() << " microseconds" << endl;
    }

  
    if ( 1 ) {
      cout << "4nu \n";
      // 4v Hypothesis

      //  3. For each 4\nu psuedoexperiment
      // - Calculate the \chi^{2}_{grid} by comparing the psuedoexperiment to the profiled grid point oscillation probabilbility equation

      // - Calculate the \chi^{2}_{3\nu{}} by comparing the psuedoexperiment to the 3\nu{} oscillation probability equation

      // - Calculate \Delta\chi^{2}_{grid toy} by taking the difference between \chi^{2}_{grid} and \chi^{2}_{3\nu{}}


      // Generating psuedo experiments
      osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
      osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
        
    
      osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments


      auto delta_chi2_loop_start = chrono::high_resolution_clock::now();
      for (int i = 0; i < NUM_TOYS; i++) {
	osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
	array_4vhyp_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
	array_4vhyp_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
	array_4vhyp_deltachisquare_grid_toy[i] = array_4vhyp_4nu_chisquare_grid[i] - array_4vhyp_4nu_chisquare_3nu[i];
      } 
      auto delta_chi2_loop_end = chrono::high_resolution_clock::now();
      // TFile *rootfile = new TFile("mendez_4nu_psuedo.root", "recreate");
      // TTree *ttree = new TTree("ttree", "ttree");
      // ttree->Branch( "chi2_grid",  &array_4nu_chisquare_grid, Form("array_4nu_chisquare_grid[%d]/D", NUM_TOYS));
      // ttree->Branch( "chi2_3v",  &array_4nu_chisquare_3nu, Form("array_4nu_chisquare_3nu[%d]/D", NUM_TOYS));
      // ttree->Branch( "delta_chi2_grid_toy",  &array_deltachisquare_grid_toy, Form("array_deltachisquare_grid_toy[%d]/D", NUM_TOYS));
      // ttree->Fill();
      // rootfile->Write();
      // rootfile->Close();
      // delete rootfile;
      auto total_delta_chi2_loop_duration = chrono::duration_cast<chrono::seconds>(delta_chi2_loop_end - delta_chi2_loop_start);
      cout << "Delta chi2 Loop duration: " << total_delta_chi2_loop_duration.count() << " microseconds" << endl;
    }    

    if ( 1 ) {
      cout << "3nu \n";
      // 4v Hypothesis

      //  3. For each 4\nu psuedoexperiment
      // - Calculate the \chi^{2}_{grid} by comparing the psuedoexperiment to the profiled grid point oscillation probabilbility equation

      // - Calculate the \chi^{2}_{3\nu{}} by comparing the psuedoexperiment to the 3\nu{} oscillation probability equation

      // - Calculate \Delta\chi^{2}_{grid toy} by taking the difference between \chi^{2}_{grid} and \chi^{2}_{3\nu{}}

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

      // TFile *rootfile = new TFile("mendez_3nu_psuedo.root", "recreate");
      // TTree *ttree = new TTree("ttree", "ttree");
      // ttree->Branch( "chi2_grid",  &array_4nu_chisquare_grid, Form("array_4nu_chisquare_grid[%d]/D", NUM_TOYS));
      // ttree->Branch( "chi2_3v",  &array_4nu_chisquare_3nu, Form("array_4nu_chisquare_3nu[%d]/D", NUM_TOYS));
      // ttree->Branch( "delta_chi2_grid_toy",  &array_deltachisquare_grid_toy, Form("array_deltachisquare_grid_toy[%d]/D", NUM_TOYS));
      // ttree->Fill();
      // rootfile->Write();
      // rootfile->Close();
      // delete rootfile;
      }}    

    auto total_cls_end = chrono::high_resolution_clock::now();
    if( 1 ) {

      cout << "Delta Chisquare Obs\n";
      osc_test->Set_oscillation_pars(0, 0.1, 0.1, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
      // osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
    
      osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
      double chi2_testA = osc_test->FCN( pars_3v_small );
      cout<<endl<<" ---> chi2_test3v "<<chi2_testA<<endl;
    
      ///////
    
      osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
      double chi2_testC = osc_test->FCN( pars_4v_grid );
      cout<<endl<<" ---> chi2_test4v "<<chi2_testC<<endl;
      ref = chi2_testC - chi2_testA;
      cout << "Delta Chi2: " << ref << endl;
        
    }

    auto total_cls_duration = chrono::duration_cast<chrono::seconds>(total_cls_end - total_cls_start);
    cout << "CLs duration: "<< total_cls_duration.count() << " microseconds" << endl;

    if (1) {
      // Create histograms for 3v/4v delta_chi2_grid_toy

      double min3 = 1e30, max3 = -1e30;
      double min4 = 1e30, max4 = -1e30;

      for (int j = 0; j < NUM_TOYS; ++j) {
	min4 = std::min(min4, array_4vhyp_deltachisquare_grid_toy[j]); max4 = std::max(max4, array_4vhyp_deltachisquare_grid_toy[j]);
	min3 = std::min(min3, array_3vhyp_deltachisquare_grid_toy[j]); max3 = std::max(max3, array_3vhyp_deltachisquare_grid_toy[j]);
      }
      TH1D *h1 = new TH1D("h1", "p3v;Delta chi2;Entries", 100, min3, max3);
      TH1D *h2 = new TH1D("h2", "p4v;Delta chi2;Entries", 100, min4, max4);
      for (int j = 0; j < NUM_TOYS; ++j) {
	// cout << array_3vhyp_deltachisquare_grid_toy[j] << ":" << array_4vhyp_deltachisquare_grid_toy[j] <<"\n";
	h1->Fill(array_3vhyp_deltachisquare_grid_toy[j]);
	h2->Fill(array_4vhyp_deltachisquare_grid_toy[j]);
      }
      
      h1->Scale(1.0 / h1->Integral());
      h2->Scale(1.0 / h2->Integral());
      // Calculate CLs from histogram bin
      int bin4v = h2->FindBin(ref);
      int bin3v = h1->FindBin(ref);
      cout << bin3v << "\n" << bin4v << "\n";
      double p3v = h1->Integral(bin3v, h1->GetNbinsX());
      double p4v = h2->Integral(bin4v, h2->GetNbinsX());
      cout <<"DEBUG\n" << p3v << "\n" << p4v <<"\n" << h1->FindBin(ref) << "\n" << h2->FindBin(ref) << "\n";
      double val_CLs = p4v / p3v;
      cout << "CLs: " << val_CLs << "\n";

      // Save Everything I need to files
      TString fname = TString::Format("./mendez_osc_analysis_validation_curve_CLs_dm2_%.4e_mue_%.4e.root", idm2, it14);
      TFile *rootfile = new TFile(fname, "recreate");
      // Create objects for cpp variables
      TParameter<double> deltachi2_ref("deltachi2_ref", ref);
      TParameter<int> psuedos("num_toys", NUM_TOYS);
      TParameter<double> cls("CLs", val_CLs);
      // Create TArray objects for all cpp arrays
      // 4v stuctures
      TVectorD arr4d(NUM_TOYS, array_4vhyp_deltachisquare_grid_toy);
      // TArrayD arr4g(NUM_TOYS, array_4vhyp_4nu_chisquare_grid);
      // TArrayD arr43(NUM_TOYS, array_4vhyp_4nu_chisquare_3nu);
      TVectorD arr3d(NUM_TOYS, array_3vhyp_deltachisquare_grid_toy);
      // TArrayD arr3g(NUM_TOYS, array_3vhyp_4nu_chisquare_grid);
      // TArrayD arr33(NUM_TOYS, array_3vhyp_4nu_chisquare_3nu);

      deltachi2_ref.Write();
      psuedos.Write();
      cls.Write();
      arr4d.Write("4v_deltachi2");
      arr3d.Write("3v_deltachi2");
      h1->Write();
      h2->Write();
      rootfile->Write();
      rootfile->Close();
      delete rootfile;

    }
  }
  if (0) {			// Binary Grid search algorithm
    cout << "Start binary search algo\n";
    auto time_start = chrono::high_resolution_clock::now();
    double counter = 0;
    double tolerance = 0.005;
    double test_point = -1;
    double val_CLs = -1;
    double previous_CLs = 0;
    double target = 0.05;
    double param_start = it14;
    double param_end = it24;
    while (val_CLs < (target - tolerance) || val_CLs > (target + tolerance)) {
      test_point = find_midpoint(param_start, param_end);
      double dm2_41_grid = idm2;
      double sin2_2theta_mue_grid = test_point;
      cout << "CLs for parameters (" << test_point << "," << idm2 << ")\n";

      double pars_4v_grid[4] ={dm2_41_grid, sin2_2theta_mue_grid, 0.0045, 0};
      double pars_3v_small[4] = {0, 0.1, 0.1, 0};
      const int NUM_TOYS = 500;
  
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
      if ( 1 ) {
	cout << "Profile \n";
	// auto profile_forloop_start = chrono::high_resolution_clock::now();
	// Grid scan
	for(int i = 0; i < steps; i++) {
	  test_sin2_theta24 = i * step_size;
	  if( sin2_2theta_mue_grid > test_sin2_theta24 ) {
	    continue; 		// Need to because I am inputing mixing ange directly
	  }	
	  double pars_4nu[4] = {dm2_41_grid, sin2_2theta_mue_grid, test_sin2_theta24, 0};// dm2, t14, t24, t34
	  double chi2 = osc_test->FCN( pars_4nu );
	  if(chi2 < chi2_min) {
	    chi2_min = chi2;
	    profiled_sin2_theta24 = test_sin2_theta24;
	  }
	}
	// auto profile_forloop_end = chrono::high_resolution_clock::now();
	cout << "Profiled sin^2 theta24: " << profiled_sin2_theta24 << endl;
	cout << "Chi Square: " << chi2_min << endl;
	// auto total_profile_duration = chrono::duration_cast<chrono::seconds>(profile_forloop_end - profile_forloop_start);
	pars_4v_grid[2] = profiled_sin2_theta24;
	// cout << "Total Profile Loop duration: " << total_profile_duration.count() << " microseconds" << endl;
      
	if ( 1 ) {
	  cout << "4nu \n";
	  // Generating psuedo experiments
	  osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);  
	  osc_test->Apply_oscillation();
	  osc_test->Set_apply_POT();// meas, CV, COV: all ready
	  //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
	  osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
	  osc_test->Set_toy_variations( NUM_TOYS );// produce NUM_TOYS pseudo experiments
	  // auto delta_chi2_loop_start = chrono::high_resolution_clock::now();
	  for (int i = 0; i < NUM_TOYS; i++) {
	    osc_test->Set_toy2fitdata(i + 1);// use the 1st pseudo experiment as the "data", which will be compared with the "pred"
	    array_4vhyp_4nu_chisquare_grid[i] = osc_test->FCN(pars_4v_grid );
	    array_4vhyp_4nu_chisquare_3nu[i] = osc_test->FCN(pars_3v_small );
	    array_4vhyp_deltachisquare_grid_toy[i] = array_4vhyp_4nu_chisquare_grid[i] - array_4vhyp_4nu_chisquare_3nu[i];
	  }
	}  
	if ( 1 ) {
	  cout << "3nu \n";
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
	}
	// auto total_cls_end = chrono::high_resolution_clock::now();
	if( 1 ) {
	  cout << "Delta Chisquare Obs\n";
	  osc_test->Set_oscillation_pars(0, 0.1, 0.1, 0);  
	  osc_test->Apply_oscillation();
	  osc_test->Set_apply_POT();// meas, CV, COV: all ready
	  osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
	  // osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
	  osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
	  double chi2_testA = osc_test->FCN( pars_3v_small );
	  cout<<endl<<" ---> chi2_test3v "<<chi2_testA<<endl;
	  osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
	  double chi2_testC = osc_test->FCN( pars_4v_grid );
	  cout<<endl<<" ---> chi2_test4v "<<chi2_testC<<endl;
	  ref = chi2_testC - chi2_testA;
	  cout << "Delta Chi2: " << ref << endl;
	}
	// auto total_cls_duration = chrono::duration_cast<chrono::seconds>(total_cls_end - total_cls_start);
	// cout << "CLs duration: "<< total_cls_duration.count() << " microseconds" << endl;
	if (1) {
	  // Create histograms for 3v/4v delta_chi2_grid_toy
	  double min3 = 1e30, max3 = -1e30;
	  double min4 = 1e30, max4 = -1e30;
	  for (int j = 0; j < NUM_TOYS; ++j) {
	    min4 = std::min(min4, array_4vhyp_deltachisquare_grid_toy[j]); max4 = std::max(max4, array_4vhyp_deltachisquare_grid_toy[j]);
	    min3 = std::min(min3, array_3vhyp_deltachisquare_grid_toy[j]); max3 = std::max(max3, array_3vhyp_deltachisquare_grid_toy[j]);
	  }
	  TH1D *h1 = new TH1D("h1", "p3v;Delta chi2;Entries", 100, min3, max3);
	  TH1D *h2 = new TH1D("h2", "p4v;Delta chi2;Entries", 100, min4, max4);
	  for (int j = 0; j < NUM_TOYS; ++j) {
	    // cout << array_3vhyp_deltachisquare_grid_toy[j] << ":" << array_4vhyp_deltachisquare_grid_toy[j] <<"\n";
	    h1->Fill(array_3vhyp_deltachisquare_grid_toy[j]);
	    h2->Fill(array_4vhyp_deltachisquare_grid_toy[j]);
	  }
	  h1->Scale(1.0 / h1->Integral());
	  h2->Scale(1.0 / h2->Integral());
	  // Calculate CLs from histogram bin
	  int bin4v = h2->FindBin(ref);
	  int bin3v = h1->FindBin(ref);
	  cout << bin3v << "\n" << bin4v << "\n";
	  double p3v = h1->Integral(bin3v, h1->GetNbinsX());
	  double p4v = h2->Integral(bin4v, h2->GetNbinsX());
	  cout <<"DEBUG\n" << p3v << "\n" << p4v <<"\n" << h1->FindBin(ref) << "\n" << h2->FindBin(ref) << "\n";
	  val_CLs = p4v / p3v;
	  cout << "CLs: " << val_CLs << "\n";
	}
      }
      // delta = abs(previous_CLs - val_CLs);
      if (val_CLs > target) {
	param_end = test_point;
      }
      else if (val_CLs < target) {
	param_start = test_point;
      }
      counter++;
      cout << "Iteration  " << counter << "\n";
      cout << "previous CLs: " << previous_CLs << "\n";
      cout << "CLs: " << val_CLs << "\n";
      previous_CLs = val_CLs;
    }// while
    double exclusion_curve_parameter = test_point;
    double exclusion_CLs = val_CLs;
    auto time_stop = chrono::high_resolution_clock::now();
    auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
    cout << "Binary Search Algo Duration: " << time_duration.count() << " microseconds" << "\n";
    cout << "exclusion_curve (sin2_2theta_mue,dm2) \n";
    cout << exclusion_curve_parameter << "," << idm2 << "\n";
    cout << "exclusion_CLs: " << exclusion_CLs << "\n" ;
  }


  // if (0) { // Iterative CLs calculator
  //   // Constants to set grid granularity
  //   const int scan_steps = 10;
  //   double profile_step_size = 0.02;
  //   int profile_steps = 100;
  //   const int ntoys = 500;

  //   // Set the CLs target
  //   double cls_target = 0.05;

  //   // Convert command line parameters to actual value range
  //   double dm2_41_grid = idm2;
  //   double sin2_2theta_mue_start = it14;
  //   double sin2_2theta_mue_end = it24;
  //   double sin2_2theta_mue_list[scan_steps];
  //   sin2_2theta_mue_list[0] = sin2_2theta_mue_start;
  //   sin2_2theta_mue_list[scan_steps - 1] = sin2_2theta_mue_end;
  //   double scan_diff = sin2_2theta_mue_end - sin2_2theta_mue_start;
  //   for (int i = 1; i<= scan_steps; i++) {
  //     sin2_2theta_mue_list[i] = sin2_2theta_mue_list[i-1] + scan_diff;
  //   }

  //   // Create the data structures to hold the previous calculation
  //   // The memory will only hold 2 sets of data max
  //   double previous_pars_4v_grid[4] = { -1, -1, -1, -1};
  //   double previous_profiled_sin2_theta_12_grid = -1;
  //   double previous_crit_delta_chi2 = -1;
  //   double previous_delta_chisquare_4v_dist[ntoys];
  //   double previous_delta_chisquare_3v_dist[ntoys];
  //   double previous_clsb = -1;
  //   double previous_cls = -1;
  //   double previous_CLs_val = -1;
    
  //   // Main Loop
  //   for (int i= 0; i< scan_steps; i++) {
  //     double pars_4v_grid[4] = {dm2, sin2_theta_mue[i], 0, 0};
  //     	// Profile
  // 	cout << "Profile \n";
  // 	// auto profile_forloop_start = chrono::high_resolution_clock::now();
  // 	// Grid scan
  // 	for(int j = 0; j < profile_steps; j++) {
  // 	  double test_sin2_theta24 = j * step_size;
  // 	  if( sin2_2theta_mue_list[i] > test_sin2_theta24 ) {
  // 	    continue; 		// Need to because I am inputing mixing ange directly
  // 	  }	
  // 	  double pars_4nu[4] = {dm2_41_grid, sin2_2theta_mue_list[i], test_sin2_theta24, 0};// dm2, t14, t24, t34
  // 	  double chi2 = osc_test->FCN( pars_4nu );
  // 	  if(chi2 < chi2_min) {
  // 	    chi2_min = chi2;
  // 	    profiled_sin2_theta24 = test_sin2_theta24;
  // 	  }
  // 	}
  // 	// auto profile_forloop_end = chrono::high_resolution_clock::now();
  // 	cout << "Profiled sin^2 theta24: " << profiled_sin2_theta24 << endl;
  // 	cout << "Chi Square: " << chi2_min << endl;

      
  //   }
  // }


  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////

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
