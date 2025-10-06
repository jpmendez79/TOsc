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

//#include <chrono> // timer
//auto time_start = chrono::high_resolution_clock::now();
//auto time_stop = chrono::high_resolution_clock::now();
//auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
//cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
//microseconds, milliseconds, seconds, minutes

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
  // double tosc_sin2_2theta_14
  // double tosc_sin2_theta_24
  // double tosc_sin2_theta_34
  //
  // which are defined in "TOsc.h"
  //
  ///////////////////////////////////////////////////////////

  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)

  double val_dm2_41         = 0;
  double val_sin2_2theta_14 = 0.36;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;

  /// standard order
  /// standard order
  
  val_dm2_41         = 0;
  val_sin2_2theta_14 = 0.236;
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
    double pars_4v_testC[4] = {1.19, 0.8, 0.015};// dm2, sin_theta14, sin_theta24, other
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
    //osc_test->Plot_user();  
    val_dm2_41         = 1.19;
    val_sin2_2theta_14 = 0.8;
    val_sin2_theta_24  = 0.015;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    //osc_test->Set_meas2fitdata();// set the real measurement as the "data", which will be compared with the "pred"
    osc_test->Set_asimov2fitdata();// set the "asimov toy sample" as the "data", which will be compared with the "pred"
  
    osc_test->Minimization_OscPars_FullCov(1.1, 0.7, 0.02, 0, "str_flag_fixpar");
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

  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  

  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

  cout<<endl;
  cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f, theta14(-it14) %d, dm2(-idm2) %d, theta24(-it24) %d",
			display, ifile, osc_test->tosc_scaleF_POT_BNB, osc_test->tosc_scaleF_POT_NuMI, it14, idm2, it24)<<endl;
  cout<<endl;

  cout<<" ---> flag_apply_oscillation_BNB  "<< osc_test->flag_apply_oscillation_BNB<<endl;
  cout<<" ---> flag_apply_oscillation_NuMI "<< osc_test->flag_apply_oscillation_NuMI<<endl;
  cout<<endl;

  cout<<" ---> flag_goodness_of_fit_CNP    "<< osc_test->flag_goodness_of_fit_CNP<<endl;
  cout<<endl;

  cout<<TString::Format(" ---> flag_syst_dirt    %d", osc_test->flag_syst_dirt)<<endl;
  cout<<TString::Format(" ---> flag_syst_mcstat  %d", osc_test->flag_syst_mcstat)<<endl;
  cout<<TString::Format(" ---> flag_syst_flux    %d", osc_test->flag_syst_flux)<<endl;
  cout<<TString::Format(" ---> flag_syst_geant   %d", osc_test->flag_syst_geant)<<endl;
  cout<<TString::Format(" ---> flag_syst_Xs      %d", osc_test->flag_syst_Xs)<<endl;
  cout<<TString::Format(" ---> flag_syst_det     %d", osc_test->flag_syst_det)<<endl;

  cout<<endl;
  cout<<" ---> Finished sucessfully"<<endl;
  
  cout<<endl;
  if( display ) {
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<endl;
    theApp.Run();
  }
  
  return 0;
}
