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
//milliseconds, minutes

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

  osc_test->scaleF_POT_BNB  = scaleF_POT_BNB;
  osc_test->scaleF_POT_NuMI = scaleF_POT_NuMI;
     
  ///////

  osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;
  
  ///////

  osc_test->flag_NuMI_nueCC_from_intnue       = Configure_Osc::flag_NuMI_nueCC_from_intnue;
  osc_test->flag_NuMI_nueCC_from_overlaynumu  = Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
  osc_test->flag_NuMI_nueCC_from_appnue       = Configure_Osc::flag_NuMI_nueCC_from_appnue;
  osc_test->flag_NuMI_nueCC_from_appnumu      = Configure_Osc::flag_NuMI_nueCC_from_appnumu;
  osc_test->flag_NuMI_nueCC_from_overlaynueNC = Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
  osc_test->flag_NuMI_nueCC_from_overlaynumuNC= Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;
  
  osc_test->flag_NuMI_numuCC_from_overlaynumu   = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue    = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue        = Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu       = Configure_Osc::flag_NuMI_numuCC_from_appnumu;
  osc_test->flag_NuMI_numuCC_from_overlaynumuNC = Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
  osc_test->flag_NuMI_numuCC_from_overlaynueNC  = Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;  
  
  osc_test->flag_NuMI_CCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue       = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  osc_test->flag_NuMI_CCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_CCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;
  
  osc_test->flag_NuMI_NCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue       = Configure_Osc::flag_NuMI_NCpi0_from_appnue;
  osc_test->flag_NuMI_NCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_NCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;


  ///////
  
  osc_test->flag_BNB_nueCC_from_intnue       = Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu  = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue       = Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu      = Configure_Osc::flag_BNB_nueCC_from_appnumu;
  osc_test->flag_BNB_nueCC_from_overlaynueNC = Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
  osc_test->flag_BNB_nueCC_from_overlaynumuNC= Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;
  
  osc_test->flag_BNB_numuCC_from_overlaynumu   = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue    = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue        = Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu       = Configure_Osc::flag_BNB_numuCC_from_appnumu;
  osc_test->flag_BNB_numuCC_from_overlaynumuNC = Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
  osc_test->flag_BNB_numuCC_from_overlaynueNC  = Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;  
  
  osc_test->flag_BNB_CCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue       = Configure_Osc::flag_BNB_CCpi0_from_appnue;
  osc_test->flag_BNB_CCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_CCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;
  
  osc_test->flag_BNB_NCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue       = Configure_Osc::flag_BNB_NCpi0_from_appnue;
  osc_test->flag_BNB_NCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_NCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;
  
  /////// set only one time
  
  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
                               Configure_Osc::default_dirtadd_file,
                               Configure_Osc::default_mcstat_file,
                               Configure_Osc::default_fluxXs_dir,
                               Configure_Osc::default_detector_dir);
  
  osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);
  
  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)
  
  double val_dm2_41         = 7.3;
  double val_sin2_2theta_14 = 0.36;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;
  
  /////////////////////////////////////////////////////////////////////////////////

  /// standard order
  // val_dm2_41         = 0;
  // val_sin2_2theta_14 = 0.1;
  // val_sin2_theta_24  = 0.1;
  
  // osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
  // osc_test->Apply_oscillation();  
  // osc_test->Set_apply_POT();// meas, CV, COV: all ready
  
  // osc_test->Set_meas2fitdata();
  
  // osc_test->Set_asimov2fitdata();

  // osc_test->FCN(const double *par);
 
  // osc_test->Minimization_OscPars_FullCov(dm2, t14, t24, t34, "str_flag_fixpar");

  // osc_test->Set_toy_variations(int num_toys);
  // osc_test->Set_toy2fitdata(int itoy)

  //
  //
  // The statistic test is implemented by the comination of the above functions
  //
  //
  
  ///////////////////////////////////////////////////////////////////////////////// one example: chi2 calculation

  if( 0 ) {
    val_dm2_41         = 0;// no-oscillation
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.1;
    val_sin2_theta_34  = 0;
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);// dm2, t14, t24, t34
    osc_test->Apply_oscillation();// Oscillation formula: Prob_oscillaion(double Etrue, double baseline, int strflag_osc);
                                  // the meaning of dm2, t14, t24, t34 are decided by the definitions in Prob_oscillaion(...)
    osc_test->Set_apply_POT();// meas, CV, COV ---> all ready

    //osc_test->Set_meas2fitdata();// set the measured data to be used in the chi2 or fitting
    osc_test->Set_asimov2fitdata();// set the Asimov data to be used in the chi2 or fitting,
                                   // which is produced by the above: Set_oscillation_pars(...), Apply_oscillation(), Set_apply_POT()

    double pars_4v[4] = {0, 0.5, 0.5, 0};// dm2, t14, t24, t34
    double chi2_result = osc_test->FCN( pars_4v );
    cout<<endl<<" ---> chi2 result: "<<chi2_result<<endl<<endl;;
  }

  ///////////////////////////////////////////////////////////////////////////////// one example: minimization

  if( 0 ) {
    val_dm2_41         = 2;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.1;
    val_sin2_theta_34  = 0;
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);// dm2, t14, t24, t34
    osc_test->Apply_oscillation();// Oscillation formula: Prob_oscillaion(double Etrue, double baseline, int strflag_osc);
    osc_test->Set_apply_POT();// meas, CV, COV ---> all ready
    
    //osc_test->Set_meas2fitdata();// set the measured data to be used in the chi2 or fitting
    osc_test->Set_asimov2fitdata();// set the Asimov data to be used in the chi2 or fitting,
                                   // which is produced by the above: Set_oscillation_pars(...), Apply_oscillation(), Set_apply_POT()

    double initial_dm2 = 1.9;
    double initial_t14 = 0.09;
    double initial_t24 = 0.09;
    osc_test->Minimization_OscPars_FullCov(initial_dm2, initial_t14, initial_t24, 0, "");

    // Minimization_OscPars_FullCov(double init_dm2_41, double init_sin2_2theta_14, double init_sin2_theta_24, double init_sin2_theta_34, TString roostr_flag_fixpar);
    // the meaning of oscillation paramters are decided by the definitions in Prob_oscillaion(...)
    // roostr_flag_fixpar: fix oscllation par(s)?
    //       if roostr_flag_fixpar contains "dm2", "t14", or "t24", then it will fix the par(s) at the initial value(s) 
    //       ---> roostr_flag_fixpar = "t14" : fix the value at init_sin2_2theta_14
    //       ---> roostr_flag_fixpar = "t14_dm2" : fix the two pars at init_sin2_2theta_14 and initial_dm2
    //       ---> roostr_flag_fixpar = "dm2_t14_t24": fix all three pars
  }
  
  ///////////////////////////////////////////////////////////////////////////////// one example: produce toys

  if( 0 ) {
    val_dm2_41         = 0;// no-oscillation
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.1;
    val_sin2_theta_34  = 0;
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);// dm2, t14, t24, t34
    osc_test->Apply_oscillation();// Oscillation formula: Prob_oscillaion(double Etrue, double baseline, int strflag_osc);
    osc_test->Set_apply_POT();// meas, CV, COV ---> all ready

    //osc_test->Set_meas2fitdata();// set the measured data to be used in the chi2 or fitting
    //osc_test->Set_asimov2fitdata();// set the Asimov data to be used in the chi2 or fitting,
                                     // which is produced by the above: Set_oscillation_pars(...), Apply_oscillation(), Set_apply_POT()

    osc_test->Set_toy_variations( 10 );// produce 10 toy samples use above oscillation pars: Set_toy_variations(int num_toys), num_toys >=1
    osc_test->Set_toy2fitdata( 3 );// set the 3rd toy sample to be used in the chi2 or fitting: Set_toy2fitdata(int itoy), itoy <=num_toys

    double pars_4v[4] = {0, 0.5, 0.5, 0};// dm2, t14, t24, t34
    double chi2_result = osc_test->FCN( pars_4v );
    cout<<endl<<" ---> chi2 result: "<<chi2_result<<endl<<endl;;
  }

  ///////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

  cout<<endl;
  cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f, dm2(-idm2) %d, theta14(-it14) %d, theta24(-it24) %d",
			display, ifile, osc_test->scaleF_POT_BNB, osc_test->scaleF_POT_NuMI, idm2, it14, it24)<<endl;

  cout<<endl;
  cout<<TString::Format(" ---> Files/Directories: ")<<endl;
  cout<<TString::Format("      ---> default_cv_file       %-50s", Configure_Osc::default_cv_file.Data())<<endl;
  cout<<TString::Format("      ---> default_dirtadd_file  %-50s", Configure_Osc::default_dirtadd_file.Data())<<endl;
  cout<<TString::Format("      ---> default_mcstat_file   %-50s", Configure_Osc::default_mcstat_file.Data())<<endl;
  cout<<TString::Format("      ---> default_fluxXs_dir    %-50s", Configure_Osc::default_fluxXs_dir.Data())<<endl;
  cout<<TString::Format("      ---> default_detector_dir  %-50s", Configure_Osc::default_detector_dir.Data())<<endl;
  cout<<TString::Format("      ---> default_eventlist_dir %-50s", Configure_Osc::default_eventlist_dir.Data())<<endl;
  
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
