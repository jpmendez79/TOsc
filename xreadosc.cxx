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

#include <chrono> // timer
//auto time_start = chrono::high_resolution_clock::now();
//auto time_stop = chrono::high_resolution_clock::now();
//auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
//cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
//microseconds, milliseconds, seconds, minutes

////////////////////////////////// SUB /////////////////////////////////////////////

void get_CL_curve(TH2 *h2_CL_input, TGraph *gh_CL_curve, int flag_index);

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

  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
			       Configure_Osc::default_dirtadd_file,
			       Configure_Osc::default_mcstat_file,
			       Configure_Osc::default_fluxXs_dir,
			       Configure_Osc::default_detector_dir);

  osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);

  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)

  double val_dm2_41         = 0;
  double val_sin2_2theta_14 = 0.36;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;

  /// standard order
  val_dm2_41         = 0;
  val_sin2_2theta_14 = 0.236;
  val_sin2_theta_24  = 0.3;
  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  osc_test->Set_meas2fitdata();
  //osc_test->Set_asimov2fitdata();

  ///////
  //osc_test->Plot_user();
  //osc_test->Minimization_OscPars_FullCov(1.1885, 0.88760418, 0.0005, 0, "str_flag_fixpar");

  ///////
  // if( 0 ) {
  //   int rows = 26*4;
  //   int num_Y = 26*2;

  //   TMatrixD matrix_gof_trans( osc_test->Get_newworld_rows(), rows );// oldworld, newworld
  //   for( int ibin=1; ibin<=52; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
  //   for( int ibin=1; ibin<=52; ibin++) matrix_gof_trans(26*7 + ibin-1, 52 + ibin-1) = 1;

  //   TMatrixD matrix_gof_trans_T = matrix_gof_trans.T(); matrix_gof_trans.T();
  //   TMatrixD matrix_gof_pred = osc_test->matrix_tosc_eff_newworld_pred * matrix_gof_trans;
  //   TMatrixD matrix_gof_data = osc_test->matrix_tosc_fitdata_newworld * matrix_gof_trans;
  //   TMatrixD matrix_gof_syst = matrix_gof_trans_T * (osc_test->matrix_tosc_eff_newworld_abs_syst_total) * matrix_gof_trans;
  //   osc_test->Exe_Goodness_of_fit( num_Y, rows-num_Y, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 501 );
  // }

  /////////////////////////////////////////////////////////// cross-check GoF
  /////////////////////////////////////////////////////////// cross-check GoF

  // if( 0 ) {
  //   int rows  = 26*14;
  //   int num_Y = 26*14;

  //   TMatrixD matrix_gof_pred = osc_test->matrix_tosc_eff_newworld_pred;
  //   TMatrixD matrix_gof_data = osc_test->matrix_tosc_fitdata_newworld;
  //   TMatrixD matrix_gof_syst = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  //   osc_test->Exe_Goodness_of_fit( num_Y, rows-num_Y, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 501 );
  // }


  // if( 0 ) {

  //   /// standard order
  //   val_dm2_41         = 0;
  //   val_sin2_2theta_14 = 0.2;
  //   val_sin2_theta_24  = 0.3;
  //   osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);
  //   osc_test->Apply_oscillation();
  //   osc_test->Set_apply_POT();// meas, CV, COV: all ready
  //   osc_test->Set_meas2fitdata();

  //   if( 1 ) {
  //     int rows = 25;
  //     int num_Y = 25;

  //     TMatrixD matrix_gof_trans( osc_test->Get_newworld_rows(), rows );// oldworld, newworld
  //     for( int ibin=1; ibin<=25; ibin++) matrix_gof_trans(26*3 + ibin-1, ibin-1) = 1;

  //     TMatrixD matrix_gof_trans_T = matrix_gof_trans.T(); matrix_gof_trans.T();
  //     TMatrixD matrix_gof_pred = osc_test->matrix_tosc_eff_newworld_pred * matrix_gof_trans;
  //     TMatrixD matrix_gof_data = osc_test->matrix_tosc_fitdata_newworld * matrix_gof_trans;
  //     TMatrixD matrix_gof_syst = matrix_gof_trans_T * (osc_test->matrix_tosc_eff_newworld_abs_syst_total) * matrix_gof_trans;
  //     osc_test->Exe_Goodness_of_fit( num_Y, rows-num_Y, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 504 );
  //   }

  // }

  /////////////////////////////////////////////////////////// Median sensitivity from frequentist CLs
  /////////////////////////////////////////////////////////// Median sensitivity from frequentist CLs

  if( 1 ) {

    cout<<endl;
    cout<<" ---> Median sensitivity from freqentist CLs "<<endl;
    cout<<endl;

    bool flag_step1_produce_toys = 0;
    bool flag_step3_CLs          = 1;

    /////// xxx (1)

    int bins_all = 26*7;
    int bins_eff = 26*7;

    TMatrixD matrix_gof_trans_eff( bins_all, bins_eff );// oldworld, newworld
    for( int ibin=1; ibin<=bins_eff;  ibin++) matrix_gof_trans_eff(ibin-1, ibin-1) = 1;

    ///////

    const int NUM_dm2 = 60;
    const int NUM_ttt = 60;

    TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -1, 2);
    TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -3, 0);

    int ittt = it14;

    double val_obj_dm2 = h1d_dm2->GetBinCenter(idm2); val_obj_dm2 = pow(10, val_obj_dm2);
    double val_obj_ttt = h1d_ttt->GetBinCenter(ittt); val_obj_ttt = pow(10, val_obj_ttt);

    ///////

    // const int NUM_dm2 = 6;
    // const int NUM_ttt = 120;

    // double array_val_dm2[NUM_dm2+1] = {0, 0.5, 2, 4, 7, 20, 40};
    // TH1D *h1d_ttt = new TH1D("h1d_tue", "h1d_tue", NUM_ttt, -3, 0);

    // int ittt = it14;

    // double val_obj_dm2 = array_val_dm2[idm2];
    // double val_obj_ttt = h1d_ttt->GetBinCenter(ittt); val_obj_ttt = pow(10, val_obj_ttt);

    // cout<<endl;
    // for(int idx=1; idx<=NUM_dm2; idx++) {
    //   cout<<TString::Format(" ---> index dm2:  %3d  %6.2f", idx, array_val_dm2[idx])<<endl;
    // }

    // cout<<endl;

    // for(int idx=1; idx<=NUM_ttt; idx++) {
    //   cout<<TString::Format(" ---> index tuu:  %3d  %6.4f", idx, pow(10, h1d_ttt->GetBinCenter(idx)))<<endl;
    // }
    // cout<<endl;

    // if( 0 ) {
    //   cout<<endl;
    //   cout<<" ---> Wilks: Asimov"<<endl;
    //   cout<<endl;


    //   /// standard order
    //   val_dm2_41         = 0;
    //   val_sin2_2theta_14 = 0.2;
    //   val_sin2_theta_24  = 0.3;
    //   osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);
    //   osc_test->Apply_oscillation();
    //   osc_test->Set_apply_POT();// meas, CV, COV: all ready

    //   //osc_test->Set_asimov2fitdata();
    //   osc_test->Set_meas2fitdata();// chi2min = 210.312885
    //                                // from /dybfs2/users/jixp/work_winxp/TOsc_nature2_after_20240522/cal_CLs_uu_BandN

    //   TH2D *h2_asimov_chi2 = new TH2D("h2_asimov_chi2", "", NUM_ttt, -3, 0, 6, 0, 6);
    //   TH2D *h2_asimov_cl = new TH2D("h2_asimov_cl", "", NUM_ttt, -3, 0, 6, 0, 6);
    //   TH2D *h2_asimov_cl2 = new TH2D("h2_asimov_cl2", "", NUM_ttt, -3, 0, 6, 0, 6);

    //   for(int jdm2=1; jdm2<=6; jdm2++) {
    // 	for(int jtuu=1; jtuu<=NUM_ttt; jtuu++) {
    // 	  double val_dm2 = array_val_dm2[jdm2];
    // 	  double val_tuu = pow(10, h1d_ttt->GetBinCenter(jtuu) );

    // 	  double pars_4v[4] = {val_dm2, val_tuu, 0, 0};
    // 	  double val_chi2 = osc_test->FCN_Pearson_FCnew( pars_4v );

    // 	  val_chi2 = val_chi2 - 210.31288;// measurement

    // 	  cout<<" ---> check "<<TString::Format("%2d %2d: %7.4f  %7.4f, chi2: %12.4f",
    // 						jdm2, jtuu, val_dm2, val_tuu,
    // 						val_chi2)<<endl;

    // 	  h2_asimov_chi2->SetBinContent(jtuu, jdm2, val_chi2);
    // 	  double val_cl = 1 - TMath::Prob(val_chi2, 2);
    // 	  h2_asimov_cl->SetBinContent(jtuu, jdm2, val_cl);
    // 	  h2_asimov_cl2->SetBinContent(jtuu, jdm2, val_cl);
    // 	}
    //   }


    //   TGraph *gh_asmiov_wilks_cl = new TGraph();
    //   gh_asmiov_wilks_cl->SetName("gh_asmiov_wilks_cl");
    //   get_CL_curve(h2_asimov_cl2, gh_asmiov_wilks_cl, 1);
    //   gh_asmiov_wilks_cl->Draw("al*");


    //   TFile *outfile = new TFile("outfile_asimov.root", "recreate");
    //   h2_asimov_chi2->Write();
    //   h2_asimov_cl->Write();
    //   gh_asmiov_wilks_cl->Write();
    //   outfile->Close();

    // }

    /////////////////////////////////////////////////////////// step_1
    /////////////////////////////////////////////////////////// step_1

    if( flag_step1_produce_toys ) {
      cout<<endl;
      cout<<" ---> Median sensitivity from freqentist CLs (step_1): pre-save toys with 3v hypothesis"<<endl;
      cout<<endl;

      /// standard order
      val_dm2_41         = 0;
      val_sin2_2theta_14 = 0.2;
      val_sin2_theta_24  = 0.3;
      osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready

      int i_toydata = 0;
      vector<double>vec_toydata_spectrum;
      TFile *outfile_toydata = new TFile("presave_3v_hypothesis_toydata_01_cv.root", "recreate");
      TTree *tree_toydata = new TTree("tree_toydata", "toydata_with_3v_hypothesis ---> i_toydata(1) for real measurement, i_toydata(2) for asimov");
      tree_toydata->Branch("i_toydata", &i_toydata, "i_toydata/I");
      tree_toydata->Branch("vec_toydata_spectrum", &vec_toydata_spectrum);

      /////// the first "toydata" is the real measurement
      osc_test->Set_meas2fitdata();
      vec_toydata_spectrum.clear();
      for(int idx=0; idx<(osc_test->matrix_tosc_fitdata_newworld.GetNcols()); idx++) {
	vec_toydata_spectrum.push_back( osc_test->matrix_tosc_fitdata_newworld(0, idx) );
      }
      i_toydata++;
      tree_toydata->Fill();

      /////// the second "toydata" is the asimov dataset
      osc_test->Set_asimov2fitdata();
      vec_toydata_spectrum.clear();
      for(int idx=0; idx<(osc_test->matrix_tosc_fitdata_newworld.GetNcols()); idx++) {
	vec_toydata_spectrum.push_back( osc_test->matrix_tosc_fitdata_newworld(0, idx) );
      }
      i_toydata++;
      tree_toydata->Fill();

      ///////

      int NUM_TOYS = 2000;// toys to be used to calculate the median result, "ccb"

      cout<<endl;
      for(int itoy=1; itoy<=NUM_TOYS; itoy++) {
	if(itoy%10==0) cout<<TString::Format(" ---> producing toy %5d / %5d", itoy, NUM_TOYS)<<endl;

	osc_test->Set_toy_variations( 1 );

	vec_toydata_spectrum.clear();
	int rows = osc_test->map_matrix_tosc_toy_pred[1].GetNcols();
	for(int idx=0; idx<rows; idx++) vec_toydata_spectrum.push_back( osc_test->map_matrix_tosc_toy_pred[1](0, idx) );
	i_toydata++;
	tree_toydata->Fill();

      }
      cout<<endl;

      ///////
      tree_toydata->Write();
      outfile_toydata->Close();

    }// step_1

    /////////////////////////////////////////////////////////// step_3'
    /////////////////////////////////////////////////////////// step_3'

    if( flag_step3_CLs ) {
      cout<<endl;
      cout<<" ---> Median sensitivity from freqentist CLs (step_3): calculate the CLs value"<<endl;
      cout<<endl;
      cout<<TString::Format(" ---> gridpoint dm2/ttt  %2d  %2d, values  %8.3f  %6.4f", idm2, ittt, val_obj_dm2, val_obj_ttt)<<endl;
      cout<<endl;

      ///
      /// CLs = CL_s+b / CL_b
      ///
      /// CL_b   = Pvalue( dchi2(3v toys) > dchi2(data) )
      /// CL_s+b = Pvalue( dchi2(4v toys) > dchi2(data) )
      ///
      /// dchi2(3v data) = chi2_4v(3v data) - chi2_3v(3v data)
      /// dchi2(4v data) = chi2_4v(4v data) - chi2_3v(4v data)
      ///
      ///

      //////////////////// read the toydata_spectrum

      const int NUM_TOYS   = 100;
      const int NUM_CYCLES = 600;

      /// total toys = num_toys * num_cycles ---> for one exclusion/sensitivity calcualtion

      map<int, TMatrixD>map_toydata_spectrum;
      map<int, double>map_prof_t24_val;
      int rows = -1;

      {
	cout<<endl;
	cout<<" ---> read the toydata_spectrum"<<endl;
	cout<<endl;

	roostr = TString::Format("input/presave_3v_hypothesis_toydata_01_cv.root");
	TFile *inputfile_toydata_cv = new TFile(roostr, "read");
	TTree *tree_toydata = (TTree*)inputfile_toydata_cv->Get("tree_toydata");

	// Declaration of leaf types
	Int_t           i_toydata;
	vector<double>  *vec_toydata_spectrum;

	// List of branches
	TBranch        *b_i_toydata;   //!
	TBranch        *b_vec_toydata_spectrum;   //!

	// Set object pointer
	vec_toydata_spectrum = 0;

	// Set branch addresses and branch pointers
	tree_toydata->SetBranchAddress("i_toydata", &i_toydata, &b_i_toydata);
	tree_toydata->SetBranchAddress("vec_toydata_spectrum", &vec_toydata_spectrum, &b_vec_toydata_spectrum);

	//
	int entries = tree_toydata->GetEntries();
	cout<<endl<<" ---> reading tree_toydata, entries "<<entries<<endl<<endl;
	tree_toydata->GetEntry(0);
	rows = vec_toydata_spectrum->size();

	for(int ientry=0; ientry<entries; ientry++) {
	  tree_toydata->GetEntry( ientry );
	  map_toydata_spectrum[i_toydata].ResizeTo(1, rows);
	  for(int idx=0; idx<rows; idx++) map_toydata_spectrum[i_toydata](0, idx) = vec_toydata_spectrum->at(idx);
	}

	delete tree_toydata;
	delete inputfile_toydata_cv;
      }

      //////////////////// for 3v: 3v_asmiov_pred, 3v_total_COV_inv, 3vToy

      map<int, TMatrixD> map_matrix_spectrum_3vToy;// index begins at 1
      TMatrixD matrix_3v_asmiov_pred(1, bins_eff);
      TMatrixD matrix_3v_total_COV_inv(bins_eff, bins_eff);

      {
	osc_test->Set_oscillation_pars(0, 0.10, 0.11, 0);
	osc_test->Apply_oscillation();
	osc_test->Set_apply_POT();// meas, CV, COV: all ready
	osc_test->Set_meas2fitdata();

	double pars_4v[4] = {0, 0.10, 0.11, 0};
	osc_test->FCN_Pearson_FCnew( pars_4v );

	{
	  cout<<endl;
	  cout<<" ---> check"<<endl;
	  cout<<" ---> rows "<<rows<<endl;
	  cout<<" ---> osc_test->matrix_tosc_chi2_COV_both_syst_stat_inv "<<osc_test->matrix_tosc_chi2_COV_both_syst_stat_inv.GetNrows()<<endl;
	  cout<<endl;
	}

	matrix_3v_asmiov_pred = osc_test->matrix_tosc_chi2_pred;
	matrix_3v_total_COV_inv = osc_test->matrix_tosc_chi2_COV_both_syst_stat_inv;

	cout<<endl<<" ---> Start producing 3v toys, total_toys: "<<NUM_TOYS*NUM_CYCLES<<endl<<endl;

	int count_3vToy = 0;
	for(int icycle=1; icycle<=NUM_CYCLES; icycle++) {
	  if(icycle%10==0) cout<<TString::Format(" ---> processing 3vToy: %4d/%4d", icycle, NUM_CYCLES)<<endl;

	  osc_test->Set_toy_variations( NUM_TOYS );
	  for(int itoy=1; itoy<=NUM_TOYS; itoy++) {
	    count_3vToy++;
	    map_matrix_spectrum_3vToy[count_3vToy].ResizeTo(1, bins_eff);
	    //map_matrix_spectrum_3vToy[count_3vToy] = osc_test->map_matrix_tosc_toy_pred[itoy];// xxx (2)
	    map_matrix_spectrum_3vToy[count_3vToy] = osc_test->map_matrix_tosc_toy_pred[itoy] * matrix_gof_trans_eff;
	  }
	}

	cout<<endl<<" ---> Finish producing 3v toys"<<endl<<endl;
      }

      //////////////////// for 4v: 4v_asmiov_pred, 4v_total_COV_inv, 4vToy

      map<int, TMatrixD> map_matrix_spectrum_4vToy;// index begins at 1
      TMatrixD matrix_4v_asmiov_pred(1, bins_eff);
      TMatrixD matrix_4v_total_COV_inv(bins_eff, bins_eff);

      {
	osc_test->Set_oscillation_pars(0, 0.10, 0.11, 0);
	osc_test->Apply_oscillation();
	osc_test->Set_apply_POT();// meas, CV, COV: all ready
	osc_test->Set_meas2fitdata();

	double pars_4v[4] = {val_obj_dm2, val_obj_ttt, 0, 0};
	osc_test->FCN_Pearson_FCnew( pars_4v );
	matrix_4v_asmiov_pred = osc_test->matrix_tosc_chi2_pred;
	matrix_4v_total_COV_inv = osc_test->matrix_tosc_chi2_COV_both_syst_stat_inv;

	cout<<endl<<" ---> Start producing 4v toys, total_toys: "<<NUM_TOYS*NUM_CYCLES<<endl<<endl;

	int count_4vToy = 0;
	for(int icycle=1; icycle<=NUM_CYCLES; icycle++) {
	  if(icycle%10==0) cout<<TString::Format(" ---> processing 4vToy: %4d/%4d", icycle, NUM_CYCLES)<<endl;

	  osc_test->Set_toy_variations( NUM_TOYS );
	  for(int itoy=1; itoy<=NUM_TOYS; itoy++) {
	    count_4vToy++;
	    map_matrix_spectrum_4vToy[count_4vToy].ResizeTo(1, bins_eff);
	    //map_matrix_spectrum_4vToy[count_4vToy] = osc_test->map_matrix_tosc_toy_pred[itoy];// xxx (3)
	    map_matrix_spectrum_4vToy[count_4vToy] = osc_test->map_matrix_tosc_toy_pred[itoy] * matrix_gof_trans_eff;
	  }
	}

	cout<<endl<<" ---> Finish producing 4v toys"<<endl<<endl;
      }

      ////////////////////
      ////////////////////
      ////////////////////

      /////// /////// /////// /////// /////// /////// ///////

      {
	int grid_idx_dm2 = idm2;
	int grid_idx_ttt = ittt;

	double data_val_dm2 = val_obj_dm2;
	double data_val_ttt = val_obj_ttt;

	vector<int>    vec_idata;
	vector<double> vec_chi2_4v_with_data;
	vector<double> vec_chi2_3v_with_data;
	vector<double> vec_dchi2_with_data;

	vector<double> vec_chi2_4v_with_4vToy;
	vector<double> vec_chi2_3v_with_4vToy;
	vector<double> vec_dchi2_with_4vToy;
	int line_dchi2_with_4vToy_GTdata = 0;

	vector<double> vec_chi2_4v_with_3vToy;
	vector<double> vec_chi2_3v_with_3vToy;
	vector<double> vec_dchi2_with_3vToy;
	int line_dchi2_with_3vToy_GTdata = 0;

	roostr = TString::Format("out_dm2_ttt_%03d_%03d.root", grid_idx_dm2, grid_idx_ttt);
	TFile *output_file = new TFile(roostr, "recreate");
	TTree *tree = new TTree("tree", "tree");

	tree->Branch("grid_idx_dm2",  &grid_idx_dm2,  "grid_idx_dm2/I" );
	tree->Branch("grid_idx_ttt",  &grid_idx_ttt,  "grid_idx_ttt/I" );

	tree->Branch("data_val_dm2",  &data_val_dm2,  "data_val_dm2/D" );
	tree->Branch("data_val_ttt",  &data_val_ttt,  "data_val_ttt/D" );

	tree->Branch("vec_idata",                &vec_idata);

	tree->Branch("vec_chi2_4v_with_data",    &vec_chi2_4v_with_data);
	tree->Branch("vec_chi2_3v_with_data",    &vec_chi2_3v_with_data);
	tree->Branch("vec_dchi2_with_data",      &vec_dchi2_with_data);

	tree->Branch("vec_chi2_4v_with_4vToy",   &vec_chi2_4v_with_4vToy);
	tree->Branch("vec_chi2_3v_with_4vToy",   &vec_chi2_3v_with_4vToy);
	tree->Branch("vec_dchi2_with_4vToy",     &vec_dchi2_with_4vToy);

	tree->Branch("vec_chi2_4v_with_3vToy",   &vec_chi2_4v_with_3vToy);
	tree->Branch("vec_chi2_3v_with_3vToy",   &vec_chi2_3v_with_3vToy);
	tree->Branch("vec_dchi2_with_3vToy",     &vec_dchi2_with_3vToy);

	// tree->Branch("line_dchi2_with_4vToy_GTdata",  &line_dchi2_with_4vToy_GTdata,  "line_dchi2_with_4vToy_GTdata/I" );
	// tree->Branch("line_dchi2_with_3vToy_GTdata",  &line_dchi2_with_3vToy_GTdata,  "line_dchi2_with_3vToy_GTdata/I" );

	///////

	cout<<endl<<" ---> Produce dchi2_with_data"<<endl<<endl;

	double user_dchi2_with_data = 0;

	for(auto it=map_toydata_spectrum.begin(); it!=map_toydata_spectrum.end(); it++) {
	  int idata = it->first;

	  TMatrixD matrix_toydata_temp(1, bins_eff);// xxx (4)
	  matrix_toydata_temp = map_toydata_spectrum[idata] * matrix_gof_trans_eff;

	  /////// calculate chi2_4v_with_data and chi2_3v_with_data
	  //TMatrixD matrix_delta_4vPred_data   = map_toydata_spectrum[idata] - matrix_4v_asmiov_pred;
	  TMatrixD matrix_delta_4vPred_data   = matrix_toydata_temp - matrix_4v_asmiov_pred;
	  TMatrixD matrix_delta_4vPred_data_T = matrix_delta_4vPred_data.T(); matrix_delta_4vPred_data.T();
	  double chi2_4v_with_data = (matrix_delta_4vPred_data * matrix_4v_total_COV_inv * matrix_delta_4vPred_data_T)(0,0);

	  //TMatrixD matrix_delta_3vPred_data   = map_toydata_spectrum[idata] - matrix_3v_asmiov_pred;
	  TMatrixD matrix_delta_3vPred_data   = matrix_toydata_temp - matrix_3v_asmiov_pred;
	  TMatrixD matrix_delta_3vPred_data_T = matrix_delta_3vPred_data.T(); matrix_delta_3vPred_data.T();
	  double chi2_3v_with_data = (matrix_delta_3vPred_data * matrix_3v_total_COV_inv * matrix_delta_3vPred_data_T)(0,0);

	  double dchi2_with_data = chi2_4v_with_data - chi2_3v_with_data;

	  if( idata==1 ) user_dchi2_with_data = dchi2_with_data;// for real measurment

	  // cout<<TString::Format(" ---> idata-%4d: chi2_4v, chi2_3v, dchi2: %10.3f  %10.3f  %10.3f",
	  // 			  idata, chi2_4v_with_data, chi2_3v_with_data, dchi2_with_data)<<endl;

	  vec_idata.push_back( idata );

	  vec_chi2_4v_with_data.push_back(chi2_4v_with_data);
	  vec_chi2_3v_with_data.push_back(chi2_3v_with_data);
	  vec_dchi2_with_data.push_back(dchi2_with_data);
	}

	///////

	cout<<endl<<" ---> Produce dchi2_with_4vToy"<<endl<<endl;

	for(auto it=map_matrix_spectrum_4vToy.begin(); it!=map_matrix_spectrum_4vToy.end(); it++) {
	  int i4vToy = it->first;

	  /////// calculate chi2_4v_with_4vToy and chi2_3v_with_4vToy
	  TMatrixD matrix_delta_4vPred_4vToy   = map_matrix_spectrum_4vToy[i4vToy] - matrix_4v_asmiov_pred;
	  TMatrixD matrix_delta_4vPred_4vToy_T = matrix_delta_4vPred_4vToy.T(); matrix_delta_4vPred_4vToy.T();
	  double chi2_4v_with_4vToy = (matrix_delta_4vPred_4vToy * matrix_4v_total_COV_inv * matrix_delta_4vPred_4vToy_T)(0,0);

	  TMatrixD matrix_delta_3vPred_4vToy   = map_matrix_spectrum_4vToy[i4vToy] - matrix_3v_asmiov_pred;
	  TMatrixD matrix_delta_3vPred_4vToy_T = matrix_delta_3vPred_4vToy.T(); matrix_delta_3vPred_4vToy.T();
	  double chi2_3v_with_4vToy = (matrix_delta_3vPred_4vToy * matrix_3v_total_COV_inv * matrix_delta_3vPred_4vToy_T)(0,0);

	  double dchi2_with_4vToy = chi2_4v_with_4vToy - chi2_3v_with_4vToy;

	  vec_chi2_4v_with_4vToy.push_back(chi2_4v_with_4vToy);
	  vec_chi2_3v_with_4vToy.push_back(chi2_3v_with_4vToy);
	  vec_dchi2_with_4vToy.push_back(dchi2_with_4vToy);

	  if( dchi2_with_4vToy >= user_dchi2_with_data ) line_dchi2_with_4vToy_GTdata++;
	}

	///////

	cout<<endl<<" ---> Produce dchi2_with_3vToy"<<endl<<endl;

	for(auto it=map_matrix_spectrum_3vToy.begin(); it!=map_matrix_spectrum_3vToy.end(); it++) {
	  int i3vToy = it->first;

	  /////// calculate chi2_4v_with_3vToy and chi2_3v_with_3vToy
	  TMatrixD matrix_delta_4vPred_3vToy   = map_matrix_spectrum_3vToy[i3vToy] - matrix_4v_asmiov_pred;
	  TMatrixD matrix_delta_4vPred_3vToy_T = matrix_delta_4vPred_3vToy.T(); matrix_delta_4vPred_3vToy.T();
	  double chi2_4v_with_3vToy = (matrix_delta_4vPred_3vToy * matrix_4v_total_COV_inv * matrix_delta_4vPred_3vToy_T)(0,0);

	  TMatrixD matrix_delta_3vPred_3vToy   = map_matrix_spectrum_3vToy[i3vToy] - matrix_3v_asmiov_pred;
	  TMatrixD matrix_delta_3vPred_3vToy_T = matrix_delta_3vPred_3vToy.T(); matrix_delta_3vPred_3vToy.T();
	  double chi2_3v_with_3vToy = (matrix_delta_3vPred_3vToy * matrix_3v_total_COV_inv * matrix_delta_3vPred_3vToy_T)(0,0);

	  double dchi2_with_3vToy = chi2_4v_with_3vToy - chi2_3v_with_3vToy;

	  vec_chi2_4v_with_3vToy.push_back(chi2_4v_with_3vToy);
	  vec_chi2_3v_with_3vToy.push_back(chi2_3v_with_3vToy);
	  vec_dchi2_with_3vToy.push_back(dchi2_with_3vToy);

	  if( dchi2_with_3vToy >= user_dchi2_with_data ) line_dchi2_with_3vToy_GTdata++;
	}

	tree->Fill();

	output_file->cd();
	tree->Write();
	output_file->Close();
      }

    }// if( flag_step3_CLs )

  }

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

////////////////////////
////////////////////////
////////////////////////

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
