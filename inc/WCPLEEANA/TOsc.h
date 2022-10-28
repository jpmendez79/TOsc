#ifndef TOsc_ana
#define TOsc_ana

#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<vector>
#include<map>
#include<set>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

/// minuit2
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

struct EventInfo {
  int           e2e_pdg;
  int           e2e_flag_FC;
  double        e2e_Etrue;
  double        e2e_Ereco;
  double        e2e_weight_xs;
  double        e2e_baseline;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////// TOsc

class TOsc {
 public:
  TOsc() {
    cout<<endl<<" ---> Hello TOsc"<<endl<<endl;

    rand = new TRandom3(0);
    
    scaleF_POT_BNB  = 1;
    scaleF_POT_NuMI = 1;
   
    ///////////////////////////
     
    dm2_41 = 0;
    sin2_2theta_14 = 0;
    sin2_theta_24  = 0;
    sin2_theta_34  = 0;

    minimization_status     = -1;
    minimization_chi2       = -1;
    minimization_dm2_41_val = -1;
    minimization_sin2_2theta_14_val = -1;
    minimization_sin2_theta_24_val  = -1;
    minimization_sin2_theta_34_val  = -1;
    minimization_dm2_41_err         = -1;
    minimization_sin2_2theta_14_err = -1;
    minimization_sin2_theta_24_err  = -1;
    minimization_sin2_theta_34_err  = -1;
    
    ///////////////////////////
    
    flag_syst_dirt   = 0;
    flag_syst_mcstat = 0;
    flag_syst_flux   = 0;
    flag_syst_geant  = 0;
    flag_syst_Xs     = 0;
    flag_syst_det    = 0;
  
    ///////////////////////////
  
    flag_NuMI_nueCC_from_intnue         = 0;
    flag_NuMI_nueCC_from_overlaynumu    = 0;
    flag_NuMI_nueCC_from_appnue         = 0;
    flag_NuMI_nueCC_from_appnumu        = 0;
    flag_NuMI_nueCC_from_overlaynueNC   = 0;
    flag_NuMI_nueCC_from_overlaynumuNC  = 0;
    
    flag_NuMI_numuCC_from_overlaynumu   = 0;
    flag_NuMI_numuCC_from_overlaynue    = 0;
    flag_NuMI_numuCC_from_appnue        = 0;
    flag_NuMI_numuCC_from_appnumu       = 0;
    flag_NuMI_numuCC_from_overlaynumuNC = 0;
    flag_NuMI_numuCC_from_overlaynueNC  = 0;
    
    flag_NuMI_CCpi0_from_overlaynumu    = 0;
    flag_NuMI_CCpi0_from_appnue         = 0;
    flag_NuMI_CCpi0_from_overlaynumuNC  = 0;
    flag_NuMI_CCpi0_from_overlaynueNC   = 0;
      
    flag_NuMI_NCpi0_from_overlaynumu    = 0;
    flag_NuMI_NCpi0_from_appnue         = 0;
    flag_NuMI_NCpi0_from_overlaynumuNC  = 0;
    flag_NuMI_NCpi0_from_overlaynueNC   = 0;
    
    ///////
  
    flag_BNB_nueCC_from_intnue       = 0;
    flag_BNB_nueCC_from_overlaynumu  = 0;
    flag_BNB_nueCC_from_appnue       = 0;
    flag_BNB_nueCC_from_appnumu      = 0;

    flag_BNB_numuCC_from_overlaynumu = 0;
    flag_BNB_numuCC_from_overlaynue  = 0;
    flag_BNB_numuCC_from_appnue      = 0;
    flag_BNB_numuCC_from_appnumu     = 0;

    flag_BNB_CCpi0_from_overlaynumu  = 0;
    flag_BNB_CCpi0_from_appnue       = 0;
  
    flag_BNB_NCpi0_from_overlaynumu  = 0;
    flag_BNB_NCpi0_from_appnue       = 0;
  
    ///////////////////////////

    NUM_TOYS = 0;
    
    ///////////////////////////
  
    default_oldworld_rows = 0;
    default_newworld_rows = 0;

  }

  ////////////////////////////////////////////////////// data members  

  TRandom3 *rand;
  
  double scaleF_POT_BNB;
  double scaleF_POT_NuMI;
  
  double dm2_41;
  double sin2_2theta_14;
  double sin2_theta_24;
  double sin2_theta_34;

  int    minimization_status;
  double minimization_chi2;
  double minimization_dm2_41_val;
  double minimization_sin2_2theta_14_val;
  double minimization_sin2_theta_24_val;
  double minimization_sin2_theta_34_val;
  double minimization_dm2_41_err;
  double minimization_sin2_2theta_14_err;
  double minimization_sin2_theta_24_err;
  double minimization_sin2_theta_34_err;

  ///////////////////////////

  bool flag_syst_dirt;
  bool flag_syst_mcstat;
  bool flag_syst_flux;
  bool flag_syst_geant;
  bool flag_syst_Xs;
  bool flag_syst_det;
  
  /////////////////////////// no specify is "CC"
  /////////////////////////// all the following true events are selected as in active volume: (cuts.h) flag_truth_inside
  /////////////////////////// flag_NuMI_nueCC_from_overlaynumu: tagged nueCC from overlay true numuCC events
  /////////////////////////// flag_NuMI_nueCC_from_overlaynueNC: tagged nueCC from overlay true nueNC events
  
  bool flag_NuMI_nueCC_from_intnue;
  bool flag_NuMI_nueCC_from_overlaynumu;
  bool flag_NuMI_nueCC_from_appnue;
  bool flag_NuMI_nueCC_from_appnumu;
  bool flag_NuMI_nueCC_from_overlaynueNC;
  bool flag_NuMI_nueCC_from_overlaynumuNC;
  
  bool flag_NuMI_numuCC_from_overlaynumu;
  bool flag_NuMI_numuCC_from_overlaynue;
  bool flag_NuMI_numuCC_from_appnue;
  bool flag_NuMI_numuCC_from_appnumu;
  bool flag_NuMI_numuCC_from_overlaynumuNC;
  bool flag_NuMI_numuCC_from_overlaynueNC;
 
  bool flag_NuMI_CCpi0_from_overlaynumu;
  bool flag_NuMI_CCpi0_from_appnue;
  bool flag_NuMI_CCpi0_from_overlaynumuNC;
  bool flag_NuMI_CCpi0_from_overlaynueNC;
  
  bool flag_NuMI_NCpi0_from_overlaynumu;
  bool flag_NuMI_NCpi0_from_appnue;
  bool flag_NuMI_NCpi0_from_overlaynumuNC;
  bool flag_NuMI_NCpi0_from_overlaynueNC;
  
  ///////
  
  bool flag_BNB_nueCC_from_intnue;
  bool flag_BNB_nueCC_from_overlaynumu;
  bool flag_BNB_nueCC_from_appnue;
  bool flag_BNB_nueCC_from_appnumu;
  bool flag_BNB_nueCC_from_overlaynueNC;
  bool flag_BNB_nueCC_from_overlaynumuNC;

  bool flag_BNB_numuCC_from_overlaynumu;
  bool flag_BNB_numuCC_from_overlaynue;
  bool flag_BNB_numuCC_from_appnue;
  bool flag_BNB_numuCC_from_appnumu;
  bool flag_BNB_numuCC_from_overlaynumuNC;
  bool flag_BNB_numuCC_from_overlaynueNC;
 
  bool flag_BNB_CCpi0_from_overlaynumu;
  bool flag_BNB_CCpi0_from_appnue;
  bool flag_BNB_CCpi0_from_overlaynumuNC;
  bool flag_BNB_CCpi0_from_overlaynueNC;
  
  bool flag_BNB_NCpi0_from_overlaynumu;
  bool flag_BNB_NCpi0_from_appnue;
  bool flag_BNB_NCpi0_from_overlaynumuNC;
  bool flag_BNB_NCpi0_from_overlaynueNC;
 
  ///////////////////////////
  
  TMatrixD matrix_transform;
  
  map<int, TH1D*>map_default_h1d_meas;
  map<int, int>map_default_h1d_meas_bins;
  map<int, double>map_default_h1d_meas_xlow;
  map<int, double>map_default_h1d_meas_xhgh;
  map<int, int>map_default_newworld_meas_bins;
  vector<double> vector_default_newworld_meas;
  TMatrixD matrix_default_newworld_meas;// assignment only once

  map<int, TH1D*>map_default_h1d_pred;
  map<int, int>map_default_h1d_pred_bins;
  map<int, double>map_default_h1d_pred_xlow;
  map<int, double>map_default_h1d_pred_xhgh;
  map<int, int>map_default_oldworld_pred_bins;
  vector<double> vector_default_oldworld_pred;
  TMatrixD matrix_default_oldworld_pred;// assignment only once
  
  TMatrixD matrix_oscillation_base_oldworld_pred;// = matrix_default_oldworld_pred - oscillation components
  TMatrixD matrix_oscillation_oldworld_pred;
  
  ///////////////////////////

  TMatrixD matrix_default_oldworld_abs_syst_addi;  // initialized, for dirt additional syst, approximation: always use the same absolute cov

  TMatrixD matrix_default_oldworld_rel_syst_flux;  // initialized
  TMatrixD matrix_default_oldworld_rel_syst_geant; // initialized
  TMatrixD matrix_default_oldworld_rel_syst_Xs;    // initialized
  TMatrixD matrix_default_oldworld_rel_syst_det;   // initialized  

  TMatrixD matrix_default_newworld_abs_syst_mcstat;// initialized, only newworld
  
  ///////////////////////////
  
  TMatrixD matrix_eff_newworld_abs_syst_addi;
  TMatrixD matrix_eff_newworld_abs_syst_mcstat;
  TMatrixD matrix_eff_newworld_abs_syst_flux;
  TMatrixD matrix_eff_newworld_abs_syst_geant;
  TMatrixD matrix_eff_newworld_abs_syst_Xs;
  TMatrixD matrix_eff_newworld_abs_syst_det;
  TMatrixD matrix_eff_newworld_abs_syst_total;
  
  TMatrixD matrix_eff_newworld_meas;
  TMatrixD matrix_eff_newworld_pred;
  TMatrixD matrix_eff_newworld_noosc;
    
  TMatrixD matrix_fitdata_newworld;

  ///////////////////////////

  int NUM_TOYS;
  map<int, TMatrixD>map_matrix_toy_pred;
  
  ///////////////////////////

  vector<double>vector_NuMI_nueCC_from_intnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo;        vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo;
  vector<double>vector_NuMI_nueCC_from_overlaynumu_scaleFPOT;   vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo;   vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo;
  vector<double>vector_NuMI_nueCC_from_appnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo;        vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo;
  vector<double>vector_NuMI_nueCC_from_appnumu_scaleFPOT;       vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_appnumu_FC_eventinfo;       vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_appnumu_PC_eventinfo;
  vector<double>vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT;  vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo;  vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo;
  vector<double>vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo; vector< vector<EventInfo> >vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo;
  
  vector<double>vector_NuMI_numuCC_from_overlaynumu_scaleFPOT;  vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo;  vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo;
  vector<double>vector_NuMI_numuCC_from_overlaynue_scaleFPOT;   vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynue_FC_eventinfo;   vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynue_PC_eventinfo;
  vector<double>vector_NuMI_numuCC_from_appnue_scaleFPOT;       vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo;       vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo;
  vector<double>vector_NuMI_numuCC_from_appnumu_scaleFPOT;      vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_appnumu_FC_eventinfo;      vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_appnumu_PC_eventinfo;
  vector<double>vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT;vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo;vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo;
  vector<double>vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo; vector< vector<EventInfo> >vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo;
  
  vector<double>vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT;   vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo;   vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo;
  vector<double>vector_NuMI_CCpi0_from_appnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo;        vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo;
  vector<double>vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo; vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo;
  vector<double>vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT;  vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo;  vector< vector<EventInfo> >vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo;
  
  vector<double>vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT;   vector< vector<EventInfo> >vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo; 
  vector<double>vector_NuMI_NCpi0_from_appnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_NuMI_NCpi0_from_appnue_eventinfo;
  vector<double>vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo;
  vector<double>vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT;  vector< vector<EventInfo> >vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo;
  
  ///////

  vector<double>vector_BNB_nueCC_from_intnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_intnue_FC_eventinfo;        vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_intnue_PC_eventinfo;
  vector<double>vector_BNB_nueCC_from_overlaynumu_scaleFPOT;   vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo;   vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo;
  vector<double>vector_BNB_nueCC_from_appnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_appnue_FC_eventinfo;        vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_appnue_PC_eventinfo;
  vector<double>vector_BNB_nueCC_from_appnumu_scaleFPOT;       vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_appnumu_FC_eventinfo;       vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_appnumu_PC_eventinfo;
  vector<double>vector_BNB_nueCC_from_overlaynueNC_scaleFPOT;  vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_overlaynueNC_FC_eventinfo;  vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_overlaynueNC_PC_eventinfo;
  vector<double>vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_overlaynumuNC_FC_eventinfo; vector< vector<EventInfo> >vector_vector_BNB_nueCC_from_overlaynumuNC_PC_eventinfo;
  
  vector<double>vector_BNB_numuCC_from_overlaynumu_scaleFPOT;  vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo;  vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo;
  vector<double>vector_BNB_numuCC_from_overlaynue_scaleFPOT;   vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynue_FC_eventinfo;   vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynue_PC_eventinfo;
  vector<double>vector_BNB_numuCC_from_appnue_scaleFPOT;       vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_appnue_FC_eventinfo;       vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_appnue_PC_eventinfo;
  vector<double>vector_BNB_numuCC_from_appnumu_scaleFPOT;      vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_appnumu_FC_eventinfo;      vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_appnumu_PC_eventinfo;
  vector<double>vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT;vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynumuNC_FC_eventinfo;vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynumuNC_PC_eventinfo;
  vector<double>vector_BNB_numuCC_from_overlaynueNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynueNC_FC_eventinfo; vector< vector<EventInfo> >vector_vector_BNB_numuCC_from_overlaynueNC_PC_eventinfo;
  
  vector<double>vector_BNB_CCpi0_from_overlaynumu_scaleFPOT;   vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo;   vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo;
  vector<double>vector_BNB_CCpi0_from_appnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_appnue_FC_eventinfo;        vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_appnue_PC_eventinfo;
  vector<double>vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_overlaynumuNC_FC_eventinfo; vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_overlaynumuNC_PC_eventinfo;
  vector<double>vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT;  vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_overlaynueNC_FC_eventinfo;  vector< vector<EventInfo> >vector_vector_BNB_CCpi0_from_overlaynueNC_PC_eventinfo;
  
  vector<double>vector_BNB_NCpi0_from_overlaynumu_scaleFPOT;   vector< vector<EventInfo> >vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo; 
  vector<double>vector_BNB_NCpi0_from_appnue_scaleFPOT;        vector< vector<EventInfo> >vector_vector_BNB_NCpi0_from_appnue_eventinfo;
  vector<double>vector_BNB_NCpi0_from_overlaynumuNC_scaleFPOT; vector< vector<EventInfo> >vector_vector_BNB_NCpi0_from_overlaynumuNC_eventinfo;
  vector<double>vector_BNB_NCpi0_from_overlaynueNC_scaleFPOT;  vector< vector<EventInfo> >vector_vector_BNB_NCpi0_from_overlaynueNC_eventinfo;
  
  ///////////////////////////

  ///////////////////////////

  ///////////////////////////

  ////////////////////////////////////////////////////// member functions
  
  void Set_default_cv_cov(TString default_cv_file, TString default_dirtadd_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir);

  void Set_oscillation_base(TString default_eventlist_dir);
  void Set_oscillation_base_subfunc(TString strfile_mcPOT, TString strfile_dataPOT, vector<double> *vec_ratioPOT, TString strfile_mc_e2e, TString str_treename, vector< vector<EventInfo> > *vec_vec_eventinfo);  
  /// matrix oscillation effect result oldworld = matrix oscillation base + matrix oscillation_ ffect;

  void Set_oscillation_base_minus(vector<double> *vec_ratioPOT, vector< vector<EventInfo> > *vec_vec_eventinfo, int pred_channel_index, TString str_osc_mode);
  void Set_oscillation_base_added(vector<double> *vec_ratioPOT, vector< vector<EventInfo> > *vec_vec_eventinfo, int pred_channel_index, TString str_osc_mode);

  void Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34) {
    dm2_41 = val_dm2_41;
    sin2_2theta_14 = val_sin2_2theta_14;
    sin2_theta_24  = val_sin2_theta_24;
    sin2_theta_34  = val_sin2_theta_34;
  }
  
  void Apply_oscillation();
  double Prob_oscillaion(double Etrue, double baseline, int strflag_osc);

  void Set_apply_POT();

  void Set_toy_variations(int num_toys);
  void Set_toy2fitdata(int itoy) {
    if( itoy > NUM_TOYS ) { cerr<<" ERROR: itoy > NUM_TOYS"<<endl; exit(1); }
    matrix_fitdata_newworld = map_matrix_toy_pred[itoy];
  }
  
  void Set_meas2fitdata()   { matrix_fitdata_newworld   = matrix_eff_newworld_meas; }
  void Set_asimov2fitdata() { matrix_fitdata_newworld   = matrix_eff_newworld_pred; }
  void Set_noosc2fitdata()  { matrix_fitdata_newworld   = matrix_eff_newworld_noosc;}
  void Set_asimov2noosc()   { matrix_eff_newworld_noosc = matrix_eff_newworld_pred; }

  void Plot_user();
  
  void Minimization_OscPars_FullCov(double init_dm2_41, double init_sin2_2theta_14, double init_sin2_theta_24, double init_sin2_theta_34, TString roostr_flag_fixpar);
  double FCN(const double *par);

  double func_CLs(double eff_d4v, double eff_d3v, double eff_dd) {
    double result = ( 1+TMath::Erf( (eff_d4v-eff_dd)/sqrt(8*fabs(eff_d4v)) ) )
      / ( 1+TMath::Erf( (eff_d3v-eff_dd)/sqrt(8*fabs(eff_d3v)) ) );
    return result;
  }
  
  ////////////////////////////////////////////////////// data members

 private:
  int default_oldworld_rows;
  int default_newworld_rows;
  
};

#endif
