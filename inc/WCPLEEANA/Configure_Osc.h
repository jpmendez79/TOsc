namespace Configure_Osc
{
  /////////////////////////// default files for spectra and covariance matrixes

  TString default_cv_file       = "./data_inputs/yyyd_BNBplusNuMI_OSC/merge.root";
  TString default_dirtadd_file  = "./data_inputs/yyyd_BNBplusNuMI_OSC/merge.root";
  TString default_mcstat_file   = "./data_inputs/yyyd_BNBplusNuMI_OSC/mc_stat/0.log";// mc_stat from no-oscillation    
  TString default_fluxXs_dir    = "./data_inputs/yyyd_BNBplusNuMI_OSC/XsFlux_edit/";// hack flux for NuMI: only one knob represents all NuMI knobs
  TString default_detector_dir  = "./data_inputs/yyyd_BNBplusNuMI_OSC/DetVar_edit/";// hack oscillation part to use the intrinsic nue
  TString default_eventlist_dir = "./data_inputs/yyyd_BNBplusNuMI_OSC/";
  
  ///////////////////////////

  bool flag_syst_dirt   = 1;
  bool flag_syst_mcstat = 1;
  bool flag_syst_flux   = 1;
  bool flag_syst_geant  = 1;
  bool flag_syst_Xs     = 1;
  bool flag_syst_det    = 1;
  
  /////////////////////////// no specify is "CC"
  /////////////////////////// all the following true events are selected as in active volume: (cuts.h) flag_truth_inside
  
  bool flag_NuMI_nueCC_from_intnue        = 1;// ####### work
  bool flag_NuMI_nueCC_from_overlaynumu   = 1;// ####### work
  bool flag_NuMI_nueCC_from_appnue        = 1;// ####### work
  bool flag_NuMI_nueCC_from_appnumu       = 0;// N/A
  bool flag_NuMI_nueCC_from_dirtnue       = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 1/557
  bool flag_NuMI_nueCC_from_dirtnumu      = 0;// approximation: ignore osc-effect. 
  bool flag_NuMI_nueCC_from_overlaynueNC  = 1;// ####### work
  bool flag_NuMI_nueCC_from_overlaynumuNC = 1;// ####### work
  
  bool flag_NuMI_numuCC_from_overlaynumu  = 1;// ####### work
  bool flag_NuMI_numuCC_from_overlaynue   = 0;// N/A
  bool flag_NuMI_numuCC_from_appnue       = 1;// ####### work
  bool flag_NuMI_numuCC_from_appnumu      = 0;// N/A
  bool flag_NuMI_numuCC_from_dirtnue      = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 693.4/129661
  bool flag_NuMI_numuCC_from_dirtnumu     = 0;// approximation: ignore osc-effect.
  bool flag_NuMI_numuCC_from_overlaynumuNC= 1;// ####### work
  bool flag_NuMI_numuCC_from_overlaynueNC = 1;// ####### work
  
  bool flag_NuMI_CCpi0_from_overlaynumu   = 1;// ####### work
  bool flag_NuMI_CCpi0_from_overlaynue    = 0;// approximation: ignore osc-effect. DocDB-36268 (NuMI, pi0-KE): nueCC/data = 4.5/255
  bool flag_NuMI_CCpi0_from_appnue        = 1;// ####### work
  bool flag_NuMI_CCpi0_from_appnumu       = 0;// approximation: ignore osc-effect. See flag_NuMI_CCpi0_from_overlaynue
  bool flag_NuMI_CCpi0_from_dirtnue       = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 10.4/7953
  bool flag_NuMI_CCpi0_from_dirtnumu      = 0;// approximation: ignore osc-effect.
  bool flag_NuMI_CCpi0_from_overlaynumuNC = 1;// ####### work
  bool flag_NuMI_CCpi0_from_overlaynueNC  = 1;// ####### work
  
  bool flag_NuMI_NCpi0_from_overlaynumu   = 1;// ####### work
  bool flag_NuMI_NCpi0_from_overlaynue    = 0;// approximation: ignore osc-effect. DocDB-36268 (NuMI, pi0-KE): nueCC/data = 34.8/874
  bool flag_NuMI_NCpi0_from_appnue        = 1;// ####### work
  bool flag_NuMI_NCpi0_from_appnumu       = 0;// approximation: ignore osc-effect. See flag_NuMI_NCpi0_from_overlaynue
  bool flag_NuMI_NCpi0_from_dirtnue       = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 188.3/5936.0
  bool flag_NuMI_NCpi0_from_dirtnumu      = 0;// approximation: ignore osc-effect. 
  bool flag_NuMI_NCpi0_from_overlaynumuNC = 1;// ####### work
  bool flag_NuMI_NCpi0_from_overlaynueNC  = 1;// ####### work
  
  ///////
  
  bool flag_BNB_nueCC_from_intnue        = 1;// ####### work
  bool flag_BNB_nueCC_from_overlaynumu   = 1;// ####### work
  bool flag_BNB_nueCC_from_appnue        = 1;// ####### work
  bool flag_BNB_nueCC_from_appnumu       = 0;
  bool flag_BNB_nueCC_from_dirtnue       = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 1/557
  bool flag_BNB_nueCC_from_dirtnumu      = 0;// approximation: ignore osc-effect.
  bool flag_BNB_nueCC_from_overlaynueNC  = 1;// ####### work
  bool flag_BNB_nueCC_from_overlaynumuNC = 1;// ####### work

  bool flag_BNB_numuCC_from_overlaynumu  = 1;// ####### work
  bool flag_BNB_numuCC_from_overlaynue   = 0; 
  bool flag_BNB_numuCC_from_appnue       = 1;// ####### work
  bool flag_BNB_numuCC_from_appnumu      = 0;
  bool flag_BNB_numuCC_from_dirtnue      = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 693.4/129661
  bool flag_BNB_numuCC_from_dirtnumu     = 0;// approximation: ignore osc-effect.
  bool flag_BNB_numuCC_from_overlaynumuNC= 1;// ####### work
  bool flag_BNB_numuCC_from_overlaynueNC = 1;// ####### work
 
  bool flag_BNB_CCpi0_from_overlaynumu   = 1;// ####### work
  bool flag_BNB_CCpi0_from_overlaynue    = 0;// approximation: ignore osc-effect. LEE PRD paper(pi0-KE): nueCC/data =  18.0/7953
  bool flag_BNB_CCpi0_from_appnue        = 1;// ####### work
  bool flag_BNB_CCpi0_from_appnumu       = 0;// approximation: ignore osc-effect. flag_BNB_CCpi0_from_overlaynue
  bool flag_BNB_CCpi0_from_dirtnue       = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 10.4/7953
  bool flag_BNB_CCpi0_from_dirtnumu      = 0;// approximation: ignore osc-effect.
  bool flag_BNB_CCpi0_from_overlaynumuNC = 1;// ####### work
  bool flag_BNB_CCpi0_from_overlaynueNC  = 1;// ####### work
  
  bool flag_BNB_NCpi0_from_overlaynumu   = 1;// ####### work
  bool flag_BNB_NCpi0_from_overlaynue    = 0;// approximation: ignore osc-effect. LEE PRD paper(pi0-KE): nueCC/data = 42.2/5936
  bool flag_BNB_NCpi0_from_appnue        = 1;// ####### work
  bool flag_BNB_NCpi0_from_appnumu       = 0;// approximation: ignore osc-effect. flag_BNB_NCpi0_from_overlaynue
  bool flag_BNB_NCpi0_from_dirtnue       = 0;// approximation: ignore osc-effect. LEE PRD paper(BNB case): dirt/data = 188.3/5936.0
  bool flag_BNB_NCpi0_from_dirtnumu      = 0;// approximation: ignore osc-effect.
  bool flag_BNB_NCpi0_from_overlaynumuNC = 1;// ####### work
  bool flag_BNB_NCpi0_from_overlaynueNC  = 1;// ####### work
 
  ///////////////////////////
  
}
