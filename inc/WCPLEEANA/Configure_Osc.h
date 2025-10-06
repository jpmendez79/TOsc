namespace Configure_Osc
{
  /////////////////////////// default files for spectra and covariance matrixes
  
  TString default_cv_file      = "/dybfs2/users/jixp/work_winxp/TOsc_nature2_after_20240522/version_to_Jesse_2025/dataset_input/merge.root";
  TString default_dirtadd_file = "/dybfs2/users/jixp/work_winxp/TOsc_nature2_after_20240522/version_to_Jesse_2025/dataset_input/merge.root";
  TString default_mcstat_file  = "/dybfs2/users/jixp/work_winxp/TOsc_nature2_after_20240522/version_to_Jesse_2025/dataset_input/0.log";// mc_stat from no oscillation
  TString default_fluxXs_dir   = "/dybfs2/users/jixp/work_winxp/TOsc_nature2_after_20240522/version_to_Jesse_2025/dataset_input/hist_rootfiles/XsFlux_edit/";// hack flux for NuMI
  TString default_detector_dir = "/dybfs2/users/jixp/work_winxp/TOsc_nature2_after_20240522/version_to_Jesse_2025/dataset_input/hist_rootfiles/DetVar_edit/";// hack oscillation: use the intrinsic
  TString default_eventlist_dir= "/dybfs2/users/jixp/work_winxp/TOsc_nature2_after_20240522/version_to_Jesse_2025/dataset_input/";
  
  /////////////////////////// By default 1. Don't modify the two values below. (X. JI)

  bool flag_apply_oscillation_BNB  = 1;
  bool flag_apply_oscillation_NuMI = 1;

  bool flag_goodness_of_fit_CNP    = 0;

  /////////////////////////// Systematics

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
