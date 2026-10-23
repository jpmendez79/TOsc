#include "stdlib.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>

#include <vector>

#include "WCPLEEANA/TOsc.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"
using namespace std;
using namespace std::chrono;


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  TString roostr = "";

  cout << endl << " ---> A Hello story ..." << endl << endl;

  int display = 0;
  if (!display) {
	gROOT->SetBatch(1);
  }
  int dummy_argc = 1;
  TApplication theApp("theApp", &dummy_argc, argv);


  double scaleF_POT_BNB = 1;
  double scaleF_POT_NuMI = 1;

  TOsc *osc_test = new TOsc();

  osc_test->tosc_scaleF_POT_BNB = scaleF_POT_BNB;
  osc_test->tosc_scaleF_POT_NuMI = scaleF_POT_NuMI;

  osc_test->flag_apply_oscillation_BNB =
	  Configure_Osc::flag_apply_oscillation_BNB;
  osc_test->flag_apply_oscillation_NuMI =
	  Configure_Osc::flag_apply_oscillation_NuMI;

  osc_test->flag_goodness_of_fit_CNP = Configure_Osc::flag_goodness_of_fit_CNP;

  osc_test->flag_syst_dirt = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det = Configure_Osc::flag_syst_det;

  osc_test->flag_NuMI_nueCC_from_intnue =
	  Configure_Osc::flag_NuMI_nueCC_from_intnue;
  osc_test->flag_NuMI_nueCC_from_overlaynumu =
	  Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
  osc_test->flag_NuMI_nueCC_from_appnue =
	  Configure_Osc::flag_NuMI_nueCC_from_appnue;
  osc_test->flag_NuMI_nueCC_from_appnumu =
	  Configure_Osc::flag_NuMI_nueCC_from_appnumu;
  osc_test->flag_NuMI_nueCC_from_overlaynueNC =
	  Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
  osc_test->flag_NuMI_nueCC_from_overlaynumuNC =
	  Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;

  osc_test->flag_NuMI_numuCC_from_overlaynumu =
	  Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue =
	  Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue =
	  Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu =
	  Configure_Osc::flag_NuMI_numuCC_from_appnumu;
  osc_test->flag_NuMI_numuCC_from_overlaynumuNC =
	  Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
  osc_test->flag_NuMI_numuCC_from_overlaynueNC =
	  Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;

  osc_test->flag_NuMI_CCpi0_from_overlaynumu =
	  Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue =
	  Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  osc_test->flag_NuMI_CCpi0_from_overlaynumuNC =
	  Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_CCpi0_from_overlaynueNC =
	  Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;

  osc_test->flag_NuMI_NCpi0_from_overlaynumu =
	  Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue =
	  Configure_Osc::flag_NuMI_NCpi0_from_appnue;
  osc_test->flag_NuMI_NCpi0_from_overlaynumuNC =
	  Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_NCpi0_from_overlaynueNC =
	  Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;

  osc_test->flag_BNB_nueCC_from_intnue =
	  Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu =
	  Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue =
	  Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu =
	  Configure_Osc::flag_BNB_nueCC_from_appnumu;
  osc_test->flag_BNB_nueCC_from_overlaynueNC =
	  Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
  osc_test->flag_BNB_nueCC_from_overlaynumuNC =
	  Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;

  osc_test->flag_BNB_numuCC_from_overlaynumu =
	  Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue =
	  Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue =
	  Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu =
	  Configure_Osc::flag_BNB_numuCC_from_appnumu;
  osc_test->flag_BNB_numuCC_from_overlaynumuNC =
	  Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
  osc_test->flag_BNB_numuCC_from_overlaynueNC =
	  Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;

  osc_test->flag_BNB_CCpi0_from_overlaynumu =
	  Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue =
	  Configure_Osc::flag_BNB_CCpi0_from_appnue;
  osc_test->flag_BNB_CCpi0_from_overlaynumuNC =
	  Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_CCpi0_from_overlaynueNC =
	  Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;

  osc_test->flag_BNB_NCpi0_from_overlaynumu =
	  Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue =
	  Configure_Osc::flag_BNB_NCpi0_from_appnue;
  osc_test->flag_BNB_NCpi0_from_overlaynumuNC =
	  Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_NCpi0_from_overlaynueNC =
	  Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;

  /////// set only one time
  /////// set only one time

  osc_test->Set_default_cv_cov(
	  Configure_Osc::default_cv_file, Configure_Osc::default_dirtadd_file,
	  Configure_Osc::default_mcstat_file, Configure_Osc::default_fluxXs_dir,
	  Configure_Osc::default_detector_dir);
  osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);

  if (argc < 5) {
	std::cerr << "Usage: " << argv[0] << " spectrum_file index init_ttt init_dm2 init_g2\n";
	return 1;
  }
  TString specpath;
  specpath = argv[1];
  int index = atoi(argv[2]);
  double init_ttt = atof(argv[3]);
  double arg_dm2 = atof(argv[4]);
  double arg_g2 = atof(argv[5]);


  //	// For testing Purposes Create an measurement to be used
  // double seed_val_dm2 = 7.3;
  // double seed_val_ttt = 0.23;
  // double seed_val_g2 = 2.5 * m_pi;
  // double pars_test_meas[3] = {seed_val_dm2, seed_val_ttt, seed_val_g2};

  // //	// generate the test measurement asimov
  // osc_test->set_oscillation_pars(pars_test_meas[0], pars_test_meas[1], 0, 0, pars_test_meas[2]);
  // osc_test->apply_oscillation();
  // osc_test->set_apply_pot(); // meas, cv, cov: all ready
  // osc_test->set_toy_variations(1);
  // osc_test->set_toy2fitdata(1);
  // osc_test->Set_asimov2fitdata(); // Set the asimov generated prediction as the test prediction


  // Load the psuedo from file
  // TString specpath =
  // "jpmendez_presave_spectrums_dm2-6.1_ttt-0.33_g2-2.5pi.root";
  // TString specpath = "jpmendez_presave.root";
  TFile *inputfile_toydata_cv = TFile::Open(specpath, "READ");

  if (!inputfile_toydata_cv || inputfile_toydata_cv->IsZombie()) {
	std::cerr << "Error: cannot open input file\n";
	return 1;
  }

  TTree *tree_toydata = (TTree *)inputfile_toydata_cv->Get("tree_toydata");

  if (!tree_toydata) {
	std::cerr << "Error: tree_toydata not found\n";
	return 1;
  }

  TTreeReader reader(tree_toydata);
  TTreeReaderValue<vector<double>> vec_toydata_spectrum(reader,
														"vec_toydata_spectrum");

  vector<vector<double>> spectrum_cache;
  while (reader.Next()) {
	spectrum_cache.push_back(*vec_toydata_spectrum);
  }
  inputfile_toydata_cv->Close();

  cout << "Cached entries: " << spectrum_cache.size() << endl;
  const int N = spectrum_cache.size();

  if (index < 0 || index >= N) {
	std::cerr << "Error: index " << index << " out of range [0, " << N - 1 << "]\n";
	return 1;
  }

  // Set psuedo to fit spectrum (FCN only uses the BNB bins 0..26*7-1)
  int bins_eff = 26 * 7;
  auto &spectrum = spectrum_cache[index];
  for (int j = 0; j < bins_eff; j++) {
	osc_test->matrix_tosc_fitdata_newworld(0, j) = spectrum[j];
  }


  // cout << meas_index <<  " " << init_ttt << " " << init_dm2_41 << " " <<
  // init_g2 << endl;

  // Diagnostics: chi2 at the generation truth, the null, and the fit start
  double pars_truth[4] = {6.1, 0.33, 0, 7.85};
  double pars_null[4] = {6.1, 0, 0, 7.85};
  double pars_start[4] = {arg_dm2, init_ttt, 0, arg_g2};
  double chi2_truth = osc_test->FCN(pars_truth);
  double chi2_null = osc_test->FCN(pars_null);
  double chi2_start = osc_test->FCN(pars_start);

  osc_test->Minimization_OscPars_FullCov(arg_dm2, init_ttt, 0, 0, arg_g2, "t24_t34", index);

  double pars_fit[4] = {osc_test->minimization_dm2_41_val,
						osc_test->minimization_sin2_2theta_14_val,
						osc_test->minimization_sin2_theta_24_val,
						osc_test->minimization_g2_val};
  double chi2_fit = osc_test->FCN(pars_fit);
  // DIAG N chi2_truth chi2_null chi2_start chi2_fit
  cout << TString::Format("DIAG %2d %8.4f %8.4f %8.4f %8.4f", index, chi2_truth,
						  chi2_null, chi2_start, chi2_fit)
	   << endl;

  if (display) {
		cout << " enter ctrl+c to end the program" << endl;
		cout << " enter ctrl+c to end the program" << endl;
		cout << endl;
		theApp.Run();
  }

  return 0;
}
