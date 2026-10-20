#include "stdlib.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
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


  // int seq = std::atoi(argv[1]);

  // Define grid used for coarse grid scan. This should remain the same as the
  // CLs calculation testing
  // Sequence Number int seq = i * (Nj * Nk) + j * Nk + k;
  // int k = seq % Nk;
  // int j = (seq / Nk) % Nj;
  // int i = seq / (Nj * Nk);
  // i = ittt
  // j = idm2
  // k = ig2
  const int NUM_dm2 = 60;
  const int NUM_ttt = 60;
  const int NUM_g2 = 60;
  const double DM2_LO = -1, DM2_HI = 2;
  const double TTT_LO = -3, TTT_HI = 0;
  const double G2_LO = 0, G2_HI = 4 * M_PI;
  int spectrum_index = -1;
  TString specpath;
  /* setup: load files, initialize, etc. (once) */
  if (argc < 3 || strcmp(argv[1], "-f") != 0) {
	cerr << "Usage: " << argv[0] << " -f <spectrum.root> ..." << endl;
	exit(1);
}
specpath = argv[2];

  // if (strcmp(argv[1], "-f") == 0) {
  //	  stringstream convert(argv[2]);
  //	  if (!(convert >> specpath)) {
  //		cerr << " Input spectrum file" << endl;
  //		exit(1);
  //		  }
  // }
	if (strcmp(argv[3], "-imeas") == 0) {
	  stringstream convert(argv[4]);
	  if (!(convert >> spectrum_index)) {
		cerr << " Input index for spectrum" << endl;
		exit(1);
		  }
		}

  // if (!inputfile_toydata_cv || inputfile_toydata_cv->iszombie()) {
  //	std::cerr << "error: cannot open input file\n";
  //	return 1;
  // }


  TFile *inputfile_toydata_cv = TFile::Open(specpath, "READ");
  TTree *tree_toydata = (TTree *)inputfile_toydata_cv->Get("tree_toydata");
  if (!tree_toydata) {
	std::cerr << "Error: tree_toydata not found\n";
	return 1;
  }

  if (spectrum_index < 0 || spectrum_index >= tree_toydata->GetEntries()) {
	cerr << "Error: -imeas " << spectrum_index << " out of range [0, "
	  << tree_toydata->GetEntries() << ")\n";
	return 1;
  }

  vector<double> spectrum;
{
	TTreeReader reader(tree_toydata);
	TTreeReaderValue<vector<double>> vec_toydata_spectrum(reader, "vec_toydata_spectrum");
	if (reader.SetEntry(spectrum_index) != TTreeReader::kEntryValid) {
		cerr << "Error: cannot read entry " << spectrum_index << "\n";
		return 1;
	}
	spectrum = *vec_toydata_spectrum;
}   // reader is destroyed here, before the file closes
inputfile_toydata_cv->Close();

  // TTreeReader reader(tree_toydata);
  // TTreeReaderValue<vector<double>> vec_toydata_spectrum(reader,
  //														"vec_toydata_spectrum");

  // vector<vector<double>> spectrum_cache;
  // while (reader.Next()) {
  //	spectrum_cache.push_back(*vec_toydata_spectrum);
  // }
  // inputfile_toydata_cv->Close();

  // cout << "Cached entries: " << spectrum_cache.size() << endl;
  // const int N = spectrum_cache.size();

for (int i = 5; i < argc; i++) {
	int seq = atoi(argv[i]);
	if (seq == -1) break;
	/* existing calculation + print for this seq */
  // Convert sequence number into bin indexes
  int ig2 = seq % NUM_g2;
  int idm2 = (seq / NUM_g2) % NUM_dm2;
  int ittt = seq / (NUM_dm2 * NUM_g2);

  auto center = [](int idx, int n, double lo, double hi) {
	return lo + (idx + 0.5) * (hi - lo) / n;
  };

  double test_dm2 = pow(10, center(idm2, NUM_dm2, DM2_LO, DM2_HI));
  double test_ttt = pow(10, center(ittt, NUM_ttt, TTT_LO, TTT_HI));
  double test_g2  =         center(ig2,  NUM_g2,  G2_LO,  G2_HI);

  // // Setup the grid
  // TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, DM2_LO, DM2_HI);
  // TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, TTT_LO, TTT_HI);
  // TH1D *h1d_g2 = new TH1D("h1d_g2", "h1d_g2", NUM_g2, G2_LO, G2_HI);
  //
  // // Convert bin indexes to actual parameters
  // double test_dm2 = h1d_dm2->GetBinCenter(idm2);
  // test_dm2 = pow(10, test_dm2);
  // double test_ttt = h1d_ttt->GetBinCenter(ittt);
  // test_ttt = pow(10, test_ttt);
  // double test_g2 = h1d_g2->GetBinCenter(ig2);


  // --------------------------------------------------
  // TOsc setup (identical in both original files) -- one-time per process,
  // shared by stage 1 and stage 2 so they operate on exactly the same state.
  // --------------------------------------------------




  // int bins_eff = 181;
  // // setting Up the bins and spectrum mask
  int bins_all = 26 * 14;
  int bins_eff = 26 * 7;

  TMatrixD matrix_gof_trans_eff( bins_all, bins_eff );// oldworld, newworld
  for( int ibin=1; ibin<=bins_eff;  ibin++) matrix_gof_trans_eff(ibin-1, ibin - 1) = 1;

  // For testing Purposes Create an measurement to be used
  // double seed_val_dm2 = 7.3;
  // double seed_val_ttt = 0.23;
  // double seed_val_g2 = 2.5 * M_PI;
  // double pars_test_meas[3] = {seed_val_dm2, seed_val_ttt, seed_val_g2};
  //
  // // Generate the test measurement Asimov
  // osc_test->Set_oscillation_pars(pars_test_meas[0], pars_test_meas[1], 0, 0, pars_test_meas[2]);
  // osc_test->Apply_oscillation();
  // osc_test->Set_apply_POT(); // meas, CV, COV: all ready
  // osc_test->Set_asimov2fitdata(); // Set the asimov generated prediction as the test prediction
  //
  // // Pull out the measurement spectrum
  // TMatrixD matrix_meas(1, bins_eff);
  // matrix_meas = osc_test->matrix_tosc_fitdata_newworld * matrix_gof_trans_eff;

  // For Testing purposes load pseudo-experiment from spectrum
  // TString specpath =
  //     "jpmendez_presave_spectrums_dm2-6.1_ttt-0.33_g2-2.5pi.root";
  // TString specpath = "jpmendez_presave_spectrums_dm2-0_ttt-0_g2-0pi.root";
  // TString specpath = "jpmendez_presave_spectrums_dm2-50_ttt-0.33_g2-2.5pi.root";

  TMatrixD matrix_meas(1, bins_eff);


  // auto &spectrum = spectrum_cache[spectrum_index];
   for (int j = 0; j < bins_eff; j++) {
	  matrix_meas(0, j) = spectrum[j];
	}

  // matrix_meas = matrix_meas * matrix_gof_trans_eff;
  // Generate the Prediction
  // Holding variables for prediction
  TMatrixD matrix_pred(1, bins_eff);
  TMatrixD matrix_total_COV_inv(bins_eff, bins_eff);
  double test_par[4] = {test_dm2, test_ttt, 0, test_g2};
  //
  osc_test->Set_oscillation_pars(test_par[0], test_par[1], test_par[2], 0, test_par[3]);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT(); // meas, CV, COV: all ready
  osc_test->FCN_Pearson_FCnew(test_par);

  // Save prediction and inverted chi2 out of tosc internal to local variables
  matrix_pred = osc_test->matrix_tosc_chi2_pred;
  matrix_total_COV_inv = osc_test->matrix_tosc_chi2_COV_both_syst_stat_inv;

  // Calculate chisquare
  TMatrixD matrix_delta = matrix_meas - matrix_pred;
  TMatrixD matrix_delta_T = matrix_delta.T();
  matrix_delta.T();
  double test_chi2 = (matrix_delta * matrix_total_COV_inv * matrix_delta_T)(0,0);



  std::ostringstream line;
  line << std::setprecision(std::numeric_limits<double>::max_digits10);

	line << spectrum_index << " " << seq << " " << test_chi2 << " " << ittt << " " << test_ttt << " "
  << idm2 << " " << test_dm2 << " " << ig2 << " " << test_g2 << "\n";

  // Save Result to file
  // std::ofstream out("output/ipv_" + std::to_string(seq) + ".txt");
  // out << line.str() << std::flush;
  // Save Result to stdout to be captured by sed in job log
  std::cout << "RESULT " << line.str() << std::flush;

}





  if (display) {
	cout << " Enter Ctrl+c to end the program" << endl;
	cout << " Enter Ctrl+c to end the program" << endl;
	cout << endl;
	theApp.Run();
  }

  return 0;
}
