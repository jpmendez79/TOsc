#include "stdlib.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>
using namespace std;

#include <vector>

#include "WCPLEEANA/TOsc.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Combines the two-stage create_dists.cxx -> writechi2obs.cxx pipeline into a
// single per-grid-point program. Stage 1 (toy Delta-chi2 generation) and stage
// 2 (observed CLs/confidence calculation) now share one process and one TOsc
// instance, so the intermediate ROOT file that used to connect them is no
// longer written at all.
//
int main(int argc, char **argv) {
  TString roostr = "";

  cout << endl << " ---> A Hello story ..." << endl << endl;

  int it14 = 0;
  int idm2 = 0;
  int inumToys = 0;
  double ig2 = 0;

  bool flag_verbose = false;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-it14") == 0) {
      stringstream convert(argv[i + 1]);
      if (!(convert >> it14)) {
        cerr << " ---> Error it14 !" << endl;
        exit(1);
      }
    }
    if (strcmp(argv[i], "-idm2") == 0) {
      stringstream convert(argv[i + 1]);
      if (!(convert >> idm2)) {
        cerr << " ---> Error idm2 !" << endl;
        exit(1);
      }
    }
    if (strcmp(argv[i], "-ig2") == 0) {
      stringstream convert(argv[i + 1]);
      if (!(convert >> ig2)) {
        cerr << " ---> Error int num g2 !" << endl;
        exit(1);
      }
    }
    if (strcmp(argv[i], "-v") == 0) {
      flag_verbose = true;
    }
    if (strcmp(argv[i], "-numToys") == 0) {
      stringstream convert(argv[i + 1]);
      if (!(convert >> inumToys)) {
        cerr << " ---> Error inumToys !" << endl;
        exit(1);
      }
    }
  }

  // Presave input path -- declared here (rather than at its point of use
  // further down) so it's available for the -v diagnostic dump before any
  // expensive setup.
  // TString xpath =
  // "input/current-xpan-presave_3v_hypothesis_toydata_01_cv.root";
  TString xpath = "jpmendez_presave_3v_hypothesis_toydata_01_cv.root";


  // --------------------------------------------------
  // Skip-if-already-done check: this is the cheapest possible point to bail
  // out, before touching TOsc, the presave file, or any ROOT object at all.
  // This is what lets GNU parallel farm the full (idm2, it14) grid across many
  // nodes without any node needing to know which grid points other nodes have
  // already finished.
  // --------------------------------------------------

  int display = 0;
  if (!display) {
    gROOT->SetBatch(1);
  }

  TApplication theApp("theApp", &argc, argv);

  // --------------------------------------------------
  // Grid values, computed exactly once and reused by both stages below.
  // Must happen before any TFile is opened: ROOT attaches freshly `new`'d
  // TH1D/TTree objects to whatever directory is currently open (gDirectory),
  // so building these histograms before the output .tmp file exists avoids
  // them being silently swept into it.
  // --------------------------------------------------
  const int NUM_dm2 = 60;
  const int NUM_ttt = 60;
  const double DM2_LO = -1, DM2_HI = 2;
  const double TTT_LO = -3, TTT_HI = 0;

  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, DM2_LO, DM2_HI);
  TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, TTT_LO, TTT_HI);

  int ittt = it14;
  double pars_3v_small[4] = {0, 0.10, 0.11, 0};
  double val_obj_dm2 = h1d_dm2->GetBinCenter(idm2);
  val_obj_dm2 = pow(10, val_obj_dm2);
  double val_obj_ttt = h1d_ttt->GetBinCenter(ittt);
  val_obj_ttt = pow(10, val_obj_ttt);
  double val_g2 = ig2 * M_PI;
  double pars_4v_grid[4] = {val_obj_dm2, val_obj_ttt, 0, val_g2};



  // --------------------------------------------------
  // TOsc setup (identical in both original files) -- one-time per process,
  // shared by stage 1 and stage 2 so they operate on exactly the same state.
  // --------------------------------------------------
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

  // --------------------------------------------------
  // Load the presave toydata file before running the expensive toy loops below,
  // so a missing/corrupt presave file fails fast instead of after burning
  // ~40,000 FCN evaluations. The number of entries here (N) is what determines
  // the size of every output vector.
  // --------------------------------------------------

  // int bins_eff = 181;
  // // setting Up the bins and spectrum mask
  int bins_all = 26 * 14;
  int bins_eff = 26 * 7;

  TMatrixD matrix_gof_trans_eff( bins_all, bins_eff );// oldworld, newworld
  for( int ibin=1; ibin<=bins_eff;  ibin++) matrix_gof_trans_eff(ibin-1, ibin - 1) = 1;

  //////////////////// for 3v: 3v_asmiov_pred, 3v_total_COV_inv, 3vToy
  map<int, TMatrixD> map_matrix_spectrum_3vToy; // index begins at 1
  TMatrixD matrix_3v_asimov_pred(1, bins_eff);
  TMatrixD matrix_3v_total_COV_inv(bins_eff, bins_eff);

  // Setting up a 3v null-osc Generation
  osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1],
                                 pars_3v_small[2], pars_3v_small[3], 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT(); // meas, CV, COV: all ready
  osc_test->Set_meas2fitdata();
  osc_test->FCN_Pearson_FCnew(pars_3v_small);


  // Save prediction and inverted chi2 out of tosc internal to local variables
  matrix_3v_asimov_pred = osc_test->matrix_tosc_chi2_pred;
  matrix_3v_total_COV_inv = osc_test->matrix_tosc_chi2_COV_both_syst_stat_inv;
  const int num_toys = inumToys;
  osc_test->Set_toy_variations(num_toys);
  for (int i = 0; i < num_toys; i++) {
    // Resize the matrix to the correct nunmber of bins
    map_matrix_spectrum_3vToy[i].ResizeTo(1, bins_eff);
    // Apply the mask to save the spectrum portion we are interested in
    map_matrix_spectrum_3vToy[i] =
      osc_test->map_matrix_tosc_toy_pred[i + 1] * matrix_gof_trans_eff;

    if (i % 100 == 0) {
      cout << "Finished 3v toy " << i << "\n";
        }
  }

  //////////////////// for 4v: 4v_asmiov_pred, 4v_total_COV_inv, 4vToy
  map<int, TMatrixD> map_matrix_spectrum_4vToy; // index begins at 1
  TMatrixD matrix_4v_asimov_pred(1, bins_eff);
  TMatrixD matrix_4v_total_COV_inv(bins_eff, bins_eff);


  // Setting up a 4v null-osc Generation
  // osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1],
  //                                pars_4v_grid[2], pars_4v_grid[3], 0);
  osc_test->Set_oscillation_pars(pars_3v_small[0], pars_3v_small[1],
                                 pars_3v_small[2], pars_3v_small[3], 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT(); // meas, CV, COV: all ready
  osc_test->Set_meas2fitdata();
  osc_test->FCN_Pearson_FCnew(pars_4v_grid);

  // Save prediction and inverted chi2 out of tosc internal to local variables
  matrix_4v_asimov_pred = osc_test->matrix_tosc_chi2_pred;
  matrix_4v_total_COV_inv = osc_test->matrix_tosc_chi2_COV_both_syst_stat_inv;
TString mat_file_name = TString::Format(
                                     "4v_asimov_inv_decay_BNB_grid_60x60_g2_dm2_ttt_%.2f_%03d_%03d.root", ig2, idm2, it14);
TFile specfile(mat_file_name, "RECREATE");
specfile.cd();


matrix_4v_asimov_pred.Write("mat_4v_asimov");
specfile.Close();

  cout << " ---> Finished successfully" << endl;

  cout << endl;
  if (display) {
    cout << " Enter Ctrl+c to end the program" << endl;
    cout << " Enter Ctrl+c to end the program" << endl;
    cout << endl;
    theApp.Run();
  }

  return 0;
}
