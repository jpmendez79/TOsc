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

  TApplication theApp("theApp", &argc, argv);

  int seq = std::atoi(argv[1]);

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

    // For testing Purposes Create an measurement to be used
  double seed_val_dm2 = 7.3;
  double seed_val_ttt = 0.23;
  double seed_val_g2 = 2.5 * M_PI;
  double pars_test_meas[3] = {seed_val_dm2, seed_val_ttt, seed_val_g2};

  // Generate the test measurement Asimov
  osc_test->Set_oscillation_pars(pars_test_meas[0], pars_test_meas[1], 0, 0, pars_test_meas[2]);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT(); // meas, CV, COV: all ready
  osc_test->Set_asimov2fitdata(); // Set the asimov generated prediction as the test prediction

  // Calculate chi2
  double test_par[4] = {test_dm2, test_ttt, 0, test_g2};
  double test_chi2 = osc_test->FCN(test_par);

  // Save Result to file
  std::ostringstream line;
  line << std::setprecision(std::numeric_limits<double>::max_digits10);

  line << seq << " " << test_chi2 << " " << ittt << " " << test_ttt << " "
  << idm2 << " " << test_dm2 << " " << ig2 << " " << test_g2 << "\n";



  std::ofstream out("output/ipv_" + std::to_string(seq) + ".txt");
  out << line.str() << std::flush;

  if (display) {
    cout << " Enter Ctrl+c to end the program" << endl;
    cout << " Enter Ctrl+c to end the program" << endl;
    cout << endl;
    theApp.Run();
  }

  return 0;
}
