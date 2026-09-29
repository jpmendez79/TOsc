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

  TApplication theApp("theApp", &argc, argv);


  // Define grid used for coarse grid scan. This should remain the same as the CLs calculation testing
  const int NUM_dm2 = 60;
  const int NUM_ttt = 60;
  const int NUM_g2 = 60;
  const double DM2_LO = -1, DM2_HI = 2;
  const double TTT_LO = -3, TTT_HI = 0;
  const double G2_LO = 0, G2_HI = 4 * M_PI;

  double test_val_dm2[NUM_dm2];
  double test_val_ttt[NUM_ttt];
  double test_val_g2[NUM_g2];

  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, DM2_LO, DM2_HI);
  TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, TTT_LO, TTT_HI);
  TH1D *h1d_g2 = new TH1D("h1d_g2", "h1d_g2", NUM_g2, G2_LO, G2_HI);

  for (int i = 0; i < NUM_dm2; i++) {
    double val_obj_dm2 = h1d_dm2->GetBinCenter(i);
    test_val_dm2[i] = pow(10, val_obj_dm2);
    double val_obj_ttt = h1d_ttt->GetBinCenter(i);
    test_val_ttt[i] = pow(10, val_obj_ttt);
    test_val_g2[i] = h1d_g2->GetBinCenter(i);
  }


  // For testing Purposes Create an measurement to be used
  double seed_val_dm2 = 7.3;
  double seed_val_ttt = 0.23;
  double seed_val_g2 = 2.5 * M_PI;
  double pars_test_meas[3] = {seed_val_dm2, seed_val_ttt, seed_val_g2};



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



  // Generate the test measurement Asimov
  osc_test->Set_oscillation_pars(pars_test_meas[0], pars_test_meas[1], 0, 0, pars_test_meas[2]);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT(); // meas, CV, COV: all ready
  osc_test->Set_asimov2fitdata(); // Set the asimov generated prediction as the test prediction


  // Variables needed for fitting
  double min_chi2 = 1e6;
  double min_ttt = -1;
  double min_dm2 = -1;
  double min_g2 = -1;
  // progress tracking
  double dimensions = NUM_dm2 * NUM_ttt * NUM_g2;
  int percent_target = 0.10 * dimensions; // Report every 10%
  int counter = 0;
  int percent_counter = 0;
  cout << "Starting IVP Coarse Grid Scan \n";
  // Record starting time
    auto start =
        high_resolution_clock::now();
  for (int idm2 = 0; idm2 < NUM_dm2; idm2++) {
    for (int ittt = 0; ittt < NUM_ttt; ittt++) {
      for (int ig2 = 0; ig2 < NUM_g2; ig2++) {
        double test_par[4] = {test_val_dm2[idm2], test_val_ttt[ittt], 0, test_val_g2[ig2]};
        double test_chi2 = osc_test->FCN(test_par);

        if (test_chi2 < min_chi2) {
          min_chi2 = test_chi2;
          min_ttt = test_par[1];
          min_dm2 = test_par[0];
          min_g2 = test_par[3];
        }
        counter++;
        if ((counter % percent_target) == 0) {
          percent_counter += 10;
          cout << percent_counter << "% Done \n";

        }
      }
    }
  }

  // Record ending time
  auto stop =
    high_resolution_clock::now();

 auto duration =
        duration_cast<microseconds>(
            stop - start);

    cout << "Time taken: "
         << duration.count()
         << " microseconds";

    cout << "Min Chi2: " << min_chi2 << endl;
    cout << "(ttt, dm2, g2)" << "(" << min_ttt << "," << min_dm2 << ","
    << min_g2 << ")" << endl;
    cout << "Actual Parameters" << endl;
    cout << "(";
    for (int i = 0; i < 3; i++) {
      cout << pars_test_meas[i];
      if (i < 2) cout << ",";
    }
    cout << endl;

  if (display) {
    cout << " Enter Ctrl+c to end the program" << endl;
    cout << " Enter Ctrl+c to end the program" << endl;
    cout << endl;
    theApp.Run();
  }

  return 0;
}
