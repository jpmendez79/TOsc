
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
  double init_ttt = -1;
  double init_dm2 = -1;
  double init_g2 = -1;
  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, DM2_LO, DM2_HI);
  TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, TTT_LO, TTT_HI);
  double chi2_min = pow(10, 6);
  vector<double>vec_toydata_spectrum;
  double val_dm2_41         = 0;
  double val_sin2_2theta_14 = 0.36;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;
  TFile *outfile_toydata = new TFile("out.root", "recreate");
  /// standard order
  val_dm2_41         = 0;
  val_sin2_2theta_14 = 0.2;
  val_sin2_theta_24  = 0.3;
  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34, 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready

  osc_test->Set_meas2fitdata();
  TMatrixD data_out = osc_test->matrix_tosc_fitdata_newworld.GetSub(0,0, 0, 26*7-1);
  int rows = data_out.GetNcols();
  TH1D hdata("hdata", "Data", rows, 1, rows);
  cout << rows << endl;
  for (int idx = 0; idx < rows; idx++) {
    hdata.SetBinContent(idx+1, data_out(0, idx) );
  }

  osc_test->Set_asimov2fitdata();
  TMatrixD dec_asimov_out = osc_test->matrix_tosc_fitdata_newworld.GetSub(0,0, 0, 26*7-1);
  TH1D hdecay("hg20", "Asimov g2=0", rows, 1, rows);
  cout << rows << endl;
  for (int idx = 0; idx < rows; idx++) {
    hdecay.SetBinContent(idx+1, dec_asimov_out(0, idx) );
  }
  hdata.Write();
  hdecay.Write();
  outfile_toydata->Close();
  // vec_toydata_spectrum.clear();
  // for(int idx=0; idx<(osc_test->matrix_tosc_fitdata_newworld.GetNcols()); idx++) {
  //   vec_toydata_spectrum.push_back( osc_test->matrix_tosc_fitdata_newworld(0, idx) );
  // }


  // Coarse Grid scan
  cout << "Grid Scan" << endl;
  osc_test->Set_meas2fitdata();
  for (int ittt = 0; ittt < NUM_ttt; ittt++) {
    // Generate the itt value once per itt
    double val_obj_ttt = h1d_ttt->GetBinCenter(ittt);
    val_obj_ttt = pow(10, val_obj_ttt);
      if (ittt % 10 == 0) cout << "ttt " << ittt << endl;
    for (int idm2 = 0; idm2 < NUM_dm2; idm2++) {
      // Generate the idm2 once per dm2
      double val_obj_dm2 = h1d_dm2->GetBinCenter(idm2);
      val_obj_dm2 = pow(10, val_obj_dm2);
      for (int ig2 = 0; ig2 < 9; ig2++) {
        // Generate the g2 value
        double val_obj_g2 = ig2 * 0.5 * M_PI;
        // Generate the 4nu prediction
        double pars_4v_grid[5] = {val_obj_dm2, val_obj_ttt, 0.0045, 0, val_obj_g2};
        double chi2 = osc_test->FCN(pars_4v_grid);
        if (chi2 < chi2_min) {
          chi2_min = chi2;
          init_ttt = val_obj_ttt;
          init_dm2 = val_obj_dm2;
          init_g2 = val_obj_g2;
          cout << "Minimum chi2 " << chi2_min << endl;
          cout << "Parameter itt: " << val_obj_ttt << endl;
          cout << "Parameter dm2: " << val_obj_dm2 << endl;
          cout << "Parameter g2: " << val_obj_g2 << endl;
        }

        //
      }
    }
  }

  cout << "Minimum chi2 " << chi2_min << endl;
  cout << "Parameter itt: " << init_ttt << endl;
  cout << "Parameter dm2: " << init_dm2 << endl;
  cout << "Parameter g2: " << init_g2 << endl;



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
