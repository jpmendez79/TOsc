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

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"
#include "TParameter.h"
#include <chrono> // timer

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A Hello story ..."<<endl<<endl;




  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;


  int it14 = 0;
  int idm2 = 0;


  const int NUM_dm2 = 60;
  const int NUM_ttt = 60;

  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, -2, 1);
  TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, -2, 0);



  // double it14 = 0;
  // double idm2 = 0;
  // double it24 = 0;


  if( !display ) {
    gROOT->SetBatch( 1 );
  }

  TApplication theApp("theApp",&argc,argv);

  /////////////////////////////////////////////////////////// Draw style


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
  /////// set only one time

  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
                               Configure_Osc::default_dirtadd_file,
                               Configure_Osc::default_mcstat_file,
                               Configure_Osc::default_fluxXs_dir,
                               Configure_Osc::default_detector_dir);
  osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);


  cout << "Read in xiangpan file \n";
  idm2 = 41;
  it14 = 30;
  // Create arrays to hold the values I need
  double pars_3v_small[4] = {0, 0.10, 0.11, 0};
  double val_obj_dm2 = h1d_dm2->GetBinCenter(idm2);
  val_obj_dm2 = pow(10, val_obj_dm2);
  double val_obj_ttt = h1d_ttt->GetBinCenter(it14);
  val_obj_ttt = pow(10, val_obj_ttt);
    double pars_4v_grid[4] ={val_obj_dm2, val_obj_ttt, 0.0045, 0};

// --------------------------------------------------
// Open input file
// --------------------------------------------------
TString xpath = "input/presave_3v_hypothesis_toydata_01_cv.root";
TFile* inputfile_toydata_cv = TFile::Open(xpath, "READ");

if (!inputfile_toydata_cv || inputfile_toydata_cv->IsZombie()) {
    std::cerr << "Error: cannot open input file\n";
    return 1;
}

// --------------------------------------------------
// Get tree
// --------------------------------------------------
TTree* tree_toydata =
    (TTree*)inputfile_toydata_cv->Get("tree_toydata");

if (!tree_toydata) {
    std::cerr << "Error: tree_toydata not found\n";
    return 1;
}

// --------------------------------------------------
// TTreeReader setup (SAFE ROOT I/O)
// --------------------------------------------------
TTreeReader reader(tree_toydata);

// This replaces SetBranchAddress entirely
TTreeReaderValue<vector<double>> vec_toydata_spectrum(
    reader, "vec_toydata_spectrum"
);

// --------------------------------------------------
// Your output containers
// --------------------------------------------------
vector<double> vec_chi2_3v;

// --------------------------------------------------
// Index map (cache)
// --------------------------------------------------
vector<vector<double>> spectrum_cache;
spectrum_cache.reserve(200);  // optional optimization

while (reader.Next()) {
    spectrum_cache.push_back(*vec_toydata_spectrum);
}

cout << "Cached entries: " << spectrum_cache.size() << endl;
// --------------------------------------------------
// Loop over entries
// --------------------------------------------------
int entry_index = 0;
	  // TMatrixD matrix_delta_3vPred_data   = matrix_toydata_temp - matrix_3v_asmiov_pred;
	  // TMatrixD matrix_delta_3vPred_data_T = matrix_delta_3vPred_data.T(); matrix_delta_3vPred_data.T();
	  // double chi2_3v_with_data = (matrix_delta_3vPred_data * matrix_3v_total_COV_inv * matrix_delta_3vPred_data_T)(0,0);

for (size_t i = 0; i < 200; i++) {

    auto& spectrum = spectrum_cache[i];
    for (int j = 0; j <= 181; j++) {
        osc_test->matrix_tosc_fitdata_newworld(0, j) = spectrum[j];
    }

    // double chi2 = osc_test->FCN(pars_3v_small);
    double chi2 = osc_test->FCN(pars_4v_grid);
cout << chi2 << endl;
    vec_chi2_3v.push_back(chi2);
}

cout << "dm2 " << val_obj_dm2 << "ttt " << val_obj_ttt << endl;


  TFile outfile("output/random_test.root", "RECREATE");
  TTree *tree = new TTree("tree", "Obs Test");
  tree->Branch("vec_chi2_3v", &vec_chi2_3v);
  tree->Fill();
  tree->Write();
  outfile.Close();

  cout << " ---> Finished successfully" << endl;

  cout<<endl;
  if( display ) {
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<endl;
    theApp.Run();
  }

  return 0;
  }
