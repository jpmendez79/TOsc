#include<iostream>
#include<sstream>
#include<cmath>
#include<filesystem>
#include "stdlib.h"
using namespace std;

#include<vector>

#include "WCPLEEANA/TOsc.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Combines the two-stage create_dists.cxx -> writechi2obs.cxx pipeline into a single
// per-grid-point program. Stage 1 (toy Delta-chi2 generation) and stage 2 (observed
// CLs/confidence calculation) now share one process and one TOsc instance, so the
// intermediate ROOT file that used to connect them is no longer written at all.
//
int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A Hello story ..."<<endl<<endl;

  int it14 = 0;
  int idm2 = 0;
  int inumXgrids = 0;
  int inumYgrids = 0;
  double iparamXmin = 0;
  double iparamXmax = 0;
  double iparamYmin = 0;
  double iparamYmax = 0;
  int  inumToys = 0;
  bool flag_verbose = false;
  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-it14")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it14 ) ) { cerr<<" ---> Error it14 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-idm2")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>idm2 ) ) { cerr<<" ---> Error idm2 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-v")==0 ) {
      flag_verbose = true;
    }
    if( strcmp(argv[i],"-numToys")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>inumToys ) ) { cerr<<" ---> Error inumToys !"<<endl; exit(1); }
    }
  }

  // Presave input path -- declared here (rather than at its point of use further
  // down) so it's available for the -v diagnostic dump before any expensive setup.
  TString xpath = "input/presave_3v_hypothesis_toydata_01_cv.root";

  // --------------------------------------------------
  // Skip-if-already-done check: this is the cheapest possible point to bail out,
  // before touching TOsc, the presave file, or any ROOT object at all. This is what
  // lets GNU parallel farm the full (idm2, it14) grid across many nodes without any
  // node needing to know which grid points other nodes have already finished.
  // --------------------------------------------------
  TString final_name = TString::Format("output/inv_decay_BNB_grid_60x60_dm2_ttt_%03d_%03d.root", idm2, it14);
  TString tmp_name = final_name + ".tmp";

  if (std::filesystem::exists(final_name.Data())) {
    cout << " ---> " << final_name << " already exists, skipping." << endl;
    return 0;
  }

  // output/ is gitignored and absent on a fresh checkout; create_directories is a
  // safe no-op if it already exists (including when raced by concurrent workers).
  std::filesystem::create_directories("output");

  int display = 0;
  if( !display ) {
    gROOT->SetBatch( 1 );
  }

  TApplication theApp("theApp",&argc,argv);

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
  const double TTT_LO = -1, TTT_HI = 1;

  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "h1d_dm2", NUM_dm2, DM2_LO, DM2_HI);
  TH1D *h1d_ttt = new TH1D("h1d_ttt", "h1d_ttt", NUM_ttt, TTT_LO, TTT_HI);

  if( flag_verbose ) {
    cout << "\n ---> Verbose configuration:\n";
    cout << "  input presave file     : " << xpath << "\n";
    cout << "  final output file      : " << final_name << "\n";
    cout << "  h1d_dm2 bounds [bins]   : [" << DM2_LO << ", " << DM2_HI << "] (" << NUM_dm2 << " bins)\n";
    cout << "  h1d_ttt bounds [bins]   : [" << TTT_LO << ", " << TTT_HI << "] (" << NUM_ttt << " bins)\n";
    cout << "  default_cv_file         : " << Configure_Osc::default_cv_file << "\n";
    cout << "  default_dirtadd_file    : " << Configure_Osc::default_dirtadd_file << "\n";
    cout << "  default_mcstat_file     : " << Configure_Osc::default_mcstat_file << "\n";
    cout << "  default_fluxXs_dir      : " << Configure_Osc::default_fluxXs_dir << "\n";
    cout << "  default_detector_dir    : " << Configure_Osc::default_detector_dir << "\n";
    cout << "  default_eventlist_dir   : " << Configure_Osc::default_eventlist_dir << "\n";
    cout << endl;
  }

  int ittt = it14;
  double pars_3v_small[4] = {0, 0.10, 0.11, 0};
  double val_obj_dm2 = h1d_dm2->GetBinCenter(idm2);
  val_obj_dm2 = pow(10, val_obj_dm2);
  double val_obj_ttt = h1d_ttt->GetBinCenter(ittt);
  val_obj_ttt = pow(10, val_obj_ttt);
  double pars_4v_grid[4] = {val_obj_dm2, val_obj_ttt, 0.0045, 0};

  // --------------------------------------------------
  // TOsc setup (identical in both original files) -- one-time per process,
  // shared by stage 1 and stage 2 so they operate on exactly the same state.
  // --------------------------------------------------
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;

  TOsc *osc_test = new TOsc();

  osc_test->tosc_scaleF_POT_BNB  = scaleF_POT_BNB;
  osc_test->tosc_scaleF_POT_NuMI = scaleF_POT_NuMI;

  osc_test->flag_apply_oscillation_BNB  = Configure_Osc::flag_apply_oscillation_BNB;
  osc_test->flag_apply_oscillation_NuMI = Configure_Osc::flag_apply_oscillation_NuMI;

  osc_test->flag_goodness_of_fit_CNP    = Configure_Osc::flag_goodness_of_fit_CNP;

  osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;

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

  // --------------------------------------------------
  // Load the presave toydata file before running the expensive toy loops below,
  // so a missing/corrupt presave file fails fast instead of after burning
  // ~40,000 FCN evaluations. The number of entries here (N) is what determines
  // the size of every output vector.
  // --------------------------------------------------
  TFile* inputfile_toydata_cv = TFile::Open(xpath, "READ");

  if (!inputfile_toydata_cv || inputfile_toydata_cv->IsZombie()) {
      std::cerr << "Error: cannot open input file\n";
      return 1;
  }

  TTree* tree_toydata = (TTree*)inputfile_toydata_cv->Get("tree_toydata");

  if (!tree_toydata) {
      std::cerr << "Error: tree_toydata not found\n";
      return 1;
  }

  TTreeReader reader(tree_toydata);
  TTreeReaderValue<vector<double>> vec_toydata_spectrum(reader, "vec_toydata_spectrum");

  vector<vector<double>> spectrum_cache;
  while (reader.Next()) {
      spectrum_cache.push_back(*vec_toydata_spectrum);
  }
  inputfile_toydata_cv->Close();

  cout << "Cached entries: " << spectrum_cache.size() << endl;
  const int N = spectrum_cache.size();

  vector<double> vec_obs_3v(N);
  vector<double> vec_obs_4v(N);
  vector<double> vec_dchi2obs(N);
  vector<double> vec_lines_3v(N);
  vector<double> vec_lines4v(N);
  vector<double> vec_confidence(N);

  // --------------------------------------------------
  // Output .tmp file, held open for the remainder of the run. No further
  // `new TH1D`/`new TTree` calls happen until the explicit outfile.cd() below,
  // so nothing gets silently attached to it in the meantime.
  // --------------------------------------------------
  TFile outfile(tmp_name, "RECREATE");

  // --------------------------------------------------
  // Stage 1 (from create_dists.cxx): generate toy Delta-chi2 distributions
  // under the 3v and 4v hypotheses. Only the Delta-chi2 = chi2_4v - chi2_3v
  // values are kept -- the individual chi2_3v/chi2_4v-per-toy vectors that
  // create_dists.cxx also wrote are never read back by stage 2, so they are
  // not stored here at all.
  //
  // Set_toy_variations() re-decomposes the systematic covariance matrix
  // (eigen-decomposition) that Apply_oscillation()/Set_apply_POT() just
  // rebuilt for the hypothesis in question. That covariance is different for
  // the 3v vs. 4v hypothesis and different at every grid point, so this work
  // is NOT cacheable across hypotheses or across grid points -- it is
  // genuinely new computation each time, not a redundancy from the old
  // two-program split.
  // --------------------------------------------------
  const int num_toys = inumToys;

  vector<double> vec_chi2_3vToy_3v;
  vector<double> vec_chi2_3vToy_4v;
  vector<double> vec_chi2_4vToy_3v;
  vector<double> vec_chi2_4vToy_4v;
  vector<double> vec_dchi2_4v;
  vector<double> vec_dchi2_3v;
  vec_dchi2_3v.reserve(num_toys);
  vec_dchi2_4v.reserve(num_toys);
  vec_chi2_3vToy_3v.reserve(num_toys);
  vec_chi2_3vToy_4v.reserve(num_toys);
  vec_chi2_4vToy_3v.reserve(num_toys);
  vec_chi2_4vToy_4v.reserve(num_toys);


  cout << "Generate 3v\n";
  osc_test->Set_oscillation_pars(0, 0.10, 0.11, 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  osc_test->Set_toy_variations(num_toys);
  for (int i = 0; i < num_toys; i++) {
    osc_test->Set_toy2fitdata(i+1);
    double chi2_3v = osc_test->FCN(pars_3v_small);
    double chi2_4v = osc_test->FCN(pars_4v_grid);
    vec_chi2_3vToy_3v.push_back(chi2_3v);
    vec_chi2_3vToy_4v.push_back(chi2_4v);
    vec_dchi2_3v.push_back(chi2_4v - chi2_3v);
  }

  cout << "Generate 4v\n";
  osc_test->Set_oscillation_pars(pars_4v_grid[0], pars_4v_grid[1], pars_4v_grid[2], pars_4v_grid[3]);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  osc_test->Set_toy_variations(num_toys);
  for (int i = 0; i < num_toys; i++) {
    osc_test->Set_toy2fitdata(i+1);
    double chi2_3v = osc_test->FCN(pars_3v_small);
    double chi2_4v = osc_test->FCN(pars_4v_grid);
    vec_chi2_4vToy_3v.push_back(chi2_3v);
    vec_chi2_4vToy_4v.push_back(chi2_4v);
    vec_dchi2_4v.push_back(chi2_4v - chi2_3v);
  }

  // Sorted in place and used directly by stage 2 below -- no intermediate
  // ROOT file write/read/copy round-trip needed since these are already
  // process-owned std::vectors.
  std::sort(vec_dchi2_3v.begin(), vec_dchi2_3v.end());
  std::sort(vec_dchi2_4v.begin(), vec_dchi2_4v.end());

  // --------------------------------------------------
  // Stage 2 (from writechi2obs.cxx): for every presave universe, compute the
  // observed Delta-chi2 and derive CLs/confidence against the toy
  // distributions generated above.
  // --------------------------------------------------
  for (int u = 0; u < N; u++) {
    auto& spectrum = spectrum_cache[u];
    for (int j = 0; j <= 181; j++) {
      osc_test->matrix_tosc_fitdata_newworld(0, j) = spectrum[j];
    }

    double obs_3v = osc_test->FCN(pars_3v_small);
    double obs_4v = osc_test->FCN(pars_4v_grid);
    double deltachi2obs = obs_4v - obs_3v;

    auto it4 = std::lower_bound(vec_dchi2_4v.data(), vec_dchi2_4v.data() + vec_dchi2_4v.size(), deltachi2obs);
    auto it3 = std::lower_bound(vec_dchi2_3v.data(), vec_dchi2_3v.data() + vec_dchi2_3v.size(), deltachi2obs);
    double count4v = (vec_dchi2_4v.data() + vec_dchi2_4v.size()) - it4;
    double count3v = (vec_dchi2_3v.data() + vec_dchi2_3v.size()) - it3;

    double cls = 0;
    if( count3v == 0 ) {
      if( count4v == 0 ) cls = 0;
      else cls = 1;
    }
    else cls = count4v / count3v;
    if( count4v>=count3v ) cls = 1;
    double confidence = 1.0 - cls;

    vec_obs_3v[u] = obs_3v;
    vec_obs_4v[u] = obs_4v;
    vec_dchi2obs[u] = deltachi2obs;
    vec_lines_3v[u] = count3v;
    vec_lines4v[u] = count4v;
    vec_confidence[u] = confidence;

    if (u%100 == 0) cout << "Finished Universe " << u << endl;
  } // Universe

  // --------------------------------------------------
  // Write final output. Explicit cd() so the new TTree attaches to outfile
  // and not to whatever directory ROOT last touched.
  // --------------------------------------------------
  outfile.cd();
  TTree outtree("tree", "CLs Grids");

  outtree.Branch("grid_idx_dm2", &idm2, "idm2/I");
  outtree.Branch("grid_idx_ttt", &ittt, "ittt/I");
  outtree.Branch("data_val_dm2", &val_obj_dm2, "val_obj_dm2/D");
  outtree.Branch("data_val_ttt", &val_obj_ttt, "val_obj_ttt/D");

  outtree.Branch("vec_chi2_4vToy_4v", &vec_chi2_4vToy_4v);
  outtree.Branch("vec_chi2_4vToy_3v", &vec_chi2_4vToy_3v);
  outtree.Branch("vec_dchi2_4v", &vec_dchi2_4v);

  outtree.Branch("vec_chi2_3vToy_4v", &vec_chi2_4vToy_4v);
  outtree.Branch("vec_chi2_3vToy_3v", &vec_chi2_4vToy_3v);
  outtree.Branch("vec_dchi2_3v", &vec_dchi2_3v);

  outtree.Branch("vec_obs_3v", &vec_obs_3v);
  outtree.Branch("vec_obs_4v", &vec_obs_4v);
  outtree.Branch("vec_dchi2obs", &vec_dchi2obs);
  outtree.Branch("vec_lines_3v", &vec_lines_3v);
  outtree.Branch("vec_lines4v", &vec_lines4v);
  outtree.Branch("vec_confidence", &vec_confidence);

  outtree.Fill();
  outtree.Write();
  outfile.Close();

  std::filesystem::rename(tmp_name.Data(), final_name.Data());
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
