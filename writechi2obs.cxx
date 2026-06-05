#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<filesystem>
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
  int it14 = 0;
  int idm2 = 0;
  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-it14")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it14 ) ) { cerr<<" ---> Error it14 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-idm2")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>idm2 ) ) { cerr<<" ---> Error idm2 !"<<endl; exit(1); }
    }
  }

TString final_name = TString::Format("output/size_1k_cls_out_dm2_ttt_%03d_%03d.root", idm2, it14);

TString tmp_name = final_name + ".tmp";

if (std::filesystem::exists(final_name.Data())) {
    return 0;
}


//
// Write Output tmp

TFile outfile(tmp_name, "RECREATE");


  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;

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


  // Create arrays to hold the values I need
  double pars_3v_small[4] = {0, 0.10, 0.11, 0};

//


// Open presaved file
// --------------------------------------------------
// Open input file
// --------------------------------------------------
TString xpath = "xiangpan-presave.root";
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
vector<double> vec_obs_3v;
vector<double> vec_obs_4v;
vector<double> vec_dchi2obs;
vector<double> vec_lines_3v;
vector<double> vec_lines4v;
vector<double> vec_confidence;
// --------------------------------------------------
// Index map (cache)
// --------------------------------------------------
vector<vector<double>> spectrum_cache;
// spectrum_cache.reserve(200);  // optional optimization

while (reader.Next()) {
    spectrum_cache.push_back(*vec_toydata_spectrum);
}

cout << "Cached entries: " << spectrum_cache.size() << endl;
  const int N = spectrum_cache.size();
  vec_obs_3v.resize(N);
vec_obs_4v.resize(N);
vec_dchi2obs.resize(N);
vec_lines_3v.resize(N);
vec_lines4v.resize(N);
vec_confidence.resize(N);

// --------------------------------------------------
// Loop over entries
// --------------------------------------------------
// Save all info needed from index file
roostr = TString::Format("output/out_dm2_ttt_%03d_%03d.root", idm2, it14);
cout << roostr << endl;
TFile file(roostr, "READ");



TTree* tree = (TTree*)file.Get("tree");

std::vector<double>* invec4v = nullptr;
std::vector<double>* invec3v = nullptr;

double val_obj_dm2;
double val_obj_ttt;

tree->SetBranchAddress("vec_dchi2_4v", &invec4v);
tree->SetBranchAddress("vec_dchi2_3v", &invec3v);
tree->SetBranchAddress("data_val_dm2", &val_obj_dm2);
tree->SetBranchAddress("data_val_ttt", &val_obj_ttt);

tree->GetEntry(0);

// // COPY OUT OF ROOT OWNERSHIP
// std::vector<double> vec4v = *invec4v;
// std::vector<double> vec3v = *invec3v;

const int size = 1000;/* number of elements to copy */

// COPY OUT OF ROOT OWNERSHIP
std::vector<double> vec4v(invec4v->begin(), invec4v->begin() + size);
std::vector<double> vec3v(invec3v->begin(), invec3v->begin() + size);

// safe to close ROOT now
file.Close();

// now pure C++
std::sort(vec4v.begin(), vec4v.end());
std::sort(vec3v.begin(), vec3v.end());

// 2000-loop here

double pars_4v_grid[4] ={val_obj_dm2, val_obj_ttt, 0.0045, 0};

      for (int u = 0; u < N; u++) {
    // Load the spectrum
    auto& spectrum = spectrum_cache[u];
    for (int j = 0; j <= 181; j++) {
      osc_test->matrix_tosc_fitdata_newworld(0, j) = spectrum[j];
    }
    // Calculate chi2_3v
    double obs_3v = osc_test->FCN(pars_3v_small);

    // Calculate chi2






    double obs_4v = osc_test->FCN(pars_4v_grid);

    // Calculate deltachi2obs
    double deltachi2obs = obs_4v - obs_3v;

    // Calculate the tail
    auto it4 = std::lower_bound(vec4v.data(), vec4v.data() + vec4v.size(), deltachi2obs);
    auto it3 = std::lower_bound(vec3v.data(), vec3v.data() + vec3v.size(), deltachi2obs);
    double count4v = (vec4v.data() + vec4v.size()) - it4;
    double count3v = (vec3v.data() + vec3v.size()) - it3;

    // Calculate the CLs
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

    if (u%100 == 0) cout << "Finished Universe " << u <<endl;

      } // Universe
      outfile.cd();
  TTree outtree("tree", "CLs Grids");
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

  return 0;
  }
