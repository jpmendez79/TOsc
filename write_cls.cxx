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
#include<algorithm>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A Hello story ..."<<endl<<endl;
  // Pre-store the 3v/4v distributions

  std::vector<std::vector<double>> vec4v_cache;
  std::vector<std::vector<double>> vec3v_cache;
  // Number of files
  vec4v_cache.resize(3600);
  vec3v_cache.resize(3600);
  int N_it14 = 60;
  int N_idm2 = 60;
  for (int idm2 = 0; idm2 < 60; ++idm2) {
    for (int it14 = 0; it14 < 60; ++it14) {

      int idx = idm2 * N_it14 + it14;

      TString roostr = TString::Format("output-6k-size/out_dm2_ttt_%03d_%03d.root", idm2+1, it14+1);

      TFile file(roostr, "READ");
      TTree* tree = (TTree*)file.Get("tree");

      std::vector<double>* invec4v = nullptr;
      std::vector<double>* invec3v = nullptr;

      tree->SetBranchAddress("vec_dchi2_4v", &invec4v);
      tree->SetBranchAddress("vec_dchi2_3v", &invec3v);

      tree->GetEntry(0);

      vec4v_cache[idx] = *invec4v;
      vec3v_cache[idx] = *invec3v;

      // file closes automatically here (scope end)
    }
  }
  for (auto& v : vec4v_cache) std::sort(v.begin(), v.end());
  for (auto& v : vec3v_cache) std::sort(v.begin(), v.end());
  cout << "Finish Presave all out*.root files\n";




  // A vectorized version
  std::vector<double> cls_grid(60 * 60);

  // input matrix
  TMatrixT<double>* inmat = nullptr;

  // cache already sorted BEFORE this loop

  TFile obsfile("input/presave_deltachisquare_obs.root", "READ");
  TTree* obstree = (TTree*)obsfile.Get("tree");
  obstree->SetBranchAddress("obs_map", &inmat);

  TFile file("cls_maps.root", "RECREATE");
  TTree tree("tree", "CLs Grids");
  tree.Branch("cls_grid", &cls_grid);

  for (int universe = 0; universe < 50; ++universe) {

    obstree->GetEntry(universe);

    const TMatrixT<double>& obs = *inmat;

    // direct pointer access (faster than operator[])
    for (int idm2 = 0; idm2 < N_idm2; ++idm2) {
      for (int it14 = 0; it14 < N_it14; ++it14) {
        double refval = (*inmat)(idm2, it14);
        int idx = idm2 * N_it14 + it14;

        const std::vector<double>& vec4v = vec4v_cache[idx];
        const std::vector<double> &vec3v = vec3v_cache[idx];
        auto it4 = std::lower_bound(vec4v.data(), vec4v.data() + vec4v.size(), refval);

        auto it3 = std::lower_bound(vec3v.data(),
                                    vec3v.data() + vec3v.size(),
                                    refval);

        double count4v = (vec4v.data() + vec4v.size()) - it4;
        double count3v = (vec3v.data() + vec3v.size()) - it3;
	    double cls = 0;
	    if( count3v == 0 ) {
	      if( count4v == 0 ) cls = 0;
	      else cls = 1;
	    }
	    else cls = count4v / count3v;
	    if( count4v>=count3v ) cls = 1;

        cls_grid[idx] = 1.0 - cls;
      }
    }

    tree.Fill();

    if (universe % 5 == 0)
      std::cout << "Finished universe " << universe << "\n";
  }

  obsfile.Close();
  tree.Write();
  file.Close();
  return 0;
  cout << " ---> Finished successfully" << endl;

  cout<<endl;


  return 0;
}
