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



  double ittt = 0;
  double idm2 = 0;
  double ig2 = 0;

	for (int i = 1; i < argc; i++) {
	if (strcmp(argv[i], "-it14") == 0) {
	  stringstream convert(argv[i + 1]);
	  if (!(convert >> ittt)) {
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
	}

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


cout<<endl;
	  cout<<" ---> Median sensitivity from freqentist CLs (step_1): pre-save toys with 3v hypothesis"<<endl;
	  cout<<endl;
	  double val_dm2_41         = 0;
	  double val_sin2_2theta_14 = 0.36;
	  double val_sin2_theta_24  = 0;
	  double val_sin2_theta_34 = 0;
	  double val_g2 = ig2;

  const int NUM_dm2 = 60;
  const int NUM_ttt = 60;
  const int NUM_g2 =  60;
  const double DM2_LO = -1, DM2_HI = 2;
  const double TTT_LO = -3, TTT_HI = 0;
  const double G2_LO = -3, G2_HI = 0;


  auto center = [](int idx, int n, double lo, double hi) {
	return lo + (idx + 0.5) * (hi - lo) / n;
  };

	/// standard order
	  val_dm2_41         = idm2;
	  val_sin2_2theta_14 = ittt;
	  val_sin2_theta_24  = 0.0;
	  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34, val_g2);
	  osc_test->Apply_oscillation();
	  osc_test->Set_apply_POT();// meas, CV, COV: all ready

	  int i_toydata = 0;
	  vector<double>vec_toydata_spectrum;
	  TFile *outfile_toydata = new TFile("jpmendez_presave_spectrums.root", "recreate");
	  TTree *tree_toydata = new TTree("tree_toydata", "toydata_with_3v_hypothesis ---> i_toydata(1) for real measurement, i_toydata(2) for asimov");
	  tree_toydata->Branch("i_toydata", &i_toydata, "i_toydata/I");
	  tree_toydata->Branch("vec_toydata_spectrum", &vec_toydata_spectrum);

  for (int seq = 0; seq < 216000; seq++) {
	if(seq%10==0) cout<<TString::Format(" ---> producing toy %5d / 216000", seq)<<endl;
	int ig2 = seq % NUM_g2;
	int idm2 = (seq / NUM_g2) % NUM_dm2;
	int ittt = seq / (NUM_dm2 * NUM_g2);

	val_dm2_41 = pow(10, center(idm2, NUM_dm2, DM2_LO, DM2_HI));
	val_sin2_2theta_14 = pow(10, center(ittt, NUM_ttt, TTT_LO, TTT_HI));
	val_g2 = (10, center(ig2, NUM_g2, G2_LO, G2_HI));

	osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34, val_g2);
	osc_test->Apply_oscillation();
	osc_test->Set_apply_POT();// meas, CV, COV: all ready

	osc_test->Set_asimov2fitdata();
	vec_toydata_spectrum.clear();
	for(int idx=0; idx<(osc_test->matrix_tosc_fitdata_newworld.GetNcols()); idx++) {
	  vec_toydata_spectrum.push_back( osc_test->matrix_tosc_fitdata_newworld(0, idx) );
	}
	i_toydata++;
	tree_toydata->Fill();

  }


	  // /////// the first "toydata" is the real measurement
	  // osc_test->Set_meas2fitdata();
	  // vec_toydata_spectrum.clear();
	  // for(int idx=0; idx<(osc_test->matrix_tosc_fitdata_newworld.GetNcols()); idx++) {
	  //	vec_toydata_spectrum.push_back( osc_test->matrix_tosc_fitdata_newworld(0, idx) );
	  // }
	  // i_toydata++;
	  // tree_toydata->Fill();

	  // /////// the second "toydata" is the asimov dataset
	  // osc_test->Set_asimov2fitdata();
	  // vec_toydata_spectrum.clear();
	  // for(int idx=0; idx<(osc_test->matrix_tosc_fitdata_newworld.GetNcols()); idx++) {
	  //	vec_toydata_spectrum.push_back( osc_test->matrix_tosc_fitdata_newworld(0, idx) );
	  // }
	  // i_toydata++;
	  // tree_toydata->Fill();

	  // ///////

	  // int NUM_TOYS = 2000;// toys to be used to calculate the median result, "ccb"

	  // cout<<endl;
	  // for(int itoy=1; itoy<=NUM_TOYS; itoy++) {
	  //	if(itoy%10==0) cout<<TString::Format(" ---> producing toy %5d / %5d", itoy, NUM_TOYS)<<endl;

	  //	osc_test->Set_toy_variations( 1 );

	  //	vec_toydata_spectrum.clear();
	  //	int rows = osc_test->map_matrix_tosc_toy_pred[1].GetNcols();
	  //	for(int idx=0; idx<rows; idx++) vec_toydata_spectrum.push_back( osc_test->map_matrix_tosc_toy_pred[1](0, idx) );
	  //	i_toydata++;
	  //	tree_toydata->Fill();

	  // }
	  // cout<<endl;

	  ///////
	  tree_toydata->Write();
	  outfile_toydata->Close();


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
