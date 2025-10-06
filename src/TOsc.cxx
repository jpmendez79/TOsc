#include "WCPLEEANA/TOsc.h"

#include "draw.icc"

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Plot_user()
{
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

double TOsc::FCN(const double *par)
{
  double chi2_final = 0;
  double fit_dm2_41         = par[0];
  double fit_sin2_2theta_14 = par[1];
  double fit_sin2_theta_24  = par[2];
        
  /////// standard order
  Set_oscillation_pars(fit_dm2_41, fit_sin2_2theta_14, fit_sin2_theta_24, 0);  
  Apply_oscillation();
  Set_apply_POT();// meas, CV, COV: all ready

  ///////
  
  TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld;
  TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred;        
  TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total;
  int rows = matrix_cov_syst_total.GetNrows();
  
  /////// modify

  /// BNB only
  // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 0, 26*7-1);
  // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 0, 26*7-1);
  // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(0, 26*7-1, 0, 26*7-1);
  // int rows = matrix_cov_syst_total.GetNrows();
 
  /// NuMI only
  // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 26*7, 26*14-1);
  // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 26*7, 26*14-1);
  // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(26*7, 26*14-1, 26*7, 26*14-1);
  // int rows = matrix_cov_syst_total.GetNrows();  

  ///////  
  
  TMatrixD matrix_cov_stat_total(rows, rows);
  TMatrixD matrix_cov_total(rows, rows);

  for(int idx=0; idx<rows; idx++) {
    double val_stat_cov = 0;        
    // double val_data = matrix_data_total(0, idx);
    double val_pred = matrix_pred_total(0, idx);
    
    // if( val_data==0 ) { val_stat_cov = val_pred/2; }
    // else {
    //   if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
    //   else val_stat_cov = val_data;
    // }

    val_stat_cov = val_pred;// Pearson's format

    if( val_stat_cov==0 ) val_stat_cov = 1e-6;
    
    matrix_cov_stat_total(idx, idx) = val_stat_cov;
  }// for(int idx=0; idx<rows; idx++)
  
  matrix_cov_total = matrix_cov_syst_total + matrix_cov_stat_total;
  TMatrixD matrix_delta = matrix_pred_total - matrix_data_total;

  /////// default method: invert takes ~n^3 calcualtion
  // TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();   
  // TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();    
  // TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
  // chi2_final = matrix_chi2(0,0);            
  
  /////// Cholesky method: decomptakes (~n^3)/6 calcualtion + (~n^2) extral cal
  chi2_final = Cal_chi2COV(matrix_delta, matrix_cov_total);


  return chi2_final;
}

///////

// double TOsc::FCN_presave_only_PRED()
// {
//   double chi2_final = 0;
//   // double fit_dm2_41         = par[0];
//   // double fit_sin2_2theta_14 = par[1];
//   // double fit_sin2_theta_24  = par[2];
        
//   /////// standard order
//   // Set_oscillation_pars(fit_dm2_41, fit_sin2_2theta_14, fit_sin2_theta_24, 0);  
//   // Apply_oscillation();
//   Set_apply_POT();// meas, CV, COV: all ready

//   ///////
  
//   TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld;
//   TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred;        
//   TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total;
//   int rows = matrix_cov_syst_total.GetNrows();
  
//   /////// modify

//   /// BNB only
//   // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 0, 26*7-1);
//   // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 0, 26*7-1);
//   // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(0, 26*7-1, 0, 26*7-1);
//   // int rows = matrix_cov_syst_total.GetNrows();
 
//   /// NuMI only
//   // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 26*7, 26*14-1);
//   // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 26*7, 26*14-1);
//   // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(26*7, 26*14-1, 26*7, 26*14-1);
//   // int rows = matrix_cov_syst_total.GetNrows();  
  
//   ///////  
  
//   TMatrixD matrix_cov_stat_total(rows, rows);
//   TMatrixD matrix_cov_total(rows, rows);

//   for(int idx=0; idx<rows; idx++) {
//     double val_stat_cov = 0;        
//     double val_data = matrix_data_total(0, idx);
//     double val_pred = matrix_pred_total(0, idx);
    
//     if( val_data==0 ) { val_stat_cov = val_pred/2; }
//     else {
//       if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
//       else val_stat_cov = val_data;
//     }

//     if( val_stat_cov==0 ) val_stat_cov = 1e-6;
    
//     matrix_cov_stat_total(idx, idx) = val_stat_cov;
//   }// for(int idx=0; idx<rows; idx++)
  
//   matrix_cov_total = matrix_cov_syst_total + matrix_cov_stat_total;
//   TMatrixD matrix_delta = matrix_pred_total - matrix_data_total;

//   /////// default method: invert takes ~n^3 calcualtion  
//   // TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();   
//   // TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();
//   // TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
//   // chi2_final = matrix_chi2(0,0);      
  
//   /////// Cholesky method: decomptakes (~n^3)/6 calcualtion + (~n^2) extral cal
//   chi2_final = Cal_chi2COV(matrix_delta, matrix_cov_total);

//   ///////

//   return chi2_final;
// }

///////

// double TOsc::FCN_presave_both_PRED_COVsyst()
// {
//   double chi2_final = 0;
//   // double fit_dm2_41         = par[0];
//   // double fit_sin2_2theta_14 = par[1];
//   // double fit_sin2_theta_24  = par[2];
        
//   /////// standard order
//   // Set_oscillation_pars(fit_dm2_41, fit_sin2_2theta_14, fit_sin2_theta_24, 0);  
//   // Apply_oscillation();
//   // Set_apply_POT();// meas, CV, COV: all ready

//   ///////
  
//   TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld;
//   TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred;        
//   TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total;
//   int rows = matrix_cov_syst_total.GetNrows();
  
//   /////// modify

//   /// BNB only
//   // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 0, 26*7-1);
//   // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 0, 26*7-1);
//   // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(0, 26*7-1, 0, 26*7-1);
//   // int rows = matrix_cov_syst_total.GetNrows();
 
//   /// NuMI only
//   // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 26*7, 26*14-1);
//   // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 26*7, 26*14-1);
//   // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(26*7, 26*14-1, 26*7, 26*14-1);
//   // int rows = matrix_cov_syst_total.GetNrows();  
  
//   ///////  
  
//   TMatrixD matrix_cov_stat_total(rows, rows);
//   TMatrixD matrix_cov_total(rows, rows);

//   for(int idx=0; idx<rows; idx++) {
//     double val_stat_cov = 0;        
//     double val_data = matrix_data_total(0, idx);
//     double val_pred = matrix_pred_total(0, idx);
    
//     if( val_data==0 ) { val_stat_cov = val_pred/2; }
//     else {
//       if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
//       else val_stat_cov = val_data;
//     }

//     if( val_stat_cov==0 ) val_stat_cov = 1e-6;
    
//     matrix_cov_stat_total(idx, idx) = val_stat_cov;
//   }// for(int idx=0; idx<rows; idx++)
  
//   matrix_cov_total = matrix_cov_syst_total + matrix_cov_stat_total;
//   TMatrixD matrix_delta = matrix_pred_total - matrix_data_total;

//   /////// default method: invert takes ~n^3 calcualtion  
//   // TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();   
//   // TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();
//   // TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
//   // chi2_final = matrix_chi2(0,0);      
  
//   /////// Cholesky method: decomptakes (~n^3)/6 calcualtion + (~n^2) extral cal
//   chi2_final = Cal_chi2COV(matrix_delta, matrix_cov_total);

//   ///////

//   return chi2_final;
// }

///////

// double TOsc::FCN_Pearson_FCnew(const double *par)
// {
//   double chi2_final = 0;
//   double fit_dm2_41         = par[0];
//   double fit_sin2_2theta_14 = par[1];
//   double fit_sin2_theta_24  = par[2];
        
//   /////// standard order
//   Set_oscillation_pars(fit_dm2_41, fit_sin2_2theta_14, fit_sin2_theta_24, 0);  
//   Apply_oscillation();
//   Set_apply_POT();// meas, CV, COV: all ready

//   ///////
  
//   TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld;
//   TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred;        
//   TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total;
//   int rows = matrix_cov_syst_total.GetNrows();
  
//   /////// modify

//   /// BNB only
//   // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 0, 26*7-1);
//   // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 0, 26*7-1);
//   // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(0, 26*7-1, 0, 26*7-1);
//   // int rows = matrix_cov_syst_total.GetNrows();
 
//   /// NuMI only
//   // TMatrixD matrix_data_total = matrix_tosc_fitdata_newworld.GetSub(0,0, 26*7, 26*14-1);
//   // TMatrixD matrix_pred_total = matrix_tosc_eff_newworld_pred.GetSub(0,0, 26*7, 26*14-1);
//   // TMatrixD matrix_cov_syst_total = matrix_tosc_eff_newworld_abs_syst_total.GetSub(26*7, 26*14-1, 26*7, 26*14-1);
//   // int rows = matrix_cov_syst_total.GetNrows();  

//   ///////  
  
//   TMatrixD matrix_cov_stat_total(rows, rows);
//   TMatrixD matrix_cov_total(rows, rows);

//   for(int idx=0; idx<rows; idx++) {
//     double val_stat_cov = 0;        
//     // double val_data = matrix_data_total(0, idx);
//     double val_pred = matrix_pred_total(0, idx);
    
//     // if( val_data==0 ) { val_stat_cov = val_pred/2; }
//     // else {
//     //   if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
//     //   else val_stat_cov = val_data;
//     // }

//     val_stat_cov = val_pred;// Pearson's format;

//     if( val_stat_cov==0 ) val_stat_cov = 1e-6;    
//     matrix_cov_stat_total(idx, idx) = val_stat_cov;

//   }// for(int idx=0; idx<rows; idx++)
  
//   matrix_cov_total = matrix_cov_syst_total + matrix_cov_stat_total;
//   TMatrixD matrix_delta = matrix_pred_total - matrix_data_total;

//   /////// default method: invert takes ~n^3 calcualtion
//   TMatrixD matrix_delta_T = matrix_delta.T(); matrix_delta.T();   
//   TMatrixD matrix_cov_total_inv = matrix_cov_total; matrix_cov_total_inv.Invert();    
//   TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
//   chi2_final = matrix_chi2(0,0);            
  
//   /// for FC new procedure
//   matrix_tosc_chi2_COV_both_syst_stat_inv.ResizeTo(rows, rows);
//   matrix_tosc_chi2_COV_both_syst_stat_inv = matrix_cov_total_inv;
//   matrix_tosc_chi2_pred.ResizeTo(1, rows);
//   matrix_tosc_chi2_pred = matrix_pred_total;

//   /////// Cholesky method: decomptakes (~n^3)/6 calcualtion + (~n^2) extral cal
//   // chi2_final = Cal_chi2COV(matrix_delta, matrix_cov_total);


//   return chi2_final;
// }

///////

void TOsc::Minimization_OscPars_FullCov(double init_dm2_41, double init_sin2_2theta_14, double init_sin2_theta_24, double init_sin2_theta_34, TString roostr_flag_fixpar)
{
  ROOT::Minuit2::Minuit2Minimizer min_osc( ROOT::Minuit2::kMigrad );
  min_osc.SetPrintLevel(2);
  min_osc.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
  min_osc.SetMaxFunctionCalls(50000);
  min_osc.SetMaxIterations(50000);
  min_osc.SetTolerance(1e-5); // tolerance*2e-3 = edm precision
  min_osc.SetPrecision(1e-18); //precision in the target function
    
  /// set fitting parameters
  ROOT::Math::Functor Chi2Functor_osc(
				      [&](const double *par) {return FCN( par );},// FCN
				      3// number of fitting parameters
				      );
    
  min_osc.SetFunction(Chi2Functor_osc);
    
  min_osc.SetVariable( 0, "dm2_41", init_dm2_41, 1e-2);
  min_osc.SetVariable( 1, "sin2_theta_14", init_sin2_2theta_14, 1e-4);
  min_osc.SetVariable( 2, "sin2_theta_24", init_sin2_theta_24, 1e-4);

  min_osc.SetLowerLimitedVariable(0, "dm2_41", init_dm2_41, 1e-3, 0);  
  min_osc.SetLimitedVariable(1, "sin2_theta_14", init_sin2_2theta_14, 1e-4, 0, 1);
  min_osc.SetLimitedVariable(2, "sin2_theta_24", init_sin2_theta_24, 1e-4, 0, 1);
    
  if( roostr_flag_fixpar.Contains("dm2") ) {
    min_osc.SetFixedVariable( 0, "dm2_41", init_dm2_41 );
  }
  if( roostr_flag_fixpar.Contains("t14") ) {
    min_osc.SetFixedVariable( 1, "sin2_theta_14", init_sin2_2theta_14 );
  }
  if( roostr_flag_fixpar.Contains("t24") ) {
    min_osc.SetFixedVariable( 2, "sin2_theta_24", init_sin2_theta_24 );
  }
    
  min_osc.Minimize();

  ///////
    
  minimization_status = min_osc.Status();
  minimization_chi2   = min_osc.MinValue();
    
  const double *par_val = min_osc.X();
  const double *par_err = min_osc.Errors();

  minimization_dm2_41_val = par_val[0];
  minimization_dm2_41_err = par_err[0];
    
  minimization_sin2_2theta_14_val = par_val[1];
  minimization_sin2_2theta_14_err = par_err[1];
  
  minimization_sin2_theta_24_val = par_val[2];
  minimization_sin2_theta_24_err = par_err[2];

  if( par_val[0]!=par_val[0] ) minimization_status = 123;
  if( par_val[1]!=par_val[1] ) minimization_status = 123;
  if( par_val[2]!=par_val[2] ) minimization_status = 123;
  if( par_err[0]!=par_err[0] ) minimization_status = 124;
  if( par_err[1]!=par_err[1] ) minimization_status = 124;
  if( par_err[2]!=par_err[2] ) minimization_status = 124;

  if( 1 ) {
    cout<<endl;
    // cout<<TString::Format(" ---> minimization, status %2d, chi2 %6.4f, dm2 %4.2f +/- %4.2f, s22t14 %5.3f +/- %5.3f",
    // 			  minimization_status, minimization_chi2,
    // 			  minimization_dm2_41_val, minimization_dm2_41_err,
    // 			  minimization_sin2_2theta_14_val, minimization_sin2_2theta_14_err
    // 			  )<<endl;
    cout<<TString::Format(" ---> best-fit, status %2d, chi2 %6.4f", minimization_status, minimization_chi2)<<endl;
  }
   
}

//////////////////////////////////////////////////////////////////////////////////////////////////// cccg

// Y constrained by X

int TOsc::Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst, int index)
{
  int rows = matrix_syst.GetNrows();
  if( num_Y+num_X!=rows ) {cout<<" ERROR: num_Y+num_X!=rows"<<endl; exit(10);}

  ///////
  
  cout<<Form(" ---> Goodness of fit: %2d", index)<<endl;

  TString roostr = "";
  int color_no = kRed;
  int color_wi = kBlue;
  int color_dd = kBlack;

  /////////////////////////////////////////////////////////////////////// noConstraint

  matrix_pred.T(); matrix_data.T();
  TMatrixD matrix_pred_Y = matrix_pred.GetSub(0, num_Y-1, 0, 0);// Yx1
  TMatrixD matrix_data_Y = matrix_data.GetSub(0, num_Y-1, 0, 0);// Yx1
  matrix_pred.T(); matrix_data.T();

  TMatrixD matrix_cov_syst_both = matrix_syst;
  TMatrixD matrix_YY = matrix_cov_syst_both.GetSub(0, num_Y-1, 0, num_Y-1);

  ///////////////////////////// goodness of fit, Pearson's format test

  /// docDB 32520, when the prediction is sufficiently low
  double array_pred_protect[11] = {0, 0.461, 0.916, 1.382, 1.833, 2.298, 2.767, 3.225, 3.669, 4.141, 4.599};

  TMatrixD matrix_goodness_cov_stat(num_Y, num_Y);

  for( int i=0; i<num_Y; i++ ) {
    double val_pred = matrix_pred_Y(i, 0);
    double val_data = matrix_data_Y(i, 0);        
    matrix_goodness_cov_stat(i,i) = val_pred;
    
    int int_data = (int)(val_data+0.1);
    if( int_data>=1 && int_data<=10 ) {
      if( val_pred<array_pred_protect[int_data] ) {
    	double numerator = pow(val_pred-val_data, 2);
        double denominator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
    	matrix_goodness_cov_stat(i,i) = numerator/denominator;
    	cout<<" --------> GoF AA: Protection Protection"<<endl;
      }
    }  

    if( flag_goodness_of_fit_CNP ) {
      double val_stat_cov = 0;
      if( val_data==0 ) { val_stat_cov = val_pred/2; }
      else {
	if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
	else val_stat_cov = val_data;
      }      
      if( val_stat_cov==0 ) val_stat_cov = 1e-6;
      matrix_goodness_cov_stat(i,i) = val_stat_cov;
    } // if( flag_goodness_of_fit_CNP )

  }// for( int i=0; i<num_Y; i++ )

  TMatrixD matrix_goodness_cov_total_noConstraint = matrix_goodness_cov_stat + matrix_YY;
  TMatrixD matrix_delta_noConstraint_T = matrix_pred_Y - matrix_data_Y;// Yx1
  TMatrixD matrix_delta_noConstraint = matrix_delta_noConstraint_T.T(); matrix_delta_noConstraint_T.T();
  TMatrixD matrix_cov_noConstraint_inv = matrix_goodness_cov_total_noConstraint; matrix_cov_noConstraint_inv.Invert();
  double val_chi2_noConstraint = (matrix_delta_noConstraint*matrix_cov_noConstraint_inv*matrix_delta_noConstraint_T)(0,0);
  double p_value_noConstraint = TMath::Prob(val_chi2_noConstraint, num_Y);

  TVectorD vector_pred_Y_noConstraint(num_Y, matrix_pred_Y.GetMatrixArray());
  TVectorD vector_data_Y_noConstraint(num_Y, matrix_data_Y.GetMatrixArray());
  double val_pred_Y_noConstraint = vector_pred_Y_noConstraint.Sum();
  double val_data_Y_noConstraint = vector_data_Y_noConstraint.Sum();
  double val_d2p_ratio_noConstraint = val_data_Y_noConstraint/val_pred_Y_noConstraint;

  cout<<TString::Format(" ---> GoF noConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, pvalue %5.3f, data %9.2f, pred %9.2f, d/p %7.3f", 
			val_chi2_noConstraint, num_Y, val_chi2_noConstraint/num_Y, p_value_noConstraint,
			val_data_Y_noConstraint, val_pred_Y_noConstraint, val_d2p_ratio_noConstraint
			)<<endl;

  //////////////////////////////////////// Plotting

  int line_width = 2;

  TF1 *f1_line1 = new TF1("f1_line1", "1", 0, 1e6);
  f1_line1->SetLineColor(color_dd);
  f1_line1->SetLineWidth(1);
  f1_line1->SetLineStyle(7);

  ///////

  roostr = TString::Format("h1d_data_Y_noConstraint_val_%d", index);
  TH1D *h1d_data_Y_noConstraint_val = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_data_Y_noConstraint_err_%d", index);
  TH1D *h1d_data_Y_noConstraint_err = new TH1D(roostr, "", num_Y, 0, num_Y);

  roostr = TString::Format("h1d_pred_Y_noConstraint_val_%d", index);
  TH1D *h1d_pred_Y_noConstraint_val = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_pred_Y_noConstraint_err_%d", index);
  TH1D *h1d_pred_Y_noConstraint_err = new TH1D(roostr, "", num_Y, 0, num_Y);

  roostr = TString::Format("h1d_d2p_Y_noConstraint_val_%d", index);
  TH1D *h1d_d2p_Y_noConstraint_val = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_d2p_Y_noConstraint_stat_err_%d", index);
  TH1D *h1d_d2p_Y_noConstraint_stat_err = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_d2p_Y_noConstraint_syst_err_%d", index);
  TH1D *h1d_d2p_Y_noConstraint_syst_err = new TH1D(roostr, "", num_Y, 0, num_Y);

  for(int idx=1; idx<=num_Y; idx++) {
    double data_val = matrix_data_Y(idx-1, 0);
    double data_err = sqrt(data_val);// statistics

    double pred_val = matrix_pred_Y(idx-1, 0);
    double pred_err = sqrt(matrix_YY(idx-1, idx-1));// systematics

    double ratio_val      = 0;
    double ratio_stat_err = 0;
    double ratio_syst_err = 0;
    if( pred_val!=0 ) {
      ratio_val = data_val/pred_val;
      ratio_stat_err = data_err/pred_val;
      ratio_syst_err = pred_err/pred_val;
    }

    // cout<<" ---> check: "<<idx<<"\t"<<pred_val<<"\t"<<pred_err/pred_val<<endl;

    h1d_data_Y_noConstraint_val->SetBinContent(idx, data_val);
    h1d_data_Y_noConstraint_err->SetBinContent(idx, data_val);
    h1d_data_Y_noConstraint_err->SetBinError(idx, data_err);

    h1d_pred_Y_noConstraint_val->SetBinContent(idx, pred_val);
    h1d_pred_Y_noConstraint_err->SetBinContent(idx, pred_val);
    h1d_pred_Y_noConstraint_err->SetBinError(idx, pred_err);

    h1d_d2p_Y_noConstraint_val->SetBinContent(idx, ratio_val);
    h1d_d2p_Y_noConstraint_stat_err->SetBinContent(idx, ratio_val);
    h1d_d2p_Y_noConstraint_stat_err->SetBinError(idx, ratio_stat_err);
    h1d_d2p_Y_noConstraint_syst_err->SetBinContent(idx, 1);
    h1d_d2p_Y_noConstraint_syst_err->SetBinError(idx, ratio_syst_err);
  }// for(int idx=1; idx<=num_Y; idx++)

  double ymax_data = h1d_data_Y_noConstraint_val->GetMaximum();
  double ymax_pred = h1d_pred_Y_noConstraint_val->GetMaximum();
  double ymax_eff  = 1e6;
  if(ymax_data>ymax_pred*1.3) ymax_eff = ymax_data*1.2;
  else ymax_eff = ymax_pred * 1.3;

  ///////
  roostr = TString::Format("canv_spectra_GoF_no_%02d", index);
  TCanvas *canv_spectra_GoF_no = new TCanvas(roostr, roostr, 1000, 950);  

  ///////
  canv_spectra_GoF_no->cd();
  TPad *pad_top_no = new TPad("pad_top_no", "pad_top_no", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_no, 0.15, 0.1, 0.1, 0.06);
  pad_top_no->Draw(); pad_top_no->cd();

  h1d_pred_Y_noConstraint_err->Draw("e2"); h1d_pred_Y_noConstraint_err->SetMinimum(0);
  h1d_pred_Y_noConstraint_err->SetMaximum(ymax_eff);
  h1d_pred_Y_noConstraint_err->SetMarkerSize(0.);
  h1d_pred_Y_noConstraint_err->SetFillColor(color_no); h1d_pred_Y_noConstraint_err->SetFillStyle(3004);
  h1d_pred_Y_noConstraint_err->SetLineColor(color_no);
  func_title_size(h1d_pred_Y_noConstraint_err, 0.06, 0.06, 0.06, 0.06);
  func_xy_title(h1d_pred_Y_noConstraint_err, "", "Events", 1);
  h1d_pred_Y_noConstraint_err->GetYaxis()->SetTitleOffset(1.2);

  h1d_pred_Y_noConstraint_val->Draw("same hist");
  h1d_pred_Y_noConstraint_val->SetLineColor(color_no);
  h1d_pred_Y_noConstraint_val->SetLineWidth(line_width);

  h1d_data_Y_noConstraint_err->Draw("same e1");
  h1d_data_Y_noConstraint_err->SetLineColor(color_dd);
  h1d_data_Y_noConstraint_err->SetLineWidth(line_width);
  h1d_data_Y_noConstraint_err->SetMarkerSize(1.2);

  h1d_pred_Y_noConstraint_err->Draw("same axis");  

  TLegend *lg_top_no = new TLegend(0.5+0.05, 0.85-0.07*4, 0.85, 0.85);
  lg_top_no->AddEntry(h1d_data_Y_noConstraint_err, TString::Format("Data: %1.0f", val_data_Y_noConstraint), "lep");
  lg_top_no->AddEntry(h1d_pred_Y_noConstraint_err, TString::Format("#color[%d]{Pred no constraint}", color_no), "lf");
  lg_top_no->AddEntry("", TString::Format("#color[%d]{Data/Pred: %4.2f}", color_no, val_d2p_ratio_noConstraint), "");
  lg_top_no->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_no, val_chi2_noConstraint, num_Y), "");
  lg_top_no->Draw();
  lg_top_no->SetBorderSize(0); lg_top_no->SetFillStyle(0); lg_top_no->SetTextSize(0.06);

  ///////
  canv_spectra_GoF_no->cd();
  TPad *pad_bot_no = new TPad("pad_bot_no", "pad_bot_no", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_no, 0.15, 0.1, 0.05, 0.3);
  pad_bot_no->Draw(); pad_bot_no->cd();

  h1d_d2p_Y_noConstraint_syst_err->Draw("e2");
  h1d_d2p_Y_noConstraint_syst_err->SetMinimum(0); h1d_d2p_Y_noConstraint_syst_err->SetMaximum(2);
  h1d_d2p_Y_noConstraint_syst_err->SetMarkerSize(0.);
  h1d_d2p_Y_noConstraint_syst_err->SetFillColor(color_no); h1d_d2p_Y_noConstraint_syst_err->SetFillStyle(3004);
  h1d_d2p_Y_noConstraint_syst_err->SetLineColor(color_no);
  func_title_size(h1d_d2p_Y_noConstraint_syst_err, 0.074, 0.074, 0.074, 0.074);
  func_xy_title(h1d_d2p_Y_noConstraint_syst_err, "Bin Index", "Data / Pred", 1);
  h1d_d2p_Y_noConstraint_syst_err->GetYaxis()->SetNdivisions(509);
  h1d_d2p_Y_noConstraint_syst_err->GetYaxis()->SetTitleOffset(0.99);

  f1_line1->Draw("same");
  
  h1d_d2p_Y_noConstraint_stat_err->Draw("same e1");
  h1d_d2p_Y_noConstraint_stat_err->SetLineColor(color_no);
  h1d_d2p_Y_noConstraint_stat_err->SetLineWidth(line_width);
  h1d_d2p_Y_noConstraint_stat_err->SetMarkerSize(1.2);
  h1d_d2p_Y_noConstraint_stat_err->SetMarkerColor(color_no);

  h1d_d2p_Y_noConstraint_syst_err->Draw("same axis");

  ///////
    
  roostr = TString::Format("canv_spectra_GoF_%06d_no.png", index); canv_spectra_GoF_no->SaveAs(roostr);
  cout<<roostr<<endl;

  ///////
  TMatrixD matrix_pred_Y_T = matrix_pred_Y.T(); matrix_pred_Y.T();
  TMatrixD matrix_data_Y_T = matrix_data_Y.T(); matrix_data_Y.T();   
  Exe_GoF_decompose(matrix_pred_Y_T, matrix_data_Y_T, matrix_goodness_cov_total_noConstraint, index, "_no");

  if( num_Y==rows ) return 1;

  /////////////////////////////////////////////////////////////////////// wiConstraint

  matrix_pred.T(); matrix_data.T();
  TMatrixD matrix_pred_X = matrix_pred.GetSub(num_Y, num_Y+num_X-1, 0, 0);// Xx1
  TMatrixD matrix_data_X = matrix_data.GetSub(num_Y, num_Y+num_X-1, 0, 0);// Xx1
  matrix_pred.T(); matrix_data.T();

  TMatrixD matrix_XX = matrix_cov_syst_both.GetSub(num_Y, num_Y+num_X-1, num_Y, num_Y+num_X-1);

  for(int ibin=1; ibin<=num_X; ibin++) {    
    double val_data = matrix_data_X(ibin-1,0);
    double val_pred = matrix_pred_X(ibin-1,0);
    double cov_stat = val_pred;

    int int_meas = (int)(val_data+0.1);    
    if( int_meas>=1 && int_meas<=10) {
      if( val_pred<array_pred_protect[int_meas] ) {
	double numerator = pow(val_pred-val_data, 2);
        double denominator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
	cov_stat = numerator/denominator;
	cout<<" --------> GoF AB: Protection Protection"<<endl;
      }
    }    

    if( flag_goodness_of_fit_CNP ) {
      double val_stat_cov = 0;
      if( val_data==0 ) { val_stat_cov = val_pred/2; }
      else {
	if( val_pred!=0 ) val_stat_cov = 3./( 1./val_data + 2./val_pred );
	else val_stat_cov = val_data;
      }      
      if( val_stat_cov==0 ) val_stat_cov = 1e-6;
      cov_stat = val_stat_cov;
    } // if( flag_goodness_of_fit_CNP )

    matrix_XX(ibin-1, ibin-1) += cov_stat;        
  }// for(int ibin=1; ibin<=num_X; ibin++)

  TMatrixD matrix_XX_inv = matrix_XX; matrix_XX_inv.Invert();

  ////////

  TMatrixD matrix_YX = matrix_cov_syst_both.GetSub(0, num_Y-1, num_Y, num_Y+num_X-1);
  TMatrixD matrix_XY = matrix_YX.T(); matrix_YX.T();
  TMatrixD matrix_Y_under_X = matrix_pred_Y + matrix_YX * matrix_XX_inv * (matrix_data_X - matrix_pred_X);
  TMatrixD matrix_YY_under_XX = matrix_YY - matrix_YX * matrix_XX_inv * matrix_XY;
  // Here, only for systetmaics uncertainty because of no stat in matrix_YY

  ////////

  TMatrixD matrix_goodness_cov_stat_wiConstraint(num_Y, num_Y);
  for( int i=0; i<num_Y; i++ ) {
    double val_pred = matrix_Y_under_X(i, 0);
    double val_data = matrix_data_Y(i, 0);    
    matrix_goodness_cov_stat_wiConstraint(i,i) = val_pred;
    
    int int_data = (int)(val_data+0.1);
    if( int_data>=1 && int_data<=10) {
      if( val_pred<array_pred_protect[int_data] ) {
	double numerator = pow(val_pred-val_data, 2);
        double denominator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
        matrix_goodness_cov_stat_wiConstraint(i,i) = numerator/denominator;
	cout<<" --------> GoF BB: Protection Protection"<<endl;
      }
    }    
  }// for( int i=0; i<num_Y; i++ )

  TMatrixD matrix_goodness_cov_total_wiConstraint = matrix_goodness_cov_stat_wiConstraint + matrix_YY_under_XX;

  ////////

  TMatrixD matrix_delta_wiConstraint_T = matrix_Y_under_X - matrix_data_Y;// Yx1
  TMatrixD matrix_delta_wiConstraint = matrix_delta_wiConstraint_T.T(); matrix_delta_wiConstraint_T.T(); 
  TMatrixD matrix_cov_wiConstraint_inv = matrix_goodness_cov_total_wiConstraint; matrix_cov_wiConstraint_inv.Invert();
  double val_chi2_wiConstraint = (matrix_delta_wiConstraint*matrix_cov_wiConstraint_inv*matrix_delta_wiConstraint_T)(0,0);
  double p_value_wiConstraint = TMath::Prob( val_chi2_wiConstraint, num_Y );

  TVectorD vector_pred_Y_wiConstraint(num_Y, matrix_Y_under_X.GetMatrixArray());
  TVectorD vector_data_Y_wiConstraint(num_Y, matrix_data_Y.GetMatrixArray());
  double val_pred_Y_wiConstraint = vector_pred_Y_wiConstraint.Sum();
  double val_data_Y_wiConstraint = vector_data_Y_wiConstraint.Sum();
  double val_d2p_ratio_wiConstraint = val_data_Y_wiConstraint/val_pred_Y_wiConstraint;

  cout<<TString::Format(" ---> GoF wiConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, pvalue %5.3f, data %9.2f, pred %9.2f, d/p %7.3f", 
			val_chi2_wiConstraint, num_Y, val_chi2_wiConstraint/num_Y, p_value_wiConstraint,
			val_data_Y_wiConstraint, val_pred_Y_wiConstraint, val_d2p_ratio_wiConstraint
			)<<endl;

  //////////////////////////////////////// Plotting

  roostr = TString::Format("h1d_data_Y_wiConstraint_val_%d", index);
  TH1D *h1d_data_Y_wiConstraint_val = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_data_Y_wiConstraint_err_%d", index);
  TH1D *h1d_data_Y_wiConstraint_err = new TH1D(roostr, "", num_Y, 0, num_Y);

  roostr = TString::Format("h1d_pred_Y_wiConstraint_val_%d", index);
  TH1D *h1d_pred_Y_wiConstraint_val = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_pred_Y_wiConstraint_err_%d", index);
  TH1D *h1d_pred_Y_wiConstraint_err = new TH1D(roostr, "", num_Y, 0, num_Y);

  roostr = TString::Format("h1d_d2p_Y_wiConstraint_val_%d", index);
  TH1D *h1d_d2p_Y_wiConstraint_val = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_d2p_Y_wiConstraint_stat_err_%d", index);
  TH1D *h1d_d2p_Y_wiConstraint_stat_err = new TH1D(roostr, "", num_Y, 0, num_Y);
  roostr = TString::Format("h1d_d2p_Y_wiConstraint_syst_err_%d", index);
  TH1D *h1d_d2p_Y_wiConstraint_syst_err = new TH1D(roostr, "", num_Y, 0, num_Y);

  for(int idx=1; idx<=num_Y; idx++) {
    double data_val = matrix_data_Y(idx-1, 0);
    double data_err = sqrt(data_val);// statistics

    double pred_val = matrix_Y_under_X(idx-1, 0);
    double pred_err = sqrt(matrix_YY_under_XX(idx-1, idx-1));// systematics

    double ratio_val      = 0;
    double ratio_stat_err = 0;
    double ratio_syst_err = 0;
    if( pred_val!=0 ) {
      ratio_val = data_val/pred_val;
      ratio_stat_err = data_err/pred_val;
      ratio_syst_err = ratio_val * (pred_err/pred_val);
    }

    double pred_val_no = matrix_pred_Y(idx-1, 0);
    if( pred_val_no!=0 ) {
      ratio_syst_err = pred_err/pred_val_no;
    }

    // cout<<" ---> check: "<<idx<<"\t"<<pred_val<<"\t"<<pred_err/pred_val<<endl;

    h1d_data_Y_wiConstraint_val->SetBinContent(idx, data_val);
    h1d_data_Y_wiConstraint_err->SetBinContent(idx, data_val);
    h1d_data_Y_wiConstraint_err->SetBinError(idx, data_err);

    h1d_pred_Y_wiConstraint_val->SetBinContent(idx, pred_val);
    h1d_pred_Y_wiConstraint_err->SetBinContent(idx, pred_val);
    h1d_pred_Y_wiConstraint_err->SetBinError(idx, pred_err);

    h1d_d2p_Y_wiConstraint_val->SetBinContent(idx, ratio_val);
    h1d_d2p_Y_wiConstraint_stat_err->SetBinContent(idx, ratio_val);
    h1d_d2p_Y_wiConstraint_stat_err->SetBinError(idx, ratio_stat_err);
    h1d_d2p_Y_wiConstraint_syst_err->SetBinContent(idx, 1);
    h1d_d2p_Y_wiConstraint_syst_err->SetBinError(idx, ratio_syst_err);
  }// for(int idx=1; idx<=num_Y; idx++)

  ///////
  roostr = TString::Format("canv_spectra_GoF_wi_%02d", index);
  TCanvas *canv_spectra_GoF_wi = new TCanvas(roostr, roostr, 1000, 950);  

  ///////
  canv_spectra_GoF_wi->cd();
  TPad *pad_top_wi = new TPad("pad_top_wi", "pad_top_wi", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_wi, 0.15, 0.1, 0.1, 0.06);
  pad_top_wi->Draw(); pad_top_wi->cd();

  h1d_pred_Y_wiConstraint_err->Draw("e2"); h1d_pred_Y_wiConstraint_err->SetMinimum(0);
  h1d_pred_Y_wiConstraint_err->SetMaximum(ymax_eff);
  h1d_pred_Y_wiConstraint_err->SetMarkerSize(0.);
  h1d_pred_Y_wiConstraint_err->SetFillColor(color_wi); h1d_pred_Y_wiConstraint_err->SetFillStyle(3005);
  h1d_pred_Y_wiConstraint_err->SetLineColor(color_wi);
  func_title_size(h1d_pred_Y_wiConstraint_err, 0.06, 0.06, 0.06, 0.06);
  func_xy_title(h1d_pred_Y_wiConstraint_err, "", "Events", 1);
  h1d_pred_Y_wiConstraint_err->GetYaxis()->SetTitleOffset(1.2);

  h1d_pred_Y_wiConstraint_val->Draw("same hist");
  h1d_pred_Y_wiConstraint_val->SetLineColor(color_wi);
  h1d_pred_Y_wiConstraint_val->SetLineWidth(line_width);

  h1d_data_Y_wiConstraint_err->Draw("same e1");
  h1d_data_Y_wiConstraint_err->SetLineColor(color_dd);
  h1d_data_Y_wiConstraint_err->SetLineWidth(line_width);
  h1d_data_Y_wiConstraint_err->SetMarkerSize(1.2);

  h1d_pred_Y_wiConstraint_err->Draw("same axis");  

  TLegend *lg_top_wi = new TLegend(0.5+0.05, 0.85-0.07*4, 0.85, 0.85);
  lg_top_wi->AddEntry(h1d_data_Y_wiConstraint_err, TString::Format("Data: %1.0f", val_data_Y_noConstraint), "lep");
  lg_top_wi->AddEntry(h1d_pred_Y_wiConstraint_err, TString::Format("#color[%d]{Pred wi constraint}", color_wi), "lf");
  lg_top_wi->AddEntry("", TString::Format("#color[%d]{Data/Pred: %4.2f}", color_wi, val_d2p_ratio_wiConstraint), "");
  lg_top_wi->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_wi, val_chi2_wiConstraint, num_Y), "");
  lg_top_wi->Draw();
  lg_top_wi->SetBorderSize(0); lg_top_wi->SetFillStyle(0); lg_top_wi->SetTextSize(0.06);

  ///////
  canv_spectra_GoF_wi->cd();
  TPad *pad_bot_wi = new TPad("pad_bot_wi", "pad_bot_wi", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_wi, 0.15, 0.1, 0.05, 0.3);
  pad_bot_wi->Draw(); pad_bot_wi->cd();

  h1d_d2p_Y_wiConstraint_syst_err->Draw("e2");
  h1d_d2p_Y_wiConstraint_syst_err->SetMinimum(0); h1d_d2p_Y_wiConstraint_syst_err->SetMaximum(2);
  h1d_d2p_Y_wiConstraint_syst_err->SetMarkerSize(0.);
  h1d_d2p_Y_wiConstraint_syst_err->SetFillColor(color_wi); h1d_d2p_Y_wiConstraint_syst_err->SetFillStyle(3005);
  h1d_d2p_Y_wiConstraint_syst_err->SetLineColor(color_wi);
  func_title_size(h1d_d2p_Y_wiConstraint_syst_err, 0.074, 0.074, 0.074, 0.074);
  func_xy_title(h1d_d2p_Y_wiConstraint_syst_err, "Bin Index", "Data / Pred", 1);
  h1d_d2p_Y_wiConstraint_syst_err->GetYaxis()->SetNdivisions(509);
  h1d_d2p_Y_wiConstraint_syst_err->GetYaxis()->SetTitleOffset(0.99);

  f1_line1->Draw("same");
  
  h1d_d2p_Y_wiConstraint_stat_err->Draw("same e1");
  h1d_d2p_Y_wiConstraint_stat_err->SetLineColor(color_wi);
  h1d_d2p_Y_wiConstraint_stat_err->SetLineWidth(line_width);
  h1d_d2p_Y_wiConstraint_stat_err->SetMarkerSize(1.2);
  h1d_d2p_Y_wiConstraint_stat_err->SetMarkerColor(color_wi);

  h1d_d2p_Y_wiConstraint_syst_err->Draw("same axis");

  roostr = TString::Format("canv_spectra_GoF_%06d_wi.png", index); canv_spectra_GoF_wi->SaveAs(roostr);

  ///////
  TMatrixD matrix_Y_under_X_T = matrix_Y_under_X.T(); matrix_Y_under_X.T();
  Exe_GoF_decompose(matrix_Y_under_X_T, matrix_data_Y_T, matrix_goodness_cov_total_wiConstraint, index, "_wi");

  //////////////////////////////////////// Plotting both

  roostr = TString::Format("canv_spectra_GoF_both_%02d", index);
  TCanvas *canv_spectra_GoF_both = new TCanvas(roostr, roostr, 1000, 950);  

  ///////
  canv_spectra_GoF_both->cd();
  TPad *pad_top_both = new TPad("pad_top_both", "pad_top_both", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_both, 0.15, 0.1, 0.1, 0.06);
  pad_top_both->Draw(); pad_top_both->cd();

  h1d_pred_Y_noConstraint_err->Draw("e2"); 
  h1d_pred_Y_noConstraint_val->Draw("same hist");
  h1d_data_Y_noConstraint_err->Draw("same e1");
 
  h1d_pred_Y_wiConstraint_err->Draw("same e2"); 
  h1d_pred_Y_wiConstraint_val->Draw("same hist");
  h1d_data_Y_wiConstraint_err->Draw("same e1");

  h1d_pred_Y_noConstraint_err->Draw("same axis");  

  TLegend *lg_top_both = new TLegend(0.5+0.05, 0.85-0.07*7, 0.85, 0.85);
  lg_top_both->AddEntry(h1d_data_Y_noConstraint_err, TString::Format("Data: %1.0f", val_data_Y_noConstraint), "lep");
  lg_top_both->AddEntry(h1d_pred_Y_noConstraint_err, TString::Format("#color[%d]{Pred no constraint}", color_no), "lf");
  lg_top_both->AddEntry("", TString::Format("#color[%d]{Data/Pred: %4.2f}", color_no, val_d2p_ratio_noConstraint), "");
  lg_top_both->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_no, val_chi2_noConstraint, num_Y), "");
  lg_top_both->AddEntry(h1d_pred_Y_wiConstraint_err, TString::Format("#color[%d]{Pred wi constraint}", color_wi), "lf");
  lg_top_both->AddEntry("", TString::Format("#color[%d]{Data/Pred: %4.2f}", color_wi, val_d2p_ratio_wiConstraint), "");
  lg_top_both->AddEntry("", TString::Format("#color[%d]{#chi^{2}/ndf: %3.2f/%d}", color_wi, val_chi2_wiConstraint, num_Y), "");
  lg_top_both->Draw();
  lg_top_both->SetBorderSize(0); lg_top_both->SetFillStyle(0); lg_top_both->SetTextSize(0.06);

  ///////
  canv_spectra_GoF_both->cd();
  TPad *pad_bot_both = new TPad("pad_bot_both", "pad_bot_both", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_both, 0.15, 0.1, 0.05, 0.3);
  pad_bot_both->Draw(); pad_bot_both->cd();

  h1d_d2p_Y_noConstraint_syst_err->Draw("e2");  
  h1d_d2p_Y_wiConstraint_syst_err->Draw("same e2");

  f1_line1->Draw("same");  

  h1d_d2p_Y_noConstraint_stat_err->Draw("same e1");
  h1d_d2p_Y_wiConstraint_stat_err->Draw("same e1");

  h1d_d2p_Y_noConstraint_syst_err->Draw("same axis");

  roostr = TString::Format("canv_spectra_GoF_%06d_both.png", index); canv_spectra_GoF_both->SaveAs(roostr);

  return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

int TOsc::Exe_GoF_decompose(TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_stat_syst_total, int index, TString postfix)
{
  int rows = matrix_stat_syst_total.GetNrows();

  TMatrixDSym DSmatrix_cov(rows, matrix_stat_syst_total.GetMatrixArray());
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();// NxN
  TVectorD vector_eigenvalue = DSmatrix_eigen.GetEigenValues();
  TMatrixD matrix_eigenvector_T = matrix_eigenvector.T(); matrix_eigenvector.T();

  TMatrixD matrix_lambda_pred = matrix_pred * matrix_eigenvector;// 1xN
  TMatrixD matrix_lambda_data = matrix_data * matrix_eigenvector;// 1xN
  TMatrixD matrix_delta_lambda = matrix_lambda_data - matrix_lambda_pred;// 1xN
  TMatrixD matrix_delta_lambda_T = matrix_delta_lambda.T(); matrix_delta_lambda.T();

  TMatrixD matrix_cov_lambda(rows, rows);
  for(int idx=0; idx<rows; idx++) matrix_cov_lambda(idx, idx) = vector_eigenvalue(idx);
  TMatrixD matrix_cov_lambda_inv = matrix_cov_lambda; matrix_cov_lambda_inv.Invert();
  double chi2_lambda = (matrix_delta_lambda * matrix_cov_lambda_inv * matrix_delta_lambda_T)(0,0);

  cout<<TString::Format(" ---> GoF decompose:    chi2 %6.2f, ndf %3d, chi2/ndf %6.2f", chi2_lambda, rows, chi2_lambda/rows)<<endl;

  //////////////////////////////////////// Plotting

  TString roostr = "";  

  roostr = TString::Format("h1d_epsilon_%d", index);     TH1D *h1d_epsilon = new TH1D(roostr+postfix, "", rows, 0, rows);  
  roostr = TString::Format("h1d_lambda_frac_%d", index); TH1D *h1d_lambda_frac = new TH1D(roostr+postfix, "", rows, 0, rows);
  
  double check_chi2 = 0;
  for(int idx=1; idx<=rows; idx++) {
    double val_eigensum   = vector_eigenvalue.Sum();
    double val_eigenvalue = vector_eigenvalue(idx-1);
    double val_epsilon = -matrix_delta_lambda(0, idx-1)/sqrt(val_eigenvalue);
    check_chi2 += val_epsilon*val_epsilon;

    h1d_epsilon->SetBinContent(idx, val_epsilon);
    h1d_lambda_frac->SetBinContent(idx, val_eigenvalue*100./val_eigensum);
  }
  cout<<TString::Format(" ---> decomp check_chi2:     %6.2f", check_chi2)<<endl;

  ///////
  roostr = TString::Format("canv_decomp_GoF_%02d", index) + postfix;
  TCanvas *canv_decomp_GoF_ = new TCanvas(roostr, roostr, 1000, 950);  

  ///////
  canv_decomp_GoF_->cd();
  TPad *pad_top_wi = new TPad("pad_top_wi", "pad_top_wi", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_wi, 0.15, 0.1, 0.1, 0.06);
  pad_top_wi->Draw(); pad_top_wi->cd();

  h1d_epsilon->Draw("p"); 
  h1d_epsilon->SetMinimum(-6); h1d_epsilon->SetMaximum(6);  
  h1d_epsilon->SetMarkerStyle(21); h1d_epsilon->SetMarkerColor(kBlack); h1d_epsilon->SetMarkerSize(1.2);  
  func_title_size(h1d_epsilon, 0.06, 0.06, 0.06, 0.06);
  func_xy_title(h1d_epsilon, "", "#epsilon value", 1);
  h1d_epsilon->GetYaxis()->SetTitleOffset(1.2);

  TF1 *f1_line0 = new TF1("f1_line0", "0", 0, 1e6);
  f1_line0->SetLineColor(kBlack);
  f1_line0->SetLineWidth(1);
  f1_line0->SetLineStyle(7);
  f1_line0->Draw("same");
     
  TH1D *h1_lambda_sigmax3 = new TH1D(TString::Format("h1_lambda_sigmax3_%d_%s", index, postfix.Data()), "", rows, 0, rows);
  TH1D *h1_lambda_sigmax2 = new TH1D(TString::Format("h1_lambda_sigmax2_%d_%s", index, postfix.Data()), "", rows, 0, rows);
  TH1D *h1_lambda_sigmax1 = new TH1D(TString::Format("h1_lambda_sigmax1_%d_%s", index, postfix.Data()), "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_lambda_sigmax3->SetBinError(ibin, 3);
    h1_lambda_sigmax2->SetBinError(ibin, 2);
    h1_lambda_sigmax1->SetBinError(ibin, 1);
  }

  h1_lambda_sigmax3->SetFillColor(kRed-9); h1_lambda_sigmax3->SetFillStyle(1001); h1_lambda_sigmax3->SetMarkerSize(0);
  h1_lambda_sigmax2->SetFillColor(kYellow); h1_lambda_sigmax2->SetFillStyle(1001); h1_lambda_sigmax2->SetMarkerSize(0);
  h1_lambda_sigmax1->SetFillColor(kGreen); h1_lambda_sigmax1->SetFillStyle(1001); h1_lambda_sigmax1->SetMarkerSize(0);

  h1_lambda_sigmax3->Draw("same e2");
  h1_lambda_sigmax2->Draw("same e2");
  h1_lambda_sigmax1->Draw("same e2");

  f1_line0->Draw("same");
  h1d_epsilon->Draw("same p"); 
  h1d_epsilon->Draw("same axis"); 

  TLegend *lg_top_ = new TLegend(0.5+0.05, 0.85-0.07*1, 0.85, 0.85);
  roostr = TString::Format("#chi^{2} = #sum #epsilon^{2}_{i} = %3.2f", check_chi2);
  lg_top_->AddEntry("", roostr, "");
  lg_top_->SetBorderSize(0); lg_top_->SetFillStyle(0); lg_top_->SetTextSize(0.06);
  lg_top_->Draw();

  ///////
  canv_decomp_GoF_->cd();
  TPad *pad_bot_wi = new TPad("pad_bot_wi", "pad_bot_wi", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_wi, 0.15, 0.1, 0.05, 0.3);
  pad_bot_wi->Draw(); pad_bot_wi->cd();
  pad_bot_wi->SetLogy();

  h1d_lambda_frac->Draw("hist");
  h1d_lambda_frac->SetMaximum(100);
  h1d_lambda_frac->SetLineColor(kBlack);
  h1d_lambda_frac->SetLineWidth(2);
  func_title_size(h1d_lambda_frac, 0.074, 0.074, 0.074, 0.074);
  func_xy_title(h1d_lambda_frac, "Decomposed Bin Index", "Error Fraction (%)", 1);
  //h1d_lambda_frac->GetYaxis()->SetNdivisions(509);
  h1d_lambda_frac->GetYaxis()->SetTitleOffset(0.99);

  ///////
  roostr = TString::Format("canv_decomp_GoF_%02d", index) + postfix + ".png";
  canv_decomp_GoF_->SaveAs(roostr);

  /////////////////////////////////    
    
  roostr = TString::Format("h2_lambda_transform_matrix_%d_%s", index, postfix.Data());
  TH2D *h2_lambda_transform_matrix = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double value = matrix_eigenvector(ibin-1, jbin-1);
      h2_lambda_transform_matrix->SetBinContent(ibin, jbin, value);
    }
  }

  roostr = TString::Format("canv_h2_lambda_transform_matrix_%d_%s", index, postfix.Data());
  TCanvas *canv_h2_lambda_transform_matrix = new TCanvas(roostr, roostr, 900, 850);
  func_canv_margin(canv_h2_lambda_transform_matrix, 0.15, 0.2,0.15,0.2);
  h2_lambda_transform_matrix->Draw("colz");
  func_xy_title(h2_lambda_transform_matrix, "Normal Bin Index", "Decomposition Bin Index", 1);
  func_title_size(h2_lambda_transform_matrix, 0.05, 0.05, 0.05, 0.05);

  cout<<endl;
  return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_toy_variations(int num_toys)
{
  //cout<<endl<<" ---> Set_toy_variations "<<endl;

  map_matrix_tosc_toy_pred.clear();
  tosc_NUM_TOYS = num_toys;

  ///////
  TMatrixDSym DSmatrix_cov(default_newworld_rows);
  for(int idx=0; idx<default_newworld_rows; idx++) {
    for(int jdx=0; jdx<default_newworld_rows; jdx++) {
      DSmatrix_cov(idx, jdx) = matrix_tosc_eff_newworld_abs_syst_total(idx, jdx);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();

  ///////
  for(int itoy=1; itoy<=num_toys; itoy++) {
    TMatrixD matrix_gaus_rand(default_newworld_rows, 1);
    int eff_count = 0;

  RANDOM_AGAIN:
    eff_count++;
    for(int idx=0; idx<default_newworld_rows; idx++) {
      if( matrix_eigenvalue(idx)>=0 ) {
	matrix_gaus_rand(idx, 0) = tosc_rand->Gaus( 0, sqrt(matrix_eigenvalue(idx)) );// systematics
      }
      else {
	matrix_gaus_rand(idx, 0) = 0;
      }
    }// for(int idx=0; idx<default_newworld_rows; idx++)

    TMatrixD matrix_variation = (matrix_eigenvector*matrix_gaus_rand).T();

    bool flag_negative = 0;

    for(int idx=0; idx<default_newworld_rows; idx++) {
      double val_with_syst = matrix_variation(0, idx) + matrix_tosc_eff_newworld_pred(0, idx);// key point

      if(val_with_syst<0) {
	if( matrix_tosc_eff_newworld_pred(0, idx)<2 ) {// hack for low statistics
	  matrix_variation(0, idx) *= (-1);
	}
	else {
	  flag_negative = 1;
	  break;
	}
      }
    }// for(int idx=0; idx<default_newworld_rows; idx++)

    if( flag_negative ) goto RANDOM_AGAIN;

    TMatrixD matrix_temp_pred(1, default_newworld_rows);
    for(int idx=0; idx<default_newworld_rows; idx++) {
      double val_with_syst = matrix_variation(0, idx) + matrix_tosc_eff_newworld_pred(0, idx);// key point
      val_with_syst = tosc_rand->PoissonD(val_with_syst);// statistics
      matrix_temp_pred(0, idx) = val_with_syst;
    }// for(int idx=0; idx<default_newworld_rows; idx++)
    map_matrix_tosc_toy_pred[itoy].Clear(); map_matrix_tosc_toy_pred[itoy].ResizeTo(1, default_newworld_rows);
    map_matrix_tosc_toy_pred[itoy] = matrix_temp_pred;
    
    // cout<<" ---> generating itoy: "<<itoy<<endl;

  }// for(int itoy=1; itoy<=num_toys; itoy++)

  ///////
  
  // if( 0 ) {// validation
  //   TPrincipal principal_cov(default_newworld_rows, "ND");
  //   double *array_pseudo_expt = new double[default_newworld_rows];
  //   for(int itoy=1; itoy<=num_toys; itoy++) {
  //     if(itoy%1000==0) cout<<TString::Format(" ---> Processing %6d", itoy)<<endl;
  //     for(int idx=0; idx<default_newworld_rows; idx++) {
  // 	array_pseudo_expt[idx] = map_matrix_tosc_toy_pred[itoy](0, idx);
  //     }
  //     principal_cov.AddRow( array_pseudo_expt );
  //   }

  //   TMatrixD *matrix_principal_cov = (TMatrixD*)principal_cov.GetCovarianceMatrix();
  //   for(int idx=0; idx<default_newworld_rows; idx++) {
  //     for(int jdx=0; jdx<idx; jdx++) {
  // 	(*matrix_principal_cov)(jdx,idx) = (*matrix_principal_cov)(idx,jdx);
  //     }      
  //   }

  //   for(int idx=0; idx<default_newworld_rows; idx++) {
  //     cout<<TString::Format("%3d cv %f, old,new: %f %f", idx+1, matrix_tosc_eff_newworld_pred(0, idx),
  // 			    sqrt(matrix_tosc_eff_newworld_abs_syst_total(idx, idx)),
  // 			    sqrt((*matrix_principal_cov)(idx,idx)))<<endl;
  //   }
  // }// if( 0 ) {// validation
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccp

void TOsc::Set_apply_POT()
{
  //////////////////////////

  //cout<<endl<<" ---> Set POT"<<endl<<endl;
  
  ////////////////////////// hack for POT scale, should know which part is BNB, and which part is NuMI
  
  matrix_tosc_eff_newworld_meas = matrix_tosc_default_newworld_meas;

  /////// hhack
  
  // tosc_scaleF_POT_BNB  = 1;  // scale BNB  run1-3 to BNB  run1-3
  // tosc_scaleF_POT_NuMI = 4;  // scale NumI run1   to NuMI run1-3
 
  // tosc_scaleF_POT_BNB  = 2;  // scale BNB  run1-3 to BNB  run1-5
  // tosc_scaleF_POT_NuMI = 5.2;// scale NumI run1   to NuMI run1-5
  
  ///////  
  for(int idx=0; idx<default_oldworld_rows; idx++) matrix_tosc_temp_oldworld_abs_syst_addi(idx, idx) = matrix_tosc_default_oldworld_abs_syst_addi(idx, idx);
  for(int idx=0; idx<default_newworld_rows; idx++) matrix_tosc_temp_newworld_abs_syst_mcstat(idx, idx) = matrix_tosc_default_newworld_abs_syst_mcstat(idx, idx);

  /////// hack, pot-scaling
  // for(int idx=0; idx<default_oldworld_rows; idx++) {
  //   if( idx<26*21 ) {// hack, oldworld
  //     matrix_tosc_oscillation_oldworld_pred(0, idx) *= tosc_scaleF_POT_BNB;
  //     matrix_tosc_temp_oldworld_abs_syst_addi(idx, idx)  *= (tosc_scaleF_POT_BNB*tosc_scaleF_POT_BNB);
  //   }
  //   else {
  //     matrix_tosc_oscillation_oldworld_pred(0, idx) *= tosc_scaleF_POT_NuMI;
  //     matrix_tosc_temp_oldworld_abs_syst_addi(idx, idx)  *= (tosc_scaleF_POT_NuMI*tosc_scaleF_POT_NuMI);
  //   }
  // }
  // /////// hack, pot-scaling
  // for(int idx=0; idx<default_newworld_rows; idx++) {
  //   if( idx<26*7 ) {// hack, oldworld
  //     matrix_tosc_eff_newworld_meas(0, idx) *= tosc_scaleF_POT_BNB;
  //     matrix_tosc_temp_newworld_abs_syst_mcstat(idx, idx) *= (tosc_scaleF_POT_BNB*tosc_scaleF_POT_BNB);
  //   }
  //   else {
  //     matrix_tosc_eff_newworld_meas(0, idx) *= tosc_scaleF_POT_NuMI;
  //     matrix_tosc_temp_newworld_abs_syst_mcstat(idx, idx) *= (tosc_scaleF_POT_NuMI*tosc_scaleF_POT_NuMI);
  //   }
  // }
  
  ///////

  ///////

  matrix_tosc_eff_newworld_pred = matrix_tosc_oscillation_oldworld_pred * matrix_tosc_transform; 
  
  user_matrix_tosc_temp_oscillation_oldworld_pred = matrix_tosc_oscillation_oldworld_pred;
  user_matrix_tosc_temp_oscillation_oldworld_pred_T.Transpose(user_matrix_tosc_temp_oscillation_oldworld_pred); 
  matrix_tosc_temp_oldworld_abs_syst = user_matrix_tosc_temp_oscillation_oldworld_pred_T * user_matrix_tosc_temp_oscillation_oldworld_pred;
  ElementMult(matrix_tosc_temp_oldworld_abs_syst, matrix_tosc_temp_oldworld_rel_syst);// key function

  /////// dirt
  for(int idx=0; idx<default_oldworld_rows; idx++) matrix_tosc_temp_oldworld_abs_syst(idx, idx) += matrix_tosc_temp_oldworld_abs_syst_addi(idx, idx);

  ///////
  
  //matrix_tosc_eff_newworld_abs_syst_total = matrix_tosc_transform_T * matrix_tosc_temp_oldworld_abs_syst * matrix_tosc_transform;

  if( 1 ) {// hhack
    int unit_block = default_oldworld_rows/3/2;// overlay, EXT, OSC ---> BNB and NuMI

    if( 1 ) {// BNB & BNB
      matrix_tosc_AA = matrix_tosc_temp_oldworld_abs_syst.GetSub(0, unit_block-1, 0, unit_block-1);
      matrix_tosc_BB = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block, 2*unit_block-1, unit_block, 2*unit_block-1);
      matrix_tosc_CC = matrix_tosc_temp_oldworld_abs_syst.GetSub(2*unit_block, 3*unit_block-1, 2*unit_block, 3*unit_block-1);
      
      matrix_tosc_AB = matrix_tosc_temp_oldworld_abs_syst.GetSub(0, unit_block-1, unit_block, 2*unit_block-1);
      matrix_tosc_AC = matrix_tosc_temp_oldworld_abs_syst.GetSub(0, unit_block-1, 2*unit_block, 3*unit_block-1);
      matrix_tosc_BC = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block, 2*unit_block-1, 2*unit_block, 3*unit_block-1);
      
      matrix_tosc_BA = matrix_tosc_AB.T(); matrix_tosc_AB.T(); 
      matrix_tosc_CA = matrix_tosc_AC.T(); matrix_tosc_AC.T(); 
      matrix_tosc_CB = matrix_tosc_BC.T(); matrix_tosc_BC.T(); 
      
      matrix_tosc_BNB2BNB.Zero();
      matrix_tosc_BNB2BNB += matrix_tosc_AA;
      matrix_tosc_BNB2BNB += matrix_tosc_BB;
      matrix_tosc_BNB2BNB += matrix_tosc_CC;
      matrix_tosc_BNB2BNB += matrix_tosc_AB;
      matrix_tosc_BNB2BNB += matrix_tosc_AC;
      matrix_tosc_BNB2BNB += matrix_tosc_BC;
      matrix_tosc_BNB2BNB += matrix_tosc_BA;
      matrix_tosc_BNB2BNB += matrix_tosc_CA;
      matrix_tosc_BNB2BNB += matrix_tosc_CB;  
    }
        
    if( 1 ) {// NuMI & NuMI
      matrix_tosc_AA = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block*3+0,            unit_block*3+unit_block-1,   unit_block*3+0,            unit_block*3+unit_block-1);
      matrix_tosc_BB = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block*3+unit_block,   unit_block*3+2*unit_block-1, unit_block*3+unit_block,   unit_block*3+2*unit_block-1);
      matrix_tosc_CC = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block*3+2*unit_block, unit_block*3+3*unit_block-1, unit_block*3+2*unit_block, unit_block*3+3*unit_block-1);
      
      matrix_tosc_AB = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block*3+0,          unit_block*3+unit_block-1,   unit_block*3+unit_block,   unit_block*3+2*unit_block-1);
      matrix_tosc_AC = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block*3+0,          unit_block*3+unit_block-1,   unit_block*3+2*unit_block, unit_block*3+3*unit_block-1);
      matrix_tosc_BC = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block*3+unit_block, unit_block*3+2*unit_block-1, unit_block*3+2*unit_block, unit_block*3+3*unit_block-1);
      
      matrix_tosc_BA = matrix_tosc_AB.T(); matrix_tosc_AB.T(); 
      matrix_tosc_CA = matrix_tosc_AC.T(); matrix_tosc_AC.T(); 
      matrix_tosc_CB = matrix_tosc_BC.T(); matrix_tosc_BC.T(); 
      
      matrix_tosc_NuMI2NuMI.Zero();
      matrix_tosc_NuMI2NuMI += matrix_tosc_AA;
      matrix_tosc_NuMI2NuMI += matrix_tosc_BB;
      matrix_tosc_NuMI2NuMI += matrix_tosc_CC;
      matrix_tosc_NuMI2NuMI += matrix_tosc_AB;
      matrix_tosc_NuMI2NuMI += matrix_tosc_AC;
      matrix_tosc_NuMI2NuMI += matrix_tosc_BC;
      matrix_tosc_NuMI2NuMI += matrix_tosc_BA;
      matrix_tosc_NuMI2NuMI += matrix_tosc_CA;
      matrix_tosc_NuMI2NuMI += matrix_tosc_CB;  
    }

    if( 1 ) {// BNB2NuMI
      matrix_tosc_AA = matrix_tosc_temp_oldworld_abs_syst.GetSub(0, unit_block-1, unit_block*3+0, unit_block*3+unit_block-1);
      matrix_tosc_BB = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block, 2*unit_block-1, unit_block*3+unit_block, unit_block*3+2*unit_block-1);
      matrix_tosc_CC = matrix_tosc_temp_oldworld_abs_syst.GetSub(2*unit_block, 3*unit_block-1, unit_block*3+2*unit_block, unit_block*3+3*unit_block-1);
      
      matrix_tosc_AB = matrix_tosc_temp_oldworld_abs_syst.GetSub(0, unit_block-1, unit_block*3+unit_block, unit_block*3+2*unit_block-1);
      matrix_tosc_AC = matrix_tosc_temp_oldworld_abs_syst.GetSub(0, unit_block-1, unit_block*3+2*unit_block, unit_block*3+3*unit_block-1);
      matrix_tosc_BC = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block, 2*unit_block-1, unit_block*3+2*unit_block, unit_block*3+3*unit_block-1);
      
      matrix_tosc_BA = matrix_tosc_temp_oldworld_abs_syst.GetSub(unit_block, 2*unit_block-1, unit_block*3+0, unit_block*3+unit_block-1);
      matrix_tosc_CA = matrix_tosc_temp_oldworld_abs_syst.GetSub(2*unit_block, 3*unit_block-1, unit_block*3+0, unit_block*3+unit_block-1);
      matrix_tosc_CB =  matrix_tosc_temp_oldworld_abs_syst.GetSub(2*unit_block, 3*unit_block-1, unit_block*3+unit_block, unit_block*3+2*unit_block-1);

      matrix_tosc_BNB2NuMI.Zero();
      matrix_tosc_BNB2NuMI += matrix_tosc_AA;
      matrix_tosc_BNB2NuMI += matrix_tosc_BB;
      matrix_tosc_BNB2NuMI += matrix_tosc_CC;
      matrix_tosc_BNB2NuMI += matrix_tosc_AB;
      matrix_tosc_BNB2NuMI += matrix_tosc_AC;
      matrix_tosc_BNB2NuMI += matrix_tosc_BC;
      matrix_tosc_BNB2NuMI += matrix_tosc_BA;
      matrix_tosc_BNB2NuMI += matrix_tosc_CA;
      matrix_tosc_BNB2NuMI += matrix_tosc_CB;  

      matrix_tosc_NuMI2BNB = matrix_tosc_BNB2NuMI.T(); matrix_tosc_BNB2NuMI.T();
    }

    //matrix_tosc_temp_total_in_POT.Zero();
    matrix_tosc_temp_total_in_POT.SetSub(0, 0, matrix_tosc_BNB2BNB);
    matrix_tosc_temp_total_in_POT.SetSub(unit_block, unit_block, matrix_tosc_NuMI2NuMI);
    matrix_tosc_temp_total_in_POT.SetSub(0, unit_block, matrix_tosc_BNB2NuMI);
    matrix_tosc_temp_total_in_POT.SetSub(unit_block, 0, matrix_tosc_NuMI2BNB);

    matrix_tosc_eff_newworld_abs_syst_total = matrix_tosc_temp_total_in_POT;    
  }
 
  /////// mc_stat
  for(int idx=0; idx<default_newworld_rows; idx++) matrix_tosc_eff_newworld_abs_syst_total(idx, idx) += matrix_tosc_temp_newworld_abs_syst_mcstat(idx, idx);

  /////// protect
  for(int idx=0; idx<default_newworld_rows; idx++) {    
    if( matrix_tosc_eff_newworld_abs_syst_total(idx, idx)==0 ) matrix_tosc_eff_newworld_abs_syst_total(idx, idx) = 1e-6;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

double TOsc::Prob_oscillaion(double Etrue, double baseline, int strflag_osc)// ooo
{  
  double prob = 0;

  // dm2*L/4E = 1.267 * dm2(eV2) * L(m) / E(MeV)  

  ///// y = 4 * x * (1-x)
  ///// x1 = ( 1 - sqrt(1-y) )/2;
  ///// x2 = ( 1 + sqrt(1-y) )/2;
  ///// if( fabs(1-y)<1e-6 ) y = 1;  

  ///////////////////////////////////////////////// dm2 vs. sin2_2Tue vs. t24
  ///////////////////////////////////////////////// dm2 vs. sin2_2Tue vs. t24
  ///////// sin2_2Tue <= sin2_T24

  // double user_sin2_2Tue = tosc_sin2_2theta_14;
  // double user_sin2_2T14 = user_sin2_2Tue/tosc_sin2_theta_24;

  // double y = user_sin2_2T14;  
  // double x = ( 1 + sqrt(1-y) )/2;// solution of >1/2
  // //double x = ( 1 - sqrt(1-y) )/2;// solution of <1/2

  // double effective_sin2_theta_14  = x;
  // double effective_cos2_theta_14  = 1 - effective_sin2_theta_14;
  // double effective_sin2_2theta_14 = y;

  ///////////////////////////////////////////////// dm2 vs. sin2_2Tee vs. t24
  ///////////////////////////////////////////////// dm2 vs. sin2_2Tee vs. t24
  
  // double user_sin2_2T14 = tosc_sin2_2theta_14;

  // double y = user_sin2_2T14;  
  // double x = ( 1 + sqrt(1-y) )/2;// solution of >1/2
  // //double x = ( 1 - sqrt(1-y) )/2;// solution of <1/2

  // double effective_sin2_theta_14  = x;
  // double effective_cos2_theta_14  = 1 - effective_sin2_theta_14;
  // double effective_sin2_2theta_14 = user_sin2_2T14;

  ///////////////////////////////////////// dm2 vs. t14 vs. t24
  ///////////////////////////////////////// dm2 vs. t14 vs. t24

  double effective_sin2_theta_14  = tosc_sin2_2theta_14;
  double effective_cos2_theta_14  = 1 - effective_sin2_theta_14;
  double effective_sin2_2theta_14 = 4 * effective_sin2_theta_14 * effective_cos2_theta_14;  
  
  /////////////////////////////////////////nominal

  double sin_Delta = sin( 1.267 * tosc_dm2_41 * baseline/Etrue );

  ///////
  double sin2_Delta = sin_Delta * sin_Delta;
  
  const int nue2nue   = 1;
  const int numu2numu = 2;
  const int numu2nue  = 3;
  const int nue2numu  = 4;
  const int nueNC     = 5;
  const int numuNC    = 6;
  
  int flag_osc = strflag_osc;    
    
  switch( flag_osc ) {
  case nue2nue:
    prob = 1 - effective_sin2_2theta_14 * sin2_Delta;
    //prob = 1;
    //prob = 1 - tosc_sin2_2theta_14 * sin2_Delta;
    break;
  case numu2numu:
    prob = 1 - 4*effective_cos2_theta_14*tosc_sin2_theta_24 * (1 - effective_cos2_theta_14*tosc_sin2_theta_24) * sin2_Delta;
    //prob = 1;
    //prob = 1 - tosc_sin2_2theta_14 * sin2_Delta;
    break;
  case numu2nue:
    prob = effective_sin2_2theta_14 * tosc_sin2_theta_24 * sin2_Delta;
    //prob = 0;
    //prob = 1;
    //prob = tosc_sin2_2theta_14 * sin2_Delta;
    break;
  case nue2numu:
    break;
  case nueNC:
    prob = 1 - effective_sin2_2theta_14 * ( 1-tosc_sin2_theta_24 ) * sin2_Delta;// theta_34 = 0
    //prob = 1;
    break;
  case numuNC:
    prob = 1 - (effective_cos2_theta_14*effective_cos2_theta_14) * (4*tosc_sin2_theta_24*(1-tosc_sin2_theta_24)) * sin2_Delta;// theta_34 = 0
    //prob = 1;
    //prob = 1 - tosc_sin2_2theta_14 * sin2_Delta;
    break;
  default:
    cerr<<"ERROR: NAN flag_osc"<<endl; exit(1);
  }  
  
  return prob;
}

///////////////////

void TOsc::Set_oscillation_base_minus(vector<double> *vec_ratioPOT, vector< vector<EventInfo> > *vec_vec_eventinfo, int pred_channel_index, TString str_osc_mode)
{
  int total_pred_chs = map_tosc_default_h1d_pred.size();
  if( pred_channel_index > total_pred_chs ) { cerr<<TString::Format(" ERROR: pred_channel_index(%d) > total_pred_chs(%d)", pred_channel_index, total_pred_chs)<<endl; exit(1); }  
  cout<<TString::Format("            ---> (self-check) work on pred_channel: %3d,  oscillation mode: %-10s", pred_channel_index, str_osc_mode.Data())<<endl;

  
  for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ ) {
    TH1D *h1d_temp = (TH1D*)map_tosc_default_h1d_pred[pred_channel_index]->Clone("h1d_temp"); h1d_temp->Reset();// define and clear
    
    for(int ievent=0; ievent<(int)vec_vec_eventinfo->at(isize).size(); ievent++) {
      EventInfo info = vec_vec_eventinfo->at(isize).at(ievent);
      h1d_temp->Fill(info.e2e_Ereco, info.e2e_weight_xs);
    }

    h1d_temp->Scale( vec_ratioPOT->at(isize) );

    ///////    
    TMatrixD matrix_temp(1, default_oldworld_rows);        
    int bin_index_base = 0;
    if( pred_channel_index > 1 ) { for(int ich=1; ich<pred_channel_index; ich++) bin_index_base += ( map_tosc_default_h1d_pred[ich]->GetNbinsX()+1 ); }
    for(int ibin=1; ibin<=h1d_temp->GetNbinsX()+1; ibin++) { matrix_temp(0, bin_index_base + ibin -1) = h1d_temp->GetBinContent(ibin); }
    matrix_tosc_oscillation_base_oldworld_pred -= matrix_temp;

    delete h1d_temp;
  }// for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ )  
}

void TOsc::Set_oscillation_base_added(vector<double> *vec_ratioPOT, vector< vector<EventInfo> > *vec_vec_eventinfo, int pred_channel_index, TString str_osc_mode)
{
  int total_pred_chs = map_tosc_default_h1d_pred.size();
  if( pred_channel_index > total_pred_chs ) { cerr<<TString::Format(" ERROR: pred_channel_index(%d) > total_pred_chs(%d)", pred_channel_index, total_pred_chs)<<endl; exit(1); }  

  
  // for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ ) {
  //   TH1D *h1d_temp = (TH1D*)map_tosc_default_h1d_pred[pred_channel_index]->Clone("h1d_temp"); h1d_temp->Reset();// define and clear
    
  //   for(int ievent=0; ievent<(int)vec_vec_eventinfo->at(isize).size(); ievent++) {
  //     EventInfo info = vec_vec_eventinfo->at(isize).at(ievent);
  //     double prob = Prob_oscillaion(info.e2e_Etrue, info.e2e_baseline, str_osc_mode);
  //     h1d_temp->Fill(info.e2e_Ereco, prob * info.e2e_weight_xs);
  //   }

  //   h1d_temp->Scale( vec_ratioPOT->at(isize) );

  //   ///////    
  //   TMatrixD matrix_temp(1, default_oldworld_rows);        
  //   int bin_index_base = 0;
  //   if( pred_channel_index > 1 ) { for(int ich=1; ich<pred_channel_index; ich++) bin_index_base += ( map_tosc_default_h1d_pred[ich]->GetNbinsX()+1 ); }
  //   for(int ibin=1; ibin<=h1d_temp->GetNbinsX()+1; ibin++) { matrix_temp(0, bin_index_base + ibin -1) = h1d_temp->GetBinContent(ibin); }
  //   matrix_tosc_oscillation_oldworld_pred += matrix_temp;
    
  //   delete h1d_temp;
  // }// for(int isize=0; isize<(int)vec_vec_eventinfo->size(); isize++ )


  /////// xxx

  int flag_osc = -1;
  if( str_osc_mode=="nue2nue" )        flag_osc = 1;
  else if( str_osc_mode=="numu2numu" ) flag_osc = 2;
  else if( str_osc_mode=="numu2nue" )  flag_osc = 3;
  else if( str_osc_mode=="nue2numu" )  flag_osc = 4;
  else if( str_osc_mode=="nueNC")      flag_osc = 5;
  else if( str_osc_mode=="numuNC")     flag_osc = 6;
  else flag_osc = -1;  
  
  if( flag_osc==-1 ) {
    cerr<<" ERROR: flag_osc==-1, str_osc_mode=="<<str_osc_mode<<" ---> Candidates: nue2nue numu2numu numu2nue nue2numu nueNC numuNC"<<endl;
    exit(1);
  }
  
  int line_user = -1;
  for( auto it_first=vec_vec_eventinfo->begin(); it_first!=vec_vec_eventinfo->end(); it_first++ ) {
    
    line_user++;
    TH1D *h1d_temp = (TH1D*)map_tosc_default_h1d_pred[pred_channel_index]->Clone("h1d_temp"); h1d_temp->Reset();// define and clear
    
    for( auto it_second=it_first->begin(); it_second!=it_first->end(); it_second++ ) {      
      EventInfo info = (*it_second);
      double prob = Prob_oscillaion(info.e2e_Etrue, info.e2e_baseline, flag_osc);
      h1d_temp->Fill(info.e2e_Ereco, prob * info.e2e_weight_xs);
    }
    
    h1d_temp->Scale( vec_ratioPOT->at(line_user) );
    
    ///////
    TMatrixD matrix_temp(1, default_oldworld_rows);        
    int bin_index_base = 0;
    if( pred_channel_index > 1 ) { for(int ich=1; ich<pred_channel_index; ich++) bin_index_base += ( map_tosc_default_h1d_pred[ich]->GetNbinsX()+1 ); }
    for(int ibin=1; ibin<=h1d_temp->GetNbinsX()+1; ibin++) { matrix_temp(0, bin_index_base + ibin -1) = h1d_temp->GetBinContent(ibin); }
    matrix_tosc_oscillation_oldworld_pred += matrix_temp;    
    
    delete h1d_temp;
    
  }  

}

///////////////////
  
void TOsc::Apply_oscillation()
{
  //cout<<endl<<" ---> Apply_oscillation"<<endl;

  matrix_tosc_oscillation_oldworld_pred = matrix_tosc_oscillation_base_oldworld_pred;

  /////////////////// same as the ones in void TOsc::Set_oscillation_base(), but use "added" instead of "minus"

  if( flag_apply_oscillation_NuMI ) {

    if( flag_NuMI_nueCC_from_intnue ) {
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo, 22, "nue2nue");// hack
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo, 23, "nue2nue");// hack
    }
  
    if( flag_NuMI_nueCC_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo, 22, "numu2numu");// hack
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo, 23, "numu2numu");// hack
    }
    if( flag_NuMI_numuCC_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo, 24, "numu2numu");// hack
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo, 25, "numu2numu");// hack
    }
    if( flag_NuMI_CCpi0_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo, 26, "numu2numu");// hack
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo, 27, "numu2numu");// hack
    }
    if( flag_NuMI_NCpi0_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo, 28, "numu2numu");// hack
    }


    if( flag_NuMI_nueCC_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo, 22, "nueNC");// hack
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo, 23, "nueNC");// hack
    }
    if( flag_NuMI_nueCC_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo, 22, "numuNC");// hack
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo, 23, "numuNC");// hack
    }
    if( flag_NuMI_numuCC_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo, 24, "nueNC");// hack
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo, 25, "nueNC");// hack
    }
    if( flag_NuMI_numuCC_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo, 24, "numuNC");// hack
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo, 25, "numuNC");// hack
    }
    if( flag_NuMI_CCpi0_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo, 26, "nueNC");// hack
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo, 27, "nueNC");// hack
    }
    if( flag_NuMI_CCpi0_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo, 26, "numuNC");// hack
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo, 27, "numuNC");// hack
    }
    if( flag_NuMI_NCpi0_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo, 28, "nueNC");// hack
    }  
    if( flag_NuMI_NCpi0_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo, 28, "numuNC");// hack
    }  

    if( flag_NuMI_nueCC_from_appnue ) {
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_appnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo, 36, "numu2nue");// hack
      Set_oscillation_base_added(&vector_NuMI_nueCC_from_appnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo, 37, "numu2nue");// hack
    }
    if( flag_NuMI_numuCC_from_appnue ) {
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_appnue_scaleFPOT, &vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo, 38, "numu2nue");// hack
      Set_oscillation_base_added(&vector_NuMI_numuCC_from_appnue_scaleFPOT, &vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo, 39, "numu2nue");// hack
    }
    if( flag_NuMI_CCpi0_from_appnue ) {
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_appnue_scaleFPOT, &vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo, 40, "numu2nue");// hack
      Set_oscillation_base_added(&vector_NuMI_CCpi0_from_appnue_scaleFPOT, &vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo, 41, "numu2nue");// hack
    }
    if( flag_NuMI_NCpi0_from_appnue ) {
      Set_oscillation_base_added(&vector_NuMI_NCpi0_from_appnue_scaleFPOT, &vector_vector_NuMI_NCpi0_from_appnue_eventinfo, 42, "numu2nue");// hack
    }
  
  }// if( flag_apply_oscillation_NuMI )

  /////////
  /////////
  /////////
  
  if( flag_apply_oscillation_BNB ) {

    if( flag_BNB_nueCC_from_intnue ) {
      Set_oscillation_base_added(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo, 1, "nue2nue");// hack
      Set_oscillation_base_added(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo, 2, "nue2nue");// hack
    }
  
    if( flag_BNB_nueCC_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo, 1, "numu2numu");// hack
      Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo, 2, "numu2numu");// hack
    }
    if( flag_BNB_numuCC_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo, 3, "numu2numu");// hack
      Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo, 4, "numu2numu");// hack
    }
    if( flag_BNB_CCpi0_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo, 5, "numu2numu");// hack
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo, 6, "numu2numu");// hack
    }
    if( flag_BNB_NCpi0_from_overlaynumu ) {
      Set_oscillation_base_added(&vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo, 7, "numu2numu");// hack
    }

    
    if( flag_BNB_nueCC_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynueNC_FC_eventinfo, 1, "nueNC");// hack
      Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynueNC_PC_eventinfo, 2, "nueNC");// hack
    }
    if( flag_BNB_nueCC_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumuNC_FC_eventinfo, 1, "numuNC");// hack
      Set_oscillation_base_added(&vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumuNC_PC_eventinfo, 2, "numuNC");// hack
    }
    if( flag_BNB_numuCC_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynueNC_FC_eventinfo, 3, "nueNC");// hack
      Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynueNC_PC_eventinfo, 4, "nueNC");// hack
    }
    if( flag_BNB_numuCC_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumuNC_FC_eventinfo, 3, "numuNC");// hack
      Set_oscillation_base_added(&vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumuNC_PC_eventinfo, 4, "numuNC");// hack
    }
    if( flag_BNB_CCpi0_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynueNC_FC_eventinfo, 5, "nueNC");// hack
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynueNC_PC_eventinfo, 6, "nueNC");// hack
    }
    if( flag_BNB_CCpi0_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumuNC_FC_eventinfo, 5, "numuNC");// hack
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumuNC_PC_eventinfo, 6, "numuNC");// hack
    }
    if( flag_BNB_NCpi0_from_overlaynueNC ) {
      Set_oscillation_base_added(&vector_BNB_NCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynueNC_eventinfo, 7, "nueNC");// hack
    }
    if( flag_BNB_NCpi0_from_overlaynumuNC ) {
      Set_oscillation_base_added(&vector_BNB_NCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynumuNC_eventinfo, 7, "numuNC");// hack
    }
 

    if( flag_BNB_nueCC_from_appnue ) {
      Set_oscillation_base_added(&vector_BNB_nueCC_from_appnue_scaleFPOT, &vector_vector_BNB_nueCC_from_appnue_FC_eventinfo, 15, "numu2nue");// hack
      Set_oscillation_base_added(&vector_BNB_nueCC_from_appnue_scaleFPOT, &vector_vector_BNB_nueCC_from_appnue_PC_eventinfo, 16, "numu2nue");// hack
    }
    if( flag_BNB_numuCC_from_appnue ) {
      Set_oscillation_base_added(&vector_BNB_numuCC_from_appnue_scaleFPOT, &vector_vector_BNB_numuCC_from_appnue_FC_eventinfo, 17, "numu2nue");// hack
      Set_oscillation_base_added(&vector_BNB_numuCC_from_appnue_scaleFPOT, &vector_vector_BNB_numuCC_from_appnue_PC_eventinfo, 18, "numu2nue");// hack
    }
    if( flag_BNB_CCpi0_from_appnue ) {
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_appnue_scaleFPOT, &vector_vector_BNB_CCpi0_from_appnue_FC_eventinfo, 19, "numu2nue");// hack
      Set_oscillation_base_added(&vector_BNB_CCpi0_from_appnue_scaleFPOT, &vector_vector_BNB_CCpi0_from_appnue_PC_eventinfo, 20, "numu2nue");// hack
    }
    if( flag_BNB_NCpi0_from_appnue ) {
      Set_oscillation_base_added(&vector_BNB_NCpi0_from_appnue_scaleFPOT, &vector_vector_BNB_NCpi0_from_appnue_eventinfo, 21, "numu2nue");// hack
    }

  }// if( flag_apply_oscillation_BNB )
    
  /////////
  /////////
  /////////

  for(int idx=0; idx<default_oldworld_rows; idx++) {

    if( matrix_tosc_oscillation_oldworld_pred(0, idx)<-1 ) {
      cout<<endl;
      cout<<"Warning: matrix_tosc_oscillation_oldworld_pred(0, idx)<-1 at index "<<idx<<endl;
      cout<<endl;      
    }

    if( matrix_tosc_oscillation_oldworld_pred(0, idx)<1e-4 ) matrix_tosc_oscillation_oldworld_pred(0, idx) = 0;// protect

  }

  ///////// winxp check with results from framework
  ///////// winxp check with results from framework
  ///////// winxp check with results from framework

  // if( 0 ) {// self-check
  //   int idx_aa = 26*0;
  //   int idx_bb = 26*21;    
  //   for(int idx=1; idx<=26*21; idx++) {
  //     cout<<TString::Format("%3d   %15.6f ---> (origin) %15.6f  (diff) %9.6f,    %15.6f  ---> (origin) %15.6f  (diff) %9.6f", idx,
  // 			    matrix_tosc_oscillation_oldworld_pred(0, idx_aa  + idx-1),  matrix_tosc_default_oldworld_pred(0, idx_aa  + idx-1),
  // 			    matrix_tosc_oscillation_oldworld_pred(0, idx_aa  + idx-1) - matrix_tosc_default_oldworld_pred(0, idx_aa  + idx-1),
  // 			    matrix_tosc_oscillation_oldworld_pred(0, idx_bb  + idx-1),  matrix_tosc_default_oldworld_pred(0, idx_bb  + idx-1),
  // 			    matrix_tosc_oscillation_oldworld_pred(0, idx_bb  + idx-1) - matrix_tosc_default_oldworld_pred(0, idx_bb  + idx-1)
  // 			    )<<endl;
  //   }
  // }// if( 0 ) {// self-check
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_oscillation_base_subfunc(TString strfile_mcPOT, TString strfile_dataPOT, vector<double> *vec_ratioPOT, TString strfile_mc_e2e, TString str_treename, vector< vector<EventInfo> > *vec_vec_eventinfo)
{
  // Declaration of leaf types
  Int_t           e2e_pdg;
  //Int_t           e2e_flag_FC;
  Double_t        e2e_Etrue;
  Double_t        e2e_Ereco;
  Double_t        e2e_weight_xs;
  Double_t        e2e_baseline;

  // List of branches
  TBranch        *b_e2e_pdg;   //!
  //TBranch        *b_e2e_flag_FC;   //!
  TBranch        *b_e2e_Etrue;   //!
  TBranch        *b_e2e_Ereco;   //!
  TBranch        *b_e2e_weight_xs;   //!
  TBranch        *b_e2e_baseline;   //!
  
  TFile *roofile_obj = new TFile(strfile_mc_e2e, "read");
  TTree *tree_obj = (TTree*)roofile_obj->Get(str_treename);
  if(!tree_obj) { cerr<<" ERROR: Set_oscillation_base_subfunc, no treename: "<<str_treename<<endl; exit(1); }
  
  // Set branch addresses and branch pointers
  tree_obj->SetBranchAddress("e2e_pdg", &e2e_pdg, &b_e2e_pdg);
  //tree_obj->SetBranchAddress("e2e_flag_FC", &e2e_flag_FC, &b_e2e_flag_FC);
  tree_obj->SetBranchAddress("e2e_Etrue", &e2e_Etrue, &b_e2e_Etrue);
  tree_obj->SetBranchAddress("e2e_Ereco", &e2e_Ereco, &b_e2e_Ereco);
  tree_obj->SetBranchAddress("e2e_weight_xs", &e2e_weight_xs, &b_e2e_weight_xs);
  tree_obj->SetBranchAddress("e2e_baseline", &e2e_baseline, &b_e2e_baseline);
      
  int entries = tree_obj->GetEntries();

  if( vec_ratioPOT!=NULL ) cout<<endl;
  cout<<TString::Format("            ---> entries %10d     %50s   --> %20s", entries, strfile_mc_e2e.Data(), str_treename.Data())<<endl;

  vector<EventInfo>vector_eventinfo;

  //if( strfile_mcPOT.Contains("numi") ) cout<<" ################ "<<strfile_mc_e2e<<endl;

  for(int ientry=0; ientry<entries; ientry++) {
    tree_obj->GetEntry(ientry);

    //if( strfile_mcPOT.Contains("numi") && strfile_mc_e2e.Contains("appnue")) e2e_baseline = 680;// hhack
    //e2e_baseline = 680;

    EventInfo eventinfo;
    eventinfo.e2e_pdg = e2e_pdg;
    //eventinfo.e2e_flag_FC = e2e_flag_FC;
    eventinfo.e2e_Etrue = e2e_Etrue;
    eventinfo.e2e_Ereco = e2e_Ereco;
    eventinfo.e2e_weight_xs = e2e_weight_xs;
    eventinfo.e2e_baseline = e2e_baseline;
    vector_eventinfo.push_back(eventinfo);
  }
    
  delete tree_obj;
  delete roofile_obj;
  
  vec_vec_eventinfo->push_back( vector_eventinfo );

  if( vec_ratioPOT!=NULL ) {
    Double_t        pot;
    TBranch        *b_pot;   //!  
    double mc_pot = 1; double data_pot = 0;

    TFile *roofile_mc = new TFile(strfile_mcPOT, "read");
    TTree *tree_mc = (TTree*)roofile_mc->Get("T");
    tree_mc->SetBranchAddress("pot", &pot, &b_pot);
    tree_mc->GetEntry(0); mc_pot = pot;
    delete tree_mc;
    delete roofile_mc;
      
    TFile *roofile_data = new TFile(strfile_dataPOT, "read");
    TTree *tree_data = (TTree*)roofile_data->Get("T");
    tree_data->SetBranchAddress("pot", &pot, &b_pot);
    tree_data->GetEntry(0); data_pot = pot;
    delete tree_data;
    delete roofile_data;
    
    vec_ratioPOT->push_back( data_pot/mc_pot );// from the output of ./convert_histo.pl
    
    cout<<"            ---> MC POT "<<mc_pot<<"\t"<<strfile_mcPOT<<endl;
    cout<<"            ---> DD POT "<<data_pot<<"\t"<<strfile_dataPOT<<endl;      
  }// if( vec_ratioPOT!=NULL )
  
}

///////////////////
  
void TOsc::Set_oscillation_base(TString evetlist_dir)
{
  /// Eventlist generator
  /// /home/xji/data0/work/work_oscillation/301_Framework_for_Osc/wcp-uboone-bdt/apps/convert_checkout_hist.cxx, winxp
  
  cout<<endl;
  cout<<" ---> Set_oscillation_base"<<endl;

  matrix_tosc_oscillation_base_oldworld_pred = matrix_tosc_default_oldworld_pred;

  TString str_dirbase = evetlist_dir;
  
  /////////////////// reference: docDB-38700
  /////////////////// reference: docDB-38700
  /////////////////// reference: docDB-38700
  ///////////////////
  ///////////////////
  ///////////////////

  if( flag_NuMI_nueCC_from_intnue ) {
    cout<<endl<<"      ---> flag_NuMI_nueCC_from_intnue"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_intrinsic_nue_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_intrinsic_nue_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_intrinsic_nue_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_intrinsic_nue_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_intrinsic_nue_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_FC_eventinfo, 22, "nue2nue");// hack
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_intnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_intnue_PC_eventinfo, 23, "nue2nue");// hack
  
  }// if( flag_NuMI_nueCC_from_intnue )

  ///////////////////
  ///////////////////  

  if( flag_NuMI_nueCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_nueCC_from_overlaynumu"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_FC_eventinfo, 22, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumu_PC_eventinfo, 23, "numu2numu");// hack
  
  }// if( flag_NuMI_nueCC_from_overlaynumu )
  
  ///////////////////  

  if( flag_NuMI_nueCC_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_NuMI_nueCC_from_overlaynueNC"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynueNC_FC_eventinfo, 22, "nueNC");// hack
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynueNC_PC_eventinfo, 23, "nueNC");// hack
  
  }// if( flag_NuMI_nueCC_from_overlaynueNC )
   
  ///////////////////  

  if( flag_NuMI_nueCC_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_NuMI_nueCC_from_overlaynumuNC"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo);
    }
 
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumuNC_FC_eventinfo, 22, "numuNC");// hack
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_nueCC_from_overlaynumuNC_PC_eventinfo, 23, "numuNC");// hack
  
  }// if( flag_NuMI_nueCC_from_overlaynumuNC )
  
  ///////////////////  

  if( flag_NuMI_numuCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_numuCC_from_overlaynumu"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_FC_eventinfo, 24, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumu_PC_eventinfo, 25, "numu2numu");// hack
  
  }// if( flag_NuMI_numuCC_from_overlaynumu )
    
  ///////////////////  

  if( flag_NuMI_numuCC_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_NuMI_numuCC_from_overlaynueNC"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynueNC_FC_eventinfo, 24, "nueNC");// hack
    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynueNC_PC_eventinfo, 25, "nueNC");// hack
  
  }// if( flag_NuMI_numuCC_from_overlaynueNC )
   
  ///////////////////  

  if( flag_NuMI_numuCC_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_NuMI_numuCC_from_overlaynumuNC"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                             strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                             strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                             strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                             strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root";   
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                             strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumuNC_FC_eventinfo, 24, "numuNC");// hack
    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_numuCC_from_overlaynumuNC_PC_eventinfo, 25, "numuNC");// hack
  
  }// if( flag_NuMI_numuCC_from_overlaynumuNC )
  
  ///////////////////  

  if( flag_NuMI_CCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_CCpi0_from_overlaynumu"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_FC_eventinfo, 26, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumu_PC_eventinfo, 27, "numu2numu");// hack
  
  }// if( flag_NuMI_CCpi0_from_overlaynumu )
      
  ///////////////////  

  if( flag_NuMI_CCpi0_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_NuMI_CCpi0_from_overlaynueNC"<<endl;        
    {// run1 fhc    
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo);
    }
    {// run1 rhc    
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo);
    }
    {// run2 fhc    
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo);
    }
    {// run2 rhc    
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo);
    }
    {// run3 rhc    
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynueNC_FC_eventinfo, 26, "nueNC");// hack
    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynueNC_PC_eventinfo, 27, "nueNC");// hack
  
  }// if( flag_NuMI_CCpi0_from_overlaynueNC )
   
  ///////////////////  

  if( flag_NuMI_CCpi0_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_NuMI_CCpi0_from_overlaynumuNC"<<endl;        
    {// run1 fhc  
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo);
    }
    {// run1 rhc  
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo);
    }
    {// run2 fhc  
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo);
    }
    {// run2 rhc  
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo);
    }
    {// run3 rhc  
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_FC_eventinfo, 26, "numuNC");// hack
    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_CCpi0_from_overlaynumuNC_PC_eventinfo, 27, "numuNC");// hack
  
  }// if( flag_NuMI_CCpi0_from_overlaynumuNC )
  
  ///////////////////  

  if( flag_NuMI_NCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_NuMI_NCpi0_from_overlaynumu"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynumu_eventinfo, 28, "numu2numu");// hack
  
  }// if( flag_NuMI_NCpi0_from_overlaynumu )
    
  ///////////////////  

  if( flag_NuMI_NCpi0_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_NuMI_NCpi0_from_overlaynueNC"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_NCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynueNC_eventinfo, 28, "nueNC");// hack
  
  }// if( flag_NuMI_NCpi0_from_overlaynueNC )
          
  ///////////////////  

  if( flag_NuMI_NCpi0_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_NuMI_NCpi0_from_overlaynumuNC"<<endl;        
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_nu_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_nu_overlay.root"; 
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_NCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_NuMI_NCpi0_from_overlaynumuNC_eventinfo, 28, "numuNC");// hack
  
  }// if( flag_NuMI_NCpi0_from_overlaynumuNC )
      
  ///////////////////  

  if( flag_NuMI_nueCC_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_NuMI_nueCC_from_appnue"<<endl;
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo);      
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo);      
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo);      
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo);      
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo);      
    }

    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_appnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_appnue_FC_eventinfo, 36, "numu2nue");// hack
    Set_oscillation_base_minus(&vector_NuMI_nueCC_from_appnue_scaleFPOT, &vector_vector_NuMI_nueCC_from_appnue_PC_eventinfo, 37, "numu2nue");// hack
    
  }// if( flag_NuMI_nueCC_from_appnue )
       
  ///////////////////  

  if( flag_NuMI_numuCC_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_NuMI_numuCC_from_appnue"<<endl;
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                      strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo);      
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                      strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo);      
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                      strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo);      
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                      strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo);      
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                      strfile_mc_e2e, str_treename, &vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo);      
    }

    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_appnue_scaleFPOT, &vector_vector_NuMI_numuCC_from_appnue_FC_eventinfo, 38, "numu2nue");// hack
    Set_oscillation_base_minus(&vector_NuMI_numuCC_from_appnue_scaleFPOT, &vector_vector_NuMI_numuCC_from_appnue_PC_eventinfo, 39, "numu2nue");// hack
    
  }// if( flag_NuMI_numuCC_from_appnue )
         
  ///////////////////  

  if( flag_NuMI_CCpi0_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_NuMI_CCpi0_from_appnue"<<endl;
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo);      
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo);      
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo);      
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo);      
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo);      
    }

    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_appnue_scaleFPOT, &vector_vector_NuMI_CCpi0_from_appnue_FC_eventinfo, 40, "numu2nue");// hack
    Set_oscillation_base_minus(&vector_NuMI_CCpi0_from_appnue_scaleFPOT, &vector_vector_NuMI_CCpi0_from_appnue_PC_eventinfo, 41, "numu2nue");// hack
    
  }// if( flag_NuMI_CCpi0_from_appnue )
           
  ///////////////////  

  if( flag_NuMI_NCpi0_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_NuMI_NCpi0_from_appnue"<<endl;
    {// run1 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_appnue_eventinfo);
    }
    {// run1 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run1_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run1_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run1_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_appnue_eventinfo);
    }
    {// run2 fhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_fhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_fhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_FHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_appnue_eventinfo);
    }
    {// run2 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run2_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run2_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run2_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_appnue_eventinfo);
    }
    {// run3 rhc
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_run3_rhc_fullosc_overlay.root";
      TString strfile_dataPOT = str_dirbase + "run3_rhc_data_numi.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_NuMI_run3_RHC_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_NuMI_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_NuMI_NCpi0_from_appnue_eventinfo);
    }

    Set_oscillation_base_minus(&vector_NuMI_NCpi0_from_appnue_scaleFPOT, &vector_vector_NuMI_NCpi0_from_appnue_eventinfo, 42, "numu2nue");// hack
    
  }// if( flag_NuMI_NCpi0_from_appnue )
         
  ///////////////////////////////////////////////////////// reference: docDB-38700
  ///////////////////////////////////////////////////////// reference: docDB-38700    
  ///////////////////////////////////////////////////////// reference: docDB-38700    
  /////////////////////////////////////////////////////////    
  /////////////////////////////////////////////////////////    
  /////////////////////////////////////////////////////////    

  if( flag_BNB_nueCC_from_intnue ) {
    cout<<endl<<"      ---> flag_BNB_nueCC_from_intnue"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_intrinsic_nue_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo);
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_intrinsic_nue_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo);
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_intrinsic_nue_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_intrinsic.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_intnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_intnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_intnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_FC_eventinfo, 1, "nue2nue");// hack
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_intnue_scaleFPOT, &vector_vector_BNB_nueCC_from_intnue_PC_eventinfo, 2, "nue2nue");// hack
    
  }// if( flag_BNB_nueCC_from_intnue )

  ///////////////////
  ///////////////////  

  if( flag_BNB_nueCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_nueCC_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_FC_eventinfo, 1, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumu_PC_eventinfo, 2, "numu2numu");// hack
    
  }// if( flag_BNB_numuCC_from_overlaynumu )
  
  ///////////////////  

  if( flag_BNB_nueCC_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_BNB_nueCC_from_overlaynueNC"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynueNC_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynueNC_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynueNC_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynueNC_FC_eventinfo, 1, "nueNC");// hack
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynueNC_PC_eventinfo, 2, "nueNC");// hack
    
  }// if( flag_BNB_nueCC_from_overlaynueNC )
  
  ///////////////////  

  if( flag_BNB_nueCC_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_BNB_nueCC_from_overlaynumuNC"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumuNC_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumuNC_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_nueCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_overlaynumuNC_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumuNC_FC_eventinfo, 1, "numuNC");// hack
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_nueCC_from_overlaynumuNC_PC_eventinfo, 2, "numuNC");// hack
    
  }// if( flag_BNB_nueCC_from_overlaynumuNC )
  
  ///////////////////  

  if( flag_BNB_numuCC_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_numuCC_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_FC_eventinfo, 3, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynumu_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumu_PC_eventinfo, 4, "numu2numu");// hack
    
  }// if( flag_BNB_numuCC_from_overlaynumu )
  
  ///////////////////  

  if( flag_BNB_numuCC_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_BNB_numuCC_from_overlaynueNC"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynueNC_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynueNC_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynueNC_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynueNC_FC_eventinfo, 3, "nueNC");// hack
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynueNC_PC_eventinfo, 4, "nueNC");// hack
    
  }// if( flag_BNB_numuCC_from_overlaynueNC )
  
  ///////////////////  

  if( flag_BNB_numuCC_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_BNB_numuCC_from_overlaynumuNC"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumuNC_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumuNC_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_numuCC_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                            strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_overlaynumuNC_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumuNC_FC_eventinfo, 3, "numuNC");// hack
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_numuCC_from_overlaynumuNC_PC_eventinfo, 4, "numuNC");// hack
    
  }// if( flag_BNB_numuCC_from_overlaynumuNC )
  
  ///////////////////  

  if( flag_BNB_CCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_CCpi0_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo);      
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo);      
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumu_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumu_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                         strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_FC_eventinfo, 5, "numu2numu");// hack
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumu_PC_eventinfo, 6, "numu2numu");// hack
    
  }// if( flag_BNB_CCpi0_from_overlaynumu )
    
  ///////////////////  

  if( flag_BNB_CCpi0_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_BNB_CCpi0_from_overlaynueNC"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynueNC_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynueNC_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynueNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynueNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynueNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                          strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynueNC_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynueNC_FC_eventinfo, 5, "nueNC");// hack
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynueNC_PC_eventinfo, 6, "nueNC");// hack
    
  }// if( flag_BNB_CCpi0_from_overlaynueNC )
  
  ///////////////////  

  if( flag_BNB_CCpi0_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_BNB_CCpi0_from_overlaynumuNC"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumuNC_PC_eventinfo);      
    }    
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumuNC_PC_eventinfo);      
    }    
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_overlaynumuNC_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumuNC_FC_eventinfo);
      str_treename = "tree_CCpi0_from_overlaynumuNC_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                           strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_overlaynumuNC_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumuNC_FC_eventinfo, 5, "numuNC");// hack
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_CCpi0_from_overlaynumuNC_PC_eventinfo, 6, "numuNC");// hack
    
  }// if( flag_BNB_CCpi0_from_overlaynumuNC )
  
  ///////////////////  

  if( flag_BNB_NCpi0_from_overlaynumu ) {
    cout<<endl<<"      ---> flag_BNB_NCpi0_from_overlaynumu"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo);
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumu";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo);
    }

    Set_oscillation_base_minus(&vector_BNB_NCpi0_from_overlaynumu_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynumu_eventinfo, 7, "numu2numu");// hack
    
  }// if( flag_BNB_NCpi0_from_overlaynumu )
      
  ///////////////////  

  if( flag_BNB_NCpi0_from_overlaynueNC ) {
    cout<<endl<<"      ---> flag_BNB_NCpi0_from_overlaynueNC xx"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynueNC_eventinfo);
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynueNC_eventinfo);
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynueNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynueNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynueNC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_BNB_NCpi0_from_overlaynueNC_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynueNC_eventinfo, 7, "nueNC");// hack
    
  }// if( flag_BNB_NCpi0_from_overlaynueNC )
         
  ///////////////////  

  if( flag_BNB_NCpi0_from_overlaynumuNC ) {
    cout<<endl<<"      ---> flag_BNB_NCpi0_from_overlaynumuNC"<<endl;    
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumuNC_eventinfo);
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumuNC_eventinfo);
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_nu_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_nu_overlay.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_overlaynumuNC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_overlaynumuNC_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_overlaynumuNC_eventinfo);
    }

    Set_oscillation_base_minus(&vector_BNB_NCpi0_from_overlaynumuNC_scaleFPOT, &vector_vector_BNB_NCpi0_from_overlaynumuNC_eventinfo, 7, "numuNC");// hack
    
  }// if( flag_BNB_NCpi0_from_overlaynumuNC )
    
  ///////////////////  

  if( flag_BNB_nueCC_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_BNB_nueCC_from_appnue"<<endl;
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_appnue_PC_eventinfo);      
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_appnue_PC_eventinfo);      
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_nueCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_nueCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_appnue_FC_eventinfo);
      str_treename = "tree_nueCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_nueCC_from_appnue_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_appnue_scaleFPOT, &vector_vector_BNB_nueCC_from_appnue_FC_eventinfo, 15, "numu2nue");// hack
    Set_oscillation_base_minus(&vector_BNB_nueCC_from_appnue_scaleFPOT, &vector_vector_BNB_nueCC_from_appnue_PC_eventinfo, 16, "numu2nue");// hack
    
  }// if( flag_BNB_nueCC_from_appnue )
    
  ///////////////////  

  if( flag_BNB_numuCC_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_BNB_numuCC_from_appnue"<<endl;
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_appnue_PC_eventinfo);      
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_appnue_PC_eventinfo);      
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_numuCC_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_numuCC_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_appnue_FC_eventinfo);
      str_treename = "tree_numuCC_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                     strfile_mc_e2e, str_treename, &vector_vector_BNB_numuCC_from_appnue_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_appnue_scaleFPOT, &vector_vector_BNB_numuCC_from_appnue_FC_eventinfo, 17, "numu2nue");// hack
    Set_oscillation_base_minus(&vector_BNB_numuCC_from_appnue_scaleFPOT, &vector_vector_BNB_numuCC_from_appnue_PC_eventinfo, 18, "numu2nue");// hack
    
  }// if( flag_BNB_numuCC_from_appnue )
    
  ///////////////////  

  if( flag_BNB_CCpi0_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_BNB_CCpi0_from_appnue"<<endl;
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_appnue_PC_eventinfo);      
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_appnue_PC_eventinfo);      
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_CCpi0_from_appnue_FC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_CCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_appnue_FC_eventinfo);
      str_treename = "tree_CCpi0_from_appnue_PC";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, NULL,                                    strfile_mc_e2e, str_treename, &vector_vector_BNB_CCpi0_from_appnue_PC_eventinfo);      
    }
    
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_appnue_scaleFPOT, &vector_vector_BNB_CCpi0_from_appnue_FC_eventinfo, 19, "numu2nue");// hack
    Set_oscillation_base_minus(&vector_BNB_CCpi0_from_appnue_scaleFPOT, &vector_vector_BNB_CCpi0_from_appnue_PC_eventinfo, 20, "numu2nue");// hack
    
  }// if( flag_BNB_CCpi0_from_appnue )
    
  ///////////////////  

  if( flag_BNB_NCpi0_from_appnue || 1 ) {
    cout<<endl<<"      ---> flag_BNB_NCpi0_from_appnue"<<endl;
    {// run1
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run1.root";
      TString strfile_dataPOT = str_dirbase + "run1_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run1_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_appnue_eventinfo);
    }
    {// run2
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run2.root";
      TString strfile_dataPOT = str_dirbase + "run2_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run2_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_appnue_eventinfo);
    }
    {// run3
      TString strfile_mcPOT   = str_dirbase + "checkout_prodgenie_bnb_numu2nue_overlay_run3.root";
      TString strfile_dataPOT = str_dirbase + "run3_data_bnb.root";
      TString strfile_mc_e2e  = str_dirbase + "roofile_obj_BNB_run3_appnue.root";
      TString str_treename = "";      
      str_treename = "tree_NCpi0_from_appnue";
      Set_oscillation_base_subfunc(strfile_mcPOT, strfile_dataPOT, &vector_BNB_NCpi0_from_appnue_scaleFPOT, strfile_mc_e2e, str_treename, &vector_vector_BNB_NCpi0_from_appnue_eventinfo);
    }
    
    Set_oscillation_base_minus(&vector_BNB_NCpi0_from_appnue_scaleFPOT, &vector_vector_BNB_NCpi0_from_appnue_eventinfo, 21, "numu2nue");// hack
    
  }// if( flag_BNB_NCpi0_from_appnue )

  
  cout<<endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////// ccc

void TOsc::Set_default_cv_cov(TString default_cv_file, TString default_dirtadd_file, TString default_mcstat_file, TString default_fluxXs_dir, TString default_detector_dir)
{
  TString roostr = "";
  
  cout<<endl<<" ---> Set_default_cv_cov"<<endl;
  
  cout<<endl;
  cout<<TString::Format("      ---> default_cv_file       %10s", default_cv_file.Data())<<endl;
  cout<<TString::Format("      ---> default_dirtadd_file  %10s", default_dirtadd_file.Data())<<endl;
  cout<<TString::Format("      ---> default_mcstat_file   %10s", default_mcstat_file.Data())<<endl;
  cout<<TString::Format("      ---> default_fluxXs_dir    %10s", default_fluxXs_dir.Data())<<endl;
  cout<<TString::Format("      ---> default_detector_dir  %10s", default_detector_dir.Data())<<endl;

  //////////////////////////////////////

  {
    TFile *roofile_default_cv_file = new TFile(default_cv_file, "read");

    ///////
    cout<<endl;
    cout<<"      ---> matrix_tosc_transform"<<endl;
    TMatrixD *temp_mat_collapse = (TMatrixD*)roofile_default_cv_file->Get("mat_collapse");  
    default_oldworld_rows = temp_mat_collapse->GetNrows();
    default_newworld_rows = temp_mat_collapse->GetNcols();
    matrix_tosc_transform.Clear(); matrix_tosc_transform.ResizeTo(default_oldworld_rows, default_newworld_rows);
    matrix_tosc_transform += (*temp_mat_collapse);
    delete temp_mat_collapse;
    
    /// hack, set LEE channel to 0
    // for(int idx=0; idx<default_oldworld_rows; idx++) {
    //   for(int jdx=0; jdx<default_newworld_rows; jdx++) {
    // 	if(idx>=26*4+11*3 && idx<26*4+11*3+26*2) matrix_tosc_transform(idx, jdx) = 0;
    //   }
    // }
    
    cout<<"      ---> default_oldworld_rows  "<<default_oldworld_rows<<endl;
    cout<<"      ---> default_newworld_rows  "<<default_newworld_rows<<endl;

    ///////

    matrix_tosc_default_newworld_meas.Clear(); matrix_tosc_default_newworld_meas.ResizeTo(1, default_newworld_rows);
    matrix_tosc_default_oldworld_pred.Clear(); matrix_tosc_default_oldworld_pred.ResizeTo(1, default_oldworld_rows);

    matrix_tosc_oscillation_base_oldworld_pred.Clear(); matrix_tosc_oscillation_base_oldworld_pred.ResizeTo(1, default_oldworld_rows);
    matrix_tosc_oscillation_oldworld_pred.Clear();      matrix_tosc_oscillation_oldworld_pred.ResizeTo(1, default_oldworld_rows);

    ///////
    
    matrix_tosc_default_oldworld_abs_syst_addi.Clear();   matrix_tosc_default_oldworld_abs_syst_addi.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  
    matrix_tosc_default_oldworld_rel_syst_flux.Clear();   matrix_tosc_default_oldworld_rel_syst_flux.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_tosc_default_oldworld_rel_syst_geant.Clear();  matrix_tosc_default_oldworld_rel_syst_geant.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_tosc_default_oldworld_rel_syst_Xs.Clear();     matrix_tosc_default_oldworld_rel_syst_Xs.ResizeTo(default_oldworld_rows, default_oldworld_rows);
    matrix_tosc_default_oldworld_rel_syst_det.Clear();    matrix_tosc_default_oldworld_rel_syst_det.ResizeTo(default_oldworld_rows, default_oldworld_rows);

    matrix_tosc_default_newworld_abs_syst_mcstat.Clear(); matrix_tosc_default_newworld_abs_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);
    
    //////
    
    matrix_tosc_eff_newworld_abs_syst_addi.Clear();   matrix_tosc_eff_newworld_abs_syst_addi.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_tosc_eff_newworld_abs_syst_mcstat.Clear(); matrix_tosc_eff_newworld_abs_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_tosc_eff_newworld_abs_syst_flux.Clear();   matrix_tosc_eff_newworld_abs_syst_flux.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_tosc_eff_newworld_abs_syst_geant.Clear();  matrix_tosc_eff_newworld_abs_syst_geant.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_tosc_eff_newworld_abs_syst_Xs.Clear();     matrix_tosc_eff_newworld_abs_syst_Xs.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_tosc_eff_newworld_abs_syst_det.Clear();    matrix_tosc_eff_newworld_abs_syst_det.ResizeTo(default_newworld_rows, default_newworld_rows);
    matrix_tosc_eff_newworld_abs_syst_total.Clear();  matrix_tosc_eff_newworld_abs_syst_total.ResizeTo(default_newworld_rows, default_newworld_rows);
	  
    matrix_tosc_eff_newworld_meas.Clear();  matrix_tosc_eff_newworld_meas.ResizeTo(1, default_newworld_rows);
    matrix_tosc_eff_newworld_pred.Clear();  matrix_tosc_eff_newworld_pred.ResizeTo(1, default_newworld_rows);
    matrix_tosc_eff_newworld_noosc.Clear(); matrix_tosc_eff_newworld_noosc.ResizeTo(1, default_newworld_rows);

    matrix_tosc_fitdata_newworld.Clear(); matrix_tosc_fitdata_newworld.ResizeTo(1, default_newworld_rows);
    
    ///////
    cout<<endl;
    cout<<"      ---> measurement"<<endl;
    for(int idx=1; idx<=10000; idx++) {
      roostr = TString::Format("hdata_obsch_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      if(h1f_temp==NULL) { delete h1f_temp; break; }
      cout<<TString::Format("      %3d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

      int bins = h1f_temp->GetNbinsX(); double xlow = h1f_temp->GetXaxis()->GetBinLowEdge(1); double xhgh = h1f_temp->GetXaxis()->GetBinUpEdge(bins);
      map_tosc_default_h1d_meas_bins[idx] = bins; map_tosc_default_h1d_meas_xlow[idx] = xlow; map_tosc_default_h1d_meas_xhgh[idx] = xhgh;                 
      for(int ibin=1; ibin<=bins+1; ibin++) vector_tosc_default_newworld_meas.push_back( h1f_temp->GetBinContent(ibin) );      
      delete h1f_temp;    
    }
    
    if( default_newworld_rows!=(int)(vector_tosc_default_newworld_meas.size()) ) { cerr<<" ---> ERROR: default_newworld_rows!=vector_tosc_default_newworld_meas"<<endl; exit(1); }
    for(int idx=0; idx<(int)(vector_tosc_default_newworld_meas.size()); idx++)  matrix_tosc_default_newworld_meas(0, idx) = vector_tosc_default_newworld_meas.at(idx);
        
    ///////
    cout<<endl;
    cout<<"      ---> prediction"<<endl;
    for(int idx=1; idx<=10000; idx++) {
      roostr = TString::Format("histo_%d", idx);
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      if(h1f_temp==NULL) { delete h1f_temp; break; }
      cout<<TString::Format("      %3d,  bins %2d, %20s", idx, h1f_temp->GetNbinsX()+1, h1f_temp->GetTitle())<<endl;

      int bins = h1f_temp->GetNbinsX(); double xlow = h1f_temp->GetXaxis()->GetBinLowEdge(1); double xhgh = h1f_temp->GetXaxis()->GetBinUpEdge(bins);
      map_tosc_default_h1d_pred_bins[idx] = bins; map_tosc_default_h1d_pred_xlow[idx] = xlow; map_tosc_default_h1d_pred_xhgh[idx] = xhgh;                 
      for(int ibin=1; ibin<=bins+1; ibin++) vector_tosc_default_oldworld_pred.push_back( h1f_temp->GetBinContent(ibin) );      
      delete h1f_temp;    
    }
    
    if( default_oldworld_rows!=(int)(vector_tosc_default_oldworld_pred.size()) ) { cerr<<" ---> ERROR: default_oldworld_rows!=vector_tosc_default_oldworld_pred"<<endl; exit(1); }
    for(int idx=0; idx<(int)(vector_tosc_default_oldworld_pred.size()); idx++)  matrix_tosc_default_oldworld_pred(0, idx) = vector_tosc_default_oldworld_pred.at(idx);
    
    ///////
    
    delete roofile_default_cv_file;
  }

  
  {
    for(auto it=map_tosc_default_h1d_meas_bins.begin(); it!=map_tosc_default_h1d_meas_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("default_h1d_meas_%d", idx);
      map_tosc_default_h1d_meas[idx] = new TH1D(roostr, roostr, it->second, map_tosc_default_h1d_meas_xlow[idx], map_tosc_default_h1d_meas_xhgh[idx]);
    }
    
    for(auto it=map_tosc_default_h1d_pred_bins.begin(); it!=map_tosc_default_h1d_pred_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("default_h1d_pred_%d", idx);
      map_tosc_default_h1d_pred[idx] = new TH1D(roostr, roostr, it->second, map_tosc_default_h1d_pred_xlow[idx], map_tosc_default_h1d_pred_xhgh[idx]);      
    }

    TFile *roofile_default_cv_file = new TFile(default_cv_file, "read");

    for(auto it=map_tosc_default_h1d_meas_bins.begin(); it!=map_tosc_default_h1d_meas_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("hdata_obsch_%d", idx);
      int bins = it->second;
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      for(int ibin=1; ibin<=bins+1; ibin++) map_tosc_default_h1d_meas[idx]->SetBinContent(ibin, h1f_temp->GetBinContent(ibin));
      delete h1f_temp;
    }
    
    for(auto it=map_tosc_default_h1d_pred_bins.begin(); it!=map_tosc_default_h1d_pred_bins.end(); it++) {
      int idx = it->first; roostr = TString::Format("histo_%d", idx);
      int bins = it->second;
      TH1F *h1f_temp = (TH1F*)roofile_default_cv_file->Get(roostr);
      for(int ibin=1; ibin<=bins+1; ibin++) map_tosc_default_h1d_pred[idx]->SetBinContent(ibin, h1f_temp->GetBinContent(ibin));
      delete h1f_temp;
    }
    
    delete roofile_default_cv_file;
  }
  
  //////////////////////////////////////

  if( flag_syst_dirt ) {
    cout<<endl;
    cout<<"      ---> Dirt: additional uncertainty, Yes"<<endl;
    
    TFile *roofile_default_dirtadd_file = new TFile(default_dirtadd_file, "read");
    TMatrixD *matrix_temp = (TMatrixD*)roofile_default_dirtadd_file->Get("cov_mat_add");
    matrix_tosc_default_oldworld_abs_syst_addi = (*matrix_temp);
    delete matrix_temp;
    delete roofile_default_dirtadd_file;
  }
  else {
    cout<<endl;
    cout<<"      ---> Dirt: additional uncertainty, No"<<endl;
  }
  
  //////////////////////////////////////

  if( flag_syst_mcstat ) {
    cout<<endl;
    cout<<"      ---> MCstat, Yes"<<endl;
    
    ifstream file_mcstat_aa(default_mcstat_file, ios::in);
    if(!file_mcstat_aa) { cerr<<" Error: No file_mcstat_aa"<<endl; exit(1); }
    int count_aa = 0;
    string str_count_aa;    
    ifstream file_check_aa(default_mcstat_file);
    while( getline(file_check_aa, str_count_aa) ) count_aa++;
    if( count_aa-1 != default_newworld_rows ) { cerr<<" Error: mcstat != default_newworld_rows"<<endl; exit(1);  }

    ifstream file_mcstat_bb(default_mcstat_file, ios::in);
    double Lee = 1; double run = 1;
    file_mcstat_bb >> Lee >> run;
    int gch = 0; int lbin = 0; double val_pred = 0; double mc_stat = 0; double nn_stat = 0;
    for(int idx=1; idx<=default_newworld_rows; idx++) {
      file_mcstat_bb >> gch >> lbin >> val_pred >> mc_stat >> nn_stat;
      matrix_tosc_default_newworld_abs_syst_mcstat(idx-1, idx-1) = mc_stat;
    }

  }
  else {
    cout<<endl;
    cout<<"      ---> MCstat, No"<<endl;
  }
  
  ////////////////////////////////////// flux, geant, Xs

  if( flag_syst_flux ) {
    cout<<endl;
    cout<<"      ---> flux, Yes"<<endl;
    
    for(int idx=1; idx<=13; idx++) {
      TFile *roofile_syst_flux = new TFile(default_fluxXs_dir+TString::Format("cov_%d.root", idx), "read");
      TMatrixD *matrix_temp_flux = (TMatrixD*)roofile_syst_flux->Get(TString::Format("frac_cov_xf_mat_%d", idx) );
      matrix_tosc_default_oldworld_rel_syst_flux += (*matrix_temp_flux);
      delete matrix_temp_flux;
      delete roofile_syst_flux;
    }
	
  }// if( flag_syst_flux )
  else {
    cout<<endl;
    cout<<"      ---> flux, No"<<endl;
  }
  
  if( flag_syst_geant ) {
    cout<<endl;
    cout<<"      ---> geant, Yes"<<endl;

    for(int idx=14; idx<=16; idx++) {
      TFile *roofile_syst_geant = new TFile(default_fluxXs_dir+TString::Format("cov_%d.root", idx), "read");
      TMatrixD *matrix_temp_geant = (TMatrixD*)roofile_syst_geant->Get(TString::Format("frac_cov_xf_mat_%d", idx) );
      matrix_tosc_default_oldworld_rel_syst_geant += (*matrix_temp_geant);
      delete matrix_temp_geant;
      delete roofile_syst_geant;
    }   
  }// if( flag_syst_geant )
  else {
    cout<<endl;
    cout<<"      ---> geant, No"<<endl;
  }
  
  if( flag_syst_Xs ) {
    cout<<endl;
    cout<<"      ---> Xs, Yes"<<endl;
    
    TFile *roofile_syst_Xs = new TFile(default_fluxXs_dir+"cov_17.root", "read");
    TMatrixD *matrix_temp_Xs = (TMatrixD*)roofile_syst_Xs->Get("frac_cov_xf_mat_17");
    matrix_tosc_default_oldworld_rel_syst_Xs = (*matrix_temp_Xs);
    delete matrix_temp_Xs;
    delete roofile_syst_Xs; 
  }// if( flag_syst_Xs )
  else {
    cout<<endl;
    cout<<"      ---> Xs, No"<<endl;
  }
  
  ////////////////////////////////////// detector

  if( flag_syst_det ) {
    cout<<endl;
    cout<<"      ---> detector, Yes"<<endl;
 
    map<int, TString>map_detectorfile_str;    
    map_detectorfile_str[1] = "cov_LYDown.root";
    map_detectorfile_str[2] = "cov_LYRayleigh.root";
    map_detectorfile_str[3] = "cov_Recomb2.root";
    map_detectorfile_str[4] = "cov_SCE.root";
    //map_detectorfile_str[5] = "cov_WMdEdx.root";
    map_detectorfile_str[6] = "cov_WMThetaXZ.root";
    map_detectorfile_str[7] = "cov_WMThetaYZ.root";
    map_detectorfile_str[8] = "cov_WMX.root";
    map_detectorfile_str[9] = "cov_WMYZ.root";
    map_detectorfile_str[10]= "cov_LYatt.root";
    
    for(auto it_map=map_detectorfile_str.begin(); it_map!=map_detectorfile_str.end(); it_map++) {
      int idx = it_map->first;
      TFile *roofile_det = new TFile(default_detector_dir + map_detectorfile_str[idx], "read");
      TMatrixD *matrix_temp_det = (TMatrixD*)roofile_det->Get(TString::Format("frac_cov_det_mat_%d", idx) );
      matrix_tosc_default_oldworld_rel_syst_det += (*matrix_temp_det);
      delete matrix_temp_det;
      delete roofile_det;
    }
  }// if( flag_syst_det )
  else {
    cout<<endl;
    cout<<"      ---> detector, No"<<endl;
  }
  
  /////////////////////////// TOsc::Set_apply_POT()
  
  matrix_tosc_temp_oldworld_abs_syst_addi.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_tosc_temp_newworld_abs_syst_mcstat.ResizeTo(default_newworld_rows, default_newworld_rows);

  matrix_tosc_temp_oldworld_rel_syst.ResizeTo(default_oldworld_rows, default_oldworld_rows);
  matrix_tosc_temp_oldworld_abs_syst.ResizeTo(default_oldworld_rows, default_oldworld_rows);

  user_matrix_tosc_temp_oscillation_oldworld_pred.ResizeTo(1, default_oldworld_rows);
  user_matrix_tosc_temp_oscillation_oldworld_pred_T.ResizeTo(default_oldworld_rows, 1);
  user_matrix_tosc_temp_cv_ij_oldworld.ResizeTo(default_oldworld_rows, default_oldworld_rows);

  matrix_tosc_temp_oldworld_rel_syst.Zero();
  matrix_tosc_temp_oldworld_rel_syst += matrix_tosc_default_oldworld_rel_syst_flux;
  matrix_tosc_temp_oldworld_rel_syst += matrix_tosc_default_oldworld_rel_syst_geant;
  matrix_tosc_temp_oldworld_rel_syst += matrix_tosc_default_oldworld_rel_syst_Xs;
  matrix_tosc_temp_oldworld_rel_syst += matrix_tosc_default_oldworld_rel_syst_det;

  ///////

  int unit_block = default_oldworld_rows/3/2;// overlay, EXT, OSC ---> BNB and NuMI

  matrix_tosc_AA.ResizeTo(unit_block, unit_block);
  matrix_tosc_BB.ResizeTo(unit_block, unit_block);
  matrix_tosc_CC.ResizeTo(unit_block, unit_block);
  
  matrix_tosc_AB.ResizeTo(unit_block, unit_block);
  matrix_tosc_AC.ResizeTo(unit_block, unit_block);
  matrix_tosc_BC.ResizeTo(unit_block, unit_block);
  
  matrix_tosc_BA.ResizeTo(unit_block, unit_block);
  matrix_tosc_CA.ResizeTo(unit_block, unit_block);
  matrix_tosc_CB.ResizeTo(unit_block, unit_block);
          
  matrix_tosc_BNB2BNB.ResizeTo(unit_block, unit_block);
  matrix_tosc_NuMI2NuMI.ResizeTo(unit_block, unit_block);
  matrix_tosc_BNB2NuMI.ResizeTo(unit_block, unit_block);
  matrix_tosc_NuMI2BNB.ResizeTo(unit_block, unit_block);

  matrix_tosc_temp_total_in_POT.ResizeTo(unit_block*2, unit_block*2);
}

