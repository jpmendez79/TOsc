// Draw the Brazil-band sensitivity plot from run_analysis.cxx output:
// median (dashed), 1sigma/2sigma containment bands, and one example
// toy-universe curve (gr_0001) overlaid.
//
// Usage:
//   root -l plot_brazil_band.C
//   root -l 'plot_brazil_band.C("cls_map_BNB_vanilla_disp_60x60.root")'
// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}

void plot_brazil_band(TString infile = "test/cls_map_BNB_vanilla_disp_60x60.root")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    cout << "*** could not open " << infile << endl;
    return;
  }

  bool has_g2 = false;
  double g2_value = 0;
  TPRegexp g2_pattern("g2_([0-9]+\\.[0-9]{2})");
  TObjArray *g2_matches = g2_pattern.MatchS(infile);
  if (g2_matches->GetLast() == 1) {
    has_g2 = true;
    g2_value = ((TObjString*)g2_matches->At(1))->GetString().Atof();
    cout << "g2 = " << g2_value << endl;
  }
  delete g2_matches;

  TGraphAsymmErrors *sigma1 = (TGraphAsymmErrors*)f->Get("sigma1");
  TGraphAsymmErrors *sigma2 = (TGraphAsymmErrors*)f->Get("sigma2");
  TGraph *gr_example = (TGraph*)f->Get("gr_0001");

  if (!sigma1 || !sigma2 || !gr_example) {
    cout << "*** missing sigma1/sigma2/gr_0001 in " << infile << endl;
    return;
  }

  int npoints = sigma1->GetN();

  // Median dashed line, read off the sigma1 central points
  TGraph *gh_median = new TGraph(npoints);
  for (int i = 0; i < npoints; i++) {
    double x, y;
    sigma1->GetPoint(i, x, y);
    gh_median->SetPoint(i, x, y);
  }

  // sigma1/sigma2 hold errors *in x* as a function of y (dm2), i.e. a
  // horizontal band. ROOT's "3"/"4" TGraph fill options instead fill a
  // vertical band (y-errors vs x), so build the bands explicitly as closed
  // polygons: low-x boundary going up in y, then high-x boundary back down.
  auto make_band = [&](TGraphAsymmErrors *sigma) {
    TGraph *band = new TGraph(2 * npoints);
    for (int i = 0; i < npoints; i++) {
      double x, y;
      sigma->GetPoint(i, x, y);
      band->SetPoint(i, x - sigma->GetErrorXlow(i), y);
    }
    for (int i = 0; i < npoints; i++) {
      double x, y;
      sigma->GetPoint(npoints - 1 - i, x, y);
      band->SetPoint(npoints + i, x + sigma->GetErrorXhigh(npoints - 1 - i), y);
    }
    return band;
  };

  TGraph *band1 = make_band(sigma1);
  TGraph *band2 = make_band(sigma2);

  ///////

  TString roostr = "canv_brazil_band";
  TCanvas *canv_brazil_band = new TCanvas(roostr, roostr, 800, 700);
  canv_brazil_band->SetLeftMargin(0.15);
  canv_brazil_band->SetRightMargin(0.1);
  canv_brazil_band->SetTopMargin(0.1);
  canv_brazil_band->SetBottomMargin(0.15);
  canv_brazil_band->SetLogx();
  canv_brazil_band->SetLogy();

  band2->SetFillColor(kYellow);
  band2->SetFillStyle(1001);
  band2->SetLineColor(kYellow);

  band1->SetFillColor(kGreen);
  band1->SetFillStyle(1001);
  band1->SetLineColor(kGreen);

  gh_median->SetLineColor(kBlack);
  gh_median->SetLineStyle(2);   // dashed
  gh_median->SetLineWidth(2);

  gr_example->SetLineColor(kBlue);
  gr_example->SetLineStyle(1);
  gr_example->SetLineWidth(2);

  band2->Draw("AF");
  band2->GetXaxis()->SetTitle("sin^{2}(2#theta_{#mu#mu})");
  band2->GetYaxis()->SetTitle("#Deltam^{2} [eV^{2}]");
  band2->GetXaxis()->CenterTitle(1);
  band2->GetYaxis()->CenterTitle(1);
  band2->GetXaxis()->SetTitleSize(0.05);
  band2->GetYaxis()->SetTitleSize(0.05);
  band2->GetXaxis()->SetLabelSize(0.04);
  band2->GetYaxis()->SetLabelSize(0.04);
  band2->GetYaxis()->SetRangeUser(0.09, 100);
  band2->GetXaxis()->SetRangeUser(0.01, 1);


  band1->Draw("F same");
  gh_median->Draw("L same");
  gr_example->Draw("L same");

  TLegend *lg = new TLegend(0.5, 0.65, 0.87, 0.87);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetTextSize(0.03);
  lg->AddEntry(gh_median, "Median sensitivity", "l");
  // lg->AddEntry(gr_example, "95% CLs Data Exclusion", "l");
  lg->Draw();

  if (has_g2) {
    TLatex *label_g2 = new TLatex(0.5, 0.89, Form("g^{2} = %.2f#pi", g2_value));
    label_g2->SetNDC();
    label_g2->SetTextSize(0.04);
    label_g2->Draw();
  }


  // Force the raster output to the requested pixel size, independent of
  // whatever size the on-screen window ended up at (window managers can
  // resize/decorate it, which otherwise crops the saved PNG).
  canv_brazil_band->SetCanvasSize(800, 700);
  canv_brazil_band->SaveAs("canv_brazil_band.png");
}
