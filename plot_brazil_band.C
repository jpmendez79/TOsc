// Draw the Brazil-band sensitivity plot from run_analysis.cxx output:
// median (dashed), 1sigma/2sigma containment bands, and one example
// toy-universe curve (gr_0001) overlaid.
//
// Usage:
//   root -l plot_brazil_band.C
//   root -l 'plot_brazil_band.C("cls_map_BNB_vanilla_disp_60x60.root")'

void plot_brazil_band(TString infile = "test/cls_map_BNB_vanilla_disp_60x60.root")
{
  gStyle->SetOptStat(0);

  TFile *f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    cout << "*** could not open " << infile << endl;
    return;
  }

  TGraphAsymmErrors *sigma1 = (TGraphAsymmErrors*)f->Get("sigma1");
  TGraphAsymmErrors *sigma2 = (TGraphAsymmErrors*)f->Get("sigma2");
  TGraph *gr_example = (TGraph*)f->Get("gr_0000");

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
  band2->GetXaxis()->SetTitle("sin^{2}2#theta");
  band2->GetYaxis()->SetTitle("#Deltam^{2} [eV^{2}]");
  band2->GetXaxis()->CenterTitle(1);
  band2->GetYaxis()->CenterTitle(1);
  band2->GetXaxis()->SetTitleSize(0.05);
  band2->GetYaxis()->SetTitleSize(0.05);
  band2->GetXaxis()->SetLabelSize(0.04);
  band2->GetYaxis()->SetLabelSize(0.04);


  band1->Draw("F same");
  gh_median->Draw("L same");
  gr_example->Draw("L same");

  TLegend *lg = new TLegend(0.5, 0.65, 0.87, 0.87);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetTextSize(0.03);
  lg->AddEntry(gh_median, "Median expected", "l");
  lg->AddEntry(band1, "1#sigma band", "f");
  lg->AddEntry(band2, "2#sigma band", "f");
  lg->AddEntry(gr_example, "gr_0001 (example universe)", "l");
  lg->Draw();

  canv_brazil_band->SaveAs("canv_brazil_band.png");
}
