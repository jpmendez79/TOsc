// Compare the MicroBooNE gr_0000 CLs exclusion curve (g^2=2.5pi) against
// the digitized IceCube (Wilks 99%) exclusion contour.
//
// Usage:
//   root -l -q plot_compare_icecube.C

void plot_compare_icecube()
{
  gStyle->SetOptStat(0);

  TString rootfile = "slice-analysis/cls_map_BNB_inv_decay_g2_2_50_numu_disp_60x60-direct.root";
  TString csvfile = "slice-analysis/digitized-icecube-g2_2_50-wilkes99.csv";

  TFile *f = TFile::Open(rootfile, "READ");
  if (!f || f->IsZombie()) {
    cout << "*** could not open " << rootfile << endl;
    return;
  }

  TGraph *gr_ub = (TGraph*)f->Get("gr_0000");
  if (!gr_ub) {
    cout << "*** missing gr_0000 in " << rootfile << endl;
    f->Close();
    return;
  }
  gr_ub = (TGraph*)gr_ub->Clone("gr_ub");
  f->Close();

  TGraph *gr_icecube = new TGraph(csvfile, "%lg, %lg");
  if (gr_icecube->GetN() == 0) {
    cout << "*** could not read " << csvfile << endl;
    return;
  }

  gr_ub->SetLineColor(kBlue);
  gr_ub->SetLineWidth(2);
  gr_ub->SetLineStyle(1);

  gr_icecube->SetLineColor(kRed);
  gr_icecube->SetLineWidth(2);
  gr_icecube->SetLineStyle(2);

  TString roostr = "canv_compare_icecube";
  TCanvas *canv_compare_icecube = new TCanvas(roostr, roostr, 800, 700);
  canv_compare_icecube->SetLeftMargin(0.15);
  canv_compare_icecube->SetRightMargin(0.1);
  canv_compare_icecube->SetTopMargin(0.1);
  canv_compare_icecube->SetBottomMargin(0.15);
  canv_compare_icecube->SetLogx();
  canv_compare_icecube->SetLogy();

  gr_ub->GetXaxis()->SetLimits(0.005, 1.2);
  gr_ub->SetMinimum(0.008);
  gr_ub->SetMaximum(120);

  gr_ub->Draw("AL");
  gr_ub->GetXaxis()->SetTitle("sin^{2}2#theta");
  gr_ub->GetYaxis()->SetTitle("#Deltam^{2} [eV^{2}]");
  gr_ub->GetXaxis()->CenterTitle(1);
  gr_ub->GetYaxis()->CenterTitle(1);
  gr_ub->GetXaxis()->SetTitleSize(0.05);
  gr_ub->GetYaxis()->SetTitleSize(0.05);
  gr_ub->GetXaxis()->SetLabelSize(0.04);
  gr_ub->GetYaxis()->SetLabelSize(0.04);

  gr_icecube->Draw("L same");

  TLegend *lg = new TLegend(0.5, 0.68, 0.87, 0.87);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetTextSize(0.03);
  lg->AddEntry(gr_ub, "MicroBooNE CLs 95%", "l");
  lg->AddEntry(gr_icecube, "IceCube Wilkes 99%", "l");
  lg->Draw();

  TLatex *label = new TLatex();
  label->SetNDC();
  label->SetTextSize(0.04);
  label->DrawLatex(0.2, 0.85, "g^{2}=2.5#pi");

  canv_compare_icecube->SaveAs("canv_compare_icecube.png");
}
