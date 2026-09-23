// Overlay the gr_0000 TGraph from each regenerated MicroBooNE CLs-map file
// (one per invisible-decay g^2 value) plus the digitized IceCube exclusion
// contour, to visually compare how the exclusion curve shifts with g^2.
//
// Usage:
//   root -l -b -q plot_direct_graphs.C

void plot_direct_graphs()
{
  gStyle->SetOptStat(0);

  struct Entry { TString file; TString label; bool is_csv; };
  std::vector<Entry> entries = {
    {"slice-analysis/cls_map_BNB_inv_decay_g2_0.00_numu_disp_60x60.root", "MicroBooNE 95% CLs g^{2}=0#pi", false},
    {"slice-analysis/cls_map_BNB_inv_decay_g2_1.00_numu_disp_60x60.root", "MicroBooNE 95% CLs g^{2}=1#pi", false},
    {"slice-analysis/cls_map_BNB_inv_decay_g2_2.00_numu_disp_60x60.root", "MicroBooNE 95% CLs g^{2}=2#pi", false},
    {"slice-analysis/cls_map_BNB_inv_decay_g2_2.50_numu_disp_60x60.root", "MicroBooNE 95% CLs g^{2}=2.50#pi", false},
    {"slice-analysis/cls_map_BNB_inv_decay_g2_3.00_numu_disp_60x60.root", "MicroBooNE 95% CLs g^{2}=3#pi", false},
    {"slice-analysis/cls_map_BNB_inv_decay_g2_4.00_numu_disp_60x60.root", "MicroBooNE 95% CLs g^{2}=4#pi", false},
    {"slice-analysis/digitized-icecube-g2_2_50-wilkes99.csv", "IceCube Wilkes 99% g^{2}=2.5#pi", true},
    {"slice-analysis/digitized-icecube-g2_2_50-wilkes95.csv", "IceCube Wilkes 95% g^{2}=2.5#pi", true},
  };

  std::vector<int> colors = {kBlue, kRed, kGreen+2, kMagenta, kOrange+1, kCyan+2, kBlack, kViolet+1};

  std::vector<TGraph*> graphs;
  std::vector<TString> labels;

  for (size_t i = 0; i < entries.size(); i++) {
    TGraph *gr_clone = nullptr;

    if (entries[i].is_csv) {
      TGraph *gr = new TGraph(entries[i].file, "%lg, %lg");
      if (gr->GetN() == 0) {
        cout << "*** could not read " << entries[i].file << endl;
        continue;
      }
      gr_clone = gr;
    } else {
      TFile *f = TFile::Open(entries[i].file, "READ");
      if (!f || f->IsZombie()) {
        cout << "*** could not open " << entries[i].file << endl;
        continue;
      }

      TGraph *gr = (TGraph*)f->Get("gr_0000");
      if (!gr) {
        cout << "*** missing gr_0000 in " << entries[i].file << endl;
        f->Close();
        continue;
      }

      gr_clone = (TGraph*)gr->Clone(Form("gr_direct_%zu", i));
      f->Close();
    }

    gr_clone->SetLineColor(colors[i % colors.size()]);
    gr_clone->SetLineWidth(2);
    gr_clone->SetLineStyle(entries[i].is_csv ? 2 : 1);

    graphs.push_back(gr_clone);
    labels.push_back(entries[i].label);
  }

  if (graphs.empty()) {
    cout << "*** no graphs loaded, nothing to plot" << endl;
    return;
  }

  TString roostr = "canv_direct_graphs";
  TCanvas *canv_direct_graphs = new TCanvas(roostr, roostr, 1000, 700);
  canv_direct_graphs->SetLeftMargin(0.12);
  canv_direct_graphs->SetRightMargin(0.37);
  canv_direct_graphs->SetTopMargin(0.1);
  canv_direct_graphs->SetBottomMargin(0.15);
  canv_direct_graphs->SetLogx();
  canv_direct_graphs->SetLogy();

  graphs[0]->GetXaxis()->SetLimits(0.005, 1.2);
  graphs[0]->SetMinimum(0.008);
  graphs[0]->SetMaximum(120);

  graphs[0]->Draw("AL");
  graphs[0]->GetXaxis()->SetTitle("sin^{2}2#theta");
  graphs[0]->GetYaxis()->SetTitle("#Deltam^{2} [eV^{2}]");
  graphs[0]->GetXaxis()->CenterTitle(1);
  graphs[0]->GetYaxis()->CenterTitle(1);
  graphs[0]->GetXaxis()->SetTitleSize(0.05);
  graphs[0]->GetYaxis()->SetTitleSize(0.05);
  graphs[0]->GetXaxis()->SetLabelSize(0.04);
  graphs[0]->GetYaxis()->SetLabelSize(0.04);

  for (size_t i = 1; i < graphs.size(); i++) {
    graphs[i]->Draw("L same");
  }

  TLegend *lg = new TLegend(0.63, 0.28, 0.99, 0.67);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->SetTextSize(0.030);
  for (size_t i = 0; i < graphs.size(); i++) {
    lg->AddEntry(graphs[i], labels[i], "l");
  }
  lg->Draw();

  canv_direct_graphs->SaveAs("canv_direct_graphs.png");
}
