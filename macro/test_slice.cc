TH1D* slice_from_2d_hist(const TH2* h, size_t idx, Int_t axis=0, Double_t scale=1.0) {
    // clone the projection, otherwise multiple calls to projection return the same object
    TH1D* hslice = axis == 0 ? (TH1D*)h->ProjectionY()->Clone() : (TH1D*)h->ProjectionX()->Clone();
    hslice->Reset();

    Int_t nbins = axis == 0 ? h->GetNbinsY() : h->GetNbinsX();
    for (Int_t i = 1; i < nbins; i++) {
        Int_t ix = idx;
        Int_t iy = i;
        if (axis != 0) {
            ix = i;
            iy = idx;
        }
        // Int_t the_bin = hslice->FindFixBin(scale * hslice->GetBinCenter(i));
        // hslice->SetBinContent(the_bin, h->GetBinContent(ix, iy));
        // hslice->SetBinError(the_bin, h->GetBinError(ix, iy));
        hslice->SetBinContent(i, h->GetBinContent(ix, iy));
        hslice->SetBinError(i, h->GetBinError(ix, iy));
    }

    return hslice;
}


TH2D* hist_from_file(const TString& fname, const TString& hname) {
    TFile* f = TFile::Open(fname, "read");
    TH2D* result = (TH2D*)f->Get(hname);
    result->SetDirectory(0);
    f->Close();
    return result;
}


void test_slice(const TString& f1, const TString& f2, const TString& var, size_t idx) {

    const TString hname(Form("h2d_%s_vs_width", var.Data()));

    TH2D* h1 = hist_from_file(f1, hname);
    TH2D* h2 = hist_from_file(f2, hname);

    TH1D* hslice1 = slice_from_2d_hist(h1, idx, 0, 1.0);
    TH1D* hslice2 = slice_from_2d_hist(h2, idx, 0, 1.1);

    hslice1->Scale(1.0 / hslice1->Integral());
    hslice2->Scale(1.0 / hslice2->Integral());

    hslice1->SetLineColor(kBlack);
    hslice2->SetLineColor(kRed);
    THStack* hs = new THStack();
    hs->Add(hslice1);
    hs->Add(hslice2);

    hs->Draw("nostack hist e");
}
