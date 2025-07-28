#include "Minuit/FCNBase.h"

data_mc_overlay() {
    const TString kDataFile("output_data_0.root");
    const TString kMCFile("output_mc_0.root");

    const TString kHistName("htxz_13");

    TFile* fd = TFile::Open(kDataFile);
    TFile* fm = TFile::Open(kMCFile);

    TH1D* hd = static_cast<TH1D*>(fd->Get(kHistName));
    TH1D* hm = static_cast<TH1D*>(fm->Get(kHistName));

    hsample->Reset();

    const UInt_t kNthrows = 10000;
    Double_t throws[kNthrows]
    for (UInt_t i = 0; i < kNthrows; i++) {
        throws[i] = hd->GetRandom();
    }

    auto fitfunc = [=](Double_t* x, Double_t* parm) {
        TH1D* hsample = static_cast<TH1D*>(hd->Clone());
        hsample->Reset();
        for (UInt_t i = 0; i < kNthrows; i++) {
            hsample->Fill(throws[i] * parm[0]);
        }

        delete hsample;
    };

    THStack* hs = new THStack();
    hd->Scale(1 / hd->Integral());
    hm->Scale(1 / hm->Integral());
    hsample->Scale(1 / hsample->Integral());
    hs->Add(hm, "histE");
    hs->Add(hd, "EP");
    hs->Add(hsample, "hist");

    hd->SetLineColor(kBlack);
    hm->SetLineColor(kRed);
    hs->Draw("nostack");
}


