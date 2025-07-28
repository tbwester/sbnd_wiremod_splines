#include <iostream>
#include <vector>

#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "THStack.h"
#include "TCanvas.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"


/*
 * Class to sample data histogram (hdata) and return chisquare of samples to MC
 * historgam (hmc)
 */
class BootstrapFCN : public ROOT::Minuit2::FCNBase {
public:
    BootstrapFCN(UInt_t _nsamples, const TH1* _hmc, const TH1* _hdata) :
        nsamples(_nsamples), hmc(static_cast<TH1*>(_hmc->Clone())),
        hdata(static_cast<TH1*>(_hdata->Clone()))
{
    samples.resize(nsamples);
    for (UInt_t i = 0; i < nsamples; i++) {
        samples.at(i) = hdata->GetRandom();
    }

    hmc->Scale((float)nsamples / hmc->Integral());

    best_chi2 = 1.0e7;
}

    virtual double operator()(const std::vector<double>&) const override;
    virtual double Up() const override { return 0.5; };

    void SetNSamples(UInt_t n) { nsamples = n; }
    UInt_t GetNSamples() const { return nsamples; }
    TH1D* GetHSample(Double_t scale, Double_t shift) const;

protected:
    TH1* hmc;
    TH1* hdata;

private:
    UInt_t nsamples;
    std::vector<Double_t> samples;

    // our own internal tracking
    Double_t best_chi2;
    TH1D* best_sample;
};


double BootstrapFCN::operator()(const std::vector<double>& par) const {
    TH1D* hsample = GetHSample(par[0], par[1]);

    double chi2 = 0.;
    for (Int_t i = 1; i <= hdata->GetNbinsX(); i++) {
        // if (hdata->GetBinCenter(i) < 2.0 || hdata->GetBinCenter(i) > 3.5) continue;
        Double_t mc = hmc->GetBinContent(i);
        Double_t dt = hsample->GetBinContent(i);
        if (mc <= 0.) continue;
        if (dt <= 0.) {
            chi2 += mc - dt;
        }
        else {
            chi2 += mc - dt + dt * TMath::Log(dt / mc);
        }
    }

    // chi2 = hsample->Chi2Test(hmc, "UW CHI2");
    // if (chi2 != chi2) chi2 = 1.0e99;
    chi2 *= 2;
    std::cout << "For scale=" << par[0] << " and shift=" << par[1] << " chisq=" << chi2 << "\n";
    delete hsample;
    return chi2;
}


TH1D* BootstrapFCN::GetHSample(Double_t scale, Double_t shift) const {
    TH1D* hsample = static_cast<TH1D*>(hdata->Clone());
    hsample->Reset();
    Double_t upper_edge = hdata->GetBinLowEdge(hdata->GetNbinsX()) 
        + hdata->GetBinWidth(hdata->GetNbinsX());

    for (UInt_t i = 0; i < nsamples; i++) {
        hsample->Fill(TMath::Max(0., TMath::Min(upper_edge - 1.0e-6,
            samples.at(i) * scale + shift
        )));
    }

    return hsample;
}


int main(int argc, char* argv[]) {
    TH1::AddDirectory(0);
    const TString kDataFile("output_data_0.root");
    const TString kMCFile("output_mc_0.root");

    const TString kHistName("htxz_13");
    const UInt_t kNsamples = 1000000;

    TFile* fd = TFile::Open(kDataFile);
    TFile* fm = TFile::Open(kMCFile);

    TH1D* hd = static_cast<TH1D*>(fd->Get(kHistName));
    TH1D* hm = static_cast<TH1D*>(fm->Get(kHistName));
    fd->Close();
    fm->Close();

    BootstrapFCN fcn(kNsamples, hm, hd);

    ROOT::Minuit2::MnUserParameters upar;
    Double_t scale = 1.0;
    Double_t shift = hm->GetMean() - hd->GetMean();
    upar.Add("scale", scale, 0.1);
    upar.SetLimits("scale", 0.0, 100.0);
    upar.Add("shift", shift, 0.1);
 
    // create MIGRAD minimizer
    ROOT::Minuit2::MnMigrad migrad(fcn, upar);
    ROOT::Minuit2::FunctionMinimum result = migrad(); 
    std::cout << "minimum: " << result << std::endl;
    scale = result.UserParameters().Params()[0];
    shift = result.UserParameters().Params()[1];

    TFile* fout = TFile::Open("~/test.root", "recreate");
    TCanvas* c = new TCanvas("c", "c");
    // c->SetCanvasSize(5, 4);
    c->Divide(1, 2);
    c->cd(1);
    THStack* hs = new THStack();
    // hd->Scale(1 / hd->Integral());
    hm->Scale(hd->Integral() / hm->Integral());
    hs->Add(hm, "histE");
    hs->Add(hd, "EP");

    TH1D* hsample = fcn.GetHSample(scale, shift);
    hsample->Scale(hd->Integral() / hsample->Integral());
    hs->Add(hsample, "hist");

    hd->SetLineColor(kBlack);
    hm->SetLineColor(kRed);
    hs->Draw("nostack");
    hs->GetYaxis()->SetTitle("Events");
    hs->GetXaxis()->SetTitle("Width");
    hs->Write();

    c->cd(2);
    TH1D* hdr = static_cast<TH1D*>(hd->Clone());
    TH1D* hsampler = static_cast<TH1D*>(hsample->Clone());
    for (Int_t i = 1; i <= hd->GetNbinsX(); i++) {
        if (hm->GetBinContent(i) == 0) {
            hdr->SetBinContent(i, 0.);
            hdr->SetBinError(i, 0.);
            hsampler->SetBinContent(i, 0.);
        }
        else { 
            hdr->SetBinContent(i, hdr->GetBinContent(i) / hm->GetBinContent(i));
            hdr->SetBinError(i, hdr->GetBinError(i) / hm->GetBinContent(i));
            hsampler->SetBinContent(i, hsampler->GetBinContent(i) / hm->GetBinContent(i));
        }
    }
    THStack* hsr = new THStack();
    hsr->Add(hdr, "EP");
    hsr->Add(hsampler, "hist");
    hsr->Draw("nostack");
    hsr->GetYaxis()->SetTitle("Ratio to MC");
    hsr->GetXaxis()->SetTitle("Width");
    hsr->Write();

    c->cd(1);
    gPad->SetLogy();
    c->Draw();
    c->Write();
}
