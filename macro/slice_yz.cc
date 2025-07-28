void hist_mean_unc(const TH1* h, Double_t tol, Double_t* result) {
    Int_t nbins = h->GetNbinsX();
    if (h->Integral() == 0) return; 

    const Int_t kNthrows = 1000;
    Double_t means[kNthrows];

    TRandom3 rng;
    for (Int_t i = 0; i < kNthrows; i++) {
        TH1* htmp = (TH1*)h->Clone();
        htmp->Reset();
        for (Int_t k = 1; k <= nbins; k++) {
            Double_t bc = h->GetBinContent(k) + rng.Gaus(0, h->GetBinError(k));
            htmp->SetBinContent(k, TMath::Max(bc, 0.));
        }
        means[i] = htmp->GetMean();
        delete htmp;
    }

    result[0] = TMath::Mean(kNthrows, means);
    result[1] = TMath::RMS(kNthrows, means);
}


void iterative_truncated_mean(TH1* h, Double_t sig_down, Double_t sig_up, Double_t tol, Double_t* result) {
    if (sig_down > sig_up) {
        fprintf(stderr, "Warning: reversing iterative truncated mean limits [%.2e,%.2e]\n", sig_down, sig_up);
        Double_t tmp = sig_down;
        sig_down = sig_up;
        sig_up = tmp;
    }

    Double_t mean = h->GetMean();
    Double_t sd = h->GetRMS();
    std::cout << mean << ", " << sd << ", " << result[0] << "\n";

    // return mean, std err of iterative truncated mean
    if (TMath::Abs(result[0] - mean) < tol) {
        hist_mean_unc(h, 1.0e-6, result);
        return;
    }
    result[0] = mean;
    result[1] = sd;
    
    // get positions to cut based on quantiles
    Double_t probs[1] = { 0.5 };
    Double_t median[1];
    h->GetQuantiles(1, median, probs);

    std::cout << "median: " << median[0] << "\n";
    std::cout << median[0] + sig_down * sd << ", " << median[0] + sig_up * sd << "\n";
    
    TH1* hnew = (TH1*)h->Clone();
    hnew->Reset();
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        if (h->GetBinLowEdge(i) > median[0] + sig_up * sd || h->GetBinLowEdge(i) + h->GetBinWidth(i) < median[0] + sig_down * sd) continue;
        hnew->SetBinContent(i, h->GetBinContent(i));
        hnew->SetBinError(i, h->GetBinError(i));
    }
    
    iterative_truncated_mean(hnew, sig_down, sig_up, tol, result);
    delete hnew;
    return;
}



void slice_yz(const TString& input_filename, const TString& hist_name, const TString& output_filename) {
    TH1::AddDirectory(0);
    // single file or file list
    const UInt_t kNdims = 7;

    const TString kLabels[kNdims] = { "x", "y", "z", "txz", "tyz", "dqdx", "width" };
    const TString kTitles[kNdims] = { "x (cm)", "y (cm)", "z (cm)", "ThetaXZ (deg)", "ThetaYZ (deg)", "dQ/dx", "Width" };
    TPRegexp re_plane_idx("hwidth([0-9])");
    TObjArray* matches = re_plane_idx.MatchS(hist_name);
    const Int_t plane_idx = std::stoi( static_cast<TObjString*>(matches->At(1))->GetString().Data());
    matches->Delete();

    TFile* fin = TFile::Open(input_filename, "read");
    THnSparseD* h = (THnSparseD*)fin->Get(hist_name);
    TH2D* proj_ntrk = (TH2D*)fin->Get(Form("hntrk_%d_%s", plane_idx, kLabels[0].Data()));
    fin->Close();
    
    TFile* fout = TFile::Open(output_filename, "recreate");

    const UInt_t kNbinsY = 4;
    const UInt_t kNbinsZ = 5;
    const Double_t kBinEdgesY[kNbinsY + 1] = { -300, -100, 0, 100, 300 };
    const Double_t kBinEdgesZ[kNbinsZ + 1] = { -100, 100, 200, 300, 400, 600 };

    for (UInt_t iy = 0; iy < kNbinsY; iy++) {
        h->GetAxis(1)->SetRangeUser(kBinEdgesY[iy], kBinEdgesY[iy + 1]);
        for (UInt_t iz = 0; iz < kNbinsZ; iz++) {
            h->GetAxis(2)->SetRangeUser(kBinEdgesZ[iz], kBinEdgesZ[iz + 1]);
            printf("Y=(%.2f, %.2f), Z=(%.2f, %.2f)\n", 
                    kBinEdgesY[iy], kBinEdgesY[iy + 1],
                    kBinEdgesZ[iz], kBinEdgesZ[iz + 1]);

            // x vs width
            TH2D* proj = h->Projection(kNdims - 1, 0);
            proj->SetTitle(Form(" Y=(%.2f, %.2f), Z=(%.2f, %.2f);%s;%s", 
                        kBinEdgesY[iy], kBinEdgesY[iy + 1],
                        kBinEdgesZ[iz], kBinEdgesZ[iz + 1],
                        kTitles[0].Data(), kTitles[kNdims - 1].Data()));
            proj->SetName(Form("h2d_%d_%d", iy, iz));
            proj->Write();

            const UInt_t nbinsx = proj->GetNbinsX();
            const UInt_t nbinsy = proj->GetNbinsY();

            std::vector<float> xs(nbinsx);
            std::vector<float> xerrs(nbinsx);
            std::vector<float> ys(nbinsx);
            std::vector<float> yerrs(nbinsx);

            TH1D* hslice = proj->ProjectionY();
            for (Int_t ii = 1; ii <= nbinsx; ii++) {
                hslice->Reset();

                xs.at(ii - 1) = proj->GetXaxis()->GetBinCenter(ii);
                xerrs.at(ii - 1) = proj->GetXaxis()->GetBinWidth(ii) / 2.0;
                for (Int_t jj = 1; jj <= nbinsy; jj++) {
                    Float_t val = proj->GetBinContent(ii, jj);
                    Int_t ntrk = proj_ntrk->GetBinContent(ii, jj);
                    if (ntrk <= 0) {
                        continue;
                    }
                    hslice->SetBinContent(jj, val);
                    Float_t hits_per_trk = val / (float)ntrk;
                    Float_t err = TMath::Sqrt(ntrk * hits_per_trk * (1. + hits_per_trk));
                    hslice->SetBinError(jj, err);
                }

                Float_t itm = 0.; // iterative truncated mean (ITM)
                Float_t itm_unc = 0.; // uncertainty of ITM
                Double_t itm_result[2];
                iterative_truncated_mean(hslice, -2, 2, 1.0e-4, itm_result);
                itm = itm_result[0];
                itm_unc = itm_result[1];
                ys.at(ii - 1) = itm;
                yerrs.at(ii - 1) = itm_unc;
            }

            TGraphErrors* g = new TGraphErrors(nbinsx, &xs[0], &ys[0], &xerrs[0], &yerrs[0]);
            g->SetTitle(Form(";%s;%s", kTitles[0].Data(), kTitles[kNdims - 1].Data()));
            g->SetName(Form("gitm_%d_%d_%s_vs_%s", iy, iz,
                        kLabels[0].Data(), kLabels[kNdims - 1].Data()));
            g->Write();
            delete hslice;

            TH2D* projyz = h->Projection(1, 2);
            projyz->SetTitle(Form(" Y=(%.2f, %.2f), Z=(%.2f, %.2f);%s;%s", 
                        kBinEdgesY[iy], kBinEdgesY[iy + 1],
                        kBinEdgesZ[iz], kBinEdgesZ[iz + 1],
                        kTitles[2].Data(), kTitles[1].Data()));
            projyz->SetName(Form("h2dyz_%d_%d", iy, iz));
            projyz->Write();
        }
    }

    fout->Close();
}
