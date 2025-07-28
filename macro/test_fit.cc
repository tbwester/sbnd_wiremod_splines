// https://root.cern/doc/master/langaus_8C.html
/*
double langaufun(double *x, double *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}
*/


void iterative_truncated_mean(Int_t n, Double_t* data, Double_t sig_down, Double_t sig_up, Double_t tol, Double_t* result) {
    if (sig_down > sig_up) {
        fprintf(stderr, "Warning: reversing iterative truncated mean limits [%.2e,%.2e]\n", sig_down, sig_up);
        Double_t tmp = sig_down;
        sig_down = sig_up;
        sig_up = tmp;
    }

    Double_t mean = TMath::Mean(n, data);
    Double_t sd = TMath::RMS(n, data);
    std::cout << mean << ", " << sd << ", " << result[0] << "\n";

    // return mean, std err of iterative truncated mean
    if (TMath::Abs(result[0] - mean) < tol) {
        result[0] = mean;
        result[1] = sd / TMath::Sqrt(n);
        return;
    }
    result[0] = mean;
    result[1] = sd;
    
    // get positions to cut based on quantiles
    Double_t probs[1] = { 0.5 };
    Double_t median[1];
    TMath::Quantiles(n, 1, data, median, probs, false);

    std::cout << "median: " << median[0] << "\n";
    std::cout << median[0] + sig_down * sd << ", " << median[0] + sig_up * sd << "\n";
    
    std::vector<Double_t> cut_data;
    cut_data.reserve(n);
    Int_t cut_count = 0;
    for (int i = 0; i < n; i++) {
        if (data[i] < median[0] + sig_down * sd || data[i] > median[0] + sig_up * sd) continue;
        cut_count++;
        cut_data.push_back(data[i]);
    }
    
    iterative_truncated_mean(cut_count, &cut_data[0], sig_down, sig_up, tol, result);
    return;
}


// return the mean & uncertainty on the mean based on bins in a histogram, accounting for each bin's error
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


Double_t langaufun(Double_t *x, Double_t *par) {
    Double_t invsq2pi = 0.398942280401;// Control constants
    // Double_t mpshift = -0.22278298l
    Double_t np = 500.0;
    Double_t sc = 5.0; // convolution extends to +-sc Gaussian sigmas
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;

    mpc=par[1];
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp-xlow)/np;

    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}


void test_fit(const TString& input_filename, const TString& hist_name, const TString& output_filename) {

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
    // THnSparseI* hi = (THnSparseI*)fin->Get("hntrk0");
    
    TFile* fout = TFile::Open(output_filename, "recreate");
    for (UInt_t i = 0; i < kNdims - 1; i++) {
        TH2D* proj = h->Projection(kNdims - 1, i);
        proj->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        proj->SetName(Form("h2d_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        proj->Write();

        TProfile* prof = proj->ProfileX();
        prof->SetName(Form("hp_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        prof->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        prof->Write();

        /*
        TH2D* proj_ntrk = hi->Projection(kNdims - 1, i);
        proj_ntrk->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        proj_ntrk->SetName(Form("h2d_ntrk_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        proj_ntrk->Write();

        TProfile* prof_ntrk = proj_ntrk->ProfileX();
        prof_ntrk->SetName(Form("hp_ntrk_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        prof_ntrk->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        prof_ntrk->Write();

        TH1D* h1d_ntrk = proj->ProjectionX();
        h1d_ntrk->SetName(Form("hs_ntrk_%s", kLabels[i].Data()));
        h1d_ntrk->SetTitle(Form(";%s;Entries", kTitles[i].Data()));
        h1d_ntrk->Write();
        */

        // number of unique tracks that contributed to each bin of hits
        TH2D* proj_ntrk = (TH2D*)fin->Get(Form("hntrk_%d_%s", plane_idx, kLabels[i].Data()));


        // make a graph of mean +- 1 sigma
        // each will use the same points along x
        const UInt_t nbinsx = proj->GetNbinsX();
        const UInt_t nbinsy = proj->GetNbinsY();
        std::vector<float> xs(nbinsx);
        std::vector<float> xerrs(nbinsx);

        // mean +- SD
        std::vector<float> ys(nbinsx);
        std::vector<float> yerrs(nbinsx);

        // laundau fit
        std::vector<float> ys2(nbinsx);
        std::vector<float> yerrs2(nbinsx);
        std::vector<float> yerrs3(nbinsx);

        // iterative truncated mean (ITM)
        std::vector<float> ys3(nbinsx);
        std::vector<float> yerrs4(nbinsx);

        // chi2 of launda fit
        std::vector<float> chisq(nbinsx);

        TH1D* h1d = proj->ProjectionX();
        TH1D* h1d_ntrk = proj_ntrk->ProjectionX();
        for (Int_t j = 1; j <= h1d->GetNbinsX(); j++) {
            Float_t ntrk = (float)h1d_ntrk->GetBinContent(j); 

            // variance of poisson-distributed number with parameter lambda_1
            // of poisson random variables with parameter lambda_2 is
            // lambda_1 * lambda_2 * (1 + lambda_2)
            Float_t hits_per_trk = h1d->GetBinContent(j) / ntrk;
            Float_t err = TMath::Sqrt(ntrk * hits_per_trk * (1. + hits_per_trk));
            h1d->SetBinError(j, err);
        }

        h1d->SetName(Form("hs_%s", kLabels[i].Data()));
        h1d->SetTitle(Form(";%s;Entries", kTitles[i].Data()));
        h1d->Write();

        TH1D* hslice = proj->ProjectionY();

        // store the fit result after each fit. If we need to retry, we'll use
        // these, since we don't expect large changes between slices
        Double_t prev_fit_parms[4];

        for (Int_t ii = 1; ii <= nbinsx; ii++) {
            xs.at(ii - 1) = proj->GetXaxis()->GetBinCenter(ii);
            xerrs.at(ii - 1) = proj->GetXaxis()->GetBinWidth(ii) / 2.0;

            hslice->Reset();

            Float_t sum = 0;
            Float_t sum2 = 0;
            Float_t max = -1;
            Float_t max_pos = -1;
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
                sum += val;
                sum2 += val * val;
                if (val > max) {
                    max = val;
                    max_pos = hslice->GetBinCenter(jj);
                }
            }

            Float_t mean = 0.;
            Float_t stddev = 0.;
            Float_t itm = 0.; // iterative truncated mean (ITM)
            Float_t itm_unc = 0.; // uncertainty of ITM
            Float_t mpv = 0.;
            Float_t mpv_unc = 0.; // ROOT fit error
            Float_t mpv_unc2 = 0.; // Landau width
            Float_t chi2 = 0.;
            int fit_status = -1;

            Double_t itm_result[2];
            iterative_truncated_mean(hslice, -2, 2, 1.0e-4, itm_result);
            itm = itm_result[0];
            itm_unc = itm_result[1];

            // set pointer outside so we can delete it after write
            TF1* best_fitfunc = 0;
            Float_t d_best_chi2 = -1;
            if (hslice->Integral() != 0) {
                mean = hslice->GetMean();
                stddev = hslice->GetStdDev();

                // fitting function, limits are set below based on CDF begin
                // fitting slightly after start of slice, to remove small
                // feature which can bias the landau fit
                const Double_t kQuantiles[2] = { 0.001, 0.999 };
                Double_t limits[2];
                hslice->GetQuantiles(2, limits, kQuantiles);

                const auto run_fit = [&](const Double_t* parms) {
                    TF1* fitfunc = new TF1("fitfunc", langaufun, limits[0], limits[1], 4);
                    for (int kk = 0; kk < 4; kk++) {
                        fitfunc->SetParameter(kk, parms[kk]);
                    }
                    fit_status = hslice->Fit(fitfunc, "RBLQE");

                    // fit failed, don't set any results
                    if (fit_status != 0 || fitfunc->GetNDF() == 0) {
                        delete fitfunc;
                        return;
                    }

                    Float_t d_this_chi2 = TMath::Abs(
                            (fitfunc->GetChisquare() / fitfunc->GetNDF()) - 1);
                    if (d_this_chi2 < d_best_chi2 || d_best_chi2 == -1) {
                        // update best fit
                        if (best_fitfunc) delete best_fitfunc;
                        best_fitfunc = fitfunc;
                        chi2 = fitfunc->GetChisquare() / fitfunc->GetNDF();
                        mpv = fitfunc->GetParameter(1);
                        mpv_unc = fitfunc->GetParError(1);
                        mpv_unc2 = fitfunc->GetParameter(0) / 2.0;
                        d_best_chi2 = d_this_chi2;

                        printf("New best fit at %s=%.2f (%d, chi2=%.2e)\n",
                                kLabels[i].Data(), xs[ii - 1], ii, chi2);

                        for (int kk = 0; kk < 4; kk++) {
                            printf(" - parm %d: in: %.2e out: %.2e\n",
                                    kk, parms[kk], fitfunc->GetParameter(kk));
                        }
                    }
                    else {
                        delete fitfunc;
                    }
                };

                const Double_t init_parms[4] = {
                    stddev, max_pos, sum / proj->GetXaxis()->GetBinWidth(0), 0.3
                };
                run_fit(init_parms);
                if (chi2 > 10 || fit_status != 0 || chi2 == 0) {
                    // fit was probably bad
                    if (ii > 1) {
                        fprintf(stderr, "Retrying fit at %s=%.2f (chi2=%.2e)\n",
                                kLabels[i].Data(), xs[ii - 1], chi2);
                        run_fit(prev_fit_parms);
                    }
                }
                if (best_fitfunc) {
                    for (int kk = 0; kk < 4; kk++) {
                        prev_fit_parms[kk] = best_fitfunc->GetParameter(kk);
                    }
                }
                else {
                    fprintf(stderr, "Fit failed at %s=%.2f\n",
                            kLabels[i].Data(), xs[ii - 1]);
                }
            }

            ys.at(ii - 1) = mean;
            yerrs.at(ii - 1) = stddev;
            ys2.at(ii - 1) = mpv;
            yerrs2.at(ii - 1) = mpv_unc;
            yerrs3.at(ii - 1) = mpv_unc2;
            chisq.at(ii - 1) = chi2;

            ys3.at(ii - 1) = itm;
            yerrs4.at(ii - 1) = itm_unc;

            // if (i == 0) {
            TH1F* hh = (TH1F*)hslice->Clone(Form("h%s_%d", kLabels[i].Data(), ii));
            hh->Write();
            // }

            if (best_fitfunc) delete best_fitfunc;
        }

        TGraphErrors* g = new TGraphErrors(nbinsx, &xs[0], &ys[0], &xerrs[0], &yerrs[0]);
        TGraphErrors* g2a = new TGraphErrors(nbinsx, &xs[0], &ys2[0], &xerrs[0], &yerrs2[0]);
        TGraphErrors* g2b = new TGraphErrors(nbinsx, &xs[0], &ys2[0], &xerrs[0], &yerrs3[0]);
        TGraphErrors* g2c = new TGraphErrors(nbinsx, &xs[0], &ys3[0], &xerrs[0], &yerrs4[0]);
        TGraph* g3 = new TGraph(nbinsx, &xs[0], &chisq[0]);

        // g->Draw("ALP");
        g->SetName(Form("g_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        g->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        g2a->SetName(Form("g2a_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        g2a->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        g2b->SetName(Form("g2b_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        g2b->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        g2c->SetName(Form("gitm_%s_vs_%s", kLabels[i].Data(), kLabels[kNdims - 1].Data()));
        g2c->SetTitle(Form(";%s;%s", kTitles[i].Data(), kTitles[kNdims - 1].Data()));
        g3->SetName(Form("g_%s_chi2", kLabels[i].Data()));
        g3->SetTitle(Form(";%s;%s", kTitles[i].Data(), "#chi^{2}/NDF"));

        g->Write();
        g2a->Write();
        g2b->Write();
        g2c->Write();
        g3->Write();

        delete hslice;
        delete h1d;
    }

    delete fout;
}
