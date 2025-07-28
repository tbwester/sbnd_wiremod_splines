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
        result[0] = mean;
        result[1] = sd / TMath::Sqrt(h->GetEntries());
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
        if (h->GetBinCenter(i) < median[0] + sig_down * sd || h->GetBinCenter(i) > median[0] + sig_up * sd) continue;
        hnew->SetBinContent(i, h->GetBinContent(i));
        hnew->SetBinError(i, h->GetBinError(i));
    }
    
    iterative_truncated_mean(hnew, sig_down, sig_up, tol, result);
    delete hnew;
    return;
}


void test_it_tr_mn() {
    TH1F* h = new TH1F("h", "", 100, -10, 10);
    h->FillRandom("gaus", 1000);

    Double_t result[2] = { -999, 0. };

    iterative_truncated_mean(h, -2, 2, 1.0e-4, result);

    printf("%.4e, %.4e\n", result[0], result[1]);
}
