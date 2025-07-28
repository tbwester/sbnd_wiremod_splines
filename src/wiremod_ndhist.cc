/*
 * Create THnSparse containing hit widths along different dimensions
 * x, y, z, ThetaXZ, ThetaYZ
 * Input: Calibration ntuples. Select T0-tagged tracks
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "Math/Vector3D.h"

#include "SCECorr.h"


std::vector<TString> filenames_from_input(const TString&, int);
TString basename_prefix(const TString&, const TString& prefix="", const TString& suffix="");
bool is_int(Float_t);
double lifetime_correction(double, double);


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;
const UInt_t kNdims = 7;
const Float_t kTrackCut = 60.; // cm

const TString kLabels[kNdims] = { "x", "y", "z", "txz", "tyz", "dqdx", "width" };
const TString kTitles[kNdims] = { "x (cm)", "y (cm)", "z (cm)", "ThetaXZ (deg)", "ThetaYZ (deg)", "dQ/dx", "Width" };
const Int_t kNbins[kNdims] = { 60, 60, 70, 30, 30, 50, 200 };
const Double_t kXmin[kNdims] = { -300, -300, -100, 0, -180, 0, 0 };
const Double_t kXmax[kNdims] = { 300, 300, 600, 180, 180, 5000, 20 };

void usage() {
    printf(" ./wiremod_ndhist [data|mc] [sce] <filename or file list>\n");
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        usage();
        exit(1);
    }
    TH1::AddDirectory(0);

    bool is_data = false;
    if (!strcmp(argv[1], "data")) {
        is_data = true;
    }
    printf("wiremod_ndhist: %s mode\n", is_data ? "DATA" : "MC");

    bool do_sce = false;
    if (!strcmp(argv[2], "sce")) {
        do_sce = true;
    }
    printf("wiremod_ndhist: SCE will %sbe applied\n", do_sce ? "" : "NOT ");

    bool do_elifetime = false;
    float e_lifetime = -1.0;
    if (argc > 4) {
        do_elifetime = true;
        e_lifetime = std::stof(argv[3]);
        printf("wiremod_ndhist: Electron lifetime correction will be applied, tau=%.2e\n", e_lifetime);
    }
    else {
        printf("wiremod_ndhist: Electron lifetime correction will NOT be applied\n");
    }

    bool do_yz = false;
    TH2F* CzyHist_sce[kNplanes][2];
    if (argc > 5) {
        do_yz = true;
        TFile* file_SCEYZ = TFile::Open(argv[4]);
        printf("wiremod_ndhist: Loading YZ nonuniformity correction histograms"
                " from %s\n", argv[4]);
        for (int l = 0; l < kNplanes; l++){
            for (int k = 0; k < 2; k++) {
                CzyHist_sce[l][k] = (TH2F*)file_SCEYZ->Get(Form("CzyHist_%i_%i",l,k));
            }
        }
        file_SCEYZ->Close();
        std::cout << CzyHist_sce[0][0] << "\n";
    }
    printf("wiremod_ndhist: YZ nonuniformity correction "
            "will %sbe applied\n", do_yz ? "" : "NOT ");

    // single file or file list
    const TString input_arg(argv[argc - 1]);
    const std::vector<TString> kFilenames = filenames_from_input(input_arg, -1);
    const TString output_filename(basename_prefix(kFilenames.at(0), "out_"));

    SCECorr* sce_corr = nullptr;
    if (do_sce) {
        printf("wiremod_ndhist: loading SCE TH3\n");
        sce_corr = new SCECorr(is_data);
        sce_corr->ReadHistograms();
        printf("wiremod_ndhist: loaded SCE TH3 complete\n");
    }

    // 1 hist per plane per TPC. We also keep track of the number of tracks in
    // each eventual projection bin using TH2Is
    THnSparseD* h[kNplanes * kNTPCs];
    TH2I* hi[kNplanes * kNTPCs * kNdims];
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        h[i] = new THnSparseD(Form("hwidth%d", i), "", kNdims, kNbins, kXmin, kXmax);
        for (unsigned j = 0; j < kNdims; j++) {
            hi[i * kNdims + j] = new TH2I(Form("hntrk_%d_%s", i, kLabels[j].Data()), "",
                    kNbins[j], kXmin[j], kXmax[j],
                    kNbins[kNdims - 1], kXmin[kNdims - 1], kXmax[kNdims - 1]);
        }
    }

    size_t nfiles_pass = 0;
    size_t nfiles_all = 0;
    size_t nevts = 0;
    size_t track_counter = 0;
    bool used = false;
    for (const auto filename : kFilenames) {
        nfiles_all++;
        TFile* f = TFile::Open(filename, "read");
        if (!f) {
            // ROOT will have thrown an error
            continue;
        }
        if (f->IsZombie()) {
            fprintf(stderr, "Could not open file %s.\n", filename.Data());
            continue;
        }

        nfiles_pass++;
        TTreeReader reader("caloskim/TrackCaloSkim", f);
        TTreeReaderValue<int> selected(reader, "trk.selected");
        TTreeReaderValue<int> whicht0(reader, "trk.whicht0");
        
        TTreeReaderValue<int> run(reader, "meta.run");
        TTreeReaderValue<int> evt(reader, "meta.evt");
        TTreeReaderValue<int> subrun(reader, "meta.subrun");
        /*
        // for hit train study
        // --------
        TTreeReaderArray<UShort_t> wire[kNplanes] = {
            { reader, "trk.hits0.h.wire" }, { reader, "trk.hits1.h.wire" }, { reader, "trk.hits2.h.wire" },
        };
        TTreeReaderArray<float> time[kNplanes] = {
            { reader, "trk.hits0.h.time" }, { reader, "trk.hits1.h.time" }, { reader, "trk.hits2.h.time" },
        };
        // --------
        // */


        // track length cut: use residual range on collection only
        TTreeReaderArray<float> rr2(reader, "trk.hits2.rr");

        TTreeReaderValue<float> trk_dirx(reader, "trk.dir.x");
        TTreeReaderValue<float> trk_diry(reader, "trk.dir.y");
        TTreeReaderValue<float> trk_dirz(reader, "trk.dir.z");

        // this is super ugly, but pointers are worse for accessing TTreeReader below
        TTreeReaderArray<unsigned short> tpc[kNplanes] = {
            { reader, "trk.hits0.h.tpc" }, { reader, "trk.hits1.h.tpc" }, { reader, "trk.hits2.h.tpc" },
        };
        TTreeReaderArray<float> goodness[kNplanes] = {
            { reader, "trk.hits0.h.goodness" }, { reader, "trk.hits1.h.goodness" }, { reader, "trk.hits2.h.goodness" },
        };
        TTreeReaderArray<float> x[kNplanes] = {
            { reader, "trk.hits0.h.sp.x" }, { reader, "trk.hits1.h.sp.x" }, { reader, "trk.hits2.h.sp.x" },
        };
        TTreeReaderArray<float> y[kNplanes] = {
            { reader, "trk.hits0.h.sp.y" }, { reader, "trk.hits1.h.sp.y" }, { reader, "trk.hits2.h.sp.y" },
        };
        TTreeReaderArray<float> z[kNplanes] = {
            { reader, "trk.hits0.h.sp.z" }, { reader, "trk.hits1.h.sp.z" }, { reader, "trk.hits2.h.sp.z" },
        };
        TTreeReaderArray<float> width[kNplanes] = {
            { reader, "trk.hits0.h.width" }, { reader, "trk.hits1.h.width" }, { reader, "trk.hits2.h.width" },
        };
        TTreeReaderArray<bool> ontraj[kNplanes] = {
            { reader, "trk.hits0.ontraj" }, { reader, "trk.hits1.ontraj" }, { reader, "trk.hits2.ontraj" },
        };
        TTreeReaderArray<float> dirx[kNplanes] = {
            { reader, "trk.hits0.dir.x" }, { reader, "trk.hits1.dir.x" }, { reader, "trk.hits2.dir.x" },
        };
        TTreeReaderArray<float> diry[kNplanes] = {
            { reader, "trk.hits0.dir.y" }, { reader, "trk.hits1.dir.y" }, { reader, "trk.hits2.dir.y" },
        };
        TTreeReaderArray<float> dirz[kNplanes] = {
            { reader, "trk.hits0.dir.z" }, { reader, "trk.hits1.dir.z" }, { reader, "trk.hits2.dir.z" },
        };
        TTreeReaderArray<float> dqdx[kNplanes] = {
            { reader, "trk.hits0.dqdx" }, { reader, "trk.hits1.dqdx" }, { reader, "trk.hits2.dqdx" },
        };

        // keep track of how many tracks go into the hit distribution
        // since, e.g., all hits from 1 track will get the same angle, the
        // errors on angle will go like sqrt(NTracks) instead of sqrt(NHits)
        std::vector<TH2I*> hflag(kNplanes * kNTPCs * kNdims, nullptr);
        for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
            for (unsigned j = 0; j < kNdims; j++) {
                hflag.at(i * kNdims + j) = (TH2I*)hi[i * kNdims + j]->Clone(Form("hflag_%d_%d", i, j));
            }
        }

        int track_idx = 0;
        while (reader.Next()) {
            track_idx++;
            if (*selected < 1) continue;

            // For CRT T0 study -- only select CRT T0-tagged tracks
            if (*whicht0 != 1) continue;

            // skip short tracks
            size_t nhits = rr2.GetSize();
            if (nhits == 0) {
                fprintf(stderr, "Warning: Selected track (idx=%d, selected=%d) with no hits? Run=%d, Subrun=%d, Evt=%d. Skipping!\n", track_idx, *selected, *run, *subrun, *evt);
                continue;
            }
            if (rr2[nhits - 1] < kTrackCut) continue;

            track_counter++;
            ROOT::Math::XYZVector trk_dir(*trk_dirx, *trk_diry, *trk_dirz);
            float trk_thxz = trk_dir.Theta() * 180. / TMath::Pi();
            float trk_thyz = trk_dir.Phi() * 180. / TMath::Pi();

            // reset track counting flags
            for (unsigned i = 0; i < kNplanes * kNTPCs * kNdims; i++) {
                hflag[i]->Reset();
            }

            for (UInt_t ip = 0; ip < kNplanes; ip++) {
                for (size_t i = 0; i < x[ip].GetSize(); i++) {
                    // skip nans
                    if (x[ip][i] != x[ip][i]) continue;
                    
                    // skip not on track
                    if (!ontraj[ip][i]) continue;

                    // goodness cut
                    if (goodness[ip][i] >= 100.) continue;

                    // hit trains have widths in increments of exactly 0.5
                    // skip hits from these
                    if (is_int(width[ip][i] * 2)) continue;

                    nevts++;

                    // for hit train study
                    /*
                    if (width[ip][i] > 10.4 && width[ip][i] < 10.7) {
                        printf("Unusually wide hit (w=%.1f, t=%.1f)! Run=%d, Evt=%d, Subrun=%d wire=%d, plane=%d TPC=%d\n",
                                width[ip][i], time[ip][i], *run, *evt, *subrun, wire[ip][i], ip, tpc[ip][i]);
                    }
                    */

                    XYZVector sp(x[ip][i], y[ip][i], z[ip][i]);
                    double dqdx_hit = dqdx[ip][i];
                    double total_correction = 1.0;
                    double sce_correction = 1.0;
                    double e_correction = 1.0;
                    double yz_correction = 1.0;

                    if (do_sce) {
                        // SCE correction
                        sp = sce_corr->WireToTrajectoryPosition(sp);
                        double pitch_sce_uncorr = sce_corr->meas_pitch(x[ip][i], y[ip][i], z[ip][i], dirx[ip][i], diry[ip][i], dirz[ip][i], ip, false);
                        double pitch_sce_corr = sce_corr->meas_pitch(x[ip][i], y[ip][i], z[ip][i], dirx[ip][i], diry[ip][i], dirz[ip][i], ip, true);
                        dqdx_hit *= pitch_sce_uncorr / pitch_sce_corr;
                        total_correction *= pitch_sce_uncorr / pitch_sce_corr;
                        sce_correction = pitch_sce_uncorr / pitch_sce_corr;
                    }

                    if (do_elifetime) {
                        dqdx_hit *= lifetime_correction(x[ip][i], e_lifetime);
                        total_correction *= lifetime_correction(x[ip][i], e_lifetime);
                        e_correction = lifetime_correction(x[ip][i], e_lifetime);
                    }

                    if (do_yz) {
                        const int nbinz=100, nbinx=40, nbiny=80;
                        const float lowz=0, highz=500, lowx=-200, highx=200, lowy=-200, highy=200;

                        int ibinx = floor((sp.X() - lowx) / (highx - lowx) * nbinx);
                        if (ibinx < 0 || ibinx >= nbinx)continue;
                        int ibiny = floor((sp.Y() - lowy) / (highy - lowy) * nbiny);
                        if (ibiny < 0 || ibiny >= nbiny)continue;
                        int ibinz = floor((sp.Z() - lowz) / (highz - lowz) * nbinz);
                        if (ibinz < 0 || ibinz >= nbinz)continue;

                        // YZ non-uniformity
                        double CF_zy = CzyHist_sce[ip][tpc[ip][i]]->GetBinContent(ibinz+1, ibiny+1);
                        dqdx_hit *= CF_zy;
                        total_correction *= CF_zy;
                        yz_correction = CF_zy;
                    }
                    printf("wiremod_ndhist: total_correction=%.6e (%.6e, %.6e, %.6e)\n",
                            total_correction, sce_correction, e_correction, yz_correction);

                    Double_t val[kNdims]  = {
                        sp.X(), sp.Y(), sp.Z(), trk_thxz, trk_thyz, dqdx_hit, width[ip][i] 
                    };

                    // select by TPC
                    unsigned hit_idx = ip + kNplanes * tpc[ip][i];
                    h[hit_idx]->Fill(val);

                    // count track up to once per bin
                    for (unsigned j = 0; j < kNdims; j++) {
                        unsigned val_idx = hit_idx * kNdims + j;
                        Int_t target_bin = hflag[val_idx]->FindFixBin(val[j], val[kNdims - 1]);
                        if (hflag[val_idx]->GetBinContent(target_bin) == 0) {
                            hi[val_idx]->Fill(val[j], val[kNdims - 1]);
                            hflag[val_idx]->Fill(val[j], val[kNdims - 1]);
                        }
                        assert(hflag[val_idx]->GetBinContent(target_bin) != 0);
                    }
                } // loop over hits
            } // loop over planes
        } // loop over events
        
        for (unsigned i = 0; i < kNplanes * kNTPCs * kNdims; i++) {
            delete hflag.at(i);
        }

        delete f;
    }

    printf("Processed %lu tracks (%lu hits)\n", track_counter, nevts);
    
    if (nfiles_pass == 0) {
        fprintf(stdout, "No files processed.\n");
        return 0;
    }
    
    TFile* fout = TFile::Open(output_filename, "recreate");
    for (unsigned i = 0; i < kNplanes * kNTPCs; i++) {
        h[i]->Write();
        for (unsigned j = 0; j < kNdims; j++) {
            hi[i * kNdims + j]->Write();
        }
    }

}


std::vector<TString> filenames_from_input(const TString& input_arg, int nmax=-1) {
    std::vector<TString> filenames;
    if (input_arg.EndsWith(".root")) {
        filenames.push_back(input_arg);
        return filenames;
    }

    // read from file
    size_t nfiles = 0;
    std::ifstream ifile(input_arg);
    std::string line;
    while (std::getline(ifile, line)) {
        nfiles++;
        fprintf(stdout, "Adding file %zu: %s...\n", nfiles, line.c_str());
        filenames.push_back(TString(line));
        if (nfiles >= nmax && nmax > 0) break;
    }
    return filenames;
}


TString basename_prefix(const TString& input, const TString& prefix, const TString& suffix) {
    // remove path from filename and return new string with prefix or suffix added before extension
    TPRegexp re(".*/(.*)");
    TObjArray* matches = re.MatchS(input);
    TString result((static_cast<TObjString*>(matches->At(1)))->String());
    matches->Delete();
    return prefix + result;
}


bool is_int(Float_t val) {
    return std::abs(roundf(val) - val) < 0.00001f;
}


double lifetime_correction(double x, double tau) {
    static const double v_drift = 156.267;
    double out = 1.;
    if(fabs(x) > 200.) return out;

    double this_tdrift = (200. - fabs(x)) / v_drift;
    out = 1. / exp(-1. * this_tdrift / tau);

    return out;
}
