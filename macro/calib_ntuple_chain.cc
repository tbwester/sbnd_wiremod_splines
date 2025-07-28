{
    TString input_dir("/pnfs/sbn/data_add/sbn_nd/poms_production/mc/MCP2025Av3/v10_04_06_01/prodgenie_corsika_proton_rockbox_sbnd/CV/calibntuples/3d/");
    std::cout << input_dir + "/*.root" << "\n";

    TChain* chain = new TChain("caloskim/TrackCaloSkim");
    auto dir = gSystem->OpenDirectory(input_dir);
    Int_t nfiles = 0;
    while (auto f = gSystem->GetDirEntry(dir)) { 
        std::cout << f << "\n";
        if (!strcmp(f, ".") || !strcmp(f, "..")) continue;
        nfiles += chain->Add(input_dir + f + "/*.root");
    }
    gSystem->FreeDirectory(dir);
    std::cout << nfiles << "\n";
}
