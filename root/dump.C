void dump() {
  //string fname = "/home/dladams/lbne/data/prodcosmics_lbne35t_milliblock_1_20150220T201346_gen_20150220T213621_g4_20150220T234543_detsim_milliblock_20150223T214501_reco_milliblock.root";
  //string fname = "/home/dladams/lbne/data/AntiMuonCutEvents_LSU_lbne35t_7_20150612T222627_gen_20150613T005121_g4_20150613T012324_detsim_20150613T015359_reco.root";
  string fname = "/home/dladams/lbne/data/prodgenie_nu_dune10kt_workspace_19_20150717T205803_gen_979df03b-1748-4ee9-97da-06c605028827_20150717T220354_g4_20150718T001808_detsim.root";
  TFile* pfile = TFile::Open(fname.c_str());
  if ( pfile == 0 || ! pfile->IsOpen() ) {
    cout << "Unable to open input file." << endl;
    return;
  }
  TObject* pobj = pfile->Get("Events");
  if ( pobj == 0 ) {
    cout << "Unable to find event tree." << endl;
    pfile->ls();
    return;
  }
  pobj->Print();
}
