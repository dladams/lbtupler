void dump() {
  string fname = "/home/dladams/lbne/data/prodcosmics_lbne35t_milliblock_1_20150220T201346_gen_20150220T213621_g4_20150220T234543_detsim_milliblock_20150223T214501_reco_milliblock.root";
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
