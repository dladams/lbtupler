{
  TFile* pfile = TFile::Open("LbTupler.root");
  LbTupler->cd();
  //gROOT->ProcessLine(".L palette.C");
  gROOT->ProcessLine(".L palette.cxx+");
  gROOT->ProcessLine(".L gettrees.cxx+");
  gROOT->ProcessLine(".L draw.cxx+");
  gROOT->ProcessLine(".L drawTracks.C");
  gStyle->SetPadRightMargin(0.14);   // For 2D plots
  if ( simtree("LbTuplerSimulation") ) {
    simtree()->SetMarkerStyle(2);
  }
}
