{
  TFile* pfile = TFile::Open("../run/LbTupler.root");
  if ( pfile == 0 || ! pfile.IsOpen() ) {
    cout << "Input file not found." << endl;
  } else {
    LbTupler->cd();
  }
  //gROOT->ProcessLine(".L palette.C");
  gROOT->ProcessLine(".L palette.cxx+");
  gROOT->ProcessLine(".L gettrees.cxx+");
  gROOT->ProcessLine(".L draw.cxx+");
  gROOT->ProcessLine(".L drawTracks.C");
  gStyle->SetPadRightMargin(0.14);   // For 2D plots
  if ( simtree("McParticleTree") ) {
    simtree()->SetMarkerStyle(2);
  }
}
