{
  gSystem->AddIncludePath("-I$BOOST_INC");
  //gSystem->AddIncludePath("-DBOOST_NO_CWCHAR");
  gSystem->AddIncludePath("-I$CPP0X_INC");
  gSystem->AddIncludePath("-I$CETLIB_INC");
  gSystem->AddIncludePath("-I$FHICLCPP_INC");
  gSystem->AddIncludePath("-I$CLHEP_INC");
  gSystem->AddIncludePath("-I$LARCORE_INC");
  gSystem->AddIncludePath("-I$LBNECODE_INC");
  gSystem->AddIncludePath("-I$LBTUPLER_INC");
  gSystem->AddDynamicPath("-L$FHICLCPP_LIB -lfhiclcpp");
  gSystem->AddLinkedLibs("$CETLIB_LIB/libcetlib.so");
  gSystem->AddLinkedLibs("$FHICLCPP_LIB/libfhiclcpp.so");
  gSystem->AddLinkedLibs("$LARCORE_LIB/libGeometry.so");
  gSystem->AddLinkedLibs("$LBNECODE_LIB/liblbne_Geometry.so");
  gSystem->AddLinkedLibs("$LBTUPLER_LIB/libLbTupler.so");
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
  gROOT->ProcessLine(".L detlar.cxx+");
  gROOT->ProcessLine(".L drawTracks.C");
  gStyle->SetPadRightMargin(0.14);   // For 2D plots
  if ( simtree("McParticleTree") ) {
    simtree()->SetMarkerStyle(2);
  }
}
