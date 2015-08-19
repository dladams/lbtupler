// det35t.C

class TpcBuilder {
public:
  double width;
  double height;
  double length;
  int ncopy = 0;
  TpcBuilder(double awidth, double aheight, double alength)
  : width(awidth), height(aheight), length(alength) { }
  void addTpc(int itpc, double xc, double yc, double zc) const {
    ++ncopy;
    ostringstream ssname;
    ssname << "TPC" << itpc;
    TGeoVolume* ptop = gGeoManager->GetTopVolume();
    TGeoMedium* pmed = gGeoManager->GetMedium("LAr");
    TGeoVolume* pvol = gGeoManager->MakeBox(ssname.str().c_str(), pmed, 0.5*width, 0.5*height, 0.5*length);
    TGeoTranslation* ptrans = new TGeoTranslation(xc, yc, zc);
    pvol->SetTransparency(75);
    pvol->SetLineColor(33);
    ptop->AddNode(pvol, ncopy, ptrans);
  }
};
TEveGeoTopNode* en = 0;

void det35t() {
  int inode = 0;
  TGeoManager* pgeo = new TGeoManager("geolardet","LAr detector");
  TGeoMaterial* pmatVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial* pmatAl = new TGeoMaterial("Al", 26.98,13,2.7);
  TGeoMaterial* pmatLAr = new TGeoMaterial("LAr", 38.95,18,1.4);

  TGeoMedium* pvacuum = new TGeoMedium("Vacuum",1, pmatVacuum);
  TGeoMedium* pAl = new TGeoMedium("Aluminium",2, pmatAl);
  TGeoMedium* pLAr = new TGeoMedium("LAr" , 3, pmatLAr);

  cout << "       pgeo: " << pgeo << endl;
  cout << "gGeoManager: " << gGeoManager << endl;
  cout << gGeoManager->GetMedium("LAr")->Print() << endl;

  //--- make the top container volume
  TGeoVolume* ptop = pgeo->MakeBox("TOP", pvacuum, 1000., 1000., 1000.);
  ptop->SetInvisible();
  pgeo->SetTopVolume(ptop);
  //pgeo->SetVisDensity(0.2);

  //TGeoVolume* pdet = new TGeoVolumeAssembly("DET");
  TpcBuilder tbld1( 28.6391, 204.564,  53.4466);
  TpcBuilder tbld2(223.939,  204.564,  53.4466);
  TpcBuilder tbld3( 28.6391,  92.8743, 51.924);
  TpcBuilder tbld4(223.939,   92.8743, 51.924);
  TpcBuilder tbld5( 28.6391, 120.794,  51.924);
  TpcBuilder tbld6(223.939,  120.794,  51.924);
  TpcBuilder tbld7( 28.6391, 204.564,  53.4466);
  TpcBuilder tbld8(223.939,  204.564,  53.4466);
  tbld1.addTpc(1, -20.8648,   18.5122,  24.6852);
  tbld2.addTpc(2, 110.49  ,   18.5122,  24.6852);
  tbld3.addTpc(3, -20.8648,  -46.4372,  77.3705);
  tbld4.addTpc(4, 110.49  ,  -46.4372,  77.3705);
  tbld5.addTpc(5, -20.8648,   60.3972,  77.3705);
  tbld6.addTpc(6, 110.49  ,   60.3972,  77.3705);
  tbld7.addTpc(7, -20.8648,   18.5122, 130.056);
  tbld8.addTpc(8, 110.49  ,   18.5122, 130.056);

  pgeo->CloseGeometry();
  TEveManager::Create();
  TGeoNode* node = gGeoManager->GetTopNode();
  en = new TEveGeoTopNode(gGeoManager, node);
  en->SetVisLevel(4);
  en->GetNode()->GetVolume()->SetVisibility(kFALSE);
  gEve->AddGlobalElement(en);
  gEve->Redraw3D(kTRUE);

   en->ExpandIntoListTreesRecursively();
   en->Save("det35t.png", "det35t");

}
