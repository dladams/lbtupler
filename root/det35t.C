class TpcBuilder {
public:
  double width;
  double height;
  double length;
  int ncopy = 0;
  TGeoManager* pgeo;
  TGeoMedium* pmed;
  TpcBuilder(TGeoManager* apgeo, TGeoMedium* apmed, double awidth, double aheight, double alength)
  : pgeo(apgeo), pmed(apmed), width(awidth), height(aheight), length(alength) { }
  void addTpc(TGeoVolume* ptop, int itpc, double xc, double yc, double zc) const {
    ++ncopy;
    ostringstream ssname;
    ssname << "TPC" << itpc;
    TGeoVolume* pvol = gGeoManager->MakeBox(ssname.str().c_str(), pmed, 0.5*width, 0.5*height, 0.5*length);
    TGeoTranslation* ptrans = new TGeoTranslation(xc, yc, zc);
    //TGeoTranslation* ptrans = new TGeoTranslation(xc+0.5*width, yc+0.5*height, zc+0.5*length);
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

  //--- make the top container volume
  TGeoVolume* ptop = pgeo->MakeBox("TOP", pvacuum, 1000., 1000., 1000.);
  ptop->SetInvisible();
  pgeo->SetTopVolume(ptop);
  pgeo->SetVisDensity(0.2);

  //TGeoVolume* pdet = new TGeoVolumeAssembly("DET");
  TpcBuilder tbld1(pgeo, pLAr,  27.175, 196.912, 53.4466);
  TpcBuilder tbld2(pgeo, pLAr, 222.475, 196.912, 53.4466);
  TpcBuilder tbld3(pgeo, pLAr,  27.175,  85.222, 51.924);
  TpcBuilder tbld4(pgeo, pLAr, 222.475,  85.222, 51.924);
  TpcBuilder tbld5(pgeo, pLAr,  27.175, 113.142, 51.92);
  TpcBuilder tbld6(pgeo, pLAr, 222.475, 113.142, 51.924);
  TpcBuilder tbld7(pgeo, pLAr,  27.175, 196.912, 53.4466); 
  TpcBuilder tbld8(pgeo, pLAr, 222.475, 196.912, 53.4466);
  tbld1.addTpc(ptop, 1, -20.8648,   18.5122,  24.6852);
  tbld2.addTpc(ptop, 2, 110.49  ,   18.5122,  24.6852);
  tbld3.addTpc(ptop, 3, -20.8648,  -46.4372,  77.3705);
  tbld4.addTpc(ptop, 4, 110.49  ,  -46.4372,  77.3705);
  tbld5.addTpc(ptop, 5, -20.8648,   60.3972,  77.3705);
  tbld6.addTpc(ptop, 6, 110.49  ,   60.3972,  77.3705);
  tbld7.addTpc(ptop, 7, -20.8648,   18.5122, 130.056);
  tbld8.addTpc(ptop, 8, 110.49  ,   18.5122, 130.056);

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
