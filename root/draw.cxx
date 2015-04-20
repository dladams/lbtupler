// draw.C

#include <string>
#include <iostream>
#include <iomanip>
#include "TDirectory.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPaletteAxis.h"
#include "palette.h"

using std::string;
using std::cout;
using std::endl;
using std::hex;
using std::dec;

int draw(std::string name, int how =0) {
  double xh1 = 0.05;
  double xh2 = 0.94;
  TObject* pobj = 0;
  gDirectory->GetObject(name.c_str(), pobj);
  if ( pobj == 0 ) {
    cout << "Object not found: " << name << endl;
    return 1;
  }
  static TCanvas* pcan = 0;
  string cname = "mycan";
  if ( how == 1 ) {
    pcan = new TCanvas(cname.c_str(), cname.c_str(), 1600, 400);
    pcan->SetLeftMargin(xh1);
    pcan->SetRightMargin(1.0-xh2);
  }
  TH2* ph2 = dynamic_cast<TH2*>(pobj);
  if ( ph2 == 0 ) {
    cout << "Object is not TH2: " << name << endl;
    return 2;
  }
  double zmin = ph2->GetMinimum();
  cout << "zmin: " << zmin << endl;
  if ( zmin < 0.0 ) palette(2);
  else palette(1);
  ph2->SetTickLength(0.010, "Y");
  ph2->SetTickLength(0.010, "Z");
  ph2->GetYaxis()->SetTitleOffset(0.50);
  ph2->Draw("colz");
  //ph2->Draw("axis same");
  // Retrieve the palette axis.
  gPad->Update();
  ph2->GetListOfFunctions()->Print(); 
  TPaletteAxis* ppalax = dynamic_cast<TPaletteAxis*>(ph2->GetListOfFunctions()->FindObject("palette")); 
  if ( ppalax != 0 ) {
    cout << "Palette axis: " << hex << long(ppalax) << dec << endl;
    ppalax->SetX1NDC(xh2);
    ppalax->SetX2NDC(xh2+0.015);
    ppalax->SetTitleOffset(0.50);
  }
  ph2->Draw("colz");
  gPad->Update();
  return 0;
}
