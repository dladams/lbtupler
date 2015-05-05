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

int draw(std::string name, int how =0, double xmin =0.0, double xmax =0.0) {
  const string myname = "draw: ";
  double xh1 = 0.05;
  double xh2 = 0.94;
  TObject* pobj = 0;
  gDirectory->GetObject(name.c_str(), pobj);
  if ( pobj == 0 ) {
    cout << myname << "Object not found: " << name << endl;
    return 1;
  }
  static TCanvas* pcan = 0;
  string cname = "mycan";
  TH2* ph2 = dynamic_cast<TH2*>(pobj);
  if ( ph2 == 0 ) {
    cout << myname << "Object is not TH2: " << name << endl;
    return 2;
  }
  if ( ph2->GetEntries() == 0 ) {
    cout << myname << "Histogram is empty: " << name << endl;
    return 3;
  }
  // Early draw to make palette available
  double zmin = ph2->GetMinimum();
  cout << "zmin: " << zmin << endl;
  if ( zmin < 0.0 ) palette(2);
  else palette(1);
  // Draw canvas and determine set palette parameters.
  double palx1 = 0.865;
  double palx2 = 0.91;
  double paltoff = 1.00;
  double axyoff = 1.0;
  double toff = 1.00;
  if ( how == 1 ) {
    // Extra-wide canvas.
    pcan = new TCanvas(cname.c_str(), cname.c_str(), 1600, 500);
    pcan->SetLeftMargin(xh1);
    pcan->SetRightMargin(1.0-xh2);
    palx1 = xh2;
    palx2 = xh2 + 0.015;
    paltoff = 0.50;
    axyoff = 0.50;
  } else if ( how == 2 ) {
    // Standard canvas.
    pcan = new TCanvas;
  }
  // Set palette parameters. First draw to make sure palette is present.
  ph2->Draw("colz");
  gPad->Update();
  TPaletteAxis* ppalax = dynamic_cast<TPaletteAxis*>(ph2->GetListOfFunctions()->FindObject("palette")); 
  if ( ppalax == 0 ) {
    cout << myname << "Unable to retrieve palette." << endl;
    return 4;
  }
  cout << "Palette axis: " << hex << long(ppalax) << dec << endl;
  ppalax->SetX1NDC(palx1);
  ppalax->SetX2NDC(palx2);
  ppalax->SetTitleOffset(paltoff);
  // Set axis parameters.
  //ph2->SetTitleOffset(toff);   // Sets x-axis label!!
  ph2->SetTickLength(0.010, "Y");
  ph2->SetTickLength(0.010, "Z");
  ph2->GetYaxis()->SetTitleOffset(axyoff);
  // Draw.
  ph2->Draw("colz");
  //ph2->Draw("axis same");
  // Retrieve the palette axis.
  gPad->Update();
  ph2->GetListOfFunctions()->Print(); 
  if ( xmax > xmin ) ph2->GetXaxis()->SetRangeUser(xmin, xmax);
  ph2->Draw("colz");
  gPad->Update();
  return 0;
}
