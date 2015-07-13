// draw.C

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "TDirectory.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPaletteAxis.h"
#include "palette.h"

using std::string;
using std::ostringstream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;

int draw(std::string name ="help", int how =0, double xmin =0.0, double xmax =0.0) {
  const string myname = "draw: ";
  if ( name == "help" ) {
    cout << myname << endl;
    cout << myname << "Usage: draw(hname, how, xmin, xmax);" << endl;
    cout << myname << "how = 0 - new standard canvas" << endl;
    cout << myname << "      1 - new stretched canvas" << endl;
    cout << myname << "     -1 - add to existing histogram/canvas" << endl;
    cout << myname << "      * - draw on current canvas" << endl;
    return 0;
  }
  double xh1 = 0.05;
  double xh2 = 0.94;
  TObject* pobj = 0;
  gDirectory->GetObject(name.c_str(), pobj);
  if ( pobj == 0 ) {
    size_t i1 = name.find('h') + 1;
    size_t i2 = name.find('_');
    if ( i1 ==1 && i2 != string::npos && i2 > i1) {
      string sevt = "event" + name.substr(i1, i2-i1);
      cout << myname << "Trying event directory " << sevt << " for " << name << endl;
      string savedir = gDirectory->GetPath();
      if ( gDirectory->cd(sevt.c_str()) ) {
        gDirectory->GetObject(name.c_str(), pobj);
      }
      gDirectory->cd(savedir.c_str());
    } else {
      cout << myname << "Name " << name << " is not in expected format: " << name << endl;
      cout << myname << "(i1=" << i1 << ", i2=" << i2 << ")" << endl;
    }
  }
  if ( pobj == 0 ) {
    cout << myname << "Object not found: " << name << endl;
    return 1;
  }
  static TCanvas* pcan = 0;
  static TH2* phdraw;
  string cname = "mycan";
  TH2* phnew = dynamic_cast<TH2*>(pobj);
  if ( phnew == 0 ) {
    cout << myname << "Object is not TH2: " << name << endl;
    return 2;
  }
  if ( phnew->GetEntries() == 0 ) {
    cout << myname << "Histogram is empty: " << name << endl;
    return 3;
  }
  // Early draw to make palette available
  double zmin = phnew->GetMinimum();
  cout << "zmin: " << zmin << endl;
  if ( zmin < 0.0 ) palette(2);
  else palette(1);
  // Draw canvas and determine set palette parameters.
  double palx1 = 0.865;
  double palx2 = 0.91;
  double paltoff = 1.00;
  double axyoff = 1.0;
  double toff = 1.00;
  string dopt = "colz";
  bool add = false;
  if ( how == 0 ) {
    // Standard canvas.
    pcan = new TCanvas;
    pcan->Print();
  } else if ( how == 1 ) {
    // Extra-wide canvas.
    pcan = new TCanvas(cname.c_str(), cname.c_str(), 1600, 500);
    pcan->SetLeftMargin(xh1);
    pcan->SetRightMargin(1.0-xh2);
    palx1 = xh2;
    palx2 = xh2 + 0.015;
    paltoff = 0.50;
    axyoff = 0.50;
  } else if ( how == -1 ) {
    if ( pcan == 0 || phdraw == 0 ) {
      cout << "There is no existing histogram!" << endl;
      return 1;
    }
    add = true;
  } else {
    if ( pcan == 0 ) {
      pcan = new TCanvas;
    } else {
      cout << "Reusing last canvas: " << pcan->GetName() << endl;
    }
  }
  static int hcount = 0;
  if ( add ) {
    cout << "Adding histogram " << phnew->GetName() << " to " << phdraw->GetName() << endl;
    // Loop over bins in the new histogram and add to corresponding bin in the
    // histogram to be displayed. Under and overflow bins are ignored.
    for ( int iy=1; iy<=phnew->GetNbinsY(); ++iy ) {
      double y = phnew->GetYaxis()->GetBinCenter(iy);
      for ( int ix=1; ix<=phnew->GetNbinsX(); ++ix ) {
        double x = phnew->GetXaxis()->GetBinCenter(ix);
        double val = phnew->GetBinContent(ix, iy);
        if ( val != 0.0 ) {
          cout << "Adding " << x << ", " << y << ": " << val << endl;
          phdraw->Fill(x, y, val);
        }
      }
    }
  } else {
    ++hcount;
    ostringstream sshname;
    sshname << "hdraw" << hcount;
    string hname = sshname.str();
    cout << "Creating histogram " << hname << endl;
    phdraw = dynamic_cast<TH2*>(phnew->Clone(hname.c_str()));
    // Set palette parameters. First draw to make sure palette is present.
    phdraw->Draw(dopt.c_str());
    gPad->Update();
    TPaletteAxis* ppalax = dynamic_cast<TPaletteAxis*>(phdraw->GetListOfFunctions()->FindObject("palette")); 
    if ( ppalax == 0 ) {
      cout << myname << "Unable to retrieve palette." << endl;
      return 4;
    }
    //cout << "Palette axis: " << hex << long(ppalax) << dec << endl;
    ppalax->SetX1NDC(palx1);
    ppalax->SetX2NDC(palx2);
    ppalax->SetTitleOffset(paltoff);
    // Set axis parameters.
    //ph2->SetTitleOffset(toff);   // Sets x-axis label!!
    phdraw->SetTickLength(0.010, "Y");
    phdraw->SetTickLength(0.010, "Z");
    phdraw->GetYaxis()->SetTitleOffset(axyoff);
  }
  // Draw.
  phdraw->Draw(dopt.c_str());
  //ph2->Draw("axis same");
  // Retrieve the palette axis.
  gPad->Update();
  //phdraw->GetListOfFunctions()->Print(); 
  if ( !add &&  xmax > xmin ) {
    phdraw->GetXaxis()->SetRangeUser(xmin, xmax);
    phdraw->Draw(dopt.c_str());
  }
  gPad->Update();
  return 0;
}
