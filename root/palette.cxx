// palette.cxx

#include "palette.h"
#include <iostream>
#include "TColor.h"

using std::cout;
using std::endl;

int palette(int ipal) {
  static Int_t  colors[50];
  static Bool_t initialized = kFALSE;
  double alpha = 1.0;
  int colout = 0;
  if ( ipal == 1 ) {
    // white->orange->red->black
    const int nRGBs = 6;
    Double_t stops[nRGBs] =       { 0.00, 0.06, 0.12, 0.24, 0.60, 1.00};
    Double_t red[nRGBs]   =       { 1.00, 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] =       { 1.00, 1.00, 0.75, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  =       { 1.00, 1.00, 0.00, 0.00, 0.10, 0.00};
    colout = 20;
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, alpha);
    return colout;
  } else if ( ipal == 2 ) {
    // blue->white->yellow->red->black
    const int nRGBs = 8;
    Double_t stops[nRGBs] = { 0.00, 0.48, 0.50, 0.53, 0.56, 0.62, 0.80, 1.00};
    Double_t red[nRGBs]   = { 0.09, 0.75, 1.00, 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] = { 0.60, 0.80, 1.00, 1.00, 0.75, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  = { 0.48, 0.93, 1.00, 1.00, 0.00, 0.00, 0.10, 0.00};
    colout = 40;
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, alpha);
    cout << TColor::GetColorPalette(0) << endl;
    return colout;
  } else if ( ipal == 11 ) {
    // white->yellow->red->black
    const int nRGBs = 6;
    Double_t stops[nRGBs] =       { 0.00, 0.06, 0.24, 0.44, 0.76, 1.00};
    Double_t red[nRGBs]   =       { 1.00, 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] =       { 1.00, 1.00, 1.00, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  =       { 1.00, 1.00, 0.20, 0.00, 0.10, 0.00};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, alpha);
    return colout;
  } else if ( ipal == 12 ) {
    // blue->white->yellow->red->black
    const int nRGBs = 7;
    Double_t stops[nRGBs] = { 0.00, 0.50, 0.53, 0.62, 0.72, 0.88, 1.00};
    Double_t red[nRGBs]   = { 0.09, 1.00, 1.00, 1.00, 1.00, 0.70, 0.00};
    Double_t green[nRGBs] = { 0.68, 1.00, 1.00, 1.00, 0.55, 0.20, 0.00};
    Double_t blue[nRGBs]  = { 0.79, 1.00, 1.00, 0.20, 0.00, 0.10, 0.00};
    colout = 40;
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, alpha);
    return colout;
  } else {
    // white->black
    const int nRGBs = 2;
    Double_t stops[nRGBs] = { 0.00, 1.00};
    Double_t red[nRGBs]   = { 1.00, 0.00};
    Double_t green[nRGBs] = { 1.00, 0.00};
    Double_t blue[nRGBs]  = { 1.00, 0.00};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255, alpha);
    return colout;
  }
  return colout;
}
