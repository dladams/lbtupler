// palette.h
//
// David Adams
// ROOT macro t define the palette used in contour plots,
// e.g. plotting 23D hists with COLZ option.
// Use phist->SetContours(n) to specify n contours.
//
// Returns the number of colors for which it is optimized.

int palette(int ipal =1);
