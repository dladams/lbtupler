void logo() {
  TCanvas* pcan = new TCanvas("lcan", "lcan", 330, 150);
  TText t1(0.01, 0.05, "D");
  TText t2(0.24, 0.05, "U");
  TText t3(0.54, 0.05, "N");
  TText t4(0.81, 0.05, "E");
  t1->SetTextSize(0.99);
  t2->SetTextSize(1.30);
  t3->SetTextSize(1.15);
  t4->SetTextSize(0.80);
  t1.Draw();
  t2.Draw();
  t3.Draw();
  t4.Draw();
  TBox box(0.68, 0.75, 0.90, 0.90);
  box->SetFillColor(0);
  box->Draw();
  pcan->Print("logo.png");
}

