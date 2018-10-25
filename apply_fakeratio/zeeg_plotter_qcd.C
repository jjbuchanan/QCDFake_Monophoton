#include <fstream>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TColor.h"
#include "boost/format.hpp"

void plot(string histname_string, Double_t leg_xoffset, Double_t leg_yoffset, TString xaxis_title, TString plotname)
{
  TString histname = TString(histname_string+"_4");
  TString histname_qcd_sidebandUp = TString(histname_string+"_5");
  TString histname_qcd_sidebandDown = TString(histname_string+"_6");
  TString histname_qcd_METUp = TString(histname_string+"_7");
  TString histname_qcd_METDown = TString(histname_string+"_8");
  TString histname_qcd_binningUp = TString(histname_string+"_9");
  TString histname_qcd_binningDown = TString(histname_string+"_10");
  TString histname_qcd_sieieLeft = TString(histname_string+"_11");
  TString histname_qcd_sieieRight = TString(histname_string+"_12");
  TString histname_qcd_templateUp = TString(histname_string+"_13");
  TString histname_qcd_templateDown = TString(histname_string+"_14");
  TString histname_qcd_unweighted = TString(histname_string+"_15");

  std::vector<TH1F*> histo_vector;
  histo_vector.clear();
  
  Float_t int_lumi = 35900.0;
  Float_t scale_factor = 0.984*1.002; // Flat pixel seed veto SF times flat pho ID SF
  
  double photon_scale_factor_unc = sqrt(pow(0.009/0.984,2)+pow(0.007/1.002,2)); // combining pix seed and pho ID components
  double electron_scale_factor_unc = 0.04; // uncertainty on the square of the electron scale factor
  
  double total_background = 0.0;
  double background_unc_sumsquares = 0.0;
  
  TCanvas *c = new TCanvas("c", "canvas",700,640);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  float t_m = 0.08; //top margin
  float b_m = 0.4; //botton margin
  float l_m = 0.09; //left margin
  float r_m = 0.05; //right margin
  c->SetTopMargin(t_m);
  c->SetBottomMargin(b_m);
  c->SetLeftMargin(l_m);
  c->SetRightMargin(r_m);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->cd();
      
  TFile *f_data = new TFile("ZeeG_data_all.root");
  TH1F* histo_data = (TH1F*)((TH1F*)f_data->Get(histname))->Clone("data_obs");
  const int nBins = histo_data->GetXaxis()->GetNbins();
  histo_data->SetBinContent(nBins, histo_data->GetBinContent(nBins)+histo_data->GetBinContent(nBins+1));
  histo_data->ClearUnderflowAndOverflow();
  Float_t int_data = histo_data->Integral();
  histo_data->SetLineWidth(3);
  histo_data->SetLineColor(kWhite);
  histo_data->SetMarkerStyle(kFullSquare);
  histo_data->SetMarkerColor(kWhite);
  histo_vector.push_back(histo_data);
  
  //Now that nBins has been specified, initialize binned systematic shift repositories to the appropriate length
  std::vector<double> syst_shiftUp_jetfake;
  syst_shiftUp_jetfake.clear();
  std::vector<double> syst_shiftDown_jetfake;
  syst_shiftDown_jetfake.clear();
  for(int i = 1; i <= nBins; i++){
    syst_shiftUp_jetfake.push_back(0);
    syst_shiftDown_jetfake.push_back(0);
  }

  TFile *f_jetfake = new TFile("ZeeG_qcd_all.root");
  TH1F* histo_jetfake = (TH1F*)((TH1F*)f_jetfake->Get(histname))->Clone("histo_jetfake");
  histo_jetfake->SetBinContent(nBins, histo_jetfake->GetBinContent(nBins)+histo_jetfake->GetBinContent(nBins+1));
  histo_jetfake->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sidebandUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sidebandUp))->Clone("histo_jetfake_sidebandUp");
  histo_jetfake_sidebandUp->SetBinContent(nBins, histo_jetfake_sidebandUp->GetBinContent(nBins)+histo_jetfake_sidebandUp->GetBinContent(nBins+1));
  histo_jetfake_sidebandUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sidebandDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sidebandDown))->Clone("histo_jetfake_sidebandDown");
  histo_jetfake_sidebandDown->SetBinContent(nBins, histo_jetfake_sidebandDown->GetBinContent(nBins)+histo_jetfake_sidebandDown->GetBinContent(nBins+1));
  histo_jetfake_sidebandDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_METUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_METUp))->Clone("histo_jetfake_METUp");
  histo_jetfake_METUp->SetBinContent(nBins, histo_jetfake_METUp->GetBinContent(nBins)+histo_jetfake_METUp->GetBinContent(nBins+1));
  histo_jetfake_METUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_METDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_METDown))->Clone("histo_jetfake_METDown");
  histo_jetfake_METDown->SetBinContent(nBins, histo_jetfake_METDown->GetBinContent(nBins)+histo_jetfake_METDown->GetBinContent(nBins+1));
  histo_jetfake_METDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_binningUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_binningUp))->Clone("histo_jetfake_binningUp");
  histo_jetfake_binningUp->SetBinContent(nBins, histo_jetfake_binningUp->GetBinContent(nBins)+histo_jetfake_binningUp->GetBinContent(nBins+1));
  histo_jetfake_binningUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_binningDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_binningDown))->Clone("histo_jetfake_binningDown");
  histo_jetfake_binningDown->SetBinContent(nBins, histo_jetfake_binningDown->GetBinContent(nBins)+histo_jetfake_binningDown->GetBinContent(nBins+1));
  histo_jetfake_binningDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sieieLeft = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sieieLeft))->Clone("histo_jetfake_sieieLeft");
  histo_jetfake_sieieLeft->SetBinContent(nBins, histo_jetfake_sieieLeft->GetBinContent(nBins)+histo_jetfake_sieieLeft->GetBinContent(nBins+1));
  histo_jetfake_sieieLeft->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sieieRight = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sieieRight))->Clone("histo_jetfake_sieieRight");
  histo_jetfake_sieieRight->SetBinContent(nBins, histo_jetfake_sieieRight->GetBinContent(nBins)+histo_jetfake_sieieRight->GetBinContent(nBins+1));
  histo_jetfake_sieieRight->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_templateUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_templateUp))->Clone("histo_jetfake_templateUp");
  histo_jetfake_templateUp->SetBinContent(nBins, histo_jetfake_templateUp->GetBinContent(nBins)+histo_jetfake_templateUp->GetBinContent(nBins+1));
  histo_jetfake_templateUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_templateDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_templateDown))->Clone("histo_jetfake_templateDown");
  histo_jetfake_templateDown->SetBinContent(nBins, histo_jetfake_templateDown->GetBinContent(nBins)+histo_jetfake_templateDown->GetBinContent(nBins+1));
  histo_jetfake_templateDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_unweighted = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_unweighted))->Clone("histo_jetfake_unweighted");
  histo_jetfake_unweighted->SetBinContent(nBins, histo_jetfake_unweighted->GetBinContent(nBins)+histo_jetfake_unweighted->GetBinContent(nBins+1));
  histo_jetfake_unweighted->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_errUp = (TH1F*)histo_jetfake->Clone("histo_jetfake_errUp");
  TH1F* histo_jetfake_errDown = (TH1F*)histo_jetfake->Clone("histo_jetfake_errDown");
  Float_t int_jetfake = histo_jetfake->Integral();
  Float_t max_int_jetfake = 0.0;
  Float_t min_int_jetfake = 0.0;
  Float_t stat_jetfake = 0.0;
  for(int i = 1; i <= nBins; i++){
    double int_bin_jetfake = histo_jetfake->GetBinContent(i);
    double int_bin_jetfake_sidebandUp = histo_jetfake_sidebandUp->GetBinContent(i);
    double int_bin_jetfake_sidebandDown = histo_jetfake_sidebandDown->GetBinContent(i);
    double int_bin_jetfake_METUp = histo_jetfake_METUp->GetBinContent(i);
    double int_bin_jetfake_METDown = histo_jetfake_METDown->GetBinContent(i);
    double int_bin_jetfake_binningUp = histo_jetfake_binningUp->GetBinContent(i);
    double int_bin_jetfake_binningDown = histo_jetfake_binningDown->GetBinContent(i);
    double int_bin_jetfake_sieieLeft = histo_jetfake_sieieLeft->GetBinContent(i);
    double int_bin_jetfake_sieieRight = histo_jetfake_sieieRight->GetBinContent(i);
    double int_bin_jetfake_templateUp = histo_jetfake_templateUp->GetBinContent(i);
    double int_bin_jetfake_templateDown = histo_jetfake_templateDown->GetBinContent(i);
    double int_bin_jetfake_unweighted = histo_jetfake_unweighted->GetBinContent(i);
    double ints_bin[] = {int_bin_jetfake, int_bin_jetfake_sidebandUp, int_bin_jetfake_sidebandDown, int_bin_jetfake_METUp, int_bin_jetfake_METDown, int_bin_jetfake_binningUp, int_bin_jetfake_binningDown, int_bin_jetfake_sieieLeft, int_bin_jetfake_sieieRight, int_bin_jetfake_templateUp, int_bin_jetfake_templateDown};
    double max_int_bin = *max_element(ints_bin, ints_bin+11);
    double min_int_bin = *min_element(ints_bin, ints_bin+11);
    histo_jetfake_errUp->SetBinContent(i, max_int_bin);
    histo_jetfake_errDown->SetBinContent(i, min_int_bin);
    max_int_jetfake += max_int_bin;
    min_int_jetfake += min_int_bin;
    syst_shiftUp_jetfake[i-1] = max_int_bin-int_bin_jetfake;
    syst_shiftDown_jetfake[i-1] = min_int_bin-int_bin_jetfake;
    double stat_bin_jetfake = 0.0;
    if (int_bin_jetfake_unweighted > 0)
      stat_bin_jetfake = int_bin_jetfake/sqrt(int_bin_jetfake_unweighted);
    histo_jetfake->SetBinError(i, stat_bin_jetfake);
    stat_jetfake += stat_bin_jetfake*stat_bin_jetfake;
  }
  Float_t syst_jetfake = TMath::Max(max_int_jetfake-int_jetfake, int_jetfake-min_int_jetfake);
  stat_jetfake = sqrt(stat_jetfake);
  Float_t err_jetfake = sqrt(syst_jetfake*syst_jetfake + stat_jetfake*stat_jetfake);
  total_background += int_jetfake;
  background_unc_sumsquares += err_jetfake*err_jetfake;
  histo_jetfake->SetFillColor(kBlue-4);
  histo_vector.push_back(histo_jetfake);
    
  // Print bin contents
  cout<<endl;
  if (histname=="Photon_Et_range_4"){
    vector<float> total_background_binned;
    total_background_binned.clear();
    for(int i = 1; i <= nBins; i ++){
      total_background_binned.push_back(0.0);
    }
    
    cout<<"$E_{T}^{\\gamma}$ &        [ 175,  200] &        [ 200,  250] &        [ 250,  300] &        [ 300,  400] &        [ 400,  600] &        [ 600, 1000] \\\\"<<endl;
    
    cout<<"\\hline"<<endl;
    
    cout<<"        jetfake ";
    for(int i = 1; i <= nBins; i++){
      cout<<"& $ "<<boost::format("%.2f")%histo_jetfake->GetBinContent(i)<<" $ ";
      total_background_binned[i-1] += histo_jetfake->GetBinContent(i);
    }
    cout<<"\\\\"<<endl;
    
    cout<<"\\hline"<<endl;
  }
  
  if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4" || histname == "h_phoRecoilMt_4")
  {
    TFile* f_ZeeG_histos;
    if(histname == "Photon_Et_range_4")
      f_ZeeG_histos = new TFile("ZeeG_histos_Pt.root","RECREATE");
    else if(histname == "h_photonic_recoil_4")
      f_ZeeG_histos = new TFile("ZeeG_histos_MET.root","RECREATE");
    else if(histname == "h_phoRecoilMt_4")
      f_ZeeG_histos = new TFile("ZeeG_histos_Mt.root","RECREATE");
    f_ZeeG_histos->cd();
    histo_jetfake->Write();
    histo_jetfake_errUp->Write();
    histo_jetfake_errDown->Write();
    f_ZeeG_histos->Close();
  }
  
  if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4"){
    for(int i = 1; i <= nBins; i++){
      double binWidth = histo_data->GetBinWidth(i);
      histo_data->SetBinContent(i,histo_data->GetBinContent(i)/binWidth);
      histo_data->SetBinError(i,histo_data->GetBinError(i)/binWidth);
      histo_jetfake->SetBinContent(i,histo_jetfake->GetBinContent(i)/binWidth);
      histo_jetfake->SetBinError(i,histo_jetfake->GetBinError(i)/binWidth);
    }
  }
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.26,0.99,0.99);
  pad1->Draw(); pad1->cd();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  TH1F *histo_allbackgrounds = (TH1F*)histo_jetfake->Clone("histo_allbackgrounds");
    
  for(int i = 1; i <= nBins; i++){
    double background = histo_allbackgrounds->GetBinContent(i);
    // Add statistical errors
    double sum_binerrors_squared = 0.0;
    sum_binerrors_squared += pow(histo_jetfake->GetBinError(i),2);
    double binerror = sqrt(sum_binerrors_squared); // Include just the statistical error
    double jetfakeerr = (fabs(syst_shiftUp_jetfake[i-1])+fabs(syst_shiftDown_jetfake[i-1]))/2.0;
    if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4"){
      double binWidth = histo_data->GetBinWidth(i);
      jetfakeerr /= binWidth;
    }
    binerror = sqrt(sum_binerrors_squared+pow(background*photon_scale_factor_unc,2)+pow(background*electron_scale_factor_unc,2)+pow(jetfakeerr,2));
    histo_allbackgrounds->SetBinError(i,binerror);
  }
  histo_allbackgrounds->SetFillColorAlpha(kGray+1,0.6);
  histo_vector.push_back(histo_allbackgrounds);
  
  if (histname == "Photon_Et_range_4"){
    for(int i = 1; i <= nBins; i++){
      float data = histo_data->GetBinContent(i);
      float binWidth = histo_data->GetBinWidth(i);
      cout<<"data, bin "<<i<<": "<<(data*binWidth)<<endl;
    }
    for(int i = 1; i <= nBins; i++){
      float background = histo_allbackgrounds->GetBinContent(i);
      float background_err = histo_allbackgrounds->GetBinError(i);
      float binWidth = histo_allbackgrounds->GetBinWidth(i);
      cout<<"background, bin "<<i<<": "<<(background*binWidth)<<" +/- "<<(background_err*binWidth)<<endl;
    }
  }

  TH1F *histo_allbackgrounds_outline = (TH1F*)histo_allbackgrounds->Clone("histo_allbackgrounds_outline");
  histo_allbackgrounds_outline->SetFillColorAlpha(kWhite,0.0);
  histo_allbackgrounds_outline->SetLineWidth(1);
  histo_vector.push_back(histo_allbackgrounds_outline);
    
  THStack *stackHisto = new THStack("stackHisto","Title");
  stackHisto->Add(histo_jetfake);
  stackHisto->SetTitle("");
  
  for(int i = 0; i < int(histo_vector.size()); i++){
    histo_vector[i]->SetStats(0);
    histo_vector[i]->SetTitle("");
    histo_vector[i]->SetLineColor(kBlack);
    histo_vector[i]->GetXaxis()->SetTitle(xaxis_title);
    histo_vector[i]->GetXaxis()->SetLabelFont(42);
    histo_vector[i]->GetXaxis()->SetLabelSize(0.06);
    histo_vector[i]->GetXaxis()->SetTitleFont(42);
    histo_vector[i]->GetXaxis()->SetTitleSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitle("Events / bin");
    if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4")
      histo_vector[i]->GetYaxis()->SetTitle("Events / GeV");
    histo_vector[i]->GetYaxis()->SetLabelFont(42);
    histo_vector[i]->GetYaxis()->SetLabelSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitleFont(42);
    histo_vector[i]->GetYaxis()->SetTitleSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitleOffset(0.9);
  }
  
  //Accommodate both the data and background plots
  double ymax_data = 0.0;
  double ymax_background = 0.0;
  for(int i = 1; i <= nBins; i++){
    double y_data = histo_data->GetBinContent(i);
    double y_error_data = histo_data->GetBinError(i);
    double y_high_data = y_data+y_error_data;
    if(y_high_data > ymax_data)
      ymax_data = y_high_data;
    double y_background = histo_allbackgrounds->GetBinContent(i);
    double y_error_background = histo_allbackgrounds->GetBinError(i);
    double y_high_background = y_background+y_error_background;
    if(y_high_background > ymax_background)
      ymax_background = y_high_background;
  }
  
  double ymin = 0.0003;
  double ymax = 1.3*ymax_data;
  if(ymax_background > ymax_data)
    ymax = 1.3*ymax_background;
  if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4"){
    pad1->SetLogy();
    // ymax *= 5;
    ymax = 5;
  }
  // else if (histname == "nJet_4"){
  //   ymax = 2.0;
  // }
  histo_data->GetYaxis()->SetRangeUser(ymin,ymax);
  if(histname == "nJet_4")
    histo_data->GetXaxis()->SetRangeUser(0,10);
  histo_data->Draw();
  stackHisto->Draw("HIST SAME");
  histo_allbackgrounds->Draw("E2 SAME");
  histo_allbackgrounds_outline->Draw("HIST SAME");
  histo_data->SetLineColor(kBlack);
  histo_data->SetMarkerColor(kBlack);
  histo_data->Draw("E0 P0 SAME");
  gPad->RedrawAxis();
  
  //Central location of leg defined to be location of leg in phoPt plot
  TLegend* leg = new TLegend(0.5+leg_xoffset,0.58075+leg_yoffset,0.885387+leg_xoffset,0.862969+leg_yoffset,"");
  leg->AddEntry(histo_data,"Data");
  leg->AddEntry(histo_jetfake,"jet#rightarrow#gamma MisID","F");
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->Draw();
  
  float lumiTextSize = 0.6;
  float lumiTextOffset = 0.2;
  float cmsTextSize = 0.75;
  TLatex *texS = new TLatex(0.60023,0.917173,"35.9 fb^{-1} (13 TeV)");
  texS->SetNDC();
  texS->SetTextFont(42);
  texS->SetTextSize(lumiTextSize*t_m);
  texS->Draw();
  TLatex *texS1 = new TLatex(0.13592,0.817173,"#bf{CMS} #it{Preliminary}");
  texS1->SetNDC();
  texS1->SetTextFont(42);
  texS1->SetTextSize(cmsTextSize*t_m);
  texS1->Draw();
  
  c->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.26);
  pad2->Draw(); pad2->cd();
  pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.35);
  
  double max_ratio = 3.5;
  
  TH1F* Ratio = (TH1F*)histo_data->Clone("Ratio");
  TH1F* Ratio_background = (TH1F*)histo_allbackgrounds->Clone("Ratio_background");
  for(int i = 1; i <= nBins; i++){
    double y_data = histo_data->GetBinContent(i);
    double y_error_data = histo_data->GetBinError(i);
    double y_background = histo_allbackgrounds->GetBinContent(i);
    double y_error_background = histo_allbackgrounds->GetBinError(i);
    double Ratiocontent = 0.0;
    double Ratioerror = max_ratio;
    double Ratioerror_background = max_ratio;
    if(y_background > 0.){
      Ratiocontent = y_data/y_background;
      Ratioerror_background = y_error_background/y_background;
      if(y_error_data > 0.)
        Ratioerror = y_error_data/y_background;
    }
    else if(y_data > 0.){
      Ratiocontent = 3.*max_ratio;
    }
    Ratio->SetBinContent(i,Ratiocontent);
    Ratio->SetBinError(i,Ratioerror);
    Ratio_background->SetBinContent(i,1);
    Ratio_background->SetBinError(i,Ratioerror_background);
  }
  
  Ratio_background->GetYaxis()->SetRangeUser(0.0,max_ratio-0.01);
  Ratio_background->GetYaxis()->SetTitle("Data/SM");
  Ratio_background->GetYaxis()->CenterTitle();
  Ratio_background->GetYaxis()->SetLabelSize(0.14);
  Ratio_background->GetYaxis()->SetTitleSize(0.15);
  Ratio_background->GetYaxis()->SetLabelFont(42);
  Ratio_background->GetYaxis()->SetTitleFont(42);
  Ratio_background->GetYaxis()->SetTitleOffset(0.30);
  Ratio_background->GetYaxis()->SetNdivisions(305);
  Ratio_background->GetXaxis()->SetTitle(xaxis_title);
  Ratio_background->GetXaxis()->SetLabelSize(0.16);
  Ratio_background->GetXaxis()->SetTitleSize(0.18);
  Ratio_background->GetXaxis()->SetLabelFont(42);
  Ratio_background->GetXaxis()->SetTitleFont(42);
  Ratio_background->GetXaxis()->SetTitleOffset(0.9);
  Ratio_background->GetXaxis()->SetTickLength(0.05);
  Ratio_background->SetStats(0);
  Ratio->SetMarkerStyle(0);
  double xmin = histo_data->GetXaxis()->GetBinLowEdge(1);
  double xmax = histo_data->GetXaxis()->GetBinUpEdge(nBins);
  if(histname == "nJet_4"){
    Ratio_background->GetXaxis()->SetRangeUser(0,10);
    xmax = 10.0;
  }
  TLine* line = new TLine(xmin,1.,xmax,1.);
  line->SetLineStyle(2);
  line->SetLineColor(kBlack);
  gStyle->SetLineStyleString(11,"3 12");
  TLine* line0 = new TLine(xmin,0.5,xmax,0.5);
  line0->SetLineStyle(11);
  line0->SetLineColor(kBlack);
  Ratio_background->Draw("E2");
  line->Draw("SAME");
  line0->Draw("SAME");
  for(int i = 1; i <= (2*max_ratio-3); i++){
    double y_coord = 1.0 + 0.5*i;
    TLine* line_i = new TLine(xmin,y_coord,xmax,y_coord);
    line_i->SetLineStyle(11);
    line_i->SetLineColor(kBlack);
    line_i->Draw("SAME");
  }
  Ratio->Draw("E0 P0 SAME");

  double background_unc = sqrt(background_unc_sumsquares);
  
  if(histname == "Photon_Et_range_4"){
    cout<<"ZeeG region"<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<"Jet faking photon: "<<int_jetfake<<" +- "<<err_jetfake<<endl;
    cout<<"Data: "<<int_data<<endl;
    cout<<"------------------------------------"<<endl;
  }

  // c->SaveAs(TString("zeeg_qcd_"+plotname+".png"));
  c->SaveAs(TString("zeeg_qcd_"+plotname+".pdf"));
  delete(c);
}

void zeeg_plotter_qcd()
{
  std::vector<string> histnames;
  histnames.clear();
  std::vector<Double_t> leg_xoffsets;
  leg_xoffsets.clear();
  std::vector<Double_t> leg_yoffsets;
  leg_yoffsets.clear();
  std::vector<TString> xaxis_titles;
  xaxis_titles.clear();
  std::vector<TString> plotnames;
  plotnames.clear();

//  histnames.push_back(TString("_4"));
//  leg_xoffsets.push_back(0.);
//  leg_yoffsets.push_back(0.);
//  xaxis_titles.push_back(TString(""));
//  plotnames.push_back(TString(""));

  histnames.push_back("Photon_Et_range");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("phoPt"));

  histnames.push_back("Photon_SCeta");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #eta"));
  plotnames.push_back(TString("phoEta"));

  histnames.push_back("Photon_SCphi");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #phi"));
  plotnames.push_back(TString("phoPhi"));

  histnames.push_back("pfMET");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("pfMET [GeV]"));
  plotnames.push_back(TString("pfMET"));

  histnames.push_back("nJet");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Number of Jets"));
  plotnames.push_back(TString("nJet"));

  histnames.push_back("h_photonic_recoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photonic Recoil [GeV]"));
  plotnames.push_back(TString("recoil"));

  histnames.push_back("h_dPhi_phoRecoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("#Delta#phi(Photon,Recoil)"));
  plotnames.push_back(TString("dPhiPhoRecoil"));

  histnames.push_back("h_leadingLeptonPt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("elePtLeading"));

  histnames.push_back("h_leadingLeptonEta");
  leg_xoffsets.push_back(0.15);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #eta"));
  plotnames.push_back(TString("eleEtaLeading"));

  histnames.push_back("h_leadingLeptonPhi");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #phi"));
  plotnames.push_back(TString("elePhiLeading"));

  histnames.push_back("h_subleadingLeptonPt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("elePtSubleading"));

  histnames.push_back("h_subleadingLeptonEta");
  leg_xoffsets.push_back(0.15);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #eta"));
  plotnames.push_back(TString("eleEtaSubleading"));

  histnames.push_back("h_subleadingLeptonPhi");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #phi"));
  plotnames.push_back(TString("elePhiSubleading"));
  
  histnames.push_back("h_dileptonPt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("dielePt"));
  
  histnames.push_back("h_dileptonM");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron Invariant Mass [GeV]"));
  plotnames.push_back(TString("dieleM"));

  histnames.push_back("h_phoPT_over_photonicRecoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} / Recoil"));
  plotnames.push_back(TString("phoPtOverRecoil"));

  histnames.push_back("h_dileptonPt_over_pfMET");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron #it{p}_{T} / pfMET"));
  plotnames.push_back(TString("dielePtOverpfMET"));

  histnames.push_back("h_min_dphijetmet");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Min. #Delta#phi(jets,MET)"));
  plotnames.push_back(TString("dPhiJetsMET"));

  histnames.push_back("h_min_dphijetrecoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Min. #Delta#phi(jets,Recoil)"));
  plotnames.push_back(TString("dphiJetsRecoil"));

  histnames.push_back("h_dPhi_leptons");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("#Delta#phi(Leading Ele,Subleading Ele)"));
  plotnames.push_back(TString("dPhiLeadingEleSubleadingEle"));

  histnames.push_back("h_phoRecoilMt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon-Recoil #it{M}_{T} [GeV]"));
  plotnames.push_back(TString("phoRecoilmT"));
  
  for(int i = 0; i < histnames.size(); i++){
    plot(histnames[i],leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
}