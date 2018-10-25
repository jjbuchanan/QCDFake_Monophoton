////postAnalyzer.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom postAnalyzer analyze
//
//To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jjbuch/LatestNtuples/ /afs/hep.wisc.edu/user/jjbuchanan/private/CMSSW_7_4_9/src/output.root -1 10000
//Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
//and storing the resulting histograms in the file output.root.
//
//To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
//root[0] TFile *f = new TFile("output.root");
//root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
//root[2] efficiency->Draw("AP")
//root[3] efficiency->SetTitle("Single photon trigger efficiency")
//root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
//root[5] efficiency->Draw("AP")
//

#define postAnalyzer_cxx
#include "postAnalyzer_WlnG_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <set>
using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
  {
    std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
    return 1;
  }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
  {
    std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
    return 1;
  }
  postAnalyzer t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void postAnalyzer::Loop(Long64_t maxEvents, int reportEvery)
{
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  double nPassing;
  nPassing = 0.0;
  double nObs_bin1, nObs_bin2, nObs_bin3, nObs_bin4, nObs_bin5;
  nObs_bin1 = nObs_bin2 = nObs_bin3 = nObs_bin4 = nObs_bin5 = 0.0;
  int nMetFilters, nTrigger, nPhoCand, nMuon, nRecoil, nDPhiPhoRec, nDPhiJetMET, nEleVeto, nLepMETmT;
  nMetFilters = nTrigger = nPhoCand = nMuon = nRecoil = nDPhiPhoRec = nDPhiJetMET = nEleVeto = nLepMETmT = 0;
  double nMetFilters_weighted, nTrigger_weighted, nPhoCand_weighted, nMuon_weighted, nRecoil_weighted, nDPhiPhoRec_weighted, nDPhiJetMET_weighted, nEleVeto_weighted, nLepMETmT_weighted;
  nMetFilters_weighted = nTrigger_weighted = nPhoCand_weighted = nMuon_weighted = nRecoil_weighted = nDPhiPhoRec_weighted = nDPhiJetMET_weighted = nEleVeto_weighted = nLepMETmT_weighted = 0.0;

  std::vector<int> phoCand1;
  phoCand1.clear();

  std::vector<int> jetveto;
  jetveto.clear();

  TFile *file = new TFile("ewk_corr_all.root");
  TH1D *ewkCorrection = (TH1D*)file->Get("znng");  cout<<"ewkCorrection histo made"<<endl;

  bool debug=true;
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Coming in: "<<std::endl;
  std::cout<<"nentries:"<<nentries<<std::endl;
  //Look at up to maxEvents events, or all if maxEvents == -1.
  Long64_t nentriesToCheck = nentries;
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;

  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  TStopwatch sw;
  sw.Start();
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
  {
    
    event_.clear();
    event_info.clear();
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //=1.0 for real data
    double event_weight=1.0;
    double EWK_corrected_weight=1.0;
    double NNLO_weight=1.0;

    TLorentzVector ele_4vec, antiele_4vec;
    TLorentzVector dilepton_4vec, dilepton_4vec_temp;

    //Reconstructed event cuts
    phoCand1   = getQcdden(175,1.4442,1);
    
    if(metFilters==1536)
    {
      nMetFilters++;
      nMetFilters_weighted++;
      if(HLTPho>>12&1 == 1)
      {
        nTrigger++;
        nTrigger_weighted++;
        if(phoCand1.size() >0)
        {
          Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCand1[0]]/TMath::CosH((*phoSCEta)[phoCand1[0]]));
          event_weight =FakeRatePt(uncorrectedPhoEt,0);
          double event_weight_sidebandUp =FakeRatePt(uncorrectedPhoEt,1);
          double event_weight_sidebandDown =FakeRatePt(uncorrectedPhoEt,2);
          double event_weight_METUp =FakeRatePt(uncorrectedPhoEt,3);
          double event_weight_METDown =FakeRatePt(uncorrectedPhoEt,4);
          double event_weight_binningUp =FakeRatePt(uncorrectedPhoEt,5);
          double event_weight_binningDown =FakeRatePt(uncorrectedPhoEt,6);
          double event_weight_sieieLeft =FakeRatePt(uncorrectedPhoEt,7);
          double event_weight_sieieRight =FakeRatePt(uncorrectedPhoEt,8);
          double event_weight_templateUp =FakeRatePt(uncorrectedPhoEt,9);
          double event_weight_templateDown =FakeRatePt(uncorrectedPhoEt,10);
          nPhoCand++;
          nPhoCand_weighted+=event_weight;
          // event_weight*=(1.013 - 0.0001168*phoEt->at(phoCand1[0]));
          if(!(phoIEta->at(phoCand1[0])==-24 && phoIPhi->at(phoCand1[0])==141) && !(phoIEta->at(phoCand1[0])==4 && phoIPhi->at(phoCand1[0])==41) && !(phoIEta->at(phoCand1[0])==5 && phoIPhi->at(phoCand1[0])==41) && !(phoIEta->at(phoCand1[0])==1 && phoIPhi->at(phoCand1[0])==81) && !(phoIEta->at(phoCand1[0])==4 && phoIPhi->at(phoCand1[0])==21))
          {
            // Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCand1[0]]/TMath::CosH((*phoSCEta)[phoCand1[0]]));
            // Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrection->GetXaxis()->FindBin(uncorrectedPhoEt));
            // EWK_corrected_weight = event_weight*(1.0+.01*EWK_percent_adjustment);
            // NNLO_weight = event_weight*EWK_corrected_weight*NNLOCorrection(phoEt->at(phoCand1[0]));

            std::vector<int> mulist = muon_veto_tightID(phoCand1[0],30.0);
            std::vector<int> looseMus = muon_veto_looseID(phoCand1[0],0,10.0);
            std::vector<int> looseEles;
            looseEles.clear();
            if(mulist.size() == 1 && looseMus.size() == 1)
            {
              nMuon++;
              nMuon_weighted+=event_weight;
              looseEles = electron_veto_looseID(phoCand1[0],mulist[0],10.0);
              jetveto = JetVetoDecision(phoCand1[0],mulist[0]);
              
              TLorentzVector lep_4vec;
              // lep_4vec.SetPtEtaPhiE(elePt->at(elelist[0]),eleEta->at(elelist[0]),elePhi->at(elelist[0]),eleEn->at(elelist[0]));
              lep_4vec.SetPtEtaPhiE(muPt->at(mulist[0]),muEta->at(mulist[0]),muPhi->at(mulist[0]),muEn->at(mulist[0]));

              Double_t lepton_pt = lep_4vec.Pt();
              TLorentzVector met_4vec;
              met_4vec.SetPtEtaPhiE(pfMET,0.,pfMETPhi,pfMET);
              TLorentzVector leptoMET_4vec = lep_4vec+met_4vec;
              Double_t leptoMET = leptoMET_4vec.Pt();
              Double_t leptoMET_phi = leptoMET_4vec.Phi();

              if(leptoMET > 170)
              {
                nRecoil++;
                nRecoil_weighted+=event_weight;
                fillHistos(0,event_weight,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                if(DeltaPhi(phoPhi->at(phoCand1[0]),leptoMET_phi)>0.5)
                {
                  nDPhiPhoRec++;
                  nDPhiPhoRec_weighted+=event_weight;
                  fillHistos(1,event_weight,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                  if(dPhiJetMET_veto(jetveto))
                  {
                    nDPhiJetMET++;
                    nDPhiJetMET_weighted+=event_weight;
                    fillHistos(2,event_weight,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                    if(looseEles.size() == 0)
                    {
                      nEleVeto++;
                      nEleVeto_weighted+=event_weight;
                      fillHistos(3,event_weight,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                      double dPhi_lepMET = DeltaPhi(muPhi->at(mulist[0]),pfMETPhi);
                      double lepMET_MT = sqrt(2*muPt->at(mulist[0])*pfMET*(1-TMath::Cos(dPhi_lepMET)));
                      if(lepMET_MT < 160)
                      {
                        nLepMETmT++;
                        nLepMETmT_weighted+=event_weight;
                        if(uncorrectedPhoEt/leptoMET < 1.4){
                          fillHistos(4,event_weight,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(5,event_weight_sidebandUp,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(6,event_weight_sidebandDown,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(7,event_weight_METUp,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(8,event_weight_METDown,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(9,event_weight_binningUp,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(10,event_weight_binningDown,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(11,event_weight_sieieLeft,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(12,event_weight_sieieRight,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(13,event_weight_templateUp,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(14,event_weight_templateDown,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          fillHistos(15,1.0,phoCand1[0],jetveto,mulist,leptoMET,leptoMET_phi);
                          
                          nPassing+=event_weight;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    tree->Fill();
    
    if (jentry%reportEvery == 0)
      {
        std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
      }
  }

  if((nentriesToCheck-1)%reportEvery != 0)
    std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
  sw.Stop();
  std::cout<<"All events checked."<<std::endl;
  //Report
  std::cout << "RealTime : " << sw.RealTime() / 60.0 << " minutes" << std::endl;
    std::cout << "CPUTime  : " << sw.CpuTime()  / 60.0 << " minutes" << std::endl;
  std::cout << std::endl;
  std::cout << "Number of events inspected: " << nTotal << std::endl;
  std::cout<<std::endl;
  cout<<"nMetFilters: "<<nMetFilters<<endl;
  cout<<"nTrigger: "<<nTrigger<<endl;
  cout<<"nPhoCand: "<<nPhoCand<<endl;
  cout<<"nMuon: "<<nMuon<<endl;
  cout<<"nRecoil: "<<nRecoil<<endl;
  cout<<"nDPhiPhoRec: "<<nDPhiPhoRec<<endl;
  cout<<"nDPhiJetMET: "<<nDPhiJetMET<<endl;
  cout<<"nEleVeto: "<<nEleVeto<<endl;
  cout<<"nLepMETmT: "<<nLepMETmT<<endl;
  cout<<endl;
  cout<<"nMetFilters_weighted: "<<nMetFilters_weighted<<endl;
  cout<<"nTrigger_weighted: "<<nTrigger_weighted<<endl;
  cout<<"nPhoCand_weighted: "<<nPhoCand_weighted<<endl;
  cout<<"nMuon_weighted: "<<nMuon_weighted<<endl;
  cout<<"nRecoil_weighted: "<<nRecoil_weighted<<endl;
  cout<<"nDPhiPhoRec_weighted: "<<nDPhiPhoRec_weighted<<endl;
  cout<<"nDPhiJetMET_weighted: "<<nDPhiJetMET_weighted<<endl;
  cout<<"nEleVeto_weighted: "<<nEleVeto_weighted<<endl;
  cout<<"nLepMETmT_weighted: "<<nLepMETmT_weighted<<endl;
}

void postAnalyzer::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ADD","ADD");
  tree->Branch("event_","std::vector<unsigned int>",&event_);
  tree->Branch("event_info","std::vector<double>",&event_info);
  fileName->cd();
  
  Float_t PtBins[7]={175.,200.,250., 300., 400., 600., 1000.0};
  Float_t MetBins[7]={170.,200.,250., 300., 400., 600., 1000.0};
  Float_t MTBins[10]={0.,200.,300.,400.,500.,600.,700.,800.,1000.,1200.};
  Float_t dPhiJetMETBins[14]={0.0,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00,3.1416};

  //h_phoIEtaIPhi = new TH2F("h_phoIEtaIPhi","Photon p_{T} > 175 GeV, E^{miss}_{T} > 140 GeV",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi->Sumw2();

  // h_dPhijetmet_min = new TH1F("h_dPhijetmet_min","h_dPhijetmet_min",40,0,3.2);h_dPhijetmet_min->Sumw2();
  //       h_dPhijetmet_min_first4 = new TH1F("h_dPhijetmet_min_first4","h_dPhijetmet_min_first4",40,0,3.2);h_dPhijetmet_min_first4->Sumw2();
  // h_HT = new TH1F("h_HT","h_HT",50,0,1000);h_HT->Sumw2();
  // h_HTMET = new TH2F("h_HTMET","h_HTMET",100,0,1000,86,140,1000);h_HTMET->Sumw2();
  // h_njetMET = new TH2F("h_njetMET","h_njetMET",20,0,20,86,140,1000);h_njetMET->Sumw2();

  //Set up the histos to be filled with method fillHistos
  for(int i=0; i<20; i++)
  {
    char ptbins[100];
    sprintf(ptbins, "_%d", i);
    std::string histname(ptbins);
    h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
    h_photon_Et_range[i] = new TH1F(("Photon_Et_range"+histname).c_str(), "Photon_Et",6,PtBins);h_photon_Et_range[i]->Sumw2();
    h_photon_SCEta[i] = new TH1F(("h_photon_SCEta"+histname).c_str(), "h_photon_SCEta",15,-1.4442,1.4442);h_photon_SCEta[i]->Sumw2();
    h_photon_SCPhi[i] = new TH1F(("h_photon_SCPhi"+histname).c_str(), "h_photon_SCPhi", 15,0,3.1416);h_photon_SCPhi[i]->Sumw2();
    h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",15,0.,700.);h_pfMET[i]->Sumw2();
    h_nJet[i] = new TH1F(("nJet"+histname).c_str(), "nJet",15,0,15);h_nJet[i]->Sumw2();
    h_leadingJetPt[i] = new TH1F(("leadingJetPt"+histname).c_str(),"leadingJetPt",15,0,300);h_leadingJetPt[i]->Sumw2();
    h_leadingJetEta[i] = new TH1F(("h_leadingJetEta"+histname).c_str(),"h_leadingJetEta",15,-1.4442,1.4442);h_leadingJetEta[i]->Sumw2();
    h_PTMET[i] = new TH1F(("PTMET"+histname).c_str(),"P_{T}/Missing E_{T}",15,0,2);h_PTMET[i]->Sumw2();
    h_leptonPt[i] = new TH1F(("h_leptonPt"+histname).c_str(),"h_leptonPt",15,30.,400.);h_leptonPt[i]->Sumw2();
    h_leptonEta[i] = new TH1F(("h_leptonEta"+histname).c_str(),"h_leptonEta",15,-2.5,2.5);h_leptonEta[i]->Sumw2();
    h_leptonPhi[i] = new TH1F(("h_leptonPhi"+histname).c_str(),"h_leptonPhi",15,0.,3.1416);h_leptonPhi[i]->Sumw2();
    h_photonic_recoil[i] = new TH1F(("h_photonic_recoil"+histname).c_str(),"h_photonic_recoil",6,MetBins);h_photonic_recoil[i]->Sumw2();
    h_dPhi_phoRecoil[i] = new TH1F(("h_dPhi_phoRecoil"+histname).c_str(),"h_dPhi_phoRecoil",15,2.0,3.1416);h_dPhi_phoRecoil[i]->Sumw2();
    h_phoPT_over_photonicRecoil[i] = new TH1F(("h_phoPT_over_photonicRecoil"+histname).c_str(),"h_phoPT_over_photonicRecoil",15,0.,2.);h_phoPT_over_photonicRecoil[i]->Sumw2();
    h_leptonPt_over_pfMET[i] = new TH1F(("h_leptonPt_over_pfMET"+histname).c_str(),"h_leptonPt_over_pfMET",20,0.,3.);h_leptonPt_over_pfMET[i]->Sumw2();
    h_lepMET_MT[i] = new TH1F(("h_lepMET_MT"+histname).c_str(),"h_lepMET_MT",15,0.,400.);h_lepMET_MT[i]->Sumw2();
    h_min_dphijetmet[i] = new TH1F(("h_min_dphijetmet"+histname).c_str(),"h_min_dphijetmet",13,dPhiJetMETBins);h_min_dphijetmet[i]->Sumw2();
    h_min_dphijetrecoil[i] = new TH1F(("h_min_dphijetrecoil"+histname).c_str(),"h_min_dphijetrecoil",10,0.,3.1416);h_min_dphijetrecoil[i]->Sumw2();
    h_dPhi_lepMET[i] = new TH1F(("h_dPhi_lepMET"+histname).c_str(),"h_dPhi_lepMET",15,0.,3.1416);h_dPhi_lepMET[i]->Sumw2();
    h_dPhi_phoRecoil_fullrange[i] = new TH1F(("h_dPhi_phoRecoil_fullrange"+histname).c_str(),"h_dPhi_phoRecoil_fullrange",20,0.,3.1416);h_dPhi_phoRecoil_fullrange[i]->Sumw2();
    h_pfMETsumEt[i] = new TH1F(("pfMETsumEt"+histname).c_str(),"pfMETsumEt",6,MetBins);
    h_METoverSqrtSumEt_extended[i] = new TH1F(("METoverSqrtSumEt_extended"+histname).c_str(),"METoverSqrtSumEt",30,0,30);
    h_METoverSqrtSumEt[i] = new TH1F(("METoverSqrtSumEt"+histname).c_str(),"METoverSqrtSumEt",30,0,10);
    h_phoIEtaIPhi[i] = new TH2F(("h_phoIEtaIPhi"+histname).c_str(),"h_phoIEtaIPhi",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi[i]->Sumw2();
    h_phoRecoilMt[i] = new TH1F(("h_phoRecoilMt"+histname).c_str(),"h_phoRecoilMt",9,MTBins);
  }
}

//Fill the sequential histos at a particular spot in the sequence
void postAnalyzer::fillHistos(int histoNumber, double event_weight,int index,std::vector<int> jets,std::vector<int> leplist,Double_t leptoMET,Double_t leptoMET_phi)
{
    Float_t uncorrectedPhoEt = ((*phoSCRawE)[index]/TMath::CosH((*phoSCEta)[index]));
    h_photon_Et_range[histoNumber]->Fill(uncorrectedPhoEt,event_weight);
    h_photon_SCEta[histoNumber]->Fill(phoSCEta->at(index),event_weight);
    h_photon_SCPhi[histoNumber]->Fill(phoSCPhi->at(index),event_weight);
    h_pfMET[histoNumber]->Fill(pfMET,event_weight);
    h_PTMET[histoNumber]->Fill(uncorrectedPhoEt/pfMET,event_weight);
    h_nJet[histoNumber]->Fill(jets.size(),event_weight);
    if(jets.size()>0){
      h_leadingJetPt[histoNumber]->Fill(jetPt->at(jets[0]),event_weight);
      h_leadingJetEta[histoNumber]->Fill(jetEta->at(jets[0]),event_weight);
      int max_njets = jets.size();
      if(jets.size() > 4)
        max_njets = 4;
      double min_dphijetmet = TMath::Pi();
      double min_dphijetrecoil = TMath::Pi();
      for(int i = 0; i < max_njets; i++)
      {
        double dphijetmet = DeltaPhi(jetPhi->at(jets[i]),pfMETPhi);
        if(dphijetmet < min_dphijetmet)
          min_dphijetmet = dphijetmet;
        double dphijetrecoil = DeltaPhi(jetPhi->at(jets[i]),leptoMET_phi);
        if(dphijetrecoil < min_dphijetrecoil)
          min_dphijetrecoil = dphijetrecoil;
      }
      h_min_dphijetmet[histoNumber]->Fill(min_dphijetmet,event_weight);
      h_min_dphijetrecoil[histoNumber]->Fill(min_dphijetrecoil,event_weight);
    }
    h_leptonPt[histoNumber]->Fill(muPt->at(leplist[0]),event_weight);
    h_leptonEta[histoNumber]->Fill(muEta->at(leplist[0]),event_weight);
    h_leptonPhi[histoNumber]->Fill(muPhi->at(leplist[0]),event_weight);
    h_photonic_recoil[histoNumber]->Fill(leptoMET,event_weight);
    double dPhi_phoRecoil = DeltaPhi(phoPhi->at(index),leptoMET_phi);
    h_dPhi_phoRecoil[histoNumber]->Fill(dPhi_phoRecoil,event_weight);
    h_dPhi_phoRecoil_fullrange[histoNumber]->Fill(dPhi_phoRecoil,event_weight);
    h_phoPT_over_photonicRecoil[histoNumber]->Fill(uncorrectedPhoEt/leptoMET,event_weight);
    h_leptonPt_over_pfMET[histoNumber]->Fill(muPt->at(leplist[0])/pfMET,event_weight);
    double dPhi_lepMET = DeltaPhi(muPhi->at(leplist[0]),pfMETPhi);
    h_dPhi_lepMET[histoNumber]->Fill(dPhi_lepMET,event_weight);
    h_lepMET_MT[histoNumber]->Fill(sqrt(2*muPt->at(leplist[0])*pfMET*(1-TMath::Cos(dPhi_lepMET))),event_weight);
    h_pfMETsumEt[histoNumber]->Fill(pfMETsumEt,event_weight);
    h_METoverSqrtSumEt_extended[histoNumber]->Fill(pfMET/sqrt(pfMETsumEt),event_weight);
    h_METoverSqrtSumEt[histoNumber]->Fill(pfMET/sqrt(pfMETsumEt),event_weight);
    // h_phoIEtaIPhi[histoNumber]->Fill(phoIPhi->at(index),phoIEta->at(index),event_weight);
    h_phoRecoilMt[histoNumber]->Fill(sqrt(2*uncorrectedPhoEt*leptoMET*(1-TMath::Cos(dPhi_phoRecoil))),event_weight);
}

void postAnalyzer::scaleHistos(int histoNumber, double scale_factor)
{
  // h_photon_Et[histoNumber]->Scale(scale_factor);
  // h_photon_Et_range[histoNumber]->Scale(scale_factor);
  // h_photon_eta[histoNumber]->Scale(scale_factor);
  // h_photon_SCEta[histoNumber]->Scale(scale_factor);
  // h_photon_phi[histoNumber]->Scale(scale_factor);
  // h_photon_SCPhi[histoNumber]->Scale(scale_factor);
  // h_pfMET[histoNumber]->Scale(scale_factor);
  // h_pfMET_300[histoNumber]->Scale(scale_factor);
  // h_dPhi[histoNumber]->Scale(scale_factor);
  // h_nJet[histoNumber]->Scale(scale_factor);
  // h_PTMET[histoNumber]->Scale(scale_factor);
}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
double postAnalyzer::DeltaPhi(double phi1, double phi2)
{
  double pi = TMath::Pi();
  double dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}


//---------------------------------------------------
// get a photon candiate based on pt eta and isolation
//----------------------------------------------------

std::vector<int> postAnalyzer::getPhoCand(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++)
    {
      //      bool kinematic = (*phoEt)[p] > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      double  uncorr_pt = (phoSCRawE->at(p)/TMath::CosH(phoSCEta->at(p)));
      bool kinematic = (phoSCRawE->at(p)/TMath::CosH(phoSCEta->at(p))) > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      bool photonId = (
           ((*phoHoverE)[p]                <  0.05   ) &&
           ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0102 ) &&
           ((*phohasPixelSeed)[p]              ==  0      ) &&
           ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.37 )  &&
           ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (1.06 + (0.014 * uncorr_pt) + (0.000019 * pow(uncorr_pt, 2.0))) )  &&
                       ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (0.28 + (0.0053 * uncorr_pt)) ) );
      
      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001 && (*phoR9)[p] < 1;
      //      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9;

      if(photonId && kinematic && noncoll){
  tmpCand.push_back(p);
      }                                                                                                                                                              
    }                                                                                                                                                                

  return tmpCand;

}






std::vector<int> postAnalyzer::getPhoCand1(double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over photons
  for(int p=0;p<nPho;p++)
    {
      Float_t uncorrectedPhoEt = ((*phoSCRawE)[p]/TMath::CosH((*phoSCEta)[p]));
      bool kinematic = fabs((*phoSCEta)[p])<phoEtaCut;
      //Short_t IsoPass;
      //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
      bool photonId = (
           ((*phoHoverE)[p]                <  0.05   ) &&
           ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.0101 ) &&
           //((*phohasPixelSeed)[p]              ==  0      ) &&
           ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 1.21 )  &&
           ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (0.65 + (0.014 * uncorrectedPhoEt) + (0.000019 * pow(uncorrectedPhoEt, 2.0))) )  &&
           ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (0.18 + (0.0053 * uncorrectedPhoEt)) ) );
      
      if(photonId && kinematic){
  tmpCand.push_back(p);
      }
    }

  return tmpCand;

}


std::vector<int> postAnalyzer::getQcdden(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons
  for(int p=0;p<nPho;p++)
  {
    Float_t uncorrectedPhoEt = ((*phoSCRawE)[p]/TMath::CosH((*phoSCEta)[p]));
    //Fail loose iso
    bool passChIsoLoose = TMath::Max( ( (*phoYuPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 3.32;
    bool passNeuIsoLoose = TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (10.910 + (0.0148* uncorrectedPhoEt) + (0.000028 * pow(uncorrectedPhoEt, 2.0)));
    bool passPhoIsoLoose = TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (3.630 + (0.0053 * uncorrectedPhoEt));

    if(!passChIsoLoose || !passNeuIsoLoose || !passPhoIsoLoose)
    {
      bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //"Very loose" ID cuts with inverted shape veto
      bool photonID = (
        ((*phoSigmaIEtaIEtaFull5x5)[p]  >  0.01040 ) &&
        ((*phoHoverE)[p]                <  0.05   ) &&
        ((*phohasPixelSeed)[p]              ==  0      ) &&
        ( TMath::Max( ( (*phoYuPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(3.32) , 0.20*uncorrectedPhoEt) )  &&
        ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(10.910 + (0.0148 * uncorrectedPhoEt) + (0.000028 * pow(uncorrectedPhoEt, 2.0))) , 0.20*uncorrectedPhoEt) )  &&
        ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < TMath::Min(5.0*(3.630+ (0.0053 * uncorrectedPhoEt)) , 0.20*uncorrectedPhoEt) )
      );

      bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;

      if(kinematic && photonID && noncoll)
      {
        tmpCand.push_back(p);
      }
    }
  }

  return tmpCand;

}

// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5
Double_t postAnalyzer::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0;

  return EffectiveArea;
}

Double_t postAnalyzer::EAchargedworst(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1064;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1026;
  return EffectiveArea;
}

Double_t postAnalyzer::EAneutral(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0597;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0807;
  return EffectiveArea;
}

Double_t postAnalyzer::EAphoton(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1210;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1107;
  return EffectiveArea;
}

std::vector<int> postAnalyzer::electron_veto_tightID(int pho_index, float elePtCut)
{
  std::vector<int> ele_cands;
  ele_cands.clear();

  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
  {
    //Make sure these get reset for every electron
    pass_SigmaIEtaIEtaFull5x5 = false;
    pass_dEtaIn = false;
    pass_dPhiIn = false;
    pass_HoverE = false;
    pass_iso = false;
    pass_ooEmooP = false;
    pass_d0 = false;
    pass_dz = false;
    pass_missingHits = false;
    pass_convVeto = false;
    //Find EA for corrected relative iso.
    if(abs(eleSCEta->at(i)) <= 1.0)
      EA = 0.1752;
    else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
      EA = 0.1862;
    else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
      EA = 0.1411;
    else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
      EA = 0.1534;
    else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
      EA = 0.1903;
    else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
      EA = 0.2243;
    else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
      EA = 0.2687;
    EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

    if(abs(eleSCEta->at(i)) <= 1.479)
    {
      pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0101;
      pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00926;
      pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.0336;
      pass_HoverE = eleHoverE->at(i) < 0.0597;
      pass_iso = EAcorrIso < 0.0354;
      pass_ooEmooP = eleEoverPInv->at(i) < 0.012;
      pass_d0 = abs(eleD0->at(i)) < 0.0111;
      pass_dz = abs(eleDz->at(i)) < 0.0466;
      pass_missingHits = eleMissHits->at(i) <= 2;
      pass_convVeto = eleConvVeto->at(i) == 1;
    }
    else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
    {
      pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0279;
      pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00724;
      pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.0918;
      pass_HoverE = eleHoverE->at(i) < 0.0615;
      pass_iso = EAcorrIso < 0.0646;
      pass_ooEmooP = eleEoverPInv->at(i) < 0.00999;
      pass_d0 = abs(eleD0->at(i)) < 0.0351;
      pass_dz = abs(eleDz->at(i)) < 0.417;
      pass_missingHits = eleMissHits->at(i) <= 1;
      pass_convVeto = eleConvVeto->at(i) == 1;
    }
      //Electron passes Loose Electron ID cuts
    // if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
    if(eleIDbit->at(i)>>1&1==1)
    {
      //Electron passes pt cut
      if(elePt->at(i) > elePtCut)
      {
        //Electron does not overlap photon
        if(dR(eleEta->at(i),elePhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
        {
          ele_cands.push_back(i);
        }
      }
    }
  }
  return ele_cands;
}

std::vector<int> postAnalyzer::muon_veto_tightID(int pho_index, float muPtCut)
{
  // bool veto_passed = true; //pass veto if no good muon found
  std::vector<int> mu_cands;
  mu_cands.clear();

  bool pass_PFMuon = true;
  bool pass_globalMuon = true;
  // bool pass_trackerMuon = true;
  bool pass_chi2ndf = false;
  bool pass_chamberHit = false;
  bool pass_matchedStations = false;
  bool pass_dxy = false;
  bool pass_dz = false;
  bool pass_pixelHits = false;
  bool pass_trackLayers = false;
  bool pass_iso = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
  {
    // pass_globalMuon = muIsGlobalMuon->at(i);
    // pass_PFMuon = muIsPFMuon->at(i);
    // pass_trackerMuon = muIsTrackerMuon->at(i);
    pass_chi2ndf = muChi2NDF->at(i) < 10.0;
    pass_chamberHit = muMuonHits->at(i) > 0;
    pass_matchedStations = muStations->at(i) > 1;
    pass_dxy = fabs(muInnerD0->at(i)) < 0.2;
    pass_dz = fabs(muInnerDz->at(i)) < 0.5;
    pass_pixelHits = muPixelHits->at(i) > 0;
    pass_trackLayers = muTrkLayers->at(i) > 5;

    muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
    tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
    pass_iso = tightIso_combinedRelative < 0.15;
    //Muon passes Tight Muon ID
    // if(pass_iso && pass_globalMuon && pass_PFMuon && pass_chi2ndf && pass_chamberHit && pass_matchedStations && pass_dxy && pass_dz && pass_pixelHits && pass_trackLayers)
    if(pass_iso && muIDbit->at(i)>>2&1==1)
    {
      //Muon passes pt cut
      if(muPt->at(i) > muPtCut)
      {
        //Muon does not overlap photon
        if(dR(muEta->at(i),muPhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
        {
          mu_cands.push_back(i);
        }
      }
    }
  }
  return mu_cands;
}

//Returns true if veto passed
//Veto failed if an electron is found that passes Loose Electron ID and elePtcut, and does not overlap the candidate photon within dR of 0.5
//Always true if no electrons with |SC eta| < 2.5, since ID always fails for |SC eta| > 2.

std::vector<int> postAnalyzer::electron_veto_looseID(int pho_index, int mu_index, float elePtCut)
{
  std::vector<int> ele_cands;
  ele_cands.clear();

  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
  {
    //Make sure these get reset for every electron
    pass_SigmaIEtaIEtaFull5x5 = false;
    pass_dEtaIn = false;
    pass_dPhiIn = false;
    pass_HoverE = false;
    pass_iso = false;
    pass_ooEmooP = false;
    pass_d0 = false;
    pass_dz = false;
    pass_missingHits = false;
    pass_convVeto = false;
    //Find EA for corrected relative iso.
    if(abs(eleSCEta->at(i)) <= 1.0)
      EA = 0.1752;
    else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
      EA = 0.1862;
    else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
      EA = 0.1411;
    else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
      EA = 0.1534;
    else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
      EA = 0.1903;
    else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
      EA = 0.2243;
    else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
      EA = 0.2687;
    EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

    if(abs(eleSCEta->at(i)) <= 1.479)
    {
      pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
      pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
      pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
      pass_HoverE = eleHoverE->at(i) < 0.104;
      pass_iso = EAcorrIso < 0.0893;
      pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
      pass_d0 = abs(eleD0->at(i)) < 0.0261;
      pass_dz = abs(eleDz->at(i)) < 0.41;
      pass_missingHits = eleMissHits->at(i) <= 2;
      pass_convVeto = eleConvVeto->at(i) == 1;
    }
    else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
    {
      pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
      pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
      pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
      pass_HoverE = eleHoverE->at(i) < 0.0897;
      pass_iso = EAcorrIso < 0.121;
      pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
      pass_d0 = abs(eleD0->at(i)) < 0.118;
      pass_dz = abs(eleDz->at(i)) < 0.822;
      pass_missingHits = eleMissHits->at(i) <= 1;
      pass_convVeto = eleConvVeto->at(i) == 1;
    }
      //Electron passes Loose Electron ID cuts
    // if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
    if(eleIDbit->at(i)>>1&1==1)
    {
      //Electron passes pt cut
      if(elePt->at(i) > elePtCut)
      {
        //Electron does not overlap photon
        if(dR(eleEta->at(i),elePhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5 && dR(eleEta->at(i),elePhi->at(i),muEta->at(mu_index),muPhi->at(mu_index)) > 0.5)
        {
          ele_cands.push_back(i);
        }
      }
    }
  }
  return ele_cands;
}



//Veto failed if a muon is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5
std::vector<int> postAnalyzer::muon_veto_looseID(int pho_index, int ele_index, float muPtCut)
{
  std::vector<int> mu_cands;
  mu_cands.clear();

  bool veto_passed = true; //pass veto if no good muon found
  bool pass_PFMuon = true;
  bool pass_globalMuon = true;
  bool pass_trackerMuon = true;
  bool pass_iso = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
  {
    // pass_PFMuon = muIsPFMuon->at(i);
    // pass_globalMuon = muIsGlobalMuon->at(i);
    // pass_trackerMuon = muIsTrackerMuon->at(i);
    muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
    tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
    pass_iso = tightIso_combinedRelative < 0.25;
    //Muon passes Loose Muon ID and PF-based combined relative, dBeta-corrected Loose Muon Isolation cuts
  //      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
        // if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon))
    if(muIDbit->at(i)>>0&1==1)
    {
      //Muon passes pt cut
      if(muPt->at(i) > muPtCut)
      {
        //Muon does not overlap photon
        if(dR(muEta->at(i),muPhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
        {
          mu_cands.push_back(i);
        }
      }
    }
  }
  return mu_cands;
}


std::vector<int> postAnalyzer::JetVetoDecision(int pho_index, int mu_index) {

  bool jetVeto=true;
  std::vector<int> jetindex;
  float value =0.0;

  for(int i = 0; i < nJet; i++)
    {

      if(0.0 < abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.5)
        value =-0.8;
      else if(2.5 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.75)
        value =-0.95;
      else if(2.75 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <3.0)
        value =-0.97;
      else if(3.00 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <5.0)
        value =-0.99;
      else
        continue;



      //      std::cout<<"Jet size: "<<nJet<<std::endl;
      double deltar = 0.0 ;
      double deltar_mu = 0.0 ;
      //      std::cout<<"Jet no:"<<i<<"coming here pujetid: "<<pfJet_pt[i]<<std::endl;
      //if(OverlapWithMuon(jetEta->at(i),jetPhi->at(i)))     continue;
      //      std::cout<<"Jet no:"<<i<<"coming here OverlapWithMuon: "<<pfJet_pt[i]<<std::endl;
      //      if(OverlapWithElectron(jetEta->at(i),jetPhi->at(i)))   continue;
        deltar= dR(jetEta->at(i),jetPhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index));
        deltar_mu = dR(jetEta->at(i),jetPhi->at(i),muEta->at(mu_index),muPhi->at(mu_index));
      if(deltar>0.4 && deltar_mu > 0.4 && jetPt->at(i) >30.0 && jetPFLooseId->at(i)==1 && jetPUID->at(i)>value)
        {
          jetindex.push_back(i);
        }

      
    }


  //  std::cout<<"Jet size: "<< jetindex.size()<<std::endl;
  //if(jetindex.size()>1)jetVeto = false;
  return jetindex;

}



bool postAnalyzer::OverlapWithMuon(double eta, double phi){
  
  bool overlap = false;
  //  std::cout<<"No of muon:"<<Muon_n<<std::endl;
  bool pass_PFMuon = false;
  bool pass_globalMuon = false;
  bool pass_trackerMuon = false;
  bool pass_iso = false;

  Float_t zero1 = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {
      // pass_PFMuon = muIsPFMuon->at(i);
      // pass_globalMuon = muIsGlobalMuon->at(i);
      // pass_trackerMuon = muIsTrackerMuon->at(i);
      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero1,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;
      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
        {
          if(muPt->at(i) > 10.)
            {
              if(dR(muEta->at(i),muPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }

  return overlap;


}



bool postAnalyzer::OverlapWithElectron(double eta, double phi){
  bool overlap = false;

  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;

  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
    {
      pass_SigmaIEtaIEtaFull5x5 = false;
      pass_dEtaIn = false;
      pass_dPhiIn = false;
      pass_HoverE = false;
      pass_iso = false;
      pass_ooEmooP = false;
      pass_d0 = false;
      pass_dz = false;
      pass_missingHits = false;
      pass_convVeto = false;

      if(abs(eleSCEta->at(i)) <= 1.0)
        EA = 0.1752;
      else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
        EA = 0.1862;
      else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
        EA = 0.1411;
      else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
        EA = 0.1534;
      else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
        EA = 0.1903;
      else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
        EA = 0.2243;
      else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
        EA = 0.2687;
      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

      if(abs(eleSCEta->at(i)) <= 1.479)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
          pass_HoverE = eleHoverE->at(i) < 0.104;
          pass_iso = EAcorrIso < 0.0893;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
          pass_d0 = abs(eleD0->at(i)) < 0.0261;
          pass_dz = abs(eleDz->at(i)) < 0.41;
          pass_missingHits = eleMissHits->at(i) <= 2;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
          pass_HoverE = eleHoverE->at(i) < 0.0897;
          pass_iso = EAcorrIso < 0.121;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
          pass_d0 = abs(eleD0->at(i)) < 0.118;
          pass_dz = abs(eleDz->at(i)) < 0.822;
          pass_missingHits = eleMissHits->at(i) <= 1;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      //Electron passes Loose Electron ID cuts
      // if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
      if(eleIDbit->at(i)>>1&1==1)
        {
          //Electron passes pt cut
          if(elePt->at(i) > 10.)
            {
              //Electron does not overlap photon
              if(dR(eleSCEta->at(i),eleSCPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }





}


double postAnalyzer::dR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}





bool postAnalyzer::dPhiJetMET_veto(std::vector<int> jets)
{
  //pass veto if no jet is found within DeltaPhi(jet,MET) < 0.5
  bool passes = false;
  int njetsMax = jets.size();
  //Only look at first four jets
  if(njetsMax > 4)
    njetsMax = 4;
  int j = 0;
  for(; j < njetsMax; j++)
    {
      //fail veto if a jet is found with DeltaPhi(jet,MET) < 0.5
      if(DeltaPhi(jetPhi->at(jets[j]),pfMETPhi) < 0.5)
  break;
    }
  if(j==njetsMax)
    passes = true;

  return passes;
}

double postAnalyzer::FakeRatePt(Float_t phoPt, int syst_number)
{
   // bin: Photon pT bin.  Index within binX: syst_number of desired systematic shift.
  //Systematics bins:
  //    0       1      2       3       4       5       6         7          8           9            10
  // Standard  sbUP  sbDown  metUP  metDown  binUP  binDown  sieieLeft  sieieRight  templateUp  templateDown
  Float_t bin1[11] = {0.0955658,0.0946783,0.096761,0.0941346,0.09363,0.0969758,0.0952309,0.0966213,0.0945601,0.0963828,0.0947487};
  Float_t bin2[11] = {0.079111,0.0779771,0.080763,0.0779608,0.0796352,0.0796462,0.0777209,0.0805415,0.0778359,0.0799271,0.0782949};
  Float_t bin3[11] = {0.0854649,0.0838512,0.0875113,0.0838365,0.0858884,0.0861928,0.0843523,0.0870061,0.084068,0.0873089,0.0836209};
  Float_t bin4[11] = {0.097131,0.0938181,0.099658,0.0945131,0.0942899,0.0992382,0.0948404,0.0989389,0.0957066,0.100559,0.0937028};
  Float_t bin5[11] = {0.119954,0.118087,0.120222,0.115721,0.11712,0.118852,0.118129,0.122956,0.117293,0.131101,0.108808};
  
  double weight = 1.0;
  if(phoPt <= 200.0)
    weight = bin1[syst_number];
  else if(200.0 < phoPt <= 250.0)
    weight = bin2[syst_number];
  else if(250.0 < phoPt <= 300.0)
    weight = bin3[syst_number];
  else if(300.0 < phoPt <= 400.0)
    weight = bin4[syst_number];
  else if(400.0 < phoPt <= 1000.0)
    weight = bin5[syst_number];
  return weight;
}
