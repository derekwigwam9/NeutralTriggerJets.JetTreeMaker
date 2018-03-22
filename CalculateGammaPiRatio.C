// 'CalculateGammaPiRatio.C'
// Derek Anderson
// 01.23.2018
//
// This takes the output of 'StJetTreeMaker' and
// plots the ratio of the pi0-spectrum to the
// gamma-spectrum.
//
// NOTE: the 1st eTtrg "bin" is defined to be
//       the entire eTtrg range.
//
// NOTE: minimum areas for different R:
//   R = 0.3 -- 0.2   R = 0.5 -- 0.65
//   R = 0.4 -- 0.35  R = 0.7 -- 1.2


#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t NBinsEt  = 4;
static const UInt_t NBinsId  = 2;
static const UInt_t NBinsPt  = 35;
static const UInt_t NPadsPt  = 2;
static const UInt_t NRegions = 2;
static const UInt_t JetType  = 0;



void CalculateGammaPiRatio(const Bool_t isInBatchMode=false) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning gamma-pi0 ratio calculation..." << endl;

  // io parameters
  const TString sOutput("pp200r9.withFitCutPlot.r03a02rm1chrg.d28m2y2018.root");
  const TString sInput("output/CollaborationMeetingJan2018/pp200r9.withOaCones.eTtrg920.r03rm1chrg.d21m1y2018.root");
  const TString sTree("JetTree");

  // trigger parameters
  const Double_t hTrgMax(0.9);
  const Double_t pi0tsp[2] = {0., 0.08};
  const Double_t gamTsp[2] = {0.2, 0.6};

  // jet parameters
  const TString  sRes("0.3");
  const Double_t rJet(0.3);
  const Double_t aJetMin(0.2);
  const Double_t pTjetMin(0.2);
  const Double_t dFrecoil(TMath::PiOver4());
  const Double_t dFdown[2] = {TMath::PiOver4(), TMath::PiOver2()};
  const Double_t dFup[2]   = {(3. * TMath::PiOver2()), (7. * TMath::PiOver4())};


  // open files and grab tree
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fInput) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  TTree *tInput = (TTree*) fInput -> Get(sTree.Data());
  if (!tInput) {
    cerr << "PANIC: couldn't grab tree!" << endl;
    return;
  }
  cout << "    Grabbed tree." << endl;


  // event leaves
  Int_t    EventIndex;
  Int_t    NJets;
  Double_t Refmult;
  Double_t TSP;
  Double_t TrgEta;
  Double_t TrgPhi;
  Double_t TrgEt;
  Double_t Rho;
  Double_t Sigma;
  Double_t Vz;
  // jet leaves
  vector<Double_t> *JetPt;
  vector<Double_t> *JetNCons;
  vector<Double_t> *JetIndex;
  vector<Double_t> *JetPtCorr;
  vector<Double_t> *JetEta;
  vector<Double_t> *JetPhi;
  vector<Double_t> *JetE;
  vector<Double_t> *JetArea;
  vector<Double_t> *JetPtOffAxisUp;
  vector<Double_t> *JetPtOffAxisDown;
  // constituent leaves
  vector<vector<Double_t>> *JetConsPt;
  vector<vector<Double_t>> *JetConsEta;
  vector<vector<Double_t>> *JetConsPhi;
  vector<vector<Double_t>> *JetConsE;

  // declare branches
  TBranch *bEventIndex;
  TBranch *bNJets;
  TBranch *bRefmult;
  TBranch *bTsp;
  TBranch *bTrgEta;
  TBranch *bTrgPhi;
  TBranch *bTrgEt;
  TBranch *bRho;
  TBranch *bSigma;
  TBranch *bVz;
  TBranch *bJetNCons;
  TBranch *bJetIndex;
  TBranch *bJetPt;
  TBranch *bJetPtCorr;
  TBranch *bJetEta;
  TBranch *bJetPhi;
  TBranch *bJetE;
  TBranch *bJetArea;
  TBranch *bJetPtOffAxisUp;
  TBranch *bJetPtOffAxisDown;
  TBranch *bJetConsPt;
  TBranch *bJetConsEta;
  TBranch *bJetConsPhi;
  TBranch *bJetConsE;

  // set branches
  tInput -> SetBranchAddress("eventIndex", &EventIndex, &bEventIndex);
  tInput -> SetBranchAddress("Refmult", &Refmult, &bRefmult);
  tInput -> SetBranchAddress("NJets", &NJets, &bNJets);
  tInput -> SetBranchAddress("TSP", &TSP, &bTsp);
  tInput -> SetBranchAddress("TrgEta", &TrgEta, &bTrgEta);
  tInput -> SetBranchAddress("TrgPhi", &TrgPhi, &bTrgPhi);
  tInput -> SetBranchAddress("TrgEt", &TrgEt, &bTrgEt);
  tInput -> SetBranchAddress("Rho", &Rho, &bRho);
  tInput -> SetBranchAddress("Sigma", &Sigma, &bSigma);
  tInput -> SetBranchAddress("Vz", &Vz, &bVz);
  tInput -> SetBranchAddress("JetIndex", &JetIndex, &bJetIndex);
  tInput -> SetBranchAddress("JetNCons", &JetNCons, &bJetNCons);
  tInput -> SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tInput -> SetBranchAddress("JetPtCorr", &JetPtCorr, &bJetPtCorr);
  tInput -> SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tInput -> SetBranchAddress("JetPhi",&JetPhi, &bJetPhi); 
  tInput -> SetBranchAddress("JetE", &JetE, &bJetE); 
  tInput -> SetBranchAddress("JetArea", &JetArea, &bJetArea);
  tInput -> SetBranchAddress("JetPtOffAxisUp", &JetPtOffAxisUp, &bJetPtOffAxisUp); 
  tInput -> SetBranchAddress("JetPtOffAxisDown", &JetPtOffAxisDown, &bJetPtOffAxisDown); 
  tInput -> SetBranchAddress("JetConsPt", &JetConsPt, &bJetConsPt);
  tInput -> SetBranchAddress("JetConsEta", &JetConsEta, &bJetConsEta);
  tInput -> SetBranchAddress("JetConsPhi", &JetConsPhi, &bJetConsPhi);
  tInput -> SetBranchAddress("JetConsE", &JetConsE, &bJetConsE);
  cout << "    Set branches." << endl;


  // define histograms
  TH1D *hPtRE[NBinsEt][NBinsId];
  TH1D *hPtUE[NBinsEt][NBinsId];

  const UInt_t   nPt(NBinsPt - 1);
  const Double_t pTbin[NBinsPt] = {-5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 30., 35., 40., 45., 50.};
  // recoil distributions
  hPtRE[0][0]  = new TH1D("hPtRecoil_pi920", "", nPt, pTbin);
  hPtRE[0][1]  = new TH1D("hPtRecoil_ga920", "", nPt, pTbin);
  hPtRE[1][0]  = new TH1D("hPtRecoil_pi911", "", nPt, pTbin);
  hPtRE[1][1]  = new TH1D("hPtRecoil_ga911", "", nPt, pTbin);
  hPtRE[2][0]  = new TH1D("hPtRecoil_pi1115", "", nPt, pTbin);
  hPtRE[2][1]  = new TH1D("hPtRecoil_ga1115", "", nPt, pTbin);
  hPtRE[3][0]  = new TH1D("hPtRecoil_pi1520", "", nPt, pTbin);
  hPtRE[3][1]  = new TH1D("hPtRecoil_ga1520", "", nPt, pTbin);
  // off-axis distributions
  hPtUE[0][0]  = new TH1D("hPtUncorrelated_pi920", "", nPt, pTbin);
  hPtUE[0][1]  = new TH1D("hPtUncorrelated_ga920", "", nPt, pTbin);
  hPtUE[1][0]  = new TH1D("hPtUncorrelated_pi911", "", nPt, pTbin);
  hPtUE[1][1]  = new TH1D("hPtUncorrelated_ga911", "", nPt, pTbin);
  hPtUE[2][0]  = new TH1D("hPtUncorrelated_pi1115", "", nPt, pTbin);
  hPtUE[2][1]  = new TH1D("hPtUncorrelated_ga1115", "", nPt, pTbin);
  hPtUE[3][0]  = new TH1D("hPtUncorrelated_pi1520", "", nPt, pTbin);
  hPtUE[3][1]  = new TH1D("hPtUncorrelated_ga1520", "", nPt, pTbin);
  // errors
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      hPtRE[iBinEt][iBinId] -> Sumw2();
      hPtUE[iBinEt][iBinId] -> Sumw2();
    }
  }
  cout << "    Defined histograms." << endl;


  // define eTtrg bins
  const Double_t eTtrgMin[NBinsEt]  = {9., 9., 11., 15.};
  const Double_t eTtrgMax[NBinsEt]  = {20., 11., 15., 20.};
  cout << "    Defined various bins." << endl;


  // no. of triggers
  UInt_t nTrgBin[NBinsEt][NBinsId];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      nTrgBin[iBinEt][iBinId] = 0;
    }
  }

  const UInt_t nEvts = tInput -> GetEntries();
  cout << "    Beginning event loop: " << nEvts << " events to process..." << endl;


  // event loop
  UInt_t bytes(0);
  UInt_t nBytes(0);
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    // load entry
    bytes = tInput -> GetEntry(iEvt);
    nBytes += bytes;
    if (bytes < 0) {
      cerr << "WARNING: issue with entry " << iEvt << "!" << endl;
      break;
    }
    else {
      if (isInBatchMode) {
        cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
      }
      else {
        cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
        if ((iEvt + 1) == nEvts) cout << endl;
      }
    }


    // trigger info
    const Double_t nJets = NJets;
    const Double_t hTrg  = TrgEta;
    const Double_t fTrg  = TrgPhi;
    const Double_t eTtrg = TrgEt;

    // trigger cuts
    const Bool_t isInTrgEtaCut = (TMath::Abs(hTrg) < hTrgMax);
    const Bool_t isInTrgEtCut  = ((eTtrg > eTtrgMin[0]) && (eTtrg < eTtrgMax[0]));
    const Bool_t isInPi0tspCut = ((TSP > pi0tsp[0]) && (TSP < pi0tsp[1]));
    const Bool_t isInGamTspCut = ((TSP > gamTsp[0]) && (TSP < gamTsp[1]));
    const Bool_t isInTspCut    = (isInPi0tspCut || isInGamTspCut);
    if (!isInTrgEtaCut || !isInTrgEtCut || !isInTspCut) continue;


    // determine trigger bin
    UInt_t eTbin(0);
    for (UInt_t iBinEt = 1; iBinEt < NBinsEt; iBinEt++) {
      const Bool_t isInEtBin = ((eTtrg >= eTtrgMin[iBinEt]) && (eTtrg < eTtrgMax[iBinEt]));
      if (isInEtBin) {
        eTbin = iBinEt;
        break;
      }
    }

    UInt_t idBin(0);
    if (isInPi0tspCut) {
      idBin = 0;
      nTrgBin[0][0]++;
      nTrgBin[eTbin][0]++;
    }
    if (isInGamTspCut) {
      idBin = 1;
      nTrgBin[0][1]++;
      nTrgBin[eTbin][1]++;
    }


    // jet loop 1
    Bool_t recoilJetIsPresent(false);
    for (UInt_t iJet = 0; iJet < nJets; iJet++) {

      // jet info
      const Double_t fJet  = JetPhi  -> at(iJet);
      const Double_t aJet  = JetArea -> at(iJet);
      const Double_t pTjet = JetPt   -> at(iJet);

      Double_t dFjet = fJet - fTrg;
      if (dFjet < 0.)
        dFjet += TMath::TwoPi();
      if (dFjet > TMath::TwoPi())
        dFjet -= TMath::TwoPi();

      // jet cuts
      const Bool_t isInAreaCut = (aJet > aJetMin);
      const Bool_t isInPtCut   = (pTjet > pTjetMin);
      const Bool_t isRecoil    = (TMath::Abs(dFjet - TMath::Pi()) < dFrecoil);
      if (!isInAreaCut || !isInPtCut || !isRecoil) {
        continue;
      }
      else {
        recoilJetIsPresent = true;
        break;
      }

    }  // end jet loop 1


    // jet loop 2
    for (UInt_t iJet = 0; iJet < nJets; iJet++) {

      // jet info
      const Double_t hJet    = JetEta           -> at(iJet);
      const Double_t fJet    = JetPhi           -> at(iJet);
      const Double_t aJet    = JetArea          -> at(iJet);
      const Double_t pTjet   = JetPt            -> at(iJet);
      const Double_t pTreco  = JetPtCorr        -> at(iJet);
      const Double_t rhoUp   = JetPtOffAxisUp   -> at(iJet);
      const Double_t rhoDown = JetPtOffAxisDown -> at(iJet);
      const Double_t rhoAvg  = (rhoUp + rhoDown) / 2.;
      const Double_t pTup    = rhoUp * aJet;
      const Double_t pTdown  = rhoDown * aJet;
      const Double_t pTavg   = rhoAvg * aJet;

      Double_t dFjet = fJet - fTrg;
      if (dFjet < 0.)
        dFjet += TMath::TwoPi();
      if (dFjet > TMath::TwoPi())
        dFjet -= TMath::TwoPi();

      // jet cuts
      const Bool_t isInAreaCut = (aJet > aJetMin);
      const Bool_t isInPtCut   = (pTjet > pTjetMin);
      if (!isInAreaCut || !isInPtCut) continue;


      // recoil jets
      const Bool_t isRecoil = (TMath::Abs(dFjet - TMath::Pi()) < dFrecoil);
      if (isRecoil) {
        hPtRE[0][idBin]     -> Fill(pTreco);
        hPtRE[eTbin][idBin] -> Fill(pTreco);
      }

      // off-axis jets
      const Bool_t isInLowerCone = ((dFjet > dFdown[0]) && (dFjet < dFdown[1]));
      const Bool_t isInUpperCone = ((dFjet > dFup[0]) && (dFjet < dFup[1]));
      const Bool_t isInUeRegion  = (isInLowerCone || isInUpperCone);
      if (isInUeRegion && recoilJetIsPresent) {
        hPtUE[0][idBin]     -> Fill(pTreco);
        hPtUE[eTbin][idBin] -> Fill(pTreco);
      }

    }  // end jet loop 2

  }  // end event loop

  cout << "    Finished event loop:\n"
       << "      " << nTrgBin[0][0] << " pi0 triggers,\n"
       << "      " << nTrgBin[0][1] << " gamma trigger."
       << endl;


  // normalize histograms
  const Double_t hBin = 2. * (1 - rJet);
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      const Double_t pTnorm = hBin * nTrgBin[iBinEt][iBinId];
      hPtRE[iBinEt][iBinId] -> Scale(1. / pTnorm);
      hPtUE[iBinEt][iBinId] -> Scale(1. / pTnorm);
      for (UInt_t iBinPt = 1; iBinPt < NBinsPt + 1; iBinPt++) {
        const Double_t pTwidth   = hPtRE[iBinEt][iBinId] -> GetBinWidth(iBinPt);
        const Double_t pTvalueRE = hPtRE[iBinEt][iBinId] -> GetBinContent(iBinPt);
        const Double_t pTvalueUE = hPtUE[iBinEt][iBinId] -> GetBinContent(iBinPt);
        const Double_t pTscaleRE = pTvalueRE / pTwidth;
        const Double_t pTscaleUE = pTvalueUE / pTwidth;
        hPtRE[iBinEt][iBinId] -> SetBinContent(iBinPt, pTscaleRE);
        hPtUE[iBinEt][iBinId] -> SetBinContent(iBinPt, pTscaleUE);
      }
    }
  }
  cout << "    Normalized histograms." << endl;


  // calculate ratios
  const Double_t weight(1.);
  const TString  baseRE("hRatioRE_et");
  const TString  baseUE("hRatioUE_et");
  const TString  eTname[NBinsEt] = {"920", "911", "1115", "1520"};

  TH1D *hRatioRE[NBinsEt];
  TH1D *hRatioUE[NBinsEt];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    TString sRatioRE(baseRE.Data());
    TString sRatioUE(baseUE.Data());
    sRatioRE += eTname[iBinEt].Data();
    sRatioUE += eTname[iBinEt].Data();

    hRatioRE[iBinEt] = new TH1D(sRatioRE.Data(), "", nPt, pTbin);
    hRatioUE[iBinEt] = new TH1D(sRatioUE.Data(), "", nPt, pTbin);
    hRatioRE[iBinEt] -> Sumw2();
    hRatioUE[iBinEt] -> Sumw2();
    hRatioRE[iBinEt] -> Divide(hPtRE[iBinEt][0], hPtRE[iBinEt][1], weight, weight);
    hRatioUE[iBinEt] -> Divide(hPtUE[iBinEt][0], hPtUE[iBinEt][1], weight, weight);
  }
  cout << "    Calculated ratios." << endl;


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fSiz(1.3);
  const UInt_t  fColRR(819);
  const UInt_t  fColRU(889);
  const UInt_t  fMarRR(20);
  const UInt_t  fMarRU(23);
  const UInt_t  fColR[NBinsId] = {859, 899};
  const UInt_t  fColU[NBinsId] = {809, 849};
  const UInt_t  fMarR[NBinsId] = {29, 20};
  const UInt_t  fMarU[NBinsId] = {22, 23};
  const Float_t fLab(0.03);
  const Float_t fLabR(0.055);
  const Float_t fSizR(0.074);
  const Float_t fOffsetX(1.);
  const Float_t fOffsetXR(0.77);
  const Float_t fOffsetY(1.07);
  const Float_t fOffsetYR(0.57);
  const Float_t fRangeY[2]  = {0.00003, 1.3};
  const Float_t fRangeYR[2] = {0.03, 13.};
  const TString sTitle("Recoil jet p_{T}^{reco}, no N_{fit} cut");
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - #rho #upoint A_{jet}");
  const TString sTitleY("(1/N^{trg}) dN^{jet}/(dp_{T}^{jet} d#eta)");
  const TString sTitleR("#pi^{0} / #gamma^{rich}");
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      // recoil jets
      hPtRE[iBinEt][iBinId] -> SetLineColor(fColR[iBinId]);
      hPtRE[iBinEt][iBinId] -> SetMarkerStyle(fMarR[iBinId]);
      hPtRE[iBinEt][iBinId] -> SetMarkerColor(fColR[iBinId]);
      hPtRE[iBinEt][iBinId] -> SetMarkerSize(fSiz);
      hPtRE[iBinEt][iBinId] -> SetTitle(sTitle.Data());
      hPtRE[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetRangeUser(fRangeY[0], fRangeY[1]);

      // off-axis jets
      hPtUE[iBinEt][iBinId] -> SetLineColor(fColU[iBinId]);
      hPtUE[iBinEt][iBinId] -> SetMarkerStyle(fMarU[iBinId]);
      hPtUE[iBinEt][iBinId] -> SetMarkerColor(fColU[iBinId]);
      hPtUE[iBinEt][iBinId] -> SetMarkerSize(fSiz);
      hPtUE[iBinEt][iBinId] -> SetTitle(sTitle.Data());
      hPtUE[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hPtUE[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hPtUE[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hPtUE[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPtUE[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hPtUE[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hPtUE[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hPtUE[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hPtUE[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hPtUE[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hPtUE[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      hPtUE[iBinEt][iBinId] -> GetYaxis() -> SetRangeUser(fRangeY[0], fRangeY[1]);

    }  // end id loop

    // recoil ratios
    hRatioRE[iBinEt] -> SetLineColor(fColRR);
    hRatioRE[iBinEt] -> SetMarkerStyle(fMarRR);
    hRatioRE[iBinEt] -> SetMarkerColor(fColRR);
    hRatioRE[iBinEt] -> SetMarkerSize(fSiz);
    hRatioRE[iBinEt] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatioRE[iBinEt] -> GetXaxis() -> SetTitleFont(fTxt);
    hRatioRE[iBinEt] -> GetXaxis() -> SetTitleSize(fSizR);
    hRatioRE[iBinEt] -> GetXaxis() -> SetTitleOffset(fOffsetXR);
    hRatioRE[iBinEt] -> GetXaxis() -> CenterTitle(fCnt);
    hRatioRE[iBinEt] -> GetXaxis() -> SetLabelSize(fLabR);
    hRatioRE[iBinEt] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatioRE[iBinEt] -> GetYaxis() -> SetTitleFont(fTxt);
    hRatioRE[iBinEt] -> GetYaxis() -> SetTitleSize(fSizR);
    hRatioRE[iBinEt] -> GetYaxis() -> SetTitleOffset(fOffsetYR);
    hRatioRE[iBinEt] -> GetYaxis() -> CenterTitle(fCnt);
    hRatioRE[iBinEt] -> GetYaxis() -> SetLabelSize(fLabR);
    hRatioRE[iBinEt] -> GetYaxis() -> SetRangeUser(fRangeYR[0], fRangeYR[1]);

    // off-axis ratios
    hRatioUE[iBinEt] -> SetLineColor(fColRU);
    hRatioUE[iBinEt] -> SetMarkerStyle(fMarRU);
    hRatioUE[iBinEt] -> SetMarkerColor(fColRU);
    hRatioUE[iBinEt] -> SetMarkerSize(fSiz);
    hRatioUE[iBinEt] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatioUE[iBinEt] -> GetXaxis() -> SetTitleFont(fTxt);
    hRatioUE[iBinEt] -> GetXaxis() -> SetTitleSize(fSizR);
    hRatioUE[iBinEt] -> GetXaxis() -> SetTitleOffset(fOffsetXR);
    hRatioUE[iBinEt] -> GetXaxis() -> CenterTitle(fCnt);
    hRatioUE[iBinEt] -> GetXaxis() -> SetLabelSize(fLabR);
    hRatioUE[iBinEt] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatioUE[iBinEt] -> GetYaxis() -> SetTitleFont(fTxt);
    hRatioUE[iBinEt] -> GetYaxis() -> SetTitleSize(fSizR);
    hRatioUE[iBinEt] -> GetYaxis() -> SetTitleOffset(fOffsetYR);
    hRatioUE[iBinEt] -> GetYaxis() -> CenterTitle(fCnt);
    hRatioUE[iBinEt] -> GetYaxis() -> SetLabelSize(fLabR);
    hRatioUE[iBinEt] -> GetYaxis() -> SetRangeUser(fRangeYR[0], fRangeYR[1]);

  }  // end eTtrg loop
  cout << "    Set styles." << endl;


  // make labels
  const UInt_t  fColT(0);
  const UInt_t  fFilT(0);
  const UInt_t  fLinT(0);
  const UInt_t  fAlign(12);
  const Float_t xyLeg1(0.1);
  const Float_t xyLeg2(0.3);
  const Float_t xyPav1(0.5);
  const Float_t xyPav2(0.7);
  const TString sSystem("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrigger("E_{T}^{trg} #in (");
  const TString sJet("anti-k_{T}, R = ");
  const TString sChrg("#bf{charged jets}");
  const TString sFull("#bf{full jets}");
  const TString sIdNameR[NBinsId] = {"#pi^{0}: recoil jets", "#gamma^{rich}: recoil jets"};
  const TString sIdNameU[NBinsId] = {"#pi^{0}: off-axis jets", "#gamma^{rich}: off-axis jets"};
  const TString sRatio[NRegions]  = {"recoil jets", "off-axis jets"};

  TLegend   *lPtJet[NBinsEt];
  TLegend   *lRatio[NBinsEt];
  TPaveText *pPtTxt[NBinsEt];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    TString sEtTrg(sTrigger.Data());
    TString sJetStuff(sJet.Data());
    sEtTrg    += eTtrgMin[iBinEt];
    sEtTrg    += ", ";
    sEtTrg    += eTtrgMax[iBinEt];
    sEtTrg    += ") GeV";
    sJetStuff += sRes.Data();

    pPtTxt[iBinEt] = new TPaveText(xyPav1, xyPav1, xyPav2, xyPav2, "NDC NB");
    pPtTxt[iBinEt] -> SetFillColor(fColT);
    pPtTxt[iBinEt] -> SetFillStyle(fFilT);
    pPtTxt[iBinEt] -> SetLineColor(fColT);
    pPtTxt[iBinEt] -> SetLineStyle(fLinT);
    pPtTxt[iBinEt] -> SetTextFont(fTxt);
    pPtTxt[iBinEt] -> SetTextAlign(fAlign);
    pPtTxt[iBinEt] -> AddText(sSystem.Data());
    pPtTxt[iBinEt] -> AddText(sEtTrg.Data());
    pPtTxt[iBinEt] -> AddText(sJetStuff.Data());
    if (JetType == 0)
      pPtTxt[iBinEt] -> AddText(sChrg.Data());
    else
      pPtTxt[iBinEt] -> AddText(sFull.Data());

    lPtJet[iBinEt] = new TLegend(xyLeg1, xyLeg1, xyLeg2, xyLeg2);
    lPtJet[iBinEt] -> SetFillColor(fColT);
    lPtJet[iBinEt] -> SetFillStyle(fFilT);
    lPtJet[iBinEt] -> SetLineColor(fColT);
    lPtJet[iBinEt] -> SetLineStyle(fLinT);
    lPtJet[iBinEt] -> SetTextFont(fTxt);
    lPtJet[iBinEt] -> AddEntry(hPtRE[0][0], sIdNameR[0].Data());
    lPtJet[iBinEt] -> AddEntry(hPtUE[0][0], sIdNameU[0].Data());
    lPtJet[iBinEt] -> AddEntry(hPtRE[0][1], sIdNameR[1].Data());
    lPtJet[iBinEt] -> AddEntry(hPtUE[0][1], sIdNameU[1].Data());

    lRatio[iBinEt] = new TLegend(xyLeg1, xyLeg1, xyLeg2, xyLeg2);
    lRatio[iBinEt] -> SetFillColor(fColT);
    lRatio[iBinEt] -> SetFillStyle(fFilT);
    lRatio[iBinEt] -> SetLineColor(fColT);
    lRatio[iBinEt] -> SetLineStyle(fLinT);
    lRatio[iBinEt] -> SetTextFont(fTxt);
    lRatio[iBinEt] -> AddEntry(hRatioRE[0], sRatio[0].Data());
    lRatio[iBinEt] -> AddEntry(hRatioUE[0], sRatio[1].Data());
  }  // end eTtrg loop
  cout << "    Made labels." << endl;


  // make line
  const UInt_t fColL(1);
  const UInt_t fLinL(2);

  TLine *lOne[NBinsEt];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    lOne[iBinEt] = new TLine(pTbin[0], 1., pTbin[NBinsPt - 1], 1.);
    lOne[iBinEt] -> SetLineColor(fColL);
    lOne[iBinEt] -> SetLineStyle(fLinL);
  }
  cout << "    Made line." << endl;


  // create directories and save histograms
  const TString sBinsEt[NBinsEt] = {"eTtrg920", "eTtrg911", "eTtrg1115", "eTtrg1520"};

  TDirectory *dBinsEt[NBinsEt];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    dBinsEt[iBinEt] = (TDirectory*) fOutput -> mkdir(sBinsEt[iBinEt].Data());
    dBinsEt[iBinEt]  -> cd();
    hRatioRE[iBinEt] -> Write();
    hRatioUE[iBinEt] -> Write();
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      hPtRE[iBinEt][iBinId] -> Write();
      hPtUE[iBinEt][iBinId] -> Write();
    }  // end id loop
  }  // end eTtrg loop
  cout << "    Made directories." << endl;


  // draw plots
  const UInt_t  widthPt(750);
  const UInt_t  heightPt(900);
  const UInt_t  widthRt(1500);
  const UInt_t  heightRt(900);
  const UInt_t  frame(0);
  const UInt_t  grid(0);
  const UInt_t  tick(1);
  const UInt_t  log(1);
  const Float_t noMargin(0.);
  const Float_t marginR(0.07);
  const Float_t marginL(0.13);
  const Float_t marginT(0.07);
  const Float_t marginB(0.13);
  const Float_t marginBR(0.17);
  const Float_t xPadPt[NPadsPt + 2]  = {0., 1., 0., 1.};
  const Float_t yPadPt[NPadsPt + 2]  = {0., 0.35, 0.35, 1.};
  const Float_t xPadRt[NBinsEt + 4]  = {0., 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.};
  const Float_t yPadRtR[NBinsEt + 4] = {0., 0.35, 0., 0.35, 0., 0.35, 0., 0.35};
  const Float_t yPadRtJ[NBinsEt + 4] = {0.35, 1., 0.35, 1., 0.35, 1., 0.35, 1.};
  const TString sPtCanvasR[NBinsEt]  = {"cPtRatio_et920", "cPtRatio_et911", "cPtRatio_et1115", "cPtRatio_et1520"};
  const TString sPtCanvasRT          = "cPtRatio_all";                                                             
  const TString sPtPads[NPadsPt]     = {"pRatio", "pJets"};                                                        
  const TString sRtPadsR[NBinsEt]    = {"pRatio920", "pRatio911", "pRatio1115", "pRatio1520"};
  const TString sRtPadsJ[NBinsEt]    = {"pJets920", "pJets911", "pJets1115", "pJets1520"};                                 

  TPad    *pPtJet[NBinsEt][NPadsPt];
  TPad    *pRtJet[NBinsEt][NPadsPt];
  TCanvas *cPtJet[NBinsEt];
  TCanvas *cRtJet;
  fOutput -> cd();
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    dBinsEt[iBinEt] -> cd();
    cPtJet[iBinEt]    = new TCanvas(sPtCanvasR[iBinEt].Data(), "", widthPt, heightPt);
    pPtJet[iBinEt][0] = new TPad(sPtPads[0].Data(), "", xPadPt[0], yPadPt[0], xPadPt[1], yPadPt[1]);
    pPtJet[iBinEt][1] = new TPad(sPtPads[1].Data(), "", xPadPt[2], yPadPt[2], xPadPt[3], yPadPt[3]);
    pPtJet[iBinEt][0] -> SetGrid(grid, grid);
    pPtJet[iBinEt][0] -> SetTicks(tick, tick);
    pPtJet[iBinEt][0] -> SetLogy(log);
    pPtJet[iBinEt][0] -> SetFrameBorderMode(frame);
    pPtJet[iBinEt][0] -> SetRightMargin(marginR);
    pPtJet[iBinEt][0] -> SetLeftMargin(marginL);
    pPtJet[iBinEt][0] -> SetTopMargin(noMargin);
    pPtJet[iBinEt][0] -> SetBottomMargin(marginBR);
    pPtJet[iBinEt][1] -> SetGrid(grid, grid);
    pPtJet[iBinEt][1] -> SetTicks(tick, tick);
    pPtJet[iBinEt][1] -> SetLogy(log);
    pPtJet[iBinEt][1] -> SetFrameBorderMode(frame);
    pPtJet[iBinEt][1] -> SetRightMargin(marginR);
    pPtJet[iBinEt][1] -> SetLeftMargin(marginL);
    pPtJet[iBinEt][1] -> SetTopMargin(marginT);
    pPtJet[iBinEt][1] -> SetBottomMargin(noMargin);
    cPtJet[iBinEt]    -> cd();
    pPtJet[iBinEt][0] -> Draw();
    pPtJet[iBinEt][1] -> Draw();
    pPtJet[iBinEt][0] -> cd();
    hRatioRE[iBinEt]  -> Draw();
    hRatioUE[iBinEt]  -> Draw("same");
    lOne[iBinEt]      -> Draw();
    lRatio[iBinEt]    -> Draw();
    pPtJet[iBinEt][1] -> cd();
    hPtRE[iBinEt][1]  -> Draw();
    hPtRE[iBinEt][0]  -> Draw("same");
    hPtUE[iBinEt][1]  -> Draw("same");
    hPtUE[iBinEt][0]  -> Draw("same");
    pPtTxt[iBinEt]    -> Draw();
    lPtJet[iBinEt]    -> Draw();
    cPtJet[iBinEt]    -> Write();
    cPtJet[iBinEt]    -> Close();
  }  // end eTtrg loop
  fOutput -> cd();
  cRtJet       = new TCanvas(sPtCanvasRT.Data(), "", widthRt, heightRt);
  pRtJet[0][0] = new TPad(sRtPadsR[0].Data(), "", xPadRt[0], yPadRtR[0], xPadRt[1], yPadRtR[1]);
  pRtJet[0][1] = new TPad(sRtPadsJ[0].Data(), "", xPadRt[0], yPadRtJ[0], xPadRt[1], yPadRtJ[1]);
  pRtJet[1][0] = new TPad(sRtPadsR[1].Data(), "", xPadRt[2], yPadRtR[2], xPadRt[3], yPadRtR[3]);
  pRtJet[1][1] = new TPad(sRtPadsJ[1].Data(), "", xPadRt[2], yPadRtJ[2], xPadRt[3], yPadRtJ[3]);
  pRtJet[2][0] = new TPad(sRtPadsR[2].Data(), "", xPadRt[4], yPadRtR[4], xPadRt[5], yPadRtR[5]);
  pRtJet[2][1] = new TPad(sRtPadsJ[2].Data(), "", xPadRt[4], yPadRtJ[4], xPadRt[5], yPadRtJ[5]);
  pRtJet[3][0] = new TPad(sRtPadsR[3].Data(), "", xPadRt[6], yPadRtR[6], xPadRt[7], yPadRtR[7]);
  pRtJet[3][1] = new TPad(sRtPadsJ[3].Data(), "", xPadRt[6], yPadRtJ[6], xPadRt[7], yPadRtJ[7]);
  pRtJet[0][0] -> SetGrid(grid, grid);
  pRtJet[0][0] -> SetTicks(tick, tick);
  pRtJet[0][0] -> SetLogy(log);
  pRtJet[0][0] -> SetFrameBorderMode(frame);
  pRtJet[0][0] -> SetRightMargin(noMargin);
  pRtJet[0][0] -> SetLeftMargin(marginL);
  pRtJet[0][0] -> SetTopMargin(noMargin);
  pRtJet[0][0] -> SetBottomMargin(marginBR);
  pRtJet[0][1] -> SetGrid(grid, grid);
  pRtJet[0][1] -> SetTicks(tick, tick);
  pRtJet[0][1] -> SetLogy(log);
  pRtJet[0][1] -> SetFrameBorderMode(frame);
  pRtJet[0][1] -> SetRightMargin(noMargin);
  pRtJet[0][1] -> SetLeftMargin(marginL);
  pRtJet[0][1] -> SetTopMargin(marginT);
  pRtJet[0][1] -> SetBottomMargin(noMargin);
  pRtJet[1][0] -> SetGrid(grid, grid);
  pRtJet[1][0] -> SetTicks(tick, tick);
  pRtJet[1][0] -> SetLogy(log);
  pRtJet[1][0] -> SetFrameBorderMode(frame);
  pRtJet[1][0] -> SetRightMargin(noMargin);
  pRtJet[1][0] -> SetLeftMargin(noMargin);
  pRtJet[1][0] -> SetTopMargin(noMargin);
  pRtJet[1][0] -> SetBottomMargin(marginBR);
  pRtJet[1][1] -> SetGrid(grid, grid);
  pRtJet[1][1] -> SetTicks(tick, tick);
  pRtJet[1][1] -> SetLogy(log);
  pRtJet[1][1] -> SetFrameBorderMode(frame);
  pRtJet[1][1] -> SetRightMargin(noMargin);
  pRtJet[1][1] -> SetLeftMargin(noMargin);
  pRtJet[1][1] -> SetTopMargin(marginT);
  pRtJet[1][1] -> SetBottomMargin(noMargin);
  pRtJet[2][0] -> SetGrid(grid, grid);
  pRtJet[2][0] -> SetTicks(tick, tick);
  pRtJet[2][0] -> SetLogy(log);
  pRtJet[2][0] -> SetFrameBorderMode(frame);
  pRtJet[2][0] -> SetRightMargin(noMargin);
  pRtJet[2][0] -> SetLeftMargin(noMargin);
  pRtJet[2][0] -> SetTopMargin(noMargin);
  pRtJet[2][0] -> SetBottomMargin(marginBR);
  pRtJet[2][1] -> SetGrid(grid, grid);
  pRtJet[2][1] -> SetTicks(tick, tick);
  pRtJet[2][1] -> SetLogy(log);
  pRtJet[2][1] -> SetFrameBorderMode(frame);
  pRtJet[2][1] -> SetRightMargin(noMargin);
  pRtJet[2][1] -> SetLeftMargin(noMargin);
  pRtJet[2][1] -> SetTopMargin(marginT);
  pRtJet[2][1] -> SetBottomMargin(noMargin);
  pRtJet[3][0] -> SetGrid(grid, grid);
  pRtJet[3][0] -> SetTicks(tick, tick);
  pRtJet[3][0] -> SetLogy(log);
  pRtJet[3][0] -> SetFrameBorderMode(frame);
  pRtJet[3][0] -> SetRightMargin(marginR);
  pRtJet[3][0] -> SetLeftMargin(noMargin);
  pRtJet[3][0] -> SetTopMargin(noMargin);
  pRtJet[3][0] -> SetBottomMargin(marginBR);
  pRtJet[3][1] -> SetGrid(grid, grid);
  pRtJet[3][1] -> SetTicks(tick, tick);
  pRtJet[3][1] -> SetLogy(log);
  pRtJet[3][1] -> SetFrameBorderMode(frame);
  pRtJet[3][1] -> SetRightMargin(marginR);
  pRtJet[3][1] -> SetLeftMargin(noMargin);
  pRtJet[3][1] -> SetTopMargin(marginT);
  pRtJet[3][1] -> SetBottomMargin(noMargin);
  cRtJet       -> cd();
  pRtJet[0][0] -> Draw();
  pRtJet[0][1] -> Draw();
  pRtJet[1][0] -> Draw();
  pRtJet[1][1] -> Draw();
  pRtJet[2][0] -> Draw();
  pRtJet[2][1] -> Draw();
  pRtJet[3][0] -> Draw();
  pRtJet[3][1] -> Draw();
  pRtJet[0][0] -> cd();
  hRatioRE[0]  -> Draw();
  hRatioUE[0]  -> Draw("same");
  lOne[0]      -> Draw();
  lRatio[1]    -> Draw();
  pRtJet[0][1] -> cd();
  hPtRE[0][1]  -> Draw();
  hPtRE[0][0]  -> Draw("same");
  hPtUE[0][1]  -> Draw("same");
  hPtUE[0][0]  -> Draw("same");
  pPtTxt[0]    -> Draw();
  lPtJet[0]    -> Draw();
  pRtJet[1][0] -> cd();
  hRatioRE[1]  -> Draw();
  hRatioUE[1]  -> Draw("same");
  lOne[1]      -> Draw();
  lRatio[1]    -> Draw();
  pRtJet[1][1] -> cd();
  hPtRE[1][1]  -> Draw();
  hPtRE[1][0]  -> Draw("same");
  hPtUE[1][1]  -> Draw("same");
  hPtUE[1][0]  -> Draw("same");
  pPtTxt[1]    -> Draw();
  lPtJet[1]    -> Draw();
  pRtJet[2][0] -> cd();
  hRatioRE[2]  -> Draw();
  hRatioUE[2]  -> Draw("same");
  lOne[2]      -> Draw();
  lRatio[2]    -> Draw();
  pRtJet[2][1] -> cd();
  hPtRE[2][1]  -> Draw();
  hPtRE[2][0]  -> Draw("same");
  hPtUE[2][1]  -> Draw("same");
  hPtUE[2][0]  -> Draw("same");
  pPtTxt[2]    -> Draw();
  lPtJet[2]    -> Draw();
  pRtJet[3][0] -> cd();
  hRatioRE[3]  -> Draw();
  hRatioUE[3]  -> Draw("same");
  lOne[3]      -> Draw();
  lRatio[3]    -> Draw();
  pRtJet[3][1] -> cd();
  hPtRE[3][1]  -> Draw();
  hPtRE[3][0]  -> Draw("same");
  hPtUE[3][1]  -> Draw("same");
  hPtUE[3][0]  -> Draw("same");
  pPtTxt[3]    -> Draw();
  lPtJet[3]    -> Draw();
  cRtJet       -> Write();
  cRtJet       -> Close();
  cout << "    Drew plots." << endl;


  // save and close
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Calculation finished!\n" << endl;

}

// End ------------------------------------------------------------------------
