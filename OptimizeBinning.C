// 'OptimizeBinning.C'
// Derek Anderson
// 03.21.2018
//
// This takes the output of 'StJetTreeMaker'
// and optimizes the binning according to
// pTreco for a given range of trigger eT.
// Needs to be compiled before running, e.g.
//   ~> .x OptimizeBinning.C++
//
// NOTE: the 1st eTtrg "bin" is defined to be
//       the entire eTtrg range.
//
// NOTE: minimum areas for different R:
//   R = 0.3 -- 0.2   R = 0.5 -- 0.65
//   R = 0.4 -- 0.35  R = 0.7 -- 1.2


#include <vector>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;


// global constants
static const UInt_t  NTypes(2);
static const UInt_t  NBinsEt(4);
static const UInt_t  NBinsId(2);
static const UInt_t  NBinsInitial(1100);
static const Float_t MaxError(0.30);
static const Float_t MinWidth(1.0);
static const Float_t MinPtJet(-10.);
static const Float_t MaxPtJet(100.);
static const Float_t MinEtTrg[NBinsEt]  = {9., 9., 11., 15.};
static const Float_t MaxEtTrg[NBinsEt]  = {20., 11., 15., 20.};
static const Float_t MinTspTrg[NBinsId] = {0., 0.2};
static const Float_t MaxTspTrg[NBinsId] = {0.08, 0.6};



void OptimizeBinning(const Bool_t isInBatchMode=false) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning binning optimization..." << endl;

  // io parameters
  const TString sOutput("optimizedBins.pp200r9.r03a02rm1full.d23m3y2018.root");
  const TString sTxtOut("optimizedBins.pp200r9.r03a02rm1full.d23m3y2018.txt");
  const TString sInput("output/CollaborationMeetingJan2018/pp200r9.withOaCones.eTtrg920.r03rm1full.d21m1y2018.root");
  const TString sTree("JetTree");

  // jet parameters
  const UInt_t   type(1);
  const TString  sRes("0.3");
  const Double_t rJet(0.3);
  const Double_t hTrgMax(0.9);
  const Double_t aJetMin(0.2);
  const Double_t pTjetMin(0.2);
  const Double_t dFrecoil(TMath::PiOver4());


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
  Int_t    EventIndex = 0;
  Int_t    NJets      = 0;
  Double_t Refmult    = 0;
  Double_t TSP        = 0;
  Double_t TrgEta     = 0;
  Double_t TrgPhi     = 0;
  Double_t TrgEt      = 0;
  Double_t Rho        = 0;
  Double_t Sigma      = 0;
  Double_t Vz         = 0;
  // jet leaves
  vector<Double_t> *JetPt            = 0;
  vector<Double_t> *JetNCons         = 0;
  vector<Double_t> *JetIndex         = 0;
  vector<Double_t> *JetPtCorr        = 0;
  vector<Double_t> *JetEta           = 0;
  vector<Double_t> *JetPhi           = 0;
  vector<Double_t> *JetE             = 0;
  vector<Double_t> *JetArea          = 0;
  vector<Double_t> *JetPtOffAxisUp   = 0;
  vector<Double_t> *JetPtOffAxisDown = 0;
  // constituent leaves
  vector<vector<Double_t> > *JetConsPt  = 0;
  vector<vector<Double_t> > *JetConsEta = 0;
  vector<vector<Double_t> > *JetConsPhi = 0;
  vector<vector<Double_t> > *JetConsE   = 0;

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


  // declare initial histograms
  TH1D *hJetPtInitial[NBinsEt][NBinsId];
  hJetPtInitial[0][0] = new TH1D("hJetPtInitial_pi920", "", NBinsInitial, MinPtJet, MaxPtJet);
  hJetPtInitial[1][0] = new TH1D("hJetPtInitial_pi911", "", NBinsInitial, MinPtJet, MaxPtJet);
  hJetPtInitial[2][0] = new TH1D("hJetPtInitial_pi1115", "", NBinsInitial, MinPtJet, MaxPtJet);
  hJetPtInitial[3][0] = new TH1D("hJetPtInitial_pi1520", "", NBinsInitial, MinPtJet, MaxPtJet);
  hJetPtInitial[0][1] = new TH1D("hJetPtInitial_ga920", "", NBinsInitial, MinPtJet, MaxPtJet);
  hJetPtInitial[1][1] = new TH1D("hJetPtInitial_ga911", "", NBinsInitial, MinPtJet, MaxPtJet);
  hJetPtInitial[2][1] = new TH1D("hJetPtInitial_ga1115", "", NBinsInitial, MinPtJet, MaxPtJet);
  hJetPtInitial[3][1] = new TH1D("hJetPtInitial_ga1520", "", NBinsInitial, MinPtJet, MaxPtJet);
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      hJetPtInitial[iBinEt][iBinId] -> Sumw2();
    }
  }
  cout << "    Declared initial histograms." << endl;


  const UInt_t nEvts = tInput -> GetEntries();
  cout << "    Beginning initial event loop: " << nEvts << " events to process..." << endl;

  // initial event loop
  Int_t  bytes(0);
  UInt_t nBytes(0);
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    // load entry
    bytes   = tInput -> GetEntry(iEvt);
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
    const Bool_t isInTrgEtCut  = ((eTtrg > MinEtTrg[0]) && (eTtrg < MaxEtTrg[0]));
    const Bool_t isInPi0tspCut = ((TSP > MinTspTrg[0]) && (TSP < MaxTspTrg[0]));
    const Bool_t isInGamTspCut = ((TSP > MinTspTrg[1]) && (TSP < MaxTspTrg[1]));
    const Bool_t isInTspCut    = (isInPi0tspCut || isInGamTspCut);
    if (!isInTrgEtaCut || !isInTrgEtCut || !isInTspCut) continue;


    // determine trigger bin
    UInt_t eTbin(0);
    for (UInt_t iBinEt = 1; iBinEt < NBinsEt; iBinEt++) {
      const Bool_t isInEtBin = ((eTtrg >= MinEtTrg[iBinEt]) && (eTtrg < MaxEtTrg[iBinEt]));
      if (isInEtBin) {
        eTbin = iBinEt;
        break;
      }
    }

    UInt_t idBin(0);
    if (isInPi0tspCut) idBin = 0;
    if (isInGamTspCut) idBin = 1;


    // jet loop
    for (UInt_t iJet = 0; iJet < nJets; iJet++) {

      // jet info
      const Double_t fJet   = JetPhi    -> at(iJet);
      const Double_t aJet   = JetArea   -> at(iJet);
      const Double_t pTjet  = JetPt     -> at(iJet);
      const Double_t pTreco = JetPtCorr -> at(iJet);

      Double_t dFjet = fJet - fTrg;
      if (dFjet < 0.)
        dFjet += TMath::TwoPi();
      if (dFjet > TMath::TwoPi())
        dFjet -= TMath::TwoPi();


      // jet cuts
      const Bool_t isInAreaCut = (aJet > aJetMin);
      const Bool_t isInPtCut   = (pTjet > pTjetMin);
      const Bool_t isRecoil    = (TMath::Abs(dFjet - TMath::Pi()) < dFrecoil);
      if (!isInAreaCut || !isInPtCut || !isRecoil) continue;

      hJetPtInitial[0][idBin]     -> Fill(pTreco);
      hJetPtInitial[eTbin][idBin] -> Fill(pTreco);

    }  // end jet loop
  }  // end initial event loop

  cout << "    Initial event loop finished." << endl;


  // declare arrays
  UInt_t  nOptimal[NBinsEt][NBinsId];
  Float_t lowEdges[NBinsEt][NBinsId][NBinsInitial];
  Float_t highEdges[NBinsEt][NBinsId][NBinsInitial];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      nOptimal[iBinEt][iBinId] = 0;
      for (UInt_t iBinPt = 0; iBinPt < NBinsInitial; iBinPt++) {
        lowEdges[iBinEt][iBinId][iBinPt]  = 0.;
        highEdges[iBinEt][iBinId][iBinPt] = 0.;
      }  // end pT loop
    }  // end id loop
  }  // end eT loop

  // optimize binning
  UInt_t iOptimal  = 0;
  UInt_t iBinStart = NBinsInitial;
  UInt_t nBinsLeft = NBinsInitial;
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      iOptimal  = 0;
      iBinStart = NBinsInitial;
      nBinsLeft = NBinsInitial;
      while (nBinsLeft > 0) {

        // initialize optimal bin
        Float_t number   = 0.;
        Float_t lowSide  = hJetPtInitial[iBinEt][iBinId] -> GetBinLowEdge(iBinStart);
        Float_t highSide = hJetPtInitial[iBinEt][iBinId] -> GetBinLowEdge(iBinStart + 1);

        // determine size of optimal bin
        for (UInt_t iBinPt = iBinStart; iBinPt > 0; iBinPt--) {
          nBinsLeft--;
          if (number != 0.) {
            Double_t error  = 1. / TMath::Sqrt(number);
            Double_t width  = highSide - lowSide;
            Bool_t   isGood = ((error < MaxError) && (width > MinWidth));
            if (isGood) {
              lowSide   = hJetPtInitial[iBinEt][iBinId] -> GetBinLowEdge(iBinPt);
              iBinStart = iBinPt - 1;
              break;
            }
            else {
              number  += hJetPtInitial[iBinEt][iBinId] -> GetBinContent(iBinPt);
              lowSide  = hJetPtInitial[iBinEt][iBinId] -> GetBinLowEdge(iBinPt);
            }
          }
          else {
            number += hJetPtInitial[iBinEt][iBinId] -> GetBinContent(iBinPt);
          }
        }  // end initial bin loop

        // fill arrays
        lowEdges[iBinEt][iBinId][iOptimal]  = lowSide;
        highEdges[iBinEt][iBinId][iOptimal] = highSide;
        nOptimal[iBinEt][iBinId]++;
        iOptimal++;

      }  // end while loop
    }  // end id loop
  }  // end eT loop
  cout << "    Optimized binning." << endl;


  // declare arrays
  Float_t lowEdgesPi0[nOptimal[0][0] + 1];
  Float_t lowEdgesPi1[nOptimal[1][0] + 1];
  Float_t lowEdgesPi2[nOptimal[2][0] + 1];
  Float_t lowEdgesPi3[nOptimal[3][0] + 1];
  Float_t lowEdgesGa0[nOptimal[0][0] + 1];
  Float_t lowEdgesGa1[nOptimal[1][0] + 1];
  Float_t lowEdgesGa2[nOptimal[2][0] + 1];
  Float_t lowEdgesGa3[nOptimal[3][0] + 1];

  // open output stream
  ofstream txtOutput(sTxtOut.Data());
  if (!txtOutput) {
    cerr << "PANIC: couldn't open output stream!" << endl;
    return;
  }

  // sort low and high edges, stream to output
  UInt_t iBinOptSort  = 0;
  UInt_t iBinOptStart = 0;
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      iBinOptSort  = 0;
      iBinOptStart = nOptimal[iBinEt][iBinId] - 1;
      for (UInt_t iBinOptimal = iBinOptStart; iBinOptimal > 0; iBinOptimal--) {

        // get bin edges
        const Float_t lo = lowEdges[iBinEt][iBinId][iBinOptimal];
        const Float_t hi = highEdges[iBinEt][iBinId][iBinOptimal];
        switch (iBinEt) {
          case 0:
            if (iBinId == 0) lowEdgesPi0[iBinOptSort] = lo;
            if (iBinId == 1) lowEdgesGa0[iBinOptSort] = lo;
            break;
          case 1:
            if (iBinId == 0) lowEdgesPi1[iBinOptSort] = lo;
            if (iBinId == 1) lowEdgesGa1[iBinOptSort] = lo;
            break;
          case 2:
            if (iBinId == 0) lowEdgesPi2[iBinOptSort] = lo;
            if (iBinId == 1) lowEdgesGa2[iBinOptSort] = lo;
            break;
          case 3:
            if (iBinId == 0) lowEdgesPi3[iBinOptSort] = lo;
            if (iBinId == 1) lowEdgesGa3[iBinOptSort] = lo;
            break;
        }

        // stream info
        txtOutput << iBinEt;
        txtOutput << " ";
        txtOutput << iBinId;
        txtOutput << " ";
        txtOutput << iBinOptSort;
        txtOutput << " ";
        txtOutput << lo;
        txtOutput << " ";
        txtOutput << hi;
        txtOutput << endl;
        iBinOptSort++;

        // for histograms
        if (iBinOptSort == nOptimal[iBinEt][iBinId] - 1) {
          const UInt_t iHistLabel  = iBinEt + (iBinId * 10);
          const UInt_t iBinMaxSort = iBinOptSort + 1;
          switch (iHistLabel) {
            case 0:
              lowEdgesPi0[iBinOptSort] = hi;
              lowEdgesPi0[iBinMaxSort] = MaxPtJet;
              break;
            case 1:
              lowEdgesPi1[iBinOptSort] = hi;
              lowEdgesPi1[iBinMaxSort] = MaxPtJet;
              break;
            case 2:
              lowEdgesPi2[iBinOptSort] = hi;
              lowEdgesPi2[iBinMaxSort] = MaxPtJet;
              break;
            case 3:
              lowEdgesPi3[iBinOptSort] = hi;
              lowEdgesPi3[iBinMaxSort] = MaxPtJet;
              break;
            case 10:
              lowEdgesGa0[iBinOptSort] = hi;
              lowEdgesGa0[iBinMaxSort] = MaxPtJet;
              break;
            case 11:
              lowEdgesGa1[iBinOptSort] = hi;
              lowEdgesGa1[iBinMaxSort] = MaxPtJet;
              break;
            case 12:
              lowEdgesGa2[iBinOptSort] = hi;
              lowEdgesGa2[iBinMaxSort] = MaxPtJet;
              break;
            case 13:
              lowEdgesGa3[iBinOptSort] = hi;
              lowEdgesGa3[iBinMaxSort] = MaxPtJet;
              break;
          }  // end switch
          txtOutput << iBinEt;
          txtOutput << " ";
          txtOutput << iBinId;
          txtOutput << " ";
          txtOutput << iBinOptSort;
          txtOutput << " ";
          txtOutput << hi;
          txtOutput << " ";
          txtOutput << MaxPtJet;
          txtOutput << endl;
        }  // end if
      }  // end bin loop
    }  // end id loop
  }  // end eT loop
  cout << "    Sorted low edges and created text file." << endl;


  // declare optimized histograms
  TH1D *hJetPtOptimal[NBinsEt][NBinsId];
  hJetPtOptimal[0][0] = new TH1D("hJetPtOptimal_pi920", "", nOptimal[0][0], lowEdgesPi0);
  hJetPtOptimal[1][0] = new TH1D("hJetPtOptimal_pi911", "", nOptimal[1][0], lowEdgesPi1);
  hJetPtOptimal[2][0] = new TH1D("hJetPtOptimal_pi1115", "", nOptimal[2][0], lowEdgesPi2);
  hJetPtOptimal[3][0] = new TH1D("hJetPtOptimal_pi1520", "", nOptimal[3][0], lowEdgesPi3);
  hJetPtOptimal[0][1] = new TH1D("hJetPtOptimal_ga920", "", nOptimal[0][1], lowEdgesGa0);
  hJetPtOptimal[1][1] = new TH1D("hJetPtOptimal_ga911", "", nOptimal[1][1], lowEdgesGa1);
  hJetPtOptimal[2][1] = new TH1D("hJetPtOptimal_ga1115", "", nOptimal[2][1], lowEdgesGa2);
  hJetPtOptimal[3][1] = new TH1D("hJetPtOptimal_ga1520", "", nOptimal[3][1], lowEdgesGa3);
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      hJetPtOptimal[iBinEt][iBinId] -> Sumw2();
    }
  }
  cout << "    Declared optimized histograms." << endl;


  UInt_t nTrg[NBinsEt][NBinsId];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      nTrg[iBinEt][iBinId] = 0;
    }
  }
  bytes = 0;
  cout << "    Beginning optimized event loop: " << nEvts << " events to process." << endl;

  // optimized event loop
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
    const Bool_t isInTrgEtCut  = ((eTtrg > MinEtTrg[0]) && (eTtrg < MaxEtTrg[0]));
    const Bool_t isInPi0tspCut = ((TSP > MinTspTrg[0]) && (TSP < MaxTspTrg[0]));
    const Bool_t isInGamTspCut = ((TSP > MinTspTrg[1]) && (TSP < MaxTspTrg[1]));
    const Bool_t isInTspCut    = (isInPi0tspCut || isInGamTspCut);
    if (!isInTrgEtaCut || !isInTrgEtCut || !isInTspCut) continue;


    // determine trigger bin
    UInt_t eTbin(0);
    for (UInt_t iBinEt = 1; iBinEt < NBinsEt; iBinEt++) {
      const Bool_t isInEtBin = ((eTtrg >= MinEtTrg[iBinEt]) && (eTtrg < MaxEtTrg[iBinEt]));
      if (isInEtBin) {
        eTbin = iBinEt;
        break;
      }
    }

    UInt_t idBin(0);
    if (isInPi0tspCut) idBin = 0;
    if (isInGamTspCut) idBin = 1;

    // increment counters
    nTrg[0][idBin]++;
    nTrg[eTbin][idBin]++;


    // jet loop
    for (UInt_t iJet = 0; iJet < nJets; iJet++) {

      // jet info
      const Double_t fJet   = JetPhi    -> at(iJet);
      const Double_t aJet   = JetArea   -> at(iJet);
      const Double_t pTjet  = JetPt     -> at(iJet);
      const Double_t pTreco = JetPtCorr -> at(iJet);

      Double_t dFjet = fJet - fTrg;
      if (dFjet < 0.)
        dFjet += TMath::TwoPi();
      if (dFjet > TMath::TwoPi())
        dFjet -= TMath::TwoPi();


      // jet cuts
      const Bool_t isInAreaCut = (aJet > aJetMin);
      const Bool_t isInPtCut   = (pTjet > pTjetMin);
      const Bool_t isRecoil    = (TMath::Abs(dFjet - TMath::Pi()) < dFrecoil);
      if (!isInAreaCut || !isInPtCut || !isRecoil) continue;

      hJetPtOptimal[0][idBin]     -> Fill(pTreco);
      hJetPtOptimal[eTbin][idBin] -> Fill(pTreco);

    }  // end jet loop
  }  // end optimized event loop

  cout << "    Optimized event loop finished." << endl;


  // normalize histograms
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      const UInt_t  nIni    = hJetPtInitial[iBinEt][iBinId] -> GetNbinsX();
      const UInt_t  nOpt    = hJetPtOptimal[iBinEt][iBinId] -> GetNbinsX();
      const Float_t trgNorm = (Float_t) nTrg[iBinEt][iBinId];
      const Float_t etaNorm = 2. * (1. - rJet);
      const Float_t normer  = trgNorm * etaNorm;
      hJetPtInitial[iBinEt][iBinId] -> Scale(1. / normer);
      hJetPtOptimal[iBinEt][iBinId] -> Scale(1. / normer);
      for (UInt_t iIni = 1; iIni <= nIni; iIni++) {
        const Double_t iniBin = hJetPtInitial[iBinEt][iBinId] -> GetBinWidth(iIni);
        const Double_t iniVal = hJetPtInitial[iBinEt][iBinId] -> GetBinContent(iIni);
        const Double_t iniErr = hJetPtInitial[iBinEt][iBinId] -> GetBinError(iIni);
        const Double_t newVal = iniVal / iniBin;
        const Double_t newErr = iniErr / iniBin;
        hJetPtInitial[iBinEt][iBinId] -> SetBinContent(iIni, newVal);
        hJetPtInitial[iBinEt][iBinId] -> SetBinError(iIni, newErr);
      }
      for (UInt_t iOpt = 1; iOpt < nOpt; iOpt++) {
        const Double_t optBin = hJetPtOptimal[iBinEt][iBinId] -> GetBinWidth(iOpt);
        const Double_t optVal = hJetPtOptimal[iBinEt][iBinId] -> GetBinContent(iOpt);
        const Double_t optErr = hJetPtOptimal[iBinEt][iBinId] -> GetBinError(iOpt);
        const Double_t newVal = optVal / optBin;
        const Double_t newErr = optErr / optBin;
        hJetPtOptimal[iBinEt][iBinId] -> SetBinContent(iOpt, newVal);
        hJetPtOptimal[iBinEt][iBinId] -> SetBinError(iOpt, newErr);
      }  // end bin loop
    }  // end id loop
  }  // end et loop
  cout << "    Normalized histograms." << endl;


  // set styles
  const UInt_t  fColIni(1);
  const UInt_t  fColOpt(2);
  const UInt_t  fMarIni(1);
  const UInt_t  fMarOpt(4);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.03);
  const Float_t fOffX(0.9);
  const Float_t fOffY(1.1);
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - #rho #upoint A^{jet}");
  const TString sTitleY("(1/N^{trg}) dN^{jet}/(dp_{T}^{reco} d#eta^{jet})");
  const TString sBase[NBinsId] = {" recoil jets, #pi^{0} (R = ", " recoil jets, #gamma^{rich} (R = "};
  const TString sType[NTypes]  = {"Charged", "Full"};
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      // make title
      TString sTitle(sType[type].Data());
      sTitle.Append(sBase[iBinId].Data());
      sTitle.Append(sRes.Data());
      sTitle.Append(")");

      // initial histograms
      hJetPtInitial[iBinEt][iBinId] -> SetLineColor(fColIni);
      hJetPtInitial[iBinEt][iBinId] -> SetMarkerColor(fColIni);
      hJetPtInitial[iBinEt][iBinId] -> SetMarkerStyle(fMarIni);
      hJetPtInitial[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hJetPtInitial[iBinEt][iBinId] -> SetTitle(sTitle.Data());
      hJetPtInitial[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetPtInitial[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetPtInitial[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hJetPtInitial[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hJetPtInitial[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hJetPtInitial[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetPtInitial[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetPtInitial[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hJetPtInitial[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hJetPtInitial[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      // optimized histograms
      hJetPtOptimal[iBinEt][iBinId] -> SetLineColor(fColOpt);
      hJetPtOptimal[iBinEt][iBinId] -> SetMarkerColor(fColOpt);
      hJetPtOptimal[iBinEt][iBinId] -> SetMarkerStyle(fMarOpt);
      hJetPtOptimal[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hJetPtOptimal[iBinEt][iBinId] -> SetTitle(sTitle.Data());
      hJetPtOptimal[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetPtOptimal[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetPtOptimal[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hJetPtOptimal[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hJetPtOptimal[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hJetPtOptimal[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetPtOptimal[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetPtOptimal[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hJetPtOptimal[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hJetPtOptimal[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);

    }  // end id loop
  }  // end et loop
  cout << "    Set styles." << endl;


  // make legends
  const UInt_t  fColLeg(0);
  const UInt_t  fLinLeg(1);
  const UInt_t  fAlign(12);
  const Float_t xyLegend[4]   = {0.1, 0.1, 0.3, 0.3};
  const TString sEntry[2]     = {"initial binning", "optimized binning"};
  const TString sTrg[NBinsEt] = {"E_{T}^{trg}#in(9,20) GeV/c", "E_{T}^{trg}#in(9,11) GeV/c", "E_{T}^{trg}#in(11,15) GeV/c", "E_{T}^{trg}#in(15,20) GeV/c"};

  TLegend *legend[NBinsEt][NBinsId];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      legend[iBinEt][iBinId] = new TLegend(xyLegend[0], xyLegend[1], xyLegend[2], xyLegend[3], sTrg[iBinEt].Data());
      legend[iBinEt][iBinId] -> SetFillColor(fColLeg);
      legend[iBinEt][iBinId] -> SetLineColor(fLinLeg);
      legend[iBinEt][iBinId] -> SetTextFont(fTxt);
      legend[iBinEt][iBinId] -> SetTextAlign(fAlign);
      legend[iBinEt][iBinId] -> AddEntry(hJetPtInitial[iBinEt][iBinId], sEntry[0].Data());
      legend[iBinEt][iBinId] -> AddEntry(hJetPtOptimal[iBinEt][iBinId], sEntry[1].Data());
    }
  }
  cout << "    Made legend." << endl;


  // save histograms
  fOutput -> cd();
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      hJetPtInitial[iBinEt][iBinId] -> Write();
      hJetPtOptimal[iBinEt][iBinId] -> Write();
    }
  }
  cout << "    Saved histograms." << endl;

  // make plots
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  log(1);

  TCanvas *canvas[NBinsEt][NBinsId];
  fOutput -> cd();
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      const UInt_t iHistLabel = iBinEt + (iBinId * 10);
      switch (iHistLabel) {
        case 0:
          canvas[iBinEt][iBinId] = new TCanvas("cPi920", "", width, height);
          break;
        case 1:
          canvas[iBinEt][iBinId] = new TCanvas("cPi911", "", width, height);
          break;
        case 2:
          canvas[iBinEt][iBinId] = new TCanvas("cPi1115", "", width, height);
          break;
        case 3:
          canvas[iBinEt][iBinId] = new TCanvas("cPi1520", "", width, height);
          break;
        case 10:
          canvas[iBinEt][iBinId] = new TCanvas("cGa920", "", width, height);
          break;
        case 11:
          canvas[iBinEt][iBinId] = new TCanvas("cGa911", "", width, height);
          break;
        case 12:
          canvas[iBinEt][iBinId] = new TCanvas("cGa1115", "", width, height);
          break;
        case 13:
          canvas[iBinEt][iBinId] = new TCanvas("cGa1520", "", width, height);
          break;
      }
      canvas[iBinEt][iBinId]        -> SetGrid(grid, grid);
      canvas[iBinEt][iBinId]        -> SetLogy(log);
      hJetPtInitial[iBinEt][iBinId] -> Draw();
      hJetPtOptimal[iBinEt][iBinId] -> Draw("same");
      legend[iBinEt][iBinId]        -> Draw();
      canvas[iBinEt][iBinId]        -> Write();
      canvas[iBinEt][iBinId]        -> Close();
    }
  }
  cout << "    Made plots." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Optimization finished!\n" << endl;

}

// End ------------------------------------------------------------------------
