// 'UncorrelatedJetStudy.C'
// Derek Anderson
// 01.21.2018
//
// This takes the output of 'StJetTreeMaker' and
// produces a series of histograms comparing
// recoil jets to so-called "uncorrelated jets."
//
// NOTE: the 1st eTtrg "bin" is defined to be
//       the entire eTtrg range.
//
// NOTE: minimum areas for different R:
//   R = 0.3 -- 0.2   R = 0.5 -- 0.65
//   R = 0.4 -- 0.35  R = 0.7 -- 1.2
//
// NOTE: 'JetType' indicates charged or full.
//   JetType = 0 -- charged
//   JetType = 1 -- full


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
static const UInt_t NBinsEt = 4;
static const UInt_t NBinsId = 2;
static const UInt_t NBinsDf = 3;
static const UInt_t NBinsPt = 35;
static const UInt_t NPadsPt = 2;
static const UInt_t NCones  = 2;
static const UInt_t JetType = 1;



void UncorrelatedJetStudy(const Bool_t isInBatchMode=false) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning uncorrelated jet study..." << endl;

  // io parameters
  const TString sOutput("pp200r9.ueStudy.eTtrg920.r03a02rm1full.d23m1y2018.root");
  const TString sInput("pp200r9.withOaCones.eTtrg920.r03rm1full.d21m1y2018.root");
  const TString sTree("JetTree");

  // trigger parameters
  const Double_t hTrgMax(0.9);
  const Double_t pi0tsp[2] = {0., 0.08};
  const Double_t gamTsp[2] = {0.2, 0.6};

  // jet parameters
  const Double_t rJet(0.3);
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
  vector<Double_t> *JetEta;
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
  TH1D *hDfJet[NBinsEt][NBinsId];
  TH1D *hDfRE[NBinsEt][NBinsId];
  TH1D *hDfUE[NBinsEt][NBinsId][NBinsDf];
  TH1D *hPtRE[NBinsEt][NBinsId];
  TH1D *hPtUE[NBinsEt][NBinsId][NBinsDf];
  TH1D *hPtOA[NBinsEt][NBinsId][NCones + 1];
  TH2D *hPtUpVsDown[NBinsEt][NBinsId];
  TH2D *hPtAvgVsRE[NBinsEt][NBinsId];

  const UInt_t   nDf(60);
  const UInt_t   nPt(NBinsPt - 1);
  const UInt_t   nPtOA(52);
  const Double_t dFbin[2]       = {0., TMath::TwoPi()};
  const Double_t pTbin[NBinsPt] = {-5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 30., 35., 40., 45., 50.};
  const Double_t pToaBin[2]     = {-1., 25.};
  // delta-phi plots
  hDfJet[0][0] = new TH1D("hDfJet_pi920", "", nDf, dFbin[0], dFbin[1]);
  hDfJet[0][1] = new TH1D("hDfJet_ga920", "", nDf, dFbin[0], dFbin[1]);
  hDfJet[1][0] = new TH1D("hDfJet_pi911", "", nDf, dFbin[0], dFbin[1]);
  hDfJet[1][1] = new TH1D("hDfJet_ga911", "", nDf, dFbin[0], dFbin[1]);
  hDfJet[2][0] = new TH1D("hDfJet_pi1115", "", nDf, dFbin[0], dFbin[1]);
  hDfJet[2][1] = new TH1D("hDfJet_ga1115", "", nDf, dFbin[0], dFbin[1]);
  hDfJet[3][0] = new TH1D("hDfJet_pi1520", "", nDf, dFbin[0], dFbin[1]);
  hDfJet[3][1] = new TH1D("hDfJet_ga1520", "", nDf, dFbin[0], dFbin[1]);
  // recoil delta-phi plots
  hDfRE[0][0] = new TH1D("hDfRecoil_pi920", "", nDf, dFbin[0], dFbin[1]);
  hDfRE[0][1] = new TH1D("hDfRecoil_ga920", "", nDf, dFbin[0], dFbin[1]);
  hDfRE[1][0] = new TH1D("hDfRecoil_pi911", "", nDf, dFbin[0], dFbin[1]);
  hDfRE[1][1] = new TH1D("hDfRecoil_ga911", "", nDf, dFbin[0], dFbin[1]);
  hDfRE[2][0] = new TH1D("hDfRecoil_pi1115", "", nDf, dFbin[0], dFbin[1]);
  hDfRE[2][1] = new TH1D("hDfRecoil_ga1115", "", nDf, dFbin[0], dFbin[1]);
  hDfRE[3][0] = new TH1D("hDfRecoil_pi1520", "", nDf, dFbin[0], dFbin[1]);
  hDfRE[3][1] = new TH1D("hDfRecoil_ga1520", "", nDf, dFbin[0], dFbin[1]);
  // uncorrelated delta-phi plots
  hDfUE[0][0][0] = new TH1D("hDfUncorrelated_pi920df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[0][0][1] = new TH1D("hDfUncorrelated_pi920df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[0][0][2] = new TH1D("hDfUncorrelated_pi920df2", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[0][1][0] = new TH1D("hDfUncorrelated_ga920df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[0][1][1] = new TH1D("hDfUncorrelated_ga920df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[0][1][2] = new TH1D("hDfUncorrelated_ga920df2", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[1][0][0] = new TH1D("hDfUncorrelated_pi911df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[1][0][1] = new TH1D("hDfUncorrelated_pi911df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[1][0][2] = new TH1D("hDfUncorrelated_pi911df2", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[1][1][0] = new TH1D("hDfUncorrelated_ga911df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[1][1][1] = new TH1D("hDfUncorrelated_ga911df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[1][1][2] = new TH1D("hDfUncorrelated_ga911df2", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[2][0][0] = new TH1D("hDfUncorrelated_pi1115df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[2][0][1] = new TH1D("hDfUncorrelated_pi1115df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[2][0][2] = new TH1D("hDfUncorrelated_pi1115df2", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[2][1][0] = new TH1D("hDfUncorrelated_ga1115df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[2][1][1] = new TH1D("hDfUncorrelated_ga1115df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[2][1][2] = new TH1D("hDfUncorrelated_ga1115df2", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[3][0][0] = new TH1D("hDfUncorrelated_pi1520df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[3][0][1] = new TH1D("hDfUncorrelated_pi1520df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[3][0][2] = new TH1D("hDfUncorrelated_pi1520df2", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[3][1][0] = new TH1D("hDfUncorrelated_ga1520df0", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[3][1][1] = new TH1D("hDfUncorrelated_ga1520df1", "", nDf, dFbin[0], dFbin[1]);
  hDfUE[3][1][2] = new TH1D("hDfUncorrelated_ga1520df2", "", nDf, dFbin[0], dFbin[1]);
  // recoil pT plots
  hPtRE[0][0]  = new TH1D("hPtRecoil_pi920", "", nPt, pTbin);
  hPtRE[0][1]  = new TH1D("hPtRecoil_ga920", "", nPt, pTbin);
  hPtRE[1][0]  = new TH1D("hPtRecoil_pi911", "", nPt, pTbin);
  hPtRE[1][1]  = new TH1D("hPtRecoil_ga911", "", nPt, pTbin);
  hPtRE[2][0]  = new TH1D("hPtRecoil_pi1115", "", nPt, pTbin);
  hPtRE[2][1]  = new TH1D("hPtRecoil_ga1115", "", nPt, pTbin);
  hPtRE[3][0]  = new TH1D("hPtRecoil_pi1520", "", nPt, pTbin);
  hPtRE[3][1]  = new TH1D("hPtRecoil_ga1520", "", nPt, pTbin);
  // uncorrelated pT plots
  hPtUE[0][0][0] = new TH1D("hPtUncorrelated_pi920df0", "", nPt, pTbin);
  hPtUE[0][0][1] = new TH1D("hPtUncorrelated_pi920df1", "", nPt, pTbin);
  hPtUE[0][0][2] = new TH1D("hPtUncorrelated_pi920df2", "", nPt, pTbin);
  hPtUE[0][1][0] = new TH1D("hPtUncorrelated_ga920df0", "", nPt, pTbin);
  hPtUE[0][1][1] = new TH1D("hPtUncorrelated_ga920df1", "", nPt, pTbin);
  hPtUE[0][1][2] = new TH1D("hPtUncorrelated_ga920df2", "", nPt, pTbin);
  hPtUE[1][0][0] = new TH1D("hPtUncorrelated_pi911df0", "", nPt, pTbin);
  hPtUE[1][0][1] = new TH1D("hPtUncorrelated_pi911df1", "", nPt, pTbin);
  hPtUE[1][0][2] = new TH1D("hPtUncorrelated_pi911df2", "", nPt, pTbin);
  hPtUE[1][1][0] = new TH1D("hPtUncorrelated_ga911df0", "", nPt, pTbin);
  hPtUE[1][1][1] = new TH1D("hPtUncorrelated_ga911df1", "", nPt, pTbin);
  hPtUE[1][1][2] = new TH1D("hPtUncorrelated_ga911df2", "", nPt, pTbin);
  hPtUE[2][0][0] = new TH1D("hPtUncorrelated_pi1115df0", "", nPt, pTbin);
  hPtUE[2][0][1] = new TH1D("hPtUncorrelated_pi1115df1", "", nPt, pTbin);
  hPtUE[2][0][2] = new TH1D("hPtUncorrelated_pi1115df2", "", nPt, pTbin);
  hPtUE[2][1][0] = new TH1D("hPtUncorrelated_ga1115df0", "", nPt, pTbin);
  hPtUE[2][1][1] = new TH1D("hPtUncorrelated_ga1115df1", "", nPt, pTbin);
  hPtUE[2][1][2] = new TH1D("hPtUncorrelated_ga1115df2", "", nPt, pTbin);
  hPtUE[3][0][0] = new TH1D("hPtUncorrelated_pi1520df0", "", nPt, pTbin);
  hPtUE[3][0][1] = new TH1D("hPtUncorrelated_pi1520df1", "", nPt, pTbin);
  hPtUE[3][0][2] = new TH1D("hPtUncorrelated_pi1520df2", "", nPt, pTbin);
  hPtUE[3][1][0] = new TH1D("hPtUncorrelated_ga1520df0", "", nPt, pTbin);
  hPtUE[3][1][1] = new TH1D("hPtUncorrelated_ga1520df1", "", nPt, pTbin);
  hPtUE[3][1][2] = new TH1D("hPtUncorrelated_ga1520df2", "", nPt, pTbin);
  // off-axis pT plots
  hPtOA[0][0][0] = new TH1D("hPtOffAxis_pi920up", "", nPt, pTbin);
  hPtOA[0][0][1] = new TH1D("hPtOffAxis_pi920down", "", nPt, pTbin);
  hPtOA[0][0][2] = new TH1D("hPtOffAxis_pi920avg", "", nPt, pTbin);
  hPtOA[0][1][0] = new TH1D("hPtOffAxis_ga920up", "", nPt, pTbin);
  hPtOA[0][1][1] = new TH1D("hPtOffAxis_ga920down", "", nPt, pTbin);
  hPtOA[0][1][2] = new TH1D("hPtOffAxis_ga920avg", "", nPt, pTbin);
  hPtOA[1][0][0] = new TH1D("hPtOffAxis_pi911up", "", nPt, pTbin);
  hPtOA[1][0][1] = new TH1D("hPtOffAxis_pi911down", "", nPt, pTbin);
  hPtOA[1][0][2] = new TH1D("hPtOffAxis_pi911avg", "", nPt, pTbin);
  hPtOA[1][1][0] = new TH1D("hPtOffAxis_ga911up", "", nPt, pTbin);
  hPtOA[1][1][1] = new TH1D("hPtOffAxis_ga911down", "", nPt, pTbin);
  hPtOA[1][1][2] = new TH1D("hPtOffAxis_ga911avg", "", nPt, pTbin);
  hPtOA[2][0][0] = new TH1D("hPtOffAxis_pi1115up", "", nPt, pTbin);
  hPtOA[2][0][1] = new TH1D("hPtOffAxis_pi1115down", "", nPt, pTbin);
  hPtOA[2][0][2] = new TH1D("hPtOffAxis_pi1115avg", "", nPt, pTbin);
  hPtOA[2][1][0] = new TH1D("hPtOffAxis_ga1115up", "", nPt, pTbin);
  hPtOA[2][1][1] = new TH1D("hPtOffAxis_ga1115down", "", nPt, pTbin);
  hPtOA[2][1][2] = new TH1D("hPtOffAxis_ga1115avg", "", nPt, pTbin);
  hPtOA[3][0][0] = new TH1D("hPtOffAxis_pi1520up", "", nPt, pTbin);
  hPtOA[3][0][1] = new TH1D("hPtOffAxis_pi1520down", "", nPt, pTbin);
  hPtOA[3][0][2] = new TH1D("hPtOffAxis_pi1520avg", "", nPt, pTbin);
  hPtOA[3][1][0] = new TH1D("hPtOffAxis_ga1520up", "", nPt, pTbin);
  hPtOA[3][1][1] = new TH1D("hPtOffAxis_ga1520down", "", nPt, pTbin);
  hPtOA[3][1][2] = new TH1D("hPtOffAxis_ga1520avg", "", nPt, pTbin);
  // off-axis 2d plots
  hPtUpVsDown[0][0] = new TH2D("hPtUpVsDown_pi920", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtUpVsDown[0][1] = new TH2D("hPtUpVsDown_ga920", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtUpVsDown[1][0] = new TH2D("hPtUpVsDown_pi911", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtUpVsDown[1][1] = new TH2D("hPtUpVsDown_ga911", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtUpVsDown[2][0] = new TH2D("hPtUpVsDown_pi1115", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtUpVsDown[2][1] = new TH2D("hPtUpVsDown_ga1115", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtUpVsDown[3][0] = new TH2D("hPtUpVsDown_pi1520", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtUpVsDown[3][1] = new TH2D("hPtUpVsDown_ga1520", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  // off-axis vs. recoil plots
  hPtAvgVsRE[0][0] = new TH2D("hPtRecoilVsAvg_pi920", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtAvgVsRE[0][1] = new TH2D("hPtRecoilVsAvg_ga920", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtAvgVsRE[1][0] = new TH2D("hPtRecoilVsAvg_pi911", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtAvgVsRE[1][1] = new TH2D("hPtRecoilVsAvg_ga911", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtAvgVsRE[2][0] = new TH2D("hPtRecoilVsAvg_pi1115", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtAvgVsRE[2][1] = new TH2D("hPtRecoilVsAvg_ga1115", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtAvgVsRE[3][0] = new TH2D("hPtRecoilVsAvg_pi1520", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  hPtAvgVsRE[3][1] = new TH2D("hPtRecoilVsAvg_ga1520", "", nPtOA, pToaBin[0], pToaBin[1], nPtOA, pToaBin[0], pToaBin[1]);
  // errors
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      hDfJet[iBinEt][iBinId]      -> Sumw2();
      hDfRE[iBinEt][iBinId]       -> Sumw2();
      hPtRE[iBinEt][iBinId]       -> Sumw2();
      hPtUpVsDown[iBinEt][iBinId] -> Sumw2();
      hPtAvgVsRE[iBinEt][iBinId]  -> Sumw2();
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hDfUE[iBinEt][iBinId][iBinDf] -> Sumw2();
        hPtUE[iBinEt][iBinId][iBinDf] -> Sumw2();
      }
      for (UInt_t iCone = 0; iCone < NCones + 1; iCone++) {
        hPtOA[iBinEt][iBinId][iCone] -> Sumw2();
      }
    }  // end id loop
  }  // end eTtrg loop
  cout << "    Defined histograms." << endl;

  // define profiles
  TProfile *pPtUpVsDown[NBinsEt][NBinsId];
  TProfile *pPtAvgVsRE[NBinsEt][NBinsId];
  // off-axis profiles
  pPtUpVsDown[0][0] = new TProfile("pPtUpVsDown_pi920", "", nPt, pTbin, "S");
  pPtUpVsDown[0][1] = new TProfile("pPtUpVsDown_ga920", "", nPt, pTbin, "S");
  pPtUpVsDown[1][0] = new TProfile("pPtUpVsDown_pi911", "", nPt, pTbin, "S");
  pPtUpVsDown[1][1] = new TProfile("pPtUpVsDown_ga911", "", nPt, pTbin, "S");
  pPtUpVsDown[2][0] = new TProfile("pPtUpVsDown_pi1115", "", nPt, pTbin, "S");
  pPtUpVsDown[2][1] = new TProfile("pPtUpVsDown_ga1115", "", nPt, pTbin, "S");
  pPtUpVsDown[3][0] = new TProfile("pPtUpVsDown_pi1520", "", nPt, pTbin, "S");
  pPtUpVsDown[3][1] = new TProfile("pPtUpVsDown_ga1520", "", nPt, pTbin, "S");
  // off-axis vs. recoil profiles
  pPtAvgVsRE[0][0] = new TProfile("pPtAvgVsRE_pi920", "", nPt, pTbin, "S");
  pPtAvgVsRE[0][1] = new TProfile("pPtAvgVsRE_ga920", "", nPt, pTbin, "S");
  pPtAvgVsRE[1][0] = new TProfile("pPtAvgVsRE_pi911", "", nPt, pTbin, "S");
  pPtAvgVsRE[1][1] = new TProfile("pPtAvgVsRE_ga911", "", nPt, pTbin, "S");
  pPtAvgVsRE[2][0] = new TProfile("pPtAvgVsRE_pi1115", "", nPt, pTbin, "S");
  pPtAvgVsRE[2][1] = new TProfile("pPtAvgVsRE_ga1115", "", nPt, pTbin, "S");
  pPtAvgVsRE[3][0] = new TProfile("pPtAvgVsRE_pi1520", "", nPt, pTbin, "S");
  pPtAvgVsRE[3][1] = new TProfile("pPtAvgVsRE_ga1520", "", nPt, pTbin, "S");
  cout << "    Defined profiles." << endl;


  // for readability
  const Double_t piOver2(TMath::PiOver2());
  const Double_t piOver4(TMath::PiOver4());
  const Double_t piOver8(piOver4 / 2.);

  // define eTtrg and dF bins
  const Double_t eTtrgMin[NBinsEt]  = {9., 9., 11., 15.};
  const Double_t eTtrgMax[NBinsEt]  = {20., 11., 15., 20.};
  const Double_t dFdownMin[NBinsDf] = {piOver4, (3. * piOver8), piOver2};
  const Double_t dFdownMax[NBinsDf] = {piOver2, (5. * piOver8), (3. * piOver4)};
  const Double_t dFupMin[NBinsDf]   = {(3. * piOver2), (11. * piOver8), (5. * piOver4)};
  const Double_t dFupMax[NBinsDf]   = {(7. * piOver4), (13. * piOver8), (3. * piOver2)};
  const TString  sDfBins[NBinsDf]   = {"UE[1]", "UE[2]", "UE[3]"};
  cout << "    Defined various bins." << endl;


  // no. of triggers and jets
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


      // delta-phi plot
      if (recoilJetIsPresent) {
        hDfJet[0][idBin]     -> Fill(dFjet);
        hDfJet[eTbin][idBin] -> Fill(dFjet);
      }

      // determine dF region
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        const Bool_t isInLowerCone = ((dFjet > dFdownMin[iBinDf]) && (dFjet < dFdownMax[iBinDf]));
        const Bool_t isInUpperCone = ((dFjet > dFupMin[iBinDf]) && (dFjet < dFupMax[iBinDf]));
        const Bool_t isInUeRegion  = (isInLowerCone || isInUpperCone);
        if (isInUeRegion && recoilJetIsPresent) {
          hDfUE[0][idBin][iBinDf]     -> Fill(dFjet);
          hDfUE[eTbin][idBin][iBinDf] -> Fill(dFjet);
          hPtUE[0][idBin][iBinDf]     -> Fill(pTreco);
          hPtUE[eTbin][idBin][iBinDf] -> Fill(pTreco);
        }
      }

      const Bool_t isRecoil = (TMath::Abs(dFjet - TMath::Pi()) < dFrecoil);
      if (isRecoil) {
        hDfRE[0][idBin]           -> Fill(dFjet);
        hDfRE[eTbin][idBin]       -> Fill(dFjet);
        hPtRE[0][idBin]           -> Fill(pTreco);
        hPtRE[eTbin][idBin]       -> Fill(pTreco);
        hPtOA[0][idBin][0]        -> Fill(pTup);
        hPtOA[eTbin][idBin][0]    -> Fill(pTup);
        hPtOA[0][idBin][1]        -> Fill(pTdown);
        hPtOA[eTbin][idBin][1]    -> Fill(pTdown);
        hPtOA[0][idBin][2]        -> Fill(pTavg);
        hPtOA[eTbin][idBin][2]    -> Fill(pTavg);
        hPtUpVsDown[0][idBin]     -> Fill(pTdown, pTup);
        hPtUpVsDown[eTbin][idBin] -> Fill(pTdown, pTup);
        pPtUpVsDown[0][idBin]     -> Fill(pTdown, pTup);
        pPtUpVsDown[eTbin][idBin] -> Fill(pTdown, pTup);
        hPtAvgVsRE[0][idBin]      -> Fill(pTreco, pTavg);
        hPtAvgVsRE[eTbin][idBin]  -> Fill(pTreco, pTavg);
        pPtAvgVsRE[0][idBin]      -> Fill(pTreco, pTavg);
        pPtAvgVsRE[eTbin][idBin]  -> Fill(pTreco, pTavg);
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

      // delta-phi histograms
      const Double_t dFwidth = hDfJet[iBinEt][iBinId] -> GetBinWidth(17);
      const Double_t dFnorm  = (hBin * dFwidth) * nTrgBin[iBinEt][iBinId];
      hDfJet[iBinEt][iBinId] -> Scale(1. / dFnorm);
      hDfRE[iBinEt][iBinId]  -> Scale(1. / dFnorm);
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hDfUE[iBinEt][iBinId][iBinDf] -> Scale(1. / dFnorm);
      }

      // pT histograms
      const Double_t pTnorm = hBin * nTrgBin[iBinEt][iBinId];
      hPtRE[iBinEt][iBinId] -> Scale(1. / pTnorm);
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hPtUE[iBinEt][iBinId][iBinDf] -> Scale(1. / pTnorm);
      }
      for (UInt_t iCone = 0; iCone < NCones + 1; iCone++) {
        hPtOA[iBinEt][iBinId][iCone] -> Scale(1. / pTnorm);
      }

    }  // end id loop
  }  // end eTtrg loop
  cout << "    Normalized histograms." << endl;


  // calculate ratios
  const Double_t weight(1.);
  const TString  baseUE("hRatioUE_");
  const TString  baseOA("hRatioOA_");
  const TString  eTname[NBinsEt] = {"920", "911", "1115", "1520"};
  const TString  idName[NBinsId] = {"pi", "ga"};
  const TString  dFname[NBinsDf] = {"df0", "df1", "df2"};

  TH1D *hRatioOA[NBinsEt][NBinsId];
  TH1D *hRatioUE[NBinsEt][NBinsId][NBinsDf];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        // make UE name
        TString sRatioUE(baseUE.Data());
        sRatioUE += idName[iBinId].Data();
        sRatioUE += eTname[iBinEt].Data();
        sRatioUE += dFname[iBinDf].Data();

        // make UE histogram
        hRatioUE[iBinEt][iBinId][iBinDf] = new TH1D(sRatioUE.Data(), "", nPt, pTbin);
        hRatioUE[iBinEt][iBinId][iBinDf] -> Sumw2();
        hRatioUE[iBinEt][iBinId][iBinDf] -> Divide(hPtRE[iBinEt][iBinId], hPtUE[iBinEt][iBinId][iBinDf], weight, weight);
      }
      // make OA name
      TString sRatioOA(baseOA.Data());
      sRatioOA += idName[iBinId].Data();
      sRatioOA += eTname[iBinEt].Data();

      // make OA histogram
      hRatioOA[iBinEt][iBinId] = new TH1D(sRatioOA.Data(), "", nPt, pTbin);
      hRatioOA[iBinEt][iBinId] -> Sumw2();
      hRatioOA[iBinEt][iBinId] -> Divide(hPtRE[iBinEt][iBinId], hPtOA[iBinEt][iBinId][2], weight, weight);
    }  // end id loop
  }  // end eTtrg loop
  cout << "    Calculated ratios." << endl;


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColP(1);
  const UInt_t  fColO(878);
  const UInt_t  fMarJ(20);
  const UInt_t  fMarD(1);
  const UInt_t  fMarP(20);
  const UInt_t  fMarO(28);
  const UInt_t  fLinD(2);
  const UInt_t  fFilR(3354);
  const UInt_t  fColR[NBinsId] = {858, 898};
  const UInt_t  fMarR[NBinsId] = {29, 20};
  const UInt_t  fColU[NBinsDf] = {798, 818, 838};
  const UInt_t  fMarU[NBinsDf] = {24, 25, 26};
  const UInt_t  fFilU[NBinsDf] = {3345, 3354, 3345};
  const Float_t fLab(0.03);
  const Float_t fLabR(0.055);
  const Float_t fSizR(0.074);
  const Float_t fOffsetX(1.);
  const Float_t fOffsetXR(0.77);
  const Float_t fOffsetY(1.07);
  const Float_t fOffsetYR(0.57);
  const TString sTitleXdf("#Delta#varphi = #varphi^{jet} - #varphi^{trg}");
  const TString sTitleYdf("(1/N^{trg}) dN^{jet}/(d#Delta#varphi d#eta)");
  const TString sTitleXpt("p_{T}^{reco} = p_{T}^{jet} - #rho #upoint A_{jet}");
  const TString sTitleYpt("(1/N^{trg}) dN^{jet})/(dp_{T}^{jet} d#eta)");
  const TString sTitleYue("RE/UE[n], RE/OA");
  const TString sTitleXud("dp_{T}^{-} = #varsigma_{T}^{-} #upoint A_{jet}");
  const TString sTitleYud("dp_{T}^{+} = #varsigma_{T}^{+} #upoint A_{jet}");
  const TString sTitleXre("p_{T}^{reco} = p_{T}^{jet} - #rho #upoint A_{jet}");
  const TString sTitleYre("dp_{T} = #frac{1}{2} (#varsigma_{T}^{-} + #varsigma_{T}^{+}) #upoint A_{jet}");
  const TString sTitleDf[NBinsId] = {"#pi^{0} trigger", "#gamma^{rich} trigger"};
  const TString sTitlePt[NBinsId] = {"#pi^{0} trigger", "#gamma^{rich} trigger"};
  const TString sTitleUd[NBinsId] = {"#pi^{0} trigger", "#gamma^{rich} trigger"};
  const TString sTitleRe[NBinsId] = {"#pi^{0} trigger", "#gamma^{rich} trigger"};
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      // delta-phi histograms
      hDfJet[iBinEt][iBinId] -> SetMarkerStyle(fMarJ);
      hDfJet[iBinEt][iBinId] -> SetTitle(sTitleDf[iBinId].Data());
      hDfJet[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hDfJet[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleXdf.Data());
      hDfJet[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hDfJet[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hDfJet[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hDfJet[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hDfJet[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleYdf.Data());
      hDfJet[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hDfJet[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hDfJet[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hDfJet[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      hDfJet[iBinEt][iBinId] -> GetZaxis() -> SetLabelSize(fLab);
      hDfRE[iBinEt][iBinId]  -> SetFillStyle(fFilR);
      hDfRE[iBinEt][iBinId]  -> SetFillColor(fColR[iBinId]);
      hDfRE[iBinEt][iBinId]  -> SetLineStyle(fLinD);
      hDfRE[iBinEt][iBinId]  -> SetLineColor(fColR[iBinId]);
      hDfRE[iBinEt][iBinId]  -> SetMarkerStyle(fMarD);
      hDfRE[iBinEt][iBinId]  -> SetMarkerColor(fColR[iBinId]);
      hDfRE[iBinEt][iBinId]  -> SetTitle(sTitleDf[iBinId].Data());
      hDfRE[iBinEt][iBinId]  -> SetTitleFont(fTxt);
      hDfRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitle(sTitleXdf.Data());
      hDfRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitleFont(fTxt);
      hDfRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hDfRE[iBinEt][iBinId]  -> GetXaxis() -> CenterTitle(fCnt);
      hDfRE[iBinEt][iBinId]  -> GetXaxis() -> SetLabelSize(fLab);
      hDfRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitle(sTitleYdf.Data());
      hDfRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitleFont(fTxt);
      hDfRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hDfRE[iBinEt][iBinId]  -> GetYaxis() -> CenterTitle(fCnt);
      hDfRE[iBinEt][iBinId]  -> GetYaxis() -> SetLabelSize(fLab);
      hDfRE[iBinEt][iBinId]  -> GetZaxis() -> SetLabelSize(fLab);
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hDfUE[iBinEt][iBinId][iBinDf] -> SetFillStyle(fFilU[iBinDf]);
        hDfUE[iBinEt][iBinId][iBinDf] -> SetFillColor(fColU[iBinDf]);
        hDfUE[iBinEt][iBinId][iBinDf] -> SetLineStyle(fLinD);
        hDfUE[iBinEt][iBinId][iBinDf] -> SetLineColor(fColU[iBinDf]);
        hDfUE[iBinEt][iBinId][iBinDf] -> SetMarkerStyle(fMarD);
        hDfUE[iBinEt][iBinId][iBinDf] -> SetMarkerColor(fColU[iBinDf]);
        hDfUE[iBinEt][iBinId][iBinDf] -> SetTitle(sTitleDf[iBinId].Data());
        hDfUE[iBinEt][iBinId][iBinDf] -> SetTitleFont(fTxt);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitle(sTitleXdf.Data());
        hDfUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleFont(fTxt);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> CenterTitle(fCnt);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetLabelSize(fLab);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitle(sTitleYdf.Data());
        hDfUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleFont(fTxt);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleOffset(fOffsetY);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> CenterTitle(fCnt);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetLabelSize(fLab);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetZaxis() -> SetLabelSize(fLab);
      }

      // pT histograms
      hPtRE[iBinEt][iBinId] -> SetLineColor(fColR[iBinId]);
      hPtRE[iBinEt][iBinId] -> SetMarkerStyle(fMarR[iBinId]);
      hPtRE[iBinEt][iBinId] -> SetMarkerColor(fColR[iBinId]);
      hPtRE[iBinEt][iBinId] -> SetTitle(sTitlePt[iBinId].Data());
      hPtRE[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleXpt.Data());
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hPtRE[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleYpt.Data());
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hPtRE[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hPtUE[iBinEt][iBinId][iBinDf]    -> SetLineColor(fColU[iBinDf]);
        hPtUE[iBinEt][iBinId][iBinDf]    -> SetMarkerStyle(fMarU[iBinDf]);
        hPtUE[iBinEt][iBinId][iBinDf]    -> SetMarkerColor(fColU[iBinDf]);
        hPtUE[iBinEt][iBinId][iBinDf]    -> SetTitle(sTitlePt[iBinId].Data());
        hPtUE[iBinEt][iBinId][iBinDf]    -> SetTitleFont(fTxt);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetXaxis() -> SetTitle(sTitleXpt.Data());
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetXaxis() -> SetTitleFont(fTxt);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetXaxis() -> CenterTitle(fCnt);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetXaxis() -> SetLabelSize(fLab);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetYaxis() -> SetTitle(sTitleYpt.Data());
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetYaxis() -> SetTitleFont(fTxt);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetYaxis() -> SetTitleOffset(fOffsetY);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetYaxis() -> CenterTitle(fCnt);
        hPtUE[iBinEt][iBinId][iBinDf]    -> GetYaxis() -> SetLabelSize(fLab);
        hRatioUE[iBinEt][iBinId][iBinDf] -> SetLineColor(fColU[iBinDf]);
        hRatioUE[iBinEt][iBinId][iBinDf] -> SetMarkerStyle(fMarU[iBinDf]);
        hRatioUE[iBinEt][iBinId][iBinDf] -> SetMarkerColor(fColU[iBinDf]);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitle(sTitleXpt.Data());
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleFont(fTxt);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleSize(fSizR);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleOffset(fOffsetXR);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> CenterTitle(fCnt);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetLabelSize(fLabR);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitle(sTitleYue.Data());
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleFont(fTxt);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleSize(fSizR);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleOffset(fOffsetYR);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> CenterTitle(fCnt);
        hRatioUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetLabelSize(fLabR);
      }
      for (UInt_t iCone = 0; iCone < NCones + 1; iCone++) {
        hPtOA[iBinEt][iBinId][iCone]    -> SetLineColor(fColO);
        hPtOA[iBinEt][iBinId][iCone]    -> SetMarkerStyle(fMarO);
        hPtOA[iBinEt][iBinId][iCone]    -> SetMarkerColor(fColO);
        hPtOA[iBinEt][iBinId][iCone]    -> SetTitle(sTitlePt[iBinId].Data());
        hPtOA[iBinEt][iBinId][iCone]    -> SetTitleFont(fTxt);
        hPtOA[iBinEt][iBinId][iCone]    -> GetXaxis() -> SetTitle(sTitleXpt.Data());
        hPtOA[iBinEt][iBinId][iCone]    -> GetXaxis() -> SetTitleFont(fTxt);
        hPtOA[iBinEt][iBinId][iCone]    -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hPtOA[iBinEt][iBinId][iCone]    -> GetXaxis() -> CenterTitle(fCnt);
        hPtOA[iBinEt][iBinId][iCone]    -> GetXaxis() -> SetLabelSize(fLab);
        hPtOA[iBinEt][iBinId][iCone]    -> GetYaxis() -> SetTitle(sTitleYpt.Data());
        hPtOA[iBinEt][iBinId][iCone]    -> GetYaxis() -> SetTitleFont(fTxt);
        hPtOA[iBinEt][iBinId][iCone]    -> GetYaxis() -> SetTitleOffset(fOffsetY);
        hPtOA[iBinEt][iBinId][iCone]    -> GetYaxis() -> CenterTitle(fCnt);
        hPtOA[iBinEt][iBinId][iCone]    -> GetYaxis() -> SetLabelSize(fLab);
      }
      hRatioOA[iBinEt][iBinId] -> SetLineColor(fColO);
      hRatioOA[iBinEt][iBinId] -> SetMarkerStyle(fMarO);
      hRatioOA[iBinEt][iBinId] -> SetMarkerColor(fColO);
      hRatioOA[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleXpt.Data());
      hRatioOA[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hRatioOA[iBinEt][iBinId] -> GetXaxis() -> SetTitleSize(fSizR);
      hRatioOA[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetXR);
      hRatioOA[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hRatioOA[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLabR);
      hRatioOA[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleYue.Data());
      hRatioOA[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hRatioOA[iBinEt][iBinId] -> GetYaxis() -> SetTitleSize(fSizR);
      hRatioOA[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetYR);
      hRatioOA[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hRatioOA[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLabR);

      // 2d pT histograms
      hPtUpVsDown[iBinEt][iBinId] -> SetTitle(sTitleUd[iBinId].Data());
      hPtUpVsDown[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleXud.Data());
      hPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleYud.Data());
      hPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      hPtUpVsDown[iBinEt][iBinId] -> GetZaxis() -> SetLabelSize(fLab);
      hPtAvgVsRE[iBinEt][iBinId]  -> SetTitle(sTitleRe[iBinId].Data());
      hPtAvgVsRE[iBinEt][iBinId]  -> SetTitleFont(fTxt);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitle(sTitleXre.Data());
      hPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitleFont(fTxt);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> CenterTitle(fCnt);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetLabelSize(fLab);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitle(sTitleYre.Data());
      hPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitleFont(fTxt);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> CenterTitle(fCnt);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetLabelSize(fLab);
      hPtAvgVsRE[iBinEt][iBinId]  -> GetZaxis() -> SetLabelSize(fLab);

      // pT profiles
      pPtUpVsDown[iBinEt][iBinId] -> SetMarkerStyle(fMarP);
      pPtUpVsDown[iBinEt][iBinId] -> SetMarkerColor(fColP);
      pPtUpVsDown[iBinEt][iBinId] -> SetTitle(sTitleUd[iBinId].Data());
      pPtUpVsDown[iBinEt][iBinId] -> SetTitleFont(fTxt);
      pPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sTitleXud.Data());
      pPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      pPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      pPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      pPtUpVsDown[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      pPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleYud.Data());
      pPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      pPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      pPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      pPtUpVsDown[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      pPtAvgVsRE[iBinEt][iBinId]  -> SetTitle(sTitleRe[iBinId].Data());
      pPtAvgVsRE[iBinEt][iBinId]  -> SetTitleFont(fTxt);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitle(sTitleXre.Data());
      pPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitleFont(fTxt);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetTitleOffset(fOffsetX);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> CenterTitle(fCnt);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetXaxis() -> SetLabelSize(fLab);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitle(sTitleYre.Data());
      pPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitleFont(fTxt);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetTitleOffset(fOffsetY);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> CenterTitle(fCnt);
      pPtAvgVsRE[iBinEt][iBinId]  -> GetYaxis() -> SetLabelSize(fLab);

    }  // end id loop
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
  const TString sOaName("OA cone");
  const TString sIdName[NBinsId] = {"#pi^{0}-jets", "#gamma^{rich}-jets"};

  TPaveText *pPtJet[NBinsEt];
  TLegend   *lPtJet[NBinsId];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    TString sEtTrg(sTrigger.Data());
    TString sJetStuff(sJet.Data());
    sEtTrg    += eTtrgMin[iBinEt];
    sEtTrg    += ", ";
    sEtTrg    += eTtrgMax[iBinEt];
    sEtTrg    += ") GeV";
    sJetStuff += rJet;

    pPtJet[iBinEt] = new TPaveText(xyPav1, xyPav1, xyPav2, xyPav2, "NDC NB");
    pPtJet[iBinEt] -> SetFillColor(fColT);
    pPtJet[iBinEt] -> SetFillStyle(fFilT);
    pPtJet[iBinEt] -> SetLineColor(fColT);
    pPtJet[iBinEt] -> SetLineStyle(fLinT);
    pPtJet[iBinEt] -> SetTextFont(fTxt);
    pPtJet[iBinEt] -> SetTextAlign(fAlign);
    pPtJet[iBinEt] -> AddText(sSystem.Data());
    pPtJet[iBinEt] -> AddText(sEtTrg.Data());
    pPtJet[iBinEt] -> AddText(sJetStuff.Data());
    if (JetType == 0)
      pPtJet[iBinEt] -> AddText(sChrg.Data());
    else
      pPtJet[iBinEt] -> AddText(sFull.Data());
  }
  for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
    lPtJet[iBinId] = new TLegend(xyLeg1, xyLeg1, xyLeg2, xyLeg2);
    lPtJet[iBinId] -> SetFillColor(fColT);
    lPtJet[iBinId] -> SetFillStyle(fFilT);
    lPtJet[iBinId] -> SetLineColor(fColT);
    lPtJet[iBinId] -> SetLineStyle(fLinT);
    lPtJet[iBinId] -> SetTextFont(fTxt);
    lPtJet[iBinId] -> AddEntry(hPtRE[0][iBinId], sIdName[iBinId].Data());
    lPtJet[iBinId] -> AddEntry(hPtUE[0][iBinId][0], sDfBins[0].Data());
    lPtJet[iBinId] -> AddEntry(hPtUE[0][iBinId][1], sDfBins[1].Data());
    lPtJet[iBinId] -> AddEntry(hPtUE[0][iBinId][2], sDfBins[2].Data());
    lPtJet[iBinId] -> AddEntry(hPtOA[0][iBinId][2], sOaName.Data());
  }
  cout << "    Made labels." << endl;


  // make line
  const UInt_t fColL(1);
  const UInt_t fLinL(2);

  TLine *lOnePt = new TLine(pTbin[0], 1., pTbin[NBinsPt - 1], 1.);
  lOnePt -> SetLineColor(fColL);
  lOnePt -> SetLineStyle(fLinL);
  cout << "    Made line." << endl;


  // create directories and save histograms
  const TString sBinsEt[NBinsEt] = {"eTtrg920", "eTtrg911", "eTtrg1115", "eTtrg1520"};
  const TString sBinsId[NBinsId]   = {"pi0", "gamma"};

  TDirectory *dBinsEt[NBinsEt];
  TDirectory *dBinsId[NBinsEt][NBinsId];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    dBinsEt[iBinEt] = (TDirectory*) fOutput -> mkdir(sBinsEt[iBinEt].Data());
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      dBinsId[iBinEt][iBinId] = (TDirectory*) dBinsEt[iBinEt] -> mkdir(sBinsId[iBinId].Data());
      dBinsId[iBinEt][iBinId] -> cd();
      hDfJet[iBinEt][iBinId]  -> Write();
      hDfRE[iBinEt][iBinId]   -> Write();
      hPtRE[iBinEt][iBinId]   -> Write();
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hDfUE[iBinEt][iBinId][iBinDf]    -> Write();
        hPtUE[iBinEt][iBinId][iBinDf]    -> Write();
        hRatioUE[iBinEt][iBinId][iBinDf] -> Write();
      }
      for (UInt_t iCone = 0; iCone < NCones + 1; iCone++) {
        hPtOA[iBinEt][iBinId][iCone] -> Write();
      }
      hRatioOA[iBinEt][iBinId]    -> Write();
      hPtUpVsDown[iBinEt][iBinId] -> Write();
      hPtAvgVsRE[iBinEt][iBinId]  -> Write();
      pPtUpVsDown[iBinEt][iBinId] -> Write();
      pPtAvgVsRE[iBinEt][iBinId]  -> Write();
    }  // end id loop
  }  // end eTtrg loop
  cout << "    Made directories." << endl;


  // draw plots
  const UInt_t  widthDf(1500);
  const UInt_t  heightDf(750);
  const UInt_t  widthPt(750);
  const UInt_t  heightPt(900);
  const UInt_t  widthUd(1500);
  const UInt_t  heightUd(750);
  const UInt_t  frame(0);
  const UInt_t  grid(0);
  const UInt_t  log(1);
  const Float_t noMargin(0.);
  const Float_t marginR(0.07);
  const Float_t marginRcolz(0.13);
  const Float_t marginL(0.13);
  const Float_t marginT(0.07);
  const Float_t marginB(0.13);
  const Float_t marginBR(0.17);
  const Float_t xPadDf[NBinsId + 2]  = {0., 0.5, 0.5, 1.};
  const Float_t yPadDf[NBinsId + 2]  = {0., 1., 0., 1.};
  const Float_t xPadPt[NPadsPt + 2]  = {0., 1., 0., 1.};
  const Float_t yPadPt[NPadsPt + 2]  = {0., 0.35, 0.35, 1.};
  const Float_t xPadUd[NBinsId + 2]  = {0., 0.5, 0.5, 1.};
  const Float_t yPadUd[NBinsId + 2]  = {0., 1., 0., 1.};
  const Float_t xPadRe[NBinsId + 2]  = {0., 0.5, 0.5, 1.};
  const Float_t yPadRe[NBinsId + 2]  = {0., 1., 0., 1.};
  const TString sDfCanvas[NBinsEt]   = {"cDfJet_et920", "cDfJet_et911", "cDfJet_et1115", "cDfJet_et1520"};
  const TString sPtCanvasP[NBinsEt]  = {"cPtJet_pi920", "cPtJet_pi911", "cPtJet_pi1115", "cPtJet_pi1520"};
  const TString sPtCanvasG[NBinsEt]  = {"cPtJet_ga920", "cPtJet_ga911", "cPtJet_ga1115", "cPtJet_ga1520"};
  const TString sUdCanvas[NBinsEt]   = {"cPtUpVsDown_et920", "cPtUpVsDown_et911", "cPtUpvsDown_et1115", "cPtUpVsDown_et1520"};
  const TString sReCanvas[NBinsEt]   = {"cPtAvgVsRe_et920", "cPtAvgVsRe_et911", "cPtAvgVsRe_et1115", "cPtAvgVsRe_et1520"};
  const TString sDfPads[NBinsId]     = {"pPi0", "pGamma"};
  const TString sPtPads[NPadsPt]     = {"pRatio", "pJets"};
  const TString sUdPads[NBinsId]     = {"pPi0", "pGamma"};
  const TString sRePads[NBinsId]     = {"pPi0", "pGamma"};

  TPad    *pDfJet[NBinsEt][NBinsId];
  TPad    *pPtPi0[NBinsEt][NPadsPt];
  TPad    *pPtGam[NBinsEt][NPadsPt];
  TPad    *pUdJet[NBinsEt][NBinsId];
  TPad    *pReJet[NBinsEt][NBinsId];
  TCanvas *cDfJet[NBinsEt];
  TCanvas *cPtPi0[NBinsEt];
  TCanvas *cPtGam[NBinsEt];
  TCanvas *cUdJet[NBinsEt];
  TCanvas *cReJet[NBinsEt];
  fOutput -> cd();
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {

    // delta-phi plots
    dBinsEt[iBinEt] -> cd();
    cDfJet[iBinEt]    = new TCanvas(sDfCanvas[iBinEt].Data(), "", widthDf, heightDf);
    pDfJet[iBinEt][0] = new TPad(sDfPads[0].Data(), "", xPadDf[0], yPadDf[0], xPadDf[1], yPadDf[1]);
    pDfJet[iBinEt][1] = new TPad(sDfPads[1].Data(), "", xPadDf[2], yPadDf[2], xPadDf[3], yPadDf[3]);
    pDfJet[iBinEt][0] -> SetGrid(grid, grid);
    pDfJet[iBinEt][0] -> SetFrameBorderMode(frame);
    pDfJet[iBinEt][0] -> SetRightMargin(noMargin);
    pDfJet[iBinEt][0] -> SetLeftMargin(marginL);
    pDfJet[iBinEt][0] -> SetTopMargin(marginT);
    pDfJet[iBinEt][0] -> SetBottomMargin(marginB);
    pDfJet[iBinEt][1] -> SetGrid(grid, grid);
    pDfJet[iBinEt][1] -> SetFrameBorderMode(frame);
    pDfJet[iBinEt][1] -> SetRightMargin(marginR);
    pDfJet[iBinEt][1] -> SetLeftMargin(noMargin);
    pDfJet[iBinEt][1] -> SetTopMargin(marginT);
    pDfJet[iBinEt][1] -> SetBottomMargin(marginB);
    cDfJet[iBinEt]    -> cd();
    pDfJet[iBinEt][0] -> Draw();
    pDfJet[iBinEt][1] -> Draw();
    pDfJet[iBinEt][0] -> cd();
    hDfJet[iBinEt][0] -> Draw();
    hDfRE[iBinEt][0]  -> Draw("SAME LF HIST");
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hDfUE[iBinEt][0][iBinDf] -> Draw("SAME LF HIST");
    }
    pDfJet[iBinEt][1] -> cd();
    hDfJet[iBinEt][1] -> Draw();
    hDfRE[iBinEt][1]  -> Draw("SAME LF HIST");
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hDfUE[iBinEt][1][iBinDf] -> Draw("SAME LF HIST");
    }
    pPtJet[iBinEt] -> Draw();
    cDfJet[iBinEt] -> Write();
    cDfJet[iBinEt] -> Close();

    // pT plots
    dBinsEt[iBinEt] -> cd();
    cPtPi0[iBinEt]    = new TCanvas(sPtCanvasP[iBinEt].Data(), "", widthPt, heightPt);
    pPtPi0[iBinEt][0] = new TPad(sPtPads[0].Data(), "", xPadPt[0], yPadPt[0], xPadPt[1], yPadPt[1]);
    pPtPi0[iBinEt][1] = new TPad(sPtPads[1].Data(), "", xPadPt[2], yPadPt[2], xPadPt[3], yPadPt[3]);
    pPtPi0[iBinEt][0] -> SetGrid(grid, grid);
    pPtPi0[iBinEt][0] -> SetLogy(log);
    pPtPi0[iBinEt][0] -> SetFrameBorderMode(frame);
    pPtPi0[iBinEt][0] -> SetRightMargin(marginR);
    pPtPi0[iBinEt][0] -> SetLeftMargin(marginL);
    pPtPi0[iBinEt][0] -> SetTopMargin(noMargin);
    pPtPi0[iBinEt][0] -> SetBottomMargin(marginBR);
    pPtPi0[iBinEt][1] -> SetGrid(grid, grid);
    pPtPi0[iBinEt][1] -> SetLogy(log);
    pPtPi0[iBinEt][1] -> SetFrameBorderMode(frame);
    pPtPi0[iBinEt][1] -> SetRightMargin(marginR);
    pPtPi0[iBinEt][1] -> SetLeftMargin(marginL);
    pPtPi0[iBinEt][1] -> SetTopMargin(marginT);
    pPtPi0[iBinEt][1] -> SetBottomMargin(noMargin);
    cPtPi0[iBinEt]    -> cd();
    pPtPi0[iBinEt][0] -> Draw();
    pPtPi0[iBinEt][1] -> Draw();
    pPtPi0[iBinEt][0] -> cd();
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      if (iBinDf == 0)
        hRatioUE[iBinEt][0][iBinDf] -> Draw();
      else
        hRatioUE[iBinEt][0][iBinDf] -> Draw("SAME");
    }
    hRatioOA[iBinEt][0] -> Draw("SAME");
    lOnePt              -> Draw();
    pPtPi0[iBinEt][1]   -> cd();
    hPtRE[iBinEt][0]    -> Draw();
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hPtUE[iBinEt][0][iBinDf] -> Draw("SAME");
    }
    hPtOA[iBinEt][0][2] -> Draw("SAME");
    lPtJet[0]           -> Draw();
    pPtJet[iBinEt]      -> Draw();
    cPtPi0[iBinEt]      -> Write();
    cPtPi0[iBinEt]      -> Close();

    dBinsEt[iBinEt] -> cd();
    cPtGam[iBinEt]    = new TCanvas(sPtCanvasG[iBinEt].Data(), "", widthPt, heightPt);
    pPtGam[iBinEt][0] = new TPad(sPtPads[0].Data(), "", xPadPt[0], yPadPt[0], xPadPt[1], yPadPt[1]);
    pPtGam[iBinEt][1] = new TPad(sPtPads[1].Data(), "", xPadPt[2], yPadPt[2], xPadPt[3], yPadPt[3]);
    pPtGam[iBinEt][0] -> SetGrid(grid, grid);
    pPtGam[iBinEt][0] -> SetLogy(log);
    pPtGam[iBinEt][0] -> SetFrameBorderMode(frame);
    pPtGam[iBinEt][0] -> SetRightMargin(marginR);
    pPtGam[iBinEt][0] -> SetLeftMargin(marginL);
    pPtGam[iBinEt][0] -> SetTopMargin(noMargin);
    pPtGam[iBinEt][0] -> SetBottomMargin(marginBR);
    pPtGam[iBinEt][1] -> SetGrid(grid, grid);
    pPtGam[iBinEt][1] -> SetLogy(log);
    pPtGam[iBinEt][1] -> SetFrameBorderMode(frame);
    pPtGam[iBinEt][1] -> SetRightMargin(marginR);
    pPtGam[iBinEt][1] -> SetLeftMargin(marginL);
    pPtGam[iBinEt][1] -> SetTopMargin(marginT);
    pPtGam[iBinEt][1] -> SetBottomMargin(noMargin);
    cPtGam[iBinEt]    -> cd();
    pPtGam[iBinEt][0] -> Draw();
    pPtGam[iBinEt][1] -> Draw();
    pPtGam[iBinEt][0] -> cd();
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      if (iBinDf == 0)
        hRatioUE[iBinEt][1][iBinDf] -> Draw();
      else
        hRatioUE[iBinEt][1][iBinDf] -> Draw("SAME");
    }
    hRatioOA[iBinEt][1] -> Draw("SAME");
    lOnePt              -> Draw();
    pPtGam[iBinEt][1]   -> cd();
    hPtRE[iBinEt][1]    -> Draw();
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hPtUE[iBinEt][1][iBinDf] -> Draw("SAME");
    }
    hPtOA[iBinEt][1][2] -> Draw("SAME");
    lPtJet[1]           -> Draw();
    pPtJet[iBinEt]      -> Draw();
    cPtGam[iBinEt]      -> Write();
    cPtGam[iBinEt]      -> Close();

    // pTup vs. pTdown
    dBinsEt[iBinEt] -> cd();
    cUdJet[iBinEt]    = new TCanvas(sUdCanvas[iBinEt].Data(), "", widthDf, heightDf);
    pUdJet[iBinEt][0] = new TPad(sUdPads[0].Data(), "", xPadUd[0], yPadUd[0], xPadUd[1], yPadUd[1]);
    pUdJet[iBinEt][1] = new TPad(sUdPads[1].Data(), "", xPadUd[2], yPadUd[2], xPadUd[3], yPadUd[3]);
    pUdJet[iBinEt][0]      -> SetGrid(grid, grid);
    pUdJet[iBinEt][0]      -> SetLogz(log);
    pUdJet[iBinEt][0]      -> SetFrameBorderMode(frame);
    pUdJet[iBinEt][0]      -> SetRightMargin(marginRcolz);
    pUdJet[iBinEt][0]      -> SetLeftMargin(marginL);
    pUdJet[iBinEt][0]      -> SetTopMargin(marginT);
    pUdJet[iBinEt][0]      -> SetBottomMargin(marginB);
    pUdJet[iBinEt][1]      -> SetGrid(grid, grid);
    pUdJet[iBinEt][1]      -> SetLogz(log);
    pUdJet[iBinEt][1]      -> SetFrameBorderMode(frame);
    pUdJet[iBinEt][1]      -> SetRightMargin(marginRcolz);
    pUdJet[iBinEt][1]      -> SetLeftMargin(marginL);
    pUdJet[iBinEt][1]      -> SetTopMargin(marginT);
    pUdJet[iBinEt][1]      -> SetBottomMargin(marginB);
    cUdJet[iBinEt]         -> cd();
    pUdJet[iBinEt][0]      -> Draw();
    pUdJet[iBinEt][1]      -> Draw();
    pUdJet[iBinEt][0]      -> cd();
    hPtUpVsDown[iBinEt][0] -> Draw("COLZ");
    pPtUpVsDown[iBinEt][0] -> Draw("SAME");
    pUdJet[iBinEt][1]      -> cd();
    hPtUpVsDown[iBinEt][1] -> Draw("COLZ");
    pPtUpVsDown[iBinEt][1] -> Draw("SAME");
    pPtJet[iBinEt]         -> Draw();
    cUdJet[iBinEt]         -> Write();
    cUdJet[iBinEt]         -> Close();

    // recoil pT vs. pTavg
    dBinsEt[iBinEt] -> cd();
    cReJet[iBinEt]    = new TCanvas(sReCanvas[iBinEt].Data(), "", widthDf, heightDf);
    pReJet[iBinEt][0] = new TPad(sRePads[0].Data(), "", xPadRe[0], yPadRe[0], xPadRe[1], yPadRe[1]);
    pReJet[iBinEt][1] = new TPad(sRePads[1].Data(), "", xPadRe[2], yPadRe[2], xPadRe[3], yPadRe[3]);
    pReJet[iBinEt][0]      -> SetGrid(grid, grid);
    pReJet[iBinEt][0]      -> SetLogz(log);
    pReJet[iBinEt][0]      -> SetFrameBorderMode(frame);
    pReJet[iBinEt][0]      -> SetRightMargin(marginRcolz);
    pReJet[iBinEt][0]      -> SetLeftMargin(marginL);
    pReJet[iBinEt][0]      -> SetTopMargin(marginT);
    pReJet[iBinEt][0]      -> SetBottomMargin(marginB);
    pReJet[iBinEt][1]      -> SetGrid(grid, grid);
    pReJet[iBinEt][1]      -> SetLogz(log);
    pReJet[iBinEt][1]      -> SetFrameBorderMode(frame);
    pReJet[iBinEt][1]      -> SetRightMargin(marginRcolz);
    pReJet[iBinEt][1]      -> SetLeftMargin(marginL);
    pReJet[iBinEt][1]      -> SetTopMargin(marginT);
    pReJet[iBinEt][1]      -> SetBottomMargin(marginB);
    cReJet[iBinEt]         -> cd();
    pReJet[iBinEt][0]      -> Draw();
    pReJet[iBinEt][1]      -> Draw();
    pReJet[iBinEt][0]      -> cd();
    hPtAvgVsRE[iBinEt][0]  -> Draw("COLZ");
    pPtAvgVsRE[iBinEt][0]  -> Draw("SAME");
    pReJet[iBinEt][1]      -> cd();
    hPtAvgVsRE[iBinEt][1]  -> Draw("COLZ");
    pPtAvgVsRE[iBinEt][1]  -> Draw("SAME");
    pPtJet[iBinEt]         -> Draw();
    cReJet[iBinEt]         -> Write();
    cReJet[iBinEt]         -> Close();

  }  // end eTtrg loop
  cout << "    Drew plots." << endl;

  // save and close
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Study finished!\n" << endl;

}

// End ------------------------------------------------------------------------
