// 'CalculateAverNumber.C'
// Derek Anderson
// 01.23.2018
//
// This takes the output of 'StJetTreeMaker' and
// calculates the average number of recoil jets
// and the average number of "uncorrelated
// jets."
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
static const UInt_t NBinsEt = 4;
static const UInt_t NBinsId = 2;
static const UInt_t NBinsDf = 3;
static const UInt_t NPadsPt = 2;
static const UInt_t JetType = 1;



void CalculateAverageNumber(const Bool_t isInBatchMode=false) {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning average number calculation..." << endl;

  // io parameters
  const TString sOutput("pp200r9.avgReVsUeNums.eTtrg920.r03a02rm1full.d31m1y2018.root");
  const TString sInput("output/CollaborationMeetingJan2018/pp200r9.withOaCones.eTtrg920.r03rm1full.d21m1y2018.root");
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
  TH1D *hNumRE[NBinsEt][NBinsId];
  TH1D *hDfUE[NBinsEt][NBinsId][NBinsDf];
  TH1D *hNumUE[NBinsEt][NBinsId][NBinsDf];
  TH2D *hNumUeVsRe[NBinsEt][NBinsId][NBinsDf];

  const UInt_t   nDf(60);
  const UInt_t   nNum(100);
  const Double_t dFbin[2] = {0., TMath::TwoPi()};
  const Double_t num[2]   = {0., (Double_t) nNum};
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
  // recoil number plots
  hNumRE[0][0] = new TH1D("hNumRecoil_pi920", "", nNum, num[0], num[1]);
  hNumRE[0][1] = new TH1D("hNumRecoil_ga920", "", nNum, num[0], num[1]);
  hNumRE[1][0] = new TH1D("hNumRecoil_pi911", "", nNum, num[0], num[1]);
  hNumRE[1][1] = new TH1D("hNumRecoil_ga911", "", nNum, num[0], num[1]);
  hNumRE[2][0] = new TH1D("hNumRecoil_pi1115", "", nNum, num[0], num[1]);
  hNumRE[2][1] = new TH1D("hNumRecoil_ga1115", "", nNum, num[0], num[1]);
  hNumRE[3][0] = new TH1D("hNumRecoil_pi1520", "", nNum, num[0], num[1]);
  hNumRE[3][1] = new TH1D("hNumRecoil_ga1520", "", nNum, num[0], num[1]);
  // uncorrelated number plots
  hNumUE[0][0][0] = new TH1D("hNumUncorrelated_pi920df0", "", nNum, num[0], num[1]);
  hNumUE[0][0][1] = new TH1D("hNumUncorrelated_pi920df1", "", nNum, num[0], num[1]);
  hNumUE[0][0][2] = new TH1D("hNumUncorrelated_pi920df2", "", nNum, num[0], num[1]);
  hNumUE[0][1][0] = new TH1D("hNumUncorrelated_ga920df0", "", nNum, num[0], num[1]);
  hNumUE[0][1][1] = new TH1D("hNumUncorrelated_ga920df1", "", nNum, num[0], num[1]);
  hNumUE[0][1][2] = new TH1D("hNumUncorrelated_ga920df2", "", nNum, num[0], num[1]);
  hNumUE[1][0][0] = new TH1D("hNumUncorrelated_pi911df0", "", nNum, num[0], num[1]);
  hNumUE[1][0][1] = new TH1D("hNumUncorrelated_pi911df1", "", nNum, num[0], num[1]);
  hNumUE[1][0][2] = new TH1D("hNumUncorrelated_pi911df2", "", nNum, num[0], num[1]);
  hNumUE[1][1][0] = new TH1D("hNumUncorrelated_ga911df0", "", nNum, num[0], num[1]);
  hNumUE[1][1][1] = new TH1D("hNumUncorrelated_ga911df1", "", nNum, num[0], num[1]);
  hNumUE[1][1][2] = new TH1D("hNumUncorrelated_ga911df2", "", nNum, num[0], num[1]);
  hNumUE[2][0][0] = new TH1D("hNumUncorrelated_pi1115df0", "", nNum, num[0], num[1]);
  hNumUE[2][0][1] = new TH1D("hNumUncorrelated_pi1115df1", "", nNum, num[0], num[1]);
  hNumUE[2][0][2] = new TH1D("hNumUncorrelated_pi1115df2", "", nNum, num[0], num[1]);
  hNumUE[2][1][0] = new TH1D("hNumUncorrelated_ga1115df0", "", nNum, num[0], num[1]);
  hNumUE[2][1][1] = new TH1D("hNumUncorrelated_ga1115df1", "", nNum, num[0], num[1]);
  hNumUE[2][1][2] = new TH1D("hNumUncorrelated_ga1115df2", "", nNum, num[0], num[1]);
  hNumUE[3][0][0] = new TH1D("hNumUncorrelated_pi1520df0", "", nNum, num[0], num[1]);
  hNumUE[3][0][1] = new TH1D("hNumUncorrelated_pi1520df1", "", nNum, num[0], num[1]);
  hNumUE[3][0][2] = new TH1D("hNumUncorrelated_pi1520df2", "", nNum, num[0], num[1]);
  hNumUE[3][1][0] = new TH1D("hNumUncorrelated_ga1520df0", "", nNum, num[0], num[1]);
  hNumUE[3][1][1] = new TH1D("hNumUncorrelated_ga1520df1", "", nNum, num[0], num[1]);
  hNumUE[3][1][2] = new TH1D("hNumUncorrelated_ga1520df2", "", nNum, num[0], num[1]);
  // 2d number plot
  hNumUeVsRe[0][0][0] = new TH2D("hNumReVsUe_pi920df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[0][0][1] = new TH2D("hNumReVsUe_pi920df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[0][0][2] = new TH2D("hNumReVsUe_pi920df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[0][1][0] = new TH2D("hNumReVsUe_ga920df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[0][1][1] = new TH2D("hNumReVsUe_ga920df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[0][1][2] = new TH2D("hNumReVsUe_ga920df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[1][0][0] = new TH2D("hNumReVsUe_pi911df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[1][0][1] = new TH2D("hNumReVsUe_pi911df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[1][0][2] = new TH2D("hNumReVsUe_pi911df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[1][1][0] = new TH2D("hNumReVsUe_ga911df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[1][1][1] = new TH2D("hNumReVsUe_ga911df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[1][1][2] = new TH2D("hNumReVsUe_ga911df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[2][0][0] = new TH2D("hNumReVsUe_pi1115df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[2][0][1] = new TH2D("hNumReVsUe_pi1115df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[2][0][2] = new TH2D("hNumReVsUe_pi1115df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[2][1][0] = new TH2D("hNumReVsUe_ga1115df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[2][1][1] = new TH2D("hNumReVsUe_ga1115df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[2][1][2] = new TH2D("hNumReVsUe_ga1115df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[3][0][0] = new TH2D("hNumReVsUe_pi1520df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[3][0][1] = new TH2D("hNumReVsUe_pi1520df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[3][0][2] = new TH2D("hNumReVsUe_pi1520df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[3][1][0] = new TH2D("hNumReVsUe_ga1520df0", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[3][1][1] = new TH2D("hNumReVsUe_ga1520df1", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  hNumUeVsRe[3][1][2] = new TH2D("hNumReVsUe_ga1520df2", "", nNum, num[0], num[1], nNum, num[0], num[1]);
  // errors
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      hDfJet[iBinEt][iBinId] -> Sumw2();
      hDfRE[iBinEt][iBinId]  -> Sumw2();
      hNumRE[iBinEt][iBinId] -> Sumw2();
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hDfUE[iBinEt][iBinId][iBinDf]      -> Sumw2();
        hNumUE[iBinEt][iBinId][iBinDf]     -> Sumw2();
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> Sumw2();
      }
    }  // end id loop
  }  // end eTtrg loop
  cout << "    Defined histograms." << endl;

  // define profiles
  TProfile *pNumUeVsRe[NBinsEt][NBinsId][NBinsDf];
  // off-axis profiles
  pNumUeVsRe[0][0][0] = new TProfile("pNumUeVsRe_pi920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[0][0][1] = new TProfile("pNumUeVsRe_pi920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[0][0][2] = new TProfile("pNumUeVsRe_pi920df2", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[0][1][0] = new TProfile("pNumUeVsRe_ga920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[0][1][1] = new TProfile("pNumUeVsRe_ga920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[0][1][2] = new TProfile("pNumUeVsRe_ga920df2", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[1][0][0] = new TProfile("pNumUeVsRe_pi920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[1][0][1] = new TProfile("pNumUeVsRe_pi920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[1][0][2] = new TProfile("pNumUeVsRe_pi920df2", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[1][1][0] = new TProfile("pNumUeVsRe_ga920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[1][1][1] = new TProfile("pNumUeVsRe_ga920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[1][1][2] = new TProfile("pNumUeVsRe_ga920df2", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[2][0][0] = new TProfile("pNumUeVsRe_pi920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[2][0][1] = new TProfile("pNumUeVsRe_pi920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[2][0][2] = new TProfile("pNumUeVsRe_pi920df2", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[2][1][0] = new TProfile("pNumUeVsRe_ga920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[2][1][1] = new TProfile("pNumUeVsRe_ga920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[2][1][2] = new TProfile("pNumUeVsRe_ga920df2", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[3][0][0] = new TProfile("pNumUeVsRe_pi920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[3][0][1] = new TProfile("pNumUeVsRe_pi920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[3][0][2] = new TProfile("pNumUeVsRe_pi920df2", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[3][1][0] = new TProfile("pNumUeVsRe_ga920df0", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[3][1][1] = new TProfile("pNumUeVsRe_ga920df1", "", nNum, num[0], num[1], "S");
  pNumUeVsRe[3][1][2] = new TProfile("pNumUeVsRe_ga920df2", "", nNum, num[0], num[1], "S");
  cout << "    Defined profiles." << endl;


  // for readability
  const Double_t piOver2(TMath::PiOver2());
  const Double_t piOver4(TMath::PiOver4());
  const Double_t piOver8(piOver4 / 2.);

  // define eTtrg and dF bins
  const Double_t eTtrgMin[NBinsEt]  = {9., 9., 11., 15.};
  const Double_t eTtrgMax[NBinsEt]  = {20., 11., 15., 20.};
  const Double_t eTtrgBins[NBinsEt] = {9., 11., 15., 20.};
  const Double_t dFdownMin[NBinsDf] = {piOver4, (3. * piOver8), piOver2};
  const Double_t dFdownMax[NBinsDf] = {piOver2, (5. * piOver8), (3. * piOver4)};
  const Double_t dFupMin[NBinsDf]   = {(3. * piOver2), (11. * piOver8), (5. * piOver4)};
  const Double_t dFupMax[NBinsDf]   = {(7. * piOver4), (13. * piOver8), (3. * piOver2)};
  cout << "    Defined various bins." << endl;


  // no. of triggers and jets
  UInt_t nTrgBin[NBinsEt][NBinsId];
  UInt_t nJetRE[NBinsEt][NBinsId];
  UInt_t nJetUE[NBinsEt][NBinsId][NBinsDf];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      nTrgBin[iBinEt][iBinId] = 0;
      nJetRE[iBinEt][iBinId]  = 0;
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        nJetUE[iBinEt][iBinId][iBinDf] = 0;
      }
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

    // reset jet counters
    for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
      for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
        nJetRE[iBinEt][iBinId] = 0;
        for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
          nJetUE[iBinEt][iBinId][iBinDf] = 0;
        }
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
      const Double_t fJet    = JetPhi  -> at(iJet);
      const Double_t aJet    = JetArea -> at(iJet);
      const Double_t pTjet   = JetPt   -> at(iJet);

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
          nJetUE[0][idBin][iBinDf]++;
          nJetUE[eTbin][idBin][iBinDf]++;
        }
      }

      const Bool_t isRecoil = (TMath::Abs(dFjet - TMath::Pi()) < dFrecoil);
      if (isRecoil) {
        hDfRE[0][idBin]     -> Fill(dFjet);
        hDfRE[eTbin][idBin] -> Fill(dFjet);
        nJetRE[0][idBin]++;
        nJetRE[eTbin][idBin]++;
      }

    }  // end jet loop 2


    // fill histograms
    hNumRE[0][idBin]     -> Fill(nJetRE[0][idBin]);
    hNumRE[eTbin][idBin] -> Fill(nJetRE[eTbin][idBin]);
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hNumUE[0][idBin][iBinDf]         -> Fill(nJetUE[0][idBin][iBinDf]);
      hNumUE[eTbin][idBin][iBinDf]     -> Fill(nJetUE[eTbin][idBin][iBinDf]);
      hNumUeVsRe[0][idBin][iBinDf]     -> Fill(nJetRE[0][idBin], nJetUE[0][idBin][iBinDf]);
      hNumUeVsRe[eTbin][idBin][iBinDf] -> Fill(nJetRE[eTbin][idBin], nJetUE[eTbin][idBin][iBinDf]);
      pNumUeVsRe[0][idBin][iBinDf]     -> Fill(nJetRE[0][idBin], nJetUE[0][idBin][iBinDf]);
      pNumUeVsRe[eTbin][idBin][iBinDf] -> Fill(nJetRE[eTbin][idBin], nJetUE[eTbin][idBin][iBinDf]);
    }

  }  // end event loop

  cout << "    Finished event loop:\n"
       << "      " << nTrgBin[0][0] << " pi0 triggers,\n"
       << "      " << nTrgBin[0][1] << " gamma trigger."
       << endl;


  // get mean no. of jets per event
  Double_t avgJetRE[NBinsEt][NBinsId];
  Double_t avgJetUE[NBinsEt][NBinsId][NBinsDf];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      const Double_t nTrg  = (Double_t) nTrgBin[iBinEt][iBinId];
      const Double_t nRE   = hDfRE[iBinEt][iBinId] -> Integral();;
      const Double_t avgRE = nRE / nTrg;
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        const Double_t nUE   = hDfUE[iBinEt][iBinId][iBinDf] -> Integral();;
        const Double_t avgUE = nUE / nTrg;
        avgJetUE[iBinEt][iBinId][iBinDf] = avgUE;
      }
      avgJetRE[iBinEt][iBinId] = avgRE;

    }  // end id loop
  }  // end eTtrg loop
  cout << "    Calculated averages." << endl;

  // calculate error on means
  Double_t errRE(0.);
  Double_t errUE(0.);
  Double_t errJetRE[NBinsEt][NBinsId];
  Double_t errJetUE[NBinsEt][NBinsId][NBinsDf];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      const Double_t intRE  = hDfRE[iBinEt][iBinId] -> IntegralAndError(1, nDf, errRE);
      const Double_t nTrg   = (Double_t) nTrgBin[iBinEt][iBinId];
      const Double_t errRES = errRE / nTrg;
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        const Double_t intUE  = hDfUE[iBinEt][iBinId][iBinDf] -> IntegralAndError(1, nDf, errUE);
        const Double_t errUES = errUE / nTrg;
        errJetUE[iBinEt][iBinId][iBinDf] = errUES;
      }
      errJetRE[iBinEt][iBinId] = errRES;

    }  // end id loop
  }  // end eTtrg loop
  cout << "    Calculated errors." << endl;



  // create avg. histograms
  TH1D *hAvgRE[NBinsId];
  TH1D *hAvgUE[NBinsId][NBinsDf];
  for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
    TString sReName("hAvgRE_");
    if (iBinId == 0) sReName += "pi";
    if (iBinId == 1) sReName += "ga";

    // RE avg.s
    hAvgRE[iBinId] = new TH1D(sReName.Data(), "", NBinsEt - 1, eTtrgBins);
    hAvgRE[iBinId] -> Sumw2();
    for (UInt_t iBinEt = 1; iBinEt < NBinsEt; iBinEt++) {
      hAvgRE[iBinId] -> SetBinContent(iBinEt, avgJetRE[iBinEt][iBinId]);
      hAvgRE[iBinId] -> SetBinError(iBinEt, errJetRE[iBinEt][iBinId]);
    }

    // UE avg.s
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      TString sUeName("hAvgUE_");
      if (iBinId == 0) sUeName += "pi";
      if (iBinId == 1) sUeName += "ga";
      sUeName += "Df";
      sUeName += iBinDf;

      hAvgUE[iBinId][iBinDf] = new TH1D(sUeName.Data(), "", NBinsEt - 1, eTtrgBins);
      hAvgUE[iBinId][iBinDf] -> Sumw2();
      for (UInt_t iBinEt = 1; iBinEt < NBinsEt; iBinEt++) {
        hAvgUE[iBinId][iBinDf] -> SetBinContent(iBinEt, avgJetUE[iBinEt][iBinId][iBinDf]);
        hAvgUE[iBinId][iBinDf] -> SetBinError(iBinEt, errJetUE[iBinEt][iBinId][iBinDf]);
      }
    }  // end dF loop
  }  // end id loop
  cout << "    Made average histograms." << endl;


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

      // number histograms
      const Double_t dNnorm = (Double_t) nTrgBin[iBinEt][iBinId];
      hNumRE[iBinEt][iBinId] -> Scale(1. / dNnorm);
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hNumUE[iBinEt][iBinId][iBinDf] -> Scale(1. / dNnorm);
      }

    }  // end id loop
  }  // end eTtrg loop
  cout << "    Normalized histograms." << endl;


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColP(1);
  const UInt_t  fMarJ(20);
  const UInt_t  fMarD(1);
  const UInt_t  fMarP(20);
  const UInt_t  fLinD(2);
  const UInt_t  fFilR(3354);
  const UInt_t  fColR[NBinsId] = {859, 899};
  const UInt_t  fMarR[NBinsId] = {29, 20};
  const UInt_t  fColU[NBinsDf] = {810, 830, 850};
  const UInt_t  fMarU[NBinsDf] = {24, 26, 28};
  const UInt_t  fFilU[NBinsDf] = {3345, 3354, 3345};
  const Float_t fLab(0.03);
  const Float_t fSizR(0.074);
  const Float_t fOffsetX(1.);
  const Float_t fOffsetY(1.07);
  const Float_t fNumX[2] = {0., 20.};
  const Float_t fDfY[2]  = {0., 1.7};
  const Float_t fNumY[2] = {0.000007, 7.};
  const TString sTitleXdf("#Delta#varphi = #varphi^{jet} - #varphi^{trg}");
  const TString sTitleXnr("N^{jet}_{RE}");
  const TString sTitleXnu("N^{jet}_{UE}");
  const TString sTitleYdf("(1/N^{trg}) dN^{jet}/(d#Delta#varphi d#eta)");
  const TString sTitleYn("arb. units");
  const TString sTitleDf[NBinsId] = {"#pi^{0} trigger", "#gamma^{rich} trigger"};
  const TString sTitleN[NBinsId]  = {"#pi^{0} trigger", "#gamma^{rich} trigger"};
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

      // make x-axis title
      TString sXaxis("#color[");
      sXaxis += fColR[iBinId];
      sXaxis += "]{";
      sXaxis += sTitleXnr.Data();
      sXaxis += "}, #color[";
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        sXaxis += fColU[iBinDf];
        sXaxis += "]{";
        sXaxis += sTitleXnu.Data();
        sXaxis += "[";
        sXaxis += iBinDf;
        if (iBinDf != NBinsDf - 1)
          sXaxis += "]}, #color[";
        else
          sXaxis += "]}";
      }

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
      hDfRE[iBinEt][iBinId]  -> GetYaxis() -> SetRangeUser(fDfY[0], fDfY[1]);
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
        hDfUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetRangeUser(fDfY[0], fDfY[1]);
        hDfUE[iBinEt][iBinId][iBinDf] -> GetZaxis() -> SetLabelSize(fLab);
      }

      // number histograms
      hNumRE[iBinEt][iBinId] -> SetLineColor(fColR[iBinId]);
      hNumRE[iBinEt][iBinId] -> SetMarkerStyle(fMarR[iBinId]);
      hNumRE[iBinEt][iBinId] -> SetMarkerColor(fColR[iBinId]);
      hNumRE[iBinEt][iBinId] -> SetTitle(sTitleN[iBinId].Data());
      hNumRE[iBinEt][iBinId] -> SetTitleFont(fTxt);
      hNumRE[iBinEt][iBinId] -> GetXaxis() -> SetTitle(sXaxis.Data());
      hNumRE[iBinEt][iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
      hNumRE[iBinEt][iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hNumRE[iBinEt][iBinId] -> GetXaxis() -> CenterTitle(fCnt);
      hNumRE[iBinEt][iBinId] -> GetXaxis() -> SetLabelSize(fLab);
      hNumRE[iBinEt][iBinId] -> GetXaxis() -> SetRangeUser(fNumX[0], fNumX[1]);
      hNumRE[iBinEt][iBinId] -> GetYaxis() -> SetTitle(sTitleYn.Data());
      hNumRE[iBinEt][iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
      hNumRE[iBinEt][iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hNumRE[iBinEt][iBinId] -> GetYaxis() -> CenterTitle(fCnt);
      hNumRE[iBinEt][iBinId] -> GetYaxis() -> SetLabelSize(fLab);
      hNumRE[iBinEt][iBinId] -> GetYaxis() -> SetRangeUser(fNumY[0], fNumY[1]);
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hNumUE[iBinEt][iBinId][iBinDf] -> SetLineColor(fColU[iBinDf]);
        hNumUE[iBinEt][iBinId][iBinDf] -> SetMarkerStyle(fMarU[iBinDf]);
        hNumUE[iBinEt][iBinId][iBinDf] -> SetMarkerColor(fColU[iBinDf]);
        hNumUE[iBinEt][iBinId][iBinDf] -> SetTitle(sTitleN[iBinId].Data());
        hNumUE[iBinEt][iBinId][iBinDf] -> SetTitleFont(fTxt);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitle(sXaxis.Data());
        hNumUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleFont(fTxt);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> CenterTitle(fCnt);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetLabelSize(fLab);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetRangeUser(fNumX[0], fNumX[1]);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitle(sTitleYn.Data());
        hNumUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleFont(fTxt);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleOffset(fOffsetY);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> CenterTitle(fCnt);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetLabelSize(fLab);
        hNumUE[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetRangeUser(fNumY[0], fNumY[1]);

        // 2d number histograms
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> SetTitle(sTitleN[iBinId].Data());
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> SetTitleFont(fTxt);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitle(sTitleXnr.Data());
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleFont(fTxt);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> CenterTitle(fCnt);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetLabelSize(fLab);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetRangeUser(fNumX[0], fNumX[1]);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitle(sTitleXnu.Data());
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleFont(fTxt);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleOffset(fOffsetY);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> CenterTitle(fCnt);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetLabelSize(fLab);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetRangeUser(fNumX[0], fNumX[1]);
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetZaxis() -> SetLabelSize(fLab);

        // number profiles
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> SetLineColor(fColP);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> SetMarkerColor(fColP);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> SetMarkerStyle(fMarP);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> SetTitle(sTitleN[iBinId].Data());
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> SetTitleFont(fTxt);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitle(sTitleXnr.Data());
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleFont(fTxt);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> CenterTitle(fCnt);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetLabelSize(fLab);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetXaxis() -> SetRangeUser(fNumX[0], fNumX[1]);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitle(sTitleXnu.Data());
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleFont(fTxt);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetTitleOffset(fOffsetY);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> CenterTitle(fCnt);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetYaxis() -> SetLabelSize(fLab);
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> GetZaxis() -> SetLabelSize(fLab);
      }

    }  // end id loop
  }  // end eTtrg loop

  const UInt_t  fColAU(878);
  const UInt_t  fMarAU(25);
  const UInt_t  fSizAU(1.5);
  const UInt_t  fColAR[NBinsId]  = {858, 898};
  const UInt_t  fMarAR[NBinsId]  = {24, 24};
  const UInt_t  fSizAR[NBinsId]  = {1.5, 1.5};
  const Float_t fOffsetYA(1.17);
  const Float_t fRangeY[2]       = {0., 1.23};
  const TString sTitleA[NBinsId] = {"#pi^{0} trigger", "#gamma^{rich} trigger"};
  const TString sTitleAX("E_{T}^{trg}");
  const TString sTitleAYR("#LTN^{jet}_{RE}#GT");
  const TString sTitleAYU("#LTN^{jet}_{UE}#GT");
  for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

    // make t-axis title
    TString sYaxis("#color[");
    sYaxis += fColAR[iBinId];
    sYaxis += "]{";
    sYaxis += sTitleAYR.Data();
    sYaxis += "}, #color[";
    sYaxis += fColAU;
    sYaxis += "]{";
    sYaxis += sTitleAYU.Data();
    sYaxis += "}";

    hAvgRE[iBinId] -> SetLineColor(fColAR[iBinId]);
    hAvgRE[iBinId] -> SetMarkerStyle(fMarAR[iBinId]);
    hAvgRE[iBinId] -> SetMarkerColor(fColAR[iBinId]);
    hAvgRE[iBinId] -> SetMarkerSize(fSizAR[iBinId]);
    hAvgRE[iBinId] -> SetTitle(sTitleA[iBinId].Data());
    hAvgRE[iBinId] -> SetTitleFont(fTxt);
    hAvgRE[iBinId] -> GetXaxis() -> SetTitle(sTitleAX.Data());
    hAvgRE[iBinId] -> GetXaxis() -> SetTitleFont(fTxt);
    hAvgRE[iBinId] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hAvgRE[iBinId] -> GetXaxis() -> CenterTitle(fCnt);
    hAvgRE[iBinId] -> GetXaxis() -> SetLabelSize(fLab);
    hAvgRE[iBinId] -> GetYaxis() -> SetTitle(sYaxis.Data());
    hAvgRE[iBinId] -> GetYaxis() -> SetTitleFont(fTxt);
    hAvgRE[iBinId] -> GetYaxis() -> SetTitleOffset(fOffsetYA);
    hAvgRE[iBinId] -> GetYaxis() -> CenterTitle(fCnt);
    hAvgRE[iBinId] -> GetYaxis() -> SetLabelSize(fLab);
    hAvgRE[iBinId] -> GetYaxis() -> SetRangeUser(fRangeY[0], fRangeY[1]);
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hAvgUE[iBinId][iBinDf] -> SetLineColor(fColAU);
      hAvgUE[iBinId][iBinDf] -> SetMarkerStyle(fMarAU);
      hAvgUE[iBinId][iBinDf] -> SetMarkerColor(fColAU);
      hAvgUE[iBinId][iBinDf] -> SetMarkerSize(fSizAU);
      hAvgUE[iBinId][iBinDf] -> SetTitle(sTitleA[iBinId].Data());
      hAvgUE[iBinId][iBinDf] -> SetTitleFont(fTxt);
      hAvgUE[iBinId][iBinDf] -> GetXaxis() -> SetTitle(sTitleAX.Data());
      hAvgUE[iBinId][iBinDf] -> GetXaxis() -> SetTitleFont(fTxt);
      hAvgUE[iBinId][iBinDf] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hAvgUE[iBinId][iBinDf] -> GetXaxis() -> CenterTitle(fCnt);
      hAvgUE[iBinId][iBinDf] -> GetXaxis() -> SetLabelSize(fLab);
      hAvgUE[iBinId][iBinDf] -> GetYaxis() -> SetTitle(sYaxis.Data());
      hAvgUE[iBinId][iBinDf] -> GetYaxis() -> SetTitleFont(fTxt);
      hAvgUE[iBinId][iBinDf] -> GetYaxis() -> SetTitleOffset(fOffsetYA);
      hAvgUE[iBinId][iBinDf] -> GetYaxis() -> CenterTitle(fCnt);
      hAvgUE[iBinId][iBinDf] -> GetYaxis() -> SetLabelSize(fLab);
      hAvgUE[iBinId][iBinDf] -> GetYaxis() -> SetRangeUser(fRangeY[0], fRangeY[1]);
    }

  }  // end id loop
  cout << "    Set styles." << endl;


  // make lines
  const UInt_t fLinLV(2);
  const UInt_t fLinLE(1);
  const UInt_t fSizLV(1);
  const UInt_t fSizLE(2);
  const UInt_t fColLU(872);
  const UInt_t fColLR[NBinsId] = {852, 902};

  TLine *lAvgEtRE[NBinsId];
  TLine *lAvgEtUE[NBinsId];
  TLine *lErrEtRE[NBinsId][2];
  TLine *lErrEtUE[NBinsId][2];
  for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
    const Double_t start = eTtrgBins[0];
    const Double_t stop  = eTtrgBins[NBinsEt - 1];
    const Double_t valRE = avgJetRE[0][iBinId];
    const Double_t errRE = errJetRE[0][iBinId];
    const Double_t valUE = avgJetUE[0][iBinId][0];
    const Double_t errUE = errJetUE[0][iBinId][0];
    const Double_t hiRE  = valRE + errRE;
    const Double_t loRE  = valRE - errRE;
    const Double_t hiUE  = valUE + errUE;
    const Double_t loUE  = valUE - errUE;

    // recoil
    lAvgEtRE[iBinId]    = new TLine(start, valRE, stop, valRE);
    lErrEtRE[iBinId][0] = new TLine(start, loRE, stop, loRE);
    lErrEtRE[iBinId][1] = new TLine(start, hiRE, stop, hiRE);
    lAvgEtRE[iBinId]    -> SetLineColor(fColLR[iBinId]);
    lAvgEtRE[iBinId]    -> SetLineStyle(fLinLV);
    lAvgEtRE[iBinId]    -> SetLineWidth(fSizLV);
    lErrEtRE[iBinId][0] -> SetLineColor(fColLR[iBinId]);
    lErrEtRE[iBinId][0] -> SetLineStyle(fLinLE);
    lErrEtRE[iBinId][0] -> SetLineWidth(fSizLE);
    lErrEtRE[iBinId][1] -> SetLineColor(fColLR[iBinId]);
    lErrEtRE[iBinId][1] -> SetLineStyle(fLinLE);
    lErrEtRE[iBinId][1] -> SetLineWidth(fSizLE);

    // recoil
    lAvgEtUE[iBinId]    = new TLine(start, valUE, stop, valUE);
    lErrEtUE[iBinId][0] = new TLine(start, loUE, stop, loUE);
    lErrEtUE[iBinId][1] = new TLine(start, hiUE, stop, hiUE);
    lAvgEtUE[iBinId]    -> SetLineColor(fColLU);
    lAvgEtUE[iBinId]    -> SetLineStyle(fLinLV);
    lAvgEtUE[iBinId]    -> SetLineWidth(fSizLV);
    lErrEtUE[iBinId][0] -> SetLineColor(fColLU);
    lErrEtUE[iBinId][0] -> SetLineStyle(fLinLE);
    lErrEtUE[iBinId][0] -> SetLineWidth(fSizLE);
    lErrEtUE[iBinId][1] -> SetLineColor(fColLU);
    lErrEtUE[iBinId][1] -> SetLineStyle(fLinLE);
    lErrEtUE[iBinId][1] -> SetLineWidth(fSizLE);
  }
  cout << "    Made lines." << endl;


  // make labels
  const UInt_t  fColT(0);
  const UInt_t  fFilT(0);
  const UInt_t  fLinT(0);
  const UInt_t  fAlign(12);
  const UInt_t  nDec(3);
  const Float_t xLeg[2]  = {0.3, 0.5};
  const Float_t yLeg[2]  = {0.7, 0.9};
  const Float_t xInfo[2] = {0.5, 0.7};
  const Float_t yInfo[2] = {0.7, 0.9};
  const Float_t xText[2] = {0.7, 0.9};
  const Float_t yText[2] = {0.7, 0.9};
  const TString sSystem("pp-collsions, #sqrt{s} = 200 GeV");
  const TString sTrigger("E_{T}^{trg} #in (");
  const TString sJet("anti-k_{T}, R = ");
  const TString sChrg("#bf{charged jets}");
  const TString sFull("#bf{full jets}");

  TPaveText *pAvgEt[NBinsId];
  TPaveText *pInfo[NBinsEt][NBinsId];
  TPaveText *pAvgNums[NBinsEt][NBinsId];
  TLegend   *lDeltaPhi[NBinsId];
  TLegend   *lNums[NBinsId];
  TLegend   *lAvg[NBinsId];
  for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {

    // average legend
    lAvg[iBinId] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[2]);
    lAvg[iBinId] -> SetFillColor(fColT);
    lAvg[iBinId] -> SetFillStyle(fFilT);
    lAvg[iBinId] -> SetLineColor(fColT);
    lAvg[iBinId] -> SetLineStyle(fLinT);
    lAvg[iBinId] -> SetTextFont(fTxt);
    lAvg[iBinId] -> AddEntry(hAvgRE[iBinId], "RE");
    lAvg[iBinId] -> AddEntry(hAvgUE[iBinId][0], "UE");

    // delta-phi legend
    lDeltaPhi[iBinId] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lDeltaPhi[iBinId] -> SetFillColor(fColT);
    lDeltaPhi[iBinId] -> SetFillStyle(fFilT);
    lDeltaPhi[iBinId] -> SetLineColor(fColT);
    lDeltaPhi[iBinId] -> SetLineStyle(fLinT);
    lDeltaPhi[iBinId] -> SetTextFont(fTxt);
    lDeltaPhi[iBinId] -> AddEntry(hDfJet[0][iBinId], "all jets");
    lDeltaPhi[iBinId] -> AddEntry(hDfRE[0][iBinId], "RE");
    lDeltaPhi[iBinId] -> AddEntry(hDfUE[0][iBinId][0], "UE[0]");
    lDeltaPhi[iBinId] -> AddEntry(hDfUE[0][iBinId][1], "UE[1]");
    lDeltaPhi[iBinId] -> AddEntry(hDfUE[0][iBinId][2], "UE[2]");

    // number legend
    lNums[iBinId] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lNums[iBinId] -> SetFillColor(fColT);
    lNums[iBinId] -> SetFillStyle(fFilT);
    lNums[iBinId] -> SetLineColor(fColT);
    lNums[iBinId] -> SetLineStyle(fLinT);
    lNums[iBinId] -> SetTextFont(fTxt);
    lNums[iBinId] -> AddEntry(hNumRE[0][iBinId], "RE");
    lNums[iBinId] -> AddEntry(hNumUE[0][iBinId][0], "UE[0]");
    lNums[iBinId] -> AddEntry(hNumUE[0][iBinId][1], "UE[1]");
    lNums[iBinId] -> AddEntry(hNumUE[0][iBinId][2], "UE[2]");


    // average info
    pAvgEt[iBinId] = new TPaveText(xInfo[0], yInfo[0], xInfo[1], yInfo[1], "NDC NB");
    pAvgEt[iBinId] -> SetFillColor(fColT);
    pAvgEt[iBinId] -> SetFillStyle(fFilT);
    pAvgEt[iBinId] -> SetLineColor(fColT);
    pAvgEt[iBinId] -> SetLineStyle(fLinT);
    pAvgEt[iBinId] -> SetTextFont(fTxt);
    pAvgEt[iBinId] -> SetTextAlign(fAlign);

    TString sRAVraw("");
    TString sUAVraw("");
    TString sRAEraw("");
    TString sUAEraw("");
    TString sRAVtxt("");
    TString sUAVtxt("");
    TString sRAEtxt("");
    TString sUAEtxt("");
    sRAVraw += avgJetRE[0][iBinId];
    sUAVraw += avgJetUE[0][iBinId][0];
    sRAEraw += errJetRE[0][iBinId];
    sUAEraw += errJetUE[0][iBinId][0];

    const UInt_t nRAVraw = sRAVraw.First(".");
    const UInt_t nUAVraw = sUAVraw.First(".");
    const UInt_t nRAEraw = sRAEraw.First(".");
    const UInt_t nUAEraw = sUAEraw.First(".");
    const UInt_t nRAVtxt = (nRAVraw + nDec) + 1;
    const UInt_t nUAVtxt = (nUAVraw + nDec) + 1;
    const UInt_t nRAEtxt = (nRAEraw + nDec) + 1;
    const UInt_t nUAEtxt = (nUAEraw + nDec) + 1;
    sRAVtxt.Append(sRAVraw.Data(), nRAVtxt);
    sUAVtxt.Append(sUAVraw.Data(), nUAVtxt);
    sRAEtxt.Append(sRAEraw.Data(), nRAEtxt);
    sUAEtxt.Append(sUAEraw.Data(), nUAEtxt);

    TString sREavg("#color[");
    TString sUEavg("#color[");
    sREavg += fColLR[iBinId];
    sREavg += "]{#LT#LTN^{jet}_{RE}#GT#GT = ";
    sREavg += sRAVtxt.Data();
    sREavg += " #pm ";
    sREavg += sRAEtxt.Data();
    sREavg += "}";
    sUEavg += fColLU;
    sUEavg += "]{#LT#LTN^{jet}_{UE}#GT#GT = ";
    sUEavg += sUAVtxt.Data();
    sUEavg += " #pm ";
    sUEavg += sUAEtxt.Data();
    sUEavg += "}";

    TString sJetInfo(sJet.Data());
    TString sJetType("");
    sJetInfo += sRes.Data();
    if (JetType == 0)
      sJetType += sChrg.Data();
    else
      sJetType += sFull.Data();

    pAvgEt[iBinId] -> AddText(sSystem.Data());
    pAvgEt[iBinId] -> AddText(sJetInfo.Data());
    pAvgEt[iBinId] -> AddText(sJetType.Data());
    pAvgEt[iBinId] -> AddText(sREavg.Data());
    pAvgEt[iBinId] -> AddText(sUEavg.Data());

    for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {

      // general info
      TString sEtTrg(sTrigger.Data());
      sEtTrg += eTtrgMin[iBinEt];
      sEtTrg += ", ";
      sEtTrg += eTtrgMax[iBinEt];
      sEtTrg += ") GeV";

      pInfo[iBinEt][iBinId] = new TPaveText(xInfo[0], yInfo[0], xInfo[1], yInfo[1], "NDC NB");
      pInfo[iBinEt][iBinId] -> SetFillColor(fColT);
      pInfo[iBinEt][iBinId] -> SetFillStyle(fFilT);
      pInfo[iBinEt][iBinId] -> SetLineColor(fColT);
      pInfo[iBinEt][iBinId] -> SetLineStyle(fLinT);
      pInfo[iBinEt][iBinId] -> SetTextFont(fTxt);
      pInfo[iBinEt][iBinId] -> SetTextAlign(fAlign);
      pInfo[iBinEt][iBinId] -> AddText(sSystem.Data());
      pInfo[iBinEt][iBinId] -> AddText(sEtTrg.Data());
      pInfo[iBinEt][iBinId] -> AddText(sJetInfo.Data());
      pInfo[iBinEt][iBinId] -> AddText(sJetType.Data());

      // average no.s
      pAvgNums[iBinEt][iBinId] = new TPaveText(xText[0], yText[0], xText[1], yText[1], "NDC NB");
      pAvgNums[iBinEt][iBinId] -> SetFillColor(fColT);
      pAvgNums[iBinEt][iBinId] -> SetFillStyle(fFilT);
      pAvgNums[iBinEt][iBinId] -> SetLineColor(fColT);
      pAvgNums[iBinEt][iBinId] -> SetLineStyle(fLinT);
      pAvgNums[iBinEt][iBinId] -> SetTextFont(fTxt);
      pAvgNums[iBinEt][iBinId] -> SetTextAlign(fAlign);

      TString sRVraw("");
      TString sREraw("");
      TString sRVtxt("");
      TString sREtxt("");
      sRVraw += avgJetRE[iBinEt][iBinId];
      sREraw += errJetRE[iBinEt][iBinId];

      const UInt_t nRVraw = sRVraw.First(".");
      const UInt_t nREraw = sREraw.First(".");
      const UInt_t nRVtxt = (nRVraw + nDec) + 1;
      const UInt_t nREtxt = (nREraw + nDec) + 1;
      sRVtxt.Append(sRVraw.Data(), nRVtxt);
      sREtxt.Append(sREraw.Data(), nREtxt);

      TString sRavg("#color[");
      sRavg += fColR[iBinId];
      sRavg += "]{#LTN^{jet}_{RE}#GT = ";
      sRavg += sRVtxt.Data();
      sRavg += " #pm ";
      sRavg += sREtxt.Data();
      sRavg += "}";
      pAvgNums[iBinEt][iBinId] -> AddText(sRavg.Data());
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {

        TString sUVraw("");
        TString sUEraw("");
        TString sUVtxt("");
        TString sUEtxt("");
        sUVraw += avgJetUE[iBinEt][iBinId][iBinDf];
        sUEraw += avgJetUE[iBinEt][iBinId][iBinDf];

        const UInt_t nUVraw = sUVraw.First(".");
        const UInt_t nUEraw = sUEraw.First(".");
        const UInt_t nUVtxt = (nUVraw + nDec) + 1;
        const UInt_t nUEtxt = (nUEraw + nDec) + 1;
        sUVtxt.Append(sUVraw.Data(), nUVtxt);
        sUEtxt.Append(sUEraw.Data(), nUEtxt);

        TString sUavg("#color[");
        sUavg += fColU[iBinDf];
        sUavg += "]{#LTN^{jet}_{UE}[";
        sUavg += iBinDf;
        sUavg += "]#GT = ";
        sUavg += sUVtxt.Data();
        sUavg += " #pm ";
        sUavg += sUEtxt.Data();
        sUavg += "}";
        pAvgNums[iBinEt][iBinId] -> AddText(sUavg.Data());

      }  // end dF loop
    }  // end eTtrg loop
  }  // end id loop
  cout << "    Made labels." << endl;


  // create directories and save histograms
  const TString sBinsEt[NBinsEt] = {"eTtrg920", "eTtrg911", "eTtrg1115", "eTtrg1520"};
  const TString sBinsId[NBinsId] = {"pi0", "gamma"};

  TDirectory *dBinsEt[NBinsEt];
  TDirectory *dBinsId[NBinsEt][NBinsId];
  for (UInt_t iBinEt = 0; iBinEt < NBinsEt; iBinEt++) {
    dBinsEt[iBinEt] = (TDirectory*) fOutput -> mkdir(sBinsEt[iBinEt].Data());
    dBinsEt[iBinEt]  -> cd();
    for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
      dBinsId[iBinEt][iBinId] = (TDirectory*) dBinsEt[iBinEt] -> mkdir(sBinsId[iBinId].Data());
      dBinsId[iBinEt][iBinId] -> cd();
      hDfJet[iBinEt][iBinId]  -> Write();
      hDfRE[iBinEt][iBinId]   -> Write();
      hNumRE[iBinEt][iBinId]  -> Write();
      for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
        hDfUE[iBinEt][iBinId][iBinDf]      -> Write();
        hNumUE[iBinEt][iBinId][iBinDf]     -> Write();
        hNumUeVsRe[iBinEt][iBinId][iBinDf] -> Write();
        pNumUeVsRe[iBinEt][iBinId][iBinDf] -> Write();
      }
    }  // end id loop
  }  // end eTtrg loop
  for (UInt_t iBinId = 0; iBinId < NBinsId; iBinId++) {
    hAvgRE[iBinId] -> Write();
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hAvgUE[iBinId][iBinDf] -> Write();
    }
  }  // end id loop
  cout << "    Made directories." << endl;


  // draw plots
  const UInt_t  widthAvg(1500);
  const UInt_t  heightAvg(750);
  const UInt_t  widthDf(1500);
  const UInt_t  heightDf(750);
  const UInt_t  widthNum(1500);
  const UInt_t  heightNum(750);
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
  const Float_t xPadAvg[NBinsId + 2] = {0., 0.5, 0.5, 1.};
  const Float_t yPadAvg[NBinsId + 2] = {0., 1., 0., 1.};
  const Float_t xPadDf[NBinsId + 2]  = {0., 0.5, 0.5, 1.};
  const Float_t yPadDf[NBinsId + 2]  = {0., 1., 0., 1.};
  const Float_t xPadN[NBinsId + 2]   = {0., 0.5, 0.5, 1.};
  const Float_t yPadN[NBinsId + 2]   = {0., 1., 0., 1.};
  const TString sAvgJets("cAvgJets");
  const TString sDfCanvas[NBinsEt]   = {"cDfJet_et920", "cDfJet_et911", "cDfJet_et1115", "cDfJet_et1520"};
  const TString sNumCanvas[NBinsEt]  = {"cNum_et920", "cNum_et911", "cNum_et1115", "cNum_et1520"};
  const TString sNumCanvas2[NBinsEt] = {"cUeVsRe_et920", "cUeVsRe_et911", "cUeVsRe_et1115", "cUeVsRe_et1520"};
  const TString sAvgPads[NBinsId]    = {"pPi0", "pGamma"};
  const TString sDfPads[NBinsId]     = {"pPi0", "pGamma"};
  const TString sNumPads[NBinsId]    = {"pPi0", "pGamma"};
  const TString sNumPads2[NBinsId]   = {"pPi0", "pGamma"};

  TPad    *pDfJet[NBinsEt][NBinsId];
  TPad    *pNumJet[NBinsEt][NBinsId];
  TPad    *pNumJet2[NBinsEt][NBinsId][NBinsDf];
  TPad    *pAvgJets[NBinsId];
  TCanvas *cDfJet[NBinsEt];
  TCanvas *cNumJet[NBinsEt];
  TCanvas *cNumJet2[NBinsEt][NBinsDf];
  TCanvas *cAvgJets;

  fOutput -> cd();
  cAvgJets    = new TCanvas(sAvgJets.Data(), "", widthAvg, heightAvg);
  pAvgJets[0] = new TPad(sAvgPads[0].Data(), "", xPadAvg[0], yPadAvg[0], xPadAvg[1], yPadAvg[1]);
  pAvgJets[1] = new TPad(sAvgPads[1].Data(), "", xPadAvg[2], yPadAvg[2], xPadAvg[3], yPadAvg[3]);
  pAvgJets[0]    -> SetGrid(grid, grid);
  pAvgJets[0]    -> SetFrameBorderMode(frame);
  pAvgJets[0]    -> SetRightMargin(marginR);
  pAvgJets[0]    -> SetLeftMargin(marginL);
  pAvgJets[0]    -> SetTopMargin(marginT);
  pAvgJets[0]    -> SetBottomMargin(marginB);
  pAvgJets[1]    -> SetGrid(grid, grid);
  pAvgJets[1]    -> SetFrameBorderMode(frame);
  pAvgJets[1]    -> SetRightMargin(marginR);
  pAvgJets[1]    -> SetLeftMargin(marginL);
  pAvgJets[1]    -> SetTopMargin(marginT);
  pAvgJets[1]    -> SetBottomMargin(marginB);
  cAvgJets       -> cd();
  pAvgJets[0]    -> Draw();
  pAvgJets[1]    -> Draw();
  pAvgJets[0]    -> cd();
  hAvgRE[0]      -> Draw();
  hAvgUE[0][0]   -> Draw("same");
  lAvgEtRE[0]    -> Draw();
  lErrEtRE[0][0] -> Draw();
  lErrEtRE[0][1] -> Draw();
  lAvgEtUE[0]    -> Draw();
  lErrEtUE[0][0] -> Draw();
  lErrEtUE[0][1] -> Draw();
  lAvg[0]        -> Draw();
  pAvgEt[0]      -> Draw();
  pAvgJets[1]    -> cd();
  hAvgRE[1]      -> Draw();
  hAvgUE[1][0]   -> Draw("same");
  lAvgEtRE[1]    -> Draw();
  lErrEtRE[1][0] -> Draw();
  lErrEtRE[1][1] -> Draw();
  lAvgEtUE[1]    -> Draw();
  lErrEtUE[1][0] -> Draw();
  lErrEtUE[1][1] -> Draw();
  lAvg[1]        -> Draw();
  pAvgEt[1]      -> Draw();
  cAvgJets       -> Write();
  cAvgJets       -> Close();

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
    lDeltaPhi[0]        -> Draw();
    pInfo[iBinEt][0]    -> Draw();
    pAvgNums[iBinEt][0] -> Draw();
    pDfJet[iBinEt][1]   -> cd();
    hDfJet[iBinEt][1]   -> Draw();
    hDfRE[iBinEt][1]    -> Draw("SAME LF HIST");
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hDfUE[iBinEt][1][iBinDf] -> Draw("SAME LF HIST");
    }
    lDeltaPhi[1]        -> Draw();
    pInfo[iBinEt][1]    -> Draw();
    pAvgNums[iBinEt][1] -> Draw();
    cDfJet[iBinEt]      -> Write();
    cDfJet[iBinEt]      -> Close();

    // number plots
    dBinsEt[iBinEt] -> cd();
    cNumJet[iBinEt]    = new TCanvas(sNumCanvas[iBinEt].Data(), "", widthNum, heightNum);
    pNumJet[iBinEt][0] = new TPad(sNumPads[0].Data(), "", xPadN[0], yPadN[0], xPadN[1], yPadN[1]);
    pNumJet[iBinEt][1] = new TPad(sNumPads[1].Data(), "", xPadN[2], yPadN[2], xPadN[3], yPadN[3]);
    pNumJet[iBinEt][0] -> SetGrid(grid, grid);
    pNumJet[iBinEt][0] -> SetLogy(log);
    pNumJet[iBinEt][0] -> SetFrameBorderMode(frame);
    pNumJet[iBinEt][0] -> SetRightMargin(marginR);
    pNumJet[iBinEt][0] -> SetLeftMargin(marginL);
    pNumJet[iBinEt][0] -> SetTopMargin(marginT);
    pNumJet[iBinEt][0] -> SetBottomMargin(marginB);
    pNumJet[iBinEt][1] -> SetGrid(grid, grid);
    pNumJet[iBinEt][1] -> SetLogy(log);
    pNumJet[iBinEt][1] -> SetFrameBorderMode(frame);
    pNumJet[iBinEt][1] -> SetRightMargin(marginR);
    pNumJet[iBinEt][1] -> SetLeftMargin(marginL);
    pNumJet[iBinEt][1] -> SetTopMargin(marginT);
    pNumJet[iBinEt][1] -> SetBottomMargin(marginB);
    cNumJet[iBinEt]    -> cd();
    pNumJet[iBinEt][0] -> Draw();
    pNumJet[iBinEt][1] -> Draw();
    pNumJet[iBinEt][0] -> cd();
    hNumRE[iBinEt][0]  -> Draw();
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hNumUE[iBinEt][0][iBinDf] -> Draw("same");
    }
    lNums[0]            -> Draw();
    pInfo[iBinEt][0]    -> Draw();
    pAvgNums[iBinEt][0] -> Draw();
    pNumJet[iBinEt][1]  -> cd();
    hNumRE[iBinEt][1]   -> Draw();
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      hNumUE[iBinEt][1][iBinDf] -> Draw("same");
    }
    lNums[1]            -> Draw();
    pInfo[iBinEt][1]    -> Draw();
    pAvgNums[iBinEt][1] -> Draw();
    cNumJet[iBinEt]     -> Write();
    cNumJet[iBinEt]     -> Close();

    // 2d number plots
    for (UInt_t iBinDf = 0; iBinDf < NBinsDf; iBinDf++) {
      TString sCanvasName(sNumCanvas2[iBinEt].Data());
      sCanvasName += "df";
      sCanvasName += iBinDf;

      dBinsEt[iBinEt] -> cd();
      cNumJet2[iBinEt][iBinDf]    = new TCanvas(sCanvasName.Data(), "", widthNum, heightNum);
      pNumJet2[iBinEt][0][iBinDf] = new TPad(sNumPads2[0].Data(), "", xPadN[0], yPadN[0], xPadN[1], yPadN[1]);
      pNumJet2[iBinEt][1][iBinDf] = new TPad(sNumPads2[1].Data(), "", xPadN[2], yPadN[2], xPadN[3], yPadN[3]);
      pNumJet2[iBinEt][0][iBinDf]   -> SetGrid(grid, grid);
      pNumJet2[iBinEt][0][iBinDf]   -> SetLogz(log);
      pNumJet2[iBinEt][0][iBinDf]   -> SetFrameBorderMode(frame);
      pNumJet2[iBinEt][0][iBinDf]   -> SetRightMargin(marginRcolz);
      pNumJet2[iBinEt][0][iBinDf]   -> SetLeftMargin(marginL);
      pNumJet2[iBinEt][0][iBinDf]   -> SetTopMargin(marginT);
      pNumJet2[iBinEt][0][iBinDf]   -> SetBottomMargin(marginB);
      pNumJet2[iBinEt][1][iBinDf]   -> SetGrid(grid, grid);
      pNumJet2[iBinEt][1][iBinDf]   -> SetLogz(log);
      pNumJet2[iBinEt][1][iBinDf]   -> SetFrameBorderMode(frame);
      pNumJet2[iBinEt][1][iBinDf]   -> SetRightMargin(marginRcolz);
      pNumJet2[iBinEt][1][iBinDf]   -> SetLeftMargin(marginL);
      pNumJet2[iBinEt][1][iBinDf]   -> SetTopMargin(marginT);
      pNumJet2[iBinEt][1][iBinDf]   -> SetBottomMargin(marginB);
      cNumJet2[iBinEt][iBinDf]      -> cd();
      pNumJet2[iBinEt][0][iBinDf]   -> Draw();
      pNumJet2[iBinEt][1][iBinDf]   -> Draw();
      pNumJet2[iBinEt][0][iBinDf]   -> cd();
      hNumUeVsRe[iBinEt][0][iBinDf] -> Draw("colz");
      pNumUeVsRe[iBinEt][0][iBinDf] -> Draw("same");
      pNumJet2[iBinEt][1][iBinDf]   -> cd();
      hNumUeVsRe[iBinEt][1][iBinDf] -> Draw("colz");
      pNumUeVsRe[iBinEt][1][iBinDf] -> Draw("same");
      cNumJet2[iBinEt][iBinDf]      -> Write();
      cNumJet2[iBinEt][iBinDf]      -> Close();
    }

  }  // end eTtrg loop
  cout << "    Drew plots." << endl;

  // save and close
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Calculation finished!\n" << endl;

}

// End ------------------------------------------------------------------------
