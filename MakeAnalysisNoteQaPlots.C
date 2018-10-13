// 'MakeAnalysisNoteQaPlots.C'
// Derek Anderson
// 10.04.2018
//
// Use this script to make various QA
// plots for the neutral-triggered
// recoil jet analysis note.


#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TDirectory.h"


using namespace std;


// global constants
static const UInt_t NTrkMax(5000);
static const UInt_t NTwrMax(5000);
static const UInt_t NMatchMax(10);
static const UInt_t NHotTwr(41);
static const UInt_t NBadRuns(45);
static const UInt_t NTrgCuts(9);
static const UInt_t NTrgBins(4);
static const UInt_t NTrkBins(4);
static const UInt_t NTrgTsp(2);
static const UInt_t NVtx(4);



void MakeAnalysisNoteQaPlots(const Bool_t isInBatchMode=false) {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Making QA plots..." << endl;


  // io parameters
  const TString sOutput("test.root");
  const TString sInput("input/pp200r9.merge.root");

  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);

  // trigger paramters
  const Int_t    adcMax(6004);
  const Double_t eStrMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t hTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  const Double_t eTbinMin[NTrgBins] = {9., 9., 11., 15.};
  const Double_t eTbinMax[NTrgBins] = {20., 11., 15., 20.};
  const Double_t tspPi0[NTrgTsp]    = {0., 0.08};
  const Double_t tspGam[NTrgTsp]    = {0.2, 0.6};

  // track parameters
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkMin(1.2);
  const Double_t pTtrkMax(20.);
  const Double_t pTbinMin[NTrkBins] = {0.2, 0.2, 2., 5.};
  const Double_t pTbinMax[NTrkBins] = {20., 2., 5., 20.};



  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fInput) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input tree
  TTree *tInput;
  fInput -> GetObject("Gfmtodst", tInput);
  if (!tInput) {
    cerr << "PANIC: couldn't grab input tree!" << endl;
    return;
  }
  cout << "    Grabbed tree." << endl;


  // declare input leaf addresses
  UInt_t   fUniqueID;
  UInt_t   fBits;
  Long64_t runNumber;
  Long64_t eventNumber;
  Int_t    trigID;
  Int_t    nGlobalTracks;
  Int_t    nPrimaryTracks;
  Int_t    refMult;
  Double_t vpdVz;
  Double_t xVertex;
  Double_t yVertex;
  Double_t zVertex;
  Double_t bbcZVertex;
  Double_t zdcCoincidenceRate;
  Double_t bbcCoincidenceRate;
  Double_t backgroundRate;
  Double_t bbcBlueBackgroundRate;
  Double_t bbcYellowBackgroundRate;
  Double_t refMultPos;
  Double_t refMultNeg;
  Double_t bTOFTrayMultiplicity;
  Int_t    nVerticies;
  Double_t MagF;
  Double_t VrtxRank;
  Int_t    FlagEvent_TrgTrkMisMtch;
  Float_t  Etsp;
  Int_t    ETwrdidT;
  Int_t    ETwradc11;
  Float_t  ETwreneT0;
  Float_t  ETwreT;
  Float_t  ETwrENET0;
  Float_t  ETwrphT;
  Float_t  ETwrPTower;
  Float_t  ETwrpidTower;
  Int_t    ETwrmoduleT;
  Float_t  EClustEneT0;
  Float_t  EClustetav1;
  Float_t  EClustphiv1;
  Float_t  EEstrpen01;
  Float_t  EEstrpen02;
  Float_t  EEstrpen03;
  Float_t  EEstrpen0;
  Float_t  EEstrpen1;
  Float_t  EEstrpen2;
  Float_t  EEstrpen3;
  Float_t  EEstrpen4;
  Float_t  EEstrpen5;
  Float_t  EEstrpen6;
  Float_t  EEstrpen7;
  Float_t  EEstrpen8;
  Float_t  EEstrpen9;
  Float_t  EEstrpen10;
  Float_t  EEstrpen11;
  Float_t  EEstrpen12;
  Float_t  EEstrpen13;
  Float_t  EEstrpen14;
  Float_t  EEstrpen15;
  Int_t    ETwrdidE;
  Float_t  EPstripenp01;
  Float_t  EPstripenp02;
  Float_t  EPstripenp03;
  Float_t  EPstripenp0;
  Float_t  EPstripenp1;
  Float_t  EPstripenp2;
  Float_t  EPstripenp3;
  Float_t  EPstripenp4;
  Float_t  EPstripenp5;
  Float_t  EPstripenp6;
  Float_t  EPstripenp7;
  Float_t  EPstripenp8;
  Float_t  EPstripenp9;
  Float_t  EPstripenp10;
  Float_t  EPstripenp11;
  Float_t  EPstripenp12;
  Float_t  EPstripenp13;
  Float_t  EPstripenp14;
  Float_t  EPstripenp15;
  Float_t  EclustEnnq1;
  Float_t  EclustEnnq20;
  Float_t  EclustEnnq19;
  Float_t  EclustEnpq1;
  Float_t  EclustEnpq20;
  Float_t  EclustEnpq19;
  Float_t  EclustEnpq21;
  Int_t    PrimaryTrackArray_;
  UInt_t   PrimaryTrackArray_fUniqueID[NTrkMax];
  UInt_t   PrimaryTrackArray_fBits[NTrkMax];
  Double_t PrimaryTrackArray_nHitsFit[NTrkMax];
  Double_t PrimaryTrackArray_nHitsPoss[NTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[NTrkMax];
  Double_t PrimaryTrackArray_pZ[NTrkMax];
  Double_t PrimaryTrackArray_pX[NTrkMax];
  Double_t PrimaryTrackArray_pY[NTrkMax];
  Double_t PrimaryTrackArray_pT[NTrkMax];
  Double_t PrimaryTrackArray_dEdx[NTrkMax];
  Double_t PrimaryTrackArray_charge[NTrkMax];
  Double_t PrimaryTrackArray_tofBeta[NTrkMax];
  Double_t PrimaryTrackArray_eta[NTrkMax];
  Double_t PrimaryTrackArray_phi[NTrkMax];
  Double_t PrimaryTrackArray_nSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_nSigPion[NTrkMax];
  Double_t PrimaryTrackArray_nSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_nSigProton[NTrkMax];
  Double_t PrimaryTrackArray_dcag[NTrkMax];
  Double_t PrimaryTrackArray_nHits[NTrkMax];
  Double_t PrimaryTrackArray_dEdxHits[NTrkMax];
  Double_t PrimaryTrackArray_firstZPoint[NTrkMax];
  Double_t PrimaryTrackArray_lastZPoint[NTrkMax];
  Double_t PrimaryTrackArray_tofSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_tofSigPion[NTrkMax];
  Double_t PrimaryTrackArray_tofSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_tofSigProton[NTrkMax];
  Double_t PrimaryTrackArray_timeOfflight[NTrkMax];
  Double_t PrimaryTrackArray_pathLength[NTrkMax];
  Int_t    PrimaryTrackArray_trkIndex[NTrkMax];
  Int_t    TowerArray_;
  UInt_t   TowerArray_fUniqueID[NTwrMax];
  UInt_t   TowerArray_fBits[NTwrMax];
  Int_t    TowerArray_TwrId[NTwrMax];
  Float_t  TowerArray_TwrEng[NTwrMax];
  Float_t  TowerArray_TwrEta[NTwrMax];
  Float_t  TowerArray_TwrPhi[NTwrMax];
  Float_t  TowerArray_TwrADC[NTwrMax];
  Float_t  TowerArray_TwrPed[NTwrMax];
  Float_t  TowerArray_TwrRMS[NTwrMax];
  Int_t    TowerArray_TwrMatchIdnex[NTwrMax];
  Int_t    TowerArray_NoOfmatchedTrk[NTwrMax];
  Float_t  TowerArray_TwrMatchP[NTwrMax];
  Float_t  TowerArray_TwrPx[NTwrMax];
  Float_t  TowerArray_TwrPy[NTwrMax];
  Float_t  TowerArray_TwrPz[NTwrMax];
  Int_t    TowerArray_fNAssocTracks[NTwrMax];
  Int_t    TowerArray_fMatchedTracksArray_[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_P[NTwrMax][NMatchMax];

  // declare input branches
  TBranch *bEventList_fUniqueID;
  TBranch *bEventList_fBits;
  TBranch *bEventList_runNumber;
  TBranch *bEventList_eventNumber;
  TBranch *bEventList_trigID;
  TBranch *bEventList_nGlobalTracks;
  TBranch *bEventList_nPrimaryTracks;
  TBranch *bEventList_refMult;
  TBranch *bEventList_vpdVz;
  TBranch *bEventList_xVertex;
  TBranch *bEventList_yVertex;
  TBranch *bEventList_zVertex;
  TBranch *bEventList_bbcZVertex;
  TBranch *bEventList_zdcCoincidenceRate;
  TBranch *bEventList_bbcCoincidenceRate;
  TBranch *bEventList_backgroundRate;
  TBranch *bEventList_bbcBlueBackgroundRate;
  TBranch *bEventList_bbcYellowBackgroundRate;
  TBranch *bEventList_refMultPos;
  TBranch *bEventList_refMultNeg;
  TBranch *bEventList_bTOFTrayMultiplicity;
  TBranch *bEventList_nVerticies;
  TBranch *bEventList_MagF;
  TBranch *bEventList_VrtxRank;
  TBranch *bEventList_FlagEvent_TrgTrkMisMtch;
  TBranch *bEventList_Etsp;
  TBranch *bEventList_ETwrdidT;
  TBranch *bEventList_ETwradc11;
  TBranch *bEventList_ETwreneT0;
  TBranch *bEventList_ETwreT;
  TBranch *bEventList_ETwrENET0;
  TBranch *bEventList_ETwrphT;
  TBranch *bEventList_ETwrPTower;
  TBranch *bEventList_ETwrpidTower;
  TBranch *bEventList_ETwrmoduleT;
  TBranch *bEventList_EClustEneT0;
  TBranch *bEventList_EClustetav1;
  TBranch *bEventList_EClustphiv1;
  TBranch *bEventList_EEstrpen01;
  TBranch *bEventList_EEstrpen02;
  TBranch *bEventList_EEstrpen03;
  TBranch *bEventList_EEstrpen0;
  TBranch *bEventList_EEstrpen1;
  TBranch *bEventList_EEstrpen2;
  TBranch *bEventList_EEstrpen3;
  TBranch *bEventList_EEstrpen4;
  TBranch *bEventList_EEstrpen5;
  TBranch *bEventList_EEstrpen6;
  TBranch *bEventList_EEstrpen7;
  TBranch *bEventList_EEstrpen8;
  TBranch *bEventList_EEstrpen9;
  TBranch *bEventList_EEstrpen10;
  TBranch *bEventList_EEstrpen11;
  TBranch *bEventList_EEstrpen12;
  TBranch *bEventList_EEstrpen13;
  TBranch *bEventList_EEstrpen14;
  TBranch *bEventList_EEstrpen15;
  TBranch *bEventList_ETwrdidE;
  TBranch *bEventList_EPstripenp01;
  TBranch *bEventList_EPstripenp02;
  TBranch *bEventList_EPstripenp03;
  TBranch *bEventList_EPstripenp0;
  TBranch *bEventList_EPstripenp1;
  TBranch *bEventList_EPstripenp2;
  TBranch *bEventList_EPstripenp3;
  TBranch *bEventList_EPstripenp4;
  TBranch *bEventList_EPstripenp5;
  TBranch *bEventList_EPstripenp6;
  TBranch *bEventList_EPstripenp7;
  TBranch *bEventList_EPstripenp8;
  TBranch *bEventList_EPstripenp9;
  TBranch *bEventList_EPstripenp10;
  TBranch *bEventList_EPstripenp11;
  TBranch *bEventList_EPstripenp12;
  TBranch *bEventList_EPstripenp13;
  TBranch *bEventList_EPstripenp14;
  TBranch *bEventList_EPstripenp15;
  TBranch *bEventList_EclustEnnq1;
  TBranch *bEventList_EclustEnnq20;
  TBranch *bEventList_EclustEnnq19;
  TBranch *bEventList_EclustEnpq1;
  TBranch *bEventList_EclustEnpq20;
  TBranch *bEventList_EclustEnpq19;
  TBranch *bEventList_EclustEnpq21;
  TBranch *bEventList_PrimaryTrackArray_;
  TBranch *bPrimaryTrackArray_fUniqueID;
  TBranch *bPrimaryTrackArray_fBits;
  TBranch *bPrimaryTrackArray_nHitsFit;
  TBranch *bPrimaryTrackArray_nHitsPoss;
  TBranch *bPrimaryTrackArray_trackFlag;
  TBranch *bPrimaryTrackArray_pZ;
  TBranch *bPrimaryTrackArray_pX;
  TBranch *bPrimaryTrackArray_pY;
  TBranch *bPrimaryTrackArray_pT;
  TBranch *bPrimaryTrackArray_dEdx;
  TBranch *bPrimaryTrackArray_charge;
  TBranch *bPrimaryTrackArray_tofBeta;
  TBranch *bPrimaryTrackArray_eta;
  TBranch *bPrimaryTrackArray_phi;
  TBranch *bPrimaryTrackArray_nSigElectron;
  TBranch *bPrimaryTrackArray_nSigPion;
  TBranch *bPrimaryTrackArray_nSigKaon;
  TBranch *bPrimaryTrackArray_nSigProton;
  TBranch *bPrimaryTrackArray_dcag;
  TBranch *bPrimaryTrackArray_nHits;
  TBranch *bPrimaryTrackArray_dEdxHits;
  TBranch *bPrimaryTrackArray_firstZPoint;
  TBranch *bPrimaryTrackArray_lastZPoint;
  TBranch *bPrimaryTrackArray_tofSigElectron;
  TBranch *bPrimaryTrackArray_tofSigPion;
  TBranch *bPrimaryTrackArray_tofSigKaon;
  TBranch *bPrimaryTrackArray_tofSigProton;
  TBranch *bPrimaryTrackArray_timeOfflight;
  TBranch *bPrimaryTrackArray_pathLength;
  TBranch *bPrimaryTrackArray_trkIndex;
  TBranch *bEventList_TowerArray_;
  TBranch *bTowerArray_fUniqueID;
  TBranch *bTowerArray_fBits;
  TBranch *bTowerArray_TwrId;
  TBranch *bTowerArray_TwrEng;
  TBranch *bTowerArray_TwrEta;
  TBranch *bTowerArray_TwrPhi;
  TBranch *bTowerArray_TwrADC;
  TBranch *bTowerArray_TwrPed;
  TBranch *bTowerArray_TwrRMS;
  TBranch *bTowerArray_TwrMatchIdnex;
  TBranch *bTowerArray_NoOfmatchedTrk;
  TBranch *bTowerArray_TwrMatchP;
  TBranch *bTowerArray_TwrPx;
  TBranch *bTowerArray_TwrPy;
  TBranch *bTowerArray_TwrPz;
  TBranch *bTowerArray_fNAssocTracks;
  TBranch *bTowerArray_fMatchedTracksArray_;
  TBranch *bTowerArray_fMatchedTracksArray_P;

  // set input branches
  tInput -> SetMakeClass(1);
  tInput -> SetBranchAddress("fUniqueID", &fUniqueID, &bEventList_fUniqueID);
  tInput -> SetBranchAddress("fBits", &fBits, &bEventList_fBits);
  tInput -> SetBranchAddress("runNumber", &runNumber, &bEventList_runNumber);
  tInput -> SetBranchAddress("eventNumber", &eventNumber, &bEventList_eventNumber);
  tInput -> SetBranchAddress("trigID", &trigID, &bEventList_trigID);
  tInput -> SetBranchAddress("nGlobalTracks", &nGlobalTracks, &bEventList_nGlobalTracks);
  tInput -> SetBranchAddress("nPrimaryTracks", &nPrimaryTracks, &bEventList_nPrimaryTracks);
  tInput -> SetBranchAddress("refMult", &refMult, &bEventList_refMult);
  tInput -> SetBranchAddress("vpdVz", &vpdVz, &bEventList_vpdVz);
  tInput -> SetBranchAddress("xVertex", &xVertex, &bEventList_xVertex);
  tInput -> SetBranchAddress("yVertex", &yVertex, &bEventList_yVertex);
  tInput -> SetBranchAddress("zVertex", &zVertex, &bEventList_zVertex);
  tInput -> SetBranchAddress("bbcZVertex", &bbcZVertex, &bEventList_bbcZVertex);
  tInput -> SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate, &bEventList_zdcCoincidenceRate);
  tInput -> SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate, &bEventList_bbcCoincidenceRate);
  tInput -> SetBranchAddress("backgroundRate", &backgroundRate, &bEventList_backgroundRate);
  tInput -> SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate, &bEventList_bbcBlueBackgroundRate);
  tInput -> SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate, &bEventList_bbcYellowBackgroundRate);
  tInput -> SetBranchAddress("refMultPos", &refMultPos, &bEventList_refMultPos);
  tInput -> SetBranchAddress("refMultNeg", &refMultNeg, &bEventList_refMultNeg);
  tInput -> SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity, &bEventList_bTOFTrayMultiplicity);
  tInput -> SetBranchAddress("nVerticies", &nVerticies, &bEventList_nVerticies);
  tInput -> SetBranchAddress("MagF", &MagF, &bEventList_MagF);
  tInput -> SetBranchAddress("VrtxRank", &VrtxRank, &bEventList_VrtxRank);
  tInput -> SetBranchAddress("FlagEvent_TrgTrkMisMtch", &FlagEvent_TrgTrkMisMtch, &bEventList_FlagEvent_TrgTrkMisMtch);
  tInput -> SetBranchAddress("Etsp", &Etsp, &bEventList_Etsp);
  tInput -> SetBranchAddress("ETwrdidT", &ETwrdidT, &bEventList_ETwrdidT);
  tInput -> SetBranchAddress("ETwradc11", &ETwradc11, &bEventList_ETwradc11);
  tInput -> SetBranchAddress("ETwreneT0", &ETwreneT0, &bEventList_ETwreneT0);
  tInput -> SetBranchAddress("ETwreT", &ETwreT, &bEventList_ETwreT);
  tInput -> SetBranchAddress("ETwrENET0", &ETwrENET0, &bEventList_ETwrENET0);
  tInput -> SetBranchAddress("ETwrphT", &ETwrphT, &bEventList_ETwrphT);
  tInput -> SetBranchAddress("ETwrPTower", &ETwrPTower, &bEventList_ETwrPTower);
  tInput -> SetBranchAddress("ETwrpidTower", &ETwrpidTower, &bEventList_ETwrpidTower);
  tInput -> SetBranchAddress("ETwrmoduleT", &ETwrmoduleT, &bEventList_ETwrmoduleT);
  tInput -> SetBranchAddress("EClustEneT0", &EClustEneT0, &bEventList_EClustEneT0);
  tInput -> SetBranchAddress("EClustetav1", &EClustetav1, &bEventList_EClustetav1);
  tInput -> SetBranchAddress("EClustphiv1", &EClustphiv1, &bEventList_EClustphiv1);
  tInput -> SetBranchAddress("EEstrpen01", &EEstrpen01, &bEventList_EEstrpen01);
  tInput -> SetBranchAddress("EEstrpen02", &EEstrpen02, &bEventList_EEstrpen02);
  tInput -> SetBranchAddress("EEstrpen03", &EEstrpen03, &bEventList_EEstrpen03);
  tInput -> SetBranchAddress("EEstrpen0", &EEstrpen0, &bEventList_EEstrpen0);
  tInput -> SetBranchAddress("EEstrpen1", &EEstrpen1, &bEventList_EEstrpen1);
  tInput -> SetBranchAddress("EEstrpen2", &EEstrpen2, &bEventList_EEstrpen2);
  tInput -> SetBranchAddress("EEstrpen3", &EEstrpen3, &bEventList_EEstrpen3);
  tInput -> SetBranchAddress("EEstrpen4", &EEstrpen4, &bEventList_EEstrpen4);
  tInput -> SetBranchAddress("EEstrpen5", &EEstrpen5, &bEventList_EEstrpen5);
  tInput -> SetBranchAddress("EEstrpen6", &EEstrpen6, &bEventList_EEstrpen6);
  tInput -> SetBranchAddress("EEstrpen7", &EEstrpen7, &bEventList_EEstrpen7);
  tInput -> SetBranchAddress("EEstrpen8", &EEstrpen8, &bEventList_EEstrpen8);
  tInput -> SetBranchAddress("EEstrpen9", &EEstrpen9, &bEventList_EEstrpen9);
  tInput -> SetBranchAddress("EEstrpen10", &EEstrpen10, &bEventList_EEstrpen10);
  tInput -> SetBranchAddress("EEstrpen11", &EEstrpen11, &bEventList_EEstrpen11);
  tInput -> SetBranchAddress("EEstrpen12", &EEstrpen12, &bEventList_EEstrpen12);
  tInput -> SetBranchAddress("EEstrpen13", &EEstrpen13, &bEventList_EEstrpen13);
  tInput -> SetBranchAddress("EEstrpen14", &EEstrpen14, &bEventList_EEstrpen14);
  tInput -> SetBranchAddress("EEstrpen15", &EEstrpen15, &bEventList_EEstrpen15);
  tInput -> SetBranchAddress("ETwrdidE", &ETwrdidE, &bEventList_ETwrdidE);
  tInput -> SetBranchAddress("EPstripenp01", &EPstripenp01, &bEventList_EPstripenp01);
  tInput -> SetBranchAddress("EPstripenp02", &EPstripenp02, &bEventList_EPstripenp02);
  tInput -> SetBranchAddress("EPstripenp03", &EPstripenp03, &bEventList_EPstripenp03);
  tInput -> SetBranchAddress("EPstripenp0", &EPstripenp0, &bEventList_EPstripenp0);
  tInput -> SetBranchAddress("EPstripenp1", &EPstripenp1, &bEventList_EPstripenp1);
  tInput -> SetBranchAddress("EPstripenp2", &EPstripenp2, &bEventList_EPstripenp2);
  tInput -> SetBranchAddress("EPstripenp3", &EPstripenp3, &bEventList_EPstripenp3);
  tInput -> SetBranchAddress("EPstripenp4", &EPstripenp4, &bEventList_EPstripenp4);
  tInput -> SetBranchAddress("EPstripenp5", &EPstripenp5, &bEventList_EPstripenp5);
  tInput -> SetBranchAddress("EPstripenp6", &EPstripenp6, &bEventList_EPstripenp6);
  tInput -> SetBranchAddress("EPstripenp7", &EPstripenp7, &bEventList_EPstripenp7);
  tInput -> SetBranchAddress("EPstripenp8", &EPstripenp8, &bEventList_EPstripenp8);
  tInput -> SetBranchAddress("EPstripenp9", &EPstripenp9, &bEventList_EPstripenp9);
  tInput -> SetBranchAddress("EPstripenp10", &EPstripenp10, &bEventList_EPstripenp10);
  tInput -> SetBranchAddress("EPstripenp11", &EPstripenp11, &bEventList_EPstripenp11);
  tInput -> SetBranchAddress("EPstripenp12", &EPstripenp12, &bEventList_EPstripenp12);
  tInput -> SetBranchAddress("EPstripenp13", &EPstripenp13, &bEventList_EPstripenp13);
  tInput -> SetBranchAddress("EPstripenp14", &EPstripenp14, &bEventList_EPstripenp14);
  tInput -> SetBranchAddress("EPstripenp15", &EPstripenp15, &bEventList_EPstripenp15);
  tInput -> SetBranchAddress("EclustEnnq1", &EclustEnnq1, &bEventList_EclustEnnq1);
  tInput -> SetBranchAddress("EclustEnnq20", &EclustEnnq20, &bEventList_EclustEnnq20);
  tInput -> SetBranchAddress("EclustEnnq19", &EclustEnnq19, &bEventList_EclustEnnq19);
  tInput -> SetBranchAddress("EclustEnpq1", &EclustEnpq1, &bEventList_EclustEnpq1);
  tInput -> SetBranchAddress("EclustEnpq20", &EclustEnpq20, &bEventList_EclustEnpq20);
  tInput -> SetBranchAddress("EclustEnpq19", &EclustEnpq19, &bEventList_EclustEnpq19);
  tInput -> SetBranchAddress("EclustEnpq21", &EclustEnpq21, &bEventList_EclustEnpq21);
  tInput -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_, &bEventList_PrimaryTrackArray_);
  tInput -> SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID, &bPrimaryTrackArray_fUniqueID);
  tInput -> SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits, &bPrimaryTrackArray_fBits);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit, &bPrimaryTrackArray_nHitsFit);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss, &bPrimaryTrackArray_nHitsPoss);
  tInput -> SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag, &bPrimaryTrackArray_trackFlag);
  tInput -> SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ, &bPrimaryTrackArray_pZ);
  tInput -> SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX, &bPrimaryTrackArray_pX);
  tInput -> SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY, &bPrimaryTrackArray_pY);
  tInput -> SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT, &bPrimaryTrackArray_pT);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx, &bPrimaryTrackArray_dEdx);
  tInput -> SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge, &bPrimaryTrackArray_charge);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta, &bPrimaryTrackArray_tofBeta);
  tInput -> SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta, &bPrimaryTrackArray_eta);
  tInput -> SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi, &bPrimaryTrackArray_phi);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron, &bPrimaryTrackArray_nSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion, &bPrimaryTrackArray_nSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon, &bPrimaryTrackArray_nSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton, &bPrimaryTrackArray_nSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag, &bPrimaryTrackArray_dcag);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits, &bPrimaryTrackArray_nHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits, &bPrimaryTrackArray_dEdxHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint, &bPrimaryTrackArray_firstZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint, &bPrimaryTrackArray_lastZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron, &bPrimaryTrackArray_tofSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion, &bPrimaryTrackArray_tofSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon, &bPrimaryTrackArray_tofSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton, &bPrimaryTrackArray_tofSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight, &bPrimaryTrackArray_timeOfflight);
  tInput -> SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength, &bPrimaryTrackArray_pathLength);
  tInput -> SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex, &bPrimaryTrackArray_trkIndex);
  tInput -> SetBranchAddress("TowerArray", &TowerArray_, &bEventList_TowerArray_);
  tInput -> SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID, &bTowerArray_fUniqueID);
  tInput -> SetBranchAddress("TowerArray.fBits", TowerArray_fBits, &bTowerArray_fBits);
  tInput -> SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId, &bTowerArray_TwrId);
  tInput -> SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng, &bTowerArray_TwrEng);
  tInput -> SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta, &bTowerArray_TwrEta);
  tInput -> SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi, &bTowerArray_TwrPhi);
  tInput -> SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC, &bTowerArray_TwrADC);
  tInput -> SetBranchAddress("TowerArray.TwrPed", TowerArray_TwrPed, &bTowerArray_TwrPed);
  tInput -> SetBranchAddress("TowerArray.TwrRMS", TowerArray_TwrRMS, &bTowerArray_TwrRMS);
  tInput -> SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex, &bTowerArray_TwrMatchIdnex);
  tInput -> SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk, &bTowerArray_NoOfmatchedTrk);
  tInput -> SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP, &bTowerArray_TwrMatchP);
  tInput -> SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx, &bTowerArray_TwrPx);
  tInput -> SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy, &bTowerArray_TwrPy);
  tInput -> SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz, &bTowerArray_TwrPz);
  tInput -> SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks, &bTowerArray_fNAssocTracks);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_, &bTowerArray_fMatchedTracksArray_);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P, &bTowerArray_fMatchedTracksArray_P);
  cout << "    Set branches." << endl;


  // define bad run and hot tower lists
  const UInt_t badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};
  const UInt_t hotTwrList[NHotTwr] = {1, 35, 141, 187, 224, 341, 424, 594, 814, 899, 900, 1046, 1128, 1132, 1244, 1382, 1388, 1405, 1588, 1766, 1773, 2066, 2160, 2253, 2281, 2284, 2301, 2303, 2306, 2590, 3007, 3495, 3840, 4043, 4047, 4053, 4057, 4121, 4442, 4569, 4617};
  cout << "    Bad run and hot tower lists defined:\n"
       << "      " << NBadRuns << " bad runs, " << NHotTwr << " hot towers."
       << endl;


  // create histograms
  TH1D *hEvtNum;
  TH1D *hEvtVz[2];
  TH1D *hEvtVr[2];
  TH1D *hTrgEt[2];
  TH1D *hTrgEta;
  TH1D *hTrgPhi;
  TH1D *hTrgTsp[NTrgTsp + 1];
  TH1D *hTrgEtBin[NTrgBins][NTrgTsp];
  TH1D *hTrkNfit[2];
  TH1D *hTrkRfit[2];
  TH1D *hTrkDca[2];
  TH1D *hTrkEta[2];
  TH1D *hTrkPt[2];
  TH1D *hTrkPtBin[NTrgTsp];
  TH1D *hTrkDfBin[NTrgTsp];
  TH2D *hEvtVyVsVx;
  TH2D *hEvtVxVsVz;
  TH2D *hEvtVyVsVz;
  TH2D *hTrgEtVsEta;
  TH2D *hTrgEtVsPhi;
  TH2D *hTrgEtaVsPhi;
  TH2D *hTrkNfitVsPt;
  TH2D *hTrkDcaVsPt;
  TH2D *hTrkPtVsDf[NTrgTsp];
  TH2D *hTrkDfVsEta[NTrgTsp];

  const UInt_t  nNum(NTrgCuts);
  const UInt_t  nVz(800);
  const UInt_t  nVr(400);
  const UInt_t  nPrim(100);
  const UInt_t  nEt(50);
  const UInt_t  nEta(40);
  const UInt_t  nPhi(30);
  const UInt_t  nTsp(100);
  const UInt_t  nFit(50);
  const UInt_t  nRat(100);
  const UInt_t  nDca(50);
  const UInt_t  nPt(275);
  const UInt_t  nPt2(110);
  const UInt_t  nPt3(55);
  const UInt_t  nDf(30);
  const Float_t num[2]  = {0., (Float_t) NTrgCuts};
  const Float_t vz[2]   = {-200., 200.};
  const Float_t vr[2]   = {-2., 2.};
  const Float_t prim[2] = {0., 100.};
  const Float_t et[2]   = {0., 50.};
  const Float_t eta[2]  = {-2., 2.};
  const Float_t phi[2]  = {-3.15, 3.15};
  const Float_t tsp[2]  = {0., 1.};
  const Float_t fit[2]  = {0., 50.};
  const Float_t rat[2]  = {0., 1.};
  const Float_t dca[2]  = {0., 5.};
  const Float_t pt[2]   = {-5., 50.};
  const Float_t df[2]   = {-1.6, 4.7};

  hEvtNum         = new TH1D("hEvtNum", "", nNum, num[0], num[1]);
  hEvtVz[0]       = new TH1D("hEvtVzAll", "", nVz, vz[0], vz[1]);
  hEvtVz[1]       = new TH1D("hEvtVzCut", "", nVz, vz[0], vz[1]);
  hEvtVr[0]       = new TH1D("hEvtVrAll", "", nVr, vr[0], vr[1]);
  hEvtVr[1]       = new TH1D("hEvtVrCut", "", nVr, vr[0], vr[1]);
  hEvtPrim        = new TH1D("hEvtPrim", "", nPrim, prim[0], prim[1]);
  hTrgEt[0]       = new TH1D("hTrgEtAll", "", nEt, et[0], et[1]);
  hTrgEt[1]       = new TH1D("hTrgEtCut", "", nEt, et[0], et[1]);
  hTrgEta         = new TH1D("hTrgEta", "", nEta, eta[0], eta[1]);
  hTrgPhi         = new TH1D("hTrgPhi", "", nPhi, phi[0], phi[1]);
  hTrgTsp[0]      = new TH1D("hTrgTspPi0", "", nTsp, tsp[0], tsp[1]);
  hTrgTsp[1]      = new TH1D("hTrgTspGam", "", nTsp, tsp[0], tsp[1]);
  hTrgTsp[2]      = new TH1D("hTrgTspAll", "", nTsp, tsp[0], tsp[1]); 
  hTrgEtBin[0][0] = new TH1D("hTrgEtPi920", "", nEt, et[0], et[1]);
  hTrgEtBin[0][1] = new TH1D("hTrgEtGa920", "", nEt, et[0], et[1]);
  hTrgEtBin[1][0] = new TH1D("hTrgEtPi911", "", nEt, et[0], et[1]);
  hTrgEtBin[1][1] = new TH1D("hTrgEtGa911", "", nEt, et[0], et[1]);
  hTrgEtBin[2][0] = new TH1D("hTrgEtPi1115", "", nEt, et[0], et[1]);
  hTrgEtBin[2][1] = new TH1D("hTrgEtGa1115", "", nEt, et[0], et[1]);
  hTrgEtBin[3][0] = new TH1D("hTrgEtPi1520", "", nEt, et[0], et[1]);
  hTrgEtBin[3][1] = new TH1D("hTrgEtGa1520", "", nEt, et[0], et[1]);
  hTrkNfit[0]     = new TH1D("hTrkNfitAll", "", nFit, fit[0], fit[1]);
  hTrkNfit[1]     = new TH1D("hTrkNfitCut", "", nFit, fit[0], fit[1]);
  hTrkRfit[0]     = new TH1D("hTrkRfitAll", "", nRat, rat[0], rat[1]);
  hTrkRfit[1]     = new TH1D("hTrkRfitCut", "", nRat, rat[0], rat[1]);
  hTrkDca[0]      = new TH1D("hTrkDcaAll", "", nDca, dca[0], dca[1]);
  hTrkDca[1]      = new TH1D("hTrkDcaCut", "", nDca, dca[0], dca[1]);
  hTrkEta[0]      = new TH1D("hTrkEtaAll", "", nEta, eta[0], eta[1]);
  hTrkEta[1]      = new TH1D("hTrkEtaCut", "", nEta, eta[0], eta[1]);
  hTrkPt[0]       = new TH1D("hTrkPtAll", "", nPt, pt[0], pt[1]);
  hTrkPt[1]       = new TH1D("hTrkPtCut", "", nPt, pt[0], pt[1]);
  hTrkPtBin[0]    = new TH1D("hTrkPtPi0", "", nPt2, pt[0], pt[1]);
  hTrkPtBin[1]    = new TH1D("hTrkPtGam", "", nPt2, pt[0], pt[1]);
  hTrkDfBin[0]    = new TH1D("hTrkDfPi0", "", nDf, df[0], df[1]);
  hTrkDfBin[1]    = new TH1D("hTrkDfGam", "", nDf, df[0], df[1]);
  hEvtVyVsVx      = new TH2D("hEvtVyVsVx", "", nVr, vr[0], vr[1], nVr, vr[0], vr[1]);
  hEvtVxVsVz      = new TH2D("hEvtVxVsVz", "", nVz, vz[0], vz[1], nVr, vr[0], vr[1]);
  hEvtVyVsVz      = new TH2D("hEvtVyVsVz", "", nVz, vz[0], vz[1], nVr, vr[0], vr[1]);
  hTrgEtVsEta     = new TH2D("hTrgEtVsEta", "", nEta, eta[0], eta[1], nEt, et[0], et[1]);
  hTrgEtVsPhi     = new TH2D("hTrgEtVsPhi", "", nPhi, phi[0], phi[1], nEt, et[0], et[1]);
  hTrgEtaVsPhi    = new TH2D("hTrgEtaVsPhi", "", nPhi, phi[0], phi[1], nEta, eta[0], eta[1]);
  hTrkNfitVsPt    = new TH2D("hTrkNfitVsPt", "", nPt3, pt[0], pt[1], nFit, fit[0], fit[1]);
  hTrkDcaVsPt     = new TH2D("hTrkDcaVsPt", "", nPt3, pt[0], pt[1], nDca, dca[0], dca[1]);
  hTrkPtVsDf[0]   = new TH2D("hTrkPtVsDfPi0", "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
  hTrkPtVsDf[1]   = new TH2D("hTrkPtVsDfGam", "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
  hTrkDfVsEta[0]  = new TH2D("hTrkDfVsEtaPi0", "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
  hTrkDfVsEta[1]  = new TH2D("hTrkDfVsEtaGam", "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
  hEvtVz[0]       -> Sumw2();
  hEvtVz[1]       -> Sumw2();
  hEvtVr[0]       -> Sumw2();
  hEvtVr[1]       -> Sumw2();
  hEvtPrim        -> Sumw2();
  hTrgEt[0]       -> Sumw2();
  hTrgEt[1]       -> Sumw2();
  hTrgEta         -> Sumw2();
  hTrgPhi         -> Sumw2();
  hTrgTsp[0]      -> Sumw2();
  hTrgTsp[1]      -> Sumw2();
  hTrgTsp[2]      -> Sumw2();
  hTrgEtBin[0][0] -> Sumw2();
  hTrgEtBin[0][1] -> Sumw2();
  hTrgEtBin[1][0] -> Sumw2();
  hTrgEtBin[1][1] -> Sumw2();
  hTrgEtBin[2][0] -> Sumw2();
  hTrgEtBin[2][1] -> Sumw2();
  hTrgEtBin[3][0] -> Sumw2();
  hTrgEtBin[3][1] -> Sumw2();
  hTrkNfit[0]     -> Sumw2();
  hTrkNfit[1]     -> Sumw2();
  hTrkRfit[0]     -> Sumw2();
  hTrkRfit[1]     -> Sumw2();
  hTrkDca[0]      -> Sumw2();
  hTrkDca[1]      -> Sumw2();
  hTrkEta[0]      -> Sumw2();
  hTrkEta[1]      -> Sumw2();
  hTrkPt[0]       -> Sumw2();
  hTrkPt[1]       -> Sumw2();
  hTrkPtBin[0]    -> Sumw2();
  hTrkPtBin[1]    -> Sumw2();
  hTrkDfBin[0]    -> Sumw2();
  hTrkDfBin[1]    -> Sumw2();
  hEvtVyVsVx      -> Sumw2();
  hEvtVxVsVz      -> Sumw2();
  hEvtVyVsVz      -> Sumw2();
  hTrgEtVsEta     -> Sumw2();
  hTrgEtVsPhi     -> Sumw2();
  hTrgEtaVsPhi    -> Sumw2();
  hTrkNfitVsPt    -> Sumw2();
  hTrkDcaVsPt     -> Sumw2();
  hTrkPtVsDf[0]   -> Sumw2();
  hTrkPtVsDf[1]   -> Sumw2();
  hTrkDfVsEta[0]  -> Sumw2();
  hTrkDfVsEta[1]  -> Sumw2();


  // no. of evts and triggers
  UInt_t nTrgCut[NTrgCuts];
  UInt_t nTrgPi0[NTrgBins];
  UInt_t nTrgGam[NTrgBins];
  for (UInt_t iCut = 0; iCut < NTrgCuts; iCut++) {
    nTrgCut[iCut] = 0;
  }
  for (UInt_t iBin = 0; iBin < NTrgBins; iBin++) {
    nTrgPi0[iBin] = 0;
    nTrgGam[iBin] = 0;
  }

  const UInt_t nEvts = tInput -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " events to process." << endl;


  // event loop
  UInt_t bytes(0);
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


    // event info
    const UInt_t   run   = runNumber;
    const Long64_t nTrks = nPrimaryTracks;
    const Long64_t nTwrs = TowerArray_;
    const Double_t xVtx  = xVertex;
    const Double_t yVtx  = yVertex;
    const Double_t zVtx  = zVertex;
    const Double_t rVtx  = TMath::Sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // filter out bad runs
    Bool_t isGoodRun = true;
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (run == badRunList[iRun]) {
        isGoodRun = false;
        break;
      }
    }
    nTrgCut[0]++;
    if (isGoodRun)
      nTrgCut[1]++;
    else
      continue;


    // fill event histograms
    hEvtVz[0]  -> Fill(zVtx);
    hEvtVr[0]  -> Fill(rVtx);
    hEvtPrim   -> Fill(nTrks);
    hEvtVyVsVx -> Fill(xVtx, yVtx);
    hEvtVxVsVz -> Fill(zVtx, xVtx);
    hEvtVyVsVz -> Fill(zVtx, yVtx);

    // vertex cuts
    const Bool_t isInRcut = (TMath::Abs(rVtx) < rVtxMax);
    const Bool_t isInZcut = (TMath::Abs(zVtx) < zVtxMax);
    if (isInRcut && isInZcut)
      nTrgCut[2]++;
    else {
      hEvtVz[1] -> Fill(zVtx);
      hEvtVr[1] -> Fill(rVtx);
      continue;
    }


    // trigger info
    const Int_t    adc    = ETwradc11;
    const UInt_t   idTrg  = ETwrdidT;
    const Double_t tspTrg = Etsp;
    const Double_t eH4    = EEstrpen4;
    const Double_t eF4    = EPstripenp4;
    const Double_t pProj  = ETwrPTower;
    const Double_t hDet   = ETwreT;
    const Double_t hPhys  = EClustetav1;
    const Double_t fTrg   = EClustphiv1;
    const Double_t eTrg   = EClustEneT0;
    const Double_t tTrg   = 2. * TMath::ATan(TMath::Exp(-1. * hPhys));
    const Double_t eTtrg  = eTrg * TMath::Sin(tTrg);

    // filter out bad towers
    Bool_t isGoodTwr = true;
    for (UInt_t iTwr = 0; iTwr < NHotTwr; iTwr++) {
      if (idTrg == hotTwrList[iTwr]) {
        isGoodTwr = false;
        break;
      }
    }
    if (isGoodTwr)
      nTrgCut[3]++;
    else
      continue;

    // trigger cuts
    const Bool_t isInAdcCut    = (adc <= adcMax);
    const Bool_t isInStrCut    = ((eH4 >= eStrMin) && (eF4 >= eStrMin));
    const Bool_t isInProjCut   = (pProj < pProjMax);
    const Bool_t isInEtaTrgCut = (TMath::Abs(hDet) < hTrgMax);
    const Bool_t isInEtCut     = ((eTtrg >= eTtrgMin) && (eTtrg < eTtrgMax));
    const Bool_t isInPi0cut    = ((tspTrg > tspPi0[0]) && (tspTrg < tspPi0[1]));
    const Bool_t isInGamCut    = ((tspTrg > tspGam[0]) && (tspTrg < tspGam[1]));
    const Bool_t isInTspCut    = (isInPi0cut || isInGamCut);

    // central strip cut
    if (isInAdcCut && isInStrCut)
      nTrgCut[4]++;
    else
      continue;

    // charged momentum cut
    if (isInProjCut)
      nTrgCut[5]++;
    else
      continue;

    // eta cut
    if (isInEtaTrgCut)
      nTrgCut[6]++;
    else
      continue;

    // eT and tsp cuts
    if (isInEtCut)
      nTrgCut[7]++;
    if (isInEtCut && isInTspCut)
      nTrgCut[8]++;


    // figure out eT bin
    UInt_t iEtBin(5);
    for (UInt_t iBin = 1; iBin < NTrgBins; iBin++) {
      const Bool_t isInBin = ((eTtrg >= eTbinMin[iBin]) && (eTtrg < eTbinMax[iBin]));
      if (isInBin) {
        iEtBin = iBin;
        break;
      }
    }

    // determine species
    if (isInEtCut && isInPi0cut) {
      hTrgTsp[0]           -> Fill(tspTrg);
      hTrgEtBin[0][0]      -> Fill(eTtrg);
      hTrgEtBin[iEtBin][0] -> Fill(eTtrg);
      nTrgPi0[0]++;
      nTrgPi0[iEtBin]++;
    }
    if (isInEtCut && isInGamCut) {
      hTrgTsp[1]           -> Fill(tspTrg);
      hTrgEtBin[0][1]      -> Fill(eTtrg);
      hTrgEtBin[iEtBin][1] -> Fill(eTtrg);
      nTrgGam[0]++;
      nTrgGam[iEtBin]++;
    }

    // fill other histograms
    hTrgEt[0] -> Fill(eTtrg);
    if (isInTspCut) {
      hTrgEta      -> Fill(hPhys);
      hTrgPhi      -> Fill(fTrg);
      hTrgEtVsEta  -> Fill(hPhys, eTtrg);
      hTrgEtVsPhi  -> Fill(fTrg, eTtrg);
      hTrgEtaVsPhi -> Fill(fTrg, hPhys);
    }

    if (isInEtCut)
      hTrgTsp[2] -> Fill(tspTrg);
    else
      hTrgEt[1] -> Fill(eTtrg);

    // leave only pi0 and gamma candidates
    if (!isInEtCut || !isInTspCut) continue;


    // track loop
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const UInt_t   nFitTrk  = PrimaryTrackArray_nHitsFit[iTrk];
      const UInt_t   nPossTrk = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFitTrk  = (Double_t) nFitTrk / (Double_t) nPossTrk;
      const Double_t dcaTrk   = PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk     = PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk     = PrimaryTrackArray_phi[iTrk];
      const Double_t pTtrk    = PrimaryTrackArray_pT[iTrk];
      const Double_t zTtrk    = pTtrk / eTtrg;

      Double_t dFtrk = fTrk - fTrg;
      if (dFtrk < (-1. * TMath::PiOver2())) dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))  dFtrk -= TMath::TwoPi();


      // track cuts
      const Bool_t isInFitCut    = (nFitTrk >= nFitMin);
      const Bool_t isInRatioCut  = (rFitTrk >= rFitMin);
      const Bool_t isInDcaCut    = (dcaTrk < dcaMax);
      const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtCut     = ((pTtrk > pTtrkMin) && (pTtrk < pTtrkMax));

      // nFit cut
      hTrkNfit[0]  -> Fill(nFitTrk);
      hTrkNfitVsPt -> Fill(pTtrk, nFitTrk);
      if (!isInFitCut)
        hTrkNfit[1] -> Fill(nFitTrk);

      // nFit / nPoss cut
      hTrkRfit[0] -> Fill(rFitTrk);
      if (!isInRatioCut)
        hTrkRfit[1] -> Fill(rFitTrk);

      // dca cut
      hTrkDca[0]  -> Fill(dcaTrk);
      hTrkDcaVsPt -> Fill(pTtrk, dcaTrk);
      if (!isInDcaCut)
        hTrkDca[1] -> Fill(dcaTrk);

      // eta cut
      hTrkEta[0] -> Fill(hTrk);
      if (!isInEtaTrkCut)
        hTrkEta[1] -> Fill(hTrk);

      // pT cut
      hTrkPt[0] -> Fill(pTtrk);
      if (!isInPtCut)
        hTrkPt[1] -> Fill(pTtrk);

      // apply cuts
      if (!isInFitCut || !isInRatioCut || !isInDcaCut || !isInEtaTrkCut) continue;


      // fill histograms
      if (isInPi0cut) {
        hTrkPtBin[0]  -> Fill(pTtrk);
        hTrkPtVsDf[0] -> Fill(dFtrk, pTtrk);
        if (isInPtCut) {
          hTrkDfBin[0]   -> Fill(dFtrk);
          hTrkDfVsEta[0] -> Fill(hTrk, dFtrk);
        }
      }
      if (isInGamCut) {
        hTrkPtBin[1]  -> Fill(pTtrk);
        hTrkPtVsDf[1] -> Fill(dFtrk, pTtrk);
        if (isInPtCut) {
          hTrkDfBin[1]   -> Fill(dFtrk);
          hTrkDfVsEta[1] -> Fill(hTrk, dFtrk);
        }
      }

    }  // end track loop

  }  // end event loop

  // report no. of triggers
  cout << "    Event loop finished:" << endl;
  for (UInt_t iCut = 0; iCut < NTrgCuts; iCut++) {
    hEvtNum -> SetBinContent(iCut + 1, nTrgCut[iCut]);
    cout << "      " << nTrgCut[iCut] << " evts. left after cut " << iCut << endl;
  }

  // report no. of pi0's and gamma's
  cout << "    Number of triggers:" << endl;
  for (UInt_t iBin = 0; iBin < NTrgBins; iBin++) {
    cout << "      eTtrg = (" << eTbinMin[iBin] << ", " << eTbinMax[iBin]
         << ") GeV: nPi0 = " << nTrgPi0[iBin] << ", nGam = " << nTrgGam[iBin]
         << endl;
  }

  // sum trigger histograms
  TH1D *hTrgEtSum = (TH1D*) hTrgEtBin[0][0] -> Clone();
  hTrgEtSum -> SetName("hTrgEtSum");
  hTrgEtSum -> Add(hTrgEtBin[0][1]);

  // create profiles
  TProfile *pTrkNfitVsPt = hTrkNfitVsPt -> ProfileX("pTrkNfitVsPt", 1, -1, "S");
  TProfile *pTrkDcaVsPt  = hTrkDcaVsPt  -> ProfileX("pTrkDcaVsPt", 1, -1, "S");


  // normalize relevant histograms
  const Double_t pTbin2 = (pt[1] - pt[0]) / (Double_t) nPt2;
  const Double_t pTbin3 = (pt[1] - pt[0]) / (Double_t) nPt3;
  const Double_t dFbin  = (df[1] - df[0]) / (Double_t) nDf;
  const Double_t etaBin = (eta[1] - eta[0]) / (Double_t) nEta;
  for (UInt_t iTsp = 0; iTsp < NTrgTsp; iTsp++) {
    // no. of triggers
    UInt_t nTrg(1);
    if (iTsp == 0) nTrg = nTrgPi0[0];
    if (iTsp == 1) nTrg = nTrgGam[0];

    // norms
    const Double_t pTnorm    = 1. / (pTbin2 * nTrg);
    const Double_t dFnorm    = 1. / (dFbin * nTrg);
    const Double_t pTdFnorm  = 1. / (pTbin3 * dFbin * nTrg);
    const Double_t etaDfNorm = 1. / (etaBin * dFbin * nTrg);

    hTrkPtBin[iTsp]   -> Scale(pTnorm);
    hTrkDfBin[iTsp]   -> Scale(dFnorm);
    hTrkPtVsDf[iTsp]  -> Scale(pTdFnorm);
    hTrkDfVsEta[iTsp] -> Scale(etaDfNorm);
  }
  cout << "    Normalized certain histograms." << endl;


  // set styles
  const UInt_t  fColAll(1);
  const UInt_t  fColFil(0);
  const UInt_t  fColCut(810);
  const UInt_t  fLinAll(1);
  const UInt_t  fLinCut(1);
  const UInt_t  fFilAll(0);
  const UInt_t  fFilCut(3472);
  const UInt_t  fMarAll(8);
  const UInt_t  fMarCut(8);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const UInt_t  fColTrg[NTrgTsp]  = {860, 810};
  const UInt_t  fColPi0[NTrgBins] = {1, 858, 848, 818};
  const UInt_t  fColGam[NTrgBins] = {1, 808, 898, 888};
  const Float_t fBar(0.6);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.3);
  const Float_t fOffZ(1.);
  const Float_t fOffL(0.007);
  const Float_t fOffB(0.2);
  const Float_t xVz2dPlot[2] = {-55., 55.};
  const Float_t xEtPlot[2]   = {7., 23.};
  const Float_t yEtPlot[2]   = {7., 37777.};
  const TString sEvt("events");
  const TString sCount("counts");
  const TString sPrim("N^{primary}");
  const TString sTrgEt("E_{T}^{trg} [GeV]");
  const TString sTrgEta("#eta^{trg}");
  const TString sTrgPhi("#varphi^{trg}");
  const TString sTrgTsp("TSP");
  const TString sTrkNfit("N^{fit}");
  const TString sTrkRfit("N^{fit} / N^{poss}");
  const TString sTrkDca("DCA [cm]");
  const TString sTrkEta("#eta^{trk}");
  const TString sTrkPt("p_{T}^{trk} [GeV/c]");
  const TString sTrkDf("#Delta#varphi^{trk}");
  const TString sTrkPtY("(1/N^{trg}) dN^{trk}/dp_{T}^{trk} [GeV/c]^{-1}");
  const TString sTrkDfY("(1/N^{trg}) dN^{trk}/d#Delta#varphi^{trk}");
  const TString sVtx[NVtx]     = {"V_{x} [cm]", "V_{y} [cm]", "V_{z} [cm]", "V_{r} [cm]"};
  const TString sCut[NTrgCuts] = {"no cuts", "bad runs removed", "V_{z}, V_{r} cuts", "bad towers removed", "e_{#eta}, e_{#varphi} cuts", "P_{proj} cut", "#eta^{trg} cut", "E_{T}^{trg} cut", "TSP cut"};

  // event no's
  hEvtNum -> SetLineColor(883);
  hEvtNum -> SetLineStyle(1);
  hEvtNum -> SetFillColor(883);
  hEvtNum -> SetFillStyle(1001);
  hEvtNum -> SetMarkerColor(883);
  hEvtNum -> SetMarkerStyle(1);
  hEvtNum -> SetBarWidth(fBar);
  hEvtNum -> SetBarOffset(fOffB);
  hEvtNum -> SetTitle("");
  hEvtNum -> SetTitleFont(fTxt);
  hEvtNum -> GetXaxis() -> SetTitle("");
  hEvtNum -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtNum -> GetXaxis() -> SetTitleSize(fTit);
  hEvtNum -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtNum -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtNum -> GetXaxis() -> SetLabelSize(fLab);
  hEvtNum -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtNum -> GetXaxis() -> CenterTitle(fCnt);
  hEvtNum -> GetYaxis() -> SetTitle(sEvt.Data());
  hEvtNum -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtNum -> GetYaxis() -> SetTitleSize(fTit);
  hEvtNum -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtNum -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtNum -> GetYaxis() -> SetLabelSize(fLab);
  hEvtNum -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iBin = 1; iBin < (NTrgCuts + 1); iBin++) {
    hEvtNum -> GetXaxis() -> SetBinLabel(iBin, sCut[iBin - 1].Data());
  }

  // vertices
  hEvtVz[0]  -> SetLineColor(fColAll);
  hEvtVz[0]  -> SetLineStyle(fLinAll);
  hEvtVz[0]  -> SetFillColor(fColFil);
  hEvtVz[0]  -> SetFillStyle(fFilAll);
  hEvtVz[0]  -> SetMarkerColor(fColAll);
  hEvtVz[0]  -> SetMarkerStyle(fMarAll);
  hEvtVz[0]  -> SetTitle("");
  hEvtVz[0]  -> SetTitleFont(fTxt);
  hEvtVz[0]  -> GetXaxis() -> SetTitle(sVtx[2].Data());
  hEvtVz[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVz[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVz[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVz[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVz[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVz[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVz[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVz[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVz[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVz[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVz[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVz[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVz[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVz[1]  -> SetLineColor(fColCut);
  hEvtVz[1]  -> SetLineStyle(fLinCut);
  hEvtVz[1]  -> SetFillColor(fColCut);
  hEvtVz[1]  -> SetFillStyle(fFilCut);
  hEvtVz[1]  -> SetMarkerColor(fColCut);
  hEvtVz[1]  -> SetMarkerStyle(fMarCut);
  hEvtVz[1]  -> SetTitle("");
  hEvtVz[1]  -> SetTitleFont(fTxt);
  hEvtVz[1]  -> GetXaxis() -> SetTitle(sVtx[2].Data());
  hEvtVz[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVz[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVz[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVz[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVz[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVz[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVz[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVz[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVz[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVz[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVz[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVz[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVz[1]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVr[0]  -> SetLineColor(fColAll);
  hEvtVr[0]  -> SetLineStyle(fLinAll);
  hEvtVr[0]  -> SetFillColor(fColFil);
  hEvtVr[0]  -> SetFillStyle(fFilAll);
  hEvtVr[0]  -> SetMarkerColor(fColAll);
  hEvtVr[0]  -> SetMarkerStyle(fMarAll);
  hEvtVr[0]  -> SetTitle("");
  hEvtVr[0]  -> SetTitleFont(fTxt);
  hEvtVr[0]  -> GetXaxis() -> SetTitle(sVtx[3].Data());
  hEvtVr[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVr[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVr[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVr[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVr[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVr[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVr[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVr[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVr[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVr[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVr[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVr[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVr[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVr[1]  -> SetLineColor(fColCut);
  hEvtVr[1]  -> SetLineStyle(fLinCut);
  hEvtVr[1]  -> SetFillColor(fColCut);
  hEvtVr[1]  -> SetFillStyle(fFilCut);
  hEvtVr[1]  -> SetMarkerColor(fColCut);
  hEvtVr[1]  -> SetMarkerStyle(fMarCut);
  hEvtVr[1]  -> SetTitle("");
  hEvtVr[1]  -> SetTitleFont(fTxt);
  hEvtVr[1]  -> GetXaxis() -> SetTitle(sVtx[3].Data());
  hEvtVr[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVr[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVr[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVr[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVr[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVr[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVr[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVr[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVr[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVr[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVr[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVr[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVr[1]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVyVsVx -> SetTitle("");
  hEvtVyVsVx -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetXaxis() -> SetTitle(sVtx[0].Data());
  hEvtVyVsVx -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVyVsVx -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVyVsVx -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVx -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVyVsVx -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVyVsVx -> GetYaxis() -> SetTitle(sVtx[1].Data());
  hEvtVyVsVx -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVyVsVx -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVyVsVx -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVx -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVyVsVx -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVyVsVx -> GetZaxis() -> SetTitle("");
  hEvtVyVsVx -> GetZaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetZaxis() -> SetTitleSize(fTit);
  hEvtVyVsVx -> GetZaxis() -> SetTitleOffset(fOffZ);
  hEvtVyVsVx -> GetZaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVx -> GetZaxis() -> SetLabelSize(fLab);
  hEvtVyVsVx -> GetZaxis() -> CenterTitle(fCnt);
  hEvtVxVsVz -> SetTitle("");
  hEvtVxVsVz -> SetTitleFont(fTxt);
  hEvtVxVsVz -> GetXaxis() -> SetTitle(sVtx[2].Data());
  hEvtVxVsVz -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVxVsVz -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVxVsVz -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVxVsVz -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVxVsVz -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVxVsVz -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVxVsVz -> GetXaxis() -> SetRangeUser(xVz2dPlot[0], xVz2dPlot[1]);
  hEvtVxVsVz -> GetYaxis() -> SetTitle(sVtx[0].Data());
  hEvtVxVsVz -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVxVsVz -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVxVsVz -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVxVsVz -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVxVsVz -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVxVsVz -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVxVsVz -> GetZaxis() -> SetTitle("");
  hEvtVxVsVz -> GetZaxis() -> SetTitleFont(fTxt);
  hEvtVxVsVz -> GetZaxis() -> SetTitleSize(fTit);
  hEvtVxVsVz -> GetZaxis() -> SetTitleOffset(fOffZ);
  hEvtVxVsVz -> GetZaxis() -> SetLabelFont(fTxt);
  hEvtVxVsVz -> GetZaxis() -> SetLabelSize(fLab);
  hEvtVxVsVz -> GetZaxis() -> CenterTitle(fCnt);
  hEvtVyVsVz -> SetTitle("");
  hEvtVyVsVz -> SetTitleFont(fTxt);
  hEvtVyVsVz -> GetXaxis() -> SetTitle(sVtx[2].Data());
  hEvtVyVsVz -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVz -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVyVsVz -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVyVsVz -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVz -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVyVsVz -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVyVsVz -> GetXaxis() -> SetRangeUser(xVz2dPlot[0], xVz2dPlot[1]);
  hEvtVyVsVz -> GetYaxis() -> SetTitle(sVtx[1].Data());
  hEvtVyVsVz -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVz -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVyVsVz -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVyVsVz -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVz -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVyVsVz -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVyVsVz -> GetZaxis() -> SetTitle("");
  hEvtVyVsVz -> GetZaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVz -> GetZaxis() -> SetTitleSize(fTit);
  hEvtVyVsVz -> GetZaxis() -> SetTitleOffset(fOffZ);
  hEvtVyVsVz -> GetZaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVz -> GetZaxis() -> SetLabelSize(fLab);
  hEvtVyVsVz -> GetZaxis() -> CenterTitle(fCnt);

  // no. of tracks
  hEvtPrim -> SetLineColor(fColAll);
  hEvtPrim -> SetLineStyle(fLinAll);
  hEvtPrim -> SetFillColor(fColFil);
  hEvtPrim -> SetFillStyle(fFilAll);
  hEvtPrim -> SetMarkerColor(fColAll);
  hEvtPrim -> SetMarkerStyle(fMarAll);
  hEvtPrim -> SetTitle("");
  hEvtPrim -> SetTitleFont(fTxt);
  hEvtPrim -> GetXaxis() -> SetTitle(sPrim.Data());
  hEvtPrim -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtPrim -> GetXaxis() -> SetTitleSize(fTit);
  hEvtPrim -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtPrim -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtPrim -> GetXaxis() -> SetLabelSize(fLab);
  hEvtPrim -> GetXaxis() -> CenterTitle(fCnt);
  hEvtPrim -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtPrim -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtPrim -> GetYaxis() -> SetTitleSize(fTit);
  hEvtPrim -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtPrim -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtPrim -> GetYaxis() -> SetLabelSize(fLab);
  hEvtPrim -> GetYaxis() -> CenterTitle(fCnt);

  // eTtrg (no tsp cuts)
  hTrgEt[0]   -> SetLineColor(fColAll);
  hTrgEt[0]   -> SetLineStyle(fLinAll);
  hTrgEt[0]   -> SetFillColor(fColFil);
  hTrgEt[0]   -> SetFillStyle(fFilAll);
  hTrgEt[0]   -> SetMarkerColor(fColAll);
  hTrgEt[0]   -> SetMarkerStyle(fMarAll);
  hTrgEt[0]   -> SetTitle("");
  hTrgEt[0]   -> SetTitleFont(fTxt);
  hTrgEt[0]   -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEt[0]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEt[0]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEt[0]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEt[0]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEt[0]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEt[0]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEt[0]   -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEt[0]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEt[0]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEt[0]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEt[0]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEt[0]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEt[0]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEt[1]   -> SetLineColor(fColCut);
  hTrgEt[1]   -> SetLineStyle(fLinCut);
  hTrgEt[1]   -> SetFillColor(fColCut);
  hTrgEt[1]   -> SetFillStyle(fFilCut);
  hTrgEt[1]   -> SetMarkerColor(fColCut);
  hTrgEt[1]   -> SetMarkerStyle(fMarCut);
  hTrgEt[1]   -> SetTitle("");
  hTrgEt[1]   -> SetTitleFont(fTxt);
  hTrgEt[1]   -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEt[1]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEt[1]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEt[1]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEt[1]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEt[1]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEt[1]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEt[1]   -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEt[1]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEt[1]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEt[1]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEt[1]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEt[1]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEt[1]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtVsEta -> SetTitle("");
  hTrgEtVsEta -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetXaxis() -> SetTitle(sTrgEta.Data());
  hTrgEtVsEta -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtVsEta -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtVsEta -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtVsEta -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtVsEta -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtVsEta -> GetYaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtVsEta -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtVsEta -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtVsEta -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtVsEta -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtVsEta -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtVsEta -> GetZaxis() -> SetTitle("");
  hTrgEtVsEta -> GetZaxis() -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetZaxis() -> SetTitleSize(fTit);
  hTrgEtVsEta -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrgEtVsEta -> GetZaxis() -> SetLabelFont(fTxt);
  hTrgEtVsEta -> GetZaxis() -> SetLabelSize(fLab);
  hTrgEtVsEta -> GetZaxis() -> CenterTitle(fCnt);
  hTrgEtVsPhi -> SetTitle("");
  hTrgEtVsPhi -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetXaxis() -> SetTitle(sTrgPhi.Data());
  hTrgEtVsPhi -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtVsPhi -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtVsPhi -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtVsPhi -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtVsPhi -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtVsPhi -> GetYaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtVsPhi -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtVsPhi -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtVsPhi -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtVsPhi -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtVsPhi -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtVsPhi -> GetZaxis() -> SetTitle("");
  hTrgEtVsPhi -> GetZaxis() -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetZaxis() -> SetTitleSize(fTit);
  hTrgEtVsPhi -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrgEtVsPhi -> GetZaxis() -> SetLabelFont(fTxt);
  hTrgEtVsPhi -> GetZaxis() -> SetLabelSize(fLab);
  hTrgEtVsPhi -> GetZaxis() -> CenterTitle(fCnt);

  // eTtrg (tsp cuts)
  for (UInt_t iBin = 1; iBin < NTrgBins; iBin++) {
    hTrgEtBin[iBin][0] -> SetLineColor(fColPi0[iBin]);
    hTrgEtBin[iBin][0] -> SetLineStyle(fLinAll);
    hTrgEtBin[iBin][0] -> SetFillColor(fColPi0[iBin]);
    hTrgEtBin[iBin][0] -> SetFillStyle(fFilAll);
    hTrgEtBin[iBin][0] -> SetMarkerColor(fColPi0[iBin]);
    hTrgEtBin[iBin][0] -> SetMarkerStyle(fMarAll);
    hTrgEtBin[iBin][0] -> SetTitle("");
    hTrgEtBin[iBin][0] -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitle(sTrgEt.Data());
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][0] -> GetXaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitle(sCount.Data());
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][0] -> GetYaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetRangeUser(yEtPlot[0], yEtPlot[1]);
    hTrgEtBin[iBin][1] -> SetLineColor(fColGam[iBin]);
    hTrgEtBin[iBin][1] -> SetLineStyle(fLinAll);
    hTrgEtBin[iBin][1] -> SetFillColor(fColGam[iBin]);
    hTrgEtBin[iBin][1] -> SetFillStyle(fFilAll);
    hTrgEtBin[iBin][1] -> SetMarkerColor(fColGam[iBin]);
    hTrgEtBin[iBin][1] -> SetMarkerStyle(fMarAll);
    hTrgEtBin[iBin][1] -> SetTitle("");
    hTrgEtBin[iBin][1] -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitle(sTrgEt.Data());
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][1] -> GetXaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitle(sCount.Data());
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][1] -> GetYaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetRangeUser(yEtPlot[0], yEtPlot[1]);
  }
  hTrgEtBin[0][0] -> SetLineColor(fColTrg[0]);
  hTrgEtBin[0][0] -> SetLineStyle(fLinAll);
  hTrgEtBin[0][0] -> SetFillColor(fColTrg[0]);
  hTrgEtBin[0][0] -> SetFillStyle(fFilAll);
  hTrgEtBin[0][0] -> SetMarkerColor(fColTrg[0]);
  hTrgEtBin[0][0] -> SetMarkerStyle(fMarAll);
  hTrgEtBin[0][0] -> SetTitle("");
  hTrgEtBin[0][0] -> SetTitleFont(fTxt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtBin[0][0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtBin[0][0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtBin[0][1] -> SetLineColor(fColTrg[1]);
  hTrgEtBin[0][1] -> SetLineStyle(fLinAll);
  hTrgEtBin[0][1] -> SetFillColor(fColTrg[1]);
  hTrgEtBin[0][1] -> SetFillStyle(fFilAll);
  hTrgEtBin[0][1] -> SetMarkerColor(fColTrg[1]);
  hTrgEtBin[0][1] -> SetMarkerStyle(fMarAll);
  hTrgEtBin[0][1] -> SetTitle("");
  hTrgEtBin[0][1] -> SetTitleFont(fTxt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtBin[0][1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtBin[0][1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][1] -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtSum       -> SetLineColor(fColAll);
  hTrgEtSum       -> SetLineStyle(fLinAll);
  hTrgEtSum       -> SetFillColor(fColAll);
  hTrgEtSum       -> SetFillStyle(fFilAll);
  hTrgEtSum       -> SetMarkerColor(fColAll);
  hTrgEtSum       -> SetMarkerStyle(fMarAll);
  hTrgEtSum       -> SetTitle("");
  hTrgEtSum       -> SetTitleFont(fTxt);
  hTrgEtSum       -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtSum       -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtSum       -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtSum       -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtSum       -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtSum       -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtSum       -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtSum       -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
  hTrgEtSum       -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEtSum       -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtSum       -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtSum       -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtSum       -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtSum       -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtSum       -> GetYaxis() -> CenterTitle(fCnt);

  // fTrg vs hTrg
  hTrgEtaVsPhi -> SetTitle("");
  hTrgEtaVsPhi -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetXaxis() -> SetTitle(sTrgPhi.Data());
  hTrgEtaVsPhi -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtaVsPhi -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtaVsPhi -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtaVsPhi -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtaVsPhi -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtaVsPhi -> GetYaxis() -> SetTitle(sTrgEta.Data());
  hTrgEtaVsPhi -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtaVsPhi -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtaVsPhi -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtaVsPhi -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtaVsPhi -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtaVsPhi -> GetZaxis() -> SetTitle("");
  hTrgEtaVsPhi -> GetZaxis() -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetZaxis() -> SetTitleSize(fTit);
  hTrgEtaVsPhi -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrgEtaVsPhi -> GetZaxis() -> SetLabelFont(fTxt);
  hTrgEtaVsPhi -> GetZaxis() -> SetLabelSize(fLab);
  hTrgEtaVsPhi -> GetZaxis() -> CenterTitle(fCnt);

  // nFit and nFit/nPoss
  hTrkNfit[0]  -> SetLineColor(fColAll);
  hTrkNfit[0]  -> SetLineStyle(fLinAll);
  hTrkNfit[0]  -> SetFillColor(fColFil);
  hTrkNfit[0]  -> SetFillStyle(fFilAll);
  hTrkNfit[0]  -> SetMarkerColor(fColAll);
  hTrkNfit[0]  -> SetMarkerStyle(fMarAll);
  hTrkNfit[0]  -> SetTitle("");
  hTrkNfit[0]  -> SetTitleFont(fTxt);
  hTrkNfit[0]  -> GetXaxis() -> SetTitle(sTrkNfit.Data());
  hTrkNfit[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkNfit[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkNfit[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkNfit[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkNfit[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkNfit[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkNfit[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkNfit[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkNfit[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkNfit[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkNfit[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkNfit[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkNfit[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkNfit[1]  -> SetLineColor(fColCut);
  hTrkNfit[1]  -> SetLineStyle(fLinCut);
  hTrkNfit[1]  -> SetFillColor(fColCut);
  hTrkNfit[1]  -> SetFillStyle(fFilCut);
  hTrkNfit[1]  -> SetMarkerColor(fColCut);
  hTrkNfit[1]  -> SetMarkerStyle(fMarCut);
  hTrkNfit[1]  -> SetTitle("");
  hTrkNfit[1]  -> SetTitleFont(fTxt);
  hTrkNfit[1]  -> GetXaxis() -> SetTitle(sTrkNfit.Data());
  hTrkNfit[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkNfit[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkNfit[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkNfit[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkNfit[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkNfit[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkNfit[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkNfit[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkNfit[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkNfit[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkNfit[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkNfit[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkNfit[1]  -> GetYaxis() -> CenterTitle(fCnt);
  pTrkNfitVsPt -> SetLineColor(fColAll);
  pTrkNfitVsPt -> SetLineStyle(fLinAll);
  pTrkNfitVsPt -> SetFillColor(fColAll);
  pTrkNfitVsPt -> SetFillStyle(fFilAll);
  pTrkNfitVsPt -> SetMarkerColor(fColAll);
  pTrkNfitVsPt -> SetMarkerStyle(fMarAll);
  pTrkNfitVsPt -> SetTitle("");
  pTrkNfitVsPt -> SetTitleFont(fTxt);
  pTrkNfitVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  pTrkNfitVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  pTrkNfitVsPt -> GetXaxis() -> SetTitleSize(fTit);
  pTrkNfitVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  pTrkNfitVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  pTrkNfitVsPt -> GetXaxis() -> SetLabelSize(fLab);
  pTrkNfitVsPt -> GetXaxis() -> CenterTitle(fCnt);
  pTrkNfitVsPt -> GetYaxis() -> SetTitle(sTrkNfit.Data());
  pTrkNfitVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  pTrkNfitVsPt -> GetYaxis() -> SetTitleSize(fTit);
  pTrkNfitVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  pTrkNfitVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  pTrkNfitVsPt -> GetYaxis() -> SetLabelSize(fLab);
  pTrkNfitVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkNfitVsPt -> SetTitle("");
  hTrkNfitVsPt -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkNfitVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetXaxis() -> SetTitleSize(fTit);
  hTrkNfitVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkNfitVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkNfitVsPt -> GetXaxis() -> SetLabelSize(fLab);
  hTrkNfitVsPt -> GetXaxis() -> CenterTitle(fCnt);
  hTrkNfitVsPt -> GetYaxis() -> SetTitle(sTrkNfit.Data());
  hTrkNfitVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetYaxis() -> SetTitleSize(fTit);
  hTrkNfitVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkNfitVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkNfitVsPt -> GetYaxis() -> SetLabelSize(fLab);
  hTrkNfitVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkNfitVsPt -> GetZaxis() -> SetTitle("");
  hTrkNfitVsPt -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetZaxis() -> SetTitleSize(fTit);
  hTrkNfitVsPt -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkNfitVsPt -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkNfitVsPt -> GetZaxis() -> SetLabelSize(fLab);
  hTrkNfitVsPt -> GetZaxis() -> CenterTitle(fCnt);

  hTrkRfit[0] -> SetLineColor(fColAll);
  hTrkRfit[0] -> SetLineStyle(fLinAll);
  hTrkRfit[0] -> SetFillColor(fColFil);
  hTrkRfit[0] -> SetFillStyle(fFilAll);
  hTrkRfit[0] -> SetMarkerColor(fColAll);
  hTrkRfit[0] -> SetMarkerStyle(fMarAll);
  hTrkRfit[0] -> SetTitle("");
  hTrkRfit[0] -> SetTitleFont(fTxt);
  hTrkRfit[0] -> GetXaxis() -> SetTitle(sTrkRfit.Data());
  hTrkRfit[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkRfit[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkRfit[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkRfit[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkRfit[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkRfit[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkRfit[0] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkRfit[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkRfit[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkRfit[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkRfit[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkRfit[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkRfit[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkRfit[1] -> SetLineColor(fColCut);
  hTrkRfit[1] -> SetLineStyle(fLinCut);
  hTrkRfit[1] -> SetFillColor(fColCut);
  hTrkRfit[1] -> SetFillStyle(fFilCut);
  hTrkRfit[1] -> SetMarkerColor(fColCut);
  hTrkRfit[1] -> SetMarkerStyle(fMarCut);
  hTrkRfit[1] -> SetTitle("");
  hTrkRfit[1] -> SetTitleFont(fTxt);
  hTrkRfit[1] -> GetXaxis() -> SetTitle(sTrkRfit.Data());
  hTrkRfit[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkRfit[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkRfit[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkRfit[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkRfit[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkRfit[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkRfit[1] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkRfit[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkRfit[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkRfit[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkRfit[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkRfit[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkRfit[1] -> GetYaxis() -> CenterTitle(fCnt);

  // dca
  hTrkDca[0]  -> SetLineColor(fColAll);
  hTrkDca[0]  -> SetLineStyle(fLinAll);
  hTrkDca[0]  -> SetFillColor(fColFil);
  hTrkDca[0]  -> SetFillStyle(fFilAll);
  hTrkDca[0]  -> SetMarkerColor(fColAll);
  hTrkDca[0]  -> SetMarkerStyle(fMarAll);
  hTrkDca[0]  -> SetTitle("");
  hTrkDca[0]  -> SetTitleFont(fTxt);
  hTrkDca[0]  -> GetXaxis() -> SetTitle(sTrkDca.Data());
  hTrkDca[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDca[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDca[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDca[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDca[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDca[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDca[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkDca[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDca[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDca[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDca[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDca[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDca[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDca[1]  -> SetLineColor(fColCut);
  hTrkDca[1]  -> SetLineStyle(fLinCut);
  hTrkDca[1]  -> SetFillColor(fColCut);
  hTrkDca[1]  -> SetFillStyle(fFilCut);
  hTrkDca[1]  -> SetMarkerColor(fColCut);
  hTrkDca[1]  -> SetMarkerStyle(fMarCut);
  hTrkDca[1]  -> SetTitle("");
  hTrkDca[1]  -> SetTitleFont(fTxt);
  hTrkDca[1]  -> GetXaxis() -> SetTitle(sTrkDca.Data());
  hTrkDca[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDca[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDca[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDca[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDca[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDca[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDca[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkDca[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDca[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDca[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDca[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDca[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDca[1]  -> GetYaxis() -> CenterTitle(fCnt);
  pTrkDcaVsPt -> SetLineColor(fColAll);
  pTrkDcaVsPt -> SetLineStyle(fLinAll);
  pTrkDcaVsPt -> SetFillColor(fColAll);
  pTrkDcaVsPt -> SetFillStyle(fFilAll);
  pTrkDcaVsPt -> SetMarkerColor(fColAll);
  pTrkDcaVsPt -> SetMarkerStyle(fMarAll);
  pTrkDcaVsPt -> SetTitle("");
  pTrkDcaVsPt -> SetTitleFont(fTxt);
  pTrkDcaVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  pTrkDcaVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  pTrkDcaVsPt -> GetXaxis() -> SetTitleSize(fTit);
  pTrkDcaVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  pTrkDcaVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  pTrkDcaVsPt -> GetXaxis() -> SetLabelSize(fLab);
  pTrkDcaVsPt -> GetXaxis() -> CenterTitle(fCnt);
  pTrkDcaVsPt -> GetYaxis() -> SetTitle(sTrkDca.Data());
  pTrkDcaVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  pTrkDcaVsPt -> GetYaxis() -> SetTitleSize(fTit);
  pTrkDcaVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  pTrkDcaVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  pTrkDcaVsPt -> GetYaxis() -> SetLabelSize(fLab);
  pTrkDcaVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDcaVsPt -> SetTitle("");
  hTrkDcaVsPt -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkDcaVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDcaVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDcaVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDcaVsPt -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDcaVsPt -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDcaVsPt -> GetYaxis() -> SetTitle(sTrkDca.Data());
  hTrkDcaVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDcaVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDcaVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDcaVsPt -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDcaVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDcaVsPt -> GetZaxis() -> SetTitle("");
  hTrkDcaVsPt -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetZaxis() -> SetTitleSize(fTit);
  hTrkDcaVsPt -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkDcaVsPt -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkDcaVsPt -> GetZaxis() -> SetLabelSize(fLab);
  hTrkDcaVsPt -> GetZaxis() -> CenterTitle(fCnt);

  // eta
  hTrkEta[0] -> SetLineColor(fColAll);
  hTrkEta[0] -> SetLineStyle(fLinAll);
  hTrkEta[0] -> SetFillColor(fColFil);
  hTrkEta[0] -> SetFillStyle(fFilAll);
  hTrkEta[0] -> SetMarkerColor(fColAll);
  hTrkEta[0] -> SetMarkerStyle(fMarAll);
  hTrkEta[0] -> SetTitle("");
  hTrkEta[0] -> SetTitleFont(fTxt);
  hTrkEta[0] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkEta[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkEta[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkEta[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkEta[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkEta[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkEta[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkEta[0] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkEta[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkEta[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkEta[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkEta[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkEta[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkEta[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkEta[1] -> SetLineColor(fColCut);
  hTrkEta[1] -> SetLineStyle(fLinCut);
  hTrkEta[1] -> SetFillColor(fColCut);
  hTrkEta[1] -> SetFillStyle(fFilCut);
  hTrkEta[1] -> SetMarkerColor(fColCut);
  hTrkEta[1] -> SetMarkerStyle(fMarCut);
  hTrkEta[1] -> SetTitle("");
  hTrkEta[1] -> SetTitleFont(fTxt);
  hTrkEta[1] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkEta[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkEta[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkEta[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkEta[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkEta[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkEta[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkEta[1] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkEta[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkEta[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkEta[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkEta[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkEta[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkEta[1] -> GetYaxis() -> CenterTitle(fCnt);

  // pTtrk
  hTrkPt[0]    -> SetLineColor(fColAll);
  hTrkPt[0]    -> SetLineStyle(fLinAll);
  hTrkPt[0]    -> SetFillColor(fColFil);
  hTrkPt[0]    -> SetFillStyle(fFilAll);
  hTrkPt[0]    -> SetMarkerColor(fColAll);
  hTrkPt[0]    -> SetMarkerStyle(fMarAll);
  hTrkPt[0]    -> SetTitle("");
  hTrkPt[0]    -> SetTitleFont(fTxt);
  hTrkPt[0]    -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPt[0]    -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPt[0]    -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPt[0]    -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPt[0]    -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPt[0]    -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPt[0]    -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPt[0]    -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkPt[0]    -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPt[0]    -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPt[0]    -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPt[0]    -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPt[0]    -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPt[0]    -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPt[1]    -> SetLineColor(fColCut);
  hTrkPt[1]    -> SetLineStyle(fLinCut);
  hTrkPt[1]    -> SetFillColor(fColCut);
  hTrkPt[1]    -> SetFillStyle(fFilCut);
  hTrkPt[1]    -> SetMarkerColor(fColCut);
  hTrkPt[1]    -> SetMarkerStyle(fMarCut);
  hTrkPt[1]    -> SetTitle("");
  hTrkPt[1]    -> SetTitleFont(fTxt);
  hTrkPt[1]    -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPt[1]    -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPt[1]    -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPt[1]    -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPt[1]    -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPt[1]    -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPt[1]    -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPt[1]    -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkPt[1]    -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPt[1]    -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPt[1]    -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPt[1]    -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPt[1]    -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPt[1]    -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtBin[0] -> SetLineColor(fColTrg[0]);
  hTrkPtBin[0] -> SetLineStyle(fLinAll);
  hTrkPtBin[0] -> SetFillColor(fColTrg[0]);
  hTrkPtBin[0] -> SetFillStyle(fFilAll);
  hTrkPtBin[0] -> SetMarkerColor(fColTrg[0]);
  hTrkPtBin[0] -> SetMarkerStyle(fMarAll);
  hTrkPtBin[0] -> SetTitle("");
  hTrkPtBin[0] -> SetTitleFont(fTxt);
  hTrkPtBin[0] -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtBin[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtBin[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtBin[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtBin[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtBin[0] -> GetYaxis() -> SetTitle(sTrkPtY.Data());
  hTrkPtBin[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtBin[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtBin[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtBin[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtBin[1] -> SetLineColor(fColTrg[1]);
  hTrkPtBin[1] -> SetLineStyle(fLinCut);
  hTrkPtBin[1] -> SetFillColor(fColTrg[1]);
  hTrkPtBin[1] -> SetFillStyle(fFilCut);
  hTrkPtBin[1] -> SetMarkerColor(fColTrg[1]);
  hTrkPtBin[1] -> SetMarkerStyle(fMarCut);
  hTrkPtBin[1] -> SetTitle("");
  hTrkPtBin[1] -> SetTitleFont(fTxt);
  hTrkPtBin[1] -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtBin[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtBin[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtBin[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtBin[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtBin[1] -> GetYaxis() -> SetTitle(sTrkPtY.Data());
  hTrkPtBin[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtBin[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtBin[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtBin[1] -> GetYaxis() -> CenterTitle(fCnt);

  // dFtrk
  hTrkDfBin[0]   -> SetLineColor(fColTrg[0]);
  hTrkDfBin[0]   -> SetLineStyle(fLinAll);
  hTrkDfBin[0]   -> SetFillColor(fColTrg[0]);
  hTrkDfBin[0]   -> SetFillStyle(fFilAll);
  hTrkDfBin[0]   -> SetMarkerColor(fColTrg[0]);
  hTrkDfBin[0]   -> SetMarkerStyle(fMarAll);
  hTrkDfBin[0]   -> SetTitle("");
  hTrkDfBin[0]   -> SetTitleFont(fTxt);
  hTrkDfBin[0]   -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfBin[0]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[0]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfBin[0]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfBin[0]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[0]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfBin[0]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfBin[0]   -> GetYaxis() -> SetTitle(sTrkDfY.Data());
  hTrkDfBin[0]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[0]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfBin[0]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfBin[0]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[0]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfBin[0]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDfBin[1]   -> SetLineColor(fColTrg[1]);
  hTrkDfBin[1]   -> SetLineStyle(fLinCut);
  hTrkDfBin[1]   -> SetFillColor(fColTrg[1]);
  hTrkDfBin[1]   -> SetFillStyle(fFilCut);
  hTrkDfBin[1]   -> SetMarkerColor(fColTrg[1]);
  hTrkDfBin[1]   -> SetMarkerStyle(fMarCut);
  hTrkDfBin[1]   -> SetTitle("");
  hTrkDfBin[1]   -> SetTitleFont(fTxt);
  hTrkDfBin[1]   -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfBin[1]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[1]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfBin[1]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfBin[1]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[1]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfBin[1]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfBin[1]   -> GetYaxis() -> SetTitle(sTrkDfY.Data());
  hTrkDfBin[1]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[1]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfBin[1]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfBin[1]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[1]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfBin[1]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[0]  -> SetLineColor(fColTrg[0]);
  hTrkPtVsDf[0]  -> SetLineStyle(fLinAll);
  hTrkPtVsDf[0]  -> SetFillColor(fColTrg[0]);
  hTrkPtVsDf[0]  -> SetFillStyle(fFilAll);
  hTrkPtVsDf[0]  -> SetMarkerColor(fColTrg[0]);
  hTrkPtVsDf[0]  -> SetMarkerStyle(fMarAll);
  hTrkPtVsDf[0]  -> SetTitle("");
  hTrkPtVsDf[0]  -> SetTitleFont(fTxt);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[1]  -> SetLineColor(fColTrg[1]);
  hTrkPtVsDf[1]  -> SetLineStyle(fLinCut);
  hTrkPtVsDf[1]  -> SetFillColor(fColTrg[1]);
  hTrkPtVsDf[1]  -> SetFillStyle(fFilCut);
  hTrkPtVsDf[1]  -> SetMarkerColor(fColTrg[1]);
  hTrkPtVsDf[1]  -> SetMarkerStyle(fMarCut);
  hTrkPtVsDf[1]  -> SetTitle("");
  hTrkPtVsDf[1]  -> SetTitleFont(fTxt);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[1]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[0] -> SetLineColor(fColTrg[0]);
  hTrkDfVsEta[0] -> SetLineStyle(fLinAll);
  hTrkDfVsEta[0] -> SetFillColor(fColTrg[0]);
  hTrkDfVsEta[0] -> SetFillStyle(fFilAll);
  hTrkDfVsEta[0] -> SetMarkerColor(fColTrg[0]);
  hTrkDfVsEta[0] -> SetMarkerStyle(fMarAll);
  hTrkDfVsEta[0] -> SetTitle("");
  hTrkDfVsEta[0] -> SetTitleFont(fTxt);
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfVsEta[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfVsEta[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[1] -> SetLineColor(fColTrg[1]);
  hTrkDfVsEta[1] -> SetLineStyle(fLinCut);
  hTrkDfVsEta[1] -> SetFillColor(fColTrg[1]);
  hTrkDfVsEta[1] -> SetFillStyle(fFilCut);
  hTrkDfVsEta[1] -> SetMarkerColor(fColTrg[1]);
  hTrkDfVsEta[1] -> SetMarkerStyle(fMarCut);
  hTrkDfVsEta[1] -> SetTitle("");
  hTrkDfVsEta[1] -> SetTitleFont(fTxt);
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfVsEta[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfVsEta[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[1] -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Styles set." << endl;


  // draw plots
  const UInt_t  width(750);
  const UInt_t  bigWidth(1500);
  const UInt_t  height(750);
  const UInt_t  bigHeight(1500);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fFrame(0);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fLogY2(0);
  const UInt_t  fLogZ(1);
  const UInt_t  fTick(1);
  const Float_t fMargin0(0.);
  const Float_t fMarginBig(0.15);
  const Float_t fMarginSmall(0.05);

  // evt no's
  TCanvas *cEvtNum = new TCanvas("cEvtNum", "", width, height);
  cEvtNum -> SetLogx(fLogX);
  cEvtNum -> SetLogy(fLogY);
  cEvtNum -> SetGrid(fGrid, fGrid);
  cEvtNum -> SetTicks(fTick, fTick);
  cEvtNum -> SetBorderMode(fMode);
  cEvtNum -> SetBorderSize(fBord);
  cEvtNum -> SetFrameBorderMode(fFrame);
  cEvtNum -> SetLeftMargin(fMarginBig);
  cEvtNum -> SetTopMargin(fMarginSmall);
  cEvtNum -> SetRightMargin(fMarginSmall);
  cEvtNum -> SetBottomMargin(fMarginBig);
  cEvtNum -> cd();
  hEvtNum -> Draw("B");
  fOutput -> cd();
  cEvtNum -> Write();
  cEvtNum -> Close();

  // vertices
  TCanvas *cEvtVz = new TCanvas("cEvtVz", "", width, height);
  cEvtVz    -> SetLogx(fLogX);
  cEvtVz    -> SetLogy(fLogY);
  cEvtVz    -> SetGrid(fGrid, fGrid);
  cEvtVz    -> SetTicks(fTick, fTick);
  cEvtVz    -> SetBorderMode(fMode);
  cEvtVz    -> SetBorderSize(fBord);
  cEvtVz    -> SetFrameBorderMode(fFrame);
  cEvtVz    -> SetLeftMargin(fMarginBig);
  cEvtVz    -> SetTopMargin(fMarginSmall);
  cEvtVz    -> SetRightMargin(fMarginSmall);
  cEvtVz    -> SetBottomMargin(fMarginBig);
  cEvtVz    -> cd();
  hEvtVz[0] -> Draw("");
  hEvtVz[1] -> Draw("hist same");
  hEvtVz[1] -> Draw("same");
  cEvtVz    -> Write();
  cEvtVz    -> Close();

  TCanvas *cEvtVyVxVsVz = new TCanvas("cEvtVyVxVsVz", "", bigWidth, height);
  TPad    *pEvtVxVsVz   = new TPad("pEvtVxVsVz", "", 0., 0., 0.5, 1.);
  TPad    *pEvtVyVsVz   = new TPad("pEvtVyVsVz", "", 0.5, 0., 1., 1.);
  pEvtVxVsVz   -> SetLogx(fLogX);
  pEvtVxVsVz   -> SetLogy(fLogY2);
  pEvtVxVsVz   -> SetLogz(fLogZ);
  pEvtVxVsVz   -> SetGrid(fGrid, fGrid);
  pEvtVxVsVz   -> SetTicks(fTick, fTick);
  pEvtVxVsVz   -> SetBorderMode(fMode);
  pEvtVxVsVz   -> SetBorderSize(fBord);
  pEvtVxVsVz   -> SetFrameBorderMode(fFrame);
  pEvtVxVsVz   -> SetLeftMargin(fMarginBig);
  pEvtVxVsVz   -> SetTopMargin(fMarginSmall);
  pEvtVxVsVz   -> SetRightMargin(fMarginBig);
  pEvtVxVsVz   -> SetBottomMargin(fMarginBig);
  pEvtVyVsVz   -> SetLogx(fLogX);
  pEvtVyVsVz   -> SetLogy(fLogY2);
  pEvtVyVsVz   -> SetLogz(fLogZ);
  pEvtVyVsVz   -> SetGrid(fGrid, fGrid);
  pEvtVyVsVz   -> SetTicks(fTick, fTick);
  pEvtVyVsVz   -> SetBorderMode(fMode);
  pEvtVyVsVz   -> SetBorderSize(fBord);
  pEvtVyVsVz   -> SetFrameBorderMode(fFrame);
  pEvtVyVsVz   -> SetLeftMargin(fMarginBig);
  pEvtVyVsVz   -> SetTopMargin(fMarginSmall);
  pEvtVyVsVz   -> SetRightMargin(fMarginBig);
  pEvtVyVsVz   -> SetBottomMargin(fMarginBig);
  cEvtVyVxVsVz -> cd();
  pEvtVxVsVz   -> Draw();
  pEvtVyVsVz   -> Draw();
  pEvtVxVsVz   -> cd();
  hEvtVxVsVz   -> Draw("colz");
  pEvtVyVsVz   -> cd();
  hEvtVyVsVz   -> Draw("colz");
  cEvtVyVxVsVz -> Write();
  cEvtVyVxVsVz -> Close();

  TCanvas *cEvtVyVsVx = new TCanvas("cEvtVyVsVx", "", width, height);
  cEvtVyVsVx -> SetLogx(fLogX);
  cEvtVyVsVx -> SetLogy(fLogY2);
  cEvtVyVsVx -> SetLogz(fLogZ);
  cEvtVyVsVx -> SetGrid(fGrid, fGrid);
  cEvtVyVsVx -> SetTicks(fTick, fTick);
  cEvtVyVsVx -> SetBorderMode(fMode);
  cEvtVyVsVx -> SetBorderSize(fBord);
  cEvtVyVsVx -> SetFrameBorderMode(fFrame);
  cEvtVyVsVx -> SetLeftMargin(fMarginBig);
  cEvtVyVsVx -> SetTopMargin(fMarginSmall);
  cEvtVyVsVx -> SetRightMargin(fMarginBig);
  cEvtVyVsVx -> SetBottomMargin(fMarginBig);
  cEvtVyVsVx -> cd();
  hEvtVyVsVx -> Draw("colz");
  cEvtVyVsVx -> Write();
  cEvtVyVsVx -> Close();

  // no. of tracks
  TCanvas *cEvtPrim = new TCanvas("cEvtPrim", "", width, height);
  cEvtPrim -> SetLogx(fLogX);
  cEvtPrim -> SetLogy(fLogY);
  cEvtPrim -> SetGrid(fGrid, fGrid);
  cEvtPrim -> SetTicks(fTick, fTick);
  cEvtPrim -> SetBorderMode(fMode);
  cEvtPrim -> SetBorderSize(fBord);
  cEvtPrim -> SetFrameBorderMode(fFrame);
  cEvtPrim -> SetLeftMargin(fMarginBig);
  cEvtPrim -> SetTopMargin(fMarginSmall);
  cEvtPrim -> SetRightMargin(fMarginSmall);
  cEvtPrim -> SetBottomMargin(fMarginBig);
  cEvtPrim -> cd();
  hEvtPrim -> Draw();
  cEvtPrim -> Write();
  cEvtPrim -> Close();

  // eTtrg (no tsp cuts)
  TCanvas *cTrgEt = new TCanvas("cTrgEt", "", width, height);
  cTrgEt    -> SetLogx(fLogX);
  cTrgEt    -> SetLogy(fLogY);
  cTrgEt    -> SetGrid(fGrid, fGrid);
  cTrgEt    -> SetTicks(fTick, fTick);
  cTrgEt    -> SetBorderMode(fMode);
  cTrgEt    -> SetBorderSize(fBord);
  cTrgEt    -> SetFrameBorderMode(fFrame);
  cTrgEt    -> SetLeftMargin(fMarginBig);
  cTrgEt    -> SetTopMargin(fMarginSmall);
  cTrgEt    -> SetRightMargin(fMarginSmall);
  cTrgEt    -> SetBottomMargin(fMarginBig);
  cTrgEt    -> cd();
  hTrgEt[0] -> Draw();
  hTrgEt[1] -> Draw("hist same");
  hTrgEt[1] -> Draw("same");
  cTrgEt    -> Write();
  cTrgEt    -> Close();

  TCanvas *cTrgEtVsEtaPhi = new TCanvas("cTrgEtVsEtaPhi", "", bigWidth, height);
  TPad    *pTrgEtVsEta    = new TPad("pTrgEtVsEta", "", 0., 0., 0.5, 1.);
  TPad    *pTrgEtVsPhi    = new TPad("pTrgEtVsPhi", "", 0.5, 0., 1., 1.);
  pTrgEtVsEta    -> SetLogx(fLogX);
  pTrgEtVsEta    -> SetLogy(fLogY2);
  pTrgEtVsEta    -> SetLogz(fLogZ);
  pTrgEtVsEta    -> SetGrid(fGrid, fGrid);
  pTrgEtVsEta    -> SetTicks(fTick, fTick);
  pTrgEtVsEta    -> SetBorderMode(fMode);
  pTrgEtVsEta    -> SetBorderSize(fBord);
  pTrgEtVsEta    -> SetFrameBorderMode(fFrame);
  pTrgEtVsEta    -> SetLeftMargin(fMarginBig);
  pTrgEtVsEta    -> SetTopMargin(fMarginSmall);
  pTrgEtVsEta    -> SetRightMargin(fMarginBig);
  pTrgEtVsEta    -> SetBottomMargin(fMarginBig);
  pTrgEtVsPhi    -> SetLogx(fLogX);
  pTrgEtVsPhi    -> SetLogy(fLogY2);
  pTrgEtVsPhi    -> SetLogz(fLogZ);
  pTrgEtVsPhi    -> SetGrid(fGrid, fGrid);
  pTrgEtVsPhi    -> SetTicks(fTick, fTick);
  pTrgEtVsPhi    -> SetBorderMode(fMode);
  pTrgEtVsPhi    -> SetBorderSize(fBord);
  pTrgEtVsPhi    -> SetFrameBorderMode(fFrame);
  pTrgEtVsPhi    -> SetLeftMargin(fMarginBig);
  pTrgEtVsPhi    -> SetTopMargin(fMarginSmall);
  pTrgEtVsPhi    -> SetRightMargin(fMarginBig);
  pTrgEtVsPhi    -> SetBottomMargin(fMarginBig);
  cTrgEtVsEtaPhi -> cd();
  pTrgEtVsEta    -> Draw();
  pTrgEtVsPhi    -> Draw();
  pTrgEtVsEta    -> cd();
  hTrgEtVsEta    -> Draw("colz");
  pTrgEtVsPhi    -> cd();
  hTrgEtVsPhi    -> Draw("colz");
  cTrgEtVsEtaPhi -> Write();
  cTrgEtVsEtaPhi -> Close();

  // eTtrg (tsp cuts)
  TCanvas *cTrgEtBin = new TCanvas("cTrgEtBin", "", bigWidth, height);
  TPad    *pTrgEtPi0 = new TPad("pTrgEtPi0", "", 0., 0., 0.5, 1.);
  TPad    *pTrgEtGam = new TPad("pTrgEtGam", "", 0.5, 0., 1., 1.);
  pTrgEtPi0       -> SetLogx(fLogX);
  pTrgEtPi0       -> SetLogy(fLogY);
  pTrgEtPi0       -> SetGrid(fGrid, fGrid);
  pTrgEtPi0       -> SetTicks(fTick, fTick);
  pTrgEtPi0       -> SetBorderMode(fMode);
  pTrgEtPi0       -> SetBorderSize(fBord);
  pTrgEtPi0       -> SetFrameBorderMode(fFrame);
  pTrgEtPi0       -> SetLeftMargin(fMarginBig);
  pTrgEtPi0       -> SetTopMargin(fMarginSmall);
  pTrgEtPi0       -> SetRightMargin(fMarginSmall);
  pTrgEtPi0       -> SetBottomMargin(fMarginBig);
  pTrgEtGam       -> SetLogx(fLogX);
  pTrgEtGam       -> SetLogy(fLogY);
  pTrgEtGam       -> SetGrid(fGrid, fGrid);
  pTrgEtGam       -> SetTicks(fTick, fTick);
  pTrgEtGam       -> SetBorderMode(fMode);
  pTrgEtGam       -> SetBorderSize(fBord);
  pTrgEtGam       -> SetFrameBorderMode(fFrame);
  pTrgEtGam       -> SetLeftMargin(fMarginBig);
  pTrgEtGam       -> SetTopMargin(fMarginSmall);
  pTrgEtGam       -> SetRightMargin(fMarginSmall);
  pTrgEtGam       -> SetBottomMargin(fMarginBig);
  cTrgEtBin       -> cd();
  pTrgEtPi0       -> Draw();
  pTrgEtGam       -> Draw();
  pTrgEtPi0       -> cd();
  hTrgEtBin[1][0] -> Draw();
  hTrgEtBin[2][0] -> Draw("same");
  hTrgEtBin[3][0] -> Draw("same");
  pTrgEtGam       -> cd();
  hTrgEtBin[1][1] -> Draw();
  hTrgEtBin[2][1] -> Draw("same");
  hTrgEtBin[3][1] -> Draw("same");
  cTrgEtBin       -> Write();
  cTrgEtBin       -> Close();

  TCanvas *cTrgEtSum = new TCanvas("cTrgEtSum", "", width, height);
  cTrgEtSum       -> SetLogx(fLogX);
  cTrgEtSum       -> SetLogy(fLogY);
  cTrgEtSum       -> SetGrid(fGrid, fGrid);
  cTrgEtSum       -> SetTicks(fTick, fTick);
  cTrgEtSum       -> SetBorderMode(fMode);
  cTrgEtSum       -> SetBorderSize(fBord);
  cTrgEtSum       -> SetFrameBorderMode(fFrame);
  cTrgEtSum       -> SetLeftMargin(fMarginBig);
  cTrgEtSum       -> SetTopMargin(fMarginSmall);
  cTrgEtSum       -> SetRightMargin(fMarginSmall);
  cTrgEtSum       -> SetBottomMargin(fMarginBig);
  hTrgEtBin[0][0] -> Draw();
  hTrgEtBin[0][1] -> Draw("same");
  hTrgEtSum       -> Draw("same");
  cTrgEtSum       -> Write();
  cTrgEtSum       -> Close();

  TCanvas *cTrgEtaVsPhi = new TCanvas("cTrgEtaVsPhi", "", width, height);
  cTrgEtaVsPhi -> SetLogx(fLogX);
  cTrgEtaVsPhi -> SetLogy(fLogY2);
  cTrgEtaVsPhi -> SetLogz(fLogZ);
  cTrgEtaVsPhi -> SetGrid(fGrid, fGrid);
  cTrgEtaVsPhi -> SetTicks(fTick, fTick);
  cTrgEtaVsPhi -> SetBorderMode(fMode);
  cTrgEtaVsPhi -> SetBorderSize(fBord);
  cTrgEtaVsPhi -> SetFrameBorderMode(fFrame);
  cTrgEtaVsPhi -> SetLeftMargin(fMarginBig);
  cTrgEtaVsPhi -> SetTopMargin(fMarginSmall);
  cTrgEtaVsPhi -> SetRightMargin(fMarginBig);
  cTrgEtaVsPhi -> SetBottomMargin(fMarginBig);
  hTrgEtaVsPhi -> Draw("colz");
  cTrgEtaVsPhi -> Write();
  cTrgEtaVsPhi -> Close();

  // nFit and nFit/nPoss
  TCanvas *cTrkNfit = new TCanvas("cTrkNfit", "", width, height);
  cTrkNfit    -> SetLogx(fLogX);
  cTrkNfit    -> SetLogy(fLogY);
  cTrkNfit    -> SetGrid(fGrid, fGrid);
  cTrkNfit    -> SetTicks(fTick, fTick);
  cTrkNfit    -> SetBorderMode(fMode);
  cTrkNfit    -> SetBorderSize(fBord);
  cTrkNfit    -> SetFrameBorderMode(fFrame);
  cTrkNfit    -> SetLeftMargin(fMarginBig);
  cTrkNfit    -> SetTopMargin(fMarginSmall);
  cTrkNfit    -> SetRightMargin(fMarginSmall);
  cTrkNfit    -> SetBottomMargin(fMarginBig);
  cTrkNfit    -> cd();
  hTrkNfit[0] -> Draw();
  hTrkNfit[1] -> Draw("hist same");
  hTrkNfit[1] -> Draw("same");
  cTrkNfit    -> Write();
  cTrkNfit    -> Close();

  TCanvas *cTrkNfitVsPt = new TCanvas("cTrkNfitVsPt", "", width, height);
  cTrkNfitVsPt -> SetLogx(fLogX);
  cTrkNfitVsPt -> SetLogy(fLogY2);
  cTrkNfitVsPt -> SetLogz(fLogZ);
  cTrkNfitVsPt -> SetGrid(fGrid, fGrid);
  cTrkNfitVsPt -> SetTicks(fTick, fTick);
  cTrkNfitVsPt -> SetBorderMode(fMode);
  cTrkNfitVsPt -> SetBorderSize(fBord);
  cTrkNfitVsPt -> SetFrameBorderMode(fFrame);
  cTrkNfitVsPt -> SetLeftMargin(fMarginBig);
  cTrkNfitVsPt -> SetTopMargin(fMarginSmall);
  cTrkNfitVsPt -> SetRightMargin(fMarginBig);
  cTrkNfitVsPt -> SetBottomMargin(fMarginBig);
  cTrkNfitVsPt -> cd();
  hTrkNfitVsPt -> Draw("colz");
  pTrkNfitVsPt -> Draw("same");
  cTrkNfitVsPt -> Write();
  cTrkNfitVsPt -> Close();

  TCanvas *cTrkRfit = new TCanvas("cTrkRfit", "", width, height);
  cTrkRfit    -> SetLogx(fLogX);
  cTrkRfit    -> SetLogy(fLogY);
  cTrkRfit    -> SetGrid(fGrid, fGrid);
  cTrkRfit    -> SetTicks(fTick, fTick);
  cTrkRfit    -> SetBorderMode(fMode);
  cTrkRfit    -> SetBorderSize(fBord);
  cTrkRfit    -> SetFrameBorderMode(fFrame);
  cTrkRfit    -> SetLeftMargin(fMarginBig);
  cTrkRfit    -> SetTopMargin(fMarginSmall);
  cTrkRfit    -> SetRightMargin(fMarginSmall);
  cTrkRfit    -> SetBottomMargin(fMarginBig);
  cTrkRfit    -> cd();
  hTrkRfit[0] -> Draw();
  hTrkRfit[1] -> Draw("hist same");
  hTrkRfit[1] -> Draw("same");
  cTrkRfit    -> Write();
  cTrkRfit    -> Close();

  // dca
  TCanvas *cTrkDca = new TCanvas("cTrkDca", "", width, height);
  cTrkDca    -> SetLogx(fLogX);
  cTrkDca    -> SetLogy(fLogY);
  cTrkDca    -> SetGrid(fGrid, fGrid);
  cTrkDca    -> SetTicks(fTick, fTick);
  cTrkDca    -> SetBorderMode(fMode);
  cTrkDca    -> SetBorderSize(fBord);
  cTrkDca    -> SetFrameBorderMode(fFrame);
  cTrkDca    -> SetLeftMargin(fMarginBig);
  cTrkDca    -> SetTopMargin(fMarginSmall);
  cTrkDca    -> SetRightMargin(fMarginSmall);
  cTrkDca    -> SetBottomMargin(fMarginBig);
  cTrkDca    -> cd();
  hTrkDca[0] -> Draw();
  hTrkDca[1] -> Draw("hist same");
  hTrkDca[1] -> Draw("same");
  cTrkDca    -> Write();
  cTrkDca    -> Close();

  TCanvas *cTrkDcaVsPt = new TCanvas("cTrkDcaVsPt", "", width, height);
  cTrkDcaVsPt -> SetLogx(fLogX);
  cTrkDcaVsPt -> SetLogy(fLogY2);
  cTrkDcaVsPt -> SetLogz(fLogZ);
  cTrkDcaVsPt -> SetGrid(fGrid, fGrid);
  cTrkDcaVsPt -> SetTicks(fTick, fTick);
  cTrkDcaVsPt -> SetBorderMode(fMode);
  cTrkDcaVsPt -> SetBorderSize(fBord);
  cTrkDcaVsPt -> SetFrameBorderMode(fFrame);
  cTrkDcaVsPt -> SetLeftMargin(fMarginBig);
  cTrkDcaVsPt -> SetTopMargin(fMarginSmall);
  cTrkDcaVsPt -> SetRightMargin(fMarginBig);
  cTrkDcaVsPt -> SetBottomMargin(fMarginBig);
  cTrkDcaVsPt -> cd();
  hTrkDcaVsPt -> Draw("colz");
  pTrkDcaVsPt -> Draw("same");
  cTrkDcaVsPt -> Write();
  cTrkDcaVsPt -> Close();

  // eta
  TCanvas *cTrkEta = new TCanvas("cTrkEta", "", width, height);
  cTrkEta    -> SetLogx(fLogX);
  cTrkEta    -> SetLogy(fLogY);
  cTrkEta    -> SetGrid(fGrid, fGrid);
  cTrkEta    -> SetTicks(fTick, fTick);
  cTrkEta    -> SetBorderMode(fMode);
  cTrkEta    -> SetBorderSize(fBord);
  cTrkEta    -> SetFrameBorderMode(fFrame);
  cTrkEta    -> SetLeftMargin(fMarginBig);
  cTrkEta    -> SetTopMargin(fMarginSmall);
  cTrkEta    -> SetRightMargin(fMarginSmall);
  cTrkEta    -> SetBottomMargin(fMarginBig);
  cTrkEta    -> cd();
  hTrkEta[0] -> Draw();
  hTrkEta[1] -> Draw("hist same");
  hTrkEta[1] -> Draw("same");
  cTrkEta    -> Write();
  cTrkEta    -> Close();

  // pTtrk
  TCanvas *cTrkPt = new TCanvas("cTrkPt", "", width, height);
  cTrkPt    -> SetLogx(fLogX);
  cTrkPt    -> SetLogy(fLogY);
  cTrkPt    -> SetGrid(fGrid, fGrid);
  cTrkPt    -> SetTicks(fTick, fTick);
  cTrkPt    -> SetBorderMode(fMode);
  cTrkPt    -> SetBorderSize(fBord);
  cTrkPt    -> SetFrameBorderMode(fFrame);
  cTrkPt    -> SetLeftMargin(fMarginBig);
  cTrkPt    -> SetTopMargin(fMarginSmall);
  cTrkPt    -> SetRightMargin(fMarginSmall);
  cTrkPt    -> SetBottomMargin(fMarginBig);
  cTrkPt    -> cd();
  hTrkPt[0] -> Draw();
  hTrkPt[1] -> Draw("hist same");
  hTrkPt[1] -> Draw("same");
  cTrkPt    -> Write();
  cTrkPt    -> Close();

  TCanvas *cTrkPtBin = new TCanvas("cTrkPtBin", "", width, height);
  cTrkPtBin    -> SetLogx(fLogX);
  cTrkPtBin    -> SetLogy(fLogY);
  cTrkPtBin    -> SetGrid(fGrid, fGrid);
  cTrkPtBin    -> SetTicks(fTick, fTick);
  cTrkPtBin    -> SetBorderMode(fMode);
  cTrkPtBin    -> SetBorderSize(fBord);
  cTrkPtBin    -> SetFrameBorderMode(fFrame);
  cTrkPtBin    -> SetLeftMargin(fMarginBig);
  cTrkPtBin    -> SetTopMargin(fMarginSmall);
  cTrkPtBin    -> SetRightMargin(fMarginSmall);
  cTrkPtBin    -> SetBottomMargin(fMarginBig);
  cTrkPtBin    -> cd();
  hTrkPtBin[0] -> Draw();
  hTrkPtBin[1] -> Draw("same");
  cTrkPtBin    -> Write();
  cTrkPtBin    -> Close();

  // dFtrk
  TCanvas *cTrkDf     = new TCanvas("cTrkDf", "", bigWidth, bigHeight);
  TPad    *pDfVsPtPi0 = new TPad("pDfVsPtPi0", "", 0., 0., 0.5, 0.5);
  TPad    *pDfVsPtGam = new TPad("pDfVsPtGam", "", 0.5, 0., 1., 0.5);
  TPad    *pDfPi0     = new TPad("pDfPi0",     "", 0., 0.5, 0.5, 1.);
  TPad    *pDfGam     = new TPad("pDfGam",     "", 0.5, 0.5, 1., 1.);
  pDfVsPtPi0    -> SetLogx(fLogX);
  pDfVsPtPi0    -> SetLogy(fLogY2);
  pDfVsPtPi0    -> SetLogz(fLogZ);
  pDfVsPtPi0    -> SetGrid(fGrid, fGrid);
  pDfVsPtPi0    -> SetTicks(fTick, fTick);
  pDfVsPtPi0    -> SetBorderMode(fMode);
  pDfVsPtPi0    -> SetBorderSize(fBord);
  pDfVsPtPi0    -> SetFrameBorderMode(fFrame);
  pDfVsPtPi0    -> SetLeftMargin(fMarginBig);
  pDfVsPtPi0    -> SetTopMargin(fMargin0);
  pDfVsPtPi0    -> SetRightMargin(fMarginBig);
  pDfVsPtPi0    -> SetBottomMargin(fMarginBig);
  pDfVsPtGam    -> SetLogx(fLogX);
  pDfVsPtGam    -> SetLogy(fLogY2);
  pDfVsPtGam    -> SetLogz(fLogZ);
  pDfVsPtGam    -> SetGrid(fGrid, fGrid);
  pDfVsPtGam    -> SetTicks(fTick, fTick);
  pDfVsPtGam    -> SetBorderMode(fMode);
  pDfVsPtGam    -> SetBorderSize(fBord);
  pDfVsPtGam    -> SetFrameBorderMode(fFrame);
  pDfVsPtGam    -> SetLeftMargin(fMarginBig);
  pDfVsPtGam    -> SetTopMargin(fMargin0);
  pDfVsPtGam    -> SetRightMargin(fMarginBig);
  pDfVsPtGam    -> SetBottomMargin(fMarginBig);
  pDfPi0        -> SetLogx(fLogX);
  pDfPi0        -> SetLogy(fLogY2);
  pDfPi0        -> SetGrid(fGrid, fGrid);
  pDfPi0        -> SetTicks(fTick, fTick);
  pDfPi0        -> SetBorderMode(fMode);
  pDfPi0        -> SetBorderSize(fBord);
  pDfPi0        -> SetFrameBorderMode(fFrame);
  pDfPi0        -> SetLeftMargin(fMarginBig);
  pDfPi0        -> SetTopMargin(fMarginSmall);
  pDfPi0        -> SetRightMargin(fMarginBig);
  pDfPi0        -> SetBottomMargin(fMargin0);
  pDfGam        -> SetLogx(fLogX);
  pDfGam        -> SetLogy(fLogY2);
  pDfGam        -> SetGrid(fGrid, fGrid);
  pDfGam        -> SetTicks(fTick, fTick);
  pDfGam        -> SetBorderMode(fMode);
  pDfGam        -> SetBorderSize(fBord);
  pDfGam        -> SetFrameBorderMode(fFrame);
  pDfGam        -> SetLeftMargin(fMarginBig);
  pDfGam        -> SetTopMargin(fMarginSmall);
  pDfGam        -> SetRightMargin(fMarginBig);
  pDfGam        -> SetBottomMargin(fMargin0);
  cTrkDf        -> cd();
  pDfVsPtPi0    -> Draw();
  pDfVsPtGam    -> Draw();
  pDfPi0        -> Draw();
  pDfGam        -> Draw();
  pDfVsPtPi0    -> cd();
  hTrkPtVsDf[0] -> Draw("colz");
  pDfVsPtGam    -> cd();
  hTrkPtVsDf[1] -> Draw("colz");
  pDfPi0        -> cd();
  hTrkDfBin[0]  -> Draw();
  pDfGam        -> cd();
  hTrkDfBin[1]  -> Draw();
  cTrkDf        -> Write();
  cTrkDf        -> Close();

  TCanvas *cTrkDfVsEta = new TCanvas("cTrkDfVsEta", "", bigWidth, height);
  TPad    *pDfVsEtaPi0 = new TPad("pDfVsEtaPi0", "", 0., 0., 0.5, 1.);
  TPad    *pDfVsEtaGam = new TPad("pDfVsEtaGam", "", 0.5, 0., 1., 1.);
  pDfVsEtaPi0    -> SetLogx(fLogX);
  pDfVsEtaPi0    -> SetLogy(fLogY2);
  pDfVsEtaPi0    -> SetLogz(fLogZ);
  pDfVsEtaPi0    -> SetGrid(fGrid, fGrid);
  pDfVsEtaPi0    -> SetTicks(fTick, fTick);
  pDfVsEtaPi0    -> SetBorderMode(fMode);
  pDfVsEtaPi0    -> SetBorderSize(fBord);
  pDfVsEtaPi0    -> SetFrameBorderMode(fFrame);
  pDfVsEtaPi0    -> SetLeftMargin(fMarginBig);
  pDfVsEtaPi0    -> SetTopMargin(fMarginSmall);
  pDfVsEtaPi0    -> SetRightMargin(fMarginBig);
  pDfVsEtaPi0    -> SetBottomMargin(fMarginBig);
  pDfVsEtaGam    -> SetLogx(fLogX);
  pDfVsEtaGam    -> SetLogy(fLogY2);
  pDfVsEtaGam    -> SetLogz(fLogZ);
  pDfVsEtaGam    -> SetGrid(fGrid, fGrid);
  pDfVsEtaGam    -> SetTicks(fTick, fTick);
  pDfVsEtaGam    -> SetBorderMode(fMode);
  pDfVsEtaGam    -> SetBorderSize(fBord);
  pDfVsEtaGam    -> SetFrameBorderMode(fFrame);
  pDfVsEtaGam    -> SetLeftMargin(fMarginBig);
  pDfVsEtaGam    -> SetTopMargin(fMarginSmall);
  pDfVsEtaGam    -> SetRightMargin(fMarginBig);
  pDfVsEtaGam    -> SetBottomMargin(fMarginBig);
  cTrkDfVsEta    -> cd();
  pDfVsEtaPi0    -> Draw();
  pDfVsEtaGam    -> Draw();
  pDfVsEtaPi0    -> cd();
  hTrkDfVsEta[0] -> Draw("colz");
  pDfVsEtaGam    -> cd();
  hTrkDfVsEta[1] -> Draw("colz");
  cTrkDfVsEta    -> Write();
  cTrkDfVsEta    -> Close();
  cout << "    Made plots." << endl;


  // close files
  fOutput         -> cd();
  hEvtNum         -> Write();
  hEvtVz[0]       -> Write();
  hEvtVz[1]       -> Write();
  hEvtVr[0]       -> Write();
  hEvtVr[1]       -> Write();
  hEvtPrim        -> Write();
  hTrgEt[0]       -> Write();
  hTrgEt[1]       -> Write();
  hTrgEta         -> Write();
  hTrgPhi         -> Write();
  hTrgTsp[0]      -> Write();
  hTrgTsp[1]      -> Write();
  hTrgTsp[2]      -> Write();
  hTrgEtBin[0][0] -> Write();
  hTrgEtBin[0][1] -> Write();
  hTrgEtBin[1][0] -> Write();
  hTrgEtBin[1][1] -> Write();
  hTrgEtBin[2][0] -> Write();
  hTrgEtBin[2][1] -> Write();
  hTrgEtBin[3][0] -> Write();
  hTrgEtBin[3][1] -> Write();
  hTrgEtSum       -> Write();
  hTrkNfit[0]     -> Write();
  hTrkNfit[1]     -> Write();
  hTrkRfit[0]     -> Write();
  hTrkRfit[1]     -> Write();
  hTrkDca[0]      -> Write();
  hTrkDca[1]      -> Write();
  hTrkEta[0]      -> Write();
  hTrkEta[1]      -> Write();
  hTrkPt[0]       -> Write();
  hTrkPt[1]       -> Write();
  hTrkPtBin[0]    -> Write();
  hTrkPtBin[1]    -> Write();
  hTrkDfBin[0]    -> Write();
  hTrkDfBin[1]    -> Write();
  hEvtVyVsVx      -> Write();
  hEvtVxVsVz      -> Write();
  hEvtVyVsVz      -> Write();
  hTrgEtVsEta     -> Write();
  hTrgEtVsPhi     -> Write();
  hTrgEtaVsPhi    -> Write();
  hTrkNfitVsPt    -> Write();
  hTrkDcaVsPt     -> Write();
  hTrkPtVsDf[0]   -> Write();
  hTrkPtVsDf[1]   -> Write();
  hTrkDfVsEta[0]  -> Write();
  hTrkDfVsEta[1]  -> Write();
  pTrkNfitVsPt    -> Write();
  pTrkDcaVsPt     -> Write();
  fOutput         -> Close();
  fInput          -> cd();
  fInput          -> Close();
  cout << "  Made QA plots!\n" << endl;

}

// End ------------------------------------------------------------------------
