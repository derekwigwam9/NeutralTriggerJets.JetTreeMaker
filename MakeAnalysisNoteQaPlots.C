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
  TH1D *hEvtVz;
  TH1D *hTrgEt;
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
  const UInt_t  nDca(30);
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
  const Float_t phi[2]  = {0., 6.3};
  const Float_t tsp[2]  = {0., 1.};
  const Float_t fit[2]  = {0., 50.};
  const Float_t rat[2]  = {0., 1.};
  const Float_t dca[2]  = {0., 3.};
  const Float_t pt[2]   = {-5., 50.};
  const Float_t df[2]   = {-1.6, 4.7};

  hEvtNum         = new TH1D("hEvtNum", "", nNum, num[0], num[1]);
  hEvtVz          = new TH1D("hEvtVz", "", nVz, vz[0], vz[1]);
  hEvtVr          = new TH1D("hEvtVr", "", nVr, vr[0], vr[1]);
  hEvtPrim        = new TH1D("hEvtPrim", "", nPrim, prim[0], prim[1]);
  hTrgEt          = new TH1D("hTrgEt", "", nEt, et[0], et[1]);
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
  hTrkPtVsDf[0]   = new TH2D("hTrkPtVsDfPi0", "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
  hTrkPtVsDf[1]   = new TH2D("hTrkPtVsDfGam", "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
  hTrkDfVsEta[0]  = new TH2D("hTrkDfVsEtaPi0", "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
  hTrkDfVsEta[1]  = new TH2D("hTrkDfVsEtaGam", "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
  hEvtVz          -> Sumw2();
  hEvtVr          -> Sumw2();
  hEvtPrim        -> Sumw2();
  hTrgEt          -> Sumw2();
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
    hEvtVz     -> Fill(zVtx);
    hEvtVr     -> Fill(rVtx);
    hEvtPrim   -> Fill(nTrks);
    hEvtVyVsVx -> Fill(xVtx, yVtx);
    hEvtVxVsVz -> Fill(zVtx, xVtx);
    hEvtVyVsVz -> Fill(zVtx, yVtx);

    // vertex cuts
    const Bool_t isInRcut = (TMath::Abs(rVtx) < rVtxMax);
    const Bool_t isInZcut = (TMath::Abs(zVtx) < zVtxMax);
    if (isInRcut && isInZcut)
      nTrgCut[2]++;
    else
      continue;


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
    hTrgEt -> Fill(eTtrg);
    if (isInTspCut) {
      hTrgEta      -> Fill(hPhys);
      hTrgPhi      -> Fill(fTrg);
      hTrgEtVsEta  -> Fill(hPhys, eTtrg);
      hTrgEtVsPhi  -> Fill(fTrg, eTtrg);
      hTrgEtaVsPhi -> Fill(fTrg, hPhys);
    }
    if (isInEtCut)
      hTrgTsp[2] -> Fill(tspTrg);

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
      hTrkNfit[0] -> Fill(nFitTrk);
      if (!isInFitCut)
        hTrkNfit[1] -> Fill(nFitTrk);

      // nFit / nPoss cut
      hTrkRfit[0] -> Fill(rFitTrk);
      if (!isInRatioCut)
        hTrkRfit[1] -> Fill(rFitTrk);

      // dca cut
      hTrkDca[0] -> Fill(dcaTrk);
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
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fBar(0.6);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.3);
  const Float_t fOffL(0.7);
  const Float_t fOffB(0.2);
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
  hEvtNum -> SetLineColor(923);
  hEvtNum -> SetLineStyle(1);
  hEvtNum -> SetFillColor(923);
  hEvtNum -> SetFillStyle(1001);
  hEvtNum -> SetMarkerColor(923);
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
  hEvtNum -> GetXaxis() -> SetLabelOffset(fLab);
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


  // draw plots
  const UInt_t width(750);
  const UInt_t height(750);
  const UInt_t fMode(0);
  const UInt_t fBord(2);
  const UInt_t fGrid(0);
  const UInt_t fFrame(0);
  const UInt_t fLogX(0);
  const UInt_t fLogY(1);
  const UInt_t fTick(1);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.05);
  const Float_t fMarginT(0.05);
  const Float_t fMarginB(0.15);

  TCanvas *cEvtNum = new TCanvas("cEvtNum", "", width, height);
  cEvtNum -> SetLogx(fLogX);
  cEvtNum -> SetLogy(fLogY);
  cEvtNum -> SetGrid(fGrid, fGrid);
  cEvtNum -> SetTicks(fTick, fTick);
  cEvtNum -> SetBorderMode(fMode);
  cEvtNum -> SetBorderSize(fBord);
  cEvtNum -> SetFrameBorderMode(fFrame);
  cEvtNum -> SetLeftMargin(fMarginL);
  cEvtNum -> SetTopMargin(fMarginT);
  cEvtNum -> SetRightMargin(fMarginR);
  cEvtNum -> SetBottomMargin(fMarginB);
  cEvtNum -> cd();
  hEvtNum -> Draw("B");
  fOutput -> cd();
  cEvtNum -> Write();
  cEvtNum -> Close();

  // close files
  fOutput         -> cd();
  hEvtNum         -> Write();
  hEvtVz          -> Write();
  hEvtVr          -> Write();
  hEvtPrim        -> Write();
  hTrgEt          -> Write();
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
  hTrkPtVsDf[0]   -> Write();
  hTrkPtVsDf[1]   -> Write();
  hTrkDfVsEta[0]  -> Write();
  hTrkDfVsEta[1]  -> Write();
  fOutput         -> Close();
  fInput          -> cd();
  fInput          -> Close();
  cout << "  Made QA plots!\n" << endl;

}

// End ------------------------------------------------------------------------
