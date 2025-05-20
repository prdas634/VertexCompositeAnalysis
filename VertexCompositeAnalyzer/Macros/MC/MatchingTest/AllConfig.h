Int_t nfiles = 20; // Can reduce to test if the script runs
vector<TString> files;

TString inputDirDataParent = "/eos/cms/store/group/phys_heavyions/prdas/";
TString inputDirData[] = {"Run2_v0_mc_wjet_trees_2/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/crab_20250302_162125/250302_152307/0000/"};

//Run2_v0_mc_wjet_trees_1/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/crab_20250103_154353/250103_144426/0000/"};

//TString inputDirDataParent = "/afs/cern.ch/work/p/prdas/";
//TString inputDirData[] = {"private/TestV0Reco/CMSSW_10_6_4_patch1/src/VertexCompositeAnalysis/VertexCompositeProducer/test/OutRootFiles/"};

TString OutputFileName = "Output_V0InJets_Matching_470to600";
//TString OutputTreeFileName = "TMVA_sig_test_new_v7_470to600_s3";

Int_t CandNo[] = {0,1,2,3,4,5,6};
TString cand_tree[] = {"Lambda","AntiLambda","Kshort","Xi","AntiXi","Omega","AntiOmega"};
TString cand_tree_ana[] = {"lambdaana","antilambdaana","kshortana","xiana","antixiana","omegaana","antiomegaana"};
Int_t cand_ID[] = {3122,3122,310,3312,3312,3334,3334};

Double_t JetR = 0.8;
//Int_t LambdaId = 3122;
//Int_t KshortId = 310;
//Double_t MinLambdaRecoMass = 1.;
//Bool_t AnaLambda = kTRUE;
//Bool_t AnaKshort = kFALSE;

//void SetAnaLambda(Bool_t c = kFALSE){AnaLambda=c;}
//void SetAnaKshort(Bool_t c = kFALSE){AnaKshort=c;}

void displayProgress(long current, long max);
void GetFiles(char const* input, vector<TString>& files, int file_limit);
Double_t dR(Double_t RefEta, Double_t RefPhi, Double_t MEta, Double_t MPhi);
Double_t dPhi(Double_t RefPhi, Double_t MPhi);
Bool_t IsSatisfyKineCuts(Double_t RefPt, Double_t RefEta);



/*Double_t DelR_RG[500][500] = {0.};

for(Int_t indx = 0; indx < 500; indx++){
    for(Int_t indy = 0; indy < 500; indy++){
        DelR_RG[indx][indy] = 999.;
    }
}*/
Double_t refdistAll = 100.;
Int_t refIndRAll = -9999;
Int_t refIndGAll = -9999;
Double_t refdist = 0.05; //decided from Dr_RG distribution
Double_t reffracpT = 0.2;
Int_t refIndR = -999;
Int_t refIndG = -999;
Double_t Dr_RG = 1000.;
Double_t FracPt_RG = 1.;
Double_t RefRMJetEta = -999.;
Double_t RefRMJetPhi = -999.;



TH1F *h1decay=  new TH1F("h1decay","h1decay",400,-200.,200.);
TH1F *h1pseudecay=  new TH1F("h1pseudecay","h1pseudecay",400,-200.,200.);

TH1F *hpdgIdL_in = new TH1F("hpdgIdL_in","hpdgIdL_in",3500,0.,3500.);

TH2F *h2LMassPdg_cand = new TH2F("h2LMassPdg_cand","h2LMassPdg_cand",2800,0.,1.4,3500,0.,3500.);
TH2F *h2LMass_matched = new TH2F("h2LMass_matched","h2LMass_matched",2800,0.,1.4,2800,0.,1.4);
TH1F *hLMass_unmatched =  new TH1F("hLMass_unmatched","hLMass_unmatched",2800,0.,1.4);
TH1F *hzVtxL=  new TH1F("hzVtxL","hzVtxL",150,-30,30);
TH1F *hLMass=  new TH1F("hLMass","hLMass",800,1.,1.4);
TH1F *hLMass_cand=  new TH1F("hLMass_cand","hLMass_cand",800,1.,1.4);
TH1F *hLMass_in=  new TH1F("hLMass_in","hLMass_in",800,1.,1.4);
TH1F *hLMass_in_PtLead=  new TH1F("hLMass_in_PtLead","hLMass_in_PtLead",800,1.,1.4);
TH1F *hLMass_out=  new TH1F("hLMass_out","hLMass_out",800,1.,1.4);

TH1F *hzVtxK=  new TH1F("hzVtxK","hzVtxK",150,-30,30);
TH1F *hKMass=  new TH1F("hKMass","hKMass",400,0.4,0.6);

TH1F *hjetN=  new TH1F("hjetN","hjetN",100,0.,100.);
TH1F *hjetEta=  new TH1F("hjetEta","hjetEta",1000,-10,10);
TH1F *hjetPt=  new TH1F("hjetPt","hjetPt",100,0,2500);
TH1F *hjetPtLead=  new TH1F("hjetPtLead","hjetPtLead",100,0,2500);
TH1F *hjetPhi=  new TH1F("hjetPhi","hjetPhi",1400,-7,7);

TH1F *hjetEta_wV0=  new TH1F("hjetEta_wV0","hjetEta_wV0",1000,-10,10);
TH1F *hjetPt_wV0=  new TH1F("hjetPt_wV0","hjetPt_wV0",100,0,2500);
TH1F *hjetPhi_wV0=  new TH1F("hjetPhi_wV0","hjetPhi_wV0",1400,-7,7);

TH1F *hjetEta_wV0_PtLead=  new TH1F("hjetEta_wV0_PtLead","hjetEta_wV0_PtLead",1000,-10,10);
TH1F *hjetPt_wV0_PtLead=  new TH1F("hjetPt_wV0_PtLead","hjetPt_wV0_PtLead",100,0,2500);
TH1F *hjetPhi_wV0_PtLead=  new TH1F("hjetPhi_wV0_PtLead","hjetPhi_wV0_PtLead",1400,-7,7);

TH1F *hjetNumDaughters=  new TH1F("hjetNumDaughters","hjetNumDaughters",500,-0.5,499.5);
TH1F *hchargedMult =  new TH1F("hchargedMult","hchargedMult",500,-0.5,499.5);
TH1F *hchargedMult_wV0 =  new TH1F("hchargedMult_wV0","hchargedMult_wV0",500,-0.5,499.5);
TH1F *hchargedMult_wV0_PtLead =  new TH1F("hchargedMult_wV0_PtLead","hchargedMult_wV0_PtLead",500,-0.5,499.5);

TH2F *h2chargedMultCandMinus2 =  new TH2F("h2chargedMultCandMinus2","h2chargedMultCandMinus2",500,-0.5,499.5,10000,0.,1.);
TH2F *h2chargedMultCandMinus =  new TH2F("h2chargedMultCandMinus","h2chargedMultCandMinus",500,-0.5,499.5,10000,0.,1.);
TH2F *h2chargedMultCand =  new TH2F("h2chargedMultCand","h2chargedMultCand",500,-0.5,499.5,10000,0.,1.);
//TH1F *hDelRV0Jet_in=  new TH1F("hDelRV0Jet_in","hDelRV0Jet_in",1100,-0.1,1.);
//TH1F *hDelRV0Jet_in_PtLead=  new TH1F("hDelRV0Jet_in_PtLead","hDelRV0Jet_in_PtLead",1100,-0.1,1.);
//TH3F *h3DelRV0Jet_in=  new TH3F("h3DelRV0Jet_in","h3DelRV0Jet_in",100,0,2500,1000,0.,1000.,1100,-0.1,1.);

//TH1F *hDelRV0Jet_out=  new TH1F("hDelRV0Jet_out","hDelRV0Jet_out",1100,-0.1,1.);
//TH2F *h2DelRV0Jet_out=  new TH2F("h2DelRV0Jet_out","h2DelRV0Jet_out",1000,0.,1000.,1100,-0.1,1.);

TH1F *hPtL=  new TH1F("hPtL","hPtL",1500,0.,1500.);
TH1F *hPtL_in=  new TH1F("hPtL_in","hPtL_in",1500,0.,1500.);
TH1F *hPtL_in_PtLead=  new TH1F("hPtL_in_PtLead","hPtL_in_PtLead",1500,0.,1500.);
TH1F *hPtL_out=  new TH1F("hPtL_out","hPtL_out",1500,0.,1500.);



TH1F *hEtaL=  new TH1F("hEtaL","hEtaL",1000,-10,10);
TH1F *hPhiL=  new TH1F("hPhiL","hPhiL",1400,-7,7);
TH1F *hEtaL_in=  new TH1F("hEtaL_in","hEtaL_in",1000,-10,10);
TH1F *hEtaL_in_PtLead=  new TH1F("hEtaL_in_PtLead","hEtaL_in_PtLead",1000,-10,10);
TH1F *hPhiL_in=  new TH1F("hPhiL_in","hPhiL_in",1400,-7,7);
TH1F *hPhiL_in_PtLead=  new TH1F("hPhiL_in_PtLead","hPhiL_in_PtLead",1400,-7,7);
TH1F *hEtaL_out=  new TH1F("hEtaL_out","hEtaL_out",1000,-10,10);
TH1F *hPhiL_out=  new TH1F("hPhiL_out","hPhiL_out",1400,-7,7);

//TH1F *hInsideJet=  new TH1F("hInsideJet","hInsideJet",20,0.,20.);
//TH1F *hDelRV0MJet=  new TH1F("hDelRV0MJet","hDelRV0MJet",1100,-0.1,1.);

TH3F *h3LMass_in=  new TH3F("h3LMass_in","h3LMass_in",500,-0.5,499.5,500,0.,2000.,800,1.,1.4);
TH3F *h3LMass_in_PtLead=  new TH3F("h3LMass_in_PtLead","h3LMass_in_PtLead",500,-0.5,499.5,500,0.,2000.,800,1.,1.4);


TH1F *h1DelR_RG = new TH1F("h1DelR_RG","h1DelR_RG",1100,-0.1,1.);
TH1F *h1FracPt_RG = new TH1F("h1FracPt_RG","h1FracPt_RG",20000,-10.,10.);

TH1F *h1DelR_RG_matched_inJ = new TH1F("h1DelR_RG_matched_inJ","h1DelR_RG_matched_inJ",1100,-0.1,1.);
TH1F *h1FracPt_RG_matched_inJ = new TH1F("h1FracPt_RG_matched_inJ","h1FracPt_RG_matched_inJ",20000,-10.,10.);

TH1F *h1DelR_RG_matched_wJ = new TH1F("h1DelR_RG_matched_wJ","h1DelR_RG_matched_wJ",1100,-0.1,1.);
TH1F *h1FracPt_RG_matched_wJ = new TH1F("h1FracPt_RG_matched_wJ","h1FracPt_RG_matched_wJ",20000,-10.,10.);

TH1F *h1DelR_RG_matched_woJ = new TH1F("h1DelR_RG_matched_woJ","h1DelR_RG_matched_woJ",1100,-0.1,1.);
TH1F *h1FracPt_RG_matched_woJ = new TH1F("h1FracPt_RG_matched_woJ","h1FracPt_RG_matched_woJ",20000,-10.,10.);

TH1F *h1DelR_RG_closest_woJ = new TH1F("h1DelR_RG_closest_woJ","h1DelR_RG_closest_woJ",1100,-0.1,1.);
TH1F *h1FracPt_RG_closest_woJ = new TH1F("h1FracPt_RG_closest_woJ","h1FracPt_RG_closest_woJ",20000,-10.,10.);

TH1F *h1DelR_RG_closest_wJ = new TH1F("h1DelR_RG_closest_wJ","h1DelR_RG_closest_wJ",1100,-0.1,1.);
TH1F *h1FracPt_RG_closest_wJ = new TH1F("h1FracPt_RG_closest_wJ","h1FracPt_RG_closest_wJ",20000,-10.,10.);

TH1F *h1DelPt_RG = new TH1F("h1DelPt_RG","h1DelPt_RG",4000,-1000.,1000.);
TH1F *h1pdgId_G = new TH1F("h1pdgId_G","h1pdgId_G",7000,-3500.,3500.);

TH1F *h1FracPt_RG_wcut = new TH1F("h1FracPt_RG_wcut","h1FracPt_RG_wcut",20000,-10.,10.);
TH1F *h1DelPt_RG_wcut = new TH1F("h1DelPt_RG_wcut","h1DelPt_RG_wcut",4000,-1000.,1000.);

TH1F *hPtL_RA =  new TH1F("hPtL_RA","hPtL_RA",1500,0.,1500.);
TH1F *hPtL_GA_wojet =  new TH1F("hPtL_GA_wojet","hPtL_GA_wojet",1500,0.,1500.);

TH1F *hPtL_GA =  new TH1F("hPtL_GA","hPtL_GA",1500,0.,1500.);
TH1F *hPtL_GA_in =  new TH1F("hPtL_GA_in","hPtL_GA_in",1500,0.,1500.);

TH1F *hPtL_GA_in_wcand =  new TH1F("hPtL_GA_in_wcand","hPtL_GA_in_wcand",1500,0.,1500.);
TH1F *hMassL_GA_in =  new TH1F("hMassL_GA_in","hMassL_GA_in",2800,0.,1.4);

TH1F *hPtL_GA_in_NchLT20 =  new TH1F("hPtL_GA_in_NchLT20","hPtL_GA_in_NchLT20",1500,0.,1500.);
TH1F *hPtL_GA_in_NchGT80 =  new TH1F("hPtL_GA_in_NchGT80","hPtL_GA_in_NchGT80",1500,0.,1500.);
TH2F *h2PtL_GA_in =  new TH2F("h2PtL_GA_in","h2PtL_GA_in",500,-0.5,499.5,1500,0.,1500.);

TH1F *hPtL_GM_in =  new TH1F("hPtL_GM_in","hPtL_GM_in",1500,0.,1500.);

TH1F *h1Pt_GA =  new TH1F("h1Pt_GA","h1Pt_GA",1500,0.,1500.);
TH1F *h1Pt_RA =  new TH1F("h1Pt_RA","h1Pt_RA",1500,0.,1500.);
TH1F *h1Mass_GA =  new TH1F("h1Mass_GA","h1Mass_GA",2800,0.,1.4);
TH1F *h1Mass_RA =  new TH1F("h1Mass_RA","h1Mass_RA",2800,0.,1.4);

TH1F *h1Pt_GM_woJ =  new TH1F("h1Pt_GM_woJ","h1Pt_GM_woJ",1500,0.,1500.);

TH1F *hMassL_GM_in =  new TH1F("hMassL_GM_in","hMassL_GM_in",2800,0.,1.4);
TH1F *hDelR_GM_in_wGJ = new TH1F("hDelR_GM_in_wGJ","hDelR_GM_in_wGJ",1100,-1.,10.);
TH1F *hMassL_GM_in_wGJ =  new TH1F("hMassL_GM_in_wGJ","hMassL_GM_in_wGJ",2800,0.,1.4);
TH1F *hDelR_GM_in_wGJ_wMV0 = new TH1F("hDelR_GM_in_wGJ_wMV0","hDelR_GM_in_wGJ_wMV0",1100,-1.,10.);

TH1F *hDelR_RGM_wMV0 = new TH1F("hDelR_RGM_wMV0","hDelR_RGM_wMV0",1100,-1.,10.);
TH1F *hDelR_RGJet = new TH1F("hDelR_RGJet","hDelR_RGJet",1100,-1.,10.);
TH2F *h2GJetPtDelR =  new TH2F("h2GJetPtDelR","h2GJetPtDelR",1100,-1.,10.,1500,0.,1500.);

TH1F *hPtL_GM_in_wGJ =  new TH1F("hPtL_GM_in_wGJ","hPtL_GM_in_wGJ",1500,0.,1500.);
TH1F *hPtL_GM_in_wGJ_NchLT20 =  new TH1F("hPtL_GM_in_wGJ_NchLT20","hPtL_GM_in_wGJ_NchLT20",1500,0.,1500.);
TH1F *hPtL_GM_in_wGJ_NchGT80 =  new TH1F("hPtL_GM_in_wGJ_NchGT80","hPtL_GM_in_wGJ_NchGT80",1500,0.,1500.);
TH2F *h2PtL_GM_in_wGJ =  new TH2F("h2PtL_GM_in_wGJ","h2PtL_GM_in_wGJ",500,-0.5,499.5,1500,0.,1500.);


TH1F *hCount_noGenJetwMV0 =  new TH1F("hCount_noGenJetwMV0","hCount_noGenJetwMV0",4,0.,4.);
TH1F *hCount_GenJetwMV0 =  new TH1F("hCount_GenJetwMV0","hCount_GenJetwMV0",4,0.,4.);
TH1F *hCount_noGenJet =  new TH1F("hCount_noGenJet","hCount_noGenJet",4,0.,4.);

TH1F *hDelPhi_RGJet = new TH1F("hDelPhi_RGJet","hDelPhi_RGJet",1400,-7.,7.);
TH1F *hDelEta_RGJet = new TH1F("hDelEta_RGJet","hDelEta_RGJet",1000,-5.,5.);

TH1F *h1GJetEta=  new TH1F("h1GJetEta","h1GJetEta",1000,-10,10);
TH1F *h1GJetPhi=  new TH1F("h1GJetPhi","h1GJetPhi",1400,-7,7);
TH1F *h1RJetEta=  new TH1F("h1RJetEta","h1RJetEta",1000,-10,10);
TH1F *h1RJetPhi=  new TH1F("h1RJetPhi","h1RJetPhi",1400,-7,7);

TH1F *h1Ncand_GA =  new TH1F("h1Ncand_GA","h1Ncand_GA",100,0.,100.);
TH1F *h1Ncand_RA =  new TH1F("h1Ncand_RA","h1Ncand_RA",100,0.,100.);
