Int_t nfiles = 1700; // Can reduce to test if the script runs
vector<TString> files;

TString inputDirDataParent = "/eos/cms/store/group/phys_heavyions/prdas/Run2_v0_wjet_trees_DataforAB/JetHT/";
TString inputDirData[] = {
  "crab_20250326_131351/250326_121419/0000/",
  "crab_20250326_132234/250326_122251/0000/",
  "crab_20250326_132655/250326_122710/0000/",
  "crab_20250326_132841/250326_122853/0000/",
  "crab_20250326_133049/250326_123102/0000/",
  "crab_20250326_135216/250326_125234/0000/",
  "crab_20250326_135300/250326_125314/0000/",
  "crab_20250326_135319/250326_125333/0000/",
  "crab_20250326_135341/250326_125353/0000/",
  "crab_20250326_141253/250326_131307/0000/",
  "crab_20250326_141313/250326_131324/0000/",
  "crab_20250326_141329/250326_131341/0000/",
  "crab_20250326_141344/250326_131356/0000/",
  "crab_20250326_141405/250326_131419/0000/",
  "crab_20250326_141430/250326_131441/0000/",
  "crab_20250326_141450/250326_131502/0000/",
  "crab_20250326_141510/250326_131523/0000/"
};


Int_t CandNo[] = {0,1,2,3,4,5,6};
TString cand_tree[] = {"Lambda","Kshort","Xi","Omega","AntiLambda","AntiXi","AntiOmega"};
TString cand_tree_ana[] = {"lambdaana","kshortana","xiana","omegaana","antilambdaana","antixiana","antiomegaana"};
Int_t cand_ID[] = {3122,310,3312,3334,3122,3312,3334};

TString OutputFileName = "TestOutput_v0InJets_forAB";

Double_t JetR = 0.8;

//Int_t GetCandNo(TString candName){CandNo = ;return CandNo;};

//Int_t LambdaId = 3122;
//Int_t KshortId = 310;
//Double_t MinLambdaRecoMass = 1.;
//Bool_t AnaLambda = kTRUE;
//Bool_t AnaKshort = kFALSE;

//void SetAnaLambda(Bool_t c = kFALSE){AnaLambda=c;}
//void SetAnaKshort(Bool_t c = kFALSE){AnaKshort=c;}

TH1F *hzVtxL=  new TH1F("hzVtxL","hzVtxL",150,-30,30);
TH1F *hLMass_in=  new TH1F("hLMass_in","hLMass_in",4000,0.,2.);
TH1F *hpdgIdL_in = new TH1F("hpdgIdL_in","hpdgIdL_in",3500,0.,3500.);

TH1F *hjetEta=  new TH1F("hjetEta","hjetEta",1000,-10,10);
TH1F *hjetPt=  new TH1F("hjetPt","hjetPt",100,0,2500);
TH1F *hjetPhi=  new TH1F("hjetPhi","hjetPhi",1400,-7,7);
TH1F *hjetNumDaughters=  new TH1F("hjetNumDaughters","hjetNumDaughters",500,-0.5,499.5);
TH1F *hchargedMult =  new TH1F("hchargedMult","hchargedMult",500,-0.5,499.5);

TH1F *hjetEta_wV0=  new TH1F("hjetEta_wV0","hjetEta_wV0",1000,-10,10);
TH1F *hjetPt_wV0=  new TH1F("hjetPt_wV0","hjetPt_wV0",100,0,2500);
TH1F *hjetPhi_wV0=  new TH1F("hjetPhi_wV0","hjetPhi_wV0",1400,-7,7);
TH1F *hchargedMult_wV0 =  new TH1F("hchargedMult_wV0","hchargedMult_wV0",500,-0.5,499.5);

TH2F *h2chargedMultCandMinus =  new TH2F("h2chargedMultCandMinus","h2chargedMultCandMinus",500,-0.5,499.5,10000,0.,1.);
TH2F *h2chargedMultCand =  new TH2F("h2chargedMultCand","h2chargedMultCand",500,-0.5,499.5,10000,0.,1.);


TH1F *hPtL_in=  new TH1F("hPtL_in","hPtL_in",1500,0.,1500.);
TH1F *hEtaL_in=  new TH1F("hEtaL_in","hEtaL_in",1000,-10,10);
TH1F *hPhiL_in=  new TH1F("hPhiL_in","hPhiL_in",1400,-7,7);

TH3F *h3LMass_in=  new TH3F("h3LMass_in","h3LMass_in",500,-0.5,499.5,400,0.,2000.,400,0.,2.);

