#include "AllConfig.h"

void displayProgress(long current, long max){
  int width = 82; // Hope the terminal is at least that wide.
  int barWidth = width - 2;
  cout << "\x1B[2K";    // Clear line
  cout << "\x1B[2000D"; // Cursor left
  cout << '[';
  for (int i = 0; i < barWidth; ++i)
    {
      if (i < barWidth * current / max)
	{
	  cout << '=';
	}
      else
	{
	  cout << ' ';
	}
    }
  cout << ']';
  cout << " " << Form("%d/%d (%5.2f%%)", (int)current, (int)max, 100.0 * current / max);
  cout.flush();
}

void GetFiles(char const* input, vector<TString>& files, int file_limit) {
  TSystemDirectory dir(input, input);
  TList *list = dir.GetListOfFiles();
  
  if (list) {
    TSystemFile *file;
    string fname;
    TIter next(list);
    while ((file = (TSystemFile*) next())) {
      fname = file->GetName();
      if (file->IsDirectory() && (fname.find(".") == std::string::npos)) {
	string newDir = std::string(input) + fname + "/";
	GetFiles(newDir.c_str(), files, file_limit);
      }
      else if ((fname.find(".root") != std::string::npos)) {
	files.push_back(std::string(input) + fname);
	// cout << files.back() << endl;
	if((int)files.size()>file_limit) return;
      }
    }
  }
  return;
}

Double_t dR(Double_t RefEta, Double_t RefPhi, Double_t MEta, Double_t MPhi){
  Double_t DelEta = TMath::Abs(MEta-RefEta);
  Double_t DelPhi = TMath::Abs(MPhi-RefPhi);
  if(DelPhi>TMath::Pi())DelPhi=TMath::Abs(DelPhi - TMath::TwoPi());
  Double_t DelR = TMath::Sqrt((DelEta*DelEta)+(DelPhi*DelPhi));
  return DelR;
}

Double_t dPhi(Double_t RefPhi, Double_t MPhi){
  Double_t DelPhi = TMath::Abs(MPhi-RefPhi);
  if(DelPhi>TMath::Pi())DelPhi=TMath::Abs(DelPhi - TMath::TwoPi());
  return DelPhi;
}

Bool_t IsSatifyKineCuts(Double_t RefPt, Double_t RefEta){
  if(RefPt < 2.0 || TMath::Abs(RefEta) > 2.4)return kFALSE;
  return kTRUE;
}

void GetObs_v7(Int_t nIpDir = 0, Bool_t DoAnaLambda = kFALSE, Bool_t DoAnaKshort = kFALSE, Bool_t IsSave = kFALSE){
  
  files.clear();
  
  GetFiles((inputDirDataParent+inputDirData[nIpDir]).Data(),files,nfiles);
  cout<<"Got "<<files.size()<< " files" << endl;
  //for(Int_t i = 0; i <= nfiles+1; i++)cout << files[i] << endl;
  
  TChain JetTree("JetTree"), LambdaTree("LambdaTree"), KshortTree("KshortTree");
  for(TString temp_path:files){
    JetTree.Add(temp_path+"/analyzerOffline/trackTree");
    LambdaTree.Add(temp_path+"/lambdaana/ParticleTree");
    KshortTree.Add(temp_path+"/kshortana/ParticleTree");
  }
  
  SetAnaLambda(DoAnaLambda);
  SetAnaKshort(DoAnaKshort);
  
  if(AnaLambda && !AnaKshort){
    Int_t nEntL = LambdaTree.GetEntries();
    cout << "\n Number of Entries L : " << nEntL << endl;
    
    Int_t nEntJ = JetTree.GetEntries();
    cout << "\n Number of Entries Jet : " << nEntJ << endl;
    
    Float_t VtxZL = 0;
    LambdaTree.SetBranchAddress("bestvtxZ",&VtxZL);
    
    Int_t jetN;
    Float_t jetPtLead;
    JetTree.SetBranchAddress("jetN",&jetN);
    JetTree.SetBranchAddress("jetPtLead",&jetPtLead);
    
    
    vector<float> *cand_massL = 0;
    vector<char> *cand_chargeL = 0;
    vector<unsigned char> *cand_statusL = 0;
    vector<float> *cand_yL = 0;
    vector<float> *cand_etaL = 0;
    vector<float> *cand_phiL = 0;
    vector<float> *cand_pTL = 0;
    vector<float> *cand_pL = 0;
    vector<int> *cand_pdgIdL = 0;
    
    vector<float> *gen_massL = 0;
    vector<float> *gen_etaL = 0;
    vector<float> *gen_phiL = 0;
    vector<float> *gen_pTL = 0;
    vector<int> *gen_pdgIdL = 0;
    
    vector<float> *cand_vtxProbL = 0;
    vector<float> *cand_vtxChi2L = 0;
    vector<float> *cand_angle2DL = 0;
    vector<float> *cand_angle3DL = 0;
    vector<float> *cand_dcaL = 0;
    vector<float> *cand_decayLength3DL = 0;
    vector<float> *cand_decayLength2DL = 0;
    vector<float> *cand_decayLengthError3DL = 0;
    vector<float> *cand_decayLengthError2DL = 0;
    vector<float> *cand_pseudoDecayLengthError3DL = 0;
    vector<float> *cand_pseudoDecayLengthError2DL = 0;
    vector<unsigned int> *cand_trkIdxL = 0;
    vector<bool> *trk_isHPL = 0;
    vector<float> *trk_nHitL = 0;
    vector<float> *trk_nChi2L = 0;
    vector<float> *trk_pTErrL = 0;
    vector<float> *trk_xyDCASignificanceL = 0;
    vector<float> *trk_zDCASignificanceL = 0;
    
    //vector<vector<unsigned int>> *cand_dauIdxL = 0;
    //vector<vector<unsigned int>> *cand_momIdxL = 0;
    
    //TBranch *L_dauIdx = 0;
    //LambdaTree.SetBranchAddress("cand_dauIdx",&cand_dauIdxL,&L_dauIdx);
    //LambdaTree.SetBranchAddress("cand_dauIdx",&cand_dauIdxL);
    //TBranch *L_momIdx = 0;
    //LambdaTree.SetBranchAddress("cand_momIdx",&cand_momIdxL,&L_momIdx);
    
    
    //vector<float> *trk_zDCASignificance = 0;
    //vector<float> *trk_xyDCASignificance = 0;
    
    
    vector<float> *jetEta = 0;
    vector<float> *jetPt = 0;
    vector<float> *jetPhi = 0;
    vector<int> *jetNumDaughters = 0;
    vector<int> *chMult = 0;
    vector<int> *genChMult = 0;
    
    vector<float> *genJetEta = 0;
    vector<float> *genJetPt = 0;
    vector<float> *genJetPhi = 0;

    
    TBranch *LMass = 0;
    LambdaTree.SetBranchAddress("cand_mass",&cand_massL,&LMass);
    TBranch *LCharge = 0;
    LambdaTree.SetBranchAddress("cand_charge",&cand_chargeL,&LCharge);
    TBranch *LStatus = 0;
    LambdaTree.SetBranchAddress("cand_status",&cand_statusL,&LStatus);
    TBranch *Ly = 0;
    LambdaTree.SetBranchAddress("cand_y",&cand_yL,&Ly);
    TBranch *LEta = 0;
    LambdaTree.SetBranchAddress("cand_eta",&cand_etaL,&LEta);
    TBranch *LPhi = 0;
    LambdaTree.SetBranchAddress("cand_phi",&cand_phiL,&LPhi);
    TBranch *LP = 0;
    LambdaTree.SetBranchAddress("cand_p",&cand_pL,&LP);
    TBranch *LPt = 0;
    LambdaTree.SetBranchAddress("cand_pT",&cand_pTL,&LPt);
    TBranch *LPdgId = 0;
    LambdaTree.SetBranchAddress("cand_pdgId",&cand_pdgIdL,&LPdgId);
    
    TBranch *LMass_G = 0;
    LambdaTree.SetBranchAddress("gen_mass",&gen_massL,&LMass_G);
    TBranch *LEta_G = 0;
    LambdaTree.SetBranchAddress("gen_eta",&gen_etaL,&LEta_G);
    TBranch *LPhi_G = 0;
    LambdaTree.SetBranchAddress("gen_phi",&gen_phiL,&LPhi_G);
    TBranch *LPt_G = 0;
    LambdaTree.SetBranchAddress("gen_pT",&gen_pTL,&LPt_G);
    TBranch *LPdgId_G = 0;
    LambdaTree.SetBranchAddress("gen_pdgId",&gen_pdgIdL,&LPdgId_G);
    
    TBranch *LvtxProb = 0;
    LambdaTree.SetBranchAddress("cand_vtxProb",&cand_vtxProbL,&LvtxProb);
    TBranch *LvtxChi2 = 0;
    LambdaTree.SetBranchAddress("cand_vtxChi2",&cand_vtxChi2L,&LvtxChi2);
    TBranch *Langle3D = 0;
    LambdaTree.SetBranchAddress("cand_angle3D",&cand_angle3DL,&Langle3D);
    TBranch *Langle2D = 0;
    LambdaTree.SetBranchAddress("cand_angle2D",&cand_angle2DL,&Langle2D);
    TBranch *Ldca = 0;
    LambdaTree.SetBranchAddress("cand_dca",&cand_dcaL,&Ldca);
    
    TBranch *LdecayLength3D = 0;
    LambdaTree.SetBranchAddress("cand_decayLength3D",&cand_decayLength3DL,&LdecayLength3D);
    TBranch *LdecayLength2D = 0;
    LambdaTree.SetBranchAddress("cand_decayLength2D",&cand_decayLength2DL,&LdecayLength2D);
    TBranch *LdecayLengthError3D = 0;
    LambdaTree.SetBranchAddress("cand_decayLengthError3D",&cand_decayLengthError3DL,&LdecayLengthError3D);
    TBranch *LdecayLengthError2D = 0;
    LambdaTree.SetBranchAddress("cand_decayLengthError2D",&cand_decayLengthError2DL,&LdecayLengthError2D);
    TBranch *LpseudoDecayLengthError3D = 0;
    LambdaTree.SetBranchAddress("cand_pseudoDecayLengthError3D",&cand_pseudoDecayLengthError3DL,&LpseudoDecayLengthError3D);
    TBranch *LpseudoDecayLengthError2D = 0;
    LambdaTree.SetBranchAddress("cand_pseudoDecayLengthError2D",&cand_pseudoDecayLengthError2DL,&LpseudoDecayLengthError2D);
    
    //TBranch *LtrkIdx = 0;
    //LambdaTree.SetBranchAddress("cand_trkIdx",&cand_trkIdxL,&LtrkIdx);
    
    TBranch *Ltrk_isHP = 0;
    LambdaTree.SetBranchAddress("trk_isHP",&trk_isHPL,&Ltrk_isHP);
    TBranch *Ltrk_nHit = 0;
    LambdaTree.SetBranchAddress("trk_nHit",&trk_nHitL,&Ltrk_nHit);
    TBranch *Ltrk_nChi2 = 0;
    LambdaTree.SetBranchAddress("trk_nChi2",&trk_nChi2L,&Ltrk_nChi2);
    TBranch *Ltrk_pTErr = 0;
    LambdaTree.SetBranchAddress("trk_pTErr",&trk_pTErrL,&Ltrk_pTErr);
    TBranch *Ltrk_xyDCASignificance = 0;
    LambdaTree.SetBranchAddress("trk_xyDCASignificance",&trk_xyDCASignificanceL,&Ltrk_xyDCASignificance);
    TBranch *Ltrk_zDCASignificance = 0;
    LambdaTree.SetBranchAddress("trk_zDCASignificance",&trk_zDCASignificanceL,&Ltrk_zDCASignificance);
    
    
    
    TBranch *jNumDaughters = 0;
    JetTree.SetBranchAddress("jetNumDaughters",&jetNumDaughters,&jNumDaughters);

    TBranch *chargedMult = 0;
    JetTree.SetBranchAddress("chargedMultiplicity",&chMult,&chargedMult);

    TBranch *jEta = 0;
    JetTree.SetBranchAddress("jetEta",&jetEta,&jEta);
    
    TBranch *jPt = 0;
    JetTree.SetBranchAddress("jetPt",&jetPt,&jPt);  
    
    TBranch *jPhi = 0;
    JetTree.SetBranchAddress("jetPhi",&jetPhi,&jPhi);

    TBranch *chargedMult_G = 0;
    JetTree.SetBranchAddress("genJetChargedMultiplicity",&genChMult,&chargedMult_G);

    TBranch *jEta_G = 0;
    JetTree.SetBranchAddress("genJetEta",&genJetEta,&jEta_G);
    
    TBranch *jPt_G = 0;
    JetTree.SetBranchAddress("genJetPt",&genJetPt,&jPt_G);

    TBranch *jPhi_G = 0;
    JetTree.SetBranchAddress("genJetPhi",&genJetPhi,&jPhi_G);
    

    
    
    TFile *fOutTree = new TFile(Form("%s.root",OutputTreeFileName.Data()),"RECREATE");
    TTree *t_sig = new TTree("TMVATree_Sig","Signal tree for TMVA");
    TTree *t_bkg = new TTree("TMVATree_Bkg","Background tree for TMVA");
    Float_t mass_L;
    char charge_L; 
    unsigned char status_L;
    Float_t y_L, eta_L, phi_L, pT_L, p_L;
    Float_t vtxProb_L, vtxChi2_L, angle2D_L,angle3D_L, dca_L;
    Float_t decayLength2D_L, decayLength3D_L, decayLengthError2D_L, decayLengthError3D_L;
    Float_t pseudoDecayLengthError2D_L, pseudoDecayLengthError3D_L;
    Int_t pdgId_L, trkIdx_L;
    Float_t trk_isHP_L, trk_nHit_L, trk_nChi2_L, trk_pTErr_L, trk_xyDCASignificance_L, trk_zDCASignificance_L;
    
    Float_t ChMultInJets_L;
    
    t_sig->Branch("ChMultInJets",&ChMultInJets_L);
    
    t_sig->Branch("mass_cand",&mass_L);
    t_sig->Branch("charge_cand",&charge_L);
    t_sig->Branch("status_cand",&status_L);
    t_sig->Branch("y_cand",&y_L);
    t_sig->Branch("eta_cand",&eta_L);
    t_sig->Branch("phi_cand",&phi_L);
    t_sig->Branch("pT_cand",&pT_L);
    t_sig->Branch("p_cand",&p_L);
    t_sig->Branch("vtxProb_cand",&vtxProb_L);
    t_sig->Branch("vtxChi2_cand",&vtxChi2_L);
    t_sig->Branch("angle2D_cand",&angle2D_L);
    t_sig->Branch("angle3D_cand",&angle3D_L);
    t_sig->Branch("dca_cand",&dca_L);
    t_sig->Branch("decayLength2D_cand",&decayLength2D_L);
    t_sig->Branch("decayLength3D_cand",&decayLength3D_L);
    t_sig->Branch("decayLengthError2D_cand",&decayLengthError2D_L);
    t_sig->Branch("decayLengthError3D_cand",&decayLengthError3D_L);
    t_sig->Branch("pseudoDecayLengthError2D_cand",&pseudoDecayLengthError2D_L);
    t_sig->Branch("pseudoDecayLengthError3D_cand",&pseudoDecayLengthError3D_L);
    t_sig->Branch("pdgId_cand",&pdgId_L);
    //t_sig->Branch("trkIdx_cand",&trkIdx_L);
    t_sig->Branch("trk_isHP_cand",&trk_isHP_L);
    t_sig->Branch("trk_nHit_cand",&trk_nHit_L);
    t_sig->Branch("trk_nChi2_cand",&trk_nChi2_L);
    t_sig->Branch("trk_pTErr_cand",&trk_pTErr_L);
    t_sig->Branch("trk_xyDCASignificance_cand",&trk_xyDCASignificance_L);
    t_sig->Branch("trk_zDCASignificance_cand",&trk_zDCASignificance_L);
    
    t_bkg->Branch("ChMultInJets",&ChMultInJets_L);
    
    t_bkg->Branch("mass_cand",&mass_L);
    t_bkg->Branch("charge_cand",&charge_L);
    t_bkg->Branch("status_cand",&status_L);
    t_bkg->Branch("y_cand",&y_L);
    t_bkg->Branch("eta_cand",&eta_L);
    t_bkg->Branch("phi_cand",&phi_L);
    t_bkg->Branch("pT_cand",&pT_L);
    t_bkg->Branch("p_cand",&p_L);
    t_bkg->Branch("vtxProb_cand",&vtxProb_L);
    t_bkg->Branch("vtxChi2_cand",&vtxChi2_L);
    t_bkg->Branch("angle2D_cand",&angle2D_L);
    t_bkg->Branch("angle3D_cand",&angle3D_L);
    t_bkg->Branch("dca_cand",&dca_L);
    t_bkg->Branch("decayLength2D_cand",&decayLength2D_L);
    t_bkg->Branch("decayLength3D_cand",&decayLength3D_L);
    t_bkg->Branch("decayLengthError2D_cand",&decayLengthError2D_L);
    t_bkg->Branch("decayLengthError3D_cand",&decayLengthError3D_L);
    t_bkg->Branch("pseudoDecayLengthError2D_cand",&pseudoDecayLengthError2D_L);
    t_bkg->Branch("pseudoDecayLengthError3D_cand",&pseudoDecayLengthError3D_L);
    t_bkg->Branch("pdgId_cand",&pdgId_L);
    //t_bkg->Branch("trkIdx_cand",&trkIdx_L);
    t_bkg->Branch("trk_isHP_cand",&trk_isHP_L);
    t_bkg->Branch("trk_nHit_cand",&trk_nHit_L);
    t_bkg->Branch("trk_nChi2_cand",&trk_nChi2_L);
    t_bkg->Branch("trk_pTErr_cand",&trk_pTErr_L);
    t_bkg->Branch("trk_xyDCASignificance_cand",&trk_xyDCASignificance_L);
    t_bkg->Branch("trk_zDCASignificance_cand",&trk_zDCASignificance_L);
    
    
    TStopwatch stopwatchL;
    stopwatchL.Start();
    //Int_t counter = 0;
    Int_t count_cand = 0;
    Int_t count_gencand = 0;
    Bool_t skipEvent = kFALSE;
    cout << "LambdaId: " << LambdaId << endl;
    for(Int_t i = 0; i < nEntL; i++){ 
    //for(Int_t i = 0; i < 20; i++){ 
      
      //displayProgress(i,nEntL);  
      LambdaTree.GetEntry(i);
      hzVtxL->Fill(VtxZL);
      Long64_t tentryL = LambdaTree.LoadTree(i);
      LMass->GetEntry(tentryL);
      //cout << "Entry " << i << ":" << endl;
      /*cout << "Entry " << i << ":" << std::endl;
	Int_t count_cand = 0;
	
	for(size_t row = 0; row < cand_dauIdxL->size(); row++){
	for(size_t col = 0; col < (*cand_dauIdxL)[row].size(); col++){
	cout << "cand_dauIdxL[" << row << "][" << col << "] = " << (*cand_dauIdxL)[row][col] << " " << endl;
	}
	}
	
	for(size_t row = 0; row < cand_momIdxL->size(); row++){
	for(size_t col = 0; col < (*cand_momIdxL)[row].size(); col++){
	cout << "cand_momIdxL[" << row << "][" << col << "] = " << (*cand_momIdxL)[row][col] << " " << endl;
	}
	}*/
      
      /*for (auto rowIt = cand_dauIdxL->begin(); rowIt != cand_dauIdxL->end(); ++rowIt) {
	count_cand++;
	for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
	cout << "colIt: " << *colIt << endl;
	}
	//cout << endl;
	}*/
      //cout << "count_cand: " << count_cand << endl;
      //cout << "ending row of dauIdx: " << cand_dauIdxL->begin() << endl;

      /*for(Int_t j = 0; j < gen_massL->size(); j++){
	//cout << "candidate charge: " << cand_chargeL->at(j) << endl;
	cout << "\n" << gen_pdgIdL->at(j) << "\t"<< gen_massL->at(j) << endl;
	//if(gen_pdgIdL->at(j) != LambdaId)continue;
	}*/
      
      count_cand += cand_massL->size();
      count_gencand += gen_massL->size();
      //cout << "candidate number: " << cand_massL->size() << ""endl;
      for(Int_t j = 0; j < cand_massL->size(); j++){
	//cout << "candidate charge: " << cand_chargeL->at(j) << endl;
	hLMass->Fill(cand_massL->at(j));
	//cout << "\n" << cand_pdgIdL->at(j) << endl;
	if(cand_pdgIdL->at(j) != LambdaId)continue;
	if(!IsSatifyKineCuts(cand_pTL->at(j),cand_etaL->at(j)))continue;
	//if(cand_massL->at(j) < MinLambdaRecoMass)continue;
	hLMass_cand->Fill(cand_massL->at(j));
	hPtL->Fill(cand_pTL->at(j));
	hEtaL->Fill(cand_etaL->at(j));
	hPhiL->Fill(cand_phiL->at(j));
	h2LMassPdg_cand->Fill(cand_massL->at(j),cand_pdgIdL->at(j));
	
	h1decay->Fill(cand_decayLengthError2DL->at(j));
	h1pseudecay->Fill(cand_pseudoDecayLengthError2DL->at(j));
	
      }        

      if(cand_massL->size() > 0){
	Int_t cnt_reco = 0;
	for(Int_t icand_reco = 0; icand_reco < cand_massL->size(); icand_reco++){
	  if(cand_pdgIdL->at(icand_reco) != LambdaId)continue;

	  if(!IsSatifyKineCuts(cand_pTL->at(icand_reco),cand_etaL->at(icand_reco)))continue;
	  //if(abs(cand_etaL->at(icand_reco)) > 2.4)continue;
	  //if(cand_pTL->at(icand_reco)< 2)continue;
	  
	  cnt_reco++;
	  h1Pt_RA -> Fill(cand_pTL->at(icand_reco));
	  h1Mass_RA -> Fill(cand_massL->at(icand_reco));
	  h1Ncand_RA->Fill(cnt_reco);
	}
      }

      if(gen_massL->size() > 0){
	Int_t cnt_gen = 0;
	for(Int_t icand_gen = 0; icand_gen < gen_massL->size(); icand_gen++){
	  if(gen_pdgIdL->at(icand_gen) != LambdaId)continue;

	  if(!IsSatifyKineCuts(gen_pTL->at(icand_gen),gen_etaL->at(icand_gen)))continue;
	  //kinematic cuts
	  //if(abs(gen_etaL->at(icand_gen)) > 2.4)continue;
	  //if(gen_pTL->at(icand_gen)< 2)continue;
	  
	  cnt_gen++;
	  h1Pt_GA -> Fill(gen_pTL->at(icand_gen));
	  h1Mass_GA -> Fill(gen_massL->at(icand_gen));
	  h1Ncand_GA->Fill(cnt_gen);
	}
      }
      
      JetTree.GetEntry(i);
      
      hjetN->Fill(jetN);
      hjetPtLead->Fill(jetPtLead);
      
      Long64_t tentry = JetTree.LoadTree(i);
      chargedMult->GetEntry(tentry);
      
      Double_t RefJetPt = 0.;
      Int_t RefJetIndex = -1;
      
      if((jetEta->size() != jetPt->size()) || (jetEta->size() != jetPhi->size())){cout << jetEta->size() << "\t" << jetPt->size() << "\t" << jetPhi->size() << endl; cout << "Mismatch in jet numbers" << endl; return;}
      for(Int_t j = 0; j < jetPt->size(); j++){
	hjetEta->Fill(jetEta->at(j));
	hjetPt->Fill(jetPt->at(j));
	hjetPhi->Fill(jetPhi->at(j));
	hjetNumDaughters->Fill(jetNumDaughters->at(j));
	hchargedMult->Fill(chMult->at(j));
	//cout << chMult->at(j) << "\t " << jetNumDaughters->at(j) << "\t" << jetPt->at(j) << endl;
	
	if(jetPt->at(j) > RefJetPt){RefJetPt = jetPt->at(j); RefJetIndex = j;}
      }
      
      //if(cand_massL->size() == 0 && jetPt->size() > 0){
      //}
      
      //Double_t DelRA[100] = {0.};
      //for(Int_t It = 0; It < 100; It++)DelRA[It] = -99.0;
      
      

      if(cand_massL->size() > 0 && jetPt->size() > 0){ // there is at least one candidate and one jet
	for(Int_t j = 0; j < cand_massL->size(); j++){ // loop over reco cand
	  if(cand_pdgIdL->at(j) != LambdaId)continue;

	  if(!IsSatifyKineCuts(cand_pTL->at(j),cand_etaL->at(j)))continue;

	  hPtL_RA->Fill(cand_pTL->at(j));
	  //if(cand_massL->at(j) < MinLambdaRecoMass)continue;

	  //finding pair of reco and gen cand with shortest distance between them
	  refdist = 0.05; // decided from Dr_RG distribution
	  reffracpT = 0.2;
	  refIndR = -999;
	  refIndG = -999;
	  for(Int_t jgen = 0; jgen < gen_massL->size(); jgen++){
	    if(gen_pdgIdL->at(jgen) != LambdaId)continue;
	    if(!IsSatifyKineCuts(gen_pTL->at(jgen),gen_etaL->at(jgen)))continue;
	    
	    Dr_RG = dR(cand_etaL->at(j),cand_phiL->at(j),gen_etaL->at(jgen),gen_phiL->at(jgen));
	    FracPt_RG = ((cand_pTL->at(j))-(gen_pTL->at(jgen)))/(gen_pTL->at(jgen));

	    h1DelR_RG->Fill(Dr_RG);
	    h1FracPt_RG->Fill(FracPt_RG);
	    //h1FracPt_RG->Fill((abs((cand_pTL->at(j))-(gen_pTL->at(jgen))))/(gen_pTL->at(jgen)));
	    h1DelPt_RG->Fill((cand_pTL->at(j))-(gen_pTL->at(jgen)));
	   
	   
	    //if(Dr_RG > 0.05)continue;                        
	    //if(Dr_RG < refdist && FracPt_RG > reffracpT && ){
	    if(Dr_RG < refdist){
	      h1FracPt_RG_wcut->Fill(FracPt_RG);
	      //h1FracPt_RG_wcut->Fill((abs((cand_pTL->at(j))-(gen_pTL->at(jgen))))/(gen_pTL->at(jgen)));
	      h1DelPt_RG_wcut->Fill((cand_pTL->at(j))-(gen_pTL->at(jgen)));
	      
	      refdist = Dr_RG;
	      refIndR = j;
	      refIndG = jgen;
	    }
	  }

	  //if(refIndR != -999 && refIndG != -999){hPtL_GM->Fill(gen_pTL->at(refIndG));}
	  Int_t InsideJet = 0;
	  Bool_t skipJet = kFALSE;
	  for(Int_t j2 = 0; j2 < jetPt->size(); j2++){
	    if(skipJet)break;
	    if(jetPt->at(j2) < 550)continue; // selecting jet with pT > 550 GeV
	    if(abs(jetEta->at(j2)) > (2.4-JetR))continue; // selecting jet with the eta acceptance < |2.4-JetR|
	    
	    Double_t DelRV0Jet = dR(jetEta->at(j2),jetPhi->at(j2),cand_etaL->at(j),cand_phiL->at(j));
	    if(DelRV0Jet <= JetR){// candidate within the jet cone
	      //DelRA[j2] = DelRV0Jet;
	      
	      skipJet = kTRUE;//if cand is found within a jet, break the jet loop at next iter
	      if(refIndR != -999 && refIndG != -999){ // candidate is matched

		RefRMJetEta = jetEta->at(j2);
		RefRMJetPhi = jetPhi->at(j2);
		
		h2LMass_matched->Fill(cand_massL->at(refIndR),gen_massL->at(refIndG));
		
		//h3LMassMult_matched->Fill(chMult->at(j2),cand_massL->at(refIndR),gen_massL->at(refIndG));
		//h3LMasspT_matched->Fill(cand_pTL->at(refIndR),cand_massL->at(refIndR),gen_massL->at(refIndG));
		
                
		ChMultInJets_L = chMult->at(j2);
		
		mass_L = cand_massL->at(refIndR);
		charge_L = cand_chargeL->at(refIndR);
		//cout << "charge: " << charge_L << endl;
		status_L = cand_statusL->at(refIndR);
		y_L = cand_yL->at(refIndR);
		eta_L = cand_etaL->at(refIndR);
		phi_L = cand_phiL->at(refIndR);
		pT_L = cand_pTL->at(refIndR);
		p_L = cand_pL->at(refIndR);
		vtxProb_L = cand_vtxProbL->at(refIndR);
		vtxChi2_L = cand_vtxChi2L->at(refIndR);
		angle2D_L = cand_angle2DL->at(refIndR);
		angle3D_L = cand_angle3DL->at(refIndR);
		dca_L = cand_dcaL->at(refIndR);
		decayLength2D_L = cand_decayLength2DL->at(refIndR);
		decayLength3D_L = cand_decayLength3DL->at(refIndR);
                
		decayLengthError2D_L = cand_decayLengthError2DL->at(refIndR);
		decayLengthError3D_L = cand_decayLengthError3DL->at(refIndR);
		pseudoDecayLengthError2D_L = cand_pseudoDecayLengthError2DL->at(refIndR);
		pseudoDecayLengthError3D_L = cand_pseudoDecayLengthError3DL->at(refIndR);
		pdgId_L = cand_pdgIdL->at(refIndR);
		//trk_isHP_L = trk_isHPL->at(refIndR);
                
		//trk_nHit_L = trk_nHitL->at(refIndR);
		//trk_nChi2_L = trk_nChi2L->at(refIndR);
		//trk_pTErr_L = trk_pTErrL->at(refIndR);
		//trk_xyDCASignificance_L = trk_xyDCASignificanceL->at(refIndR);
		//trk_zDCASignificance_L = trk_zDCASignificanceL->at(refIndR);
                
		t_sig->Fill();

		/*if(genJetPt->size() == 0){
		  cout << "No gen jet found!" << endl;
		  skipEvent = kTRUE;
		  }*/
		//if((!skipEvent) && cand_vtxProbL->at(refIndR) > 0.05 && cand_angle3DL->at(refIndR) < 0.005){
		if(cand_vtxProbL->at(refIndR) < 0.05 || cand_angle3DL->at(refIndR) > 0.005)continue;
		hPtL_GM_in->Fill(gen_pTL->at(refIndG));

		hMassL_GM_in->Fill(cand_massL->at(refIndR));
		
		//cout << "refIndR: " << refIndR << endl;


	      }
	      else{                    
		hLMass_unmatched->Fill(cand_massL->at(j));
		
		ChMultInJets_L = chMult->at(j2);
		
		mass_L = cand_massL->at(j);
		charge_L = cand_chargeL->at(j);
		status_L = cand_statusL->at(j);
		y_L = cand_yL->at(j);
		eta_L = cand_etaL->at(j);
		phi_L = cand_phiL->at(j);
		pT_L = cand_pTL->at(j);
		p_L = cand_pL->at(j);
		vtxProb_L = cand_vtxProbL->at(j);
		vtxChi2_L = cand_vtxChi2L->at(j);
		angle2D_L = cand_angle2DL->at(j);
		angle3D_L = cand_angle3DL->at(j);
		dca_L = cand_dcaL->at(j);
		decayLength2D_L = cand_decayLength2DL->at(j);
		decayLength3D_L = cand_decayLength3DL->at(j);
                
		decayLengthError2D_L = cand_decayLengthError2DL->at(j);
		decayLengthError3D_L = cand_decayLengthError3DL->at(j);
		pseudoDecayLengthError2D_L = cand_pseudoDecayLengthError2DL->at(j);
		pseudoDecayLengthError3D_L = cand_pseudoDecayLengthError3DL->at(j);
		pdgId_L = cand_pdgIdL->at(j);
                
		//trk_isHP_L = trk_isHPL->at(j);
		//trk_nHit_L = trk_nHitL->at(j);
		//trk_nChi2_L = trk_nChi2L->at(j);
		//trk_pTErr_L = trk_pTErrL->at(j);
		//trk_xyDCASignificance_L = trk_xyDCASignificanceL->at(j);
		//trk_zDCASignificance_L = trk_zDCASignificanceL->at(j);
                
		
		t_bkg->Fill();
                
	      }
	      
	      
	      
	      
	      
	      
	      //hDelRV0Jet_in->Fill(DelRV0Jet);
	      //h3DelRV0Jet_in->Fill(jetPt->at(j2),cand_pTL->at(j),DelRV0Jet);
	      hLMass_in->Fill(cand_massL->at(j));
	      hPtL_in->Fill(cand_pTL->at(j));
	      hEtaL_in->Fill(cand_etaL->at(j));
	      hPhiL_in->Fill(cand_phiL->at(j));
	      h3LMass_in->Fill(chMult->at(j2),cand_pTL->at(j),cand_massL->at(j));
	      hchargedMult_wV0->Fill(chMult->at(j2));
	      hjetEta_wV0->Fill(jetEta->at(j2));
	      hjetPt_wV0->Fill(jetPt->at(j2));
	      hjetPhi_wV0->Fill(jetPhi->at(j2));
	      
	      h2chargedMultCandMinus->Fill(chMult->at(j2)-1,1.0/(chMult->at(j2)-1));
	      h2chargedMultCand->Fill(chMult->at(j2),1.0/chMult->at(j2));
	      
	      InsideJet++;
	      
	      if(j2 != RefJetIndex)continue;
	      //if(jetPt->at(j2) < 550)continue;
	      if(abs(jetEta->at(j2)) > (2.4-JetR))continue;
	      //hDelRV0Jet_in_PtLead->Fill(DelRV0Jet);
	      hLMass_in_PtLead->Fill(cand_massL->at(j));
	      hPtL_in_PtLead->Fill(cand_pTL->at(j));
	      hEtaL_in_PtLead->Fill(cand_etaL->at(j));
	      hPhiL_in_PtLead->Fill(cand_phiL->at(j));
	      h3LMass_in_PtLead->Fill(chMult->at(j2),cand_pTL->at(j),cand_massL->at(j));
	      hchargedMult_wV0_PtLead->Fill(chMult->at(j2));
	      hjetEta_wV0_PtLead->Fill(jetEta->at(j2));
	      hjetPt_wV0_PtLead->Fill(jetPt->at(j2));
	      hjetPhi_wV0_PtLead->Fill(jetPhi->at(j2));
              
	    }
	  }
	  if(!InsideJet){
	    //hDelRV0Jet_out->Fill(DelRV0Jet);
	    //h2DelRV0Jet_out->Fill(cand_pTL->at(j),DelRV0Jet);
	    hLMass_out->Fill(cand_massL->at(j));
	    hPtL_out->Fill(cand_pTL->at(j));
	    hEtaL_out->Fill(cand_etaL->at(j));
	    hPhiL_out->Fill(cand_phiL->at(j));
	  }
	  if(InsideJet > 1){
	    //counter++; 
	    //cout << "InsideJet: " << InsideJet << endl; 
	    //hInsideJet->Fill(InsideJet);
	    //for(Int_t It = 0; It < 100; It++){
	    //  if(DelRA[It] != -99.0)hDelRV0MJet->Fill(DelRA[It]);
	    //}
	  }
	
	  if(!skipJet)continue;
	  if(refIndR == -999 || refIndG == -999)continue;

	  if(cand_vtxProbL->at(refIndR) < 0.05 || cand_angle3DL->at(refIndR) > 0.005)continue;
	  
	  Int_t count_noGenJet = 0;
	  Int_t count_noGenJetwMV0 = 0;
	  Int_t count_GenJetwMV0 = 0;
	  if(genJetPt->size() > 0){
	    for(Int_t j2 = 0; j2 < genJetPt->size(); j2++){
	      if(genJetPt->at(j2) < 550)continue; // selecting jet with pT > 550 GeV
	      if(abs(genJetEta->at(j2)) > (2.4-JetR))continue; // selecting jet with the eta acceptance < |2.4-JetR|
	      Double_t DelRV0Jet_G = dR(genJetEta->at(j2),genJetPhi->at(j2),gen_etaL->at(refIndG),gen_phiL->at(refIndG));

	      hDelR_GM_in_wGJ->Fill(DelRV0Jet_G);
	      
	      if(DelRV0Jet_G <= JetR){
		count_GenJetwMV0++;
		hCount_GenJetwMV0->Fill(count_GenJetwMV0);
		
		hPtL_GM_in_wGJ->Fill(gen_pTL->at(refIndG));

		hMassL_GM_in_wGJ->Fill(cand_massL->at(refIndR));
		hDelR_GM_in_wGJ_wMV0->Fill(DelRV0Jet_G);

		
		  
		if(genChMult->at(j2) < 20)hPtL_GM_in_wGJ_NchLT20->Fill(gen_pTL->at(refIndG));
		if(genChMult->at(j2) > 80)hPtL_GM_in_wGJ_NchGT80->Fill(gen_pTL->at(refIndG));
		h2PtL_GM_in_wGJ->Fill(genChMult->at(j2),gen_pTL->at(refIndG));
		count_GenJetwMV0 = 0;

		if(RefRMJetEta == -999. || RefRMJetPhi == -999.)continue;
		Double_t DelR_Jet_RG_wMV0 = dR(genJetEta->at(j2),genJetPhi->at(j2),RefRMJetEta,RefRMJetPhi);
		hDelR_RGM_wMV0->Fill(DelR_Jet_RG_wMV0);
	      }
	      else{
		count_noGenJetwMV0++;
		hCount_noGenJetwMV0->Fill(count_noGenJetwMV0);
		count_noGenJetwMV0 = 0;
	      }
	      
	    }
	  }
	  else{
	    count_noGenJet++;
	    hCount_noGenJet->Fill(count_noGenJet);
	    count_noGenJet = 0;
	    //cout << "No gen jet found!" << endl;
	  }
	  
	  
	}
	//cout << "\nAn event with both jet and V0 particle found! No. of jets: " << jetPt->size() << " and no. of Lambdas: " << cand_massL->size() << endl;
	//break;



      }

      /*if(!skipEvent){
	
	}*/


      if(genJetPt->size() > 0 && jetPt->size() > 0){
	for(Int_t jRec = 0; jRec < jetPt->size(); jRec++){
	  for(Int_t j2 = 0; j2 < genJetPt->size(); j2++){
	    if(abs(jetEta->at(jRec)) > (2.4-JetR))continue;
	    if(abs(genJetEta->at(j2)) > (2.4-JetR))continue;
	    if(jetPt->at(jRec) < 550)continue;	    
	    Double_t DelR_Jet_RG = dR(genJetEta->at(j2),genJetPhi->at(j2),jetEta->at(jRec),jetPhi->at(jRec));
	    Double_t DelPhi_Jet_RG = dPhi(genJetPhi->at(j2),jetPhi->at(jRec));
	    Double_t DelEta_Jet_RG = genJetEta->at(j2) - jetEta->at(jRec);
	    hDelR_RGJet->Fill(DelR_Jet_RG);
	    hDelPhi_RGJet->Fill(DelPhi_Jet_RG);
	    hDelEta_RGJet->Fill(DelEta_Jet_RG);

	    
	    h1GJetEta->Fill(genJetEta->at(j2));
	    h1GJetPhi->Fill(genJetPhi->at(j2));

	    h1RJetEta->Fill(jetEta->at(jRec));
	    h1RJetPhi->Fill(jetPhi->at(jRec));
	    
	    //if(abs(jetEta->at(jRec)) > (2.4-JetR))continue;
	    h2GJetPtDelR->Fill(DelR_Jet_RG,genJetPt->at(j2));
	  }
	}
      }

      if(gen_massL->size() > 0){ // there is at least one gen candidate
	for(Int_t jG = 0; jG < gen_massL->size(); jG++){ // loop over gen cand
	  h1pdgId_G->Fill(gen_pdgIdL->at(jG));
	  if(gen_pdgIdL->at(jG) != LambdaId)continue;
	  if(!IsSatifyKineCuts(gen_pTL->at(jG),gen_etaL->at(jG)))continue;
	  //cout << gen_pdgIdL->at(jG) << endl;
	  hPtL_GA_wojet->Fill(gen_pTL->at(jG));
	}
      }
      
      if(gen_massL->size() > 0 && genJetPt->size() > 0){ // there is at least one candidate and one jet
	for(Int_t jG = 0; jG < gen_massL->size(); jG++){ // loop over gen cand
	  h1pdgId_G->Fill(gen_pdgIdL->at(jG));
	  if(gen_pdgIdL->at(jG) != LambdaId)continue;
	  if(!IsSatifyKineCuts(gen_pTL->at(jG),gen_etaL->at(jG)))continue;
	  //cout << gen_pdgIdL->at(jG) << endl;
	  hPtL_GA->Fill(gen_pTL->at(jG));
	  Bool_t skipGenJet = kFALSE;
	  for(Int_t j2 = 0; j2 < genJetPt->size(); j2++){
	    if(skipGenJet)break;
	    if(genJetPt->at(j2) < 550)continue; // selecting jet with pT > 550 GeV
	    if(abs(genJetEta->at(j2)) > (2.4-JetR))continue; // selecting jet with the eta acceptance < |2.4-JetR|
	    
	    Double_t DelRV0Jet = dR(genJetEta->at(j2),genJetPhi->at(j2),gen_etaL->at(jG),gen_phiL->at(jG));
	    if(DelRV0Jet <= JetR){// candidate within the jet cone
	      skipGenJet = kTRUE;
	      hPtL_GA_in->Fill(gen_pTL->at(jG));
	      if(cand_massL->size() > 0)hPtL_GA_in_wcand->Fill(gen_pTL->at(jG));
	      hMassL_GA_in->Fill(gen_massL->at(jG));
	      if(genChMult->at(j2) < 20)hPtL_GA_in_NchLT20->Fill(gen_pTL->at(jG));
	      if(genChMult->at(j2) > 80)hPtL_GA_in_NchGT80->Fill(gen_pTL->at(jG));
	      h2PtL_GA_in->Fill(genChMult->at(j2),gen_pTL->at(jG));
	    }
	  }
	}
      }
      
      //if(RefRMJetEta == -999. || RefRMJetPhi == -999.)continue;
      //Double_t DelR_Jet_RG_wMV0 = dR(genJetEta->at(j2),genJetPhi->at(j2),RefRMJetEta,RefRMJetPhi);
      //hDelR_RGM_wMV0->Fill(DelR_Jet_RG_wMV0);
      
    }
    cout << "candidates: " << count_cand << "\t gen candidates: " << count_gencand << endl;
    t_sig->Write();
    t_bkg->Write();
    fOutTree->Close();
    
    //cout << hchargedMult->GetEntries() << "\t" << hchargedMult->GetMean() << endl;
    //cout << "\n counter: " << counter << endl;
    
    stopwatchL.Stop();
    cout << "\nReal time: " << stopwatchL.RealTime() << "\t CPU time: " << stopwatchL.CpuTime() << endl;
    
    JetTree.ResetBranchAddresses();
    LambdaTree.ResetBranchAddresses();
    
    
    TCanvas *cL = new TCanvas();
    //hLMass->Draw("e"); 
    hLMass->SetTitle(""); hLMass->GetXaxis()->SetTitle("p #pi^{-} + #bar{p} #pi^{+} invariant mass [GeV]"); hLMass->GetYaxis()->SetTitle("No. of #Lambda candidates"); //hLMass->GetXaxis()->SetRangeUser(1.08,1.16);// hLMass->SetStats(0);
    
    
  }
  
  if(AnaKshort && !AnaLambda){
    Int_t nEntK = KshortTree.GetEntries();
    cout << "\n Number of Entries K : " << nEntK << endl;
    
    Float_t VtxZK = 0;
    KshortTree.SetBranchAddress("bestvtxZ",&VtxZK);
    
    vector<float> *cand_massK = 0;
    
    TBranch *KMass = 0;
    KshortTree.SetBranchAddress("cand_mass",&cand_massK,&KMass);
    
    for(Int_t i = 0; i < nEntK; i++){  
      displayProgress(i,nEntK); 
      KshortTree.GetEntry(i);
      hzVtxK->Fill(VtxZK);
      Long64_t tentryK = KshortTree.LoadTree(i);
      KMass->GetEntry(tentryK);
      
      for(Int_t j = 0; j < cand_massK->size(); j++){
	hKMass->Fill(cand_massK->at(j));
      }
    }
    KshortTree.ResetBranchAddresses();
    
    TCanvas *cK = new TCanvas();
    //hKMass->Draw("e"); 
    hKMass->SetTitle(""); hKMass->GetXaxis()->SetTitle("#pi^{+} #pi^{-} invariant mass [GeV]"); hKMass->GetYaxis()->SetTitle("No. of K^{0}_{S} candidates"); //hKMass->GetXaxis()->SetRangeUser(0.4,0.6);// hKMass->SetStats(0);
  }
  
  if(AnaLambda && AnaKshort){
    Int_t nEntL = LambdaTree.GetEntries();
    cout << "\n Number of Entries L : " << nEntL << endl;
    Int_t nEntK = KshortTree.GetEntries();
    cout << "\n Number of Entries K : " << nEntK << endl;
    if(nEntL != nEntK){cout << "Mismatch in no. of entries! Fatal Error" << endl; return;} 
    
    Float_t VtxZL = 0;
    LambdaTree.SetBranchAddress("bestvtxZ",&VtxZL);
    Float_t VtxZK = 0; 
    KshortTree.SetBranchAddress("bestvtxZ",&VtxZK);
    
    vector<float> *cand_massL = 0;
    vector<float> *cand_massK = 0;
    
    TBranch *LMass = 0;
    LambdaTree.SetBranchAddress("cand_mass",&cand_massL,&LMass);
    TBranch *KMass = 0;
    KshortTree.SetBranchAddress("cand_mass",&cand_massK,&KMass);
    
    TStopwatch stopwatch;
    stopwatch.Start();
    for(Int_t i = 0; i < nEntL; i++){ 
      displayProgress(i,nEntL);  
      
      LambdaTree.GetEntry(i);
      hzVtxL->Fill(VtxZL);
      Long64_t tentryL = LambdaTree.LoadTree(i);
      LMass->GetEntry(tentryL);
      
      KshortTree.GetEntry(i);
      hzVtxK->Fill(VtxZK);
      Long64_t tentryK = KshortTree.LoadTree(i);
      KMass->GetEntry(tentryK);
      
      /*for(Int_t j = 0; j < cand_massL->size(); j++){
	hLMass->Fill(cand_massL->at(j));
	}
	for(Int_t j = 0; j < cand_massK->size(); j++){
	  hKMass->Fill(cand_massK->at(j));
	  }*/
      
      
      Int_t NoL = cand_massL->size();
      Int_t NoK = cand_massK->size();
      Int_t NoMax = 0;
      Bool_t IsLgtK = kFALSE;
      if(NoL>NoK){
	NoMax = NoL;
	IsLgtK = kTRUE;
      }
      else{
	NoMax = NoK;
      }
      //cout << "NoL: " << NoL << "\t NoK: " << NoK << "\t NoMax: " << NoMax << endl;
      
      for(Int_t j = 0; j < NoMax; j++){
	if(IsLgtK){
	  //cout << "Hello1" << endl;
	  hLMass->Fill(cand_massL->at(j));
	  if(j >= NoK)continue;
	  //cout << "Hello2" << endl;
	  hKMass->Fill(cand_massK->at(j));
	}
	else{
	  //cout << "Hello3" << endl;
	  hKMass->Fill(cand_massK->at(j));
	  if(j >= NoL)continue;
	  //cout << "Hello4" << endl;
	  hLMass->Fill(cand_massL->at(j));
	}
      }
      
    }
    cout << endl;
    stopwatch.Stop();
    cout << "Real time: " << stopwatch.RealTime() << "\t CPU time: " << stopwatch.CpuTime() << endl;
    LambdaTree.ResetBranchAddresses();
    KshortTree.ResetBranchAddresses();
    
    TCanvas *cL = new TCanvas();
    //hLMass->Draw("e"); 
    hLMass->SetTitle(""); hLMass->GetXaxis()->SetTitle("p #pi^{-} + #bar{p} #pi^{+} invariant mass [GeV]"); hLMass->GetYaxis()->SetTitle("No. of #Lambda candidates"); //hLMass->GetXaxis()->SetRangeUser(1.08,1.16);// hLMass->SetStats(0);
    TCanvas *cK = new TCanvas();
    //hKMass->Draw("e"); 
    hKMass->SetTitle(""); hKMass->GetXaxis()->SetTitle("#pi^{+} #pi^{-} invariant mass [GeV]"); hKMass->GetYaxis()->SetTitle("No. of K^{0}_{S} candidates"); //hKMass->GetXaxis()->SetRangeUser(0.4,0.6);// hKMass->SetStats(0);
  }    
  
  
  if(IsSave){
    TFile *fO = new TFile(Form("%s.root",OutputFileName.Data()),"RECREATE");
    
    h1decay->Write();
    h1pseudecay->Write();
      
    hzVtxL->Write("bestVtxZL");
    hLMass->Write("InvMassL");
    hzVtxK->Write("bestVtxZK");
    hKMass->Write("InvMassK");
    
    h2LMassPdg_cand->Write("h2LMassPdg_cand");
    h2LMass_matched->Write();
    hLMass_unmatched->Write();
    hLMass_cand->Write("InvMassL_cand");
    hLMass_in->Write("InvMassL_in");
    hLMass_in_PtLead->Write("InvMassL_in_PtLead");
    hLMass_out->Write("InvMassL_out");
    
    hjetN->Write();
    hjetEta->Write();
    hjetPhi->Write();
    hjetPt->Write();
    hjetPtLead->Write();
    
    hjetEta_wV0->Write();
    hjetPhi_wV0->Write();
    hjetPt_wV0->Write();
    
    hjetEta_wV0_PtLead->Write();
    hjetPhi_wV0_PtLead->Write();
    hjetPt_wV0_PtLead->Write();
    
    
    hjetNumDaughters->Write();
    hchargedMult->Write();
    hchargedMult_wV0->Write();
    hchargedMult_wV0_PtLead->Write();
    
    h2chargedMultCandMinus->Write();
    h2chargedMultCand->Write();
    
    //hDelRV0Jet_in->Write();
    //hDelRV0Jet_in_PtLead->Write();
    //hDelRV0Jet_out->Write();
    //h3DelRV0Jet_in->Write();
    //h2DelRV0Jet_out->Write();
    
    hPtL->Write();
    hPtL_in->Write();
    hPtL_in_PtLead->Write();
    hPtL_out->Write();
    
    hEtaL->Write();
    hEtaL_in->Write();
    hEtaL_in_PtLead->Write();
    hEtaL_out->Write();
    
    hPhiL->Write();
    hPhiL_in->Write();
    hPhiL_in_PtLead->Write();
    hPhiL_out->Write();  
    
    //hInsideJet->Write();    
    //hDelRV0MJet->Write();  
    
    h3LMass_in->Write();        
    h3LMass_in_PtLead->Write();        

    h1DelR_RG->Write();
    h1FracPt_RG->Write();
    h1DelPt_RG->Write();
    h1pdgId_G->Write();

    h1FracPt_RG_wcut->Write();
    h1DelPt_RG_wcut->Write();

    hPtL_RA->Write();
    hPtL_GA_wojet->Write();
    
    hPtL_GA->Write();
    hPtL_GA_in->Write();
    hPtL_GA_in_NchLT20->Write();
    hPtL_GA_in_NchGT80->Write();
    h2PtL_GA_in->Write();

    hPtL_GA_in_wcand->Write();
    hMassL_GA_in->Write();
    
    hPtL_GM_in->Write();
    hPtL_GM_in_wGJ->Write();
    hPtL_GM_in_wGJ_NchLT20->Write();
    hPtL_GM_in_wGJ_NchGT80->Write();
    h2PtL_GM_in_wGJ->Write();

    hMassL_GM_in->Write();
    hMassL_GM_in_wGJ->Write();
    hDelR_GM_in_wGJ->Write();
    hDelR_GM_in_wGJ_wMV0->Write();

    hDelR_RGM_wMV0->Write();
    
    h2GJetPtDelR->Write();
    
    hCount_noGenJetwMV0->Write();
    hCount_GenJetwMV0->Write();
    hCount_noGenJet->Write();

    hDelR_RGJet->Write();
    hDelPhi_RGJet->Write();
    hDelEta_RGJet->Write();
    
    h1RJetEta->Write();
    h1RJetPhi->Write();
    h1GJetEta->Write();
    h1GJetPhi->Write();

    h1Pt_RA->Write();
    h1Pt_GA->Write();
    h1Mass_RA->Write();
    h1Mass_GA->Write();

    h1Ncand_RA->Write();
    h1Ncand_GA->Write();
    
    fO->Close();
  }
  cout << endl;
  
}
