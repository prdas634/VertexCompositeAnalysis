#include "MyConfigAll_v2.h"

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

Bool_t IsSatisfyKineCuts(Double_t RefPt, Double_t RefEta){
  if(RefPt < 1.0 || TMath::Abs(RefEta) > 2.4)return kFALSE;
  return kTRUE;
}


void GetObs_v2(Int_t nIpDir = 0, Int_t iCand = 0, Bool_t IsSave = kFALSE){
  
  files.clear();
  
  GetFiles((inputDirDataParent+inputDirData[nIpDir]).Data(),files,nfiles);
  cout<<"Got "<<files.size()<< " files" << endl;
  //for(Int_t i = 0; i <= nfiles+1; i++)cout << files[i] << endl;
  
  TChain JetTree("JetTree"), CandTree(Form("%sTree",cand_tree[iCand].Data()));
  for(TString temp_path:files){
    JetTree.Add(temp_path+"/analyzerOffline/trackTree");
    CandTree.Add(temp_path+"/"+cand_tree_ana[iCand]+"/ParticleTree");
  }
  
  Int_t nEntL = CandTree.GetEntries();
  cout << "\n Number of Entries L : " << nEntL << endl;
  
  Int_t nEntJ = JetTree.GetEntries();
  cout << "\n Number of Entries Jet : " << nEntJ << endl;

  Float_t VtxZL = 0;
  CandTree.SetBranchAddress("bestvtxZ",&VtxZL);

  Int_t jetN;
  Float_t jetPtLead;
  JetTree.SetBranchAddress("jetN",&jetN);
  JetTree.SetBranchAddress("jetPtLead",&jetPtLead);
  
  vector<float> *cand_massL = 0;
  vector<float> *cand_etaL = 0;
  vector<float> *cand_phiL = 0;
  vector<float> *cand_pTL = 0;
  vector<int> *cand_pdgIdL = 0;

  vector<float> *jetEta = 0;
  vector<float> *jetPt = 0;
  vector<float> *jetPhi = 0;
  vector<int> *jetNumDaughters = 0;
  vector<int> *chMult = 0;
  
  TBranch *LMass = 0;
  CandTree.SetBranchAddress("cand_mass",&cand_massL,&LMass);
  TBranch *LEta = 0;
  CandTree.SetBranchAddress("cand_eta",&cand_etaL,&LEta);
  TBranch *LPhi = 0;
  CandTree.SetBranchAddress("cand_phi",&cand_phiL,&LPhi);
  TBranch *LPt = 0;
  CandTree.SetBranchAddress("cand_pT",&cand_pTL,&LPt);
  TBranch *LPdgId = 0;
  CandTree.SetBranchAddress("cand_pdgId",&cand_pdgIdL,&LPdgId);
  
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
  

  TStopwatch stopwatch;
  stopwatch.Start();
  Int_t counter = 0;
  for(Int_t i = 0; i < nEntL; i++){ 
  //for(Int_t i = 0; i < 1000; i++){ 

    //displayProgress(i,nEntL);  
    CandTree.GetEntry(i);
    hzVtxL->Fill(VtxZL);
    Long64_t tentryL = CandTree.LoadTree(i);
    LMass->GetEntry(tentryL);

    JetTree.GetEntry(i);
    
    Long64_t tentry = JetTree.LoadTree(i);
    chargedMult->GetEntry(tentry);
  
    if((jetEta->size() != jetPt->size()) || (jetEta->size() != jetPhi->size())){cout << jetEta->size() << "\t" << jetPt->size() << "\t" << jetPhi->size() << endl; cout << "Mismatch in jet numbers" << endl; return;}
    for(Int_t j = 0; j < jetPt->size(); j++){
      if(jetPt->at(j) < 550)continue;
      if(abs(jetEta->at(j)) > (2.4-JetR))continue;      
      hjetEta->Fill(jetEta->at(j));
      hjetPt->Fill(jetPt->at(j));
      hjetPhi->Fill(jetPhi->at(j));
      hjetNumDaughters->Fill(jetNumDaughters->at(j));
      hchargedMult->Fill(chMult->at(j));
      //cout << chMult->at(j) << "\t " << jetNumDaughters->at(j) << "\t" << jetPt->at(j) << endl;
    }
  
    if(cand_massL->size() > 0 && jetPt->size() > 0){
      //cout << "Here I am" << endl;
      //counter = 0;
      for(Int_t j = 0; j < cand_massL->size(); j++){
	if(cand_pdgIdL->at(j) != cand_ID[iCand])continue;
	//if(cand_massL->at(j) < MinLambdaRecoMass)continue;
	if(!IsSatisfyKineCuts(cand_pTL->at(j),cand_etaL->at(j)))continue;
	//cout << "got candidate\t";
	Bool_t skipJet = kFALSE;
	for(Int_t j2 = 0; j2 < jetPt->size(); j2++){
	  if(skipJet)break;//the cand is already found within a jet, so break the jet loop
	  if(jetPt->at(j2) < 550)continue;
	  if(abs(jetEta->at(j2)) > (2.4-JetR))continue;
	  Double_t DelRV0Jet = dR(jetEta->at(j2),jetPhi->at(j2),cand_etaL->at(j),cand_phiL->at(j));
	  if(DelRV0Jet <= JetR){
	    skipJet = kTRUE;//the cand is found within the jet
	    counter++;
	    hLMass_in->Fill(cand_massL->at(j));
	    hpdgIdL_in->Fill(cand_pdgIdL->at(j));
	    hPtL_in->Fill(cand_pTL->at(j));
	    hEtaL_in->Fill(cand_etaL->at(j));
	    hPhiL_in->Fill(cand_phiL->at(j));
	    h3LMass_in->Fill(chMult->at(j2),cand_pTL->at(j),cand_massL->at(j));
	    hchargedMult_wV0->Fill(chMult->at(j2));
	    hjetEta_wV0->Fill(jetEta->at(j2));
	    hjetPt_wV0->Fill(jetPt->at(j2));
	    hjetPhi_wV0->Fill(jetPhi->at(j2));

	    h2chargedMultCandMinus2->Fill(chMult->at(j2)-2,1.0/(chMult->at(j2)-2));
	    h2chargedMultCandMinus->Fill(chMult->at(j2)-1,1.0/(chMult->at(j2)-1));
	    h2chargedMultCand->Fill(chMult->at(j2),1.0/chMult->at(j2));
	  }
	}
      }
    }

  }

  cout << "No of candidates: " << counter << endl;
  
  stopwatch.Stop();

  cout << "Real time: " << stopwatch.RealTime() << "\t CPU time: " << stopwatch.CpuTime() << endl;
  
  JetTree.ResetBranchAddresses();
  CandTree.ResetBranchAddresses();
  
  
  
  if(IsSave){
    TFile *fO = new TFile(Form("%s_%s.root",OutputFileName.Data(),cand_tree[iCand].Data()),"RECREATE");
    hzVtxL->Write("bestVtxZL");    
    hLMass_in->Write("InvMassL_in");
    hpdgIdL_in->Write();
    hjetEta->Write();
    hjetPhi->Write();
    hjetPt->Write();
    hjetNumDaughters->Write();
    hchargedMult->Write();
    
    hjetEta_wV0->Write();
    hjetPhi_wV0->Write();
    hjetPt_wV0->Write();
    hchargedMult_wV0->Write();

    h2chargedMultCandMinus2->Write();
    h2chargedMultCandMinus->Write();
    h2chargedMultCand->Write();
    
    hPtL_in->Write();
    hEtaL_in->Write();
    hPhiL_in->Write();
        
    h3LMass_in->Write();        
        
    fO->Close();
  }
  
  cout << endl;
    
}
