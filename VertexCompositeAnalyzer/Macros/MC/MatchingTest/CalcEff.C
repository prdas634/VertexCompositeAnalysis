#include "AllConfig.h"

void CalcEff(Int_t nIpDir = 0, Int_t iCand = 0, Bool_t IsSave = kFALSE){
    files.clear();
  
    GetFiles((inputDirDataParent+inputDirData[nIpDir]).Data(),files,nfiles);
    cout<<"Got "<<files.size()<< " files" << endl;

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

    vector<float> *cand_vtxProbL = 0;
    vector<float> *cand_angle3DL = 0;

    vector<float> *gen_massL = 0;
    vector<float> *gen_etaL = 0;
    vector<float> *gen_phiL = 0;
    vector<float> *gen_pTL = 0;
    vector<int> *gen_pdgIdL = 0;

    vector<float> *jetEta = 0;
    vector<float> *jetPt = 0;
    vector<float> *jetPhi = 0;
    vector<int> *jetNumDaughters = 0;
    vector<int> *chMult = 0;

    vector<float> *genJetEta = 0;
    vector<float> *genJetPt = 0;
    vector<float> *genJetPhi = 0;
    vector<int> *genChMult = 0;    

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

    TBranch *LvtxProb = 0;
    CandTree.SetBranchAddress("cand_vtxProb",&cand_vtxProbL,&LvtxProb);
    TBranch *Langle3D = 0;
    CandTree.SetBranchAddress("cand_angle3D",&cand_angle3DL,&Langle3D);

    TBranch *LMass_G = 0;
    CandTree.SetBranchAddress("gen_mass",&gen_massL,&LMass_G);
    TBranch *LEta_G = 0;
    CandTree.SetBranchAddress("gen_eta",&gen_etaL,&LEta_G);
    TBranch *LPhi_G = 0;
    CandTree.SetBranchAddress("gen_phi",&gen_phiL,&LPhi_G);
    TBranch *LPt_G = 0;
    CandTree.SetBranchAddress("gen_pT",&gen_pTL,&LPt_G);
    TBranch *LPdgId_G = 0;
    CandTree.SetBranchAddress("gen_pdgId",&gen_pdgIdL,&LPdgId_G);
  
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

    TBranch *jEta_G = 0;
    JetTree.SetBranchAddress("genJetEta",&genJetEta,&jEta_G); 
    TBranch *jPt_G = 0;
    JetTree.SetBranchAddress("genJetPt",&genJetPt,&jPt_G);
    TBranch *jPhi_G = 0;
    JetTree.SetBranchAddress("genJetPhi",&genJetPhi,&jPhi_G);
    TBranch *chargedMult_G = 0;
    JetTree.SetBranchAddress("genJetChargedMultiplicity",&genChMult,&chargedMult_G);

    TStopwatch stopwatch;
    stopwatch.Start();
    Int_t count_cand = 0;
    Int_t count_gencand = 0;
    Bool_t skipEvent = kFALSE;
    Int_t counter = 0;
    for(Int_t i = 0; i < nEntL; i++){ 
    //for(Int_t i = 0; i < 500; i++){
    
        //displayProgress(i,nEntL);  
        CandTree.GetEntry(i);
        hzVtxL->Fill(VtxZL);
        Long64_t tentryL = CandTree.LoadTree(i);
        LMass->GetEntry(tentryL);

        JetTree.GetEntry(i);
    
        Long64_t tentry = JetTree.LoadTree(i);
        chargedMult->GetEntry(tentry);    
        
        if((jetEta->size() != jetPt->size()) || (jetEta->size() != jetPhi->size())){cout << jetEta->size() << "\t" << jetPt->size() << "\t" << jetPhi->size() << endl; cout << "Mismatch in jet numbers" << endl; return;}
        
        hjetN->Fill(jetN);
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

        count_cand += cand_massL->size();//reco cand
        count_gencand += gen_massL->size();//gen cand
        
        if(cand_massL->size() > 0){
	        for(Int_t icand_reco = 0; icand_reco < cand_massL->size(); icand_reco++){
            hLMass->Fill(cand_massL->at(icand_reco));
            h2LMassPdg_cand->Fill(cand_massL->at(icand_reco),cand_pdgIdL->at(icand_reco));
	          if(cand_pdgIdL->at(icand_reco) != cand_ID[iCand])continue;
	          if(!IsSatisfyKineCuts(cand_pTL->at(icand_reco),cand_etaL->at(icand_reco)))continue;
            if(cand_vtxProbL->at(icand_reco) < 0.05 || cand_angle3DL->at(icand_reco) > 0.005)continue;
            hLMass_cand->Fill(cand_massL->at(icand_reco));
	          hEtaL->Fill(cand_etaL->at(icand_reco));
	          hPhiL->Fill(cand_phiL->at(icand_reco));
	          h1Pt_RA -> Fill(cand_pTL->at(icand_reco));//pT dist of reco cand irrespective of existing jet
	          h1Mass_RA -> Fill(cand_massL->at(icand_reco));
	          h1Ncand_RA->Fill(1.0);
          }
        }


        
        if(cand_massL->size() > 0 && gen_massL->size() > 0){
          for(Int_t j = 0; j < cand_massL->size(); j++){//loop over reco cand while there is at least one jet reconstructed in the event
            if(cand_pdgIdL->at(j) != cand_ID[iCand])continue;
            if(!IsSatisfyKineCuts(cand_pTL->at(j),cand_etaL->at(j)))continue;
            if(cand_vtxProbL->at(j) < 0.05 || cand_angle3DL->at(j) > 0.005)continue;

            //finding pair of reco and gen cand with shortest distance between them
            refdistAll = 100.; 
            refIndRAll = -9999;
	          refIndGAll = -9999;
            refIndR = -999;
	          refIndG = -999;
            Dr_RG = 1000.;

            for(Int_t jgen = 0; jgen < gen_massL->size(); jgen++){
              if(gen_pdgIdL->at(jgen) != cand_ID[iCand])continue;
	            if(!IsSatisfyKineCuts(gen_pTL->at(jgen),gen_etaL->at(jgen)))continue;

              Dr_RG = dR(cand_etaL->at(j),cand_phiL->at(j),gen_etaL->at(jgen),gen_phiL->at(jgen));
              FracPt_RG = ((cand_pTL->at(j))-(gen_pTL->at(jgen)))/(gen_pTL->at(jgen));

              if(Dr_RG < refdistAll){	                
                refdistAll = Dr_RG; 
                refIndRAll = j;
	              refIndGAll = jgen;
              }
            }
            
            if(refIndRAll != -9999 && refIndGAll != -9999){ // closest gen and reco candidates 
              h1DelR_RG_closest_woJ->Fill(dR(cand_etaL->at(refIndRAll),cand_phiL->at(refIndRAll),gen_etaL->at(refIndGAll),gen_phiL->at(refIndGAll)));  
              h1FracPt_RG_closest_woJ->Fill(((cand_pTL->at(refIndRAll))-(gen_pTL->at(refIndGAll)))/(gen_pTL->at(refIndGAll)));
            }

            if(refdistAll <= refdist){
	            refIndR = refIndRAll;
	            refIndG = refIndGAll;
            }

            if(refIndR != -999 && refIndG != -999){ // candidate is matched
              h1DelR_RG_matched_woJ->Fill(dR(cand_etaL->at(refIndR),cand_phiL->at(refIndR),gen_etaL->at(refIndG),gen_phiL->at(refIndG)));  
              h1FracPt_RG_matched_woJ->Fill(((cand_pTL->at(refIndR))-(gen_pTL->at(refIndG)))/(gen_pTL->at(refIndG)));
              h1Pt_GM_woJ -> Fill(gen_pTL->at(refIndG));//pT dist of matched gen cand irrespective of existing jet
            }
            
          }
        }

        if(gen_massL->size() > 0){
	        for(Int_t icand_gen = 0; icand_gen < gen_massL->size(); icand_gen++){
	          if(gen_pdgIdL->at(icand_gen) != cand_ID[iCand])continue;
            if(!IsSatisfyKineCuts(gen_pTL->at(icand_gen),gen_etaL->at(icand_gen)))continue;
	          h1Pt_GA -> Fill(gen_pTL->at(icand_gen));//pT dist of gen cand irrespective of existing jet
	          h1Mass_GA -> Fill(gen_massL->at(icand_gen));
	          h1Ncand_GA->Fill(1.0);
	        }
        }

        if(cand_massL->size() > 0 && jetPt->size() > 0){
          for(Int_t j = 0; j < cand_massL->size(); j++){//loop over reco cand while there is at least one jet reconstructed in the event
            if(cand_pdgIdL->at(j) != cand_ID[iCand])continue;
            if(!IsSatisfyKineCuts(cand_pTL->at(j),cand_etaL->at(j)))continue;
            if(cand_vtxProbL->at(j) < 0.05 || cand_angle3DL->at(j) > 0.005)continue;
            hPtL_RA->Fill(cand_pTL->at(j));//pT dist of reco cand with cuts while there is at least one jet reconstructed in the event             

            //finding pair of reco and gen cand with shortest distance between them
            refdistAll = 100.; 
            refIndRAll = -9999;
	          refIndGAll = -9999;
            refIndR = -999;
	          refIndG = -999;
            Dr_RG = 1000.;

            for(Int_t jgen = 0; jgen < gen_massL->size(); jgen++){
              if(gen_pdgIdL->at(jgen) != cand_ID[iCand])continue;
	            if(!IsSatisfyKineCuts(gen_pTL->at(jgen),gen_etaL->at(jgen)))continue;

              Dr_RG = dR(cand_etaL->at(j),cand_phiL->at(j),gen_etaL->at(jgen),gen_phiL->at(jgen));
              FracPt_RG = ((cand_pTL->at(j))-(gen_pTL->at(jgen)))/(gen_pTL->at(jgen));

	            h1DelR_RG->Fill(Dr_RG);
	            h1FracPt_RG->Fill(FracPt_RG);
	            h1DelPt_RG->Fill((cand_pTL->at(j))-(gen_pTL->at(jgen)));

              if(Dr_RG < refdistAll){
                h1FracPt_RG_wcut->Fill(FracPt_RG);
                h1DelPt_RG_wcut->Fill((cand_pTL->at(j))-(gen_pTL->at(jgen)));    
	                
                refdistAll = Dr_RG; 
                refIndRAll = j;
	              refIndGAll = jgen;
              }
            }
            
            if(refIndRAll != -9999 && refIndGAll != -9999){ // closest gen and reco candidates 
              h1DelR_RG_closest_wJ->Fill(dR(cand_etaL->at(refIndRAll),cand_phiL->at(refIndRAll),gen_etaL->at(refIndGAll),gen_phiL->at(refIndGAll)));  
              h1FracPt_RG_closest_wJ->Fill(((cand_pTL->at(refIndRAll))-(gen_pTL->at(refIndGAll)))/(gen_pTL->at(refIndGAll)));
            }

            if(refdistAll <= refdist){
	            refIndR = refIndRAll;
	            refIndG = refIndGAll;
            }

            if(refIndR != -999 && refIndG != -999){ // candidate is matched
              h1DelR_RG_matched_wJ->Fill(dR(cand_etaL->at(refIndR),cand_phiL->at(refIndR),gen_etaL->at(refIndG),gen_phiL->at(refIndG)));  
              h1FracPt_RG_matched_wJ->Fill(((cand_pTL->at(refIndR))-(gen_pTL->at(refIndG)))/(gen_pTL->at(refIndG)));    
              
              Int_t InsideJet = 0;   
              Bool_t skipJet = kFALSE;
	            for(Int_t j2 = 0; j2 < jetPt->size(); j2++){//loop over reco jets starting with leading pT jet
                if(skipJet)break; //if a jet with matched cand found, then break
	              if(jetPt->at(j2) < 550)continue; // selecting jet with pT > 550 GeV
	              if(abs(jetEta->at(j2)) > (2.4-JetR))continue; // selecting jet with the eta acceptance < |2.4-JetR|

                Double_t DelRV0Jet = dR(jetEta->at(j2),jetPhi->at(j2),cand_etaL->at(j),cand_phiL->at(j));
                if(DelRV0Jet <= JetR){// candidate within the jet cone
                  skipJet = kTRUE;//if matched cand is found within a jet, break the jet loop at next iter
                  counter++;
                  RefRMJetEta = jetEta->at(j2);
                  RefRMJetPhi = jetPhi->at(j2);
	                h2LMass_matched->Fill(cand_massL->at(refIndR),gen_massL->at(refIndG));
                  hPtL_GM_in->Fill(gen_pTL->at(refIndG));
		              hMassL_GM_in->Fill(cand_massL->at(refIndR));
                  	                
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

                  h1DelR_RG_matched_inJ->Fill(dR(cand_etaL->at(refIndR),cand_phiL->at(refIndR),gen_etaL->at(refIndG),gen_phiL->at(refIndG)));  
                  h1FracPt_RG_matched_inJ->Fill(((cand_pTL->at(refIndR))-(gen_pTL->at(refIndG)))/(gen_pTL->at(refIndG)));
	     
	                InsideJet++;
                }
              }
              if(!skipJet)continue;
              Bool_t skipGenJet = kFALSE;
	            if(genJetPt->size() > 0){
	              for(Int_t j2 = 0; j2 < genJetPt->size(); j2++){
                  if(skipGenJet)break;
	                if(genJetPt->at(j2) < 550)continue; // selecting jet with pT > 550 GeV
	                if(abs(genJetEta->at(j2)) > (2.4-JetR))continue; // selecting jet with the eta acceptance < |2.4-JetR|
	                Double_t DelRV0Jet_G = dR(genJetEta->at(j2),genJetPhi->at(j2),gen_etaL->at(refIndG),gen_phiL->at(refIndG));

          	      //hDelR_GM_in_wGJ->Fill(DelRV0Jet_G);
	      
	                if(DelRV0Jet_G <= JetR){
                		hCount_GenJetwMV0->Fill(1.0);
		
		                hPtL_GM_in_wGJ->Fill(gen_pTL->at(refIndG));
                		hMassL_GM_in_wGJ->Fill(cand_massL->at(refIndR));
  		              hDelR_GM_in_wGJ_wMV0->Fill(DelRV0Jet_G);

                		if(genChMult->at(j2) < 20)hPtL_GM_in_wGJ_NchLT20->Fill(gen_pTL->at(refIndG));
	  	              if(genChMult->at(j2) > 80)hPtL_GM_in_wGJ_NchGT80->Fill(gen_pTL->at(refIndG));
	                	h2PtL_GM_in_wGJ->Fill(genChMult->at(j2),gen_pTL->at(refIndG));

	              	  //if(RefRMJetEta == -999. || RefRMJetPhi == -999.)continue;
              		  Double_t DelR_Jet_RG_wMV0 = dR(genJetEta->at(j2),genJetPhi->at(j2),RefRMJetEta,RefRMJetPhi);
	              	  hDelR_RGM_wMV0->Fill(DelR_Jet_RG_wMV0);
                    skipGenJet = kTRUE;                    
      	          }
	                else{
	            	    hCount_noGenJetwMV0->Fill(1.0);
	                }
      	        }
	            }
	            else{
	              hCount_noGenJet->Fill(1.0);
	              //cout << "No gen jet found!" << endl;
	            }
            }
          }
        }

        if(genJetPt->size() > 0 && jetPt->size() > 0){
	        for(Int_t jRec = 0; jRec < jetPt->size(); jRec++){
       	    if(jetPt->at(jRec) < 550)continue;	    
            if(abs(jetEta->at(jRec)) > (2.4-JetR))continue;
	          
            for(Int_t j2 = 0; j2 < genJetPt->size(); j2++){
	            if(abs(genJetEta->at(j2)) > (2.4-JetR))continue;

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

        if(gen_massL->size() > 0 && genJetPt->size() > 0){ // there is at least one candidate and one jet in gen level
	        for(Int_t jG = 0; jG < gen_massL->size(); jG++){ // loop over gen cand
	          h1pdgId_G->Fill(gen_pdgIdL->at(jG));
	          if(gen_pdgIdL->at(jG) != cand_ID[iCand])continue;
	          if(!IsSatisfyKineCuts(gen_pTL->at(jG),gen_etaL->at(jG)))continue;
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

    }

    cout << "No of candidates: " << counter << endl;
    stopwatch.Stop();
    cout << "Real time: " << stopwatch.RealTime() << "\t CPU time: " << stopwatch.CpuTime() << endl;
    JetTree.ResetBranchAddresses();
    CandTree.ResetBranchAddresses();

    if(IsSave){
    TFile *fO = new TFile(Form("%s.root",OutputFileName.Data()),"RECREATE");
       
    hzVtxL->Write("bestVtxZL");
    hLMass->Write("InvMassL");
    
    hpdgIdL_in->Write();
    h2LMassPdg_cand->Write("h2LMassPdg_cand");
    h2LMass_matched->Write();
    hLMass_cand->Write("InvMassL_cand");
    hLMass_in->Write("InvMassL_in");
    
    hjetN->Write();
    hjetEta->Write();
    hjetPhi->Write();
    hjetPt->Write();
    
    hjetEta_wV0->Write();
    hjetPhi_wV0->Write();
    hjetPt_wV0->Write();    
    
    hjetNumDaughters->Write();
    hchargedMult->Write();
    hchargedMult_wV0->Write();

    h2chargedMultCandMinus2->Write();    
    h2chargedMultCandMinus->Write();
    h2chargedMultCand->Write();
    
    //hDelRV0Jet_in->Write();
    //hDelRV0Jet_in_PtLead->Write();
    //hDelRV0Jet_out->Write();
    //h3DelRV0Jet_in->Write();
    //h2DelRV0Jet_out->Write();
    
    hPtL_in->Write();
    
    hEtaL->Write();
    hEtaL_in->Write();
    
    hPhiL->Write();
    hPhiL_in->Write();
    
    //hInsideJet->Write();    
    //hDelRV0MJet->Write();  
    
    h3LMass_in->Write();        

    h1DelR_RG->Write();
    h1FracPt_RG->Write();
    h1DelR_RG_matched_inJ->Write();
    h1FracPt_RG_matched_inJ->Write();
    h1DelR_RG_matched_wJ->Write();
    h1FracPt_RG_matched_wJ->Write();
    h1DelR_RG_matched_woJ->Write();
    h1FracPt_RG_matched_woJ->Write();
    h1DelR_RG_closest_woJ->Write();
    h1FracPt_RG_closest_woJ->Write();
    h1DelR_RG_closest_wJ->Write();
    h1FracPt_RG_closest_wJ->Write();

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
    
    h1Pt_GM_woJ->Write();
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