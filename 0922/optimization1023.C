#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "isPassZee.h"
#include "isPassZmumu.h"


using namespace std;



TTree *tree;
TFile *f;

void optimization1023(){

  bool isSignal=0,isBG=0,isDYHT=0;
//string  masspoint[4]={"1000","2000","3000","4500"};
string  masspoint[13]={"600","800","1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};

TH1F* th1[6];
  th1[0] =new TH1F("FATjetPRmass","FATjetPRmass",600,0,600); 
  th1[1] =new TH1F("FATjetSDmass","FATjetSDmass",600,0,600); 
  th1[2] =new TH1F("FATjetPRmassL2L3Corr","FATjetPRmassL2L3Corr",600,0,600); 
  th1[3] =new TH1F("FATjetSDmassPuppiL2L3Corr","FATjetSDmassPuppiL2L3Corr",600,0,600); 
  th1[4] =new TH1F("FATjetPRmassPuppiL2L3Corr","FATjetPRmassPuppiL2L3Corr",600,0,600); 
  th1[5] =new TH1F("FATjetATLASmassL2L3Corr","FATjetATLASmassL2L3Corr",600,0,600); 
  
  for(int i=0;i<6;i++)th1[i]->Sumw2();

  for (int massP=0;massP<13;massP++){
    double scaleNTotal=0;
    
    //DY100-200
    // double scaleF=0,xsecF=139.4;isDYHT=1;
    // TString endfix =Form("DYHT100-%s",masspoint[massP].data());;
    // for(int w=1;w<90;w++){
    //   f = TFile::Open(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root",w));if (!f || !f->IsOpen())continue;
    //   TDirectory * dir = (TDirectory*)f->Get(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root:/tree",w));   dir->GetObject("treeMaker",tree);

    //DY200-400
 // double scaleF=0,xsecF=42.75;isDYHT=1;
 //     TString endfix =Form("DYHT200-%s",masspoint[massP].data());;
 //     for(int w=1;w<45;w++){
 //       f = TFile::Open(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root",w));if (!f || !f->IsOpen())continue;
 //       TDirectory * dir = (TDirectory*)f->Get(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root:/tree",w));   dir->GetObject("treeMaker",tree);


    //DY400-600
    // double scaleF=0,xsecF=5.497;isDYHT=1;
    // TString endfix =Form("DYHT400-%s",masspoint[massP].data());;
    // for(int w=1;w<45;w++){
    //    f = TFile::Open(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root",w));if (!f || !f->IsOpen())continue;
    //    TDirectory * dir = (TDirectory*)f->Get(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root:/tree",w));   dir->GetObject("treeMaker",tree);


    //DY600-inf
     // double scaleF=0,xsecF=2.21;isDYHT=1;
     // TString endfix =Form("DYHT600-%s",masspoint[massP].data());;
     // for(int w=1;w<48;w++){
     //    f = TFile::Open(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root",w));if (!f || !f->IsOpen())continue;
     //    TDirectory * dir = (TDirectory*)f->Get(Form("/data7/syu/jet_CMSSW747/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_25ns/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/results/NCUGlobalTuples_%d.root:/tree",w));   dir->GetObject("treeMaker",tree);



  //Signal
     isSignal=1; double scaleF=0; double xsecF=1;
     TString endfix =Form("signal-%s",masspoint[massP].data());
       for(int w=1;w<2;w++){
     f = TFile::Open(Form("/data7/syu/jet_CMSSW747/ZprimeToZhToZlephbb_25ns/ZprimeToZhToZlephbb_narrow_M-%s_13TeV-madgraph.root",masspoint[massP].data()));    if (!f || !f->IsOpen())continue;
        TDirectory * dir = (TDirectory*)f->Get(Form("/data7/syu/jet_CMSSW747/ZprimeToZhToZlephbb_25ns/ZprimeToZhToZlephbb_narrow_M-%s_13TeV-madgraph.root:/tree",masspoint[massP].data()));    dir->GetObject("treeMaker",tree);

      cout<<"massP="<<massP<<",w="<<w<<endl;

  TreeReader data(tree);
  //data.Print();
  Long64_t nLepton[3]={0};
  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  ofstream fout;
  fout.open("ele_Eiko.txt");

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal ++;
    
    Int_t nGenPar        = data.GetInt("nGenPar");
    Int_t* genParId      = data.GetPtrInt("genParId");
    Int_t* genParSt      = data.GetPtrInt("genParSt");
    Int_t* genMomParId   = data.GetPtrInt("genMomParId");
    Int_t* genDa1      = data.GetPtrInt("genDa1");
    Int_t* genDa2      = data.GetPtrInt("genDa2");
    TClonesArray* genParP4 = (TClonesArray*) data.GetPtrTObject("genParP4");

    bool hasElectron=false;

    for(int ig=0; ig < nGenPar; ig++){

        int pid = abs(genParId[ig]);
        if(pid!=11)continue;
        int momId = abs(genMomParId[ig]);
        if(
            momId!=23 &&
            momId!=9000001 &&
            momId!=pid)
        continue;
        hasElectron=true;
        if(hasElectron)break;

    }

    if(!hasElectron)continue;
    nPass[0]++;


  
    TLorentzVector* thisB ,* thatB;
    bool isThisB=0,isThatB=0;
    vector<int> thisBArray;

    for(int ig=0; ig < nGenPar; ig++){
      if(isThisB && isThatB)break;
        int pid = abs(genParId[ig]);
        int momId = abs(genMomParId[ig]);
        if(
            momId!=25 &&
            momId!=9000001 &&
            momId!=pid)
        continue;
	if (pid==5)thisBArray.push_back(ig);
        if (pid==5 && !isThisB){
	  thisB= (TLorentzVector*)genParP4->At(ig);
	  isThisB=1;
	  continue;
	}
	if(pid==5 && isThisB && !isThatB){
	  thatB= (TLorentzVector*)genParP4->At(ig);
          isThatB=1;
	}
	

    }
    //if (thisBArray.size()!=2)cout<<"thisBArray.size()="<<thisBArray.size()<<" , entry="<<jEntry<<endl;

     
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
 	std::string thisTrig= trigName[it];
 	bool results = trigResult[it];

	// std::cout << thisTrig << " : " << results << std::endl;
	
 	if( (thisTrig.find("HLT_Ele105")!= std::string::npos && results==1)
	    ||
	    (thisTrig.find("HLT_Mu45")!= std::string::npos && results==1)
	    )
 	  {
 	    passTrigger=true;
 	    break;
 	  }


      }


    if(!passTrigger)continue;
    nPass[2]++;

    

    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[3]++;

    Int_t nEle         = data.GetInt("nEle");
    Int_t run          = data.GetInt("runId");
    Int_t lumi         = data.GetInt("lumiSection");
    Int_t event        = data.GetInt("eventId");
    vector<bool> &passHEEPID = *((vector<bool>*) data.GetPtr("eleIsPassHEEPNoIso"));
    TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
    Float_t* eleSCEta         = data.GetPtrFloat("eleScEta");
    Float_t* eleMiniIso       = data.GetPtrFloat("eleMiniIso");
    Int_t*   eleCharge        = data.GetPtrInt("eleCharge");
	TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");

    bool findEPair=false;
    TLorentzVector l4_Z(0,0,0,0);
    std::vector<int> myElectrons;

    // select_electrons(data, myElectrons);

    // select good electrons

    for(int ie=0; ie< nEle; ie++)
      {

    	TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(ie);

    	if(fabs(thisEle->Eta())>2.5)continue;

    	if(! (fabs(eleSCEta[ie])<1.442 || fabs(eleSCEta[ie])>1.566))continue;
    	
    	if(thisEle->Pt() < 35)continue;

    	if(!passHEEPID[ie])continue;
    	
    	if(eleMiniIso[ie]>0.1)continue;

    	myElectrons.push_back(ie);
      }

	vector<int> Zee;
	TLorentzVector * checkThisEle,* checkThatEle;
	findEPair=isPassZee(data,Zee);
	//out<<jEntry<<"=";
	if(findEPair){
	checkThisEle=(TLorentzVector*)eleP4->At(Zee[0]), 
	checkThatEle=(TLorentzVector*)eleP4->At(Zee[1]);	
	}
	
	if(!findEPair){
		
		findEPair=isPassZmumu(data,Zee);
		if(findEPair){
			checkThisEle=(TLorentzVector*)muP4->At(Zee[0]), 
			checkThatEle=(TLorentzVector*)muP4->At(Zee[1]);		
		}
		
	}
	
	if(!findEPair)continue;
	l4_Z=(*checkThisEle+*checkThatEle);
	
    /*
    for(unsigned int i=0; i< myElectrons.size(); i++)
      {
	int ie = myElectrons[i];
	TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(ie);
	if(thisEle->Pt()<115)continue;
		
	for(unsigned int j=0; j< i; j++)
	  {
	    int je= myElectrons[j];

	    if(eleCharge[ie]*eleCharge[je]>0)continue;

	    TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(je);

	    Float_t mll  = (*thisEle+*thatEle).M();
	    Float_t ptll = (*thisEle+*thatEle).Pt();
	    

	    if(mll<70 || mll>110)continue;
	    if(ptll<200)continue;

	    if(!findEPair){
	      l4_Z=(*thisEle+*thatEle);
	      checkThisEle = thisEle;
              checkThatEle = thatEle;
	    }
	    findEPair=true;
	  }	
      }

    if(!findEPair)
      continue;
    */
	nPass[4]++;
     
    //cout<<"a"<<jEntry<<endl;

    Int_t nJet         = data.GetInt("FATnJet");
    TClonesArray* jetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  jetSDmass = data.GetPtrFloat("FATjetSDmass");
	Float_t*  jetPRmass = data.GetPtrFloat("FATjetPRmass");
	Float_t*  jetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
	Float_t*  jetSDmassPuppiL2L3Corr = data.GetPtrFloat("FATjetSDmassPuppiL2L3Corr");
	Float_t*  jetPRmassPuppiL2L3Corr = data.GetPtrFloat("FATjetPRmassPuppiL2L3Corr");
	Float_t*  jetATLASmassL2L3Corr = data.GetPtrFloat("FATjetATLASmassL2L3Corr");

    TLorentzVector l4_leadingJet(0,0,0,0);
    bool findAJet=false;
    bool goodJet=0;
    
	


    for(int ij=0; ij<nJet; ij++)
      {
	
	TLorentzVector* thisJet = (TLorentzVector*)jetP4->At(ij);

	//if(jetSDmass[ij]<50 || jetSDmass[ij]>110)continue;
	//cout<<thisJet->DeltaR(*checkThisEle)<<endl;
        if (thisJet->DeltaR(*checkThisEle)<0.8)continue;
        if (thisJet->DeltaR(*checkThatEle)<0.8)continue;
        if (thisJet->Pt()<200 || thisJet->Eta()>2.4)continue;

        bool checkMyEle=0;
        for (unsigned int k=0;k<myElectrons.size();k++){
          TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(myElectrons[k]);
		  if(thisEle->Pt()<115)continue;
	  if (thisJet->DeltaR(*thisEle)<0.8)checkMyEle=1;
	}
	if (checkMyEle)continue;
        if (jetSDmass[ij]<20||jetSDmass[ij]>220)continue;
		if (jetPRmass[ij]<=-999)continue;
		if (jetPRmassL2L3Corr[ij]<=-999)continue;
		if (jetSDmassPuppiL2L3Corr[ij]<=-999)continue;
		if (jetPRmassPuppiL2L3Corr[ij]<=-999)continue;
		if (jetATLASmassL2L3Corr[ij]<=-999)continue;
	
	  TLorentzVector* thisGraviton = (TLorentzVector*)jetP4->At(ij);
	  *thisGraviton=*thisGraviton+l4_Z;
	  //cout<<<<endl;
	  int temp ;
	  stringstream convert(masspoint[massP]);
	  convert >> temp;
	  //istringstream (masspoint[massP] ) >> temp;
	  double temp2=temp*0.85;
	  //cout<< temp2<<endl;
	  if(thisGraviton->M()<temp2 )continue;
	
	if(!findAJet){
	  l4_leadingJet = *thisJet;
          th1[0]->Fill(jetPRmass[ij]);
		  th1[1]->Fill(jetSDmass[ij]);
		  th1[2]->Fill(jetPRmassL2L3Corr[ij]);
		  th1[3]->Fill(jetSDmassPuppiL2L3Corr[ij]);
		  th1[4]->Fill(jetPRmassPuppiL2L3Corr[ij]);
		  th1[5]->Fill(jetATLASmassL2L3Corr[ij]);
          
	}
	//cout<<"b"<<endl;
	findAJet=true;
      }
    
    
    if(!findAJet)
      continue;
    nPass[5]++;
    //if(!goodJet)l4_leadingJet.Print();



    Float_t MGrav = (l4_leadingJet + l4_Z).M();
    if(MGrav<400)continue;
    nPass[6]++;

    fout << run << " " << lumi << " " << event << endl;
    

  } // end of loop over entries

  fout.close();
  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0; i<3; i++)
    std::cout << "nLepton[" << i << "]    = " << nLepton[i] << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;
  //th1->Draw();
  scaleF+=nPass[0];
  scaleNTotal+=nTotal;
  }
  TFile* outFile = new TFile(Form("root_files_1023/%s.root",endfix.Data()),"recreate");       
  
  
  if(!isSignal){
    
    for(int i=0;i<6;i++)th1[i]->Scale(5000*xsecF*1.23/scaleNTotal);
  }
  for(int i=0;i<6;i++)th1[i]->Write();
  
  outFile->Close();
  for(int i=0;i<6;i++)th1[i]->Reset();
  


  
  }

    }//sig mass point



