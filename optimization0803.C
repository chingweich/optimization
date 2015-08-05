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

#define  nWidth 5
#define  nBmin 10 

using namespace std;

double  getEff(TH1F *h,int bmin,int bmax);

TTree *tree;
TFile *f;

void optimization0803(){

  

  TH1F * th1 =new TH1F("HMass","HMass",600,0,600); 
  TH1F * th2 =new TH1F("Pt of good","Pt of good",2000,0,2000); 
  TH1F * th3 =new TH1F("Pt of bad","P of bad",2000,0,2000);
  TH1F * thr =new TH1F("Dr of good","deltaR of good",50,0,3);
  TH1F * thr2 =new TH1F("Dr of bad","deltaR of bad",50,0,3);

  //DY bkg
  // TString endfix ="DYBkg";
  // for(int w=1;w<363;w++){
  //    f = TFile::Open(Form("/data7/khurana/NCUGlobalTuples/SPRING15/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ForEIKO/150729_202330/0000/NCUGlobalTuples_%d.root",w));
  //     if (!f || !f->IsOpen())continue;
  //    TDirectory * dir = (TDirectory*)f->Get(Form("/data7/khurana/NCUGlobalTuples/SPRING15/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ForEIKO/150729_202330/0000/NCUGlobalTuples_%d.root:/tree",w));
  //  dir->GetObject("treeMaker",tree);

  //BulkGravitonZlepZqq  noCleaning_BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph.root
    string  masspoint[13]={"600","800","1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
   // string  crabInf[13]={"150729_212416","150729_212459","150729_211555","150729_211640","150729_211727","150729_211727","150729_211855","150729_211940","150729_212023","150729_212110","150729_212159","150729_212245","150729_212331"};
    for (int massP=0;massP<13;massP++){
    TString endfix =Form("BulkGravitonZlepZqq-%s",masspoint[massP].data());
    for(int w=1;w<2;w++){
      f = TFile::Open(Form("/data2/syu/13TeV/BulkGravitonZlepZqq/noCleaning_BulkGravToZZToZlepZhad_narrow_M-%s_13TeV-madgraph.root",masspoint[massP].data()));
      if (!f || !f->IsOpen())continue;
      TDirectory * dir = (TDirectory*)f->Get(Form("/data2/syu/13TeV/BulkGravitonZlepZqq/noCleaning_BulkGravToZZToZlepZhad_narrow_M-%s_13TeV-madgraph.root:/tree",masspoint[massP].data()));
      dir->GetObject("treeMaker",tree);


  //Signal

   // string  masspoint[13]={"600","800","1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
   // string  crabInf[13]={"150729_212416","150729_212459","150729_211555","150729_211640","150729_211727","150729_211727","150729_211855","150729_211940","150729_212023","150729_212110","150729_212159","150729_212245","150729_212331"};
   // for (int massP=0;massP<13;massP++){
   // TString endfix =Form("signal-%s",,masspoint[massP].data());
   // for(int w=1;w<3;w++){
   //   f = TFile::Open(Form("/data7/khurana/NCUGlobalTuples/SPRING15/ZPrimeSignal/ZprimeToZhToZlephbb_narrow_M-%s_13TeV-madgraph/crab_ZprimeToZhToZlephbb_narrow_M-%s_13TeV-madgraph/%s/0000/NCUGlobalTuples_%d.root",masspoint[massP].data(),masspoint[massP].data(),crabInf[massP].data(),w));
   // if (!f || !f->IsOpen())continue;
   // TDirectory * dir = (TDirectory*)f->Get(Form("/data7/khurana/NCUGlobalTuples/SPRING15/ZPrimeSignal/ZprimeToZhToZlephbb_narrow_M-%s_13TeV-madgraph/crab_ZprimeToZhToZlephbb_narrow_M-%s_13TeV-madgraph/%s/0000/NCUGlobalTuples_%d.root:/tree",masspoint[massP].data(),masspoint[massP].data(),crabInf[massP].data(),w));
   // dir->GetObject("treeMaker",tree);

  cout<<w<<endl;

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

    // bool hasHadron=false;

    // for(int ig=0; ig < nGenPar; ig++){

    //     int pid = abs(genParId[ig]);
    //     //if(pid!=11)continue;
    //     int momId = abs(genMomParId[ig]);
    //     if(
    //         momId!=25 &&
    //         momId!=9000001 &&
    //         momId!=pid)
    //     continue;
    //     hasHadron=true;
    //     if(hasHadron)break;


    // }
    


    // if(!hasHadron)continue;
    // nPass[1]++;
  
    TLorentzVector* thisB ,* thatB;
    bool isThisB=0,isThatB=0;
    vector<int> thisBArray;

    for(int ig=0; ig < nGenPar; ig++){
      if(isThisB && isThatB)continue;
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
    	
    	if(thisEle->Pt() < 115)continue;

    	if(!passHEEPID[ie])continue;
    	
    	if(eleMiniIso[ie]>0.1)continue;

    	myElectrons.push_back(ie);
      }

    TLorentzVector * checkThisEle, * checkThatEle;

    for(unsigned int i=0; i< myElectrons.size(); i++)
      {
	int ie = myElectrons[i];
	TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(ie);

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
    nPass[4]++;

    //cout<<"a"<<jEntry<<endl;

    Int_t nJet         = data.GetInt("FATnJet");
    TClonesArray* jetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  jetSDmass = data.GetPtrFloat("FATjetSDmass");

    TLorentzVector l4_leadingJet(0,0,0,0);
    bool findAJet=false;
    bool goodJet=0;



    for(int ij=0; ij<nJet; ij++)
      {
	
	TLorentzVector* thisJet = (TLorentzVector*)jetP4->At(ij);

	//if(jetSDmass[ij]<50 || jetSDmass[ij]>110)continue;
        if (thisJet->DeltaR(*checkThisEle)<0.8)continue;
        if (thisJet->DeltaR(*checkThatEle)<0.8)continue;
        if (thisJet->Pt()<20 || thisJet->Eta()>2.4)continue;
        bool checkMyEle=0;
        for (unsigned int k=0;k<myElectrons.size();k++){
          TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(myElectrons[k]);
	  if (thisJet->DeltaR(*thisEle)<0.8)checkMyEle=1;
	}
	if (checkMyEle)continue;
        if (jetSDmass[ij]<20||jetSDmass[ij]>220)continue;
	if(!findAJet){
	  l4_leadingJet = *thisJet;
          th1->Fill(jetSDmass[ij]);
          if(jetSDmass[ij]<20){
	    th3->Fill(thisJet->Pt());
	    //thr2->Fill(thisB->DeltaR(*thatB));
	  }
	  else {
	    //thr->Fill(thisB->DeltaR(*thatB));
	    th2->Fill(thisJet->Pt());
	  }
	  //if (thisJet->DeltaR(*thisB)<0.8)goodJet=1;
          //else if (thisJet->DeltaR(*thatB)<0.8)goodJet=1;
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

  }
  TFile* outFile = new TFile(Form("root_files/%s.root",endfix.Data()),"recreate");       
  
  th1->Write();
  th2->Write();
  th3->Write();
  thr->Write();
  thr2->Write();
  outFile->Close();

  int width [nWidth]={30,35,50,45,50};
  int bmin[nBmin]={91,93,95,97,99,101,103,105,107,109};
  double eff[nWidth][nBmin];
  for(int i=0;i<nWidth;i++){
    for(int j=0;j<nBmin;j++){
      //eff[i][j]=getEff(th1,bmin[j],bmin[j]+width[i]);
      //cout<<"range={"<<bmin[j]<<","<<bmin[j]+width[i]<<"} and efficiency="<<eff[i][j]<<endl;
    }
  }

  
}

   }//sig mass point


double  getEff(TH1F *h,int bmin,int bmax){
  double denom=h->Integral();
  double nomin=h->Integral(bmin,bmax);
  return nomin/denom;
}
