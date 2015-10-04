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
#include <string>
#include "TGaxis.h"
#include "TPad.h"

TFile *f,*f2,*f3;
TCanvas* c1,*c2;


void setLeg(TLegend *leg){
leg->SetFillColor(18);
leg->SetFillStyle(0);
leg->SetTextSize(0.02);
leg->SetBorderSize(2);
}

void drawSameOP(){
	TH1F * thn[5],* thn2[5];
	c1 = new TCanvas("c1","",1360,768);
	//string  masspoint[4]={"1000","2000","3000","4500"};
	string  masspoint[13]={"600","800","1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
  for (int massP=0;massP<13;massP++){
	  c1->cd();
	  TString fin =Form("root_files/signal-%s.root",masspoint[massP].data()),
      fin2 = Form("root_files/DY-%s.root",masspoint[massP].data());
      //fin3=Form("root_files/DYHTall-%s.root",masspoint[massP].data());
	  f= TFile::Open(fin.Data());
      f2= TFile::Open(fin2.Data());
	  //f3= TFile::Open(fin3.Data());
	  TH1F * th1[6];
      th1[0]	  = (TH1F*)f->FindObjectAny("FATjetPRmass");
	  th1[1]	  = (TH1F*)f->FindObjectAny("FATjetSDmass");
	  th1[2]	  = (TH1F*)f->FindObjectAny("FATjetPRmassL2L3Corr");
	  th1[3]	  = (TH1F*)f->FindObjectAny("FATjetSDmassPuppiL2L3Corr");
	  th1[4]	  = (TH1F*)f->FindObjectAny("FATjetPRmassPuppiL2L3Corr");
	  th1[5]	  = (TH1F*)f->FindObjectAny("FATjetATLASmassL2L3Corr");
      TH1F * th2 = (TH1F*)f2->FindObjectAny("FATjetPRmass");
	  //TH1F * th3 = (TH1F*)f2->FindObjectAny("HMass");
	  
	  gStyle->SetOptStat(0000000000);
	  double x1NDC = 0.6522;
      double y1NDC = 0.7835;
      double x2NDC = 0.9122;
      double y2NDC = 0.9860;
	  
      
	  th1[0]->SetTitle(Form("Z' mass %s Gev",masspoint[massP].data()));
	  
	  
	  for(int i=0;i<6;i++)th1[i]->Sumw2();
	  for(int i=0;i<6;i++)th1[i]->Scale(1/th1[i]->Integral());
	  th2->Sumw2();
	  th2->Scale(1/th2->Integral());
	  
	  th1[0]->SetYTitle("normalized to 1");
	  
	  
	  int bmin=0,bmax=0;

    for (int k=1;k<3000;k++){
		bmin=k;
		if (th1[0]->GetBinContent(k)!=0) break;
	}

	for (int k=3000;k>0;k--){
		bmax=k;
		if (th1[0]->GetBinContent(k)!=0) break;
	}
	if((bmin-300)<0)bmin=0;
	else bmin=bmin-300;
	//bmax=bmax+bmax/10;
	  bmax=250;
	  
	  int rebinNum=5;
	  
	  for(int i=0;i<6;i++)th1[i]->Rebin(rebinNum);
	  th2->Rebin(rebinNum);
	  
	  
	
	
	  
	  for(int i=0;i<6;i++)th1[i]->GetXaxis()->SetRangeUser(bmin,bmax);
	  th2->GetXaxis()->SetRangeUser(bmin,bmax);
	  
	  
	  
	  
	  
	  
	  th2->SetLineColor(46);
	  th2->SetLineStyle(7);
	  
	  
	  
	  //th1->SetMaximum(th1->GetMaximum()>th3->GetMaximum()?th1->GetMaximum()*1.1:th3->GetMaximum()*1.1);
	  th1[2]->SetMaximum(th1[0]->GetMaximum()>th2->GetMaximum()?th1[0]->GetMaximum()*1.2:th2->GetMaximum()*1.2);
	  th1[2]->SetXTitle("mass of AK8 jet [GeV]"); 
	  th1[2]->SetYTitle("Arbitrari Unit");
	  th1[2]->SetTitle(Form("Z' mass %s GeV",masspoint[massP].data()));
	  
	  
   if(massP==2){
	   thn[0]=th1[0];
	   thn2[0]=(TH1F*)th1[2]->Clone("");
   }
   if(massP==4){
	   thn[1]=th1[0];
	   thn2[1]=(TH1F*)th1[2]->Clone("");
   }
   if(massP==7){
	   thn[2]=th1[0];
	   thn2[2]=(TH1F*)th1[2]->Clone("");
   }
   if(massP==9){
	   thn[3]=th1[0];
	   thn2[3]=(TH1F*)th1[2]->Clone("");
   }
   if(massP==11){
	   thn[4]=th1[0];
	   thn2[4]=(TH1F*)th1[2]->Clone("");
   }
   
	  for(int i=0;i<6;i++)th1[i]->SetLineStyle(i+1);
	  for(int i=0;i<6;i++)th1[i]->SetLineColor(i+1);
	  th1[2]->SetFillStyle(1001);
	  th1[2]->SetFillColor(3);
	  c1->cd();
	  th1[2]->Draw("hist,e");
	  for(int i=0;i<6;i++){
		  if(i==2)continue;
		  th1[i]->Draw("same,hist,e");
		  
	  }
	  th2->Draw("same,hist,e");

	  
      
  
 
   TLegend* leg ;
   if (massP==0)leg=new TLegend(0.591452,0.662447,0.790645,0.883966);
  else leg=new TLegend(0.191452,0.662447,0.390645,0.883966);
  setLeg(leg);
   leg->AddEntry(th1[0],"Uncorrected pruned");
   leg->AddEntry(th1[1],"Uncorrected softdrop");
   leg->AddEntry(th1[2],"L2L3-corrected pruned");
   leg->AddEntry(th1[3],"PUPPI+L2L3-corrected pruned");
   leg->AddEntry(th1[4],"PUPPI+L2L3-corrected softdrop");
   leg->AddEntry(th1[5],"L2L3-corrected trimmed mass");
   leg->AddEntry(th2,"DY");
   
   leg->Draw("same");
   if(massP==0)c1->Print("pdf/drawSameOPv7DYHT.pdf(");
   else if (massP==12)c1->Print("pdf/drawSameOPv7DYHT.pdf)");
   else c1->Print("pdf/drawSameOPv7DYHT.pdf");
   
   
   
   
	  
  }
	
	c2 = new TCanvas("c1","",1360,768);
	c2->SetGridx(1);

  c2->SetGridy(1);
	thn[0]->SetLineColor(kRed);
	thn2[0]->SetLineColor(kRed);
	thn[1]->SetLineColor( kBlue);
	thn2[1]->SetLineColor( kBlue);
	thn[2]->SetLineColor(kViolet-3);
	thn2[2]->SetLineColor(kViolet-3);
	thn[3]->SetLineColor(kGreen+2);
	thn2[3]->SetLineColor(kGreen+2);
	thn[4]->SetLineColor(kAzure+1);
	thn2[4]->SetLineColor(kAzure+1);
	
	thn[0]->SetTitle("Uncorrected pruned");
	thn2[0]->SetTitle("L2L3-corrected pruned");
	
	TLegend* leg ;
	leg=new TLegend(0.191452,0.662447,0.390645,0.883966);
	leg->AddEntry(thn[0],"Z' mass 1000 GeV");
	leg->AddEntry(thn[1],"Z' mass 1400 GeV");
	leg->AddEntry(thn[2],"Z' mass 2000 GeV");
	leg->AddEntry(thn[3],"Z' mass 3000 GeV");
	leg->AddEntry(thn[4],"Z' mass 4000 GeV");
	
	thn[0]->SetMaximum(thn[0]->GetMaximum()*1.2);
	thn[0]->Draw("C,hist");
	for(int i=1;i<5;i++)thn[i]->Draw("same,C,hist");
	leg->Draw("same");
	c2->Print("1002/1.png");
	
	thn2[0]->SetMaximum(thn2[0]->GetMaximum()*1);
	thn2[0]->Draw("C,hist");
	for(int i=1;i<5;i++)thn2[i]->Draw("same,C,hist");
	leg->Draw("same");
	c2->Print("1002/2.png");
}