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
TCanvas* c1;


void setLeg(TLegend *leg){
leg->SetFillColor(18);
leg->SetFillStyle(0);
leg->SetTextSize(0.03);
leg->SetBorderSize(2);
}

void drawSameOP(){
	c1 = new TCanvas("c1","",1360,768);
	string  masspoint[13]={"600","800","1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
  for (int massP=0;massP<1;massP++){
	  
	  TString fin =Form("root_files/signal-%s.root",masspoint[massP].data()),
      fin2 = Form("root_files/BulkGravitonZlepZqq-%s.root",masspoint[massP].data()),
      fin3=Form("root_files/DYtest-%s.root",masspoint[massP].data());
	  f= TFile::Open(fin.Data());
      f2= TFile::Open(fin2.Data());
	  f3= TFile::Open(fin3.Data());
	  TH1F * th1 = (TH1F*)f->FindObjectAny("HMass");
      TH1F * th2 = (TH1F*)f3->FindObjectAny("HMass");
	  TH1F * th3 = (TH1F*)f2->FindObjectAny("HMass");
	  
	  gStyle->SetOptStat(0000000000);
	  double x1NDC = 0.6522;
      double y1NDC = 0.7835;
      double x2NDC = 0.9122;
      double y2NDC = 0.9860;
	  
      
	  th1->SetTitle(Form("Z' mass %s Gev",masspoint[massP].data()));
	  
	  
	  th1->Sumw2();
	  th1->Scale(1/th1->Integral());
	  th2->Sumw2();
	  th2->Scale(1/th2->Integral());
	  th3->Sumw2();
	  th3->Scale(1/th3->Integral());
	  th1->SetYTitle("normalized to 1");
	  th1->Rebin(5);
	  th2->Rebin(5);
	  th3->Rebin(5);
	  th1->GetXaxis()->SetRangeUser(20,220);
	  th2->GetXaxis()->SetRangeUser(20,220);
	  th3->GetXaxis()->SetRangeUser(20,220);
	  
	  
	  
	  
	  th2->SetLineColor(46);
	  th3->SetLineColor(30);
	  th1->SetMaximum(th1->GetMaximum()>th3->GetMaximum()?th1->GetMaximum()*1.1:th3->GetMaximum()*1.1);
	  th1->Draw("hist,e");
	  th3->Draw("same,hist");
	  th2->Draw("same,hist");

	  
      
  
 
   TLegend* leg ;
  leg=new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  setLeg(leg);
   leg->AddEntry(th1,"signal");
   leg->AddEntry(th2,"DY");
   leg->AddEntry(th3,"BG");
   leg->Draw("same");
   //if(massP==0)c1->Print("pdf/drawSameOPv5.pdf(");
   //else if (massP==12)c1->Print("pdf/drawSameOPv5.pdf)");
   //else c1->Print("pdf/drawSameOPv5.pdf");
	  
  }
	
	
	
}