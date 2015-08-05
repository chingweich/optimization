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
#include <string>

#define  nWidth 12
#define  nBmin 41 

double  getEff(TH1F *h,int bmin,int bmax);
double getSign(double d1,double d2);
double  getErr(TH1F *h,int bmin,int bmax);

using namespace std;

TFile *f,*f2;
TCanvas* c1;

void setLeg(TLegend *leg){
leg->SetFillColor(18);
leg->SetFillStyle(0);
leg->SetTextSize(0.03);
leg->SetBorderSize(2);
}

void makeEff(TString fin){

  c1 = new TCanvas("c1","",1360,768);
  f= TFile::Open(fin.Data());
  TH1F * th1 = (TH1F*)f->FindObjectAny("HMass");
  int width [nWidth]={30,35,50,45,50};
  int bmin[nBmin]={91,93,95,97,99,101,103,105,107,109};
  double eff[nWidth][nBmin];
  for(int i=0;i<nWidth;i++){
    for(int j=0;j<nBmin;j++){
      eff[i][j]=getEff(th1,bmin[j],bmin[j]+width[i]);
      cout<<"range={"<<bmin[j]<<","<<bmin[j]+width[i]<<"} and efficiency="<<eff[i][j]<<endl;
    }
  }
}

void makeEff(){
  c1 = new TCanvas("c1","",1360,768);
  string  masspoint[13]={"600","800","1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
  for (int massP=0;massP<13;massP++){
  TString fin =Form("root_files/signal-%s.root",masspoint[massP].data())
  ,fin2 = Form("root_files/BulkGravitonZlepZqq-%s.root",masspoint[massP].data());
  f= TFile::Open(fin.Data());
  f2= TFile::Open(fin2.Data());
  TH1F * th1 = (TH1F*)f->FindObjectAny("HMass");
  TH1F * th2 = (TH1F*)f2->FindObjectAny("HMass");
  int width [nWidth]={20,25,30,35,40,45,50,55,60,65,70,75};
  int bmin[nBmin]={81,82,83,84,85,86,87,88,89,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120};
  double eff[nWidth][nBmin],eff2[nWidth][nBmin],err[nWidth][nBmin],err2[nWidth][nBmin];
  double sign[nWidth][nBmin],signCP[nWidth][nBmin],signErr[nWidth][nBmin];
  for(int i=0;i<nWidth;i++){
    for(int j=0;j<nBmin;j++){
      eff[i][j]=getEff(th1,bmin[j],bmin[j]+width[i]);
      eff2[i][j]=getEff(th2,bmin[j],bmin[j]+width[i]);
	  err[i][j]=getErr(th1,bmin[j],bmin[j]+width[i]);
      err2[i][j]=getErr(th2,bmin[j],bmin[j]+width[i]);
      sign[i][j]=getSign(eff[i][j],eff2[i][j]);
      signCP[i][j]=getSign(eff[i][j],eff2[i][j]);
	  signErr[i][j]=(1/(1+sqrt(eff2[i][j])))*sqrt(err[i][j]*err[i][j]+ (eff[i][j]*eff[i][j]*err2[i][j]*err2[i][j])/(4*eff2[i][j]*(1+sqrt(eff2[i][j]))*(1+sqrt(eff2[i][j]))));
      //cout<<"range={"<<bmin[j]<<","<<bmin[j]+width[i]<<"} and efficiency="<<eff[i][j]<<endl;
      //cout<<"range=["<<bmin[j]<<","<<bmin[j]+width[i]<<"] and significant="<<sign[i][j]<<endl;
    }
  }
  vector<int> signI,signJ;
  for(int ij=0;ij<nWidth*nBmin;ij++){
    double tempSign=0;
    int tempI=0,tempJ=0;
    for(int i=0;i<nWidth;i++){
      for(int j=0;j<nBmin;j++){
        if(signCP[i][j]>tempSign){
          //cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
          tempSign=signCP[i][j];
	  tempI=i;
	  tempJ=j;
	}
      }
    }
    signI.push_back(tempI);
    signJ.push_back(tempJ);
    //cout<<signCP[tempI][tempJ]<<endl;
    signCP[tempI][tempJ]=0;
  }
  
  int numberBin=15;
  int numberBinFull=nWidth*nBmin;
  TH1F * th1Sign=new TH1F("Sign","Sign",numberBin,0,numberBin);
  TH1F * th1Eff=new TH1F("th1Eff","th1Eff",numberBin,0,numberBin);
  TH1F * th1Eff2=new TH1F("th1Eff2","th1Eff2",numberBin,0,numberBin);
  TH1F * th1SignFull=new TH1F("SignFull","SignFull",numberBinFull,0,numberBinFull);
  TH1F * th1EffFull=new TH1F("th1EffFull","th1EffFull",numberBinFull,0,numberBinFull);
  TH1F * th1Eff2Full=new TH1F("th1Eff2Full","th1Eff2Full",numberBinFull,0,numberBinFull);
  for(int i=0;i<numberBin;i++){
	  th1Sign->SetBinContent(i+1,sign[signI[i]][signJ[i]]);
	  th1Sign->SetBinError(i+1,signErr[signI[i]][signJ[i]]);
	  th1Eff->SetBinContent(i+1,eff[signI[i]][signJ[i]]);
	  th1Eff2->SetBinContent(i+1,eff2[signI[i]][signJ[i]]);
	  th1Eff->SetBinError(i+1,err[signI[i]][signJ[i]]);
	  th1Eff2->SetBinError(i+1,err2[signI[i]][signJ[i]]);
      th1Sign->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJ[i]],bmin[signJ[i]]+width[signI[i]]));
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  
  
  for(unsigned int i=0;i<signI.size();i++){
	  th1SignFull->SetBinContent(i+1,sign[signI[i]][signJ[i]]);
	  th1SignFull->SetBinError(i+1,signErr[signI[i]][signJ[i]]);
	  th1EffFull->SetBinContent(i+1,eff[signI[i]][signJ[i]]);
	  th1Eff2Full->SetBinContent(i+1,eff2[signI[i]][signJ[i]]);
	  th1EffFull->SetBinError(i+1,err[signI[i]][signJ[i]]);
	  th1Eff2Full->SetBinError(i+1,err2[signI[i]][signJ[i]]);
      //th1SignFull->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJ[i]],bmin[signJ[i]]+width[signI[i]]));
      th1SignFull->GetXaxis()->SetBinLabel(i+1,"");
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  gStyle->SetOptStat(0000000000);
  double x1NDC = 0.6522;
  double y1NDC = 0.7835;
  double x2NDC = 0.9122;
  double y2NDC = 0.9860;
  
  th1Sign->SetTitle(Form("largest 15 Sign. windows,%sGev",masspoint[massP].data()));
  th1Sign->SetMinimum(0);
  th1Sign->SetMaximum(1);
  th1Sign->Draw();
  th1Eff->SetLineColor(2);
  th1Eff2->SetLineColor(3);
  th1Eff->Draw("same");
  th1Eff2->Draw("same");
  TLegend* leg ;
  leg=new TLegend(x1NDC,y1NDC,x2NDC,y2NDC);
  setLeg(leg);
  leg->AddEntry(th1Sign,"significance");
  leg->AddEntry(th1Eff,"signal efficiency");
  leg->AddEntry(th1Eff2,"Bkg. efficiency");
  leg->Draw("same");
  
  if(massP==0)c1->Print("pdf/signv2.pdf(");
  //else if(massP==12)c1->Print("pdf/signv2.pdf)");
  else c1->Print("pdf/signv2.pdf");
  
  th1SignFull->SetTitle(Form("%s",masspoint[massP].data()));
  th1SignFull->SetMinimum(0);
  th1SignFull->SetMaximum(1);
  th1SignFull->Draw("Hist");
  th1EffFull->SetLineColor(2);
  th1Eff2Full->SetLineColor(3);
  th1EffFull->Draw("same");
  th1Eff2Full->Draw("same");
  leg->Draw("same");
  
  if(massP==12)c1->Print("pdf/signv2.pdf)");
  else c1->Print("pdf/signv2.pdf");
  
  }

}

double  getEff(TH1F *h,int bmin,int bmax){
  double denom=h->Integral();
  double nomin=h->Integral(bmin,bmax);
  return nomin/denom;
}

double  getErr(TH1F *h,int bmin,int bmax){
  double temp=0;
  double denom=h->Integral();
  //double nomin=h->Integral(bmin,bmax);
  for(int i=bmin;i<bmax+1;i++){
	  temp+=h->GetBinError(i)*h->GetBinError(i);
  }
  temp=sqrt(temp);
  return temp/denom;
}

double getSign(double d1,double d2){
  return d1/(1+sqrt(d2));
}
