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
#include "TGaxis.h"
#include "TPad.h"
#include <fstream>

#define  nWidth 9
#define  nBmin 5

double  getEff(TH1F *h,int bmin,int bmax);
double getSign(double d1,double d2);
double  getErr(TH1F *h,int bmin,int bmax);
double  getErr2(TH1F *h,int bmin,int bmax);

using namespace std;

TFile *f,*f2,*f3;
TCanvas* c1;

void setLeg(TLegend *leg){
leg->SetFillColor(18);
leg->SetFillStyle(0);
leg->SetTextSize(0.03);
leg->SetBorderSize(2);
}
/*
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
*/
void makeEffPRN(){
  c1 = new TCanvas("c1","",1360,768);
  string output="signPRNv1BG";
  int twikiSign[13][nWidth][nBmin];
  int twikiSignNum[13][nWidth][nBmin];
  double twikiSignValue[13][nWidth][nBmin];
  double twikiSignNumValue[13][nWidth][nBmin];
  double twikiEffValue[13][nWidth][nBmin];
  //double twikiEffNumValue[13][nWidth][nBmin];
  //double twikiWidth[13][nWidth][nBmin];
  //double twikiBmin[13][nWidth][nBmin];
  string  masspoint[13]={"600","800","1000","1200","1400","1600","1800","2000","2500","3000","3500","4000","4500"};
  int width [nWidth]={20,25,30,35,40,45,50,55,60};
  int bmin[nBmin]={90,95,100,105,110};
  for (int massP=0;massP<13;massP++){
  TString fin =Form("root_files/signal-%s.root",masspoint[massP].data()),
    fin2 = Form("root_files/BulkGravitonZlepZqq-%s.root",masspoint[massP].data());
  // fin2=Form("root_files/DYHTall-%s.root",masspoint[massP].data());
  f= TFile::Open(fin.Data());
  f2= TFile::Open(fin2.Data());
  TH1F * th1 = (TH1F*)f->FindObjectAny("HMass");
  TH1F * th2 = (TH1F*)f2->FindObjectAny("HMass");
  
  TString fin3 = Form("root_files/BulkGravitonZlepZqq-%s.root",masspoint[massP].data());
  f3= TFile::Open(fin3.Data());
  TH1F * th3 = (TH1F*)f3->FindObjectAny("HMass");
  //cout<<"before"<<th2->GetEntries()<<endl;
  th2->Sumw2();
  //th2->Add(th3);
  //cout<<"after"<<th2->GetEntries()<<endl;
  double eff[nWidth][nBmin],eff2[nWidth][nBmin],err[nWidth][nBmin],err2[nWidth][nBmin],err2Num[nWidth][nBmin];
  double sign[nWidth][nBmin],signCP[nWidth][nBmin],signErr[nWidth][nBmin];
  double signNum[nWidth][nBmin],signNumCP[nWidth][nBmin],signNumErr[nWidth][nBmin];
  for(int i=0;i<nWidth;i++){
    for(int j=0;j<nBmin;j++){
      eff[i][j]=getEff(th1,bmin[j],bmin[j]+width[i]);
      eff2[i][j]=getEff(th2,bmin[j],bmin[j]+width[i]);
	  //if(eff2[i][j]==0)cout<<massP<<","<<eff2[i][j]<<endl;
	  err[i][j]=getErr(th1,bmin[j],bmin[j]+width[i]);
      err2[i][j]=getErr(th2,bmin[j],bmin[j]+width[i]);
	  err2Num[i][j]=getErr2(th2,bmin[j],bmin[j]+width[i]);
	  //if(massP==0)cout<<"method1="<<err2[i][j]<<",method2="<<err2Num[i][j]/th2->Integral()<<endl;
      sign[i][j]=getSign(eff[i][j],eff2[i][j]);
	  signNum[i][j]=getSign(eff[i][j],eff2[i][j]*th2->Integral());
      signCP[i][j]=sign[i][j];
      signNumCP[i][j]=signNum[i][j];
	  signErr[i][j]=(1/(1+sqrt(eff2[i][j])))*sqrt(err[i][j]*err[i][j]+ (eff[i][j]*eff[i][j]*err2[i][j]*err2[i][j])/(4*eff2[i][j]*(1+sqrt(eff2[i][j]))*(1+sqrt(eff2[i][j]))));
      //signNumErr[i][j]=(1/(1+sqrt(eff2[i][j]*th2->Integral())))*sqrt(err[i][j]*err[i][j]+ (eff[i][j]*eff[i][j]*err2[i][j]*err2[i][j]*th2->Integral())/(4*eff2[i][j]*(1+sqrt(eff2[i][j]*th2->Integral()))*(1+sqrt(eff2[i][j]*th2->Integral()))));
	  signNumErr[i][j]=(1/(1+sqrt(eff2[i][j]*th2->Integral())))*sqrt(err[i][j]*err[i][j]+ (eff[i][j]*eff[i][j]*err2Num[i][j]*err2Num[i][j])/(4*eff2[i][j]*th2->Integral()*(1+sqrt(eff2[i][j]*th2->Integral()))*(1+sqrt(eff2[i][j]*th2->Integral()))));
		if(eff2[i][j]<1e-7){
			signErr[i][j]=err[i][j]/(1+sqrt(eff2[i][j]));
			signNumErr[i][j]=err[i][j]/(1+sqrt(eff2[i][j]*th2->Integral()));
		}
	  //cout<<"range={"<<bmin[j]<<","<<bmin[j]+width[i]<<"} and efficiency="<<eff[i][j]<<endl;
      //cout<<"range=["<<bmin[j]<<","<<bmin[j]+width[i]<<"] and significant="<<sign[i][j]<<endl;
	 //if(massP==0)cout<<eff[i][j]<<","<<eff2[i][j]<<endl;
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

    //cout<<signCP[tempI][tempJ]<<endl;
    signCP[tempI][tempJ]=0;
	if(width[tempI]+bmin[tempJ]>150){
		
		continue;
	}
	//if(massP==0)cout<<"tempI"<<tempI<<",tempJ"<<tempJ<<endl;
	signI.push_back(tempI);
    signJ.push_back(tempJ);
  }
  
  vector<int> signINum,signJNum;
  for(int ij=0;ij<nWidth*nBmin;ij++){
    double tempSign=0;
    int tempI=0,tempJ=0;
    for(int i=0;i<nWidth;i++){
      for(int j=0;j<nBmin;j++){
        if(signNumCP[i][j]>tempSign){
          //cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
          tempSign=signNumCP[i][j];
	  tempI=i;
	  tempJ=j;
	}
      }
    }
    
    //cout<<signCP[tempI][tempJ]<<endl;
    signNumCP[tempI][tempJ]=0;
    if(width[tempI]+bmin[tempJ]>150){
		//cout<<"!";
		continue;
	}
	signINum.push_back(tempI);
    signJNum.push_back(tempJ);
  }
  
  int numberBin=15;
  int numberBinFull=35;
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
  
  //cout<<"signI.size="<<signI.size()<<endl;
  for(unsigned int i=0;i<signI.size();i++){
	  twikiSign[massP][signI[i]][signJ[i]]=i+1;
	  
	  th1SignFull->SetBinContent(i+1,sign[signI[i]][signJ[i]]);
	  th1SignFull->SetBinError(i+1,signErr[signI[i]][signJ[i]]);
	  th1EffFull->SetBinContent(i+1,eff[signI[i]][signJ[i]]);
	  th1Eff2Full->SetBinContent(i+1,eff2[signI[i]][signJ[i]]);
	  th1EffFull->SetBinError(i+1,err[signI[i]][signJ[i]]);
	  th1Eff2Full->SetBinError(i+1,err2[signI[i]][signJ[i]]);
      //th1SignFull->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJ[i]],bmin[signJ[i]]+width[signI[i]]));
      th1SignFull->GetXaxis()->SetBinLabel(i+1,"");
	  twikiSignValue[massP][signI[i]][signJ[i]]=sign[signI[i]][signJ[i]];
	  twikiEffValue[massP][signI[i]][signJ[i]]=eff[signI[i]][signJ[i]];
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  gStyle->SetOptStat(0000000000);
  double x1NDC = 0.6522;
  double y1NDC = 0.7835;
  double x2NDC = 0.9122;
  double y2NDC = 0.9860;
  
  th1Sign->SetTitle(Form("largest 15 Sign.(eff),%sGev",masspoint[massP].data()));
  th1Sign->SetMinimum(0);
  th1Sign->SetMaximum(1);
  th1Sign->SetLineColor(4);
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
  
  if(massP==0)c1->Print(Form("pdfPR/%s.pdf(",output.data()));
  //else if(massP==12)c1->Print("pdf/signv2.pdf)");
  else c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  th1SignFull->SetTitle(Form("all windows(eff),%sGev",masspoint[massP].data()));
  th1SignFull->SetMinimum(0);
  th1SignFull->SetMaximum(1);
  th1SignFull->SetLineColor(4);
  th1SignFull->Draw("Hist");
  th1EffFull->SetLineColor(2);
  th1Eff2Full->SetLineColor(3);
  th1EffFull->Draw("same");
  th1Eff2Full->Draw("same");
  leg->Draw("same");
  
  c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  for(int i=0;i<numberBin;i++){
	  th1Sign->SetBinContent(i+1,signNum[signINum[i]][signJNum[i]]);
	  th1Sign->SetBinError(i+1,signNumErr[signINum[i]][signJNum[i]]);
	  th1Eff->SetBinContent(i+1,eff[signINum[i]][signJNum[i]]);
	  th1Eff2->SetBinContent(i+1,eff2[signINum[i]][signJNum[i]]);
	  th1Eff->SetBinError(i+1,err[signINum[i]][signJNum[i]]);
	  th1Eff2->SetBinError(i+1,err2[signINum[i]][signJNum[i]]);
      th1Eff->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJNum[i]],bmin[signJNum[i]]+width[signINum[i]]));
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  
  
  for(unsigned int i=0;i<signINum.size();i++){
	  twikiSignNum[massP][signINum[i]][signJNum[i]]=i+1;
	  
	  th1SignFull->SetBinContent(i+1,signNum[signINum[i]][signJNum[i]]);
	  th1SignFull->SetBinError(i+1,signNumErr[signINum[i]][signJNum[i]]);
	  th1EffFull->SetBinContent(i+1,eff[signINum[i]][signJNum[i]]);
	  th1Eff2Full->SetBinContent(i+1,eff2[signINum[i]][signJNum[i]]);
	  th1EffFull->SetBinError(i+1,err[signINum[i]][signJNum[i]]);
	  th1Eff2Full->SetBinError(i+1,err2[signINum[i]][signJNum[i]]);
      //th1SignFull->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",bmin[signJ[i]],bmin[signJ[i]]+width[signI[i]]));
      th1SignFull->GetXaxis()->SetBinLabel(i+1,"");
      th1EffFull->GetXaxis()->SetBinLabel(i+1,"");
      th1Eff2Full->GetXaxis()->SetBinLabel(i+1,"");
	  twikiSignNumValue[massP][signINum[i]][signJNum[i]]=signNum[signINum[i]][signJNum[i]];
	  //twikiEffNumValue[massP][signINum[i]][signJNum[i]]=effNum[signINum[i]][signJNum[i]];
    //cout<<"range=["<<bmin[signJ[i]]<<","<<bmin[signJ[i]]+width[signI[i]]<<"] and significant="<<sign[signI[i]][signJ[i]]<<endl;
  }
  
  th1Eff->SetMinimum(0);
  th1Eff->SetMaximum(1);
  th1EffFull->SetMinimum(0);
  th1EffFull->SetMaximum(1);
  th1Eff->SetYTitle("efficiency");
  th1EffFull->SetYTitle("efficiency");
  th1Eff->SetTitle(Form("largest 15 Sign.(Num),%sGev",masspoint[massP].data()));
  th1Eff->Draw();
  th1Eff2->Draw("same");
  
  Float_t rightmax = 2*th1Sign->GetBinContent(1);
  //cout<<rightmax<<endl;
  Float_t scale = gPad->GetUymax()/rightmax;
  //hint1->SetLineColor(kRed);
  th1Sign->Scale(scale);
  //hint1->Draw("same");
  //draw an axis on the right side
  c1->Update();
  th1Sign->Draw("same");
  //leg->Clear();
  //leg->AddEntry(th1Eff,"signal efficiency");
  //leg->AddEntry(th1Eff2,"Bkg. efficiency");
  //leg->SetY1(0.8335);
  leg->Draw("same");
  
  TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis->SetTitle("significance");
   axis->SetTitleColor(4);
   axis->SetLineColor(4);
   axis->SetLabelColor(4);
   axis->Draw();
  
  c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  
  
  
  th1EffFull->SetTitle(Form("all windows(Num),%sGev",masspoint[massP].data()));
  th1EffFull->Draw();
  th1Eff2Full->Draw("same");
  c1->Update();
  rightmax = 1.1*th1SignFull->GetBinContent(1);
  scale = gPad->GetUymax()/rightmax;
  th1SignFull->Scale(scale);
  
  th1SignFull->Draw("Hist,same");
  leg->Draw("same");
  
  TGaxis *axis2 = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
   axis2->SetTitle("significance");
   axis2->SetTitleColor(4);
   axis2->SetLineColor(4);
   axis2->SetLabelColor(4);
   axis2->Draw();
  
  if(massP==12)c1->Print(Form("pdfPR/%s.pdf)",output.data()));
  else c1->Print(Form("pdfPR/%s.pdf",output.data()));
  
  }
  ofstream myfile;
  myfile.open ("txt/twikiOPPRN.txt");
  myfile<<"|*windowRange*|";
  for (int massP=0;massP<13;massP++)myfile<<"*"<<masspoint[massP].data()<<"Num*|";
  //myfile<<"*";
  for (int massP=0;massP<13;massP++)myfile<<"*"<<masspoint[massP].data()<<"Eff*|";
  myfile<<"*avg.rank(eff)*|*avg.rank(Num)*|*avg.rank(total)*|"<<endl;
  for(int j=0;j<nBmin;j++){
	for(int i=0;i<nWidth;i++){
		if(width[i]+bmin[j]>150)continue;
		myfile<<"|"<<bmin[j]<<"to"<<bmin[j]+width[i]<<"|";
	    double temp1=0,temp2=0;
		for (int massP=0;massP<13;massP++){
			myfile<<twikiSignNum[massP][i][j]<<"("<<twikiSignNumValue[massP][i][j]<<")|";
			temp1+=twikiSign[massP][i][j];
			temp2+=twikiSignNum[massP][i][j];
		}
		for (int massP=0;massP<13;massP++){
			myfile<<twikiSign[massP][i][j]<<"("<<twikiSignValue[massP][i][j]<<")|";
		}
		myfile<<temp1/13<<"|"<<temp2/13<<"|"<<(temp1+temp2)/26<<"|"<<endl;
		//<<"to"<<bmin[i]+width[j]<<"|"<<endl;
		//cout<<j<<","<<i<<endl;
	}
  }
  myfile<<endl;
  
  myfile<<"|*rank(signalEff)*|";
  for (int massP=0;massP<13;massP++)myfile<<"*"<<masspoint[massP].data()<<"Eff*|";
  myfile<<endl;
  vector<int> signINum,signJNum;
 //for(int ij=0;ij<nWidth*nBmin;ij++){
    for(int ij=0;ij<15;ij++){
    myfile<<"|"<<ij+1<<"|";
	for (int massP=0;massP<13;massP++){
		double tempSign=1000;
    	int tempI=0,tempJ=0;
		//for(int i=0;i<nWidth;i++)for(int j=0;j<nBmin;j++)if(width[tempI]+bmin[tempJ]>150)twikiSign[massP][tempI][tempJ]=10000;
		
		for(int i=0;i<nWidth;i++){
			for(int j=0;j<nBmin;j++){
				if(twikiSign[massP][i][j]<tempSign){
					if(width[i]+bmin[j]>150)continue;
					//cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
					tempSign=twikiSign[massP][i][j];
					tempI=i;
					tempJ=j;
					//cout<<tempI<<","<<tempJ<<","<<tempSign<<endl;
				}
			}
		}
		//signINum.push_back(tempI);
		//signJNum.push_back(tempJ);
		//cout<<signCP[tempI][tempJ]<<endl;
	//cout<<"twikiSign[massP][tempI][tempJ]="<<twikiSign[massP][tempI][tempJ]<<"|"<<bmin[tempJ]<<"to"<<bmin[tempJ]+width[tempI]<<"|"<<endl;
		twikiSign[massP][tempI][tempJ]=10000;
		//if(width[tempI]+bmin[tempJ]>150)continue;
		myfile<<bmin[tempJ]<<"to"<<bmin[tempJ]+width[tempI]<<"("<<setprecision(3)<<twikiEffValue[massP][tempI][tempJ]<<")|";
		//cout<<twikiSign[massP][tempI][tempJ];
		
	}
	myfile<<endl;
  }
  myfile<<endl;
  
  myfile<<"|*rank(signalEff)*|";
  for (int massP=0;massP<13;massP++)myfile<<"*"<<masspoint[massP].data()<<"Num*|";
  myfile<<endl;
  //vector<int> signINum,signJNum;
 //for(int ij=0;ij<nWidth*nBmin;ij++){
    for(int ij=0;ij<15;ij++){
    myfile<<"|"<<ij+1<<"|";
	for (int massP=0;massP<13;massP++){
		double tempSign=10000;
    	int tempI=0,tempJ=0;
		//for(int i=0;i<nWidth;i++)for(int j=0;j<nBmin;j++)if(width[tempI]+bmin[tempJ]>150)twikiSignNum[massP][tempI][tempJ]=10000;
		
		for(int i=0;i<nWidth;i++){
			for(int j=0;j<nBmin;j++){
				if(twikiSignNum[massP][i][j]<tempSign){
					if(width[i]+bmin[j]>150)continue;
					//cout<<"i="<<i<<",j="<<j<<",signcp="<<signCP[i][j]<<",temp="<<tempSign<<endl;
					tempSign=twikiSignNum[massP][i][j];
					tempI=i;
					tempJ=j;
					//cout<<tempI<<","<<tempJ<<","<<tempSign<<endl;
				}
			}
		}
		//signINum.push_back(tempI);
		//signJNum.push_back(tempJ);
		//cout<<signCP[tempI][tempJ]<<endl;
		twikiSignNum[massP][tempI][tempJ]=10000;
		//if(width[tempI]+bmin[tempJ]>150)continue;
		myfile<<bmin[tempJ]<<"to"<<bmin[tempJ]+width[tempI]<<"("<<setprecision(3)<<twikiEffValue[massP][tempI][tempJ]<<")|";
		//cout<<twikiSign[massP][tempI][tempJ];
		
	}
	myfile<<endl;
  }
  myfile<<endl;
}

double  getEff(TH1F *h,int bmin,int bmax){
	if(h->Integral()==0){
		cout<<h->Integral();
		return 0;
	}
  double denom=h->Integral();//cout<<"de="<<denom;
  double nomin=h->Integral(bmin,bmax);//cout<<"no="<<nomin<<endl;

  
  return nomin/denom;
}

double  getErr(TH1F *h,int bmin,int bmax){
	if(h->Integral()==0)return 0;
  double temp=0,temp2=0;
  double denom=h->Integral();
  //double nomin=h->Integral(bmin,bmax);
  for(int i=1;i<h->GetNbinsX();i++){
		if(  i>bmin && i<bmax+1){
			temp+=h->GetBinError(i)*h->GetBinError(i);
		}
		else 
			temp2+=h->GetBinError(i)*h->GetBinError(i);
	}
  
  temp=temp/(h->Integral()*h->Integral())-2*h->Integral(bmin,bmax)*temp/(h->Integral()*h->Integral()*h->Integral())+h->Integral(bmin,bmax)*h->Integral(bmin,bmax)*(temp+temp2)/(h->Integral()*h->Integral()*h->Integral()*h->Integral());
  
  return sqrt(temp);
}

double  getErr2(TH1F *h,int bmin,int bmax){
  double temp=0;
  double denom=h->Integral();
  //double nomin=h->Integral(bmin,bmax);
  for(int i=bmin;i<bmax+1;i++){
	  temp+=h->GetBinError(i)*h->GetBinError(i);
  }
  temp=sqrt(temp);
  return temp;
}

double getSign(double d1,double d2){
  return d1/(1+sqrt(d2));
}
