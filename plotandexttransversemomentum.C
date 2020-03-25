#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TRandom.h>

#include <TCanvas.h>
#include <TList.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TASImage.h>
#include <TLine.h>
#include <TBox.h>
#include <TAxis.h>
#include "esa_tools.C"
#include "DebugClassesMultESA2016.C"


#include "my_settings.C"
//#include <TGxis.h>
#include <iostream>

//--------------------------------
void CreateDir(const Char_t* dirName)
{
   TString pwd(gSystem->pwd());
   gSystem->cd(pwd.Data());

   if(gSystem->cd(dirName)) {
   gSystem->cd(pwd.Data());
   } else {
   gSystem->mkdir(dirName, kTRUE); // kTRUE means recursive
   }   
}   


//root -l plotandextlowpt("EsBin7","png");
//void plotandextlowpt( const Char_t * type, TString suffix){
void plotandexttransversemomentum( TString suffix="PNG"){
//void plotandextlowpt( const Char_t * type="EsBin0", const Char_t * inFN1="inel/MeanPtVsNm_EsBin0.root",const Char_t * outDir="EsBin0", TString suffix="png",const char *outf="MeanPtVsNm_inelf.root"){


const Char_t * inFN1="SOpp13TeV_MERGE_newbinsSOPCSO_LHC15fpass2.root";
const Char_t * outDir="Extrapolation";
const char *outf="SOpp13TeV_MERGE_newbinsSOPCSO_LHC15fpass2_ext.root";

CreateDir(outDir);
CreateDir(Form("%s",outDir));


TFile *f1 = TFile::Open(inFN1, "read");//outputdedx
if(!f1) {
        Printf("FATAL: File 1 doesn't exist");
        return;
}


Double_t
LevyTsallis_Char(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.); 
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) /(n*C); 
  Double_t part5 = TMath::Power(part4, -n); 
  return pt * norm * part3 * part5;
}


TF1 *
LevyTsallisChar(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Char, 0.1, 20.0, 4); //20
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-6, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-6, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}

//pion =139.57 MeVs  =0.13957 GeVs
//TF1* fpions05= LevyTsallis("fpions05",0.139,8.0,0.07,4631.0);  
TF1* fTsallisc= LevyTsallisChar("fTsallisc",0.139,8.0,0.07,4631.0);
//h->Fit(fpions05);

TH1D* hPtmvsNm[nMultbins];
TH1D *hPtmvsNmnew[nMultbins];
TH1D *hNorm[nbinesntrk-initbin+1];

Printf("getting pt spectra");
for(Int_t i_so=0; i_so<1; i_so++){
 for(Int_t i=1; i<nMultbins; i++){
   hNorm[i] = 0;
   hNorm[i] =(TH1D* )f1->Get(Form("hSoM%d",i));
   if(i_so==0){
   hPtmvsNm[i]=(TH1D *)f1->Get(Form("hMcOutMult%d",i));
   }else{
   hPtmvsNm[i]=(TH1D *)f1->Get(Form("hMcOutMult%dES%d",i_ntrk, i_so-1));
   }
           
   Double_t nev = hNorm[i]->GetEntries();
                cout<<"nev="<<nev<<endl;
 
//   hPtmvsNmnew[i]=(TH1D *)hPtmvsNm[i]->Clone(Form("hPtmvsNmnew%s",i));
  // hPtmvsNm->Fit(fTsallisc);

   Double_t totalpt = 0; 
                Double_t totalyield = 0; 
                Double_t e_totalpt = 0; 
                Double_t e_totalyield = 0; 

   for(Int_t bin = binStart; bin <= binStop; ++bin){
                       Double_t pt =  hPtmvsNm[[i]->GetBinCenter(bin);
                       cout<<"  low edge, pt="<<hPtmvsNm[[i]->GetBinLowEdge(bin)<<"   pt="<<pt<<endl;
                       Double_t dpt =  hPtIn[i_ntrk-initbin]->GetBinWidth(bin);
                       Double_t deta = 0.6;
                       Double_t yield = hPtIn[i_ntrk-initbin]->GetBinContent(bin);
                       Double_t e_yield = hPtIn[i_ntrk-initbin]->GetBinError(bin);
   }


 }
}
Printf("choose");




//TF1 *myfit = (TF1*) hPtmvsNm->GetFunction("fTsallisc");
//Double_t chis=myfit->GetChisquare();
//Int_t ndf=myfit->GetNDF();


//Double_t a=myfit->Eval(1.5);
//Double_t b=myfit->Eval(2.5);
/*
Printf("a=%f, b=%f, chis=%f,ndf=%d",a,b, chis,ndf);
Double_t cont=0;
Double_t cont2=0;

for(Int_t i=0; i<=hPtmvsNm->GetNbinsX(); i++){
  cont=hPtmvsNm->GetBinContent(i);
  cont2=hPtvsNm->GetBinContent(i);
  if(i==2){
   hPtmvsNmnew->SetBinContent(i,a);
   hPtvsNmnew->SetBinContent(i,a2);
  }
  else if(i==3){
   hPtmvsNmnew->SetBinContent(i,b);
   hPtvsNmnew->SetBinContent(i,b2);
  }
  else{
    hPtmvsNmnew->SetBinContent(i,cont);
    hPtvsNmnew->SetBinContent(i,cont2);
  }
}
*/

  TCanvas* cmeanpt;
  TLegend* legmeanpt;
  cmanpt= new TCanvas("cmeanpt%","THmeanpt",200,10,900,600);
  hPtmvsNm[4]->SetMarkerStyle(24);
//  hPtmvsNmnew[4]->SetMarkerStyle(7);
  hPtmvsNm[4]->SetMarkerColor(2);
//  hPtmvsNmnew[4]->SetMarkerColor(4);
  hPtmvsNm[4]->Draw();
//  hPtmvsNmnew[4]->Draw("sames");


/*

 TFile *fout = TFile::Open( outf,"RECREATE");
 fout->cd();
 hPtmvsNm->Write();
 myfit->Write();
 hPtmvsNmnew->Write();
 hPtvsNm->Write();
 myfit2->Write();
 hPtvsNmnew->Write();
 fout->Close();
*/
}
