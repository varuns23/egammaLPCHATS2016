#include <iostream>
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGlobalFunc.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include  "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include  "RooFitResult.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TPave.h"
#include "TStyle.h"
#include "TROOT.h"
#include "RooGenericPdf.h"

using namespace RooFit ;

void fit(){
  
  gStyle->Print();
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasDefH(500); //Height of canvas                                                                                                                       
  gStyle->SetCanvasDefW(550); //Width of canvas                                                                                                                        
  gStyle->SetErrorX(0.);
  gStyle->SetMarkerStyle(20);

  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

  // Margins:                                                                                                                                                          
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.17); //was 0.16                                                                                                                           
  gStyle->SetPadRightMargin(0.05);//was 0.02                                                                                                                           

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.40); //was 1.25                                                                                                                            
  // For the axis labels:                                                                                                                                              
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

  RooRealVar x("x","M_{ep} (in GeV)",60,120) ;
  TFile *inputfile = new TFile("ELE_FakeRate.root","READ");
  inputfile->cd();

  TH1F* heg1 = dynamic_cast<TH1F*>(inputfile->Get("E-E Invariant Mass5060"));
  RooDataHist dh1("dh1","dh1",x,Import(*heg1)) ;
  RooPlot* frame = x.frame(Title("M_{ep}")) ;
  dh1.plotOn(frame,Name("dh1")) ;
  RooRealVar alpha1("alpha1","alpha1",0.001,-0.01,0.01);
  RooRealVar beta1("beta1","beta1",0.01,-0.1,0.1);//200
  RooRealVar mG("a1","a1",88.93,88,92);
  RooRealVar sG("a1","a1",4.038,4,5);
  RooGaussian gaus("g","g",x,mG,sG);
  RooRealVar ac("a","a",1,-10,10);
  RooRealVar bc("b","b",1,-10,10);
  RooRealVar pc("p","p",91.18);
  RooRealVar gc("g","g",1,-10,10);
  RooExponential shape1("shape1","shape1",x,beta1);
  RooRealVar cs("cs","cs",50); //RooRealVar cs("cs","cs",60);  
  RooRealVar as("as","as",.1); //RooRealVar as("as","as",.01);
  RooRealVar bs("bs","bs",-0.04 ); // RooRealVar bs("bs","bs",.01 );
  RooRealVar nsig1("nsig1","signal events1",10,0.,100000.);
  RooRealVar nbkg1("nbkg1","signal background events1",10,0.,100000.);
  
  RooRealVar mBW("mBW","mean BW",91.1875) ;
  RooRealVar sBW("sBW","sigma BW",2.4952) ;
  RooBreitWigner brietwigner("BW","BW",x,mBW,sBW) ;

  RooRealVar mCB1("mCB1", "mCB" ,0.0,-0.01,0.01) ;
  RooRealVar sCB1("sCB1", "sCB1" ,1.64,1,2.5 )  ;
  RooRealVar nCB1("nCB1","", 100);
  RooRealVar alphaCB1("alphaCB1","", 0.6,.2,1); // RooRealVar alphaCB1("alphaCB1","", 0.6,.2,1);
  RooCBShape cball1("cball1","cball1",x,mCB1,sCB1, alphaCB1, nCB1);


  RooGenericPdf sigmoid("sigmoid","sigmoid","(cs/(1+exp(-as*(x-65)))*exp(bs*x))",RooArgSet(x,as,bs,cs)) ;
  RooFFTConvPdf bwcball1("BWxCB","Breit-Wigner (X) CB",x,brietwigner, cball1) ;
  RooAddPdf  model1("model1","BW X CB and SIG",RooArgList(bwcball1,sigmoid),RooArgList(nsig1,nbkg1)) ;
  
 
  //  frame->Draw() ;
  // c1->SaveAs("photonelectron1(win80100).png")  ;*/
  
  
  model1.fitTo(dh1,Extended(),Range(60,120));                                                                                
  model1.plotOn(frame,Name("model1"),LineColor(kRed));
  model1.plotOn(frame,Components(sigmoid),LineStyle(kDashed)) ;
  model1.plotOn(frame,Components(bwcball1),LineColor(kGreen)) ;
  //model2.plotOn(frame,Components(cball2),LineColor(kGreen)) ;
  //  model.plotOn(frame,Components(gaus),LineColor(kGreen)) ;
  model1.paramOn(frame,Layout(0.62,0.90),Format("NEU",AutoPrecision(1))) ;
  frame->getAttText()->SetTextSize(0.026) ;
  TCanvas* c2 = new TCanvas("DATA","DATA",800,800) ;
  gPad->SetLeftMargin(0.19) ; frame->GetYaxis()->SetTitleOffset(1.8) ;
  gPad->SetRightMargin(0.10) ; frame->GetXaxis()->SetTitleOffset(1.2) ;
  
  Double_t chi2 = frame->chiSquare("model1", "dh1", 9);
  std::cout<<"Chi Square=:"<<chi2<<std::endl;

  //  RooCurve* func = (RooCurve*) frame->findObject("sigmoid") ;
  x.setRange("signal",80,100) ;
  RooAbsReal* fracSig = bwcball1.createIntegral(x, x, "signal") ;
  Double_t nsig_frac  = nsig1.getVal() * fracSig->getVal();
  RooAbsReal* fracBkg = sigmoid.createIntegral(x, x, "background") ;
  Double_t nbkg_frac  = nbkg1.getVal() * fracBkg->getVal();

  cout<< " Signal = " << nsig_frac << endl;
  cout<< " BKGD = " << nbkg_frac << endl;
  
  frame->Draw();
  
}
