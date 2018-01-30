#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TLine.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TLatex.h"
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
//#include <unistd.h>
#include <getopt.h>
using namespace std;


int ctob(int c = 3){
  //cout << c <<endl;
  int remain = 0;    
  int b = 0;         
  int count = 0;     
  while(c>0){        
    remain = c%2;    
    b = b+remain*pow(10,count);
    c=c/2;                     
    count++;                   
  }                            
  return b;                    
}    



Double_t fitf(Double_t *v, Double_t *par){
  Double_t arg = 0;                       
  if (par[2] != 0) arg = (v[0] - par[1])/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg) + par[3];
  return fitval;                                             
}                                                            

Double_t linear(Double_t *v, Double_t *par){

  Double_t fitval = -v[0]+par[0] ;
  return fitval;
}
Double_t inverse(Double_t *v, Double_t *par){

  Double_t fitval = par[0]/v[0] ;
  return fitval;
}



Double_t horizontal(Double_t *v, Double_t *par){

  Double_t fitval=par[0];
  return fitval;
}



TGraph* GetGraph(vector<Double_t> Px, vector<Double_t> Py){
  const Int_t nPy = Py.size();

  Double_t A[nPy];
  Double_t B[nPy];

  for(Int_t j=0;j<nPy;j++){
    A[j]=Px.at(j);
    B[j]=Py.at(j);
  }
  TGraph* fGraph = new TGraph(nPy,A,B);
  fGraph->SetFillStyle(0);
  fGraph->SetFillColor(0);
  return fGraph;
}


int main(int argc, char** argv)
{
  if(argc != 2) {
    cout << "Usage: " << argv[0] << " [any number]" << endl;
    return 1;
  }

  string filename1 = "LC.pdf";

  string forcename = "";
  int DataList = 1;
  vector<string> DataName;
  vector<string> LegName;

  DataName.push_back("./data/CAST.txt");
  DataName.push_back("./data/sikivie.txt"); // 
  DataName.push_back("./data/proposal.txt"); //
  DataName.push_back("./data/qcdaxion.txt"); //
  DataName.push_back("./data/prelim.txt"); //

  DataList = DataName.size();
  vector<Int_t> Color;
  Color.push_back(2);
  Color.push_back(3);
  Color.push_back(4);
  Color.push_back(6);
  Color.push_back(7);
  Color.push_back(8);
  Color.push_back(9);
  Color.push_back(12);
  Color.push_back(28);
  Color.push_back(32);
  Color.push_back(37);
  Color.push_back(38);
  Color.push_back(41);
  Color.push_back(44);
  Color.push_back(46);
  Color.push_back(49);


  vector<Int_t> Style;
  Style.push_back(20);
  Style.push_back(24);
  Style.push_back(21);
  Style.push_back(25);
  Style.push_back(22);
  Style.push_back(26);
  Style.push_back(23);
  Style.push_back(32);
  Style.push_back(33);
  Style.push_back(27);
  Style.push_back(34);
  Style.push_back(28);
  Style.push_back(29);
  Style.push_back(30);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetFrameBorderMode(0);
  TCanvas *c1 = new TCanvas("c1", "Human-Readable Name");
  c1->SetFillColor(0);  

  const int input =DataList;

  ifstream fin[input];
  double r1, r2;
  vector<double> Lambda[input];
  vector<double> Coupling[input];
  TGraph *g[input];

  for(int i = 0;i<input;i++){
    cout << DataName.at(i) << endl;
    fin[i].open(DataName.at(i));
    if(fin[i].is_open()){
      while(!fin[i].eof()){
	fin[i] >> r1 >> r2;
	Lambda[i].push_back(r1);
	Coupling[i].push_back(r2);
	cout << r1 << " " << r2 << endl;
      }
    }
    Lambda[i].pop_back();
    Coupling[i].pop_back();
    g[i] = GetGraph(Lambda[i],Coupling[i]);
  }

  c1->SetLogx();
  c1->SetLogy();
  TMultiGraph *mg = new TMultiGraph();
  for(int i=0;i<input-1;i++){
    g[i]->GetXaxis()->SetTitle("lambda(m)");
    g[i]->GetYaxis()->SetTitle("Log(gg)");
    g[i]->GetYaxis()->SetTitleOffset(1.3);
    g[i]->SetLineColor(Color.at(i));
    g[i]->SetMarkerStyle(Style.at(i));
    g[i]->SetMarkerSize(1);
    g[i]->SetLineWidth(3);
    mg->Add(g[i],"l");
  }
  int i = input-1;
  g[i]->GetXaxis()->SetTitle("lambda(m)");
  g[i]->GetYaxis()->SetTitle("Log(gg)");
  g[i]->GetYaxis()->SetTitleOffset(1.3);
  g[i]->SetLineColor(Color.at(i));
  g[i]->SetMarkerStyle(22);
  g[i]->SetMarkerSize(1);
  g[i]->SetLineWidth(3);
  mg->Add(g[i],"p");

  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTopMargin(0.05);


  mg->SetTitle(";m_{a}(eV);g_{a#gamma#gamma}(GeV^{-1})");
  mg->Draw("a");
  mg->GetXaxis()->SetLimits(1e-11,1e-6);
  mg->SetMinimum(1e-19);
  mg->SetMaximum(1e-2);

  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleOffset(1.1);
  mg->GetYaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleOffset(1.1);
  /*
  TF1 *f1 = new TF1("f1",inverse,1e-6,2e-1,1);
  f1->SetParameter(0,2e-7);
  TGaxis *ax1 = new TGaxis(2e-1,max,1e-6,max,"f1",510,"G+");
  ax1->SetTitle("m(eV)");
  ax1->SetLabelSize(0.);
  ax1->SetLabelOffset(-0.04);
  ax1->SetLabelFont(40);
  ax1->SetTitleFont(42);
  ax1->SetTitleSize(0.);
  ax1->SetTitleOffset(1.2);
  ax1->Draw();
  */

  TLatex latex;

  latex.SetTextSize(0.05);
  latex.DrawLatex(1e-9,8e-10,"CAST");
  //latex.SetTextSize(0.04);
  latex.DrawLatex(1e-9,1e-17,"Sikivie");
  //latex.SetTextSize(0.035);
  latex.DrawLatex(1e-10,1e-5,"Preliminary Data");
  latex.DrawLatex(1e-9,1e-14,"Proposal");
  latex.DrawLatex(9e-8,1e-16,"QCD Axion");


  c1->Print(filename1.c_str());
  return 1;
}


