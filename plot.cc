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
  if(argc != 5) {
    cout << "Usage: " << argv[0] << " [index_f] [position] [max] [min]" << endl;
    return 1;
  }

  int index_f = atoi(argv[1]);
  int position = atoi(argv[2]);
  int max = atoi(argv[3]);
  int min = atoi(argv[4]);
  string ind = to_string(index_f);
  string pos = to_string(position);

  string file1 = "./potential_";
  string file2 = "./coupling_";
  string filename1 = file1 + ind + "_" +pos+".txt";
  string filename2 = file2 + ind + "_" +pos+".txt";
  string filename3 = file2 + ind + "_" +pos+".pdf";

  string forcename = "";
  int DataList = 1;
  vector<string> DataName;
  vector<string> LegName;
  DataName.push_back(filename2);
  if(index_f ==2){
    forcename = "2";
    DataName.push_back("./data/Eot-Wash.2.txt"); // Heckel PRL111,151802
    DataName.push_back("./data/UVA.2.txt"); //Ritter, PRD42, 977
    DataName.push_back("./data/Kotler.2.txt");
  }else if(index_f == 3){
    forcename = "3";
    DataName.push_back("./data/Terrano.3.txt");
    DataName.push_back("./data/Eot-Wash.3.txt");
    DataName.push_back("./data/NTHU.3.txt");
    DataName.push_back("./data/Kotler.3.txt");
    DataName.push_back("./data/axion.txt");
  }else if(index_f == 4){
    forcename = "4+5";
    DataName.push_back("./data/stellar.4.txt");
  }else if(index_f == 6){
    forcename = "6+7";
  }else if(index_f == 8){
    forcename = "8";
  }else if(index_f == 9){
    forcename = "9+10";
    DataName.push_back("./data/Hammond.9.txt");
    DataName.push_back("./data/Ni.9.txt");
    DataName.push_back("./data/Eot.9.txt");
    DataName.push_back("./data/Terrano.9.txt");
    DataName.push_back("./data/quax.9.txt");
    DataName.push_back("./data/axion.9.txt");
    DataName.push_back("./data/stellar.9.txt");
  }else if(index_f == 11){
    forcename = "11";
    DataName.push_back("./data/Eot-Wash.11.txt");

  }else if(index_f == 12){
    forcename = "12+13";
  }else if(index_f == 14){
    forcename = "14";
  }else if(index_f == 15){
    forcename = "15";
    DataName.push_back("./data/darkphoton.3.txt");
  }else if(index_f == 16){
    forcename = "16";
  }
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
  TMultiGraph *mg = new TMultiGraph();
  for(int i=0;i<input;i++){
    g[i]->GetXaxis()->SetTitle("lambda(m)");
    g[i]->GetYaxis()->SetTitle("Log(gg)");
    g[i]->GetYaxis()->SetTitleOffset(1.3);
    g[i]->SetLineColor(Color.at(i));
    g[i]->SetMarkerStyle(Style.at(i));
    g[i]->SetMarkerSize(1);
    g[i]->SetLineWidth(3);
    mg->Add(g[i],"l");
  }
  

  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.08);
  gPad->SetTopMargin(0.05);


  mg->SetTitle(Form(";#lambda(m);Log(f_{%s})",forcename.c_str()));
  mg->Draw("a");
  mg->GetXaxis()->SetLimits(0.000001,0.2);
  mg->SetMinimum(min);
  mg->SetMaximum(max);

  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleOffset(1.1);
  mg->GetYaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleOffset(1.1);

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
  TLatex latex;
  if(ind == "2"){
    latex.SetTextSize(0.05);
    latex.DrawLatex(2e-4,-30,"case 4");
    latex.SetTextSize(0.04);
    latex.DrawLatex(2e-6,-16,"Kotler, PRL115, 081801");
    latex.SetTextSize(0.035);
    latex.DrawLatex(3e-3,-28,"Ritter, PRD42, 977");
    latex.DrawLatex(6e-4,-38,"Heckel, PRL111,151802");
  }else if(ind == "3"){
    latex.SetTextSize(0.05);
    latex.DrawLatex(3e-5,-29,"Axion");
    latex.DrawLatex(1.3e-4,-7,"This Work");
    latex.SetTextSize(0.04);
    latex.DrawLatex(2e-6,-5,"Kotler, PRL115, 081801");
    latex.DrawLatex(2e-3,-18,"Terrano, PRL115,201801");
    latex.SetTextSize(0.035);
    latex.DrawLatex(2e-3,-3,"Ni, Physica 194B-196B, 153");
    latex.DrawLatex(4e-3,-9,"Heckel, PRL111,151802");
  }else if(ind == "4"){
    latex.DrawLatex(1e-3,-21,"This Work");
    latex.DrawLatex(3e-4,-31,"stellar cooling");
  }else if(ind == "9"){
    latex.SetTextSize(0.04);
    latex.DrawLatex(2.25e-5,-34.5,"Axion");
    latex.SetTextSize(0.05);
    latex.DrawLatex(1.4e-4,-20,"This Work");
    latex.SetTextSize(0.04);
    latex.DrawLatex(3e-5,-28,"Terrano, PRL115,201801");
    latex.SetTextSize(0.035);
    latex.DrawLatex(3e-2,-27,"Ni,PRL82 2439");
    latex.DrawLatex(1e-2,-30,"Crescinia, PLB773, 677");
    latex.DrawLatex(2e-6,-24,"Hodel, PRL106,041801");
    latex.DrawLatex(3e-3,-22,"Hammond, PRL98,081101");
    latex.DrawLatex(2e-4,-31.2,"stellar cooling");
  }else if(ind == "11"){
    latex.SetTextSize(0.05);
    latex.DrawLatex(2e-4,-18,"This Work");
    latex.DrawLatex(1e-4,-26,"Heckel, PRL111,151802");
  }else if(ind == "15"){
    latex.SetTextSize(0.05);
    latex.DrawLatex(2e-2,-2,"This Work");
    latex.SetTextSize(0.04);
    latex.SetTextAngle(15);
    latex.DrawLatex(1e-4,-19.,"Dark Photon");
  }

  //latex.SetTextSize(0.05);
  //latex.DrawLatex(1.4e-4,-25,"case 1");

  /*

  //  TLegend *leg = new TLegend(0.5,0.8,0.92,0.95);
  //  TLegend *leg = new TLegend(0.14,0.15,0.6,0.23);
  TLegend *leg = new TLegend(0.15,0.15,0.7,0.23);
  leg->AddEntry(g[0],"10 fT LANL SERF, gap = 1 cm; case 1","l");
  //  leg->AddEntry(g[3],"1 fT/#sqrt{Hz}, gap = 0.1 mm; case 1","l");
  //leg->AddEntry(g[1],"Lesile, et al, PRD89,114022","l");
  //leg->AddEntry(g[2],"Chu, et al, PRD91,102006","l");
  //  leg->AddEntry(g[3],"Yan, et al, PRL110,082003(g_{A}^{n}g_{V}^{N})","l");
  leg->SetTextFont(132);
  leg->SetBorderSize(0);
  leg->Draw();
  */
  c1->Print(filename3.c_str());
  return 1;
}


