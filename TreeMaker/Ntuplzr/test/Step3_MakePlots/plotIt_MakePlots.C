#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TChain.h"
#include "TKey.h"
#include "TTree.h"
#include "Riostream.h"

#pragma once

#define TITLE_FONTSIZE 26
#define LABEL_FONTSIZE 18

#define LEFT_MARGIN 0.17
#define RIGHT_MARGIN 0.03
#define TOP_MARGIN 0.05
#define BOTTOM_MARGIN 0.13

//////////////////////
//TDRStyle////////////
//////////////////////

TStyle* createMyStyle(bool logy=0) {
  TStyle *myStyle = new TStyle("myStyle", "myStyle");

  TGaxis::SetExponentOffset(-0.07, -0.01, "y");

  // For the canvas:
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetCanvasColor(kWhite);
  myStyle->SetCanvasDefH(800); //Height of canvas
  myStyle->SetCanvasDefW(800); //Width of canvas
  myStyle->SetCanvasDefX(0);   //POsition on screen
  myStyle->SetCanvasDefY(0);

  // For the Pad:
  myStyle->SetPadBorderMode(0);
  myStyle->SetPadColor(kWhite);
  myStyle->SetPadGridX(false);
  myStyle->SetPadGridY(false);
  myStyle->SetGridColor(0);
  myStyle->SetGridStyle(3);
  myStyle->SetGridWidth(1);

  // For the frame:
  myStyle->SetFrameBorderMode(0);
  myStyle->SetFrameBorderSize(1);
  myStyle->SetFrameFillColor(0);
  myStyle->SetFrameFillStyle(0);
  myStyle->SetFrameLineColor(1);
  myStyle->SetFrameLineStyle(1);
  myStyle->SetFrameLineWidth(1);

  // For the histo:
  myStyle->SetHistLineStyle(1);
  myStyle->SetHistLineWidth(2);
  myStyle->SetEndErrorSize(2);

  //For the fit/function:
  myStyle->SetFitFormat("5.4g");
  myStyle->SetFuncColor(2);
  myStyle->SetFuncStyle(1);
  myStyle->SetFuncWidth(1);

  // For the statistics box:
  myStyle->SetOptFile(0);
  myStyle->SetStatColor(kWhite);
  //myStyle->SetStatFont(43);
  //myStyle->SetStatFontSize(0.025);
  myStyle->SetStatTextColor(1);
  myStyle->SetStatFormat("6.4g");
  myStyle->SetStatBorderSize(1);
  myStyle->SetStatH(0.12);
  myStyle->SetStatW(0.3);
  myStyle->SetStatY(0.92);
  myStyle->SetStatX(0.94);

  //For the date:
  myStyle->SetOptDate(0);

  // Margins:
  myStyle->SetPadTopMargin(TOP_MARGIN);
  myStyle->SetPadBottomMargin(BOTTOM_MARGIN);
  myStyle->SetPadLeftMargin(LEFT_MARGIN);
  myStyle->SetPadRightMargin(RIGHT_MARGIN);

  // For the Global title:
  myStyle->SetOptTitle(0);
  myStyle->SetTitleFont(63);
  myStyle->SetTitleColor(1);
  myStyle->SetTitleTextColor(1);
  myStyle->SetTitleFillColor(10);
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleAlign(33); 
  myStyle->SetTitleX(1);
  myStyle->SetTitleFontSize(TITLE_FONTSIZE);

  // For the axis titles:

  myStyle->SetTitleColor(1, "XYZ");
  myStyle->SetTitleFont(43, "XYZ");
  myStyle->SetTitleSize(TITLE_FONTSIZE, "XYZ");
  myStyle->SetTitleYOffset(1.75); 
  myStyle->SetTitleXOffset(1.5);

  myStyle->SetLabelColor(1, "XYZ");
  myStyle->SetLabelFont(43, "XYZ");
  myStyle->SetLabelOffset(0.01, "YZ");
  myStyle->SetLabelOffset(0.015, "X");
  myStyle->SetLabelSize(LABEL_FONTSIZE, "XYZ");

  myStyle->SetAxisColor(1, "XYZ");
  myStyle->SetStripDecimals(kTRUE);
  myStyle->SetTickLength(0.03, "XYZ");
  myStyle->SetNdivisions(510, "XYZ");
  myStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  myStyle->SetPadTickY(1);

  myStyle->SetOptLogx(0);
  myStyle->SetOptLogy(logy);
  myStyle->SetOptLogz(0);

  myStyle->SetHatchesSpacing(1.3);
  myStyle->SetHatchesLineWidth(1);

  myStyle->cd();

  return myStyle;
}


//////////////////////////
//////legend style////////
//////////////////////////

void leg_myStyle(TLegend *leg,
    //int text_size = 0.035,
    //int text_font = 43,
    //int text_align = 22,
    int text_align = 12,
    int fill_style = 1,
    int fill_color = 10,
    int line_color = 0,
    int line_width = 0,
    int border_size = 1) {

  //leg->SetTextSize(text_size);
  //leg->SetTextFont(text_font);
  leg->SetTextAlign(text_align);
  leg->SetFillStyle(fill_style);
  leg->SetFillColor(fill_color);
  leg->SetLineColor(line_color);
  leg->SetLineWidth(line_width);
  leg->SetBorderSize(border_size);
}

//////////////////////////
//////TH1F style////////
//////////////////////////

void h_myStyle(TH1 *h,
	       int line_color=1,
	       int line_style=1, 
	       int fill_color=50,
	       int fill_style=1001,
	       float y_min=-1111.,
	       float y_max=-1111.,
	       int ndivx=510,
	       int ndivy=510,
	       int marker_style=20,
	       int marker_color=1,
	       float marker_size=1.2,
	       int optstat=0){
  //    TString xtitle="") {

  h->SetLineWidth(3);
  h->SetLineColor(line_color);
  h->SetLineStyle(line_style);
  h->SetFillColor(fill_color);
  h->SetFillStyle(fill_style);
  h->SetMaximum(y_max);
  h->SetMinimum(y_min);
  h->GetXaxis()->SetNdivisions(ndivx);
  h->GetYaxis()->SetNdivisions(ndivy);
  h->GetYaxis()->SetTitleOffset(2.);

  h->SetMarkerStyle(marker_style);
  h->SetMarkerColor(marker_color);
  h->SetMarkerSize(marker_size);
  h->SetStats(optstat);
  //h->SetBinErrorOption(TH1::kPoisson2);
  
  //if (xtitle.Length() > 0)
  //  h->GetXaxis()->SetTitle(xtitle);
}


//////////////////////////
//////graph style/////////
//////////////////////////

void g_myStyle(TGraphAsymmErrors *h,
	       int line_color=1,
	       int line_style=1, 
	       int fill_color=50,
	       int fill_style=1001,
	       float y_min=-1111.,
	       float y_max=-1111.,
	       int ndivx=510,
	       int ndivy=510,
	       int marker_style=20,
	       int marker_color=1,
	       float marker_size=1.2){
  //    TString xtitle="") {

  h->SetLineWidth(3);
  h->SetLineColor(line_color);
  h->SetLineStyle(line_style);
  h->SetFillColor(fill_color);
  h->SetFillStyle(fill_style);
  h->SetMaximum(y_max);
  h->SetMinimum(y_min);
  h->GetXaxis()->SetNdivisions(ndivx);
  h->GetYaxis()->SetNdivisions(ndivy);
  h->GetYaxis()->SetTitleOffset(2.);

  h->SetMarkerStyle(marker_style);
  h->SetMarkerColor(marker_color);
  h->SetMarkerSize(marker_size);
  //h->SetStats(optstat);
  //h->SetBinErrorOption(TH1::kPoisson2);
  
  //if (xtitle.Length() > 0)
  //  h->GetXaxis()->SetTitle(xtitle);
}


/////////////////////////////
//////String txt/style///////
/////////////////////////////

void cms_myStyle(TString pileup){
  //std::string status = "Phase-2 Simulation Preliminary";
  std::string status = "Phase-2 Simulation";
  TPaveText* pt_exp = new TPaveText(LEFT_MARGIN, 1 - 0.5 * TOP_MARGIN, 1 - RIGHT_MARGIN, 1, "brNDC");
  pt_exp->SetFillStyle(0);
  pt_exp->SetBorderSize(0);
  pt_exp->SetMargin(0);
  pt_exp->SetTextFont(62);
  pt_exp->SetTextSize(0.75 * TOP_MARGIN);
  pt_exp->SetTextAlign(13);
  TString d = TString::Format("CMS #font[52]{#scale[0.76]{%s}}", status.c_str());
  pt_exp->AddText(d);
  pt_exp->Draw();

  //TString lumi_s = "3 ab^{-1} (14 TeV";
  TString lumi_s = "(14 TeV";
  if (pileup.Length() > 0) lumi_s = lumi_s + TString::Format(", %s PU)", pileup.Data());
  else lumi_s = lumi_s + ")";
  TPaveText* pt_lumi = new TPaveText(LEFT_MARGIN, 1 - 0.5 * TOP_MARGIN, 1 - RIGHT_MARGIN, 1, "brNDC");
  pt_lumi->SetFillStyle(0);
  pt_lumi->SetBorderSize(0);
  pt_lumi->SetMargin(0);
  pt_lumi->SetTextFont(42);
  pt_lumi->SetTextSize(0.6 * TOP_MARGIN);
  pt_lumi->SetTextAlign(33);
  pt_lumi->AddText(lumi_s);
  pt_lumi->Draw();
}


/////////////////////////////////////////////////////////////////
//////Efficiency plotting for a given num and denominator////////
/////////////////////////////////////////////////////////////////


//---------------------------------------------------------------
int Eff(TString outdir, TH1F* h1, TH1F* h2, TString title="", TString pu="", bool sig=1)
//---------------------------------------------------------------
{
  TString sName = TString::Format("%s",h1->GetName());
  std::string sName2 = (std::string) TString::Format("%s",h2->GetName());

  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  g->Divide(h1,h2);  
  
  g->GetHistogram()->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  g->GetHistogram()->GetYaxis()->SetTitle(title+" Efficiency");
  float y_min;
  if(sig) y_min = 0.1;
  else y_min = 0.00005;
  g_myStyle(g,1,1,1,2,y_min,1.1,510,510,22,1,1.2);
  TCanvas* cn = new TCanvas("cn","canvas");
  cn->cd();
  g->Draw("AP");
  g->SetName("Eff"+sName);

  if(sName2.find("pt")!=-1)
    {
      g->GetHistogram()->GetXaxis()->SetRangeUser(0,350);
    }
  cms_myStyle(pu);
  if(sName2.find("Gen")!=-1)
    {
      cn->SaveAs(outdir+"Eff"+sName+"Gen.C");
      cn->SaveAs(outdir+"Eff"+sName+"Gen.pdf");
      cn->SaveAs(outdir+"Eff"+sName+"Gen.png");
    }
  else 
    {
      cn->SaveAs(outdir+"Eff"+sName+".C");
      cn->SaveAs(outdir+"Eff"+sName+".pdf");
      cn->SaveAs(outdir+"Eff"+sName+".png");
    }
  delete cn;
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////
//////gets efficiency plots for all IDs/particles on separate canvases///////////////
/////////////////////////////////////////////////////////////////////////////////////

int plotEffSeparate(TFile *inFile,TString outdir = "HistogramsTogether/", TString pu="", bool sig=1)
{
 TString nums1[4]={"Loose","Medium","Tight","Good"};
 TString den1[1]={"All"};
 TString match[2]={"Matched",""}; //signal = 0, bkg =1;
 int m_id;
 if(sig) m_id = 0;
 else m_id = 1;
 TString histos1[4]={"pt","eta", "phi","nvtx"};
 TString part1[2]={"Elec","Muon"};
 for (int nh=0; nh< 4 ; nh++)
   for (int npart=0; npart< 1 ; npart++)
     for (int id=0; id< 4 ; id++)
       {
	 TString NumHisto = histos1[nh]+nums1[id]+match[m_id]+part1[npart];
	 TString DenHisto = histos1[nh]+den1[0]+match[m_id]+part1[npart];
	 std::cout<<NumHisto<< "  " << DenHisto << std::endl;
	 TH1F *h1 = (TH1F*) inFile->Get(NumHisto); 
	 TH1F *h2 = (TH1F*) inFile->Get(DenHisto); 
	 std::cout<<h1->GetName() << std::endl;
	 std::cout<<h2->GetName() << std::endl;
	 Eff(outdir, h1, h2,nums1[id]+part1[npart],pu,sig);
	 if(nh==0) std::cout<<nums1[id]<<" Efficiency: " <<h1->GetEntries()/h2->GetEntries() <<std::endl;
       }

 //overall efficiency wrt pt only for Loose/Tight/Iso ID in EB and EE sepearately after genmatching
 TString eta[2]={"_EB","_EE"};
 for (int nh=0; nh< 1 ; nh++)
   for (int npart=0; npart< 1 ; npart++)
     for (int id=0; id< 4 ; id++)
       for (int ieta=0; ieta< 2 ; ieta++)
       {
	 TString NumHisto = histos1[nh]+nums1[id]+match[m_id]+part1[npart]+eta[ieta];
	 TString DenHisto = histos1[nh]+den1[0]+match[m_id]+part1[npart]+eta[ieta];
	 std::cout<<NumHisto<< "  " << DenHisto << std::endl;
	 TH1F *h1 = (TH1F*) inFile->Get(NumHisto); 
	 TH1F *h2 = (TH1F*) inFile->Get(DenHisto); 
	 std::cout<<h1->GetName() << std::endl;
	 std::cout<<h2->GetName() << std::endl;
	 Eff(outdir, h1, h2,nums1[id]+part1[npart], pu,sig);
       }
 return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
///getting graphs with efficiencies for a given num and denominator with ////////////////
///////////////////their names to be plotted on the same canvas//////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------
int Eff2(TH1F* h1, TH1F* h2, TGraphAsymmErrors *g, int id_color, float y_min)
//---------------------------------------------------------------
{
  TString sName = TString::Format("%s",h1->GetName());
  //std::string sName2 = (std::string) TString::Format("%s",h2->GetName());
  //g->Divide(h1,h2);  
  g->GetHistogram()->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  g->GetHistogram()->GetYaxis()->SetTitle("ID Efficiency");
  g_myStyle(g,id_color,1,1,2,y_min,1.1,510,510,22,id_color,1.2);
  g->SetName("Eff"+sName);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////
//////gets efficiency plots for all IDs on the same canvas///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

int plotEffAll(TFile *inFile, TString outdir = "HistogramsTogether/",TString pu="", bool sig=1)
//---------------------------------------------------------------
{  
  const int nhistos1= 4;
  const int nids    = 4;
  const int nparts  = 2;
  
  TString nums1[nids]={"Loose","Medium","Tight","Good"};
  TString den1[1]={"All"};
  TString match[2]={"Matched",""}; //signal = 0, bkg =1;
  TString histos1[nhistos1]={"pt","eta", "phi","nvtx"};
  TString part1[nparts]={"Elec","Muon"};
  
  //my_style = createMyStyle(0);
  //my_style->cd();
  
  //TGraphAsymmErrors *g[nhistos1][nids] = new TGraphAsymmErrors();
  TGraphAsymmErrors *g[nhistos1][nids];
  char name[100];


  for (int npart=0; npart< 2 ; npart++)
  //for (int npart=0; npart< 1 ; npart++)
    for (int nh=0; nh< nhistos1 ; nh++)
    //for (int nh=0; nh< 1 ; nh++)
      {
	TCanvas* cn = new TCanvas("cn","canvas");
	TLegend* leg = new TLegend(0.40,0.15, 0.70, 0.30, NULL,"brNDC");
	leg_myStyle(leg);
	for (int id=0; id< nids ; id++)
	//for (int id=0; id< 1 ; id++)
	  {
	    if(npart==1 && id==1 ) continue;
	    //sprintf(name,"ptJet%", i);
	    TString NumHisto, DenHisto;
	    int m_id;
	    if(sig) m_id = 0;
	    else m_id = 1;
	    NumHisto = histos1[nh]+nums1[id]+match[m_id]+part1[npart]; //ptLooseElec for bkg; ptLooseMatchedElec for sig
	    DenHisto = histos1[nh]+den1[0]+match[m_id]+part1[npart];//ptAllElec for bkg; ptAllMatchedElec for sig
	    //std::cout<<NumHisto<< "  " << DenHisto << std::endl;
	    TH1F *h1 = (TH1F*) inFile->Get(NumHisto); 
	    TH1F *h2 = (TH1F*) inFile->Get(DenHisto); 
	    std::cout<<"Num: "<<h1->GetName() << "    Nbins: "<< h1->GetNbinsX()<<std::endl;
	    std::cout<<"Den: "<<h2->GetName() << "    Nbins: "<< h2->GetNbinsX()<<std::endl;
	    float y_min;
	    if(sig) y_min = 0.1;
	    else y_min = 0.00005;
	    g[nh][id] = new TGraphAsymmErrors(h1,h2);
	    Eff2(h1, h2,g[nh][id], id+1, y_min);
	    cn->cd();
	    if(id==0)
	      g[nh][id]->Draw("AP");
	    else 
	      g[nh][id]->Draw("P");

	    if(nh==0)
	      g[nh][id]->GetHistogram()->GetXaxis()->SetRangeUser(0,350);
	    

	
	    TString sName = TString::Format("%s",h1->GetName());
	    leg->AddEntry("Eff"+sName, nums1[id]+" ID","pl");
	    
	    //if(nh==0) std::cout<<nums1[id]<<" Efficiency: " <<h1->GetEntries()/h2->GetEntries() <<std::endl;
	}
	cms_myStyle(pu);
	leg->Draw();
	cn->SaveAs(outdir+"Eff"+histos1[nh]+part1[npart]+".C");
	cn->SaveAs(outdir+"Eff"+histos1[nh]+part1[npart]+".pdf");
	cn->SaveAs(outdir+"Eff"+histos1[nh]+part1[npart]+".png");
	delete cn;
      }
  return 0;
}

////////////////////////////////////////////////////////
/////plots the TH1F histogram///////////////////////////
////////////////////////////////////////////////////////
//---------------------------------------------------------------
int plot1D(TString outdir, TH1F* histo, bool optstat, TString pu)
//---------------------------------------------------------------
{

  TString sName = TString::Format("%s",histo->GetName());
  //if (sName.BeginsWith("Gen"))
  //  h_myStyle(histo,46,1,46,2,-1111,-1111,510,510,22,46,1.2,optstat);
  //else
    h_myStyle(histo,38,1,38,2,-1111,-1111,510,510,22,38,1.2,optstat);

  TCanvas* cn = new TCanvas("cn","canvas");
  cn->cd();
  if (sName.EndsWith("Phi")) histo->SetMinimum(0);
  //histo->DrawCopy("hist");
  histo->Draw("hist");
  cms_myStyle(pu);

  cn->SaveAs(outdir+sName+".C");
  cn->SaveAs(outdir+sName+".pdf");
  cn->SaveAs(outdir+sName+".png");

  delete cn;

  return 0;
}

////////////////////////////////////////////////////////
///// for plotting the Cumulative plot of isolation////
////////////////////////////////////////////////////////

//---------------------------------------------------------------
int plot_isoEff(TString outdir, TH1F* histo, bool optstat, TString pu="")
//---------------------------------------------------------------
{

  TString sName = TString::Format("%s",histo->GetName());
  h_myStyle(histo,1,1,1,0,-1111,1.2,510,510,22,1,1.2,optstat);

  TCanvas* cn = new TCanvas("cn","canvas");
  cn->cd();
  if (sName.EndsWith("Phi")) histo->SetMinimum(0);
  histo->Draw("hist,c");
  cms_myStyle(pu);

  cn->SaveAs(outdir+"Eff"+sName+".C");
  cn->SaveAs(outdir+"Eff"+sName+".pdf");
  cn->SaveAs(outdir+"Eff"+sName+".png");

  delete cn;

  return 0;
}

/////////////////////////////////////////////////////////
///// gets and plots the Cumulative plot of isolation////
/////////////////////////////////////////////////////////

int iso_eff(TH1F *h1, TString outdir = "HistogramsTogether/" , bool optstat= 0, TString pu="")
{
  double sum;
  TH1F *h2=(TH1F*)h1->Clone();
  const int nbins =h2->GetNbinsX();
  for(int i=1; i<=nbins; i++)
    {
      sum=h1->Integral(0,i,"");
      sum=(sum/h1->GetEntries());
      //cout<<sum<<endl;
      h2->SetBinContent(i,sum);
    }
  plot_isoEff(outdir, h2 , optstat, pu);
  return 0;
}



//---------------------------------------------------------------
int plotIt_Eff(TString fiName = "histos.root", TString outdir = "HistogramsTogether/" , TString pu = "",bool sig =0, bool optstat= 0)
//---------------------------------------------------------------
{  
  TFile *inFile = TFile::Open(fiName);
  std::cout<<"opened the root file successfully"<<std::endl;
  gROOT->ProcessLine(".! mkdir -p "+outdir);

  TString outdir_eff;
  if(sig)
    outdir_eff = outdir+"/Eff/";
  else 
    outdir_eff = outdir+"/Eff_bkg/";

  gROOT->ProcessLine(".! mkdir -p "+outdir_eff);
  gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+outdir);
  gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+outdir_eff);
  
  TStyle* my_style = createMyStyle(1);
  my_style->cd();
  // gROOT->SetBatch(true); 
  
  TIter nextkey(inFile->GetListOfKeys());
  
  TKey *key, *oldkey = 0;
  while ((key = (TKey*)nextkey())) {
    //keep only the highest cycle number for each key
    //if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;
    
    TObject *obj = key->ReadObj();
    if (obj->IsA()->InheritsFrom(TH1F::Class()) || obj->IsA()->InheritsFrom(TH1I::Class())) {
      TH1F* histo = (TH1F*)obj;
      if((TString::Format("%s",histo->GetName()).BeginsWith("rel_")))
      	{
      	  iso_eff(histo, outdir_eff, optstat, pu);
      	}
      plot1D(outdir_eff, histo, optstat,pu);
    }
    else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      cout << "Found subdirectory " << obj->GetName() << endl;
      //gSystem->MakeDirectory(outdir_eff+"/"+obj->GetName());
      //plot1D(outdir_eff, histo, optstat);
    } // end of IF a TDriectory 
  }

  //overall efficiency for Loose/Tight/Iso ID after genmatching
  if(sig)
    my_style = createMyStyle(0); // logy? 
  else
    my_style = createMyStyle(1); // logy?
  my_style->cd();
  //plotEffSeparate(inFile,outdir_eff,sig);
  plotEffAll(inFile, outdir_eff,pu, sig);
  
  inFile->Close();
  delete inFile;
  
  return 0;
}



//---------------------------------------------------------------
int plotIt_ROC(TString fiName1 = "histos.root", TString fiName2 = "histos.root",TString outdir = "HistogramsTogether/", TString pu="")
//---------------------------------------------------------------
{ 
  gROOT->ProcessLine(".! mkdir -p "+outdir);
  gSystem->Exec("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+outdir);

  TStyle* my_style = createMyStyle(0);
  my_style->cd();


  TFile *inFile1 = TFile::Open(fiName1);
  TFile *inFile2 = TFile::Open(fiName2);
  std::cout<<"opened the root files successfully"<<std::endl;

  TString label;
  if(fiName1.BeginsWith("../Step2_PostAN/DYToLL"))
    label = "DY";
  TString ID[3] = {"Loose","Tight","Medium"};
  TString iEB[2] = {"EB","EE"};
  TString match[2]= {"Matched",""};
  TString part1[2]={"Elec","Muon"};
  for (int id=0; id< 3 ; id++)
    {  
      for (int npart = 0;npart<2;npart++)
	{ 
	  if(npart==1)
	    { 
	      if  (id==2) continue;
	      TH1F *h1,*h2;
	      TH2D *AccRejPlot = new TH2D("AccRejPlot","", 120, 0, 1.2, 120, 0 ,1.2);
	      TH2D *AccRejPlot2 = new TH2D("AccRejPlot2","", 120, 0, 1.2, 120, 0 ,1.2);
	      TString h1Histo = "rel_iso"+ID[id]+match[0]+part1[npart]; // Signal
	      TString h2Histo = "rel_iso"+ID[id]+match[1]+part1[npart]; // QCD
	      std::cout<<h1Histo<<"," <<h2Histo<< std::endl;
	      h1 = (TH1F*) inFile1->Get(h1Histo);
	      h2 = (TH1F*) inFile2->Get(h2Histo);
	      std::cout<<"got the histograms: "<<h1->GetName()<< ", " << h2->GetName()<<std::endl;
	      TGraph *g;
	      int i = 0;
	      double gx[10], gy[10];
	      for ( float j = 0.0 ; j < 1.0 ; j += 0.1)
		{
		  int bin = h1->FindBin(j);
		  double y = 1 - h2->Integral(0,bin-1)/(double)h2->GetEntries();
		  double x = h1->Integral(0,bin)/(double)h1->GetEntries();
		  
		  std::cout<<i<<std::endl;
		  gx[i]=x;
		  gy[i]=y;
		  
		  AccRejPlot->Fill(x,y,j);
		  AccRejPlot2->Fill(x,y);
		  i++;
		}
	      
	      g = new TGraph(10, gx, gy);
	      std::cout<<i<<std::endl;
	      AccRejPlot->Divide(AccRejPlot2);
	      AccRejPlot->GetZaxis()->SetTitle("Isolation Value");
	      AccRejPlot->GetYaxis()->SetTitle("ID Rejection");
	      AccRejPlot->GetXaxis()->SetTitle("ID Efficiency");
	      AccRejPlot->GetZaxis()->SetRangeUser(0.0, 1.0);
	      AccRejPlot->SetStats(0);
	      //AccRejPlot->GetYaxis()->SetRangeUser(0.1, 1.02);
	      //AccRejPlot->GetXaxis()->SetRangeUser(0.30, 1.01);
	      
	      AccRejPlot->GetYaxis()->SetRangeUser(0.3, 1.02);
	      AccRejPlot->GetXaxis()->SetRangeUser(0.9, 1.01);
	      
	      TCanvas* cn = new TCanvas("cn","canvas");
	      cn->cd();
	      AccRejPlot->Draw("text");
	      //AccRejPlot->Draw("");
	      g->Draw("PL");
	      cms_myStyle(pu);
	      cn->SaveAs(outdir+ID[id]+part1[npart]+"_Pt20_"+label+"_QCD.pdf");
	      cn->SaveAs(outdir+ID[id]+part1[npart]+"_Pt20_"+label+"_QCD.png");
	    }
	  else 
	    {
	      for(int isB=0;isB<2;isB++)
		{
		  TH1F *h1,*h2;
		  TH2D *AccRejPlot = new TH2D("AccRejPlot","", 120, 0, 1.2, 120, 0 ,1.2);
		  TH2D *AccRejPlot2 = new TH2D("AccRejPlot2","", 120, 0, 1.2, 120, 0 ,1.2);
		  
		  TString h1Histo = "rel_iso"+ID[id]+match[0]+part1[npart]+"_"+iEB[isB]; // Signal
		  TString h2Histo = "rel_iso"+ID[id]+match[1]+part1[npart]+"_"+iEB[isB]; // QCD
		  h1 = (TH1F*) inFile1->Get(h1Histo);
		  h2 = (TH1F*) inFile2->Get(h2Histo);
		  std::cout<<h1Histo<<"," <<h2Histo<< std::endl;
		  //std::cout<<"got the histograms"<<std::endl;
		  std::cout<<"got the histograms: "<<h1->GetName()<< ", " << h2->GetName()<<std::endl;
		  TGraph *g;
		  int i = 0;
		  double gx[20], gy[20];
		  for ( float j = 0.0 ; j < 1 ; j += 0.05)
		    {
		      int bin = h1->FindBin(j);
		      double y = 1 - h2->Integral(0,bin-1)/(double)h2->GetEntries();
		      double x = h1->Integral(0,bin)/(double)h1->GetEntries();
		      
		      std::cout<<i<<std::endl;
		      gx[i]=x;
		      gy[i]=y;
		      
		      AccRejPlot->Fill(x,y,j);
		      AccRejPlot2->Fill(x,y);
		      i++;
		    }
		  
		  g = new TGraph(20, gx, gy);
		  std::cout<<i<<std::endl;
		  AccRejPlot->Divide(AccRejPlot2);
		  AccRejPlot->GetZaxis()->SetTitle("Isolation Value");
		  AccRejPlot->GetYaxis()->SetTitle("ID Rejection");
		  AccRejPlot->GetXaxis()->SetTitle("ID Efficiency");
		  AccRejPlot->GetZaxis()->SetRangeUser(0.0, 1.0);
		  AccRejPlot->SetStats(0);
		  //AccRejPlot->GetYaxis()->SetRangeUser(0.1, 1.02);
		  //AccRejPlot->GetXaxis()->SetRangeUser(0.30, 1.01);
		  
		  AccRejPlot->GetYaxis()->SetRangeUser(0.3, 1.02);
		  if(isB==0)
		AccRejPlot->GetXaxis()->SetRangeUser(0.9, 1.01);
		  if(isB==1)
		    AccRejPlot->GetXaxis()->SetRangeUser(0.4, 1.01);
		  
		  TCanvas* cn = new TCanvas("cn","canvas");
		  cn->cd();
		  AccRejPlot->Draw("text");
		  //AccRejPlot->Draw("");
		  g->Draw("PL");
		  cms_myStyle(pu);
		  
		  cn->SaveAs(outdir+ID[id]+part1[npart]+"_Pt20_"+iEB[isB]+label+"_QCD.pdf");
		  cn->SaveAs(outdir+ID[id]+part1[npart]+"_Pt20_"+iEB[isB]+label+"_QCD.png");
		  //cn->SaveAs(outdir+"MediumElec_Pt20"+iEB[isB]+"_"+label+"QCD.pdf");
		  //cn->SaveAs(outdir+"MediumElec_Pt20"+iEB[isB]+"_"+label+"QCD.png");
		}
	    }
	}
    }
  inFile1->Close();
  inFile2->Close();
  delete inFile1;
  delete inFile2;
  
  return 0;
}
    


void plotIt_MakePlots()
{
  plotIt_Eff("../Step2_PostAN/DYToLL_M-50_14TeV_TuneCP5_pythia8.root", "/eos/user/s/sandhya/www/RTB/Iter1/","200",1,0);
  plotIt_Eff("../Step2_PostAN/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8.root", "/eos/user/s/sandhya/www/RTB/Iter1/","200",0,0);
  plotIt_ROC("../Step2_PostAN/DYToLL_M-50_14TeV_TuneCP5_pythia8.root","../Step2_PostAN/QCD_Pt-15To7000_TuneCP5_Flat_14TeV-pythia8.root", "/eos/user/s/sandhya/www/RTB/Iter1/","200");
}
