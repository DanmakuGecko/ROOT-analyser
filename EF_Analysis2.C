#define EF_Analysis2_cxx
#include "EF_Analysis2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <iostream>
#include <TH1F.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include <iostream>
#include "TCanvas.h"
#include "TROOT.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooPlot.h"
#include <RooPlot.h> // Inclure l'en-tête de la classe RooPlot
#include <TColor.h> // Inclure l'en-tête pour la classe TColor de ROOT
#include "TColor.h"



#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <RooDataSet.h>

using namespace RooFit;





void EF_Analysis2::Loop()
{
//   In a ROOT session, you can do:
//      root> .L EF_Analysis2.C
//      root> EF_Analysis2 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   //cout << "IN1" << endl ;
   if (fChain == 0) return;
   //cout <<"IN2" << endl ;

   Long64_t nentries = fChain->GetEntriesFast();
   cout << "nentries = " << nentries << endl ;	
   Long64_t nbytes = 0, nb = 0;
   //return ;
   
   // ne lisons que ce qui est nÃ©cÃ©ssaire
   fChain->SetBranchStatus("*",0);  // disable all branches	
   fChain->SetBranchStatus("B_PX",1);  // activate branchname
   fChain->SetBranchStatus("Topo_xgb",1);
   fChain->SetBranchStatus("B_PY",1);
   fChain->SetBranchStatus("B_PZ",1);
   fChain->SetBranchStatus("B_PE",1);
   fChain->SetBranchStatus("B_M_pipKS",1);
   fChain->SetBranchStatus("B_M_ppiKS",1);
   fChain->SetBranchStatus("B_M_pipiKS",1);
   fChain->SetBranchStatus("B_M_KKKS",1);
   fChain->SetBranchStatus("B_MM",1);
   fChain->SetBranchStatus("B_M_pipiKS",1);
   fChain->SetBranchStatus("h2_PROBNNp",1);
   fChain->SetBranchStatus("h1_PROBNNp",1);
   fChain->SetBranchStatus("h2_PROBNNK",1);
   fChain->SetBranchStatus("h1_PROBNNK",1);
   fChain->SetBranchStatus("h2_PROBNNpi",1);
   fChain->SetBranchStatus("h1_PROBNNpi",1);
   

   
   fChain->SetBranchStatus("m13Sq_ppiKS",1);
   fChain->SetBranchStatus("m12Sq_ppiKS",1);
   fChain->SetBranchStatus("m23Sq_ppiKS",1);
   
   
   fChain->SetBranchStatus("m13Sq_pipKS",1);
   fChain->SetBranchStatus("m12Sq_pipKS",1);
   fChain->SetBranchStatus("m23Sq_pipKS",1);
   
   TH1I * histogram1 = new TH1I ("Masse invariante", "Candidates pipKS Topo > 0.95", 250, 5090, 7000);
   TH1I * histogram2 = new TH1I ("Masse invariante", "candidates pipKS", 67, 5400, 6200);
   TH1I * histogram3 = new TH1I ("Masse invariante", "pks", 150, 1300, 5650);

   
   

   
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=0; jentry<10000000;jentry++) {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);
   nbytes += nb;
   
   
   
   if (Topo_xgb > 0.95) {
       histogram1->Fill(B_M_pipKS);
       histogram1->Fill(B_M_ppiKS);
 
   }
   
   if (Topo_xgb > 0.95 && abs(B_M_pipiKS-5279.)>15. && abs(B_M_KKKS-5279.)>15.) {
	if ((B_M_pipKS>5400.)) {
  		//histogram1->Fill(sqrt(abs(B_PE*B_PE-(B_PX*B_PX+B_PY*B_PY+B_PZ*B_PZ))));
  		//if (h2_PROBNNp>0.95 && h1_PROBNNpi>0.95) {//if (abs(B_M_pipiKS-5279.62)>30) {
  		if ((h2_PROBNNp>0.95 && h2_PROBNNpi<0.1 && h2_PROBNNK<0.2) && (h1_PROBNNpi>0.95 && h1_PROBNNp<0.1 && h1_PROBNNK<0.2)) {
  		
  		histogram2->Fill(B_M_pipKS);
  		histogram3->Fill(sqrt(m23Sq_pipKS));
  		}
  	}
  	if ((B_M_ppiKS>5400.)) {
  		if ((h1_PROBNNp>0.95 && h1_PROBNNpi<0.1 && h1_PROBNNK<0.2) && (h2_PROBNNpi>0.95 && h2_PROBNNp<0.1 && h2_PROBNNK<0.2)) {//if (abs(B_M_pipiKS-5279.62)>40) {
  		//if (h1_PROBNNp>0.95 && h2_PROBNNpi>0.95)
  		
  		histogram2->Fill(B_M_ppiKS);
  		histogram3->Fill(sqrt(m13Sq_ppiKS));
  		}
  	}
  	
      //cout << "Topo_xgb =" << Topo_xgb << endl ;

   }

   //if ((jentry % 10) == 0) cout << "jentry = " << jentry << endl;

   
   
   
   
   if ((jentry%10000)==0) cout << "jentry = " << jentry << endl ;
      ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout << "h1_MM =" << h1_MM << endl ;  
      // if (Cut(ientry) < 0) continue;
   }
   TCanvas* canvas1 = new TCanvas("canvas1", "Histogram1", 800, 600);
   canvas1->Update();
   canvas1->Draw();
   histogram1->Draw("EP");

   TCanvas* canvas2 = new TCanvas("canvas2", "Histogram2", 800, 600);
   canvas2->Update();
   canvas2->Draw();
   histogram2->Draw("EP");
   
   TCanvas* canvas3 = new TCanvas("canvas3", "Histogram3", 800, 600);
   canvas3->Update();
   canvas3->Draw();
   histogram3->Draw();

}


void runArgusModel()
{
   RooRealVar x("x","MeV",5400.,6200.) ;
   RooRealVar mean_Lb("mean_Lb","mean of Lb distribution",5620.,5610.,5630.) ;
   RooRealVar sigma_Lb("sigma_Lb","width of Lb distribution",15.28) ;
   RooRealVar alpha_Lb("alpha_Lb", "alpha_Lb", 2.);
   RooRealVar n_Lb("n_Lb", "n_Lb", 1.);

   RooRealVar mean_Xib("mean_Xib","mean of Xib distribution",5791.) ;
   RooRealVar massdiffrence_XibLb("massdiffrence_XibLb", "massdiffrence_XibLb", 171.);

   RooRealVar mean_PRLb("mean_PRLb","mean of Lb distribution",5450.5) ;
   RooRealVar sigma_PRLb("sigma_PRLb","width of Lb distribution",30.) ;
   RooRealVar alpha_PRLb("alpha_PRLb", "alpha_PRLb", 0.35);

   RooFormulaVar mean_PRXib("mean_PRXib","mean of Xib distribution","mean_Xib-massdiffrence_XibLb",RooArgList(mean_Xib, massdiffrence_XibLb)) ;
   RooRealVar sigma_PRXib("sigma_PRXib","width of Xib distribution",30.) ;

         RooRealVar a1("a1","a1",5.);
         RooRealVar a2("a2","a2",5.);
         RooRealVar a3("a3","a3",3.);
         RooBernstein bg_bern("bg_bern","background",x,RooArgList(a1,a2,a3));



      

   RooRealVar frac("frac","fraction of Xib w.r.t. Lb", 0.05, 0.01, 0.3) ; 

   RooRealVar nLb("nLb", "#signal Lb events", 1000., 0., 10000);
   RooRealVar nXib("nXib", "#signal Xib events", 90., 0., 10000);
   RooRealVar nPRLb("nPRLb", "#PRLb", 270., 0., 10000);
   RooRealVar nPRXib("nPRXib", "#PRXib", 1., 0., 10000);
   

   

   RooRealVar c1("c1","c1", -0.5, -1, 1); 
   RooChebychev  cheb_poly("cheb_poly","cheb_poly", x, RooArgList(c1));
   
   
   


   
   
   RooRealVar nsig("nsig", "#signal events", 9500, 0., 10000);
   RooRealVar ncomb("ncomb", "#background events", 900, 0., 10000);
   
   




   RooCBShape lambdab("lambdab", "crystal ball PDF", x, mean_Lb, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape Xib("Xib", "crystal ball PDF", x, mean_Xib, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape PRLb("PRLb", "crystal ball PDF", x, mean_PRLb, sigma_PRLb, alpha_PRLb, n_Lb);
   RooCBShape PRXib("PRXib", "crystal ball PDF", x, mean_PRXib, sigma_PRLb, alpha_PRLb, n_Lb);

   //RooAddPdf signal_model("sgnal_model", "Lb+Xib", {lambdab, Xib}, {nLb, nXib});
   RooAddPdf model("model","all",{bg_bern}, {ncomb});
   
   //RooAddPdf model1("model1","all",{lambdab}, {nLb});
   //RooAddPdf model2("model2","all",{Xib}, {nXib});
   //RooAddPdf model3("model3","all",{PRLb}, {nPRLb});
   //RooAddPdf model4("model4","all",{cheb_poly}, {ncomb});
   //RooAddPdf model5("model5","all",{PRXib}, {nPRXib});


   RooDataSet* data = model.generate(x,10000) ;
   //RooDataSet* data2 = model1.generate(x,10000) ;
   //RooDataSet* data1 = model5.generate(x,10000) ;

   model.fitTo(*data ) ;
   //model1.fitTo(*data2 ) ;
   //model5.fitTo(*data1 ) ;

   //RooDataSet* data3 = model2.generate(x,10000) ;
   //RooDataSet* data4 = model3.generate(x,10000) ;
   //RooDataSet* data5 = model4.generate(x,10000) ;
   //model4.fitTo(*data5 ) ;

   //model2.fitTo(*data3 ) ;
   //model3.fitTo(*data4 ) ;


   RooPlot* xframe = x.frame() ;
   
   data->plotOn(xframe) ;
   model.plotOn(xframe) ;
   
   //data2->plotOn(xframe) ;
   //lambdab.plotOn(xframe, LineColor(kRed)) ;
   
   //data3->plotOn(xframe) ;
   //model2.plotOn(xframe, LineColor(kGreen)) ;

   //data4->plotOn(xframe) ;
   //model3.plotOn(xframe, LineColor(kPink)) ;

   //data5->plotOn(xframe) ;
   //model4.plotOn(xframe, LineColor(kPink)) ;
   
   //data1->plotOn(xframe) ;
   //model5.plotOn(xframe, LineColor(kYellow)) ;
   
   
   xframe->Draw() ;
}

void EF_Analysis2::Loop2()
{
//   In a ROOT session, you can do:
//      root> .L EF_Analysis2.C
//      root> EF_Analysis2 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   //cout << "IN1" << endl ;
   if (fChain == 0) return;
   //cout <<"IN2" << endl ;

   Long64_t nentries = fChain->GetEntriesFast();
   cout << "nentries = " << nentries << endl ;	
   Long64_t nbytes = 0, nb = 0;
   //return ;
   
   // ne lisons que ce qui est nÃ©cÃ©ssaire
   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("B_PX",1);  // activate branchname
   fChain->SetBranchStatus("Topo_xgb",1);
   fChain->SetBranchStatus("B_PY",1);
   fChain->SetBranchStatus("B_PZ",1);
   fChain->SetBranchStatus("B_PE",1);
   fChain->SetBranchStatus("B_M_pipKS",1);
   fChain->SetBranchStatus("B_M_KKKS",1);
   fChain->SetBranchStatus("B_M_piKKS",1);
   fChain->SetBranchStatus("B_M_KpiKS",1);
   fChain->SetBranchStatus("B_M_ppiKS",1);
   fChain->SetBranchStatus("B_MM",1);
   fChain->SetBranchStatus("B_M_pipiKS",1);
   fChain->SetBranchStatus("h2_PROBNNp",1);
   fChain->SetBranchStatus("h1_PROBNNp",1);
   fChain->SetBranchStatus("h2_PROBNNK",1);
   fChain->SetBranchStatus("h1_PROBNNK",1);
   fChain->SetBranchStatus("h2_PROBNNpi",1);
   fChain->SetBranchStatus("h1_PROBNNpi",1);
   

   
   fChain->SetBranchStatus("m13Sq_ppiKS",1);
   fChain->SetBranchStatus("m12Sq_ppiKS",1);
   fChain->SetBranchStatus("m23Sq_ppiKS",1);
   
   
   fChain->SetBranchStatus("m13Sq_pipKS",1);
   fChain->SetBranchStatus("m12Sq_pipKS",1);
   fChain->SetBranchStatus("m23Sq_pipKS",1);
   
   TH1I * histogram1 = new TH1I ("Masse invariante", "Masse invariante B_MM Topo>0.9", 250, 4800, 6000);
   TH1I * histogram2 = new TH1I ("Masse invariante", "candidates pipiKS",100, 5100, 5400);
   TH1I * histogram3 = new TH1I ("Masse invariante", "candidates piKKS", 100, 5100, 5480);
   TH1I * histogram4 = new TH1I ("Masse invariante", "candidates KKKS", 100, 5100, 5400);
   
   

   
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=0; jentry<10000000;jentry++) {
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);
   nbytes += nb;

   if (Topo_xgb > 0.95) {
	if ((B_M_pipKS>5400.)) {
  		histogram1->Fill(sqrt(abs(B_PE*B_PE-(B_PX*B_PX+B_PY*B_PY+B_PZ*B_PZ))));
  		//if (h2_PROBNNp>0.95 && h1_PROBNNpi>0.95) {//if (abs(B_M_pipiKS-5279.62)>30) {
  		if ((h2_PROBNNp>0.95 && h2_PROBNNpi<0.1 && h2_PROBNNK<0.2) && (h1_PROBNNpi>0.95 && h1_PROBNNp<0.1 && h1_PROBNNK<0.2)) {
  		histogram2->Fill(B_M_pipiKS);
  		histogram3->Fill(B_M_piKKS);
  		histogram4->Fill(B_M_KKKS);
  		}
  	}
  	if ((B_M_ppiKS>5400.)) {
  		if ((h1_PROBNNp>0.95 && h1_PROBNNpi<0.1 && h1_PROBNNK<0.2) && (h2_PROBNNpi>0.95 && h2_PROBNNp<0.1 && h2_PROBNNK<0.2)) {//if (abs(B_M_pipiKS-5279.62)>40) {
  		//if (h1_PROBNNp>0.95 && h2_PROBNNpi>0.95)
  		histogram2->Fill(B_M_pipiKS);
  		histogram3->Fill(B_M_KpiKS);
  		histogram4->Fill(B_M_KKKS);
  		}
  	}
  	
      //cout << "Topo_xgb =" << Topo_xgb << endl ;

   }

   //if ((jentry % 10) == 0) cout << "jentry = " << jentry << endl;

   
   
   
   
   if ((jentry%10000)==0) cout << "jentry = " << jentry << endl ;
      ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout << "h1_MM =" << h1_MM << endl ;  
      // if (Cut(ientry) < 0) continue;
   }
   TCanvas* canvas1 = new TCanvas("canvas1", "Histogram1", 800, 600);
   histogram1->Draw();
   canvas1->Update();
   canvas1->Draw();
   TCanvas* canvas2 = new TCanvas("canvas2", "Histogram2", 800, 600);
   canvas2->Update();
   canvas2->Draw();
   histogram2->Draw();
   
   TCanvas* canvas3 = new TCanvas("canvas3", "Histogram3", 800, 600);
   canvas3->Update();
   canvas3->Draw();
   histogram3->Draw();
   
   TCanvas* canvas4 = new TCanvas("canvas4", "Histogram4", 800, 600);
   canvas4->Update();
   canvas4->Draw();
   histogram4->Draw();

}







#include <TFile.h>
#include <TTree.h>

/*

void EF_Analysis2::copytree()
{

    

    Long64_t nentries = fChain->GetEntriesFast();
    cout << "Nombre de candidats total : " << nentries << endl;
    Long64_t nbytes = 0, nb = 0;

    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("m23Sq_pipKS", 1);
    fChain->SetBranchStatus("B_MM", 1);
    fChain->SetBranchStatus("m13Sq_ppiKS", 1);
    fChain->SetBranchStatus("h2_PROBNNp", 1);
    fChain->SetBranchStatus("h1_PROBNNpi", 1);
    fChain->SetBranchStatus("h1_PROBNNp", 1);
    fChain->SetBranchStatus("h2_PROBNNpi", 1);
    fChain->SetBranchStatus("Topo_xgb", 1);
    fChain->SetBranchStatus("B_M_pipiKS", 1);
    fChain->SetBranchStatus("h2_PROBNNp", 1);
    fChain->SetBranchStatus("h1_PROBNNp", 1);
    fChain->SetBranchStatus("h2_PROBNNK", 1);
    fChain->SetBranchStatus("h1_PROBNNK", 1);
    fChain->SetBranchStatus("h2_PROBNNpi", 1);
    fChain->SetBranchStatus("h1_PROBNNpi", 1);
    fChain->SetBranchStatus("B_M_pipKS", 1);
    fChain->SetBranchStatus("B_M_ppiKS", 1);

    TH1I *histogram = new TH1I("Invariant mass", "Candidates for ...", 80, 5400, 6200);

    RooRealVar x("x", "Invariant mass", 5400, 6200);
    RooDataHist dataHist("dataHist", "Data histogram", x, histogram);

    TTree* outputTree0 = new TTree("outputTree0", "Output Tree0");
    Float_t B_M_pipKS_out, B_M_ppiKS_out;
    outputTree0->Branch("B_M_pipKS", &B_M_pipKS_out, "B_M_pipKS/F");

    Long64_t jentry = 0;
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry);
    
    for (Long64_t jentry=0; jentry<10000000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
        
   if (Topo_xgb > 0.95 && abs(B_M_pipiKS-5279.)>15. && abs(B_M_KKKS-5279.)>15.) {
	if ((B_M_pipKS>5399. && B_M_pipKS<6201.)) {
  		//histogram1->Fill(sqrt(abs(B_PE*B_PE-(B_PX*B_PX+B_PY*B_PY+B_PZ*B_PZ))));
  		//if (h2_PROBNNp>0.95 && h1_PROBNNpi>0.95) {//if (abs(B_M_pipiKS-5279.62)>30) {
  		if ((h2_PROBNNp>0.95 && h2_PROBNNpi<0.1 && h2_PROBNNK<0.2) && (h1_PROBNNpi>0.95 && h1_PROBNNp<0.1 && h1_PROBNNK<0.2)) {

                    
                    B_M_pipKS_out = B_M_pipKS;
                    outputTree0->Fill();
                }
            }
            
  	if ((B_M_ppiKS>5399. && B_M_ppiKS<6201.)) {
  		if ((h1_PROBNNp>0.95 && h1_PROBNNpi<0.1 && h1_PROBNNK<0.2) && (h2_PROBNNpi>0.95 && h2_PROBNNp<0.1 && h2_PROBNNK<0.2)) {

                     
                    B_M_pipKS_out = B_M_ppiKS;
                    outputTree0->Fill();
                }
            }
        }
        
        if ((jentry % 10000) == 0) cout << "jentry = " << jentry << endl;
        ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
    }



    TFile outputFile("output0.root", "RECREATE");
    outputTree0->Write();
    outputFile.Close();
}

*/
void roodataset()
{
 //LIRE LES FICHIERS 
 //-----------------

 TFile *file = TFile::Open("pKKS.root", "READ");
 TTree* tree = (TTree*)file->Get("outputTree0");
 
 //Déclarer les variables : 

 Float_t B_M_KpKS; 
 TBranch* b_B_M_KpKS;
 

 tree->SetBranchAddress("B_M_KpKS", &B_M_KpKS, &b_B_M_KpKS); //On lie la branchye a la variable

 RooRealVar x("B_M_KpKS", "B_M_KpKS", 5400., 6150.); 
 int nbins=50;

 RooDataSet data("data", "data", RooArgSet(x)); //on va utiliser RooArgSet pour définir notre dataset
 Long64_t nbEntrees = tree->GetEntries();
 for (Long64_t i = 0; i < nbEntrees; i++) {

      b_B_M_KpKS->GetEntry(i);

      cout<<B_M_KpKS<<endl;
      //if (abs(B_M_pipKS-5447.)>1. && abs(B_M_pipKS-5980.)>2. ){   //Nouvelle condition pour enelver l'artéfact
          x.setVal(B_M_KpKS);
          data.add(RooArgSet(x));
   //}
   } 



     RooRealVar mean_Lb("mean_Lb","mean of Lb distribution",5618.2) ;
   RooRealVar sigma_Lb("sigma_Lb","width of Lb distribution",15.28) ;
   RooRealVar alpha_Lb("alpha_Lb", "alpha_Lb", 2.);
   RooRealVar n_Lb("n_Lb", "n_Lb", 1.);

   RooRealVar mean_Xib("mean_Xib","mean of Xib distribution",5791.) ;
   RooRealVar massdiffrence_XibLb("massdiffrence_XibLb", "massdiffrence_XibLb", 171.);

   RooRealVar mean_PRLb("mean_PRLb","mean of Lb distribution",5450.5) ;
   RooRealVar sigma_PRLb("sigma_PRLb","width of Lb distribution",30.) ;
   RooRealVar alpha_PRLb("alpha_PRLb", "alpha_PRLb", 0.35);

   RooFormulaVar mean_PRXib("mean_PRXib","mean of Xib distribution","mean_Xib-massdiffrence_XibLb",RooArgList(mean_Xib, massdiffrence_XibLb)) ;
   RooRealVar sigma_PRXib("sigma_PRXib","width of Xib distribution",30.) ;


         RooRealVar a1("a1","a1",5.);
         RooRealVar a2("a2","a2",5.);
         RooRealVar a3("a3","a3",3.);
   

         RooBernstein bg_bern("bg_bern","background",x,RooArgList(a1,a2,a3));


      

   RooRealVar frac("frac","fraction of Xib w.r.t. Lb", 0.05, 0.01, 0.3) ; 

   RooRealVar nLb("nLb", "#signal Lb events", 1000., 0., 10000);
   RooRealVar nXib("nXib", "#signal Xib events", 90., 0., 10000);
   RooFormulaVar nPRLb("nPRXib","nPRXib","nLb/(400/66)", RooArgList(nLb));
   //RooRealVar nPRLb("nPRLb", "#PRLb", 270., 0., 10000);
   RooFormulaVar nPRXib("nPRXib", "#PRXib","nXib/(400/66)",RooArgList(nXib) );
   

   

   //RooRealVar c1("c1","c1", -0.5, -1.03, 1.03); 
   //RooChebychev  cheb_poly("cheb_poly","cheb_poly", x, RooArgList(c1));


   
   
   RooRealVar nsig("nsig", "#signal events", 9500, 0., 10000);
   RooRealVar ncomb("ncomb", "#background events", 900, 0., 10000);
   
   




   RooCBShape lambdab("lambdab", "crystal ball PDF", x, mean_Lb, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape Xib("Xib", "crystal ball PDF", x, mean_Xib, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape PRLb("PRLb", "crystal ball PDF", x, mean_PRLb, sigma_PRLb, alpha_PRLb, n_Lb);
   RooCBShape PRXib("PRXib", "crystal ball PDF", x, mean_PRXib, sigma_PRLb, alpha_PRLb, n_Lb);

   //RooAddPdf signal_model("sgnal_model", "Lb+Xib", {lambdab, Xib}, {nLb, nXib});
   RooAddPdf model("model","all",{lambdab, Xib, PRLb,PRXib,bg_bern}, {nLb, nXib, nPRLb,nPRXib,ncomb});
   

 
 

 RooFitResult* result = model.fitTo(data, RooFit::Save()); //On fit les données du fichier avec notre modèle 

 result->Print();

RooPlot* frame = x.frame();
data.plotOn(frame, RooFit::Binning(nbins));
model.plotOn(frame);

// Ajoute une légende au graphique
TLegend* legend = new TLegend(0.6,0.65,0.9,0.9);
// Ajoute des entrées à la légende pour chaque composante du modèle
legend->AddEntry(frame->getObject(0), "Données", "lep"); // Ajoute une entrée pour les données
legend->AddEntry(frame->getObject(1), "Modèle complet", "l"); // Ajoute une entrée pour le modèle complet

// Ajoute des entrées pour chaque composante individuelle du modèle
model.plotOn(frame, Components(lambdab), LineStyle(kDashed), LineColor(kRed));
legend->AddEntry(frame->getObject(2), "Composante lambdab", "l");

model.plotOn(frame, Components(Xib), LineStyle(kDashed), LineColor(kBlack));
legend->AddEntry(frame->getObject(3), "Composante Xib", "l");

model.plotOn(frame, Components(PRLb), LineStyle(kDashed), LineColor(kGreen));
legend->AddEntry(frame->getObject(4), "Composante PRLb", "l");

model.plotOn(frame, Components(PRXib), LineStyle(kDashed), LineColor(kOrange));
legend->AddEntry(frame->getObject(5), "Composante PRXib", "l");

model.plotOn(frame, Components(bg_bern), LineStyle(kDashed), LineColor(kYellow));
legend->AddEntry(frame->getObject(6), "Composante bg_bern", "l");

// Dessine le graphique avec la légende
frame->Draw();
legend->Draw();
}


void H0()

{
    //Création du vecteur pour les nll
    
    vector<double_t> nlls; //Vecteur pour les ordonnées
    //double nXib = 0;

		

	   //LIRE LES FICHIERS 
	   //-----------------

	   TFile *file = TFile::Open("pKKS.root", "READ");
	   TTree* tree = (TTree*)file->Get("outputTree0");
	   
	   //Déclarer les variables : 

	   Float_t B_M_KpKS;  
	   TBranch* b_B_M_KpKS;

	   tree->SetBranchAddress("B_M_KpKS", &B_M_KpKS, &b_B_M_KpKS); //On lie la branche à la variable

	   RooRealVar x("B_M_KpKS", "B_M_KpKS", 5400., 6200.);
   


	   RooDataSet data("data", "data", RooArgSet(x)); //on va utiliser RooArgSet pour définir notre dataset
	   Long64_t nbEntrees = tree->GetEntries();
	   
	   for (Long64_t i = 0; i < nbEntrees; i++) {

	      b_B_M_KpKS->GetEntry(i);
	      
	      //cout<<B_M_pipKS<<endl;
	      //if (abs(B_M_KpKS-5447)>1 && abs(B_M_KpKS-5980)>2 ){
		      x.setVal(B_M_KpKS);
		      data.add(RooArgSet(x));
	   //}
	   }
	   
	   
	   
		  
	    //Destruction des artéfacts :

	    
	    
			  
		  
  RooRealVar mean_Lb("mean_Lb","mean of Lb distribution",5618.2) ;
   RooRealVar sigma_Lb("sigma_Lb","width of Lb distribution",15.28) ;
   RooRealVar alpha_Lb("alpha_Lb", "alpha_Lb", 2.);
   RooRealVar n_Lb("n_Lb", "n_Lb", 1.);

   RooRealVar mean_Xib("mean_Xib","mean of Xib distribution",5791.) ;
   RooRealVar massdiffrence_XibLb("massdiffrence_XibLb", "massdiffrence_XibLb", 171.);

   RooRealVar mean_PRLb("mean_PRLb","mean of Lb distribution",5450.5) ;
   RooRealVar sigma_PRLb("sigma_PRLb","width of Lb distribution",30.) ;
   RooRealVar alpha_PRLb("alpha_PRLb", "alpha_PRLb", 0.35);

   RooFormulaVar mean_PRXib("mean_PRXib","mean of Xib distribution","mean_Xib-massdiffrence_XibLb",RooArgList(mean_Xib, massdiffrence_XibLb)) ;
   RooRealVar sigma_PRXib("sigma_PRXib","width of Xib distribution",30.) ;




         RooRealVar a1("a1","a1",5.);
         RooRealVar a2("a2","a2",5.);
         RooRealVar a3("a3","a3",3.);
   

         RooBernstein bg_bern("bg_bern","background",x,RooArgList(a1,a2,a3));



      

   RooRealVar frac("frac","fraction of Xib w.r.t. Lb", 0.05, 0.01, 0.3) ; 

   RooRealVar nLb("nLb", "#signal Lb events", 1000., 0., 10000);
   RooRealVar nXib("nXib", "#signal Xib events", 0.);
   RooFormulaVar nPRLb("nPRXib","nPRXib","nLb/(400/66)", RooArgList(nLb));
   //RooRealVar nPRLb("nPRLb", "#PRLb", 270., 0., 10000);
   RooFormulaVar nPRXib("nPRXib", "#PRXib","nXib/(400/66)",RooArgList(nXib) );
   

   

   //RooRealVar c1("c1","c1", -0.5, -1.03, 1.03); 
   //RooChebychev  cheb_poly("cheb_poly","cheb_poly", x, RooArgList(c1));


   
   
   RooRealVar nsig("nsig", "#signal events", 9500, 0., 10000);
   RooRealVar ncomb("ncomb", "#background events", 900, 0., 10000);
   
   




   RooCBShape lambdab("lambdab", "crystal ball PDF", x, mean_Lb, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape Xib("Xib", "crystal ball PDF", x, mean_Xib, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape PRLb("PRLb", "crystal ball PDF", x, mean_PRLb, sigma_PRLb, alpha_PRLb, n_Lb);
   RooCBShape PRXib("PRXib", "crystal ball PDF", x, mean_PRXib, sigma_PRLb, alpha_PRLb, n_Lb);

   //RooAddPdf signal_model("sgnal_model", "Lb+Xib", {lambdab, Xib}, {nLb, nXib});
   RooAddPdf model("model","all",{lambdab, Xib, PRLb,PRXib,bg_bern}, {nLb, nXib, nPRLb,nPRXib,ncomb});
   


	   //Calcul de la log-vraisemblance
	   
	   RooAbsReal* nll = model.createNLL(data, RooFit::Extended(true)); //variable liée a nll
	   RooMinimizer(*nll).migrad() ;
	   Double_t nllValue = -nll->getVal();
	   

	   
	    
	   RooFitResult* result = model.fitTo(data, RooFit::Save()); //On fit les données du fichier avec notre modèle 

	   result->Print();
	   
	   
	   
	   RooPlot* frame = x.frame();
	   data.plotOn(frame, RooFit::Binning(66));
	   model.plotOn(frame);
	  
	   frame->Draw();
	   
	   cout << "nll =  " << nllValue << endl;
	   nlls.push_back(nllValue);
	   

}

void H1()

{
    //Création du vecteur pour les nll
    
    vector<double_t> nlls; //Vecteur pour les ordonnées

		

	   //LIRE LES FICHIERS 
	   //-----------------

	   TFile *file = TFile::Open("pKKS.root", "READ");
	   TTree* tree = (TTree*)file->Get("outputTree0");
	   
	   //Déclarer les variables : 

	   Float_t B_M_KpKS;  
	   TBranch* b_B_M_KpKS;

	   tree->SetBranchAddress("B_M_KpKS", &B_M_KpKS, &b_B_M_KpKS); //On lie la branche à la variable

	   RooRealVar x("B_M_KpKS", "B_M_KpKS", 5400., 6200.);
   


	   RooDataSet data("data", "data", RooArgSet(x)); //on va utiliser RooArgSet pour définir notre dataset
	   Long64_t nbEntrees = tree->GetEntries();
	   
	   for (Long64_t i = 0; i < nbEntrees; i++) {

	      b_B_M_KpKS->GetEntry(i);
	      
	      //cout<<B_M_pipKS<<endl;
	      //if (abs(B_M_KpKS-5447)>1 && abs(B_M_KpKS-5980)>2 ){
		      x.setVal(B_M_KpKS);
		      data.add(RooArgSet(x));
	   //}
	   }
	   
	   
	   
		  
	    //Destruction des artéfacts :

	    
	    
			  
		  
  RooRealVar mean_Lb("mean_Lb","mean of Lb distribution",5618.2) ;
   RooRealVar sigma_Lb("sigma_Lb","width of Lb distribution",15.28) ;
   RooRealVar alpha_Lb("alpha_Lb", "alpha_Lb", 2.);
   RooRealVar n_Lb("n_Lb", "n_Lb", 1.);

   RooRealVar mean_Xib("mean_Xib","mean of Xib distribution",5791.) ;
   RooRealVar massdiffrence_XibLb("massdiffrence_XibLb", "massdiffrence_XibLb", 171.);

   RooRealVar mean_PRLb("mean_PRLb","mean of Lb distribution",5450.5) ;
   RooRealVar sigma_PRLb("sigma_PRLb","width of Lb distribution",30.) ;
   RooRealVar alpha_PRLb("alpha_PRLb", "alpha_PRLb", 0.35);

   RooFormulaVar mean_PRXib("mean_PRXib","mean of Xib distribution","mean_Xib-massdiffrence_XibLb",RooArgList(mean_Xib, massdiffrence_XibLb)) ;
   RooRealVar sigma_PRXib("sigma_PRXib","width of Xib distribution",30.) ;


         RooRealVar a1("a1","a1",5.);
         RooRealVar a2("a2","a2",5.);
         RooRealVar a3("a3","a3",3.);
   

         RooBernstein bg_bern("bg_bern","background",x,RooArgList(a1,a2,a3));





      

   RooRealVar frac("frac","fraction of Xib w.r.t. Lb", 0.05, 0.01, 0.3) ; 

   RooRealVar nLb("nLb", "#signal Lb events", 1000., 0., 10000);
   RooRealVar nXib("nXib", "#signal Xib events", 90., 0., 10000);
   RooFormulaVar nPRLb("nPRXib","nPRXib","nLb/(400/66)", RooArgList(nLb));
   //RooRealVar nPRLb("nPRLb", "#PRLb", 270., 0., 10000);
   RooFormulaVar nPRXib("nPRXib", "#PRXib","nXib/(400/66)",RooArgList(nXib) );
   

   

   //RooRealVar c1("c1","c1", -0.5, -1.03, 1.03); 
   //RooChebychev  cheb_poly("cheb_poly","cheb_poly", x, RooArgList(c1));


   
   
   RooRealVar nsig("nsig", "#signal events", 9500, 0., 10000);
   RooRealVar ncomb("ncomb", "#background events", 900, 0., 10000);
   
   




   RooCBShape lambdab("lambdab", "crystal ball PDF", x, mean_Lb, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape Xib("Xib", "crystal ball PDF", x, mean_Xib, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape PRLb("PRLb", "crystal ball PDF", x, mean_PRLb, sigma_PRLb, alpha_PRLb, n_Lb);
   RooCBShape PRXib("PRXib", "crystal ball PDF", x, mean_PRXib, sigma_PRLb, alpha_PRLb, n_Lb);

   //RooAddPdf signal_model("sgnal_model", "Lb+Xib", {lambdab, Xib}, {nLb, nXib});
   RooAddPdf model("model","all",{lambdab, Xib, PRLb,PRXib,bg_bern}, {nLb, nXib, nPRLb,nPRXib,ncomb});
	   //Calcul de la log-vraisemblance
	   
	   RooAbsReal* nll = model.createNLL(data, RooFit::Extended(true)); //variable liée a nll
	   RooMinimizer(*nll).migrad() ;
	   Double_t nllValue = -nll->getVal();
	   

	   
	    
	   RooFitResult* result = model.fitTo(data, RooFit::Save()); //On fit les données du fichier avec notre modèle 

	   result->Print();
	   
	   
	   
	   RooPlot* frame = x.frame();
	   data.plotOn(frame, RooFit::Binning(66));
	   model.plotOn(frame);
	  
	   frame->Draw();
	   
	   cout << "nll =  " << nllValue << endl;
	   nlls.push_back(nllValue);
	   
   

}

#include <TFile.h>
#include <TTree.h>


void copySelectedBranches()
{
    
    // Ouvrir le fichier ROOT en lecture
    TFile* inputFile = new TFile("B2KShh-Collision2018-Stripping34-DD-B2pipiKS-TrigandPresel-XGBs.root", "READ");

    // Créer un nouveau fichier ROOT en écriture
    TFile* outputFile = new TFile("output.root", "RECREATE");

    // Obtenir l'arbre d'origine à partir du fichier d'entrée
    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get("B2pipiKS")); 

    // Cloner l'arbre d'origine dans le fichier de sortie en ne copiant que les branches nécessaires
    TTree* outputTree = inputTree->CloneTree(0);

    // Ajouter les branches nécessaires à l'arbre de sortie
    outputTree->SetBranchStatus("*", 0); // Désactive toutes les branches par défaut
    outputTree->SetBranchStatus("B_M_pipKS", 1);


    // Copier les entrées sélectionnées dans l'arbre de sortie
    Long64_t numEntries = inputTree->GetEntries();
    for (Long64_t iEntry = 0; iEntry < numEntries; ++iEntry) {
        inputTree->GetEntry(iEntry);
        
            outputTree->Fill();
       if ((iEntry % 10000) == 0) cout << "jentry = " << iEntry << endl;

    }

    // Écrire l'arbre sélectionné dans le fichier de sortie
    outputFile->Write();
    outputFile->Close();

    // Fermer le fichier d'entrée
    inputFile->Close();
}



void EF_Analysis2::copytree()
{

    
    Long64_t nentries = fChain->GetEntriesFast();
    cout << "Nombre de candidats total : " << nentries << endl;
    Long64_t nbytes = 0, nb = 0;

   fChain->SetBranchStatus("*",0);  // disable all branches	
   fChain->SetBranchStatus("B_PX",1);  // activate branchname
   fChain->SetBranchStatus("Topo_xgb",1);
   fChain->SetBranchStatus("B_PY",1);
   fChain->SetBranchStatus("B_PZ",1);
   fChain->SetBranchStatus("B_PE",1);
   fChain->SetBranchStatus("B_M_KpKS",1);
   fChain->SetBranchStatus("B_M_pKKS",1);
   fChain->SetBranchStatus("B_M_pipiKS",1);
   fChain->SetBranchStatus("B_M_KKKS",1);
   fChain->SetBranchStatus("B_MM",1);
   fChain->SetBranchStatus("B_M_pipiKS",1);
   fChain->SetBranchStatus("h2_PROBNNp",1);
   fChain->SetBranchStatus("h1_PROBNNp",1);
   fChain->SetBranchStatus("h2_PROBNNK",1);
   fChain->SetBranchStatus("h1_PROBNNK",1);
   fChain->SetBranchStatus("h2_PROBNNpi",1);
   fChain->SetBranchStatus("h1_PROBNNpi",1);
   

   
   fChain->SetBranchStatus("m13Sq_ppiKS",1);
   fChain->SetBranchStatus("m12Sq_ppiKS",1);
   fChain->SetBranchStatus("m23Sq_ppiKS",1);
   
   
   fChain->SetBranchStatus("m13Sq_pipKS",1);
   fChain->SetBranchStatus("m12Sq_pipKS",1);
   fChain->SetBranchStatus("m23Sq_pipKS",1);

     TH1D *histogram1 = new TH1D("histogram1", "histogram1", 80, 5400, 6200);

    RooRealVar x("x", "Invariant mass", 5400, 6200);
    RooDataHist dataHist("dataHist", "Data histogram", x, histogram1);

    
    TTree* outputTree0 = new TTree("outputTree0", "Output Tree0");
    Float_t B_M_KpKS_out, B_M_pKKS_out;
    outputTree0->Branch("B_M_KpKS", &B_M_KpKS_out, "B_M_KpKS/F");

    Long64_t jentry = 0;
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry);
    
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        
   if (Topo_xgb > 0.95 && abs(B_M_pipiKS-5279.)>15. && abs(B_M_KKKS-5279.)>15.) {
	if ((B_M_KpKS>5400. && B_M_KpKS<6200.)) {
  		
  		if ((h2_PROBNNp>0.95 && h2_PROBNNpi<0.1 && h2_PROBNNK<0.2) && (h1_PROBNNK>0.95 && h1_PROBNNp<0.1 && h1_PROBNNpi<0.2)) {

                    
                    B_M_KpKS_out = B_M_KpKS;
                    outputTree0->Fill();
                    histogram1->Fill(B_M_KpKS);
                }
            }
            
  	if ((B_M_pKKS>5400. && B_M_pKKS<6200.)) {
  		if ((h1_PROBNNp>0.95 && h1_PROBNNpi<0.1 && h1_PROBNNK<0.2) && (h2_PROBNNK>0.95 && h2_PROBNNp<0.1 && h2_PROBNNpi<0.2)) {

                     
                    B_M_KpKS_out = B_M_pKKS;
                    outputTree0->Fill();
                    histogram1->Fill(B_M_pKKS);
                }
            }
        }
        
        if ((jentry % 10000) == 0) cout << "jentry = " << jentry << endl;
        ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
    }



    TFile outputFile("pKKS.root", "RECREATE");
    outputTree0->Write();
    outputFile.Close();
   TCanvas* canvas1 = new TCanvas("canvas1", "histogram1", 800, 600);
   canvas1->Update();
   canvas1->Draw();
   histogram1->Draw("EP");
}



void H3()

{
    //Création du vecteur pour les nll
    
    vector<double_t> nlls; //Vecteur pour les ordonnées
    


		

	   //LIRE LES FICHIERS 
	   //-----------------

	   TFile *file = TFile::Open("pKKS.root", "READ");
	   TTree* tree = (TTree*)file->Get("outputTree0");
	   
	   //Déclarer les variables : 

	   Float_t B_M_KpKS;  
	   TBranch* b_B_M_KpKS;

	   tree->SetBranchAddress("B_M_KpKS", &B_M_KpKS, &b_B_M_KpKS); //On lie la branche à la variable

	   RooRealVar x("B_M_KpKS", "B_M_KpKS", 5400., 6200.);
   


	   RooDataSet data("data", "data", RooArgSet(x)); //on va utiliser RooArgSet pour définir notre dataset
	   Long64_t nbEntrees = tree->GetEntries();
	   
	   for (Long64_t i = 0; i < nbEntrees; i++) {

	      b_B_M_KpKS->GetEntry(i);
	      
	      //cout<<B_M_pipKS<<endl;
	      //if (abs(B_M_KpKS-5447)>1 && abs(B_M_KpKS-5980)>2 ){
		      x.setVal(B_M_KpKS);
		      data.add(RooArgSet(x));
	   //}
	   }
	    
	    
			  
		  
RooRealVar mean_Lb("mean_Lb","mean of Lb distribution",5618.2) ;
   RooRealVar sigma_Lb("sigma_Lb","width of Lb distribution",15.28) ;
   RooRealVar alpha_Lb("alpha_Lb", "alpha_Lb", 2.);
   RooRealVar n_Lb("n_Lb", "n_Lb", 1.);

   RooRealVar mean_Xib("mean_Xib","mean of Xib distribution",5791.) ;
   RooRealVar massdiffrence_XibLb("massdiffrence_XibLb", "massdiffrence_XibLb", 171.);

   RooRealVar mean_PRLb("mean_PRLb","mean of Lb distribution",5450.5) ;
   RooRealVar sigma_PRLb("sigma_PRLb","width of Lb distribution",30.) ;
   RooRealVar alpha_PRLb("alpha_PRLb", "alpha_PRLb", 0.35);

   RooFormulaVar mean_PRXib("mean_PRXib","mean of Xib distribution","mean_Xib-massdiffrence_XibLb",RooArgList(mean_Xib, massdiffrence_XibLb)) ;
   RooRealVar sigma_PRXib("sigma_PRXib","width of Xib distribution",30.) ;


         RooRealVar a1("a1","a1",5.);
         RooRealVar a2("a2","a2",5.);
         RooRealVar a3("a3","a3",3.);
   

         RooBernstein bg_bern("bg_bern","background",x,RooArgList(a1,a2,a3));





      

   RooRealVar frac("frac","fraction of Xib w.r.t. Lb", 0.05, 0.01, 0.3) ; 

   RooRealVar nLb("nLb", "#signal Lb events", 0.);
   RooRealVar nXib("nXib", "#signal Xib events", 90., 0., 10000);
   RooFormulaVar nPRLb("nPRXib","nPRXib","nLb/(400/66)", RooArgList(nLb));
   //RooRealVar nPRLb("nPRLb", "#PRLb", 270., 0., 10000);
   RooFormulaVar nPRXib("nPRXib", "#PRXib","nXib/(400/66)",RooArgList(nXib) );
   

   

   //RooRealVar c1("c1","c1", -0.5, -1.03, 1.03); 
   //RooChebychev  cheb_poly("cheb_poly","cheb_poly", x, RooArgList(c1));


   
   
   RooRealVar nsig("nsig", "#signal events", 9500, 0., 10000);
   RooRealVar ncomb("ncomb", "#background events", 900, 0., 10000);
   
   




   RooCBShape lambdab("lambdab", "crystal ball PDF", x, mean_Lb, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape Xib("Xib", "crystal ball PDF", x, mean_Xib, sigma_Lb, alpha_Lb, n_Lb);
   RooCBShape PRLb("PRLb", "crystal ball PDF", x, mean_PRLb, sigma_PRLb, alpha_PRLb, n_Lb);
   RooCBShape PRXib("PRXib", "crystal ball PDF", x, mean_PRXib, sigma_PRLb, alpha_PRLb, n_Lb);

   //RooAddPdf signal_model("sgnal_model", "Lb+Xib", {lambdab, Xib}, {nLb, nXib});
   RooAddPdf model("model","all",{lambdab, Xib, PRLb,PRXib,bg_bern}, {nLb, nXib, nPRLb,nPRXib,ncomb});
   
	   //Calcul de la log-vraisemblance
	   
	   RooAbsReal* nll = model.createNLL(data, RooFit::Extended(true)); //variable liée a nll
	   RooMinimizer(*nll).migrad() ;
	   Double_t nllValue = -nll->getVal();
	   

	   
	    
	   RooFitResult* result = model.fitTo(data, RooFit::Save()); //On fit les données du fichier avec notre modèle 

	   result->Print();
	   
	   
	   
	   RooPlot* frame = x.frame();
	   data.plotOn(frame, RooFit::Binning(66));
	   model.plotOn(frame);
	  
	   frame->Draw();
	   
	   cout << "nll =  " << nllValue << endl;
	   nlls.push_back(nllValue);

}







