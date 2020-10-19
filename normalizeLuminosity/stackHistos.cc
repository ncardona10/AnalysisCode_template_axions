#include "PhenoAnalyzer.h"
#include <bits/stdc++.h>
#include "./Plots/MyHistograms.h"

using namespace std;

int main(int argc, char const *argv[])
{
  cout<<"creating output file"<<endl;
  TFile *HistoOutputFile = new TFile("/home/n.cardonac/AnalysisCode/PhenoAnalyzer_SUSY_VBF_Higgsino/test_stack.root", "RECREATE");


  TString inputFolder = "/home/n.cardonac/RunPhenoCodes/higgsino_nomjj/";

  vector<TString> folderNames = {
                                //  "m_n2_100_c1_75_n1_50",
                                //  "m_n2_400_c1_385_n1_370",
                                 "vv",
                                 "w+jets",
                                 "z+jets",
                                 "ttbar",
                                 };

  map<TString, TFile *> sfiles;

  cout<<"reading input files..."<<endl;

  for (int i = 0; (unsigned)i < folderNames.size(); i++)
  {
    sfiles[folderNames[i]] = new TFile(inputFolder + folderNames[i] + "/" + folderNames[i] + "_nocut-1.root","READ");
    cout<<"reading file: "<< inputFolder + folderNames[i] + "/" + folderNames[i] + "_nocut-1.root" <<endl;
  }

  vector<TString> dirNames = {"nLeptons",
			      "Jets",	
                              "MET",
                              "BJets",
                              "VBF",
                              "single_e_met_bjets_vbf",
                              "single_mu_met_bjets_vbf"};

  //create ouput folders in output root
  cout<<"creating output subdirectories"<<endl;
  HistoOutputFile->cd();
  map<TString, TDirectory *> outputDirs;
  for (int i = 0; (unsigned)i < dirNames.size(); i++)
  { 
    cout<<"subdirectories names "<< dirNames[i]<<endl;
    outputDirs[dirNames[i]] = HistoOutputFile->mkdir(dirNames[i]);
  }

  //histogram names
  vector<TString> histoNames = {"ptelectron",
                                "ptmuon",
                                "pttau",
                                "ptjet",
                                "etaelectron",
                                "etamuon",
                                "etatau",
                                "etajet",
                                "phielectron",
                                "phimuon",
                                "phitau",
                                "phijet",
                                "Mtelectron",
                                "Mtmuon",
                                "Mttau",
                                "Mtjet",
                                "Mjj",
                                "MET",
                                "DEta"
                                };

  /*

  for every folder( ww,wz,...):
    for every directory (nleptons,...):
      for every histogram in folderhistogram in folder:
        get histogram
        change name to folder
      plot 


  */

  for (int directory_i = 0; (unsigned)directory_i < dirNames.size(); directory_i++)
  {
    for (int histo_i = 0; (unsigned)histo_i < histoNames.size(); histo_i++)
    {

      cout << dirNames[directory_i] << " " << histoNames[histo_i] << endl;
      TObjArray histos;

      for (int folder_i = 0; (unsigned)folder_i < folderNames.size(); folder_i++)
      {
        cout<<"changing to "<< folderNames[folder_i]<< " directory and "<<dirNames[directory_i]<<" subfolder"  <<endl;
        sfiles[folderNames[folder_i]]->cd(dirNames[directory_i]);

        cout<<"getting histogram info"<<endl;
        TH1F *histo = (TH1F *)sfiles[folderNames[folder_i]]->Get(dirNames[directory_i] + "/" + histoNames[histo_i]);


        cout<<"setting title and name"<<endl;
        histo->SetTitle(folderNames[folder_i]);
        
        histo->SetName(folderNames[folder_i]);

        cout<<"adding to histogram list"<<endl;

        if(histoNames[histo_i] == "Mjj"){
          histo->Rebin(5);
        }

        histos.AddLast(histo);
      }


      outputDirs[dirNames[directory_i]]->cd();

      TCanvas *cl = new TCanvas(dirNames[directory_i] + " " + histoNames[histo_i],dirNames[directory_i] + " " +  histoNames[histo_i], 600, 500);

      // cl->Divide(2,2); //create subplots if needed

      string histoNameString = (string)histoNames[histo_i];
      Draw_Stacked(histos, (TPad *)cl->cd(0), true, histoNameString);
      
      //if((folderNames[i] == "m_n2_100_c1_80_n1_60")){
        
      //} 
      cl->Write();
    }
  }


  HistoOutputFile->Close();

  for (int i = 0; (unsigned)i < folderNames.size(); i++)
  {
    sfiles[folderNames[i]]->Close();
  }
  cout << "DONE" << endl;
  return 0;
}
