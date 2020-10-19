/*                             __
                           .--()°"."
 Author: Nathalia Cardona "|, . ,"
                           !_-(_\
*/

#include "PhenoAnalyzer.h"
#include "LeptonCounter.h"
#include "Cuts.h"
#include "Physics.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <bits/stdc++.h>
#include <utility>

using namespace std;

void writeCsv(int count, string path, string cut)
{
  ofstream outfile;
  outfile.open("/home/n.cardonac/AnalysisCode_Axions/PhenoAnalyzer/counts.csv", ios_base::app); // append instead of overwrite
  outfile << path << "," << cut << "," << count << "\n";
}

int main(int argc, char *argv[])
{
  cout << "Starting phenoanalyzer..." << endl;

  // Importing Delphes data
  TChain chain("Delphes");
  chain.Add(argv[1]);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

  // output file manager
  TFile *HistoOutputFile = new TFile(argv[2], "RECREATE");

  // directory to store the histograms
  TDirectory *noCutsDirectory = HistoOutputFile->mkdir("noCuts");
  TDirectory *vbfCutsDirectory = HistoOutputFile->mkdir("vbfCuts");
  TDirectory *photonCuts_mppDirectory = HistoOutputFile->mkdir("photonCuts_mpp");
  TDirectory *photonCuts_ptDirectory = HistoOutputFile->mkdir("photonCuts_pt");

  // get tree info
  vector<string> branches = {
      "Electron",
      "Muon",
      "Jet",
      "MissingET",
      "Photon"
      };

  map<string, TClonesArray *> branchDict;

  // create a dictionary with the branches
  for (int i = 0; (unsigned)i < branches.size(); i++)
  {
    TClonesArray *branch = treeReader->UseBranch(branches[i].c_str());
    branchDict[branches[i]] = branch;
  }

  // boolean mask to avoid over computation
  vector<int> cutsArr;

  // vector of pairs for jet pair that max mjj
  vector<pair<Jet *, Jet *>> maxMjjJetPairs;

  cout<<"creating pairs"<<endl;

  for (int i = 0; (unsigned)i < treeReader->GetEntries(); i++)
  {
    cutsArr.push_back(1);
    maxMjjJetPairs.push_back(getMaxMjjJetPair(treeReader, branchDict, i));
  }

  int nEvents;

  // write number of events to csv
  nEvents = (int)treeReader->GetEntries();
  writeCsv(nEvents, string(argv[1]), "C0");

  // open output file
  HistoOutputFile->cd();

  // ---------------------------------No cuts--------------------------------------------

  noCutsDirectory->cd();
  cout << "No cuts" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr,maxMjjJetPairs, noFilter);
  cout << "No cuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "noCuts");

  // ---------------------------------vbf Cuts--------------------------------------------

  vbfCutsDirectory->cd();
  cout << "vbfCuts" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr,maxMjjJetPairs, vbfCut);
  cout << "vbfCuts done." << endl;

  writeCsv(nEvents, string(argv[1]), "vbfCuts");

  // ---------------------------------photon Cuts mpp-----------------------------------------

  photonCuts_mppDirectory->cd();
  cout << "photonCuts_mpp" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr,maxMjjJetPairs, photonCuts_mpp);
  cout << "photonCuts_mpp done." << endl;

  writeCsv(nEvents, string(argv[1]), "photonCuts_mpp");

  // ---------------------------------photon Cuts pt-----------------------------------------

  photonCuts_ptDirectory->cd();
  cout << "photonCuts_pt" << endl;
  nEvents = ptEtaPhiMjjMt(treeReader, branchDict, cutsArr,maxMjjJetPairs, photonCuts_pt);
  cout << "photonCuts_pt done." << endl;

  writeCsv(nEvents, string(argv[1]), "photonCuts_pt");



  // ------------------------------------------------------------------------------------

  // close output file
  cout << "closing output file" << endl;
  HistoOutputFile->Close();

  cout << "DONE." << endl;
}
