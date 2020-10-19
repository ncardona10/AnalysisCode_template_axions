/*                             __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Counts the number of leptons in different PT ranges.                        
*/

#ifndef LEPTONCOUNTER_H
#define LEPTONCOUNTER_H

#include "./ROOTFunctions.h"
#include "./DelphesFunctions.h"
#include "../plots/MyHistograms.h"
#include "./Physics.h"
#include "./Overlaps.h"
#include <string>
#include <map>
#include <vector>
#include <set>

void fillHisto(TH1 *histo, float value)
{
  if (value > 0)
  {
    histo->Fill(value);
  }
}

bool inSet(int val, set<int> theSet)
{
  return theSet.count(val) > 0;
}

int ptEtaPhiMjjMt(ExRootTreeReader *treeReader,
                  map<string, TClonesArray *> branchDict,
                  vector<int> &cutsArr,
                  vector<pair<Jet *, Jet *>> &jetPairs,
                  bool (*filter)(ExRootTreeReader *,
                                 map<string, TClonesArray *>,
                                 int,
                                 vector<int> &,
                                 vector<pair<Jet *, Jet *>> &))
{

  int numEvents = 0;

  vector<string> variables = {"pt", "eta", "phi", "Mt","leadingPT", "dphi_MET_", "Mt_leading"};
  vector<string> particleTypes = {"electron",
                                  "muon",
                                  "tau",
                                  "jet",
                                  "photon"};

  // create histograms
  map<string, TH1 *> histos;
  for (int i = 0; (unsigned)i < variables.size(); i++)
  {
    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      int bins = 100;
      float x_min = 0.0;
      float x_max = 15.0;

      if (variables[i].compare("pt") == 0 || variables[i].compare("leadingPT") == 0)
      {
        x_max = 2000;
        bins = 1000;
      }
      else if (variables[i].compare("phi") == 0)
      {
        x_max = TMath::Pi();
        x_min = -TMath::Pi();
      }
      else if (variables[i].compare("dphi_MET_") == 0)
      {
        x_max = TMath::Pi()*2;
      }
      else if (variables[i].compare("Mt") == 0 || variables[i].compare("Mt_leading") == 0)
      {
        x_max = 1000;
      }
      else
      {
        x_min = -5;
        x_max = 5;
      }
      histos[variables[i] + particleTypes[j]] = blankHistogram(particleTypes[j] + " " + variables[i],
                                                               variables[i] + particleTypes[j],
                                                               bins, x_min, x_max); // check the histogram limits & bins
    }
  }
  
  histos["mass"] = blankHistogram("Mjj", "Mjj", 1000, 0, 10000);
  histos["DEta"] = blankHistogram("DEta", "DEta", 100, 0, 10);
  histos["MET"] = blankHistogram("MET", "MET", 100, 0, 1000);
  histos["ST"] = blankHistogram("ST", "ST", 1000, 0, 5000);

  Long64_t numberOfEntries = treeReader->GetEntries();

  for (Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    treeReader->ReadEntry(entry);

    if (filter(treeReader, branchDict, entry, cutsArr, jetPairs))
    {
      numEvents += 1;

      // MET stuff
      MissingET *METPointer = (MissingET *)branchDict["MissingET"]->At(0);
      double MET = METPointer->MET;

      // taus
      float maxPt_tau = -1.0;
      Jet *leadingTau;
      for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
      {
        Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
        if (jet->TauTag == 1)
        {
          // its a tau!
          if ((abs(jet->Eta) < 2.4) && (jet->PT > 20.0))
          {
            histos["pttau"]->Fill(jet->PT);
            histos["etatau"]->Fill(jet->Eta);
            histos["phitau"]->Fill(normalizedDphi(jet->Phi));

            double mtval = mt(jet->PT, MET, jet->Phi - METPointer->Phi);
            histos["Mttau"]->Fill(mtval);

            //leading pt
            if(jet->PT > maxPt_tau)
            {
              maxPt_tau = jet->PT;
              leadingTau = jet;
            }
          }
        }
      }
      if (maxPt_tau != -1.0)
      {
        double mtval_leading = mt(leadingTau->PT, MET, leadingTau->Phi - METPointer->Phi);
        histos["Mt_leadingtau"]->Fill(mtval_leading);
        histos["leadingPTtau"]->Fill(maxPt_tau);
        float dphi = abs(normalizedDphi(leadingTau->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_tau"]->Fill(dphi);
      }


      // electrons
      float maxPt_electrons = -1.0;
      Electron *leadingElectron;
      for (int leaf = 0; leaf < branchDict["Electron"]->GetEntries(); leaf++)
      {
        Electron *electron = (Electron *)branchDict["Electron"]->At(leaf);

        if ((abs(electron->Eta) < 2.4) && (electron->PT > 8.0))
        {
          histos["ptelectron"]->Fill(electron->PT);
          histos["etaelectron"]->Fill(electron->Eta);
          histos["phielectron"]->Fill(normalizedDphi(electron->Phi));

          double mtval = mt(electron->PT, MET, electron->Phi - METPointer->Phi);
          histos["Mtelectron"]->Fill(mtval);

          //leading pt
          if(electron->PT > maxPt_electrons)
          {
            maxPt_electrons = electron->PT;
            leadingElectron = electron;
          }
        }
      }
      if (maxPt_electrons != -1.0)
      {
        double mtval_leading = mt(leadingElectron->PT, MET, leadingElectron->Phi - METPointer->Phi);
        histos["Mt_leadingelectron"]->Fill(mtval_leading);
        histos["leadingPTelectron"]->Fill(maxPt_electrons);
        float dphi = abs(normalizedDphi(leadingElectron->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_electron"]->Fill(dphi);
      }

      // muons
      float maxPt_muons = -1.0;
      Muon *leadingMuon;
      for (int leaf = 0; leaf < branchDict["Muon"]->GetEntries(); leaf++)
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(leaf);
        if ((abs(muon->Eta) < 2.4) && (muon->PT > 5.0))
        {
          histos["ptmuon"]->Fill(muon->PT);
          histos["etamuon"]->Fill(muon->Eta);
          histos["phimuon"]->Fill(normalizedDphi(muon->Phi));

          double mtval = mt(muon->PT, MET, muon->Phi - METPointer->Phi);
          histos["Mtmuon"]->Fill(mtval);

          //leading pt
          if(muon->PT > maxPt_muons)
          {
            maxPt_muons = muon->PT;
            leadingMuon = muon;
          }
        }
      }
      if (maxPt_muons != -1.0)
      {
        double mtval_leading = mt(leadingMuon->PT, MET, leadingMuon->Phi - METPointer->Phi);
        histos["Mt_leadingmuon"]->Fill(mtval_leading);
        histos["leadingPTmuon"]->Fill(maxPt_muons);
        float dphi = abs(normalizedDphi(leadingMuon->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_muon"]->Fill(dphi);
      }

      //jets: jets minimum cuts are applied in jets cuts
      float maxPt_jet = -1.0;
      Jet *leadingJet;
      for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
      {
        Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
        if (!jet->TauTag)
        {
          histos["ptjet"]->Fill(jet->PT);
          histos["etajet"]->Fill(jet->Eta);
          histos["phijet"]->Fill(normalizedDphi(jet->Phi));

          //Doesnt make sense but needed to keep code symmetry
          double mtval = mt(jet->PT, MET, jet->Phi - METPointer->Phi);
          histos["Mtjet"]->Fill(mtval);

          //leading pt
          if(jet->PT > maxPt_jet)
          {
            maxPt_jet = jet->PT;
            leadingJet = jet;
          }
        }
      }
      if (maxPt_jet != -1.0)
      {
        double mtval_leading = mt(leadingJet->PT, MET, leadingJet->Phi - METPointer->Phi);
        histos["Mt_leadingjet"]->Fill(mtval_leading);
        histos["leadingPTjet"]->Fill(maxPt_jet);
        float dphi = abs(normalizedDphi(leadingJet->Phi)-normalizedDphi(METPointer->Phi));
        histos["dphi_MET_jet"]->Fill(dphi);        
      }

      //Photons
      Photon *photon;
      for (int leaf = 0; leaf < branchDict["Photon"]->GetEntries(); leaf++)
      {
        photon = (Photon *)branchDict["Photon"]->At(leaf);
        histos["ptphoton"]->Fill(photon->PT);
        histos["etaphoton"]->Fill(photon->Eta);
        histos["phiphoton"]->Fill(normalizedDphi(photon->Phi));

        // Doesnt make sense but needed to keep code symmetry
        // double mtval = mt(photon->PT, MET, photon->Phi - METPointer->Phi);
        histos["Mtphoton"]->Fill(1.0);

        histos["leadingPTphoton"]->Fill(1.0);
        histos["dphi_MET_photon"]->Fill(1.0);
        histos["Mt_leadingphoton"]->Fill(1.0);
 
      }

      // Mjj and deta
      // min2JetsNotTau guarantees not null pair
      if (min2JetsNotTau(treeReader, branchDict, entry))
      {
        double mass = mjj(jetPairs[entry].first, jetPairs[entry].second);
        histos["mass"]->Fill(mass);

        float deta = deltaEta(jetPairs[entry].first, jetPairs[entry].second);
        histos["DEta"]->Fill(deta);
      }

      //MET
      histos["MET"]->Fill(MET);
      float ST = MET + leadingMuon->PT + leadingElectron->PT + leadingTau->PT;
      histos["ST"]->Fill(ST);
    }
  }

  for (int i = 0; (unsigned)i < variables.size(); i++)
  {
    for (int j = 0; (unsigned)j < particleTypes.size(); j++)
    {
      histos[variables[i] + particleTypes[j]]->Write();
    }
  }

  histos["DEta"]->Write();
  histos["mass"]->Write();
  histos["MET"]->Write();
  histos["ST"]->Write();

  return numEvents;
}

#endif
