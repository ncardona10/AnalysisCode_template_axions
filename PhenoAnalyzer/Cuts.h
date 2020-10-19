/*
                               __
                           .--()°'.'
 Author: Nathalia Cardona '|, . ,'
                           !_-(_\
 Cuts for single particles
*/

#include "./ROOTFunctions.h"
#include "./DelphesFunctions.h"
#include "./Physics.h"
#include "./Overlaps.h"
#include "./LeptonCounter.h"

/*
  Apply no filter to the events. Use to the get the no cuts distributions
*/

bool noFilter(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr,
              vector<pair<Jet *, Jet *>> &jetPairs)
{
  return true;
}

// ----------------------------------------------------------------------------------------------

/*
  Verify met conditions over the events
  Conditions:
    met > 150.0 GeV
*/

bool met(ExRootTreeReader *treeReader,
         map<string, TClonesArray *> branchDict,
         int entry,
         vector<int> &cutsArr,
         vector<pair<Jet *, Jet *>> &jetPairs)
{
  // met>150

  if (cutsArr[entry])
  {
    treeReader->ReadEntry(entry);

    bool metBool = met(treeReader, branchDict, entry) > 150;

    cutsArr[entry] = metBool;
    return metBool;
  }

  return false;
}

// -----------------------------------------------------------------------------------------------

/*
  Verify bjets conditions over the jets collection
  Conditions:
    N(bjets) = 0
    Reconstruction cuts for bjets
      pT(bjets) > 30.0 GeV
      |eta(bjets)| < 2.5
*/

bool bjets(ExRootTreeReader *treeReader,
           map<string, TClonesArray *> branchDict,
           int entry,
           vector<int> &cutsArr,
           vector<pair<Jet *, Jet *>> &jetPairs)
{
  // # bjets = 0

  treeReader->ReadEntry(entry);

  bool bJetsBool = true;

  if (cutsArr[entry])
  {

    for (int leaf = 0; leaf < branchDict["Jet"]->GetEntries(); leaf++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(leaf);
      if (jet->BTag == 1 && jet->PT > 30.0 && abs(jet->Eta) < 2.5)
      {
        bJetsBool = false;
        break;
      }
    }
  }
  else
  {
    bJetsBool = false;
  }
  cutsArr[entry] = bJetsBool;
  return bJetsBool;
}

// -----------------------------------------------------------------------------------------------

/*
  Verify ISR jet conditions over the leading jet
  Conditions:
    pT(leading) > 100.0 GeV
    |eta(leading)| < 5.0
*/

bool jetsCuts(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr,
              vector<pair<Jet *, Jet *>> &jetPairs)
{
  treeReader->ReadEntry(entry);
  if (cutsArr[entry])
  {
    bool ans = false;
    float max_pt = -1;

    Jet *leadingJet;
    for (int i = 0; i < branchDict["Jet"]->GetEntries(); i++)
    {
      Jet *jet = (Jet *)branchDict["Jet"]->At(i);
      if (jet->PT > max_pt)
      {
        max_pt = jet->PT;
        leadingJet = jet;
      }
    }

    // Check that at least one jet exists
    if (max_pt != -1)
    {
      ans = (leadingJet->PT > 100) && abs(leadingJet->Eta) < 5.0;
    }

    cutsArr[entry] = ans;

    return ans;
  }
  return false;
}

// -----------------------------------------------------------------------------------------------

/*
  Find the leading and subleading photons
  Return: Pair of photons (leading, subleading)
  Note: 
    To obtain the leading from the pair use
      pairOfPhotons.first
    To obtain the subleading from the pair use
      pairOfPhotons.second
  
*/

pair<Photon *, Photon *> getLeadAndSubLeadPhotons(ExRootTreeReader *treeReader,
                                                  map<string, TClonesArray *> branchDict,
                                                  int entry)
{
  treeReader->ReadEntry(entry);

  // Create the leading photon
  Photon *leadingPhoton;
  // Create the subleading photon
  Photon *subleadingPhoton;

  float leadingPT = -1.0;
  float subleadingPT = -1.0;

  // Find the leading jet
  Photon *photon;
  for (int i = 0; i < branchDict["Photon"]->GetEntries(); i++)
  {
    photon = (Photon *)branchDict["Photon"]->At(i);

    if (photon->PT > leadingPT)
    {
      subleadingPhoton = leadingPhoton;
      subleadingPT = leadingPT;

      leadingPhoton = photon;
      leadingPT = photon->PT;
    }

    else if (photon->PT > subleadingPT && photon->PT <= leadingPT)
    {
      subleadingPhoton = photon;
      subleadingPT = photon->PT;
    }
  }

  pair<Photon *, Photon *> pairOfPhotons(leadingPhoton, subleadingPhoton);

  return pairOfPhotons;
}

// -----------------------------------------------------------------------------------------------

/*
  Verify if exist at least one pair of jets that satisfies the VBF conditions and the reconstruction conditions
  Jets contidions:
    pT > 20.0 GeV
    |eta| < 5.0
  VBF Conditions:
    m(j,j) > 750.0 GeV
    eta(j1) * eta(j2) < 0.0
    | DEta(j1, j2)| > 3.6
*/

bool vbfCut(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry,
            vector<int> &cutsArr,
            vector<pair<Jet *, Jet *>> &jetPairs)

{
  bool ans = 0;

  treeReader->ReadEntry(entry);

  if (cutsArr[entry])
  {
    bool min2JetsBool = min2JetsNotTau(treeReader, branchDict, entry);

    if (min2JetsBool)
    {
      Jet *jet1;
      Jet *jet2;

      bool isVBF = false;

      for (int i = 0; (i < (branchDict["Jet"]->GetEntries() - 1)) && !isVBF; i++)
      {
        jet1 = (Jet *)branchDict["Jet"]->At(i);
        if ((jet1->TauTag == 0) && (jet1->PT > 30.0) && (abs(jet1->Eta) < 5.0))
        {
          for (int j = i + 1; (j < branchDict["Jet"]->GetEntries()) && !isVBF; j++)
          {
            jet2 = (Jet *)branchDict["Jet"]->At(j);
            if ((jet2->TauTag == 0) && (jet2->PT > 30.0) && (abs(jet2->Eta) < 5.0))
            {
              bool mjjBool = mjj(jet1, jet2) > 750.0;
              bool etaMultBool = (jet1->Eta * jet2->Eta) < 0;
              bool deltaEtaBool = deltaEta(jet1, jet2) > 3.6;

              isVBF = mjjBool
                      && etaMultBool
                      && deltaEtaBool;
            }
          }
        }
      }
      ans = isVBF;
    }
  }

  cutsArr[entry] = ans;
  return ans;
}

// -----------------------------------------------------------------------------------------------

/*
  Verify pT condition over the leading and subleading photons
  Conditions:
    pT(leading photon) > 300.0 GeV
    pT(subleading photon) > 10.0 GeV
*/

bool photonCuts_pt(ExRootTreeReader *treeReader,
                   map<string, TClonesArray *> branchDict,
                   int entry,
                   vector<int> &cutsArr,
                   vector<pair<Jet *, Jet *>> &jetPairs)
{
  treeReader->ReadEntry(entry);
  bool ans = false;
  if (cutsArr[entry])
  {
    bool min2PhotonsBool = minNPhotons(treeReader, branchDict, entry, 2);
    if (min2PhotonsBool)
    {
      pair<Photon *, Photon *> leadingPhotons = getLeadAndSubLeadPhotons(treeReader, branchDict, entry);

      // verify Pt conditions
      ans = (leadingPhotons.first->PT > 300.0)
          && (leadingPhotons.second->PT > 30.0) 
          && (abs(leadingPhotons.first->Eta) < 2.5) 
          && (abs(leadingPhotons.second->Eta) < 2.5) 
          ;
    }
  }
  cutsArr[entry] = ans;
  return ans;
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool photonCuts_mpp(ExRootTreeReader *treeReader,
                    map<string, TClonesArray *> branchDict,
                    int entry,
                    vector<int> &cutsArr,
                    vector<pair<Jet *, Jet *>> &jetPairs)
{
  treeReader->ReadEntry(entry);
  bool ans = false;
  if (cutsArr[entry])
  {
    bool min2PhotonsBool = minNPhotons(treeReader, branchDict, entry, 2);
    if (min2PhotonsBool)
    {
      const double mpp = 500.0;
      Photon *photon1;
      Photon *photon2;

      bool reqConditionMpp = false;

      for (int i = 0; (i < branchDict["Photon"]->GetEntries() - 1) && !reqConditionMpp; i++)
      {
        photon1 = (Photon *)branchDict["Photon"]->At(i);
        for (int j = i + 1; (j < branchDict["Photon"]->GetEntries()) && !reqConditionMpp; j++)
        {
          photon2 = (Photon *)branchDict["Photon"]->At(j);
          reqConditionMpp = mphoton_photon(photon1, photon2) > mpp;
        }
      }
      ans = reqConditionMpp;
    }
  }
  cutsArr[entry] = ans;
  return ans;
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

//Verify if exist at least one electron and tau that satisfy the minimum pt and |eta| cuts
bool leptonsCuts(ExRootTreeReader *treeReader,
                 map<string, TClonesArray *> branchDict,
                 int entry,
                 vector<int> &cutsArr,
                 vector<pair<Jet *, Jet *>> &jetPairs)
{
  treeReader->ReadEntry(entry);
  if (cutsArr[entry])
  {

    bool ansElec = false;
    bool ansTau = false;
    bool ans = false;

    int cont = 0;

    while (!ansElec && cont < branchDict["Electron"]->GetEntries())
    {
      // Get a lepton
      Electron *elec = (Electron *)branchDict["Electron"]->At(cont);

      // To create a muon use:
      //Muon *muon = (Muon *)branchDict["Muon"]->At(cont);

      if (elec->PT > 5.0 && abs(elec->Eta) < 2.4)
      {
        ansElec = true;
      }
      cont++;
    }

    //Reset counter
    cont = 0;

    while (!ansTau && cont < branchDict["Jet"]->GetEntries())
    {
      // Get a jet
      Jet *jet = (Jet *)branchDict["Jet"]->At(cont);

      // Verify if the jet is a tau
      if (jet->TauTag == 1)
      {
        if (jet->PT >= 20.0 && abs(jet->Eta) < 2.4)
        {
          ansTau = true;
        }
      }
      cont++;
    }

    ans = ansElec && ansTau;

    cutsArr[entry] = ans;
    return ans;
  }
  return false;
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

//Baseline to create the lepton channels. Conside only an specific number for every lepton
bool nParticle(ExRootTreeReader *treeReader,
               map<string, TClonesArray *> branchDict,
               int entry,
               int n_electrons,
               int n_muon,
               int n_tau,
               vector<int> &cutsArr,
               bool checkElecPT,
               bool checkMuonPT,
               bool checkTauPT)
{

  /*
    must comply with VBF cuts  & cuts
    n_electrons electron with:
      pt>8 
      abs(eta) < 2.4
    n_tau taus with:
      pt>20
      abs(eta)<2.4
    n_muon muons with:
      pt>5
      abs(eta)<2.4
  */

  treeReader->ReadEntry(entry);

  if (cutsArr[entry])
  {
    // verify electron condition
    int nElectrons = 0;
    int i = 0;
    vector<TLorentzVector> particleCharacteristics;

    while (nElectrons <= n_electrons && i < branchDict["Electron"]->GetEntries())
    {
      Electron *elec = (Electron *)branchDict["Electron"]->At(i);
      if (elec->PT >= 8 && abs(elec->Eta) < 2.4)
      {
        if (checkElecPT)
        {
          if (elec->PT <= 40)
          {
            nElectrons++;
            particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
          }
        }

        else
        {
          nElectrons++;
          particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
        }
      }
      i++;
    }

    if (nElectrons == n_electrons)
    {
      //verify number of muons
      int nMuons = 0;
      int i = 0;

      particleCharacteristics.clear();

      while (nMuons <= n_muon && i < branchDict["Muon"]->GetEntries())
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(i);
        if (muon->PT >= 5 && abs(muon->Eta) < 2.4)
        {

          if (checkMuonPT)
          {
            if (muon->PT <= 40)
            {
              nMuons++;
              particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
            }
          }
          else
          {
            nMuons++;
            particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
          }
        }
        i++;
      }

      if (nMuons == n_muon)
      {
        //verify number of taus
        int nTaus = 0;
        int i = 0;
        particleCharacteristics.clear();

        while (nTaus <= n_tau && i < branchDict["Jet"]->GetEntries())
        {
          Jet *jet = (Jet *)branchDict["Jet"]->At(i);
          if (jet->TauTag == 1)
          {
            if (checkTauPT)
            {
              if (jet->PT >= 20.0 && abs(jet->Eta) < 2.4)
              {
                nTaus++;
                particleCharacteristics.push_back(createTLorentzVector(jet->PT, jet->Eta, jet->Mass, jet->Phi));
              }
            }
          }
          i++;
        }
        return nTaus == n_tau;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

/*
-------------------------------- SINGLE PARTICLE CHANNELS ----------------------------------------
*/

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool mono_e(ExRootTreeReader *treeReader,
            map<string, TClonesArray *> branchDict,
            int entry,
            vector<int> &cutsArr,
            vector<pair<Jet *, Jet *>> &jetPairs)
{
  return nParticle(treeReader, branchDict, entry, 1, 0, 0, cutsArr, true, false, false);
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool mono_mu(ExRootTreeReader *treeReader,
             map<string, TClonesArray *> branchDict,
             int entry,
             vector<int> &cutsArr,
             vector<pair<Jet *, Jet *>> &jetPairs)
{
  return nParticle(treeReader, branchDict, entry, 0, 1, 0, cutsArr, false, true, false);
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool mono_tau(ExRootTreeReader *treeReader,
              map<string, TClonesArray *> branchDict,
              int entry,
              vector<int> &cutsArr,
              vector<pair<Jet *, Jet *>> &jetPairs)
{
  return nParticle(treeReader, branchDict, entry, 0, 0, 1, cutsArr, false, false, true);
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

//Baseline to create the lepton channels. Consider at least an specific number for every lepton
bool nParticle_atLeast(ExRootTreeReader *treeReader,
                       map<string, TClonesArray *> branchDict,
                       int entry,
                       int n_electrons,
                       int n_muon,
                       int n_tau,
                       vector<int> &cutsArr,
                       bool checkElecPT,
                       bool checkMuonPT,
                       bool checkTauPT)
{

  /*
    n_electrons electron with:
      pt>8 
      abs(eta) < 2.4
    n_tau taus with:
      pt>20
      abs(eta)<2.4
    n_muon muons with:
      pt>5
      abs(eta)<2.4
  */

  treeReader->ReadEntry(entry);

  if (cutsArr[entry])
  {
    // verify electron condition
    int nElectrons = 0;
    int i = 0;
    vector<TLorentzVector> particleCharacteristics;

    while (nElectrons <= n_electrons && i < branchDict["Electron"]->GetEntries())
    {
      Electron *elec = (Electron *)branchDict["Electron"]->At(i);
      if (elec->PT >= 8 && abs(elec->Eta) < 2.4)
      {
        if (checkElecPT)
        {
          if (elec->PT <= 40)
          {
            nElectrons++;
            particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
          }
        }

        else
        {
          nElectrons++;
          particleCharacteristics.push_back(createTLorentzVector(elec->PT, elec->Eta, 0.000510998902, elec->Phi));
        }
      }
      i++;
    }

    if (nElectrons >= n_electrons)
    {
      //verify number of muons
      int nMuons = 0;
      int i = 0;

      particleCharacteristics.clear();

      while (nMuons <= n_muon && i < branchDict["Muon"]->GetEntries())
      {
        Muon *muon = (Muon *)branchDict["Muon"]->At(i);
        if (muon->PT >= 5 && abs(muon->Eta) < 2.4)
        {

          if (checkMuonPT)
          {
            if (muon->PT <= 40)
            {
              nMuons++;
              particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
            }
          }
          else
          {
            nMuons++;
            particleCharacteristics.push_back(createTLorentzVector(muon->PT, muon->Eta, 0.1056583715, muon->Phi));
          }
        }
        i++;
      }

      if (nMuons >= n_muon)
      {
        //verify number of taus
        int nTaus = 0;
        int i = 0;
        particleCharacteristics.clear();

        while (nTaus <= n_tau && i < branchDict["Jet"]->GetEntries())
        {
          Jet *jet = (Jet *)branchDict["Jet"]->At(i);
          if (jet->TauTag == 1)
          {
            if (checkTauPT)
            {
              if (jet->PT >= 20.0 && abs(jet->Eta) < 2.4)
              {
                nTaus++;
                particleCharacteristics.push_back(createTLorentzVector(jet->PT, jet->Eta, jet->Mass, jet->Phi));
              }
            }
          }
          i++;
        }
        return nTaus >= n_tau;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

/*
---------------------------------- MULTILEPTONS CHANNEL ----------------------------------------------
*/

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool mu_mu_e_atLeast(ExRootTreeReader *treeReader,
                     map<string, TClonesArray *> branchDict,
                     int entry,
                     vector<int> &cutsArr,
                     vector<pair<Jet *, Jet *>> &jetPairs)
{
  return nParticle_atLeast(treeReader, branchDict, entry, 1, 2, 0, cutsArr, false, false, false);
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool e_e_e_atLeast(ExRootTreeReader *treeReader,
                   map<string, TClonesArray *> branchDict,
                   int entry,
                   vector<int> &cutsArr,
                   vector<pair<Jet *, Jet *>> &jetPairs)
{
  return nParticle_atLeast(treeReader, branchDict, entry, 3, 0, 0, cutsArr, false, false, false);
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool mu_mu_mu_atLeast(ExRootTreeReader *treeReader,
                      map<string, TClonesArray *> branchDict,
                      int entry,
                      vector<int> &cutsArr,
                      vector<pair<Jet *, Jet *>> &jetPairs)
{
  return nParticle_atLeast(treeReader, branchDict, entry, 0, 3, 0, cutsArr, false, false, false);
}

// -----------------------------------------------------------------------------------------------

/*
  
*/

bool e_e_mu_atLeast(ExRootTreeReader *treeReader,
                    map<string, TClonesArray *> branchDict,
                    int entry,
                    vector<int> &cutsArr,
                    vector<pair<Jet *, Jet *>> &jetPairs)
{
  return nParticle_atLeast(treeReader, branchDict, entry, 2, 1, 0, cutsArr, false, false, false);
}

// -----------------------------------------------------------------------------------------------