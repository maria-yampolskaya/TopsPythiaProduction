//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: VjetsPythia8_main
//
// Purpose:  This is the main code for generating events with pythia 8, as well as
//           storing relevant output variables in Root Trees and histograms. Jets are
//           reconstructed with FastJet. The goal is to produce datasets with all the
//           relevant variables for doing theoretical predictions and sensitivity
//           studies.
//
//    Note: The main() method is used to generate events. It should not be modified
//          by users. General flags setting, histogram definitions and object
//          declarations should be made in the MyAnalysis::init(). The analysis to be
//          run in the loop must be written in void MyAnalysis::analyze(Event& event).
//          Note that the analysis here is only to calculate the physics quantities
//          to store in the Ntuple. Finally, the Pythia settings needed for a specific
//          production are set in VjetsPythia8.cmnd.
//
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


//*************************************************************************************
// Include headers
//*************************************************************************************

#include "VjetsPythia8.h"


using namespace Pythia8;


//*************************************************************************************






//*************************************************************************************
// Initialization Code
//*************************************************************************************
void MyAnalysis::init()
{


  // Initialize counters and other global variables
  // ----------------------------------------------

  // Note: These must be define in VjetsPythia8.h to be in-scope for all functions


  // Debug flag
  // ..........

  debug = false;


  // Number of events
  // ................

  nEvt = 0;



  // Book Histograms
  // ---------------

  //int nbins = 100;
  //lepton_px = new TH1F("Lepton", "Lepton X-Momentum", nbins, 0.0, 100);



  // Book Ntuples
  // ------------

  tree->Branch("top_pt",&top_pt);
  tree->Branch("top_eta",&top_eta);
  tree->Branch("top_phi",&top_phi);
  tree->Branch("top_E",&top_E);

  tree->Branch("neutrino_pt",&neutrino_pt);
  tree->Branch("neutrino_eta",&neutrino_eta);
  tree->Branch("neutrino_phi",&neutrino_phi);
  tree->Branch("neutrino_E",&neutrino_E);

  tree->Branch("muon_pt",&muon_pt);
  tree->Branch("muon_eta",&muon_eta);
  tree->Branch("muon_phi",&muon_phi);
  tree->Branch("muon_E",&muon_E);
  tree->Branch("muon_charge",&muon_charge);

  tree->Branch("electron_pt",&electron_pt);
  tree->Branch("electron_eta",&electron_eta);
  tree->Branch("electron_phi",&electron_phi);
  tree->Branch("electron_E",&electron_E);
  tree->Branch("electron_charge",&electron_charge);

  tree->Branch("boson_pt",&boson_pt);
  tree->Branch("boson_eta",&boson_eta);
  tree->Branch("boson_phi",&boson_phi);
  tree->Branch("boson_E",&boson_E);
  tree->Branch("boson_ID",&boson_ID);

  tree->Branch("lightjet_pt",&lightjet_pt);
  tree->Branch("lightjet_eta",&lightjet_eta);
  tree->Branch("lightjet_phi",&lightjet_phi);
  tree->Branch("lightjet_E",&lightjet_E);

  tree->Branch("bjet_pt",&bjet_pt);
  tree->Branch("bjet_eta",&bjet_eta);
  tree->Branch("bjet_phi",&bjet_phi);
  tree->Branch("bjet_E",&bjet_E);

  tree->Branch("lightpartonjet_pt",&lightpartonjet_pt);
  tree->Branch("lightpartonjet_eta",&lightpartonjet_eta);
  tree->Branch("lightpartonjet_phi",&lightpartonjet_phi);
  tree->Branch("lightpartonjet_E",&lightpartonjet_E);

  tree->Branch("bpartonjet_pt",&bpartonjet_pt);
  tree->Branch("bpartonjet_eta",&bpartonjet_eta);
  tree->Branch("bpartonjet_phi",&bpartonjet_phi);
  tree->Branch("bpartonjet_E",&bpartonjet_E);

  tree->Branch("bquark_pt",&bquark_pt);
  tree->Branch("bquark_eta",&bquark_eta);
  tree->Branch("bquark_phi",&bquark_phi);
  tree->Branch("bquark_E",&bquark_E);

  tree->Branch("nTop",&nTop);
  tree->Branch("nNeutrino",&nNeutrino);
  tree->Branch("nMuon",&nMuon);
  tree->Branch("nElectron",&nElectron);
  tree->Branch("nLightjet",&nLightjet);
  tree->Branch("nBjet",&nBjet);
  tree->Branch("nLightpartonjet",&nLightpartonjet);
  tree->Branch("nBpartonjet",&nBpartonjet);
  tree->Branch("nBoson",&nBoson);
  tree->Branch("mEt",&mEt);
  tree->Branch("mEt_phi",&mEt_phi);

}

//*************************************************************************************






//*************************************************************************************
// Analysis Code
//*************************************************************************************
void MyAnalysis::analyze(Event& event, Event& partonevent)
{

  // Declare an Analysis Utilities class object
  // ------------------------------------------

  // Note: To be able to access the functions define there

  ANA_utils myUtils;


  // Fill truth particle and jets information
  // ----------------------------------------

  // Set pointers
  // ............

  // Note: In each case, we first need to set a pointer to the vector of TruthPart containing the relevant particles,
  //       and then we call the function using the pointer.

  p_Top_Coll = &Top_Coll;
  p_Vecboson_Coll = &Vecboson_Coll;
  p_LeptonBare_Coll = &LeptonBare_Coll;
  p_Neutrino_Coll = &Neutrino_Coll;
  p_TruthJets_Coll = &TruthJets_Coll;
  p_PartonJets_Coll = &PartonJets_Coll;

  // Get Tops
  // ........

  myUtils.Get_Tops(event, p_Top_Coll);

  // Get Vector bosons
  // .................

  myUtils.Get_VectorBosons(event, p_Vecboson_Coll);


  // Bare leptons
  // ............

  // Note: Need to find the indices of the vector bosons found above. Don't call the function of none are found.

  std::vector<int> vecbosindex;
  for (int vb_i = 0; vb_i < Vecboson_Coll.size(); vb_i++)
  {
    vecbosindex.push_back( (Vecboson_Coll[vb_i]).Index() );
  }

  if (vecbosindex.size() > 0) myUtils.Get_BarePromptLepton(event, vecbosindex, p_LeptonBare_Coll, p_Neutrino_Coll);


  // True Jets
  // .........

  // Note 1: A list of stable particles not to be clustered in jets must first be defined for truth jets.

  // Note 2: A different function is called for final state particle jets and for pre-hadronization parton jets


  std::vector<int> skippart;
  for (int i_part = 0; i_part < LeptonBare_Coll.size(); i_part++) skippart.push_back((LeptonBare_Coll[i_part]).Index());

  // Particle jets
  myUtils.TrueJetsReco(event, skippart, p_TruthJets_Coll);


  // Parton jets
  myUtils.PartonJetsReco(event, partonevent, p_PartonJets_Coll);

  // Fill ntuples and histograms
  // ---------------------------

  nBoson = p_Vecboson_Coll->size();
  if (nBoson != 0) {
    for (size_t i = 0; i < p_Vecboson_Coll->size(); i++){
      boson_pt.push_back((Vecboson_Coll[i]).Pt());
      boson_eta.push_back((Vecboson_Coll[i]).Eta());
      boson_phi.push_back((Vecboson_Coll[i]).Phi());
      boson_E.push_back((Vecboson_Coll[i]).E());
      boson_ID.push_back((Vecboson_Coll[i]).Pdgid());
    }
  }
  else if (nBoson == 0) {
    boson_pt.push_back(0);
    boson_eta.push_back(0);
    boson_phi.push_back(0);
    boson_E.push_back(0);
    boson_ID.push_back(0);
  }

  nLightjet = 0;
  nBjet = 0;
  if (p_TruthJets_Coll->size() != 0) {
    for (size_t i = 0; i < p_TruthJets_Coll->size(); i++){
      if ((TruthJets_Coll[i]).BQTag() || (TruthJets_Coll[i]).BHTag()) {
        bjet_pt.push_back((TruthJets_Coll[i]).Pt());
        bjet_eta.push_back((TruthJets_Coll[i]).Eta());
        bjet_phi.push_back((TruthJets_Coll[i]).Phi());
        bjet_E.push_back((TruthJets_Coll[i]).E());
        nBjet += 1;
      }
      else {
        lightjet_pt.push_back((TruthJets_Coll[i]).Pt());
        lightjet_eta.push_back((TruthJets_Coll[i]).Eta());
        lightjet_phi.push_back((TruthJets_Coll[i]).Phi());
        lightjet_E.push_back((TruthJets_Coll[i]).E());
        nLightjet += 1;
      }
    }
  }
  if (nLightjet == 0) {
    lightjet_pt.push_back(0);
    lightjet_eta.push_back(0);
    lightjet_phi.push_back(0);
    lightjet_E.push_back(0);
  }
  if (nBjet == 0) {
    bjet_pt.push_back(0);
    bjet_eta.push_back(0);
    bjet_phi.push_back(0);
    bjet_E.push_back(0);
  }

  nLightpartonjet = 0;
  nBpartonjet = 0;
  if (p_PartonJets_Coll->size() != 0) {
    for (size_t i = 0; i < p_PartonJets_Coll->size(); i++){
      if ((PartonJets_Coll[i]).BQTag() || (PartonJets_Coll[i]).BHTag()) {
        bpartonjet_pt.push_back((PartonJets_Coll[i]).Pt());
        bpartonjet_eta.push_back((PartonJets_Coll[i]).Eta());
        bpartonjet_phi.push_back((PartonJets_Coll[i]).Phi());
        bpartonjet_E.push_back((PartonJets_Coll[i]).E());
        nBpartonjet += 1;
      }
      else {
        lightpartonjet_pt.push_back((PartonJets_Coll[i]).Pt());
        lightpartonjet_eta.push_back((PartonJets_Coll[i]).Eta());
        lightpartonjet_phi.push_back((PartonJets_Coll[i]).Phi());
        lightpartonjet_E.push_back((PartonJets_Coll[i]).E());
        nLightpartonjet += 1;
      }
    }
  }
  if (nLightpartonjet == 0) {
    lightpartonjet_pt.push_back(0);
    lightpartonjet_eta.push_back(0);
    lightpartonjet_phi.push_back(0);
    lightpartonjet_E.push_back(0);
  }
  if (nBpartonjet == 0) {
    bpartonjet_pt.push_back(0);
    bpartonjet_eta.push_back(0);
    bpartonjet_phi.push_back(0);
    bpartonjet_E.push_back(0);
  }

  nElectron = 0;
  nMuon = 0;
  if (p_LeptonBare_Coll->size() != 0) {
    for (size_t i = 0; i < p_LeptonBare_Coll->size(); i++){
      if ((LeptonBare_Coll[i]).Pdgid() == 11 || (LeptonBare_Coll[i]).Pdgid() == -11) {
        electron_pt.push_back((LeptonBare_Coll[i]).Pt());
        electron_eta.push_back((LeptonBare_Coll[i]).Eta());
        electron_phi.push_back((LeptonBare_Coll[i]).Phi());
        electron_E.push_back((LeptonBare_Coll[i]).E());
        electron_charge.push_back((LeptonBare_Coll[i]).Charge());
        nElectron += 1;
      }
      if ((LeptonBare_Coll[i]).Pdgid() == 13 || (LeptonBare_Coll[i]).Pdgid() == -13) {
        muon_pt.push_back((LeptonBare_Coll[i]).Pt());
        muon_eta.push_back((LeptonBare_Coll[i]).Eta());
        muon_phi.push_back((LeptonBare_Coll[i]).Phi());
        muon_E.push_back((LeptonBare_Coll[i]).E());
        muon_charge.push_back((LeptonBare_Coll[i]).Charge());
        nMuon += 1;
      }
    }
  }
  if (nElectron == 0) {
    electron_pt.push_back(0);
    electron_eta.push_back(0);
    electron_phi.push_back(0);
    electron_E.push_back(0);
    electron_charge.push_back(0);
  }
  if (nMuon == 0) {
    muon_pt.push_back(0);
    muon_eta.push_back(0);
    muon_phi.push_back(0);
    muon_E.push_back(0);
    muon_charge.push_back(0);
  }


  nTop = p_Top_Coll->size();
  if (nTop != 0) {
    for (size_t i = 0; i < p_Top_Coll->size(); i++){
      top_pt.push_back((Top_Coll[i]).Pt());
      top_eta.push_back((Top_Coll[i]).Eta());
      top_phi.push_back((Top_Coll[i]).Phi());
      top_E.push_back((Top_Coll[i]).E());
    }
  }
  else if (nTop == 0) {
    top_pt.push_back(0);
    top_eta.push_back(0);
    top_phi.push_back(0);
    top_E.push_back(0);
  }

  nNeutrino = p_Neutrino_Coll->size();
  if (nNeutrino != 0) {
    for (size_t i = 0; i < p_Neutrino_Coll->size(); i++){
      neutrino_pt.push_back((Neutrino_Coll[i]).Pt());
      neutrino_eta.push_back((Neutrino_Coll[i]).Eta());
      neutrino_phi.push_back((Neutrino_Coll[i]).Phi());
      neutrino_E.push_back((Neutrino_Coll[i]).E());
    }
  }
  else if (nNeutrino == 0) {
    neutrino_pt.push_back(0);
    neutrino_eta.push_back(0);
    neutrino_phi.push_back(0);
    neutrino_E.push_back(0);
  }

  float sqsum;
  float sum_y;
  float sum_x;

  for (auto i : Neutrino_Coll) {
    sqsum += pow(i.Px(), 2.0) + pow(i.Py(), 2.0);
    sum_y += i.Py();
    sum_x += i.Px();
  }

  mEt = sqrt(sqsum);
  mEt_phi = atan(sum_y/sum_x);

  tree->Fill();



  // Clear event-based vectors
  // -------------------------

  top_pt.clear();
  top_eta.clear();
  top_phi.clear();
  top_E.clear();

  neutrino_pt.clear();
  neutrino_eta.clear();
  neutrino_phi.clear();
  neutrino_E.clear();

  muon_pt.clear();
  muon_eta.clear();
  muon_phi.clear();
  muon_E.clear();
  muon_charge.clear();

  electron_pt.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_E.clear();
  electron_charge.clear();

  bjet_pt.clear();
  bjet_eta.clear();
  bjet_phi.clear();
  bjet_E.clear();

  lightjet_pt.clear();
  lightjet_eta.clear();
  lightjet_phi.clear();
  lightjet_E.clear();

  bpartonjet_pt.clear();
  bpartonjet_eta.clear();
  bpartonjet_phi.clear();
  bpartonjet_E.clear();

  lightpartonjet_pt.clear();
  lightpartonjet_eta.clear();
  lightpartonjet_phi.clear();
  lightpartonjet_E.clear();

  bquark_pt.clear();
  bquark_eta.clear();
  bquark_phi.clear();
  bquark_E.clear();

  boson_pt.clear();
  boson_eta.clear();
  boson_phi.clear();
  boson_E.clear();
  boson_ID.clear();

  Top_Coll.clear();
  Vecboson_Coll.clear();
  LeptonBare_Coll.clear();
  Neutrino_Coll.clear();
  TruthJets_Coll.clear();
  PartonJets_Coll.clear();
  skippart.clear();


}

//*************************************************************************************






//*************************************************************************************
// Finishing Code
//*************************************************************************************
void MyAnalysis::finish()
{

  // Normalize histograms
  // --------------------

  //Double_t norm = lepton_px->GetEntries();
  //lepton_px->Scale(1/norm);


  // Print histograms
  // ----------------


  // Write trees

  tree->Write();

  // Final Print Statements if needed
  // --------------------------------


}

//*************************************************************************************






//*************************************************************************************
// Main Code
//*************************************************************************************
int main(int argc, char* argv[])
{


  //===========================================================================
  // Initialization Phase
  //===========================================================================


  // Safety checks to make sure the code is run properly
  // ---------------------------------------------------

  // Check that correct number of command-line arguments
  // ...................................................

  if (argc != 2)
  {
    cerr << " Unexpected number of command-line arguments. \n"
    << " You are expected to provide a file name and nothing else. \n"
    << " Program stopped! " << endl;
    return 1;
  }


  // Check that the provided file name corresponds to an existing file
  // .................................................................

  ifstream is(argv[1]);
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
    << " Program stopped! " << endl;
    return 1;
  }


  // Confirm that external file will be used for settings
  // ....................................................

  cout << " PYTHIA settings will be read from file " << argv[1] << endl;



  // Pythia 8 initialization
  // -----------------------

  // Declare a pythia object
  // .......................

  Pythia pythia;


  // Specification of Pythia settings
  // ................................

  // Note: These settings are specified in the input file

  pythia.readFile(argv[1]);


  // Initialize Pythia
  // .................

  pythia.init();



  // Prepare input file
  // ------------------

  TFile *myfile = TFile::Open("outfile.root","recreate");



  // Analysis initialization
  // -----------------------

  // Declare user analysis class and initialize it
  // .............................................

  MyAnalysis myAnalysis;
  myAnalysis.init();


  // Read in number of event
  // .......................

  int nEvent = pythia.mode("Main:numberOfEvents");


  // Some global variables
  // ---------------------

  int nListEvts = 2;


  // Declare Event Variables
  // -----------------------

  // An event record for parton level particles
  // ..........................................

  // Note: Partons will be taken at the end of the parton evolution
  //       i.e. just before the hadronization.

  Event partonLevelEvent;
  partonLevelEvent.init("Parton Level event record", &pythia.particleData);



  //===========================================================================
  // Loop to Generate Events
  //===========================================================================


  // Begin event loop. Generate event. Skip if error
  // -----------------------------------------------

  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythia.next()) continue;


    // Get the parton-level event
    // --------------------------

    // Declare an Analysis Utilities class object
    // ------------------------------------------

    // Note: To be able to access the functions define there
    
    ANA_utils myUtilsMain;



    myUtilsMain.getPartonLevelEvent(pythia.event, partonLevelEvent);


    // Display some event info
    // -----------------------

    // List first few events
    // .....................

    if (iEvent < nListEvts)
    {
      pythia.event.list();
      partonLevelEvent.list();
    }


    // User Analysis of current event
    // ------------------------------

    myAnalysis.analyze(pythia.event, partonLevelEvent);



  }   // End of event loop.


  //===========================================================================
  // Control Output and Run Information
  //===========================================================================


  // Pythia Statistics display
  // -------------------------

  pythia.stat();


  // Write root output
  // -----------------

  //T->Write();


  // User finishing
  // --------------

  myAnalysis.finish();


  // Write info to file
  // ------------------

  myfile->Write();
  myfile->Close();


  return 0;
}

//*************************************************************************************
