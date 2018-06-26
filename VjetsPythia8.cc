//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: VjetsPythia8_main
//
// Purpose:  This is the main production code of W/Z+jets events using the pythia 8
//           generator. Jets are reconstructed with FastJet. Systematic variations
//           are applied. Root is used for the output. The goal is to use this to
//           produce datasets for which all relevant quantities are stored in a root
//           Ntuple and in validation histograms.
//
//    Note: The main() method is used to generate events. It should not be modified
//          by users. General flags setting, histogram definitions and object
//          declarations should be made in the MyAnalysis::init(). The analysis to be
//          run in the loop must be written in void MyAnalysis::analyze(Event& event).
//          Note that the analysis here is only to calculate the physics quantities
//          to store in the Ntuple. The histogram/ntuple booking is done in the init()
//          function. They are filled in analyze(), and write to file in finish().
//          Event statistics must also be output from void MyAnalysis::finish().
//          Finally, the Pythia settings needed for a specific production are set in
//          VjetsPythia8.cmnd.

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


  // Packages initialization
  // -----------------------

     // FastJet Initialization
     // ......................

  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, Rparam,
                                      recombScheme, strategy);

     // Fastjet input
     // .............

  std::vector <fastjet::PseudoJet> fjInputs;



  // Initialize counters and other global variables
  // ----------------------------------------------

     // Note: These must be define in VjetsPythia8.h to be in-scope for all functions


     // Number of events
     // ................

  nEvt = 0;


  // Book Histograms
  // ---------------
  int nbins = 100;
  lepton_px = new TH1F("Lepton", "Lepton X-Momentum", nbins, 0.0, 100);


  // Book Trees
  // ------------
  lepton_tree = new TTree("LeptonTree","Lepton Data");
  lepton_tree->Branch("px",&px);
  lepton_tree->Branch("py",&py);
  lepton_tree->Branch("pz",&pz);
  lepton_tree->Branch("pt",&pt);
  lepton_tree->Branch("v",&v);


}

//*************************************************************************************






//*************************************************************************************
// Analysis Code
//*************************************************************************************
void MyAnalysis::analyze(Event& event)
{

  // Define an Analysis Utilities class object to be able to access the functions define there
  // -----------------------------------------------------------------------------------------

  ANA_utils myUtils;
  TruthPart mytest;

  mytest.Set(0.01, 1.02, 2.03, 3.04, 4.05, 5.06, 6.07, 7.08, 1.1, 12, 3, 33, -11);

  // Find kinematic information about vector bosons and leptons
  // ----------------------------------------------------------

     // Born electrons
     // ..............

  p_Lepton_Born = &Lepton_Born;  // Set the pointer to the vector of TLorentzVector containing the kinematic of each electrons
  TLorentzVector *p_v = &v;

  // Fill ntuples and histograms
  // ------------
  for (size_t i = 0; i < p_Lepton_Born->size(); i++) {
    px = (Lepton_Born[i]).Px();
    py = (Lepton_Born[i]).Py();
    pz = (Lepton_Born[i]).Pz();
    pt = (Lepton_Born[i]).Pt();

    p_v->SetPxPyPzE(px, py, pz, pt);

    lepton_px->Fill(px);
    lepton_tree->Fill();
  }

  Lepton_Born.clear();            // Clear the vector from previous event information

  p_Electron_Born = &Electron_Born;
  Electron_Born.clear();
  v.Clear();

  myUtils.Born_Welectron(event, p_Lepton_Born);
  //  myUtils.Born_Welectron2(event, Electron_Born);
  //  myUtils.Born_Welectron2(event);
  //myUtils.Born_Welectron2(Electron_Born);

  //std::cout << "Pt = " << (Lepton_Born[0]).Pt() << ", Px = " << (Lepton_Born[0]).Px() << ", Py = " << (Lepton_Born[0]).Py() << ", Pz = " << (Lepton_Born[0]).Pz() << std::endl;
  //std::cout << "Pt = " << (Electron_Born[0]).Pt() << ", Px = " << (Electron_Born[0]).Px() << ", Py = " << (Electron_Born[0]).Py() << ", Pz = " << (Electron_Born[0]).Pz() << std::endl;
  //std::cout << "Pt = " << mytest.Pt() << ", Px = " << mytest.Px() << ", Py = " << mytest.Py() << ", Pz = " << mytest.Pz() << std::endl;



}

//*************************************************************************************






//*************************************************************************************
// Finishing Code
//*************************************************************************************
void MyAnalysis::finish()
{

  // Normalize histograms
  // --------------------
  Double_t norm = lepton_px->GetEntries();
  lepton_px->Scale(1/norm);

  // Print histograms
  // ----------------
  lepton_px->Draw();
  lepton_tree->Write();
  //lepton_tree->StartViewer();
  lepton_tree->Print();

  // Close Ntuple files
  // ------------------
  //delete lepton_file;


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
//  lepton_file = new TFile("LeptonFile.root","CREATE","Lepton Data");
  TFile *lepton_file = new TFile("LeptonFile.root","RECREATE");

  // Analysis initialization
  // -----------------------

    // Declare user analysis class and initialize it
    // .............................................

  MyAnalysis myAnalysis;
  myAnalysis.init();


    // Read in number of event
    // .......................

  int nEvent = pythia.mode("Main:numberOfEvents");



//===========================================================================
// Loop to Generate Eevents
//===========================================================================


  // Begin event loop. Generate event. Skip if error
  // -----------------------------------------------

  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
    {
      if (!pythia.next()) continue;



   // User Analysis of current event
   // ------------------------------

      myAnalysis.analyze(pythia.event);



    }   // End of event loop.


//===========================================================================
// Control Output and Run Information
//===========================================================================


  // Pythia Statistics display
  // -------------------------

  pythia.stat();

  // User finishing
  // --------------

  myAnalysis.finish();
  lepton_file->Write();
  lepton_file->Close();

  // Done
  // ----

  return 0;
}

//*************************************************************************************
