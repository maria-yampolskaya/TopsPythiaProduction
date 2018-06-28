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

  int nbins = 100;
  lepton_px = new TH1F("Lepton", "Lepton X-Momentum", nbins, 0.0, 100);


  
  // Book Ntuples
  // ------------

  lepton_ntuple = new TNtuple("Lepton", "Lepton Data", "px:py:pz:E");

}

//************************************************************************************* 






//*************************************************************************************
// Analysis Code
//*************************************************************************************
void MyAnalysis::analyze(Event& event)
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
  p_TruthJetsColl = &TruthJetsColl;   


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

        // Note: A list of stable particles not to be clustered in jets must first be defined

  std::vector<int> skippart;
  for (int i_part = 0; i_part < LeptonBare_Coll.size(); i_part++) skippart.push_back((LeptonBare_Coll[i_part]).Index());

  myUtils.TrueJetsReco(event, skippart, p_TruthJetsColl);



  // Fill ntuples and histograms
  // ---------------------------
  
  for (size_t i = 0; i < p_LeptonBare_Coll->size(); i++)
    {
      lepton_ntuple->Fill((LeptonBare_Coll[i]).Px(),(LeptonBare_Coll[i]).Py(),(LeptonBare_Coll[i]).Pz(),(LeptonBare_Coll[i]).E());
      lepton_px->Fill((LeptonBare_Coll[i]).Px());
    }

  
  // Clear event-based vectors
  // -------------------------

  Top_Coll.clear();              
  Vecboson_Coll.clear();         
  LeptonBare_Coll.clear();           
  Neutrino_Coll.clear();           
  TruthJetsColl.clear();
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

  lepton_px->Draw();

  
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


  // Declare Event Variables
  // -----------------------



//=========================================================================== 
// Loop to Generate Events
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
