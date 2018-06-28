//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: VjetsPythia8.h
// 
// Purpose:  Header file for the Pythia8 production and analysis
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// Dependencies (#includes)
// ========================

// Pythia include
// --------------

#include "Pythia8/Pythia.h"


// Packages include
// ----------------

// FastJet
// .......

  // Note: The FastJet3.h header enables automatic initialisation of fastjet::PseudoJet 
  //       objects from Pythia8 Particle and Vec4 objects, as well as advanced features
  //       such as access to (a copy of) the original Pythia 8 Particle directly from
  //       the PseudoJet, and fastjet selectors that make use of the Particle properties.
  //       See the extensive comments in the header file for further details and examples.

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "Pythia8Plugins/FastJet3.h"


// Root
// ....

#include "TH1.h"
#include "TVirtualPad.h"     // Interactive graphics.
#include "TApplication.h"
#include "TFile.h"           // Saving file.
#include "TTree.h"
#include "TNtuple.h"
#include "TLorentzVector.h" 


// My include
// ----------

#include "ANA_utils.h"
#include "TruthPart.h"
#include "TruthJets.h"


// C/C++ include
// -------------

#include <vector>
#include <string>
#include <iostream>

using namespace Pythia8;




// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
// class MyAnakysis
//
//   Note: Main analysis class definition
//  
// -----------------------------------------------------------------------------------
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

class MyAnalysis
{

 public:


  // Constructor and destructor
  // --------------------------

  MyAnalysis() {}


  // Declare functions
  // -----------------

     // Note: Utility functions are defined in ANA_utils.cc, while the main generator
     //       functions and the analysis are defined in VjetsPythia.cc.


  // Initialization actions
  // ......................
  
  void init();

  
  // Analysis of each new event
  // ..........................
  
  void analyze(Event& event);

  
  // Show final results
  // ..................
  
  void finish();



  
  // Declare variables and objects that span init - analyze - finish
  // ---------------------------------------------------------------
  
private:

  
  // Global variables
  // ................
  
  int nEvt;
  int nEventAccept;
  int vetoCount[4];

  bool firstEvent;

  bool debug;

  // Lepton and jet level
  // ....................

  string level;

  
  // Declare Objects and Vectors of Objects
  // --------------------------------------

  std::vector<TruthPart> Top_Coll;
  std::vector<TruthPart>* p_Top_Coll;

  std::vector<TruthPart> Vecboson_Coll;
  std::vector<TruthPart>* p_Vecboson_Coll;

  std::vector<TruthPart> LeptonBare_Coll;
  std::vector<TruthPart>* p_LeptonBare_Coll;

  std::vector<TruthPart> Neutrino_Coll;
  std::vector<TruthPart>* p_Neutrino_Coll;

  std::vector<TruthJets> TruthJetsColl;
  std::vector<TruthJets>* p_TruthJetsColl;

  
  // Histograms
  // ..........

  TH1* lepton_px;



  // Others
  // ......
  
  TNtuple *lepton_ntuple;
  TTree *myTree = new TTree("T","ev1 Tree");


};




