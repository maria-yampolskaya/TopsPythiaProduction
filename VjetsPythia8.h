//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: VjetsPythia8.h
//
// Purpose:  Header file for the Pythia8 production and analysis
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//#ifndef VjetsPythia8_h
//#define VjetsPythia8_h


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
#include "TLorentzVector.h"


// My include
// ----------

#include "ANA_utils.h"
#include "TruthPart.h"


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

     // Note: Utility functions are defined in VjetsPythia8_utils.cc, the main generator
     //       functions are defined in VjetsPythia_main.cc, and the analysis is defined
     //       in VjetsPythia_ana.cc.


  // Initialization actions
  // ......................

  void init();


  // Analysis of each new event
  // ..........................

  void analyze(Event& event);


  // Show final results
  // ..................

  void finish();


  //double compute_weight(const int NB, const double x_sect, const double inst_lum);


  // Declare variables and objects that span init - analyze - finish
  // ---------------------------------------------------------------

private:

  // Global variables
  // ................

  int nEvt;
  int nEventAccept;
  int vetoCount[4];

  bool firstEvent;


  // Lepton and jet level
  // ....................

  string level;

  // Declare Objects and Vectors of Objects
  // --------------------------------------

  std::vector<TLorentzVector> Lepton_Born;
  std::vector<TLorentzVector>* p_Lepton_Born;

  std::vector<TruthPart> Electron_Born;
  std::vector<TruthPart>* p_Electron_Born;



  // Histograms
  // ..........
  TH1* lepton_px;



  // Others
  // ......

  double dSigmaBin[5];
  //TFile *lepton_file;
  TTree *lepton_tree;
  Double_t px, py, pz, pt;
  TLorentzVector v;


};
