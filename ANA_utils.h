//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: ANA_utils.h
// 
// Purpose:  Header file for analysis
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// Dependencies (#includes)
// ------------------------

  // C++
  // ...

#include <cmath>
#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <fstream>

  // Root
  // ....

#include "TLorentzVector.h" 


// FastJet
// .......

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "Pythia8Plugins/FastJet3.h"


  // My Includes
  // ...........

#include "TruthPart.h"
#include "TruthJets.h"


#include "Pythia8/Pythia.h"

using namespace std;
using namespace Pythia8;





// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
// class ANA_utils:
//
//       Functions to be used in the analysis part of the generation process in order
//       to store the relevant information in the ntuples, ttrees, and histograms.
//  
// -----------------------------------------------------------------------------------
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

class ANA_utils
{

public :

  // Constructor and destructor
  // --------------------------

  ANA_utils(){}
  
  ~ANA_utils(){}


  
  // Declare functions
  // -----------------

     // To compute alpgen weight
     // ........................

  double compute_weight(const int NB, const double x_sect, const double inst_lum);


     // Fill truth particle and truth jet information
     // .............................................

  void Fill_TruthPart(Pythia8::Event event, int index, TruthPart* p_TruthPart);

  void Get_Tops(Pythia8::Event event, std::vector<TruthPart>* p_Top_Coll);

  void Get_VectorBosons(Pythia8::Event event, std::vector<TruthPart>* p_VecBoson_Coll);

  void Get_BarePromptLepton(Pythia8::Event event, std::vector<int> vecboson_index, std::vector<TruthPart>* p_BareLept_Coll, std::vector<TruthPart>* p_Neutrino_Coll);
    
  void Bare_Welectron4Vec(Pythia8::Event event, std::vector<TLorentzVector>* p_Born_Coll);

  void Bare_WelectronTruePart(Pythia8::Event event, std::vector<TruthPart>* p_BornElec_Coll);

  void TrueJetsReco(Pythia8::Event event, std::vector<int> partskipped, std::vector<TruthJets>* p_TruthJets_Coll);


};


