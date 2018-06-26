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


  // My Includes
  // ...........

#include "TruthPart.h"


#include "Pythia8/Pythia.h"

using namespace std;
using namespace Pythia8;





// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
// class ANARb : public RootTuple_bu
//
//   Note: Analysis class definition
//
// -----------------------------------------------------------------------------------
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

class ANA_utils
{

public :

  // Constructor and destructor
  // --------------------------


//============================================================================================================
//
// Constructor and destructor
//
//============================================================================================================

  ANA_utils(){}

  ~ANA_utils(){}

//------------------------------------------------------------------------------------------------------------




  // Declare functions
  // -----------------

     // Note: Utility functions are defined in ANARb_utils.cpp, while the
     //       main analysis function is defined in ANARb_Muon.cpp or ANARb_Elec.cpp






     // To compute alpgen weight
     // ........................

  double compute_weight(const int NB, const double x_sect, const double inst_lum);

  void Born_Welectron(Pythia8::Event event, std::vector<TLorentzVector>* p_Born_Coll);

  //void Born_Welectron2(Pythia8::Event event, std::vector<TruthPart> p_BornElec_Coll);
  void Born_Welectron2(std::vector<TruthPart> p_BornElec_Coll);
  //void Born_Welectron2(Pythia8::Event event);

  void find_bquarks(Pythia8::Event event, std::vector<TruthPart>* p_bquarks);






//private:

  // Define internal variables
  // -------------------------

//  int _nmax;



//public:

//  Float_t _WEIGHT; // weight for histos


};
