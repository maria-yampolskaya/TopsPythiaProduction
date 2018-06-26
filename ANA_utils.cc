//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: ANA_utils.hpp
//
// Purpose: Provide functions to fill physics objects with proper corrections and
//          systematics. Utility functions for various kinematic calculations too.
//
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// Dependencies (#includes)
// ------------------------

//#include "VjetsPythia8.h"
#include "ANA_utils.h"

   // C++
   // ...

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cstring>
#include <stdio.h>
#include <ios>
#include <fstream>
#include <cstdio>


using namespace std;


//*************************************************************************************************************
//
//    Event Handling Functions:
//
//         -compute_weight
//
//*************************************************************************************************************



// ============================================================================================================
//double MyAnalysis::compute_weight(const int NB, const double x_sect, const double instLum)
double ANA_utils::compute_weight(const int NB, const double x_sect, const double instLum)
//
// To compute the Alpgen weigth allowing to combine the various Np samples
//
// ============================================================================================================
{
  double weight;

  weight = x_sect/NB*instLum;

  return weight;

}

//------------------------------------------------------------------------------------------------------------






//*************************************************************************************************************
//
//    Truth Level Analysis Functions:
//
//         -Born_Welectron
//
//*************************************************************************************************************



// ============================================================================================================
void ANA_utils::Born_Welectron(Pythia8::Event event, std::vector<TLorentzVector>* p_Born_Coll)
//
// Find the electron at Born from the decay of a W
//
// ============================================================================================================
{


  int idxW = -1;

  for (int i = event.size() - 1; i > 0; i--)
    {
      if (event[i].idAbs() == 24)
	{
	  idxW = i;
	  break;
	}
    }
    if (idxW == -1)
      {
	cout << "Error: Could not find W" << endl;
	return;
      }


    // Find the electron from the W decay
    int idxElec = idxW;
    while(true)
      {
	int daughter = event[idxElec].daughter1();
	if   (daughter == 0) break;
	else                 idxElec = daughter;
      }
    if (event[idxElec].idAbs() != 11 || !event[idxElec].isFinal())
      {
	cout << "Error: Found incorrect decay product of the W" << endl;
	return;
      }


    double elecPt  = event[idxElec].pT();
    double elecEta = event[idxElec].eta();
    double elecPhi = event[idxElec].phi();
    double elecM = event[idxElec].m();

    TLorentzVector temp_elec;
    temp_elec.SetPtEtaPhiM(elecPt,elecEta,elecPhi,elecM);
    p_Born_Coll->push_back(temp_elec);


    TruthPart mytest;

    mytest.Set(0.01, 1.02, 2.03, 3.04, 4.05, 5.06, 6.07, 7.08, 1.1, 12, 3, 33, -11);

    //std::cout << "Pt = " << mytest.Pt() << ", Px = " << mytest.Px() << ", Py = " << mytest.Py() << ", Pz = " << mytest.Pz() << std::endl;


}

//------------------------------------------------------------------------------------------------------------





// ============================================================================================================
//void Born_Welectron2(Pythia8::Event event)
//void Born_Welectron2(Pythia8::Event event, std::vector<TruthPart> p_BornElec_Coll)
void Born_Welectron2(std::vector<TruthPart> p_BornElec_Coll)
//
// Find the electron at Born from the decay of a W
//
// ============================================================================================================
{

  /*
  int idxW = -1;

  for (int i = event.size() - 1; i > 0; i--)
    {
      if (event[i].idAbs() == 24)
	{
	  idxW = i;
	  break;
	}
    }
    if (idxW == -1)
      {
	cout << "Error: Could not find W" << endl;
	return;
      }


    // Find the electron from the W decay
    int idxElec = idxW;
    while(true)
      {
	int daughter = event[idxElec].daughter1();
	if   (daughter == 0) break;
	else                 idxElec = daughter;
      }
    if (event[idxElec].idAbs() != 11 || !event[idxElec].isFinal())
      {
	cout << "Error: Found incorrect decay product of the W" << endl;
	return;
      }


    double elecPt  = event[idxElec].pT();
    double elecEta = event[idxElec].eta();
    double elecPhi = event[idxElec].phi();
    double elecM = event[idxElec].m();


    TLorentzVector temp_elec;
    temp_elec.SetPtEtaPhiM(elecPt,elecEta,elecPhi,elecM);

      double charge = event[idxElec].charge();
      int status = event[idxElec].status();
      int statusHepMC = event[idxElec].statusHepMC();
      int pdgid = event[idxElec].id();
      int part_index = idxElec;

    TruthPart temp_truth;
  */
    //temp_truth.Set(temp_elec.E(), temp_elec.Px(), temp_elec.Py(), temp_elec.Pz(), elecM, elecPt, elecEta, elecPhi, charge, part_index, status, statusHepMC, pdgid);
    //p_BornElec_Coll->push_back(temp_truth);

}

//------------------------------------------------------------------------------------------------------------

// ============================================================================================================
void ANA_utils::find_bquarks(Pythia8::Event event, std::vector<TruthPart>* p_bquarks)
//
// Find the bottom quarks from the initial processes
//
// ============================================================================================================
{

  std::vector<int> bquark_indices;

  for (int i = 1; i < event.size(); i++) {
      int mother = event[i].mother1();
      if ((event[i].idAbs() == 5) && (mother == 0 || event[mother].idAbs() == 6 || event[mother].idAbs() == 23 || event[mother].idAbs() == 24 || event[mother].idAbs() == 25)) {
	        bquark_indices.push_back(i);
	     }
    }

  for (auto i : bquark_indices) {
      TruthPart temp_b;
      temp_b.Set_M(event[i].m());
      temp_b.Set_Phi(event[i].phi());
      temp_b.Set_Eta(event[i].eta());
      temp_b.Set_Pt(event[i].pT());

      p_bquarks->push_back(temp_b);
    }

}
