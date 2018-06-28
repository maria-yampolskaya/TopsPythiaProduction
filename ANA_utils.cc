//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: ANA_utils.cc
// 
// Purpose: Provide functions to fill physics objects with proper corrections and
//          systematics. Utility functions for various kinematic calculations too.
//
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// Dependencies (#includes)
// ------------------------

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
double ANA_utils::compute_weight(const int NB, const double x_sect, const double instLum) 
//
// To compute the event weigth allowing to combine various samples
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
//         -Fill_TruthPart: Fill a TruthParticle object.
//
//         -Get_Tops: Fill a collection of TruthParticle with the top and anti-tops in the event.
//
//         -Get_VectorBosons: Fill a collection of TruthParticle with all vector bosons in the event.
//
//         -Bare_Welectron4Vec: Fill a Collection of TLorentzVector with bare level electrons from W
//
//         -Bare_WelectronTruePart: Fill a Collection of TruthParticle with bare level electrons from W
//
//         -Get_BarePromptLepton: Fill a collection of TruthParticle with prompt leptons.
//
//         -TrueJetsReco: Reconstruct jets, and fill a Collection of TruthJets
//
//*************************************************************************************************************



// ============================================================================================================
void ANA_utils::Fill_TruthPart(Pythia8::Event event, int index, TruthPart* p_TruthPart)
//
//  Fill a TruthPart object from the particle "index" in the Pythia event record.
//  
//     Note: The object is filled via the pointer set in the function which calls Fill_TruthPart
//  
// ============================================================================================================
{

  // Fill a temporary 4-vector from the particle kinematic
  // -----------------------------------------------------

  double partPt  = event[index].pT();
  double partEta = event[index].eta();
  double partPhi = event[index].phi();
  double partM = event[index].m();

    
  TLorentzVector temp_4vec;
  temp_4vec.SetPtEtaPhiM(partPt,partEta,partPhi,partM);


  // Get the internal Pythia information on this particle
  // ----------------------------------------------------
  
  double charge = event[index].charge();
  int status = event[index].status();
  int statusHepMC = event[index].statusHepMC();
  int pdgid = event[index].id();
  int part_index = index;


  // Fill the TruthPart object
  // -------------------------
  
  p_TruthPart->Set(temp_4vec.E(), temp_4vec.Px(), temp_4vec.Py(), temp_4vec.Pz(), partM, partPt, partEta, partPhi, charge, part_index, status, statusHepMC, pdgid);
  
}

// ============================================================================================================






// ============================================================================================================
void ANA_utils::Get_Tops(Pythia8::Event event, std::vector<TruthPart>* p_Top_Coll)
//
// Find the top and the anti-top in an event using the last entry in the event record before the top decay. Fill
// a TruthPart object with each of the top and anti-top and store them in a collection. Make sure that if there
// is one single top, only that particle will be stored in the collection and that the function won't crash.  
// 
//
// ============================================================================================================
{
  // Find the top particles
  // ----------------------

     // Note: Start from the last entry in the event record and move up the list until a top is found. That is
     //       the last one before it decays.
  
  int top = -1;
  int antitop = -1;

  for (int i = event.size() - 1; i > 0; i--)
    {
      if (event[i].id() == 6 && top == -1) top = i;
      if (event[i].id() == -6 && antitop == -1) antitop = i;
    }

  
  
  // Fill truth particle object and store in the collection for each top
  // -------------------------------------------------------------------

      // Note: Need to define a pointer to the top TruthPart in order for the function Fill_TruthPart to
      //       fill it properly.
  
  TruthPart temp_top;
  TruthPart* p_temp_top;
  TruthPart temp_antitop;
  TruthPart* p_temp_antitop;

  p_temp_top = &temp_top;   
  p_temp_antitop = &temp_antitop; 

  if (top != -1) Fill_TruthPart(event, top, p_temp_top);
  if (antitop != -1) Fill_TruthPart(event, antitop, p_temp_antitop);


    // Store the top particles in the top collection
    // .............................................

  if (top != -1 ) p_Top_Coll->push_back(temp_top);
  if (antitop != -1 )p_Top_Coll->push_back(temp_antitop);
  
}

// ============================================================================================================






// ============================================================================================================
void ANA_utils::Get_VectorBosons(Pythia8::Event event, std::vector<TruthPart>* p_VecBoson_Coll)
//
// Find the W and Z bosons, fill a TruthPart object for each vector bosons and store them in a collection. Note
// that carbon copies and chains of intermediate steps are also stored in the even record, giving many more W
// or Z in the event record than there is in the event process. The last copy is the one that carry the right
// energy and momentum. To find it, the Pythia function Particle::iBotCopyID() looks for copy of the same
// particle down stream. If there is no further copy, the function returns the index of the particle tested.
// The vector bosons to keep are those for which there is no further copy in the "descendence".
//
//    Note 1: This is a different approach for doing the same thing as what is done to find the top and the
//            anti-top. It has been shown that they yield the same result, although the current strategy seems
//            more robust. Both are kept to provide different examples.
//
//    Note 2: For further comments on how this function works, see Get_Tops.
//
// ============================================================================================================
{

  // Find the W or Z bosons
  // ----------------------

  std::vector<int> vecboson_index;
  
  for (int i = event.size() - 1; i > 0; i--)
    {
      if (event[i].idAbs() == 24 || event[i].idAbs() == 23)
	{
	  if (event[i].iBotCopyId()==i) vecboson_index.push_back(i);
	}
    }
    

  // Fill truth particle object and store in the collection for each vector boson 
  // ----------------------------------------------------------------------------

  if (vecboson_index.size() == 0) return;
  
  for (int i_vb = 0; i_vb < vecboson_index.size(); i_vb++)
    {

      TruthPart temp_vb;
      TruthPart* p_temp_vb;


      p_temp_vb = &temp_vb;   


      Fill_TruthPart(event, vecboson_index[i_vb], p_temp_vb);

      
    // Store the truth particle the vector boson collection
    // ....................................................

      p_VecBoson_Coll->push_back(temp_vb);

    }

} 

//------------------------------------------------------------------------------------------------------------






// ============================================================================================================
void ANA_utils::Bare_Welectron4Vec(Pythia8::Event event, std::vector<TLorentzVector>* p_Born_Coll)
//
// Find the electron at Bare level from the decay of a W, and store the 4-vector in the empty collection.
//
//   Note: This function is very simple and works only for events with a single W decaying to an electron. It
//         is kept as a simple example that can be used. Get_BarePromptLepton is a much more general function.
//
// ============================================================================================================
{

  // Find the W boson
  // ----------------
  
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

    
  // Find the bare electron from the W decay
  // ---------------------------------------

  int idxElec = idxW;
  while(true)
    {
      int daughter = event[idxElec].daughter1();
      if   (daughter == 0) break;
      else                 idxElec = daughter;
    }

  
  // Check if the last W "descendent" is a stable electron
  // -----------------------------------------------------
  
  if (event[idxElec].idAbs() != 11 || !event[idxElec].isFinal())
    {
      cout << "Error: Found incorrect decay product of the W" << endl;
      return;
    }
 

  // Fill 4-vector
  // -------------
  
  double elecPt  = event[idxElec].pT();
  double elecEta = event[idxElec].eta();
  double elecPhi = event[idxElec].phi();
  double elecM = event[idxElec].m();
  
  TLorentzVector temp_elec;
  temp_elec.SetPtEtaPhiM(elecPt,elecEta,elecPhi,elecM);

    
    // Store the 4 vector in the bare electron collection
    // ..................................................
    
  p_Born_Coll->push_back(temp_elec);


} 

//------------------------------------------------------------------------------------------------------------






// ============================================================================================================
void ANA_utils::Bare_WelectronTruePart(Pythia8::Event event, std::vector<TruthPart>* p_BornElec_Coll)
//
// Find the electron at Bare level from the decay of a W, fill a true particle objects with the information about
// this electron and store it in the empty collection.
//
//    Note: While the function if more sophisticated than Bare_Welectron4Vec, it is still less general than
//          Get_BarePromptLepton which should be the favored function to find bare leptons. This one is
//          nevertheless kept for reference and for testing of other functions.   
//
// ============================================================================================================
{

  // Find the W or Z bosons
  // ----------------------

  std::vector<int> vecboson_index;
  
  for (int i = event.size() - 1; i > 0; i--)
    {
      if (event[i].idAbs() == 24 || event[i].idAbs() == 23)
	{
	  if (event[i].iBotCopyId()==i) vecboson_index.push_back(i);
	}
    }
    
  // Find the electron from the W decay
  // ----------------------------------
  
  int idxElec = vecboson_index[0];
  while(true)
    {
      int daughter = event[idxElec].daughter1();
      if   (daughter == 0) break;
      else                 idxElec = daughter;
    }

    

  // Check if the last W "descendent" is a stable electron
  // -----------------------------------------------------

  if (event[idxElec].idAbs() != 11 || !event[idxElec].isFinal())
    {
      cout << "Error: Found incorrect decay product of the W" << endl;
      return;
    }
 
      
  // Fill truth particle object and store in the collection
  // ------------------------------------------------------

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
  
  temp_truth.Set(temp_elec.E(), temp_elec.Px(), temp_elec.Py(), temp_elec.Pz(), elecM, elecPt, elecEta, elecPhi, charge, part_index, status, statusHepMC, pdgid);


    // Store the 4 vector in the bare electron collection
    // ..................................................

  p_BornElec_Coll->push_back(temp_truth);


} 

//------------------------------------------------------------------------------------------------------------






// ============================================================================================================
void ANA_utils::Get_BarePromptLepton(Pythia8::Event event, std::vector<int> vecboson_index, std::vector<TruthPart>* p_BareLept_Coll, std::vector<TruthPart>* p_Neutrino_Coll)
//
// Find all leptons at Bare level coming from the decay of a W or a Z. Fill a true particle objects collection
// with the information about these leptons. This function is valid for both electrons and muons, and for events
// with single or multiple vector bosons (e.g. diboson, ttbar, etc.). The StableDaughterList vector keeps all
// the stable particles descending from a vector boson. It therefore has all the photons for dressing of leptons
// if needed. These photons are not stored for the moment but can easily be.
//
//    Note 1: Because of the QED emission of the charged leptons, and of the carbon copies of the particles stored
//            in the Pythia event record, the immidiate daughters of a vector bosons immediately before it decays
//            is only rarely the stable particles to be used at bare level. It is therefore important to follow
//            the entire decay chain until such stable particles are reached. The number of generation to cover
//            varies for different vector bosons in one event, and even for different leptons. For example, if one
//            Z lepton emits two photons, but not the other, the second lepton lepton will be stable from the first
//            generation of the Z daugthers, while the first one will only be stable at the third generation, and
//            the two photons will respectively be stable at the second and third generations. To this end, a
//            complete enquiry of each daughters at each generation is required. After a stress test of 1000 events,
//            there was no need to further than the 6th generation. The code therefore stops at the 7th generation.
//  
//    Note 2: It happens in a few events that there are two more charged leptons than expected. This is the result
//            of a photon conversion into a pair of lepton-anti-lepton. These are kept as stable leptons assigned
//            to the decay of a vector boson. A cut on the number of charged leptons expected in the event will
//            take care of these. It happens both in the predictions and in the data.
//  
// ============================================================================================================
{

    
  // Find the stable lepton from the vector boson decays
  // ---------------------------------------------------

  std::vector<int> lepton_index;
  std::vector<int> neutrino_index;

  
     // Loop over all vector bosons found in the event
     // ..............................................

        // Note: This require using the function Get_VectorBosons, and passing the index of each vector boson as
        //       input to this function
  
  for (int i=0; i< vecboson_index.size(); i++)
    {
      int idxVb = vecboson_index[i];      
      

        // Define vectors to keep the indices of each particles in each generations
        // ........................................................................
      
      std::vector<int> StableDaughterList;

      std::vector<int> gen1;
      std::vector<int> gen2;
      std::vector<int> gen3;
      std::vector<int> gen4;
      std::vector<int> gen5;
      std::vector<int> gen6;
      std::vector<int> gen7;


        // The first generation is straighforwardly the daughters of the vector boson
        // ..........................................................................
      
      gen1.push_back(event[idxVb].daughter1());
      gen1.push_back(event[idxVb].daughter2());


        // If the vector bosons have not been forced to be stable, checked if the daughters are stable
        // ...........................................................................................

            // Note: Any unstable daughter will have their own daughters populating the next generation
      
      if (gen1.size() != 0)
	{
	  if ( event[gen1[0]].isFinal()  )  StableDaughterList.push_back(gen1[0]);

	  else
	    {

	      for (int d1 = 0; d1 <  (event[gen1[0]].daughterList()).size(); d1++)
		{
		  gen2.push_back((event[gen1[0]].daughterList())[d1]);
		}
	    }

	  if ( event[gen1[1]].isFinal()  )  StableDaughterList.push_back(gen1[1]);

	  else
	    {
	      for (int d2 = 0; d2 <  (event[gen1[1]].daughterList()).size(); d2++)
		{
		  gen2.push_back((event[gen1[1]].daughterList())[d2]);
		}

	    }

	  
        // Check if the particles in the second generation are stable, if not moved to the next
	// ....................................................................................  

	  if (gen2.size() != 0)
	    {

      
	      for (int i_gen2 = 0; i_gen2 < gen2.size(); i_gen2++)
		{
		  if ( event[gen2[i_gen2]].isFinal()  ) StableDaughterList.push_back(gen2[i_gen2]);
		  else
		    {
		      for (int d3 = 0; d3 <  (event[gen2[i_gen2]].daughterList()).size(); d3++)
			{
			  gen3.push_back((event[gen2[i_gen2]].daughterList())[d3]);
			}
		      
		    }
		}

	      
        // Check if the particles in the third generation are stable, if not moved to the next
	// ....................................................................................  

	      if (gen3.size() != 0)
		{

		  for (int i_gen3 = 0; i_gen3 < gen3.size(); i_gen3++)
		    {
		      if ( event[gen3[i_gen3]].isFinal()  ) StableDaughterList.push_back(gen3[i_gen3]);
		      else
			{
			  for (int d4 = 0; d4 <  (event[gen3[i_gen3]].daughterList()).size(); d4++)
			    {
			      gen4.push_back((event[gen3[i_gen3]].daughterList())[d4]);
			    }

			}
		    }


        // Check if the particles in the forth generation are stable, if not moved to the next
	// ....................................................................................  
		  
		  if (gen4.size() != 0)
		    {
		      for (int i_gen4 = 0; i_gen4 < gen4.size(); i_gen4++)
			{
			  if ( event[gen4[i_gen4]].isFinal()  ) StableDaughterList.push_back(gen4[i_gen4]);
			  else
			    {
			      for (int d5 = 0; d5 <  (event[gen4[i_gen4]].daughterList()).size(); d5++)
				{
				  gen5.push_back((event[gen4[i_gen4]].daughterList())[d5]);
				}
			      
			    }
			}


	// Check if the particles in the fith generation are stable, if not moved to the next
	// ....................................................................................  
		      
		      if (gen5.size() != 0)
			{
			  for (int i_gen5 = 0; i_gen5 < gen5.size(); i_gen5++)
			    {
			      if ( event[gen5[i_gen5]].isFinal()  ) StableDaughterList.push_back(gen5[i_gen5]);
			      else
				{
				  for (int d6 = 0; d6 <  (event[gen5[i_gen5]].daughterList()).size(); d6++)
				    {
				      gen6.push_back((event[gen5[i_gen5]].daughterList())[d6]);
				    }

				}
			    }

			  
        // Check if the particles in the sixth generation are stable, if not moved to the next
	// ....................................................................................  

			  if (gen6.size() != 0)
			    {
			      for (int i_gen6 = 0; i_gen6 < gen6.size(); i_gen6++)
				{
				  if ( event[gen6[i_gen6]].isFinal()  ) StableDaughterList.push_back(gen6[i_gen6]);
				  else
				    {
				      for (int d7 = 0; d7 <  (event[gen6[i_gen6]].daughterList()).size(); d7++)
					{
					  gen7.push_back((event[gen6[i_gen6]].daughterList())[d7]);
					}
				      
				    }
				}

			    } //end gen6!=0
			} //end gen5!=0
		    } // end if gen4!=0
		} // end if gen3!=0
	    } // end if gen2!=0
	} // end if gen1!=0


        // Of all the stable particles, find those that are leptons (electrons or muons)
        // .............................................................................
      
      for (int stab_i = 0; stab_i < StableDaughterList.size(); stab_i++)
	{
	  int temp_id = StableDaughterList[stab_i];
	  
	  if ( event[temp_id].idAbs() == 11 || event[temp_id].idAbs() == 13 ) lepton_index.push_back(StableDaughterList[stab_i]);
	  if ( event[temp_id].idAbs() == 12 || event[temp_id].idAbs() == 14  || event[temp_id].idAbs() == 16 ) neutrino_index.push_back(StableDaughterList[stab_i]);
	}

    } // end loop over vector bosons



  // Fill truth particle object and store in the collection for each stable charged lepton 
  // -------------------------------------------------------------------------------------

  if (lepton_index.size() == 0)
    {
  
      for (int i_lep = 0; i_lep < lepton_index.size(); i_lep++)
	{

	  TruthPart temp_lep;
	  TruthPart* p_temp_lep;


	  p_temp_lep = &temp_lep;   


	  Fill_TruthPart(event, lepton_index[i_lep], p_temp_lep);

      
    // Store the 4 vector in the bare electron collection
    // ..................................................

	  p_BareLept_Coll->push_back(temp_lep);

	}
    }
  
  
  // Fill truth particle object and store in the collection for each neutrino
  // ------------------------------------------------------------------------

  if (neutrino_index.size() == 0) return;
  
  for (int i_neu = 0; i_neu < neutrino_index.size(); i_neu++)
    {

      TruthPart temp_neu;
      TruthPart* p_temp_neu;


      p_temp_neu = &temp_neu;   


      Fill_TruthPart(event, neutrino_index[i_neu], p_temp_neu);

      
    // Store the 4 vector in the bare electron collection
    // ..................................................

      p_Neutrino_Coll->push_back(temp_neu);

    }
  

} 

//------------------------------------------------------------------------------------------------------------






// ============================================================================================================
void ANA_utils::TrueJetsReco(Pythia8::Event event, std::vector<int> partskipped, std::vector<TruthJets>* p_TruthJets_Coll)
//
// Reconstruct the truth-level jets using FastJets, sort them in pT, fill a TruthJets objects with the relevant
// variables for each jet, and store it in the empty collection.
//
// ============================================================================================================
{

  // Build Jets
  // ----------
  
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

  
     // Reset Fastjet input
     // ...................
  
  fjInputs.resize(0);

  
     // Loop over event record to decide what to pass to FastJet
     // ........................................................
  
  for (int i = 0; i < event.size(); ++i)
    {

      // Final state only
      if (!event[i].isFinal())        continue;


      // No neutrinos
      if (event[i].idAbs() == 12 || event[i].idAbs() == 14 ||
          event[i].idAbs() == 16)     continue;


      // Only |eta| < 4.9
      if (fabs(event[i].eta()) > 4.9) continue;


      // Do not include the stable particles listed in the partskipped vector

          // Note: This is used to not include the decay product of the W, Z, H and top for example

      bool skip_i = false;
      for (int j = 0; j < partskipped.size(); j++)
	{
	  if (partskipped[j] == i) skip_i = true;
	}
      if (skip_i == true) continue;
      

      // Store as input to Fastjet
      
      fjInputs.push_back( fastjet::PseudoJet( event[i].px(),
        event[i].py(), event[i].pz(), event[i].e() ) );
    }

     // Check that there is some input to Fastjet
     // .........................................
  
  if (fjInputs.size() == 0)
    {
      cout << "Error: event with no final state particles" << endl;
      return;
    }


     // Run Fastjet algorithm
     // .....................
  
  vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
  fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);


     // Extract inclusive jets sorted by pT (note minimum pT of 30.0 GeV)
     // .................................................................
  
  inclusiveJets = clustSeq.inclusive_jets(30.0);
  sortedJets    = sorted_by_pt(inclusiveJets);



  // Fill the jet object and store it in the true jet collection
  // -----------------------------------------------------------

     // Loop over all jets produced by Fastjet
     // ......................................
  
  for (int i_jet = 0; i_jet < sortedJets.size(); i_jet++)
    {
      bool btag;

     // Check if the jet is a B-jet
     // ...........................
      
      btag = false;


      // Fill a TruthJets object
      // .......................
      
      TruthJets temp_truthjet;
      temp_truthjet.Set((sortedJets[i_jet]).e(), (sortedJets[i_jet]).px(), (sortedJets[i_jet]).py(), (sortedJets[i_jet]).pz(), (sortedJets[i_jet]).m(), (sortedJets[i_jet]).pt(), (sortedJets[i_jet]).eta(), (sortedJets[i_jet]).rap(), (sortedJets[i_jet]).phi(), btag);

      // Put in the jet collection
      // .........................
      
      p_TruthJets_Coll->push_back(temp_truthjet);
      
    }



  // Close objects
  // -------------
  
  delete jetDef;

  
}
