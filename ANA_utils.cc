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
//         -compute_weight: Weight event for lumi
//         -getPartonLevelEvent: fill a pythia Event object with particles at parton level
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






// ============================================================================================================
void  ANA_utils::getPartonLevelEvent( Event& event, Event& partonLevelEvent) 
//
// A generic routine to extract the particles that existed right before the hadronization machinery is invoked
// by Pythia. This is useful to get parton-level quantities, but after the entire shower history, and not
// randomly between the hard process and the end-point of the evolution.
//
// ============================================================================================================
{


  partonLevelEvent.reset();
  

  // Loop over the entire event to select parton level particles just before hadronization
  // -------------------------------------------------------------------------------------
  
  for (int i = 0; i < event.size(); ++i)
    {
      bool accept = false;


      // Only partons after full evolution and before hadronization are kept
      // ...................................................................
      
      if (event[i].isFinalPartonLevel()) accept = true;


      // Don't keep neutrinos and charged leptons because they come from hard process
      // ............................................................................

      int idAbs = event[i].idAbs();
      if (idAbs >10 && idAbs < 17) accept = false;

      
      // Reject particles outside calorimeter acceptance
      // ...............................................

      if (event[i].eta() > 4.9) accept = false;
    


      if (accept == true)
	{
	  int iNew = partonLevelEvent.append( event[i] );

	  
  // Set copied properties more appropriately
  // ----------------------------------------

	  // Note: Set a positive status, original location as "mother", and with no daughters.
	  
	  partonLevelEvent[iNew].statusPos();
	  partonLevelEvent[iNew].mothers( i, i);
	  partonLevelEvent[iNew].daughters( 0, 0);
	}
    }
}

//--------------------------------------------------------------------------






//*************************************************************************************************************
//
//    Truth Level Analysis Functions:
//
//         -Fill_TruthPart: Fill a TruthParticle object.
//
//         -Get_Tops: Fill a collection of TruthParticle with the top and anti-tops in the event.
//
//         -Get_VectorBosons: Fill a collection of TruthParticle with all prompt vector bosons in an event.
//
//         -Bare_Welectron4Vec: Fill a Collection of TLorentzVector with bare level electrons from W
//
//         -Bare_WelectronTruePart: Fill a Collection of TruthParticle with bare level electrons from W
//
//         -Get_BarePromptLepton: Fill a collection of TruthParticle with prompt leptons and another one for neutrinos.
//
//         -TrueJetsReco: Reconstruct jets, and fill a Collection of TruthJets, and deal with b-tagging.
//
//         -Get_BottomQuarks: Fill a collection of TruthParticle with the bottom and anti-bottoms in the event.
//
//         -Get_BottomHadrons: Fill a collection of TruthParticle with the B-hadrons from promt b-quarks.
//
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

// ============================================================================================================







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

// ============================================================================================================






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

// ============================================================================================================






// ============================================================================================================
void ANA_utils::Get_BarePromptLepton(Pythia8::Event event, std::vector<int> vecboson_index, std::vector<TruthPart>* p_BareLept_Coll, std::vector<TruthPart>* p_Neutrino_Coll)
//
// Find all leptons at Bare level coming from the decay of a W or a Z. Fill a true particle objects collection
// with the information about these leptons. This function is valid for both electrons and muons, and for events
// with single or multiple vector bosons (e.g. diboson, ttbar, etc.). The StableDaughterList vector keeps all
// the stable particles descending from a vector boson. It therefore has all the photons for dressing of leptons
// if needed. These photons are not stored for the moment but can easily be. Neutrinos are also stored. 
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

  if (lepton_index.size() != 0)
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

// ============================================================================================================






// ============================================================================================================
void ANA_utils::TrueJetsReco(Pythia8::Event event, std::vector<int> partskipped, std::vector<TruthJets>* p_TruthJets_Coll)
//
// Reconstruct the truth-level jets using FastJets, sort them in pT, fill a TruthJets objects with the relevant
// variables for each jet, and store it in the empty collection. A b-tagger is also implemented in this function.
// To this end, first, the b-quarks or B-hadrons are found, then they are used to creat ghost particles. These
// ghosts are clustered with the jets, but they do not impact jet reconstruction. In a last step, the constituents
// of each jets are investigated to see if they have a ghost or not. If they do, they will be b-tagged.
//
//    Note 1: A dR b-quark/b-hadron assignment is also implemented for comparison, but the ghost approach is prefered.
//
//    Note 2: The "partskipped" vectors given the particle event index of the stable particles not to be clustered
//            in the jets, such as final state leptons coming from W and Z for example.
//
//    The code contains further notes below.   
//
// ============================================================================================================
{

  // Find b-quarks and b-hadrons in the event
  // ----------------------------------------

      // Note: There is a very strong correlation between the kinematic of the b-quark and the kinematic of the b-hadron
      //       following the b-quark hadronization. When printing the pT and eta of the b-quarks and the b-hadrons, one 
      //       would see that these quantities have very similar values for b-quarks and b-hadrons. The correlation is 
      //       however much weaker with the bjets. 
  
  std::vector<TruthPart> BottomQuark_Coll;
  std::vector<TruthPart>* p_BottomQuark_Coll = &BottomQuark_Coll;

  Get_BottomQuarks(event, p_BottomQuark_Coll);


  std::vector<TruthPart> BottomHadron_Coll;
  std::vector<TruthPart>* p_BottomHadron_Coll = &BottomHadron_Coll;

  Get_BottomHadrons(event, BottomQuark_Coll, p_BottomHadron_Coll);

  
  // Create ghost particles with each b-quark or b-hadron
  // ----------------------------------------------------

      // Note: A ghost particle is a particle that exactly has the direction of a given particle, but has its energy to almost 0.
      //       Such ghost particle are useful for b-quarks or b-hadrons association to jets, because the ghost will be clustered
      //       in the jet, but the energy and direction of the jet will not be affected at all by the ghost. By looking to which 
      //       jet the ghost has been clustered, we'll know to which jet the b-quark or b-hadron can be assigned.

  std::vector<TLorentzVector> Ghost_Bquarks;
  std::vector<TLorentzVector> Ghost_Bhadrons;
  

      // Create ghosts from b-quarks
      // ...........................
  
  for (int i_bq = 0; i_bq<BottomQuark_Coll.size(); i_bq++)
    {
      TLorentzVector ghost_bqrk;
      ghost_bqrk.SetPtEtaPhiM( ((BottomQuark_Coll[i_bq]).Pt()/1000000.), (BottomQuark_Coll[i_bq]).Eta(), (BottomQuark_Coll[i_bq]).Phi(), ((BottomQuark_Coll[i_bq]).M()/1000000.));
      
      Ghost_Bquarks.push_back(ghost_bqrk);
    }


      // Create ghosts from b-hadrons
      // ............................
  
  for (int i_bh = 0; i_bh<BottomHadron_Coll.size(); i_bh++)
    {
      TLorentzVector ghost_bhad;
      ghost_bhad.SetPtEtaPhiM( ((BottomHadron_Coll[i_bh]).Pt()/1000000.), (BottomHadron_Coll[i_bh]).Eta(), (BottomHadron_Coll[i_bh]).Phi(), ((BottomHadron_Coll[i_bh]).M()/1000000.));
      
      Ghost_Bhadrons.push_back(ghost_bhad);
    }

  
  // Build Inclusive Jets
  // --------------------

     // FastJet Initialization
     // ......................

  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
				      recombScheme, strategy);
  
     // Fastjet input
     // .............
  
  std::vector <fastjet::PseudoJet> fjInputs;

  
     // Reset Fastjet input
     // ...................
  
  fjInputs.resize(0);
  int index=0;

  
     // Loop over event record and specify which particles to pass to FastJet
     // .....................................................................
  
  for (int i = 0; i < event.size(); ++i)
    {

      // Final state only
      if (!event[i].isFinal())        continue;


      // No neutrinos
      if (event[i].idAbs() == 12 || event[i].idAbs() == 14 || event[i].idAbs() == 16)     continue;


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
      

      // Store as input to Fastjet and set a unique identifier for each input particle
      // .............................................................................

      fastjet::PseudoJet particle( event[i].px(),event[i].py(), event[i].pz(), event[i].e() );

      particle.set_user_index(index);
      fjInputs.push_back( particle );
      index++;
    }


     // Add the ghost to the input particles
     // ....................................

         // Note: In contrary to the stable particles constituting the jets, the index value used for
         //       ghost particles is negative, with 0 < ghost_index < -10 for b-quark ghosts, and
         //       with -9 < ghost_index < -100 for b-hadron ghosts.
  
  int ghost_index = -1;
  
  for (int i_gbq = 0; i_gbq<Ghost_Bquarks.size(); i_gbq++)
    {
      fastjet::PseudoJet ghost_particle( (Ghost_Bquarks[i_gbq]).Px(), (Ghost_Bquarks[i_gbq]).Py(), (Ghost_Bquarks[i_gbq]).Pz(), (Ghost_Bquarks[i_gbq]).E() );

      ghost_particle.set_user_index(ghost_index);
      fjInputs.push_back( ghost_particle );
      ghost_index--;
    }


  for (int i_gbh = 0; i_gbh<Ghost_Bhadrons.size(); i_gbh++)
    {
      fastjet::PseudoJet ghost_particle( (Ghost_Bhadrons[i_gbh]).Px(), (Ghost_Bhadrons[i_gbh]).Py(), (Ghost_Bhadrons[i_gbh]).Pz(), (Ghost_Bhadrons[i_gbh]).E() );
      
      ghost_index = ghost_index - 10;
      ghost_particle.set_user_index(ghost_index);
      fjInputs.push_back( ghost_particle );
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


     // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
     // .................................................................
  
  inclusiveJets = clustSeq.inclusive_jets(20.0);
  sortedJets    = sorted_by_pt(inclusiveJets);



  // Fill the jet object and store it in the true jet collection
  // -----------------------------------------------------------

     // Loop over all jets produced by Fastjet
     // ......................................


  for (int i_jet = 0; i_jet < sortedJets.size(); i_jet++)
    {

     
     // Get the constituents of the jet
     // ...............................

      std::vector<fastjet::PseudoJet> constituents = (sortedJets[i_jet]).constituents();
      
  

     // Check if the jet is a B-jet
     // ...........................

        // Note: We have two methods for btagging:
        //
        //            1- Assign a b-quark or a b-hadron to a jet if its angular distance to this jet is smaller than 0.4
        //
        //            2- Include a b-quark or b-hadron ghost (particle pointing in the same direction but with almost no energy and mass)
        //               in the jet clustering and determine in which jets are the ghosts.
        //
        //       These methods are only approximative because the b-quarks hadronize and the b-hadrons decay, and the product could be far
        //       from each others, have low pT or end up in different jets. It is frequent that a b-quark or a b-hadron doesn't get assigned
        //       to any jet. If the jet threshold is low or the jet size is large, the b-quark or b-hadron can be assigned to many jets, all
        //       close in distance compared to the jet size paramater. These problems do not affect the ghost approach which, by construction
        //       always assign each quark or each hadron to a single jet. It is however possible that two b-quarks or b-hadrons get assigned
        //       to the same jet. It is more frequent for larger jet size parameter (about <0.1% for jets of R=0.1, 0.5% for jets of R=0.4 and
        //       about 5% for jets of R=0.8). Also  the assigned jet could be the ghost itself (so no real jet assignment), or it could be  
        //       assign to a jet that does not contain the bulk of the original b-quark or b-hadron energy.


        // dR Approach
      
      bool bqtag = false;
      bool bhtag = false;

      for (auto bottom : BottomQuark_Coll)
	{
	  if (sqrt(pow(bottom.Eta() - (sortedJets[i_jet]).eta(),2) + pow(bottom.Phi() - (sortedJets[i_jet]).phi(),2) <= 0.4))
	    {
	      bqtag = true;
	    }
	}

      for (auto bottom : BottomHadron_Coll)
	{
	  if (sqrt(pow(bottom.Eta() - (sortedJets[i_jet]).eta(),2) + pow(bottom.Phi() - (sortedJets[i_jet]).phi(),2) <= 0.4))
	    {
	      bhtag = true;
	    }
	}
	

      
        // Ghost Approach

      bool bq_ghosttag = false;
      bool bh_ghosttag = false;
      int ghost_mult = 0;
      
      for (int i_const = 0; i_const<constituents.size(); i_const++)
	{
	  if ( (constituents[i_const]).user_index() < 0 && (constituents[i_const]).user_index() > -10 )
	    {
	      bq_ghosttag = true;
	      ghost_mult++;
	    }

	  if ( (constituents[i_const]).user_index() < -9  )
	    {
	      bh_ghosttag = true;
	      ghost_mult++;
	    }
	    
	}

      
      int npart_injet = constituents.size() - ghost_mult;


      // Fill a TruthJets object
      // .......................
      
      TruthJets temp_truthjet;
      temp_truthjet.Set((sortedJets[i_jet]).e(), (sortedJets[i_jet]).px(), (sortedJets[i_jet]).py(), (sortedJets[i_jet]).pz(), (sortedJets[i_jet]).m(), (sortedJets[i_jet]).pt(), (sortedJets[i_jet]).eta(), (sortedJets[i_jet]).rap(), (sortedJets[i_jet]).phi(), bq_ghosttag, bh_ghosttag, npart_injet);


      // Put in the jet collection
      // .........................
      
      p_TruthJets_Coll->push_back(temp_truthjet);
      
    }


  // Close objects
  // -------------
  
  delete jetDef;

  
}

// ============================================================================================================






// ============================================================================================================
void ANA_utils::PartonJetsReco(Pythia8::Event event, Pythia8::Event partonevent, std::vector<TruthJets>* p_PartonJets_Coll)
//
// Reconstruct the parton-level jets using FastJets, sort them in pT, fill a TruthJets objects with the relevant
// variables for each jet, and store it in the empty collection. A b-tagger is also implemented in this function.
// To this end, the entire event record is used to find the b-quarks after all QCD radiation. It is then used to
// create ghost particles. Then the parton-only event record is used to build the jets. The ghosts are also included
// in this clustering without impacting the final parton-level jets, just before hadronization. As for the Truth jets
// defined in TruthJetsReco, in the last step the constituents of each jets are investigated to see if they have a
// ghost or not. If they do, they will be b-tagged.
//
//
//    Note: The "partskipped" vectors used in the TruthJetsReco function is useless here because the event record
//          is not the same. It is however easier to get the lepton to ignore here because they come from the
//          hard interaction and therefore have a status value between 21 and 29. We can simply ignore them. 
//
//    More explanations can be found in the function TrueJetsReco. 
//
// ============================================================================================================
{

  
  // Find b-quarks in the total event
  // --------------------------------

  std::vector<TruthPart> BottomQuark_Coll;
  std::vector<TruthPart>* p_BottomQuark_Coll = &BottomQuark_Coll;

  Get_BottomQuarks(event, p_BottomQuark_Coll);


  
  // Create ghost particles with each b-quark
  // ----------------------------------------

  std::vector<TLorentzVector> Ghost_Bquarks;
  
  for (int i_bq = 0; i_bq<BottomQuark_Coll.size(); i_bq++)
    {
      TLorentzVector ghost_bqrk;
      ghost_bqrk.SetPtEtaPhiM( ((BottomQuark_Coll[i_bq]).Pt()/1000000.), (BottomQuark_Coll[i_bq]).Eta(), (BottomQuark_Coll[i_bq]).Phi(), ((BottomQuark_Coll[i_bq]).M()/1000000.));
      
      Ghost_Bquarks.push_back(ghost_bqrk);
    }


  
  // Build Inclusive Jets with FastJet
  // ---------------------------------

     // FastJet Initialization
     // ......................

  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam,
				      recombScheme, strategy);
  
     // Fastjet input
     // .............
  
  std::vector <fastjet::PseudoJet> fjInputs;

  
     // Reset Fastjet input
     // ...................
  
  fjInputs.resize(0);
  int index=0;

  
     // Loop over the parton-event record and set the fastjet input particles
     // .....................................................................

        // Note: There is no need to remove unwanted particles because they have already be removed from the parton event record

  
  for (int i = 0; i < partonevent.size(); ++i)
    {

      // Store as input to Fastjet and set a unique identifier for each input particle
      // .............................................................................

      fastjet::PseudoJet particle( partonevent[i].px(), partonevent[i].py(), partonevent[i].pz(), partonevent[i].e() );

      particle.set_user_index(index);
      fjInputs.push_back( particle );
      index++;
    }


     // Add the ghost to the input particles
     // ....................................

         // Note: In contrary to the stable particles constituting the jets, the index value used for
         //       ghost particles is negative, with 0 < ghost_index < -10 for b-quark ghosts.
  
  int ghost_index = -1;
  
  for (int i_gbq = 0; i_gbq<Ghost_Bquarks.size(); i_gbq++)
    {
      fastjet::PseudoJet ghost_particle( (Ghost_Bquarks[i_gbq]).Px(), (Ghost_Bquarks[i_gbq]).Py(), (Ghost_Bquarks[i_gbq]).Pz(), (Ghost_Bquarks[i_gbq]).E() );

      ghost_particle.set_user_index(ghost_index);
      fjInputs.push_back( ghost_particle );
      ghost_index--;
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


     // Extract inclusive jets sorted by pT (note minimum pT of 20.0 GeV)
     // .................................................................
  
  inclusiveJets = clustSeq.inclusive_jets(20.0);
  sortedJets    = sorted_by_pt(inclusiveJets);



  // Fill the jet object and store it in the true jet collection
  // -----------------------------------------------------------

     // Loop over all jets produced by Fastjet
     // ......................................


  for (int i_jet = 0; i_jet < sortedJets.size(); i_jet++)
    {

     
     // Get the constituents of the jet
     // ...............................

      std::vector<fastjet::PseudoJet> constituents = (sortedJets[i_jet]).constituents();
  
  

     // Check if the jet is a B-jet
     // ...........................

      bool bq_ghosttag = false;
      int ghost_mult = 0;
      
      for (int i_const = 0; i_const<constituents.size(); i_const++)
	{
	  if ( (constituents[i_const]).user_index() < 0 && (constituents[i_const]).user_index() > -10 )
	    {
	      bq_ghosttag = true;
	      ghost_mult++;
	    }
	}

      
      int npart_injet = constituents.size() - ghost_mult;


      // Fill a TruthJets object
      // .......................
      
      TruthJets temp_truthjet;
      temp_truthjet.Set((sortedJets[i_jet]).e(), (sortedJets[i_jet]).px(), (sortedJets[i_jet]).py(), (sortedJets[i_jet]).pz(), (sortedJets[i_jet]).m(), (sortedJets[i_jet]).pt(), (sortedJets[i_jet]).eta(), (sortedJets[i_jet]).rap(), (sortedJets[i_jet]).phi(), bq_ghosttag, false, npart_injet);


      // Put in the jet collection
      // .........................
      
      p_PartonJets_Coll->push_back(temp_truthjet);
      
    }


  // Close objects
  // -------------
  
  delete jetDef;

  
}

// ============================================================================================================






// ============================================================================================================
void ANA_utils::Get_BottomQuarks(Pythia8::Event event, std::vector<TruthPart>* p_BQuarks_Coll)
//
// Find the bottom quarks that have come from a top, a Z, a W, or a Higgs. Fill a TruthPart object for each of
// the b-quarks that come from the above particles, just before it decays or hadronizes. Store it in the
// b-qaurk collection.
//  
//    Note: The b-quark must absolutely come from one of the aforementioned particles, but it should be the
//          the last in line before decaying or hadronizing, i.e. it must be after any gluon or photon
//          radiation, as well as not a carbon copy of the right b-quark. The strategy is therefore to find
//          the b-quarks immediately following the top/Z/W/H, and following their generations of descendants,
//          until there is no more b-quark in the next generation. A test showed that about 0.1% of the events
//          had up to the 13th generation, on a 1000 events test. The code therefore follow up to the 15th
//          generation of gluon/photon FSR or carbon copies, to make sure that there essentially no chances to
//          pick the wrong b-quarks. For most events, there are between 3 and 7 generations.   
//
// ============================================================================================================
{

  // Find the indices of the b-quarks directly coming from tops, W, Z, or H
  // ----------------------------------------------------------------------
  
  std::vector<int> bquark_uplevelindices;

  for (int i = 1; i < event.size(); i++)
    {
      if ((event[i].idAbs() == 5) &&
	  (event[event[i].mother1()].idAbs() == 6 || event[event[i].mother1()].idAbs() == 24 || event[event[i].mother1()].idAbs() == 23 || event[event[i].mother1()].idAbs() == 25) )
	{
	  bquark_uplevelindices.push_back(i);
	}
    }


  
  // Find the indices of the b-quarks just before they decay or hadronize
  // --------------------------------------------------------------------

  std::vector<int> bquark_indices;

  
        // Loop over all b-quarks directly coming from top/W/Z/H
        // .....................................................

            // Note: An integer will be used to keep the particle index of the last b-quark before decay or hadronization

  for (int i=0; i< bquark_uplevelindices.size(); i++)
    {
      
      
      int lastbquark_index=bquark_uplevelindices[i];


        // Get all the daughters of the considered b-quark and loop over them
        // ..................................................................

            // Note: Define an identifier to keep track of the generation that follows the last b-quark before decay and hadronization. That
            //       generatio is the first one that follows a top/W/Z/H decaying to b-quarks that do not contain any b-quark in it.

      std::vector<int> gen2 = event[lastbquark_index].daughterList();
      int topgen =1;
      if (gen2.size()>0) topgen = 2;
      for (int i2 = 0; i2<gen2.size();i2++)
	{

	// Check for b-quarks in the generation
	// ....................................
	  
	    // Note: If there are no b-quarks in this generation, there is nothing else to do because we already know which is the last b-quark index.
	    //       If there is a b-quark, move to the next generation by looking at its daughters.

	    // Note: This process is repeated over 15 generations, to make sure that we really get the last b-quark of the chain. The constraint on the
	    //       for loop (i<genX.size() )) guarantees that it is not possible to get further than this last generation.
	  
	  if (event[gen2[i2]].idAbs() != 5) continue;
	  else
	    {
	      lastbquark_index=gen2[i2];
	      std::vector<int> gen3 = event[lastbquark_index].daughterList();
	      if (gen3.size()>0) topgen = 3;
	      for (int i3 = 0; i3<gen3.size();i3++)
		{
		  
		  if (event[gen3[i3]].idAbs() !=5) continue;
		  else
		    {
		      lastbquark_index=gen3[i3];
		      std::vector<int> gen4 = event[lastbquark_index].daughterList();
		      if (gen4.size()>0) topgen = 4;
		      for (int i4 = 0; i4<gen4.size();i4++)
			{
			  if (event[gen4[i4]].idAbs() !=5) continue;
			  else
			    {
			      lastbquark_index=gen4[i4];
			      std::vector<int> gen5 = event[lastbquark_index].daughterList();
			      if (gen5.size()>0) topgen = 5;
			      for (int i5 = 0; i5<gen5.size();i5++)
				{
				  if (event[gen5[i5]].idAbs() !=5) continue;
				  else
				    {
				      lastbquark_index=gen5[i5];
				      std::vector<int> gen6 = event[lastbquark_index].daughterList();
				      if (gen6.size()>0) topgen = 6;
				      for (int i6 = 0; i6<gen6.size();i6++)
					{
					  if (event[gen6[i6]].idAbs() !=5) continue;
					  else
					    {
					      lastbquark_index=gen6[i6];
					      std::vector<int> gen7 = event[lastbquark_index].daughterList();
					      if (gen7.size()>0) topgen = 7;
					      for (int i7 = 0; i7<gen7.size();i7++)
						{
						  if (event[gen7[i7]].idAbs() !=5) continue;
						  else
						    {
						      lastbquark_index=gen7[i7];
						      std::vector<int> gen8 = event[lastbquark_index].daughterList();
						      if (gen8.size()>0) topgen = 8;
						      for (int i8 = 0; i8<gen8.size();i8++)
							{
							  if (event[gen8[i8]].idAbs() !=5) continue;
							  else
							    {
							      lastbquark_index=gen8[i8];
							      std::vector<int> gen9 = event[lastbquark_index].daughterList();
							      if (gen9.size()>0) topgen = 9;
							      for (int i9 = 0; i9<gen9.size();i9++)
								{
								  if (event[gen9[i9]].idAbs() !=5) continue;
								  else
								    {
								      lastbquark_index=gen9[i9];
								      std::vector<int> gen10 = event[lastbquark_index].daughterList();
								      if (gen10.size()>0) topgen = 10;
								      for (int i10 = 0; i10<gen10.size();i10++)
									{
									  if (event[gen10[i10]].idAbs() !=5) continue;
									  else
									    {
									      lastbquark_index=gen10[i10];
									      std::vector<int> gen11 = event[lastbquark_index].daughterList();
									      if (gen11.size()>0) topgen = 11;
									      for (int i11 = 0; i11<gen11.size();i11++)
										{
										  if (event[gen11[i11]].idAbs() !=5) continue;
										  else
										    {
										      lastbquark_index=gen11[i11];
										      std::vector<int> gen12 = event[lastbquark_index].daughterList();
										      if (gen12.size()>0) topgen = 12;
										      for (int i12 = 0; i12<gen12.size();i12++)
											{
											  if (event[gen12[i12]].idAbs() !=5) continue;
											  else
											    {
											      lastbquark_index=gen12[i12];
											      std::vector<int> gen13 = event[lastbquark_index].daughterList();
											      if (gen13.size()>0) topgen = 13;
											      for (int i13 = 0; i13<gen13.size();i13++)
												{
												  if (event[gen13[i13]].idAbs() !=5) continue;
												  else
												    {
												      lastbquark_index=gen13[i13];
												      std::vector<int> gen14 = event[lastbquark_index].daughterList();
												      if (gen14.size()>0) topgen = 14;
												      for (int i14 = 0; i14<gen14.size();i14++)
													{
													  if (event[gen14[i14]].idAbs() !=5) continue;
													  else
													    {
													      lastbquark_index=gen14[i14];
													      std::vector<int> gen15 = event[lastbquark_index].daughterList();
													      if (gen15.size()>0) topgen = 15;
													    }
													} // End Generation 14
												    }
												} // End Generation 13
											    }
											} // End Generation 12
										    }
										} // End Generation 11
									    }
									} // End Generation 10
								    }
								} // End Generation 9
							    }
							} // End Generation 8
						    }
						} // End Generation 7
					    }
					} // End Generation 6
				    }
				} // End Generation 5
			    }
			} // End Generation 4
		    }
		} // End Generation 3
	    }
	} // End Generation 2
      

        //  Push the last index in the vector of b-quarks indices for which the b-quark will now decay or hadronize
        //  .......................................................................................................
      
      bquark_indices.push_back(lastbquark_index);

    }


  
  // Fill a TruthPart object for each b-quark, and push it back into the b-quark collection
  // -------------------------------------------------------------------------------------
  
  for (auto i : bquark_indices)
    {
      
      TruthPart temp_b;
      TruthPart* p_temp_b = &temp_b;

      Fill_TruthPart(event, i, p_temp_b);
      p_BQuarks_Coll->push_back(temp_b);
    }


}

// ============================================================================================================






// ============================================================================================================
void ANA_utils::Get_BottomHadrons(Pythia8::Event event, std::vector<TruthPart> BQuark_Coll, std::vector<TruthPart>* p_BHadrons_Coll)
//
// Starting from a collection of b-quarks, this function returns a collection of B-hadrons that come from these
// b-quarks.  
//
// ============================================================================================================
{

  // Create a list of all B-hadrons PdgIds
  // -------------------------------------
  
  int BhadronsPdgId[88] = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 515, 525, 531, 10531, 533, 10533, 20553, 535, 541, 10541, 543, 10543, 20543, 545, 551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 20553, 30553, 100553, 110553, 120553, 130553, 200553, 210553, 220553, 300553, 9000553, 9010553, 555, 10555, 20555, 100555, 110555, 120555, 200555, 557, 100557, 5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 5312, 5322, 5314, 5324, 5332, 5334, 5142, 5242, 5412, 5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 5522, 5514, 5524, 5532, 5534, 5542, 5544, 5554};


  // Loop over all the bquarks coming from top/W/Z/H
  // -----------------------------------------------

  for (int i_bq=0; i_bq<BQuark_Coll.size(); i_bq++)
    {


  // Find the list of daughters of this b-quark
  // ------------------------------------------

      int bquark_index = (BQuark_Coll[i_bq]).Index();
      std::vector<int> bquark_daughterlist = event[bquark_index].daughterList();

      
  // Find the particle index of the B-hadrons following the b-quark hadronization
  // ----------------------------------------------------------------------------    

      // Loop over the daughters of the b-quark
      // ......................................

      for (int i_bd = 0; i_bd<bquark_daughterlist.size(); i_bd++)
	{
	  int temp_index = bquark_daughterlist[i_bd];

	  
      // Check if the particle is in the list of possible B-hadrons and if yes, push this particle back in the B-Hadron Collection
      // .........................................................................................................................

	  int temp_id = event[temp_index].id();

	  for (auto i_list : BhadronsPdgId)
	    {

	      if ( i_list == fabs(temp_id) )
		{
		  TruthPart temp_bhad;
		  TruthPart* p_temp_bhad = &temp_bhad;

		  Fill_TruthPart(event, temp_index, p_temp_bhad);
		  p_BHadrons_Coll->push_back(temp_bhad);
		}
	    }
	}
    }
}

// ============================================================================================================
