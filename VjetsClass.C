#define VjetsClass_cxx
#include "VjetsClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>
#include "TLorentzVector.h"

void VjetsClass::Loop()
{
  //   In a ROOT session, you can do:
  //      root> .L VjetsClass.C
  //      root> VjetsClass t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;
  Int_t nentries = Int_t(fChain->GetEntry());

  TFile *myfile = TFile::Open("histograms.root","recreate");
  //TH1F *lightjet_hist  = new TH1F("lightjet_hist", "Invariant Mass of the Most Energetic Light Jets", 100, 0, 37);
  TH1F *Wminus_neut  = new TH1F("Wminus_neut", "Transverse Mass of W- Boson, using neutrino data", 50, 0, 100);
  TH1F *Wplus_neut  = new TH1F("Wplus_neut", "Transverse Mass of W+ Boson, using neutrino data", 50, 0, 100);

  TH1F *Wminus_Met  = new TH1F("Wminus_Met", "Transverse Mass of W- Boson, using missing Et", 50, 0, 100);
  TH1F *Wplus_Met  = new TH1F("Wplus_Met", "Transverse Mass of W+ Boson, using missing Et", 50, 0, 100);

  // Iterate through the entries
  for (Int_t ientry=0; ientry<nentries;ientry++) {
    GetEntry(ientry);
    // Show(ientry);

    // First cuts: bjet and mEt
    bool bjet_satisfied = false;
    bool mEt_satisfied = false;

    for (size_t iBjet = 0; iBjet < nBjet; iBjet++) {
      if ((bjet_E->at(iBjet) >= 25.0) && (-2.5 < bjet_eta->at(iBjet)) && (bjet_eta->at(iBjet) < 2.5)) {
        bjet_satisfied = true;
        break;
      }
    }

    if (Met >= 25.0) {
      mEt_satisfied = true;
    }

    // Initial cut
    if (bjet_satisfied && mEt_satisfied) {
      // Semilepton checks
      bool electron_satisfied = false;
      bool lightjet_satisfied = false;

      int dilepton_nElect = 0;
      // Loop through electrons to check for at least 1 electron with 25 GeV, then exactly 2 electroncs with 25 GeV and certain nu
      for (size_t i = 0; i < nElectron; i++) {
        if (electron_E->at(i) >= 25.0) {
          electron_satisfied = true;
        }

        if ((electron_E->at(i) >= 25.0) && (-2.5 < electron_eta->at(i)) && (electron_eta->at(i) < 2.5)) {
          dilepton_nElect += 1;
        }
      }

      int semilepton_nLightjet = 0;
      // Loop through light jets
      for (size_t i = 0; i < nLightjet; i++) {
        if (lightjet_E->at(i) >= 30.0) {
          semilepton_nLightjet += 1;
          if (semilepton_nLightjet == 2) {
            lightjet_satisfied = true;
            break;
          }
        }
      }

      // Semilepton channel
      if (electron_satisfied && lightjet_satisfied) {/*

        // Find the two most energetic light jets
        std::vector<int> energetic_check; // vector of indices of light jets to checks

        // Check only above 30 GeV
        for (size_t i = 0; i < nLightjet; i++) {
          if (lightjet_E->at(i) >= 30.0) {
            energetic_check.push_back(i);
          }
        }

        // Compare energies
        int most_energetic = 0;
        int secondmost_energetic;
        for (auto i : energetic_check) {
          if (i == 0) {
            continue;
          }
          if (i == 1) {
            if (lightjet_E->at(0) <= lightjet_E->at(1)) {
              most_energetic = 1;
              secondmost_energetic = 0;
            }
            else if (lightjet_E->at(0) > lightjet_E->at(1)) {
              secondmost_energetic = 1;
            }
          }

          if (lightjet_E->at(i) > lightjet_E->at(secondmost_energetic)) {
            secondmost_energetic = i;
            if (lightjet_E->at(i) > lightjet_E->at(most_energetic)) {
              secondmost_energetic = most_energetic;
              most_energetic = i;
            }
          }
        }

        // Create Lorentz vectors to calculate invariant mass of highest energy light jets
        TLorentzVector lightjet1, lightjet2;
        lightjet1.SetPtEtaPhiE(lightjet_pt->at(most_energetic), lightjet_eta->at(most_energetic), lightjet_phi->at(most_energetic), lightjet_E->at(most_energetic));
        lightjet2.SetPtEtaPhiE(lightjet_pt->at(secondmost_energetic), lightjet_eta->at(secondmost_energetic), lightjet_phi->at(secondmost_energetic), lightjet_E->at(secondmost_energetic));

        lightjet_hist->Fill(lightjet1.M());
        lightjet_hist->Fill(lightjet2.M());

        energetic_check.clear(); */
      }

      // Dilepton channel
      if (dilepton_nElect == 2) {
        for (size_t i = 0; i < nElectron; i++) {
          TLorentzVector elect;
          elect.SetPtEtaPhiE(electron_pt->at(i), electron_eta->at(i), electron_phi->at(i), electron_E->at(i));

          for (size_t j = 0; j < nNeutrino; j++) {
            TLorentzVector neut;
            neut.SetPtEtaPhiE(neutrino_pt->at(j), neutrino_eta->at(j), neutrino_phi->at(j), neutrino_E->at(j));

            // W- boson
            if (electron_charge->at(i) < 0 && neutrino_PdgId->at(j) < 0) {
              Wminus_neut->Fill(sqrt(2.0*(elect.Px() + elect.Py())*(neut.Px() + neut.Py())*(1 - cos(elect.Phi() - neut.Phi()))));
              Wminus_Met->Fill(sqrt(2.0*(elect.Px() + elect.Py())*(Met/2.0)*(1 - cos(elect.Phi() - Met_phi))));
              /*
              std::cout << "Missing energy: " << Met << '\n';
              std::cout << "Transverse according to Met: " << sqrt(abs(2.0*(elect.Px() + elect.Py())*(Met/2.0)*(1 - cos(elect.Phi() - Met_phi)))) << '\n';
              std::cout << "Transverse according to Neutrino: " << sqrt(abs(2.0*(elect.Px() + elect.Py())*(neut.Px() + neut.Py())*(1 - cos(elect.Phi() - neut.Phi())))) << '\n';
              */
            }

            // W+ boson
            if (electron_charge->at(i) > 0 && neutrino_PdgId->at(j) > 0) {
              Wplus_neut->Fill(sqrt(2.0*(elect.Px() + elect.Py())*(neut.Px() + neut.Py())*(1 - cos(elect.Phi() - neut.Phi()))));
              Wplus_Met->Fill(sqrt(2.0*(elect.Px() + elect.Py())*(Met/2.0)*(1 - cos(elect.Phi() - Met_phi))));
            }
          }
        }
        /*
        for (size_t i = 0; i < nMuon; i++) {
          TLorentzVector muon;
          muon.SetPtEtaPhiE(muon_pt->at(i), muon_eta->at(i), muon_phi->at(i), muon_E->at(i));
          lepton_hist->Fill(sqrt(pow(muon_E->at(i), 2) - pow(muon.Pz(),2)));
        }
        */
      }
    }
  }

  //lightjet_hist->Write();
  Wminus_neut->Write();
  Wplus_neut->Write();
  Wminus_Met->Write();
  Wplus_Met->Write();
}
