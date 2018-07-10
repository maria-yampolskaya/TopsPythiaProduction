//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 10 09:39:17 2018 by ROOT version 6.15/01
// from TTree ParticleTree/Particle Data
// found on file: outfile.root
//////////////////////////////////////////////////////////

#ifndef VjetsClass_h
#define VjetsClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class VjetsClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *top_pt;
   vector<float>   *top_eta;
   vector<float>   *top_phi;
   vector<float>   *top_E;
   vector<float>   *neutrino_pt;
   vector<float>   *neutrino_eta;
   vector<float>   *neutrino_phi;
   vector<float>   *neutrino_E;
   vector<int>     *neutrino_PdgId;
   vector<float>   *muon_pt;
   vector<float>   *muon_eta;
   vector<float>   *muon_phi;
   vector<float>   *muon_E;
   vector<float>   *muon_charge;
   vector<float>   *electron_pt;
   vector<float>   *electron_eta;
   vector<float>   *electron_phi;
   vector<float>   *electron_E;
   vector<float>   *electron_charge;
   vector<float>   *boson_pt;
   vector<float>   *boson_eta;
   vector<float>   *boson_phi;
   vector<float>   *boson_E;
   vector<int>     *boson_ID;
   vector<float>   *lightjet_pt;
   vector<float>   *lightjet_eta;
   vector<float>   *lightjet_phi;
   vector<float>   *lightjet_E;
   vector<int>     *lightjet_nPart;
   vector<float>   *bjet_pt;
   vector<float>   *bjet_eta;
   vector<float>   *bjet_phi;
   vector<float>   *bjet_E;
   vector<int>     *bjet_nPart;
   vector<float>   *lightpartonjet_pt;
   vector<float>   *lightpartonjet_eta;
   vector<float>   *lightpartonjet_phi;
   vector<float>   *lightpartonjet_E;
   vector<float>   *bpartonjet_pt;
   vector<float>   *bpartonjet_eta;
   vector<float>   *bpartonjet_phi;
   vector<float>   *bpartonjet_E;
   Int_t           nTop;
   Int_t           nNeutrino;
   Int_t           nMuon;
   Int_t           nElectron;
   Int_t           nLightjet;
   Int_t           nBjet;
   Int_t           nLightpartonjet;
   Int_t           nBpartonjet;
   Int_t           nBoson;
   Float_t         Met;
   Float_t         Met_phi;

   // List of branches
   TBranch        *b_top_pt;   //!
   TBranch        *b_top_eta;   //!
   TBranch        *b_top_phi;   //!
   TBranch        *b_top_E;   //!
   TBranch        *b_neutrino_pt;   //!
   TBranch        *b_neutrino_eta;   //!
   TBranch        *b_neutrino_phi;   //!
   TBranch        *b_neutrino_E;   //!
   TBranch        *b_neutrino_PdgId;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_E;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_electron_pt;   //!
   TBranch        *b_electron_eta;   //!
   TBranch        *b_electron_phi;   //!
   TBranch        *b_electron_E;   //!
   TBranch        *b_electron_charge;   //!
   TBranch        *b_boson_pt;   //!
   TBranch        *b_boson_eta;   //!
   TBranch        *b_boson_phi;   //!
   TBranch        *b_boson_E;   //!
   TBranch        *b_boson_ID;   //!
   TBranch        *b_lightjet_pt;   //!
   TBranch        *b_lightjet_eta;   //!
   TBranch        *b_lightjet_phi;   //!
   TBranch        *b_lightjet_E;   //!
   TBranch        *b_lightjet_nPart;   //!
   TBranch        *b_bjet_pt;   //!
   TBranch        *b_bjet_eta;   //!
   TBranch        *b_bjet_phi;   //!
   TBranch        *b_bjet_E;   //!
   TBranch        *b_bjet_nPart;   //!
   TBranch        *b_lightpartonjet_pt;   //!
   TBranch        *b_lightpartonjet_eta;   //!
   TBranch        *b_lightpartonjet_phi;   //!
   TBranch        *b_lightpartonjet_E;   //!
   TBranch        *b_bpartonjet_pt;   //!
   TBranch        *b_bpartonjet_eta;   //!
   TBranch        *b_bpartonjet_phi;   //!
   TBranch        *b_bpartonjet_E;   //!
   TBranch        *b_nTop;   //!
   TBranch        *b_nNeutrino;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_nLightjet;   //!
   TBranch        *b_nBjet;   //!
   TBranch        *b_nLightpartonjet;   //!
   TBranch        *b_nBpartonjet;   //!
   TBranch        *b_nBoson;   //!
   TBranch        *b_Met;   //!
   TBranch        *b_Met_phi;   //!

   VjetsClass(TTree *tree=0);
   virtual ~VjetsClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef VjetsClass_cxx
VjetsClass::VjetsClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outfile.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outfile.root");
      }
      f->GetObject("ParticleTree",tree);

   }
   Init(tree);
}

VjetsClass::~VjetsClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t VjetsClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t VjetsClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void VjetsClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   top_pt = 0;
   top_eta = 0;
   top_phi = 0;
   top_E = 0;
   neutrino_pt = 0;
   neutrino_eta = 0;
   neutrino_phi = 0;
   neutrino_E = 0;
   neutrino_PdgId = 0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_E = 0;
   muon_charge = 0;
   electron_pt = 0;
   electron_eta = 0;
   electron_phi = 0;
   electron_E = 0;
   electron_charge = 0;
   boson_pt = 0;
   boson_eta = 0;
   boson_phi = 0;
   boson_E = 0;
   boson_ID = 0;
   lightjet_pt = 0;
   lightjet_eta = 0;
   lightjet_phi = 0;
   lightjet_E = 0;
   lightjet_nPart = 0;
   bjet_pt = 0;
   bjet_eta = 0;
   bjet_phi = 0;
   bjet_E = 0;
   bjet_nPart = 0;
   lightpartonjet_pt = 0;
   lightpartonjet_eta = 0;
   lightpartonjet_phi = 0;
   lightpartonjet_E = 0;
   bpartonjet_pt = 0;
   bpartonjet_eta = 0;
   bpartonjet_phi = 0;
   bpartonjet_E = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("top_pt", &top_pt, &b_top_pt);
   fChain->SetBranchAddress("top_eta", &top_eta, &b_top_eta);
   fChain->SetBranchAddress("top_phi", &top_phi, &b_top_phi);
   fChain->SetBranchAddress("top_E", &top_E, &b_top_E);
   fChain->SetBranchAddress("neutrino_pt", &neutrino_pt, &b_neutrino_pt);
   fChain->SetBranchAddress("neutrino_eta", &neutrino_eta, &b_neutrino_eta);
   fChain->SetBranchAddress("neutrino_phi", &neutrino_phi, &b_neutrino_phi);
   fChain->SetBranchAddress("neutrino_E", &neutrino_E, &b_neutrino_E);
   fChain->SetBranchAddress("neutrino_PdgId", &neutrino_PdgId, &b_neutrino_PdgId);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_E", &muon_E, &b_muon_E);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("electron_pt", &electron_pt, &b_electron_pt);
   fChain->SetBranchAddress("electron_eta", &electron_eta, &b_electron_eta);
   fChain->SetBranchAddress("electron_phi", &electron_phi, &b_electron_phi);
   fChain->SetBranchAddress("electron_E", &electron_E, &b_electron_E);
   fChain->SetBranchAddress("electron_charge", &electron_charge, &b_electron_charge);
   fChain->SetBranchAddress("boson_pt", &boson_pt, &b_boson_pt);
   fChain->SetBranchAddress("boson_eta", &boson_eta, &b_boson_eta);
   fChain->SetBranchAddress("boson_phi", &boson_phi, &b_boson_phi);
   fChain->SetBranchAddress("boson_E", &boson_E, &b_boson_E);
   fChain->SetBranchAddress("boson_ID", &boson_ID, &b_boson_ID);
   fChain->SetBranchAddress("lightjet_pt", &lightjet_pt, &b_lightjet_pt);
   fChain->SetBranchAddress("lightjet_eta", &lightjet_eta, &b_lightjet_eta);
   fChain->SetBranchAddress("lightjet_phi", &lightjet_phi, &b_lightjet_phi);
   fChain->SetBranchAddress("lightjet_E", &lightjet_E, &b_lightjet_E);
   fChain->SetBranchAddress("lightjet_nPart", &lightjet_nPart, &b_lightjet_nPart);
   fChain->SetBranchAddress("bjet_pt", &bjet_pt, &b_bjet_pt);
   fChain->SetBranchAddress("bjet_eta", &bjet_eta, &b_bjet_eta);
   fChain->SetBranchAddress("bjet_phi", &bjet_phi, &b_bjet_phi);
   fChain->SetBranchAddress("bjet_E", &bjet_E, &b_bjet_E);
   fChain->SetBranchAddress("bjet_nPart", &bjet_nPart, &b_bjet_nPart);
   fChain->SetBranchAddress("lightpartonjet_pt", &lightpartonjet_pt, &b_lightpartonjet_pt);
   fChain->SetBranchAddress("lightpartonjet_eta", &lightpartonjet_eta, &b_lightpartonjet_eta);
   fChain->SetBranchAddress("lightpartonjet_phi", &lightpartonjet_phi, &b_lightpartonjet_phi);
   fChain->SetBranchAddress("lightpartonjet_E", &lightpartonjet_E, &b_lightpartonjet_E);
   fChain->SetBranchAddress("bpartonjet_pt", &bpartonjet_pt, &b_bpartonjet_pt);
   fChain->SetBranchAddress("bpartonjet_eta", &bpartonjet_eta, &b_bpartonjet_eta);
   fChain->SetBranchAddress("bpartonjet_phi", &bpartonjet_phi, &b_bpartonjet_phi);
   fChain->SetBranchAddress("bpartonjet_E", &bpartonjet_E, &b_bpartonjet_E);
   fChain->SetBranchAddress("nTop", &nTop, &b_nTop);
   fChain->SetBranchAddress("nNeutrino", &nNeutrino, &b_nNeutrino);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("nLightjet", &nLightjet, &b_nLightjet);
   fChain->SetBranchAddress("nBjet", &nBjet, &b_nBjet);
   fChain->SetBranchAddress("nLightpartonjet", &nLightpartonjet, &b_nLightpartonjet);
   fChain->SetBranchAddress("nBpartonjet", &nBpartonjet, &b_nBpartonjet);
   fChain->SetBranchAddress("nBoson", &nBoson, &b_nBoson);
   fChain->SetBranchAddress("Met", &Met, &b_Met);
   fChain->SetBranchAddress("Met_phi", &Met_phi, &b_Met_phi);
   Notify();
}

Bool_t VjetsClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void VjetsClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t VjetsClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef VjetsClass_cxx
