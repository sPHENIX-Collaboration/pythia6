/**
 \file factory.cpp
 
 \author Thomas Burton 
 \date 10/9/12
 \copyright 2012 BNL.
 */

#include "factory.h"

#include <TBranch.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TVector3.h>

#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/erhic/Kinematics.h>

#include "pythia_commons.h"
#include "pythia_erhic.h"

// =====================================================================
// =====================================================================
Factory::Factory()
: mEvent(new erhic::EventPythia()) {
}

// =====================================================================
// =====================================================================
Factory::~Factory() {
   if(mEvent) {
      delete mEvent;
      mEvent = NULL;
   } // if
}

// =====================================================================
// Populate the stored event with the current contents of
// PYTHIA and return a pointer to the event.
// Do not delete this pointer!
// =====================================================================
erhic::EventPythia* Factory::Create() {
   mEvent->Clear("");
   // See the PYTHIA manual for more details about the meaning
   // of the various msti(), pari(), vint() variables.
   mEvent->SetGenEvent(__pythia6_MOD_genevent);
   mEvent->SetProcess(msti(1));
   mEvent->SetNucleon(msti(12));
   mEvent->SetTargetParton(msti(16));
   mEvent->SetBeamParton(msti(15));
   mEvent->SetTargetPartonX(pari(34));
   mEvent->SetBeamPartonX(pari(33));
   mEvent->SetBeamPartonTheta(pari(53));
   mEvent->SetLeptonPhi(vint(313));
   mEvent->SetF1(__pythia6_MOD_f1);
   mEvent->SetF2(__pythia6_MOD_f2);
   mEvent->SetSigmaRad(__pythia6_MOD_sigma_rad);
   mEvent->SetHardS(pari(14));
   mEvent->SetHardT(pari(15));
   mEvent->SetHardU(pari(16));
   mEvent->SetHardQ2(pari(22));
   mEvent->SetHardPt2(pari(18));
   mEvent->SetSigRadCor(__pythia6_MOD_sigradcor);
   mEvent->SetEBrems(__pythia6_MOD_ebrems);
   mEvent->SetPhotonFlux(vint(319));
   mEvent->SetTrueY(vint(309));
   mEvent->SetTrueQ2(vint(307));
   mEvent->SetTrueX(__pythia6_MOD_truex);
   mEvent->SetTrueW2(__pythia6_MOD_truew2);
   mEvent->SetTrueNu(__pythia6_MOD_truenu);
   mEvent->SetR(__pythia6_MOD_r);
   static std::stringstream stream;
   const int ntracks = pyjets_.n;
   // Loop with indices from [1, N] as per Fortran convention.
   // Accessor functions k(), p(), v() take care of correct indexing
   // of the equivalent C++ arrays.
   for(int i(1); i <= ntracks; ++i) {
      erhic::ParticleMC* track = new erhic::ParticleMC;
      track->SetIndex(i);
      track->SetStatus(k(i, 1));
      track->SetId(k(i, 2));
      track->SetParentIndex(k(i, 3));
      track->SetChild1Index(k(i, 4));
      track->SetChildNIndex(k(i, 5));
      track->Set4Vector(TLorentzVector(p(i, 1), p(i, 2),
                                       p(i, 3), p(i, 4)));
      track->SetM(p(i, 5));
      track->SetVertex(TVector3(v(i, 1), v(i, 2), v(i, 3)));
      track->SetEvent(mEvent);
      mEvent->AddLast(track);
   } // for
   // If a radiative photon was generated append it to the track record.
   if(__pythia6_MOD_ebrems > 0.) {
      erhic::ParticleMC* track = new erhic::ParticleMC;
      track->SetIndex(mEvent->GetNTracks() + 1);
      track->SetStatus(55);
      track->SetId(22);
      track->SetParentIndex(1);
      track->SetChild1Index(0);
      track->SetChildNIndex(0);
      track->Set4Vector(TLorentzVector(__pythia6_MOD_radgamp[0],
                                       __pythia6_MOD_radgamp[1],
                                       __pythia6_MOD_radgamp[2],
                                       __pythia6_MOD_radgame));
      track->SetVertex(TVector3(0., 0., 0.));
      track->SetEvent(mEvent);
      mEvent->AddLast(track);
   } // if
   // After adding all tracks, compute event kinematics.
   // Compute using all three methods: electron, Jacquet-Blondel
   // and double angle.
   erhic::DisKinematics* nm =
   erhic::LeptonKinematicsComputer(*mEvent).Calculate();
   erhic::DisKinematics* jb =
   erhic::JacquetBlondelComputer(*mEvent).Calculate();
   erhic::DisKinematics* da =
   erhic::DoubleAngleComputer(*mEvent).Calculate();
   if(nm) {
      mEvent->SetLeptonKinematics(*nm);
      delete nm;
   } // if
   if(jb) {
      mEvent->SetJacquetBlondelKinematics(*jb);
      delete jb;
   } // if
   if(da) {
      mEvent->SetDoubleAngleKinematics(*da);
      delete da;
   } // if
   // Set particle properties that depend on other particles.
   // This has to come last as some of these properties depend
   // on the computed event kinematics.
   for(unsigned i(0); i < mEvent->GetNTracks(); ++i) {
      mEvent->GetTrack(i)->ComputeEventDependentQuantities(*mEvent);
   } // for
   return mEvent;
}

// =====================================================================
// Returns a string with the full (including namespace) class name
// of the event type produced.
// This is important for use with ROOT TTree to ensure the correct
// event type in branches.
// =====================================================================
std::string Factory::EventName() const {
   return mEvent->Class()->GetName();
}

// =====================================================================
// Add a branch named "name" for the event type generated
// by this factory to a ROOT TTree.
// Returns a pointer to the branch, or NULL in the case of an error.
// =====================================================================
TBranch* Factory::Branch(TTree& tree, const std::string& name) {
   return tree.Branch(name.c_str(), &mEvent);
}
