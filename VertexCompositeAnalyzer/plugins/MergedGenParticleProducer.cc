#include "MergedGenParticleProducer.hh"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "HepPDT/ParticleID.hh"

MergedGenParticleProducerRice::MergedGenParticleProducerRice(const edm::ParameterSet& config)
{
  input_pruned_ = consumes<edm::View<reco::GenParticle>>(config.getParameter<edm::InputTag>("inputPruned"));
  input_packed_ = consumes<edm::View<pat::PackedGenParticle>>(config.getParameter<edm::InputTag>("inputPacked"));

  produces<reco::GenParticleCollection>();
}

void MergedGenParticleProducerRice::produce(edm::Event& event, const edm::EventSetup& setup)
{
  // Need a ref to the product now for creating the mother/daughter refs
  auto ref = event.getRefBeforePut<reco::GenParticleCollection>();

  // Get the input collections
  edm::Handle<edm::View<reco::GenParticle> > pruned_handle;
  event.getByToken(input_pruned_, pruned_handle);

  edm::Handle<edm::View<pat::PackedGenParticle> > packed_handle;
  event.getByToken(input_packed_, packed_handle);

  size_t pruned_handle_size = pruned_handle.isValid() ? pruned_handle->size() : 0;

  // First determine which packed particles are also still in the pruned collection
  // so that we can skip them later
  std::map<pat::PackedGenParticle const*, reco::GenParticle const*> st1_dup_map;

  // Also map pointers in the original pruned collection to their index in the vector.
  // This index will be the same in the merged collection.
  std::map<reco::Candidate const*, std::size_t> pruned_idx_map;

  for (unsigned int i = 0; i < pruned_handle_size; ++i) {
    reco::GenParticle const& src = pruned_handle->at(i);
    pruned_idx_map[&src] = i;
    if (src.status() != 1) continue;

    // Convert the pruned GenParticle into a PackedGenParticle then do an exact
    // floating point comparison of the pt, eta and phi
    // This of course relies on the PackedGenParticle constructor being identical
    // between the CMSSW version this sample was produced with and the one we're
    // analysing with
    pat::PackedGenParticle pks(src, reco::GenParticleRef());
    unsigned found_matches = 0;
    for (unsigned j = 0; j < packed_handle->size(); ++j) {
      pat::PackedGenParticle const& pk = packed_handle->at(j);
      if ( pks.pdgId() != pk.pdgId() or pks.p4() != pk.p4() ) continue;
      ++found_matches;
      st1_dup_map[&pk] = &src;
    }
    if (found_matches > 1) {
      edm::LogWarning("MergedGenParticleProducerRice") << "Found multiple packed matches for: " << i << "\t" << src.pdgId() << "\t" << src.pt() << "\t" << src.y() << "\n";
    }
    else if (found_matches == 0 && std::abs(src.y()) < 6.0) {
      edm::LogWarning("MergedGenParticleProducerRice") << "unmatched status 1: " << i << "\t" << src.pdgId() << "\t" << src.pt() << "\t" << src.y() << "\n";
    }
  }

  // Fix by Markus
  // check for photons from pruned (light) hadrons
  unsigned int nPhotonsFromPrunedHadron = 0;
  for (unsigned int j = 0; j < packed_handle->size(); ++j) {
    pat::PackedGenParticle const& pk = packed_handle->at(j);
    if (isPhotonFromPrunedHadron(pk)) ++nPhotonsFromPrunedHadron;
  }

  // At this point we know what the size of the merged GenParticle will be so we can create it
  const unsigned int n = pruned_handle_size + (packed_handle->size() - st1_dup_map.size()) + nPhotonsFromPrunedHadron;
  auto cands = std::unique_ptr<reco::GenParticleCollection>(new reco::GenParticleCollection(n));

  // First copy in all the pruned candidates
  for (unsigned i = 0; i < pruned_handle_size; ++i) {
    reco::GenParticle const& old_cand = pruned_handle->at(i);
    reco::GenParticle & new_cand = cands->at(i);
    new_cand = reco::GenParticle(pruned_handle->at(i));
    // Update the mother and daughter refs to this new merged collection
    new_cand.resetMothers(ref.id());
    new_cand.resetDaughters(ref.id());
    for (unsigned m = 0; m < old_cand.numberOfMothers(); ++m) {
      new_cand.addMother(reco::GenParticleRef(ref, pruned_idx_map.at(old_cand.mother(m))));
    }
    for (unsigned d = 0; d < old_cand.numberOfDaughters(); ++d) {
      new_cand.addDaughter(reco::GenParticleRef(ref, pruned_idx_map.at(old_cand.daughter(d))));
    }
  }

  // Now copy in the packed candidates that are not already in the pruned
  for (unsigned i = 0, idx = pruned_handle_size; i < packed_handle->size(); ++i) {
    pat::PackedGenParticle const& pk = packed_handle->at(i);
    if (st1_dup_map.count(&pk)) continue;
    reco::GenParticle & new_cand = cands->at(idx);
    new_cand = reco::GenParticle(pk.charge(), pk.p4(), pk.vertex(), pk.pdgId(), 1, true);
    reco::GenStatusFlags& new_status = new_cand.statusFlags();
    new_status = pk.statusFlags();

    // Insert dummy pi0 mothers for orphaned photons
    if (isPhotonFromPrunedHadron(pk)) {
      ++idx;
      reco::GenParticle & dummy_mother = cands->at(idx);
      dummy_mother = reco::GenParticle(0, pk.p4(), pk.vertex(), 111, 2, true);
      for (unsigned m = 0; m < pk.numberOfMothers(); ++m) {
        new_cand.addMother(reco::GenParticleRef(ref, idx));
        // Since the packed candidates drop the vertex position we'll take this from the mother
        if (m == 0) {
          dummy_mother.setP4(pk.mother(m)->p4());
          dummy_mother.setVertex(pk.mother(m)->vertex());
          new_cand.setVertex(pk.mother(m)->vertex());
        }
        // Should then add the GenParticle as a daughter of its dummy mother
        dummy_mother.addDaughter(reco::GenParticleRef(ref, idx-1));
      }
    }
    // Connect to mother from pruned particles
    reco::GenParticle & daughter = cands->at(idx);
    for (unsigned m = 0; m < pk.numberOfMothers(); ++m) {
      if (!pruned_handle.isValid()) break;
      daughter.addMother(reco::GenParticleRef(ref, pruned_idx_map.at(pk.mother(m))));
      // Since the packed candidates drop the vertex position we'll take this from the mother
      if (m == 0) {
        daughter.setVertex(pk.mother(m)->vertex());
      }
      // Should then add this GenParticle as a daughter of its mother
      cands->at(pruned_idx_map.at(pk.mother(m))).addDaughter(reco::GenParticleRef(ref, idx));
    }
    ++idx;
  }

  event.put(std::move(cands));
}

bool MergedGenParticleProducerRice::isPhotonFromPrunedHadron(const pat::PackedGenParticle& pk) const
{
  if (pk.pdgId() == 22 and pk.statusFlags().isDirectHadronDecayProduct()) {
    // no mother
    if (pk.numberOfMothers() == 0) return true;
    // miniaod mother not compatible with the status flag
    HepPDT::ParticleID motherid(pk.mother(0)->pdgId());
    if (not (motherid.isHadron() and pk.mother(0)->status() == 2)) return true;
  }
  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MergedGenParticleProducerRice);

