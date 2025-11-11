#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TEfficiency.h"



class efficiencyMC : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit efficiencyMC(const edm::ParameterSet&);
      ~efficiencyMC();

      edm::ConsumesCollector iC = consumesCollector();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::ParameterSet parameters;

      //
      // --- Tokens and Handles
      //


      // GenParticles
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken;
      edm::Handle<edm::View<reco::GenParticle> > gens;

      // updated muons
      edm::EDGetTokenT<edm::View<Run3ScoutingMuon> > muonsToken;
      edm::Handle<edm::View<Run3ScoutingMuon> > muons;

      edm::EDGetTokenT<edm::View<Run3ScoutingVertex> > svsToken;
      edm::Handle<edm::View<Run3ScoutingVertex> > svs;

      edm::EDGetTokenT<edm::View<Run3ScoutingVertex> > primaryVerticesToken;
      edm::Handle<edm::View<Run3ScoutingVertex> > primaryVertices;

      //
      // --- Variables
      //

      bool isData = false;

      // --- HLT
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      triggerExpression::Data triggerCache_;

      //
      // --- Output
      //
      std::string output_filename;

      TH1F* histo_pt_gen;
      TH1F* histo_eta_gen;
      TH1F* histo_phi_gen;
      TH1F* histo_lxy_gen;
      TH2F* histo_pt_vs_lxy_gen;
      
      TH1F* histo_pt_gen_matched;
      TH1F* histo_eta_gen_matched;
      TH1F* histo_phi_gen_matched;
      TH1F* histo_lxy_gen_matched;

      TH2F* histo_pt_comp_matched;

      TH1F* histo_pt_scouting;
      TH1F* histo_eta_scouting;
      TH1F* histo_phi_scouting;
      TH1F* histo_lxy_scouting;

      TH1F* histo_pt_scouting_matched;
      TH1F* histo_eta_scouting_matched;
      TH1F* histo_phi_scouting_matched;
      TH1F* histo_lxy_scouting_matched;

      TProfile* histo_pt_resolution;
      TProfile* histo_pt_resolution_vs_abseta;
      TProfile* histo_pt_resolution_vs_lxy;      

      TEfficiency *reco_efficiency_pt_scouting;
      TEfficiency *reco_efficiency_lxy_scouting;

      TEfficiency *efficiency_trg_DoubleMuon;
      TEfficiency *efficiency_minpt_DoubleMuon;
      TEfficiency *efficiency_minlxy_DoubleMuon;
      TEfficiency *efficiency_maxlxy_DoubleMuon;
      TEfficiency *efficiency_npv_DoubleMuon;

      TH1F *histo_npv;
      TH1F *event_counts;
      TFile *file_out;


};

// Constructor
efficiencyMC::efficiencyMC(const edm::ParameterSet& iConfig) :
   triggerCache_(triggerExpression::Data(iConfig.getParameterSet("triggerConfiguration"), consumesCollector()))
{

   usesResource("TFileService");

   parameters = iConfig;

   event_counts = new TH1F("event_counts", "", 1, 0, 1);

   histo_pt_gen = new TH1F("histo_pt_gen", ";Generated muon p_{T} (GeV); Number of muons", 100, 0, 50);
   histo_eta_gen = new TH1F("histo_eta_gen", ";Generated muon #eta; Number of muons", 30, -4, 4);
   histo_phi_gen = new TH1F("histo_phi_gen", ";Generated muon #phi (rad); Number of muons", 30, -3.2, 3.2);
   histo_lxy_gen = new TH1F("histo_lxy_gen", ";Generated muon l_{xy} (cm); Number of muons", 30, 0.0, 30.0);
   histo_pt_gen_matched = new TH1F("histo_pt_gen_matched", ";Generated muon p_{T} (GeV); Number of muons", 100, 0, 50);
   histo_eta_gen_matched = new TH1F("histo_eta_gen_matched", ";Generated muon #eta; Number of muons", 30, -4, 4);
   histo_phi_gen_matched = new TH1F("histo_phi_gen_matched", ";Generated muon #phi (rad); Number of muons", 30, -3.2, 3.2);
   histo_lxy_gen_matched = new TH1F("histo_lxy_gen_matched", ";Generated muon l_{xy} (cm); Number of muons", 30, 0.0, 30.0);

   histo_pt_vs_lxy_gen = new TH2F("histo_pt_vs_lxy_gen", ";Generated muon p_{T} (GeV); Generated muon l_{xy} (cm)", 100, 0, 50, 80, 0, 80);

   histo_pt_comp_matched = new TH2F("histo_pt_comp_matched", ";Generated muon p_{T} (GeV); Reconstructed muon p_{T} (GeV)", 100, 0, 100, 100, 0, 100);

   histo_pt_scouting = new TH1F("histo_pt_scouting", ";Scouting muon p_{T} (GeV); Number of muons", 100, 0, 50);
   histo_eta_scouting = new TH1F("histo_eta_scouting", ";Scouting muon #eta; Number of muons", 30, -4, 4);
   histo_phi_scouting = new TH1F("histo_phi_scouting", ";Scouting muon #phi (rad); Number of muons", 30, -3.2, 3.2);
   histo_lxy_scouting = new TH1F("histo_lxy_scouting", ";Scouting SV (cm); Number of vertices", 80, 0, 80);

   histo_pt_scouting_matched = new TH1F("histo_pt_scouting_matched", ";Scouting muon p_{T} (GeV); Number of muons", 100, 0, 50);
   histo_eta_scouting_matched = new TH1F("histo_eta_scouting_matched", ";Scouting muon #eta; Number of muons", 30, -4, 4);
   histo_phi_scouting_matched = new TH1F("histo_phi_scouting_matched", ";Scouting muon #phi (rad); Number of muons", 30, -3.2, 3.2);
   histo_lxy_scouting_matched = new TH1F("histo_lxy_scouting_matched", ";Scouting SV (cm); Number of vertices", 80, 0, 80);

   histo_pt_resolution = new TProfile("histo_pt_resolution",
					";Generated muon p_{T} (GeV); (p_{T}^{scouting} - p_{T}^{gen}) / p_{T}^{gen}",
					100, 0, 50);

   histo_pt_resolution_vs_abseta = new TProfile("histo_pt_resolution_vs_abseta",
					";Generated muon |#eta|; (p_{T}^{scouting} - p_{T}^{gen}) / p_{T}^{gen}",
					24, 0, 2.4);

   histo_pt_resolution_vs_lxy = new TProfile("histo_pt_resolution_vs_lxy",
					";Generated muon l_{xy} (cm); (p_{T}^{scouting} - p_{T}^{gen}) / p_{T}^{gen}",
					80, 0, 80);

   
   reco_efficiency_pt_scouting = new TEfficiency("reco_efficiency_pt_scouting", ";Generated p_{T} (GeV); HLT efficiency", 100, 0, 50);
   reco_efficiency_lxy_scouting = new TEfficiency("reco_efficiency_lxy_scouting", ";Generated l_{xy} (cm); HLT efficiency", 80, 0, 80);

   efficiency_trg_DoubleMuon = new TEfficiency("efficiency_trg_DoubleMuon", "DoubleMuon efficiency; Leading p_{T} (GeV) ; Subleading p_{T} (GeV)", 20, 0, 20, 20, 0, 20);
   efficiency_minpt_DoubleMuon = new TEfficiency("efficiency_minpt_DoubleMuon", ";Subleading p_{T} (GeV); HLT efficiency", 100, 0, 50);
   efficiency_minlxy_DoubleMuon = new TEfficiency("efficiency_minlxy_DoubleMuon", ";Generated min l_{xy} (cm); HLT efficiency", 80, 0, 80);
   efficiency_maxlxy_DoubleMuon = new TEfficiency("efficiency_maxlxy_DoubleMuon", ";Generated max l_{xy} (cm); HLT efficiency", 80, 0, 80);
   efficiency_npv_DoubleMuon = new TEfficiency("efficiency_npv_DoubleMuon", ";Number of Primary Vertices; HLT efficiency", 50, 0, 50);

   histo_npv = new TH1F("histo_npv", ";Number of Primary Vertices; Number of events", 50, 0, 50);

   isData = parameters.getParameter<bool>("isData");

   genToken = consumes<edm::View<reco::GenParticle> >  (parameters.getParameter<edm::InputTag>("generatedParticles"));
   muonsToken = consumes<edm::View<Run3ScoutingMuon> >  (parameters.getParameter<edm::InputTag>("muonPacker"));
   svsToken = consumes<edm::View<Run3ScoutingVertex> >  (parameters.getParameter<edm::InputTag>("svPacker"));
   primaryVerticesToken = consumes<edm::View<Run3ScoutingVertex> >  (parameters.getParameter<edm::InputTag>("primaryVertexPacker"));
}


// Destructor
efficiencyMC::~efficiencyMC() {
}


// beginJob (Before first event)
void efficiencyMC::beginJob() {
    std::cout << "Begin Job" << std::endl;

    // Init the file and the TTree
    output_filename = parameters.getParameter<std::string>("nameOfOutput");
    file_out = new TFile(output_filename.c_str(), "RECREATE");

    // Initialize histograms for counters
    event_counts = new TH1F("event_counts", "Event Counters", 2, 0, 2); // Bin 0: Total events, Bin 1: Trigger-passed events
}

// endJob (After event loop has finished)
void efficiencyMC::endJob() {

    std::cout << "End Job" << std::endl;
    file_out->cd();

    // Write the counters histogram
    event_counts->Write();

    // Write other histograms and efficiencies
    histo_pt_gen->Write();
    histo_eta_gen->Write();
    histo_phi_gen->Write();
    histo_lxy_gen->Write();
    histo_pt_gen_matched->Write();
    histo_eta_gen_matched->Write();
    histo_phi_gen_matched->Write();
    histo_lxy_gen_matched->Write();
    histo_pt_vs_lxy_gen->Write();

    histo_pt_comp_matched->Write();

    histo_pt_scouting->Write();
    histo_eta_scouting->Write();
    histo_phi_scouting->Write();
    histo_lxy_scouting->Write();

    histo_pt_scouting_matched->Write();
    histo_eta_scouting_matched->Write();
    histo_phi_scouting_matched->Write();
    histo_lxy_scouting_matched->Write();

    histo_pt_resolution->Write();
    histo_pt_resolution_vs_abseta->Write();
    histo_pt_resolution_vs_lxy->Write();
    
    reco_efficiency_pt_scouting->Write();
    reco_efficiency_lxy_scouting->Write();

    efficiency_trg_DoubleMuon->Write();
    efficiency_minpt_DoubleMuon->Write();
    efficiency_minlxy_DoubleMuon->Write();
    efficiency_maxlxy_DoubleMuon->Write();
    efficiency_npv_DoubleMuon->Write();

    histo_npv->Write();

    file_out->Close();
}


// fillDescriptions
void efficiencyMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

// Analyze (per event)
void efficiencyMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    event_counts->Fill(0); // Fill bin 0 for total events

    bool validMuons = iEvent.getByToken(muonsToken, muons);
    bool validSV = iEvent.getByToken(svsToken, svs);
    bool validPV = iEvent.getByToken(primaryVerticesToken, primaryVertices);
    iEvent.getByToken(genToken, gens);

    // Count primary vertices
    int nPrimaryVertices = 0;
    if (validPV) {
        nPrimaryVertices = primaryVertices->size();
        histo_npv->Fill(nPrimaryVertices);
    }

    // // Debug output
    // if (validMuons) {
    //     std::cout << "Event has " << muons->size() << " scouting muons" << std::endl;
    // } else {
    //     std::cout << "No valid scouting muons in this event" << std::endl;
    // }

    // Fill histograms

    if (validMuons) {
        for (const auto& mu : *muons){
            histo_pt_scouting->Fill(mu.pt());
            histo_eta_scouting->Fill(mu.eta());
            histo_phi_scouting->Fill(mu.phi());
        }
    }

    if (validSV) {
        for (const auto& sv : *svs){
            histo_lxy_scouting->Fill(sqrt(sv.x()*sv.x() + sv.y()*sv.y()));
        }
    }

    // Fill efficiencies (wrt generation)
    std::vector<float> muon_pt;
    std::vector<float> muon_lxy;
    std::vector<unsigned int> matched_scouting;
    
    for (const auto& gp : *gens) {
        // remove particles outside of the detector range
        if (fabs(gp.eta()) > 2.4)
            continue;
        // remove non-muons
        if (fabs(gp.pdgId()) != 13)
            continue;
        // remove non-final state particles
        if (gp.status() != 1)
            continue;
            
        // Get mother particle
        reco::GenParticleRef mref;
        reco::GenParticle m;
        if (gp.mother()->pdgId() == gp.pdgId()) {
            mref = gp.motherRef();
            m = *mref;
            while (m.pdgId() == m.mother()->pdgId()) {
                mref = m.motherRef();
                m = *mref;
            }
        } else {
            m = gp;
        }

        //// Supported mother particles
        // 1023: dark photon from HAHM
        // 6000211: LL scalar from b-hadron decay
        // 333: phi (BPH)
        if (fabs(m.mother()->pdgId()) != 1023 && fabs(m.mother()->pdgId()) != 6000211 && fabs(m.mother()->pdgId()) != 333)
            continue;

        // Get Pt
        muon_pt.push_back(gp.pt());

        // Get lxy distance
        double vx = m.vx();
        double vy = m.vy();
        double vz = m.vz();
        double lxy = sqrt(vx*vx + vy*vy);

        // Get lxy
        muon_lxy.push_back(lxy);

        histo_pt_gen->Fill(gp.pt());
        histo_eta_gen->Fill(gp.eta());
        histo_phi_gen->Fill(gp.phi());
        histo_lxy_gen->Fill(lxy);
        histo_pt_vs_lxy_gen->Fill(gp.pt(), lxy);


        //
        // Match with scouting muons
        float delR_scouting = 999.;
        Run3ScoutingMuon best_scouting;
        unsigned int best_idx_scouting = 999;
        if (validMuons) {
            for (unsigned int i = 0; i < muons->size(); i++){
                if (std::find(matched_scouting.begin(), matched_scouting.end(), i) != matched_scouting.end())
                    continue;
                const auto& mu = muons->at(i);
                float _delPhi = deltaPhi(gp.phi(), mu.phi());
                float _delR = sqrt(_delPhi*_delPhi + (mu.eta() - gp.eta())*(mu.eta() - gp.eta()));
                float _delPt = fabs(1./gp.pt() - 1./mu.pt());
                if (_delR < delR_scouting) {
                    best_idx_scouting = i;
                    best_scouting = mu;
                    delR_scouting = _delR;
                }
            }
            // Fill efficiency
            if (delR_scouting < 0.1) {
                //std::cout << "Matched muon found! delR = " << delR_scouting << ", gen_pt = " << gp.pt() << ", scout_pt = " << best_scouting.pt() << std::endl;
                matched_scouting.push_back(best_idx_scouting);
                // Fill matched gen muon histograms
                histo_pt_gen_matched->Fill(gp.pt());
                histo_eta_gen_matched->Fill(gp.eta());
                histo_phi_gen_matched->Fill(gp.phi());
                histo_lxy_gen_matched->Fill(lxy);
                histo_pt_comp_matched->Fill(gp.pt(), best_scouting.pt());
                // Fill matched scouting muon histograms
                histo_pt_scouting_matched->Fill(best_scouting.pt());
                histo_eta_scouting_matched->Fill(best_scouting.eta());
                histo_phi_scouting_matched->Fill(best_scouting.phi());

		// Fill pt resolution histogram
		float pt_rel_diff = (best_scouting.pt() - gp.pt()) / gp.pt();
		histo_pt_resolution->Fill(gp.pt(), pt_rel_diff);
		histo_pt_resolution_vs_abseta->Fill(fabs(gp.eta()), pt_rel_diff);
		histo_pt_resolution_vs_lxy->Fill(lxy, pt_rel_diff);
		
		
                // Find the vertex associated with the matched scouting muon
                double lxy_scouting = -1;
                if (validSV && !best_scouting.vtxIndx().empty()) {
                    int vtxIdx = best_scouting.vtxIndx().at(0);
                    if (vtxIdx >= 0 && vtxIdx < (int)svs->size()) {
                        const auto& vtx = svs->at(vtxIdx);
                        double lxy_scouting = sqrt(vtx.x()*vtx.x() + vtx.y()*vtx.y());
                        histo_lxy_scouting_matched->Fill(lxy_scouting);
                    }
                }

                // Fill reconstruction efficiency histograms
                reco_efficiency_pt_scouting->Fill(true, gp.pt());
                reco_efficiency_lxy_scouting->Fill(true, lxy);
            } else {
                reco_efficiency_pt_scouting->Fill(false, gp.pt());
                reco_efficiency_lxy_scouting->Fill(false, lxy);
            }
        }

    }
    // Trigger evaluation
    bool passDiMuHLT = false;

    // Skip trigger evaluation if running on data
    if (!isData) {
        if (muon_pt.size() > 1) {
            std::sort( std::begin(muon_pt), std::end(muon_pt), [&](int i1, int i2){ return i1 > i2; });
            double minlxy = *std::min_element(muon_lxy.begin(), muon_lxy.end());
            double maxlxy = *std::max_element(muon_lxy.begin(), muon_lxy.end());
            //std::cout << minlxy << " " << maxlxy << std::endl;

            if (triggerCache_.setEvent(iEvent, iSetup)){
                const auto& vts_dimu(triggerExpression::parse("DST_PFScouting_DoubleMuonNoVtx_v5"));
                if (vts_dimu){
                    vts_dimu->init(triggerCache_);
                    passDiMuHLT = (*vts_dimu)(triggerCache_);
                }
            } 
            if (passDiMuHLT) {
                event_counts->Fill(1); // Fill bin 1 for trigger-passed events
            }
            
            // Debug: Print info when trigger fails but we have good reco muons
            if (!passDiMuHLT && validMuons && muons->size() >= 2) {
                std::cout << "=== TRIGGER FAILED DEBUG ===" << std::endl;
                std::cout << "Gen muons: " << muon_pt.size() << " with pt: ";
                for (float pt : muon_pt) std::cout << pt << " ";
                std::cout << std::endl;
                
                std::cout << "Reco muons: " << muons->size() << " with pt: ";
                for (const auto& mu : *muons) std::cout << mu.pt() << " ";
                std::cout << std::endl;
                
                // Check if muons meet basic pt requirements
                int nMuonsAbove3 = 0;
                int nMuonsAbove5 = 0;
                for (const auto& mu : *muons) {
                    if (mu.pt() > 3.0) nMuonsAbove3++;
                    if (mu.pt() > 5.0) nMuonsAbove5++;
                }
                std::cout << "Muons above 3 GeV: " << nMuonsAbove3 << ", above 5 GeV: " << nMuonsAbove5 << std::endl;
                std::cout << "Min lxy: " << minlxy << ", Max lxy: " << maxlxy << std::endl;
                std::cout << "============================" << std::endl;
            }
            
            efficiency_trg_DoubleMuon->Fill(passDiMuHLT, muon_pt.at(0), muon_pt.at(1));
            efficiency_minpt_DoubleMuon->Fill(passDiMuHLT, muon_pt.at(1));
            efficiency_minlxy_DoubleMuon->Fill(passDiMuHLT, minlxy);
            efficiency_maxlxy_DoubleMuon->Fill(passDiMuHLT, maxlxy);
            efficiency_npv_DoubleMuon->Fill(passDiMuHLT, nPrimaryVertices);
        }
    }

}


DEFINE_FWK_MODULE(efficiencyMC);
