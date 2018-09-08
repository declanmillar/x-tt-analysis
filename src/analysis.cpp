#include "analysis.hpp"
#include "trim.hpp"
#include "progress-bar.hpp"
#include "bool-to-string.hpp"
#include "lester_mt2_bisect.h"
#include "neutrino-weighter.hpp"
#include "atlas-style.hpp"
#include "solve-poly.hpp"
#include "kinematic-reconstructer.hpp"
#include "highest-pt.hpp"
#include "match-bjets-to-leps.hpp"
#include "get-parameter.hpp"
#include <exception>
#include <regex>

void Analysis::Run()
{
    this->Loop();
    this->PostLoop();
}

void Analysis::EachEvent(double weight)
{
    if (m_debug) cout << "Starting EachEvent ...\n";
    if (m_debug) cout << "Event weight = " << weight << "\n";

    // hard process particles
    this->GetHardParticles();

    pair<TLorentzVector, TLorentzVector> p_l_truth = make_pair(m_hardLepP->P4(), m_hardLepM->P4());
    pair<TLorentzVector, TLorentzVector> p_b_truth = make_pair(m_hardB->P4(), m_hardBbar->P4());
    pair<TLorentzVector, TLorentzVector> p_t_truth = make_pair(m_hardTop->P4(), m_hardTbar->P4());
    pair<TLorentzVector, TLorentzVector> p_v_truth = make_pair(m_hardNu->P4(), m_hardNuBar->P4());
    TLorentzVector p_ttbar_truth = p_t_truth.first + p_t_truth.second;
    TLorentzVector p_miss_truth = p_v_truth.first + p_v_truth.second;
    double mass_ttbar_truth = p_ttbar_truth.M() / 1000;
    double pT_top_truth = p_t_truth.first.Pt();
    double pT_tbar_truth = p_t_truth.second.Pt();
    double ET_miss_truth = p_miss_truth.Pt();
    double dR_l1b1 = p_l_truth.first.DeltaR(p_b_truth.first);
    double dR_l2b2 = p_l_truth.second.DeltaR(p_b_truth.second);
    double dR_t1b1 = p_t_truth.first.DeltaR(p_b_truth.first);
    double dR_t2b2 = p_t_truth.second.DeltaR(p_b_truth.second);
    double dR_t1l1 = p_t_truth.first.DeltaR(p_l_truth.first);
    double dR_t2l2 = p_t_truth.second.DeltaR(p_l_truth.second);
    double dR_t1t2 = p_t_truth.first.DeltaR(p_t_truth.second);
    vector<double> eff_values = {mass_ttbar_truth, pT_top_truth, pT_tbar_truth};
    m_mass_ttbar_truth->push_back(mass_ttbar_truth);
    m_dR_lb_truth->push_back(dR_l1b1);
    
    double HT_truth = p_l_truth.first.Pt() + p_l_truth.second.Pt() + p_b_truth.first.Pt() + p_b_truth.second.Pt();
    double HTmet_truth = HT_truth + ET_miss_truth;
    TLorentzVector p_bbll_truth = p_l_truth.first + p_l_truth.second + p_b_truth.first + p_b_truth.second;
    double mass_bbll_truth = p_bbll_truth.M();
    double pT_bbll_truth = p_bbll_truth.Pt();
    double KT_truth = sqrt(mass_bbll_truth * mass_bbll_truth + pT_bbll_truth * pT_bbll_truth) + ET_miss_truth;
    HT_truth = HT_truth / 1000;
    HTmet_truth = HTmet_truth / 1000;
    KT_truth = KT_truth / 1000;
    
    // double deltaEta_tt_truth = p_t_truth.first.DeltaEta(p_t_truth.first);
    // h_deltaEta_tt_truth->Fill(deltaEta_tt_truth, weight);
    
    double deltaPhi_tt_truth = p_t_truth.first.DeltaPhi(p_t_truth.first);
    h_deltaPhi_tt_truth->Fill(deltaPhi_tt_truth, weight);
    
    // double dRmax = min(dR_l1b1, dR_l2b2) / 2;
    // cout << mass_ttbar_truth << " " << dR_l1b1 << "\n";
    
    h_dR_lb_truth->Fill(dR_l1b1, weight);
    h_dR_lb_truth->Fill(dR_l2b2, weight);
    h_dR_t1l1_truth->Fill(dR_t1l1, weight);
    h_dR_t2l2_truth->Fill(dR_t2l2, weight);
    h_dR_t1b1_truth->Fill(dR_t1b1, weight);
    h_dR_t2b2_truth->Fill(dR_t2b2, weight);
    h_dR_t1t2_truth->Fill(dR_t1t2, weight);
    
    h2_dR_lb_mtt_truth->Fill(mass_ttbar_truth, dR_l1b1, weight);
    h2_dR_lb_mtt_truth->Fill(mass_ttbar_truth, dR_l2b2, weight);
    h2_dR_t1l1_mtt_truth->Fill(mass_ttbar_truth, dR_t1l1, weight);
    h2_dR_t2l2_mtt_truth->Fill(mass_ttbar_truth, dR_t2l2, weight);
    h2_dR_t1b1_mtt_truth->Fill(mass_ttbar_truth, dR_t1b1, weight);
    h2_dR_t2b2_mtt_truth->Fill(mass_ttbar_truth, dR_t2b2, weight);
    h2_dR_l1b1_pTl_truth->Fill(p_l_truth.first.Pt(), dR_t1b1, weight);
    h2_dR_l2b2_pTl_truth->Fill(p_l_truth.second.Pt(), dR_t2b2, weight);

    h_pT_top_truth->Fill(pT_top_truth, weight);
    h_eta_top_truth->Fill(p_t_truth.first.Eta(), weight);
    h_y_top_truth->Fill(p_t_truth.first.Rapidity(), weight);
    h_phi_top_truth->Fill(p_t_truth.first.Phi() / m_pi, weight);
    h_mass_top_truth->Fill(p_t_truth.first.M(), weight);

    h_pT_tbar_truth->Fill(pT_tbar_truth, weight);
    h_eta_tbar_truth->Fill(p_t_truth.second.Eta(), weight);
    h_y_tbar_truth->Fill(p_t_truth.second.Rapidity(), weight);
    h_phi_tbar_truth->Fill(p_t_truth.second.Phi() / m_pi, weight);
    h_mass_tbar_truth->Fill(p_t_truth.second.M(), weight);

    h_pT_ttbar_truth->Fill(p_ttbar_truth.Pt(), weight);
    h_eta_ttbar_truth->Fill(p_ttbar_truth.Eta(), weight);
    h_phi_ttbar_truth->Fill(p_ttbar_truth.Phi() / m_pi, weight);
    h_mass_ttbar_truth->Fill(mass_ttbar_truth, weight);
    h_y_ttbar_truth->Fill(p_ttbar_truth.Rapidity(), weight);
    
    // get truth 4-momenta in reconstructed truth ttbar frame
    auto pttbar_top_truth = p_t_truth.first, pttbar_tbar_truth = p_t_truth.second;
    // , pttbar_ttbar = p_ttbar;
    // auto pttbar_b1 = p_b1, pttbar_b2 = p_b2, pttbar_v1 = p_v1, pttbar_v2 = p_v2;
    // auto pttbar_l = p_l;
    TVector3 v_ttbar_truth = - (p_ttbar_truth).BoostVector();
    pttbar_top_truth.Boost(v_ttbar_truth);
    pttbar_tbar_truth.Boost(v_ttbar_truth);
    // pttbar_ttbar.Boost(v_ttbar);
    // pttbar_b1.Boost(v_ttbar);
    // pttbar_b2.Boost(v_ttbar);
    // pttbar_l.first.Boost(v_ttbar);
    // pttbar_l.second.Boost(v_ttbar);
    // pttbar_v1.Boost(v_ttbar);
    // pttbar_v2.Boost(v_ttbar);

    // get truth momenta in reconstructed truth top frame
    // auto ptop_top = p_top, ptop_tbar = p_tbar, ptop_ttbar = p_ttbar;
    // auto ptop_b1 = p_b1, ptop_b2 = p_b2, ptop_v1 = p_v1, ptop_v2 = p_v2;
    auto ptop_l_truth = p_l_truth;
    TVector3 v_top_truth = - (p_t_truth.first).BoostVector();
    // ptop_top.Boost(v_top);
    // ptop_tbar.Boost(v_top);
    // ptop_ttbar.Boost(v_top);
    // ptop_b1.Boost(v_top);
    // ptop_b2.Boost(v_top);
    ptop_l_truth.first.Boost(v_top_truth);
    // ptop_l.second.Boost(v_top);
    // ptop_v1.Boost(v_top);
    // ptop_v2.Boost(v_top);
    
    // get truth momenta in reconstructed truth tbar frame
    // auto ptbar_top = p_top, ptbar_tbar = p_tbar, ptbar_ttbar = p_ttbar;
    // auto ptbar_b1 = p_b1, ptbar_b2 = p_b2, ptbar_v1 = p_v1, ptbar_v2 = p_v2;
    auto ptbar_l_truth = p_l_truth;
    TVector3 v_tbar_truth = - (p_t_truth.second).BoostVector();
    // ptbar_top.Boost(v_tbar);
    // ptbar_tbar.Boost(v_tbar);
    // ptbar_ttbar.Boost(v_tbar);
    // ptbar_b1.Boost(v_tbar);
    // ptbar_b2.Boost(v_tbar);
    // ptbar_l.first.Boost(v_tbar);
    ptbar_l_truth.second.Boost(v_tbar_truth);
    // ptbar_v1.Boost(v_tbar);
    // ptbar_v2.Boost(v_tbar);
    
    double cosTheta_tl1_truth = cos(ptop_l_truth.first.Angle(pttbar_top_truth.Vect()));
    h_cosTheta1_truth->Fill(cosTheta_tl1_truth, weight);
    double cosTheta_tl2_truth = cos(ptbar_l_truth.second.Angle(pttbar_top_truth.Vect()));
    h_cosTheta2_truth->Fill(cosTheta_tl2_truth, weight);
    
    auto deltaY_top_truth = abs(p_t_truth.first.Rapidity()) - abs(p_t_truth.second.Rapidity());
    h_deltaY_top_truth->Fill(deltaY_top_truth, weight);
    
    h_HT_truth->Fill(HT_truth, weight);
    h_HTmet_truth->Fill(HTmet_truth, weight);
    h_KT_truth->Fill(KT_truth, weight);
    
    h_ETmiss_truth->Fill(ET_miss_truth, weight);
    h2_mtt_truth_ETmiss_truth->Fill(mass_ttbar_truth, p_miss_truth.Pt(), weight);
    
    this->GetTruthParticles();

    h_n_truthElectrons->Fill(m_truthElectrons->size(), weight);
    h_n_truthMuons->Fill(m_truthMuons->size(), weight);
    h_n_truthBquarks->Fill(m_truthBquarks->size(), weight);
    
    this->SelectElectrons();
    this->SelectMuons();
    this->SelectJets();
    
    this->FillPurities(1);
    // this->FillTaggedHistograms();
    int n_tagged = 0;
    for (auto tag : *m_muon_truth_tags) if (tag) ++n_tagged;
    h_n_muons_truth_tagged->Fill(n_tagged, 1);
    n_tagged = 0;
    for (auto tag : *m_electron_truth_tags) if (tag) ++n_tagged;
    h_n_electrons_truth_tagged->Fill(n_tagged, 1);
    n_tagged = 0;
    
    int ijet = 0;
    for (auto tag : *m_jet_truth_tags) {
        Jet *jet = m_jets->at(ijet);
        if (tag) {
            int nTracks = 0;
            for (int k = 0; k < jet->Constituents.GetEntriesFast(); ++k) {
                TObject* object = jet->Constituents.At(k);
                if (object == 0) continue;
                if (object->IsA() == Track::Class()) ++nTracks;
            }
            h_ntracks_truth_tagged_jets->Fill(nTracks);
            ++n_tagged;
        }
        ++ijet;
    }
        
    h_n_jets_truth_tagged->Fill(n_tagged, 1);

    this->IsolateElectrons();
    this->IsolateMuons();

    // objects pass selection
    h_n_selElectrons->Fill(m_electrons->size(), weight);
    h_n_selMuons->Fill(m_muons->size(), weight);
    h_n_selJets->Fill(m_jets->size(), weight);
    
    this->FillPurities(2);
    n_tagged = 0;
    for (auto tag : *m_jet_truth_tags) if (tag) ++n_tagged;
    h_n_selJets_truth_tagged->Fill(n_tagged, 1);
    n_tagged = 0;
    for (auto tag : *m_electron_truth_tags) if (tag) ++n_tagged;
    h_n_selElectrons_truth_tagged->Fill(n_tagged, 1);
    n_tagged = 0;
    for (auto tag : *m_muon_truth_tags) if (tag) ++n_tagged;
    h_n_selMuons_truth_tagged->Fill(n_tagged, 1);
    
    this->OverlapRemoval();
    
    h_n_uniqueElectrons->Fill(m_electrons->size(), weight);
    h_n_uniqueMuons->Fill(m_muons->size(), weight);
    h_n_uniqueJets->Fill(m_jets->size(), weight);
    
    this->FillPurities(3);
    n_tagged = 0;
    for (auto tag : *m_jet_truth_tags) if (tag) ++n_tagged;
    h_n_uniqueJets_truth_tagged->Fill(n_tagged, 1);
    n_tagged = 0;
    for (auto tag : *m_electron_truth_tags) if (tag) ++n_tagged;
    h_n_uniqueElectrons_truth_tagged->Fill(n_tagged, 1);
    n_tagged = 0;
    for (auto tag : *m_muon_truth_tags) if (tag) ++n_tagged;
    h_n_uniqueMuons_truth_tagged->Fill(n_tagged, 1);
    
    // cuts
    UpdateCutflow(c_events, true);
    
    // leptons
    if (!this->ExactlyTwoLeptons()) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_2l_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_2l_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);

    this->AssignChannel();

    if (!this->OppositeCharge()) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_oc_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_oc_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);

    pair<TLorentzVector, TLorentzVector> p_l = this->GetLeptonMomenta();
    if (!this->SufficientMll(p_l)) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_mll_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_mll_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);

    if (!this->OutsideZmassWindow(p_l)) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_mZ_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_mZ_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);

    MissingET* missingET = (MissingET*) b_MissingET->At(0);
    double ET_miss = missingET->MET;
    h_ETmiss->Fill(ET_miss, weight);

    if (!this->SufficientMET(ET_miss)) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_ETmiss_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_ETmiss_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);

    TLorentzVector p_miss;
    p_miss.SetPtEtaPhiM(ET_miss, missingET->Eta, missingET->Phi, 0.0);
    
    if (!this->SufficientJets()) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_2j_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_2j_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);

    if (!this->SufficientBtags()) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_2b_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_2b_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);

    // Get jet momenta and split by b-tag
    vector<TLorentzVector> p_j, p_b, p_q;
    for (int i = 0; i < m_jets->size(); ++i) {
        Jet *jet = (Jet*) m_jets->at(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        // if (jet->BTag & (1 << i))
        if (jet->BTag > 0) {
            if (m_debug) cout << "b-jet pT = " << p.Pt() << "\n";
            p_b.push_back(p);
            h_pT_bjets->Fill(p.Pt(), weight);
            h_eta_bjets->Fill(p.Eta(), weight);
        }
        else {
            p_q.push_back(p);
            h_pT_qjets->Fill(p.Pt(), weight);
            h_eta_qjets->Fill(p.Eta(), weight);
        }
        p_j.push_back(p);
        h_pT_jets->Fill(p.Pt(), weight);
        h_eta_jets->Fill(p.Eta(), weight);
    }

    if (m_debug) {
        cout << "p_l+   = "; p_l.first.Print();
        cout << "p_l-   = "; p_l.second.Print();
        for (int i = 0; i < p_b.size(); ++i) {
            cout << "p_b" << i << "   = ";
            p_b.at(i).Print();
        }
        for (int i = 0; i < p_q.size(); ++i) {
            cout << "p_q" << i << "   = ";
            p_q.at(i).Print();
        }
        cout << "p_miss = "; p_miss.Print();
    }
    
    // select leading jets (hypothetical b-jets)
    pair<TLorentzVector, TLorentzVector> p_b_hypo;
    if (m_bTags >= 2) {
        if (m_debug) cout << "Selecting two b-tagged jets with highest pT ...\n";
        p_b_hypo = TwoHighestPt(p_b);
    }
    else if (m_bTags == 0) {
        if (m_debug) cout << "Selecting two jets with highest pT ...\n";
        p_b_hypo = TwoHighestPt(p_q);
    }
    else if (m_bTags == 1) {
        if (m_debug) cout << "Selecting b-tagged jet and jet with highest pT\n";
        p_b_hypo.first = p_b.at(0);
        p_b_hypo.second = HighestPt(p_q);
    }
    else {
        cout << "ERROR: invalid b-tags\n";
    }
    if (m_debug) cout << "Collected hypothetical b-quark jets.\n";
    
    // HT in arXiv:1604.05538v2 is defined as the scalar sum of the pT of the two leading jets and leptons
    double HT = p_l.first.Pt() + p_l.second.Pt() + p_b_hypo.first.Pt() + p_b_hypo.second.Pt();
    
    if (!this->SufficientHT(HT)) {
        this->FillCutsEfficiencies(eff_values, 0);
        h_eff_cut_HT_mass_ttbar_truth->Fill(mass_ttbar_truth, 0);
        return;
    }
    h_eff_cut_HT_mass_ttbar_truth->Fill(mass_ttbar_truth, 1);
    
    this->FillCutsEfficiencies(eff_values, 1);
    h_n_passElectrons->Fill(m_electrons->size(), weight);
    h_n_passMuons->Fill(m_muons->size(), weight);
    h_n_passJets->Fill(m_jets->size(), weight);
    h_n_passBjets->Fill(m_bTags, weight);
    
    if (m_debug) cout << "Passed all cuts\n";
    
    double HTmet = HT + ET_miss;
    
    double HTjMET = p_l.first.Pt() + p_l.second.Pt() + ET_miss;
    for (auto& p : p_j) HTjMET += p.Pt();
    
    TLorentzVector p_bbll = p_l.first + p_l.second + p_b_hypo.first + p_b_hypo.second;
    double mass_bbll = p_bbll.M();
    double pT_bbll = p_bbll.Pt();
    double KT = sqrt(mass_bbll * mass_bbll + pT_bbll * pT_bbll) + ET_miss;

    TLorentzVector p_vis = p_l.first + p_l.second;
    for (auto& p : p_j) p_vis += p;
    double mass_vis = p_vis.M();
    double pT_vis = p_vis.Pt();
    double KTj = sqrt(mass_vis * mass_vis + pT_vis * pT_vis) + ET_miss;
    
    mass_bbll = mass_bbll / 1000;
    HT = HT / 1000;
    HTmet = HTmet / 1000;
    KT = KT / 1000;
    mass_vis = mass_vis / 1000;
    HTjMET = HT / 1000;
    KTj = KT / 1000;
    
    if (m_debug) cout << "filling transverse histograms\n";

    h_mass_bbll->Fill(mass_bbll, weight);
    h_HT->Fill(HT, weight);
    h_HTmet->Fill(HTmet, weight);
    h_KT->Fill(KT, weight);
    h_mass_vis->Fill(mass_vis, weight);
    h_HTjMET->Fill(HTjMET, weight);
    h_KTj->Fill(KTj, weight);

    if (m_debug) cout << "filled transverse histograms\n";

    auto deltaAbsEta_l = abs(p_l.first.Eta()) - abs(p_l.second.Eta());

    TLorentzVector p_ll = p_l.first + p_l.second;
    auto y_ll = p_ll.Rapidity();
    auto pll_l = p_l.first;
    TVector3 v_ll = - (p_ll).BoostVector();
    pll_l.Boost(v_ll);
    double cosTheta_l = pll_l.CosTheta();
    double cosThetaStar_l = int(y_ll / abs(y_ll)) * cosTheta_l;

    double deltaPhi_ll = p_l.first.DeltaPhi(p_l.second) / m_pi;

    h_cosTheta_l->Fill(cosTheta_l, weight);
    h_cosThetaStar_l->Fill(cosThetaStar_l, weight);

    h_pT_l1->Fill(p_l.first.Pt(), weight);
    h_pT_l2->Fill(p_l.second.Pt(), weight);
    h_eta_l1->Fill(p_l.first.Eta(), weight);
    h_eta_l2->Fill(p_l.second.Eta(), weight);

    h_deltaPhi_ll->Fill(deltaPhi_ll, weight);
    if (deltaPhi_ll > 0.5) {
        h_HT_DphiF->Fill(HT, weight);
        h_KT_DphiF->Fill(KT, weight);
    }
    if (deltaPhi_ll < 0.5) {
        h_HT_DphiB->Fill(HT, weight);
        h_KT_DphiB->Fill(KT, weight);
    }

    h_cosThetaStar_l->Fill(cosThetaStar_l, weight);
    if (cosThetaStar_l > 0) {
        h_HT_lF->Fill(HT, weight);
        h_KT_lF->Fill(KT, weight);
    }
    if (cosThetaStar_l < 0) {
        h_HT_lB->Fill(HT, weight);
        h_KT_lB->Fill(KT, weight);
    }

    h_deltaEta_l->Fill(deltaAbsEta_l, weight);
    if (deltaAbsEta_l > 0) {
        h_HT_lCF->Fill(HT, weight);
        h_KT_lCF->Fill(KT, weight);
    }
    if (deltaAbsEta_l < 0) {
        h_HT_lCB->Fill(HT, weight);
        h_KT_lCB->Fill(KT, weight);
    }
    
    h2_HT_deltaPhi->Fill(HT, deltaPhi_ll, weight);
    h2_HTmet_deltaPhi->Fill(HTmet, deltaPhi_ll, weight);
    h2_mvis_deltaPhi->Fill(mass_vis, deltaPhi_ll, weight);
    h2_KT_deltaPhi->Fill(KT, deltaPhi_ll, weight);

    // reconstruction
    if (m_reconstruction == "KIN" or m_reconstruction == "NuW") {
        if (m_debug) cout << "Starting reconstruction ...\n";

        TLorentzVector p_top, p_tbar, p_ttbar, p_b1, p_b2, p_v1, p_v2;

        if (m_debug) cout << "Matching b-jets to leptons ...\n";
        auto p_b_match = MatchBjetsToLeps(p_l, p_b_hypo);

        if (m_reconstruction == "KIN") {
            if (m_debug) cout << "Setting up kinematic reconstruction ...\n";
            KinematicReconstructer KIN = KinematicReconstructer(m_mass_b, m_mass_W);
            bool isSolution = KIN.Reconstruct(p_l, p_b_match, p_miss, m_mass_top);
            if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_match, p_miss, m_mass_top - 0.5);
            if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_match, p_miss, m_mass_top + 0.5);
            if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_match, p_miss, m_mass_top - 1.0);
            if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_match, p_miss, m_mass_top + 1.0);
            if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_match, p_miss, m_mass_top - 1.5);
            if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_match, p_miss, m_mass_top + 1.5);

            // std::pair<TLorentzVector, TLorentzVector> p_b_rmatch = std::make_pair(p_b_match.second, p_b_match.first);
            // if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_rmatch, p_miss, m_mass_top);
            // if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_rmatch, p_miss, m_mass_top - 0.5);
            // if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_rmatch, p_miss, m_mass_top + 0.5);
            // if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_rmatch, p_miss, m_mass_top - 1.0);
            // if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_rmatch, p_miss, m_mass_top + 1.0);
            // if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_rmatch, p_miss, m_mass_top - 1.5);
            // if (!isSolution) isSolution = KIN.Reconstruct(p_l, p_b_rmatch, p_miss, m_mass_top + 1.5);

            if (isSolution) {
                this->UpdateCutflow(c_validSolution, true);
                this->FillRecoEfficiencies(eff_values, 1);
                p_top = KIN.GetTop();
                p_tbar = KIN.GetTbar();
                p_ttbar = KIN.GetTtbar();
                p_b1 = KIN.GetB();
                p_b2 = KIN.GetBbar();
                p_v1 = KIN.GetNu();
                p_v2 = KIN.GetNubar();
            }
            else {
                this->UpdateCutflow(c_validSolution, false);
                this->FillRecoEfficiencies(eff_values, 0);
                return;
            }
            if (m_debug) cout << "Finished kinematic reconstruction\n";
        }
        else if (m_reconstruction == "NuW") {
            if (m_debug) cout << "Setting up NuW reconstruction ...\n";
            NeutrinoWeighter nuW = NeutrinoWeighter(1, p_l.first.Pt() + p_l.first.Phi(), m_mass_b); // 2nd argument is random seed same for specific event
            double weight_max = nuW.Reconstruct(p_l.first, p_l.second, p_b_match.first, p_b_match.second, p_miss.Px(), p_miss.Py(), p_miss.Phi());
            // if (weight_max <= 0.0) weight_max = nuW.Reconstruct(p_l.first, p_l.second, p_b_match.second, p_b_match.first, p_miss.Px(), p_miss.Py(), p_miss.Phi());
            if (weight_max > 0.0) {
                this->UpdateCutflow(c_validSolution, true);
                this->FillRecoEfficiencies(eff_values, 1);
                p_top = nuW.GetTop();
                p_tbar = nuW.GetTbar();
                p_ttbar = nuW.GetTtbar();
                p_b1 = nuW.GetB();
                p_b2 = nuW.GetBbar();
                p_v1 = nuW.GetNu();
                p_v2 = nuW.GetNubar();
            }
            else {
                this->UpdateCutflow(c_validSolution, false);
                this->FillRecoEfficiencies(eff_values, 0);
                return;
            }
            if (m_debug) cout << "Finished NuW reconstruction\n";
        }
        else {
            cout << "error: Reconstruction = {KIN, NuW}\n";
        }
        
        h_dR_l->Fill(min(p_l.first.DeltaR(p_l_truth.first), p_l.first.DeltaR(p_l_truth.first)));
        h_dR_l->Fill(min(p_l.second.DeltaR(p_l_truth.second), p_l.second.DeltaR(p_l_truth.second)));

        TLorentzVector p_W1 = p_l.first + p_v1;
        TLorentzVector p_W2 = p_l.second + p_v2;

        auto mass_ttbar = p_ttbar.M();
        auto y_ttbar = p_ttbar.Rapidity();
        auto cosTheta_ttbar = cos(p_top.Angle(p_tbar.Vect()));
        auto DeltaY_top = abs(p_top.Rapidity()) - abs(p_tbar.Rapidity());
        double deltaPhi_tt = p_top.DeltaPhi(p_tbar);
        h_deltaPhi_tt->Fill(deltaPhi_tt, weight);

        // get 4-momenta in reconstructed ttbar frame
        auto pttbar_top = p_top, pttbar_tbar = p_tbar, pttbar_ttbar = p_ttbar;
        auto pttbar_b1 = p_b1, pttbar_b2 = p_b2, pttbar_v1 = p_v1, pttbar_v2 = p_v2;
        auto pttbar_l = p_l;
        TVector3 v_ttbar = - (p_ttbar).BoostVector();
        pttbar_top.Boost(v_ttbar);
        pttbar_tbar.Boost(v_ttbar);
        pttbar_ttbar.Boost(v_ttbar);
        pttbar_b1.Boost(v_ttbar);
        pttbar_b2.Boost(v_ttbar);
        pttbar_l.first.Boost(v_ttbar);
        pttbar_l.second.Boost(v_ttbar);
        pttbar_v1.Boost(v_ttbar);
        pttbar_v2.Boost(v_ttbar);

        // get momenta in reconstructed top frame
        auto ptop_top = p_top, ptt = p_tbar, ptop_ttbar = p_ttbar;
        auto ptop_b1 = p_b1, ptop_b2 = p_b2, ptop_v1 = p_v1, ptop_v2 = p_v2;
        auto ptop_l = p_l;
        TVector3 v_top = - (p_top).BoostVector();
        ptop_top.Boost(v_top);
        ptt.Boost(v_top);
        ptop_ttbar.Boost(v_top);
        ptop_b1.Boost(v_top);
        ptop_b2.Boost(v_top);
        ptop_l.first.Boost(v_top);
        ptop_l.second.Boost(v_top);
        ptop_v1.Boost(v_top);
        ptop_v2.Boost(v_top);
        
        // get momenta in reconstructed tbar frame
        auto ptbar_top = p_top, ptbar_tbar = p_tbar, ptbar_ttbar = p_ttbar;
        auto ptbar_b1 = p_b1, ptbar_b2 = p_b2, ptbar_v1 = p_v1, ptbar_v2 = p_v2;
        auto ptbar_l = p_l;
        TVector3 v_tbar = - (p_tbar).BoostVector();
        ptbar_top.Boost(v_tbar);
        ptbar_tbar.Boost(v_tbar);
        ptbar_ttbar.Boost(v_tbar);
        ptbar_b1.Boost(v_tbar);
        ptbar_b2.Boost(v_tbar);
        ptbar_l.first.Boost(v_tbar);
        ptbar_l.second.Boost(v_tbar);
        ptbar_v1.Boost(v_tbar);
        ptbar_v2.Boost(v_tbar);

        double cosTheta = pttbar_top.CosTheta();
        double cosThetaStar = int(y_ttbar / abs(y_ttbar)) * cosTheta;

        double cosPhi = cos(ptop_l.first.Phi() - ptbar_l.second.Phi());

        double cosTheta_tl1 = cos(ptop_l.first.Angle(pttbar_top.Vect()));
        double cosTheta_tl2 = cos(ptbar_l.second.Angle(pttbar_tbar.Vect()));
        double cos1cos2 = cosTheta_tl1 * cosTheta_tl2;

        if (m_debug) cout << "Filling reconstruction histograms ...\n";

        mass_ttbar = mass_ttbar / 1000;

        h_E_top->Fill(p_top.E(), weight);
        h_pT_top->Fill(p_top.Pt(), weight);
        h_eta_top->Fill(p_top.Eta(), weight);
        h_y_top->Fill(p_top.Rapidity(), weight);
        h_phi_top->Fill(p_top.Phi() / m_pi, weight);
        h_mass_top->Fill(p_top.M(), weight);

        h_E_tbar->Fill(p_tbar.E(), weight);
        h_pT_tbar->Fill(p_tbar.Pt(), weight);
        h_eta_tbar->Fill(p_tbar.Eta(), weight);
        h_y_tbar->Fill(p_tbar.Rapidity(), weight);
        h_phi_tbar->Fill(p_tbar.Phi() / m_pi, weight);
        h_mass_tbar->Fill(p_tbar.M(), weight);

        h_pT_ttbar->Fill(p_ttbar.Pt(), weight);
        h_eta_ttbar->Fill(p_ttbar.Eta(), weight);
        h_phi_ttbar->Fill(p_ttbar.Phi() / m_pi, weight);
        h_mass_ttbar->Fill(mass_ttbar, weight);
        h_y_ttbar->Fill(y_ttbar, weight);

        double dR_top = p_top.DeltaR(p_t_truth.first);
        h_dR_top->Fill(dR_top, weight);

        double dR_tbar = p_tbar.DeltaR(p_t_truth.second);
        h_dR_tbar->Fill(dR_tbar, weight);

        double dR_ttbar = p_ttbar.DeltaR(p_ttbar_truth);
        h_dR_ttbar->Fill(dR_ttbar, weight);
        
        if (dR_top < 0.3 and dR_tbar < 0.3) h_reco_quality->Fill(0.5, 1);
        else h_reco_quality->Fill(0.5, 0);

        h_perf_pT_top->Fill((p_t_truth.first.Pt() - p_top.Pt()) / p_t_truth.first.Pt(), weight);
        h_perf_eta_top->Fill((p_t_truth.first.Eta() - p_top.Eta()) / p_t_truth.first.Eta(), weight);
        h_perf_phi_top->Fill((p_t_truth.first.Phi() - p_top.Phi()) / p_t_truth.first.Phi(), weight);
        h_perf_mass_top->Fill((p_t_truth.first.M() - p_top.M()) / p_top.M(), weight);

        h_perf_pT_tbar->Fill((p_t_truth.second.Pt() - p_tbar.Pt()) / p_t_truth.second.Pt(), weight);
        h_perf_eta_tbar->Fill((p_t_truth.second.Eta() - p_tbar.Eta()) / p_t_truth.second.Eta(), weight);
        h_perf_phi_tbar->Fill((p_t_truth.second.Phi() - p_tbar.Phi()) / p_t_truth.second.Phi(), weight);
        h_perf_mass_tbar->Fill((p_t_truth.second.M() - p_tbar.M()) / p_tbar.M(), weight);

        h_perf_pT_ttbar->Fill((p_ttbar_truth.Pt() - p_ttbar.Pt()) / p_ttbar_truth.Pt(), weight);
        h_perf_eta_ttbar->Fill((p_ttbar_truth.Eta() - p_ttbar.Eta()) / p_ttbar_truth.Eta(), weight);
        h_perf_phi_ttbar->Fill((p_ttbar_truth.Phi() - p_ttbar.Phi()) / p_ttbar_truth.Phi(), weight);
        double perf_mass_ttbar = (mass_ttbar_truth - mass_ttbar) / mass_ttbar_truth;
        h_perf_mass_ttbar->Fill(perf_mass_ttbar, weight);
        
        h_res_pT_tops->Fill(p_t_truth.first.Pt() - p_top.Pt(), weight);
        h_res_y_tops->Fill(p_t_truth.first.Rapidity() - p_top.Rapidity(), weight);
        h_res_phi_tops->Fill(p_t_truth.first.Phi() - p_top.Phi(), weight);
        h_res_eta_tops->Fill(p_t_truth.first.Eta() - p_top.Eta(), weight);
        h_res_pT_tops->Fill(p_t_truth.second.Pt() - p_tbar.Pt(), weight);
        h_res_y_tops->Fill(p_t_truth.second.Rapidity() - p_tbar.Rapidity(), weight);
        h_res_phi_tops->Fill(p_t_truth.second.Phi() - p_tbar.Phi(), weight);
        h_res_eta_tops->Fill(p_t_truth.second.Eta() - p_tbar.Eta(), weight);
        
        h_res_pT_ttbar->Fill(p_ttbar_truth.Pt() - p_ttbar.Pt(), weight);
        h_res_y_ttbar->Fill(p_ttbar_truth.Rapidity() - p_ttbar.Rapidity(), weight);
        h_res_mass_ttbar->Fill(1000* (mass_ttbar_truth - mass_ttbar), weight);
        h_res_eta_ttbar->Fill(p_ttbar_truth.Eta() - p_ttbar.Eta(), weight);

        h2_perf_mass_ttbar->Fill(mass_ttbar_truth, perf_mass_ttbar, weight);
        h2_perf_mass_ttbar_pTtop->Fill(p_t_truth.first.Pt(), perf_mass_ttbar, weight);
        h2_perf_mass_ttbar_pTtbar->Fill(p_t_truth.second.Pt(), perf_mass_ttbar, weight);

        h2_pT_top_TvR->Fill(p_top.Pt(), p_t_truth.first.Pt(), weight);
        h2_eta_top_TvR->Fill(p_top.Eta(), p_t_truth.first.Eta(), weight);
        h2_phi_top_TvR->Fill(p_top.Phi(), p_t_truth.first.Phi(), weight);
        h2_mass_top_TvR->Fill(p_top.M(), p_t_truth.first.M(), weight);

        h2_pT_tbar_TvR->Fill(p_tbar.Pt(), p_t_truth.second.Pt(), weight);
        h2_eta_tbar_TvR->Fill(p_tbar.Eta(), p_t_truth.second.Eta(), weight);
        h2_phi_tbar_TvR->Fill(p_tbar.Phi(), p_t_truth.second.Phi(), weight);
        h2_mass_tbar_TvR->Fill(p_tbar.M(), p_t_truth.second.M(), weight);

        h2_pT_ttbar_TvR->Fill(p_ttbar.Pt(), p_ttbar_truth.Pt(), weight);
        h2_eta_ttbar_TvR->Fill(p_ttbar.Eta(), p_ttbar_truth.Eta(), weight);
        h2_phi_ttbar_TvR->Fill(p_ttbar.Phi(), p_ttbar_truth.Phi(), weight);
        h2_mass_ttbar_TvR->Fill(mass_ttbar, mass_ttbar_truth, weight);

        h_mass_W1->Fill(p_W1.M(), weight);
        h_mass_W2->Fill(p_W2.M(), weight);

        if (deltaPhi_ll > 0.5) h_mtt_DphiF->Fill(mass_ttbar, weight);
        if (deltaPhi_ll < 0.5) h_mtt_DphiB->Fill(mass_ttbar, weight);

        h_cosPhi->Fill(cosPhi, weight);
        if (cosPhi > 0) h_mtt_cPhiF->Fill(mass_ttbar, weight);
        if (cosPhi < 0) h_mtt_cPhiB->Fill(mass_ttbar, weight);

        h_cosTheta_ttbar->Fill(cosTheta_ttbar, weight);
        h_cosTheta->Fill(cosTheta, weight);
        h_cosThetaStar->Fill(cosThetaStar, weight);
        if (cosThetaStar > 0) h_mtt_tF->Fill(mass_ttbar, weight);
        if (cosThetaStar < 0) h_mtt_tB->Fill(mass_ttbar, weight);

        h_cosThetaStar_l->Fill(cosThetaStar_l, weight);
        if (cosThetaStar_l > 0) h_mtt_lF->Fill(mass_ttbar, weight);
        if (cosThetaStar_l < 0) h_mtt_lB->Fill(mass_ttbar, weight);

        h_deltaEta_l->Fill(deltaAbsEta_l, weight);
        if (deltaAbsEta_l > 0) h_mtt_lCF->Fill(mass_ttbar, weight);
        if (deltaAbsEta_l < 0) h_mtt_lCB->Fill(mass_ttbar, weight);

        h_deltaY_top->Fill(DeltaY_top, weight);
        if (DeltaY_top > 0) h_mtt_tCF->Fill(mass_ttbar, weight);
        if (DeltaY_top < 0) h_mtt_tCB->Fill(mass_ttbar, weight);

        h_cosTheta1->Fill(cosTheta_tl1, weight);
        h_cosTheta2->Fill(cosTheta_tl2, weight);
        h_cos1cos2->Fill(cos1cos2, weight);

        if (cosTheta_tl1 > 0) h_mtt_c1F->Fill(mass_ttbar, weight);
        if (cosTheta_tl1 < 0) h_mtt_c1B->Fill(mass_ttbar, weight);

        if (cosTheta_tl2 > 0) h_mtt_c2F->Fill(mass_ttbar, weight);
        if (cosTheta_tl2 < 0) h_mtt_c2B->Fill(mass_ttbar, weight);

        if (cos1cos2 > 0) h_mtt_c1c2F->Fill(mass_ttbar, weight);
        if (cos1cos2 < 0) h_mtt_c1c2B->Fill(mass_ttbar, weight);

        h2_mtt_cosThetaStar->Fill(mass_ttbar, cosThetaStar, weight);
        h2_mtt_deltaPhi->Fill(mass_ttbar, deltaPhi_ll, weight);
        h2_mtt_cosTheta1->Fill(mass_ttbar, cosTheta_tl1, weight);
        h2_mtt_cosTheta2->Fill(mass_ttbar, cosTheta_tl2, weight);
        h2_mtt_cos1cos2->Fill(mass_ttbar, cos1cos2, weight);
    }
    else {
        if (m_debug) cout << "Skipping reconstruction\n";
    }

    if (m_debug) cout << "Finished filling histograms\n";
    this->CleanupEvent();
}

void Analysis::CleanupEvent()
{
    if (m_debug) cout << "Cleaning up event ...\n";
    delete m_electrons;
    delete m_muons;
    delete m_jets;
    delete m_truthElectrons;
    delete m_truthMuons;
    delete m_truthBquarks;
    delete m_electron_truth_tags;
    delete m_muon_truth_tags;
    delete m_lepton_truth_tags;
    delete m_jet_truth_tags;
    
    m_hardTop = nullptr;
    m_hardTbar = nullptr;
    m_hardB = nullptr;
    m_hardBbar = nullptr;
    m_hardLepP = nullptr;
    m_hardLepM = nullptr;
    m_hardNu = nullptr;
    m_hardNuBar = nullptr;
    if (m_debug) cout << "Finished cleaning up event\n";
}

    

void Analysis::EveryEvent(double weight)
{
    // runs for every event with no event selection or cuts

    if (m_debug) cout << "Starting EveryEvent ...\n";

    h_n_electrons->Fill(b_Electron->GetEntries(), weight);
    h_n_muons->Fill(b_Muon->GetEntries(), weight);
    h_n_jets->Fill(b_Jet->GetEntries(), weight);

    if (m_debug) cout << "Fetching all jets ...\n";
    vector<TLorentzVector> p_j;
    int bTags = 0;
    for (int i = 0; i < b_Jet->GetEntries(); ++i) {
        Jet *jet = (Jet*) b_Jet->At(i);
        if (jet->BTag > 0) bTags++;
        TLorentzVector p;
        p.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        h_pT_alljets->Fill(p.Pt(), weight);
        h_eta_alljets->Fill(p.Eta(), weight);
        p_j.push_back(p);
    }

    h_n_bJets->Fill(bTags, weight);

    if (m_debug) cout << "Fetching all electrons ...\n";
    vector<TLorentzVector> p_el;
    for (int i = 0; i < b_Electron->GetEntries(); ++i) {
        Electron *electron = (Electron*) b_Electron->At(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, 0.0);
        h_pT_allel->Fill(p.Pt(), weight);
        h_eta_allel->Fill(p.Eta(), weight);
        p_el.push_back(p);
    }

    if (m_debug) cout << "Fetching all muons ...\n";
    vector<TLorentzVector> p_mu;
    for (int i = 0; i < b_Muon->GetEntries(); ++i) {
        Muon *muon = (Muon*) b_Muon->At(i);
        TLorentzVector p;
        p.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0.0);
        h_pT_allmu->Fill(p.Pt(), weight);
        h_eta_allmu->Fill(p.Eta(), weight);
        p_mu.push_back(p);
    }

    if (m_debug) cout << "Fetching missing ET ...\n";
    MissingET *missingET = (MissingET*) b_MissingET->At(0);
    TLorentzVector p_miss;
    double ET_miss = missingET->MET;

    if (m_debug) cout << "Calculating HT ...\n";
    double HT = 0;
    for (auto p : p_j) HT += p.Pt();
    for (auto p : p_el) HT += p.Pt();
    for (auto p : p_mu) HT += p.Pt();
    HT += ET_miss;
    HT = HT / 1000;
    h_HT_all->Fill(HT, weight);

    if (m_debug) cout << "Calculating pvis ...\n";
    TLorentzVector pvis(0, 0, 0, 0);
    for (auto p : p_j) pvis += p;
    for (auto p : p_el) pvis += p;
    for (auto p : p_mu) pvis += p;

    if (m_debug) cout << "Calculating mass_vis ...\n";
    double mass_vis = pvis.M();
    mass_vis = mass_vis / 1000;
    h_mass_vis_all->Fill(mass_vis, weight);

    if (m_debug) cout << "Calculating KT ...\n";
    double pT_vis = pvis.Pt();
    double KT = sqrt(mass_vis * mass_vis + pT_vis * pT_vis) + ET_miss;
    KT = KT / 1000;
    h_KT_all->Fill(KT, weight);

    if (m_debug) cout << "Finished EveryEvent.\n";
}

void Analysis::SetupInputFile()
{
    if (m_debug) cout << "Setting up single input file...\n";
    m_input = new vector< tuple<string, int> >;
    m_input->push_back(make_tuple(m_dataDirectory + m_inputFileName, 0));
    
    // processfilename, proc_id, n_proc, cross_section, uncertainty, weight
    m_processes = new vector< tuple<string, int, int, double, double, double> >;
    if (m_xSec) m_processes->push_back(make_tuple(m_dataDirectory + m_processfilename, 0, 1, -999, -999, -999));
    else m_processes->push_back(make_tuple("", 0, 1, 1.0, 0.0, 1.0));
    
    if (m_debug) cout << "Finished setting up inputfile.\n";
}

void Analysis::SetupInputFiles()
{
    m_input = new vector< tuple<string, int> >;
    m_processes = new vector< tuple<string, int, int, double, double, double> >;
    string filename;

    string E = to_string(m_energy);

    vector<string> initials = { "gg", "qq", "dd", "uu" };

    vector<string> finals;
    size_t pos = m_process.find("-tt");
    if (pos == std::string::npos) {
        cout << "error: final state doesn't contain -tt\n";
        exit(false);
    }
    string final_state = m_process.substr(pos);
    if (boost::contains(final_state, "-tt-bbllvv")) finals = { "-tt-bbeevv", "-tt-bbmumuvv", "-tt-bbemuvv", "-tt-bbmuevv" };
    else finals = { final_state };

    int proc_id = 0;

    for (auto fin : finals) {
        for (auto initial : initials) {

            // check initial state has been specified for analysis
            if (!boost::contains(m_process, initial)) continue;

            // get model name (if only QCD processes - set to Standard Model)
            string model = "";
            if (initial == "gg" or initial == "qq") model = "SM";
            else model = m_model;

            // infer intermediate EW sector bosons based on model and intial states
            string intermediates = "";
            if (initial == "uu" or initial == "dd") {
                intermediates = intermediates + "-AZ";
                if (m_model != "SM") intermediates = intermediates + "X";
                // intermediates = "-X";
            }

            // get options
            string options = "";
            options = m_options;

            // combine for file name
            filename = initial + intermediates + fin + "_" + model + "_" + E + "TeV" + "_" + m_pdf + options;

            cout << "Inputs:         " << m_dataDirectory << filename << "_*_delphes.root ...\n";

            // loop over all matching files (e.g. *_01.root and *_02.root)
            boost::filesystem::directory_iterator end_itr; // Default ctor yields past-the-end
            int nfiles = 0;

            if (m_useMassSlices) cout << "WARNING: I AM USING MASS SLICES\n";

            int end = 1;
            if (m_useMassSlices) end =  m_energy;

            for (int j = 0; j < end; ++j) {
                string range = "";
                if (m_useMassSlices) range = "_" + to_string(j) + "-" + to_string(j + 1);

                int nfiles_per_slice = 0;
                for (boost::filesystem::directory_iterator i(m_dataDirectory); i != end_itr; ++i) {

                    string file = i->path().filename().string();

                    if (i->path().extension() != ".root") continue;
                    if (!boost::contains(file, filename)) continue;
                    if (boost::contains(file, "KIN")) continue;
                    if (boost::contains(file, "NuW")) continue;
                    if (boost::contains(file, "TRN")) continue;
                    if (!boost::contains(file, "_delphes")) continue;
                    if (!boost::contains(file, range)) continue;
                    if (!boost::filesystem::is_regular_file(i->status())) continue;

                    regex reg(filename + "_[0-9]+-[0-9]+_[0-9]+_delphes");
                    if (!m_useMassSlices and regex_search(file, reg)) continue;
                    if (m_useMassSlices and !regex_search(file, reg)) continue;

                    nfiles++;
                    nfiles_per_slice++;
                    tuple<string, int> input = make_tuple(m_dataDirectory + file, proc_id);
                    m_input->push_back(input);
                    if (m_debug)
                    {
                        if (nfiles < 10) cout << "found " << nfiles << ":        " << get<0>(input);
                        else if (nfiles < 100) cout << "found " << nfiles << ":       " << get<0>(input);
                        else cout << "found " << nfiles << ":      " << get<0>(input);
                        cout << ", process: " << get<1>(input) << "\n";
                    }
                }

                if (m_useMassSlices and nfiles_per_slice == 0) {
                    cout << "no files in energy range " << range << " [TeV]\n";
                    continue;
                }
                string proc_filename = "";
                if (m_xSec) {
                    string proc_filename = m_dataDirectory + initial + intermediates + "-tt-bbllvv" + "_" + model + "_" + E + "TeV" + "_" + m_pdf + options + range + ".txt";
                    cout << "Process file:    " << proc_filename << " ...\n";
                }
                cout << "No. files:        " << nfiles << "\n";
                tuple< string, int, int, double, double, double > process = make_tuple(proc_filename, proc_id, nfiles, 1.0, 0.0, 1.0);
                m_processes->push_back(process);
                // proc_id++;
            }
        }
        std::sort(m_input->begin(), m_input->end());
    }

    // check some input files have been specified
    if (m_input->size() < 1) {
        cout << "error: no input files found.";
        exit(false);
    }

    // check all input files exist
    for (auto input : *m_input) {
        struct stat buffer;
        bool exists = stat((get<0>(input)).c_str(), &buffer) == 0;
        if (exists == false) {
            cout << "error: no " << get<0>(input) << "\n";
            exit(exists);
        }
    }
    if (m_xSec) {
        for (auto process : *m_processes) {
            struct stat buffer;
            bool exists = stat((get<0>(process)).c_str(), &buffer) == 0;
            if (exists == false) {
                cout << "error: no " << get<0>(process) << "\n";
                exit(exists);
            }
        }
    }
}

bool Analysis::SetupOutputFile()
{
    if (m_debug) cout << "Starting SetupOutputFile ...\n";

    string E = to_string(m_energy) + "TeV";
    string L = to_string(m_luminosity) + "fb-1";

    m_outputFileName = m_dataDirectory;
    if (m_inputFileName != "") m_outputFileName += m_inputFileName.substr(0,m_inputFileName.size() - 5);
    else m_outputFileName = m_dataDirectory + m_process + "_" + m_model + "_" + E + "_" + m_pdf + m_options + "_delphes";
    if (m_useMassSlices) m_outputFileName += "_sliced";
    m_outputFileName += "_" + m_reconstruction;
    m_outputFileName += "_b" + to_string(m_minBtags);
    m_outputFileName += m_tag;
    if (m_luminosity > 0) m_outputFileName += L;
    m_outputFileName += ".root";
    cout << "Output file:      " << m_outputFileName << "\n";
    m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
    return m_outputFile->IsOpen();
}


void Analysis::PostLoop()
{
    this->PrintCutflow();
    this->MakeDistributions();
    m_outputFile->Close();
    cout << "\nOUTPUT\n";
    cout << m_outputFileName << "\n";
}


void Analysis::Asymmetry(const string& name, const string& title, const string& xtitle, TH1D* h1, TH1D* h2)
{
    TH1D* h_numerator = (TH1D*) h1->Clone(name.data());
    TH1D* h_denominator = (TH1D*) h1->Clone();
    h_numerator->SetTitle(title.data());
    h_numerator->Add(h2, -1);
    h_denominator->Add(h2, 1);
    h_numerator->Divide(h_denominator);
    delete h_denominator;
    if (m_luminosity > 0) this->AsymmetryUncertainty(h_numerator, h1, h2);
    h_numerator->GetYaxis()->SetTitle(h_numerator->GetTitle());
    h_numerator->GetXaxis()->SetTitle(xtitle.data());
    h_numerator->Write();
}


void Analysis::AsymmetryUncertainty(TH1D* hA, TH1D* h1, TH1D* h2)
{
    double A, dA, N, N1, N2;
    for (int i = 1; i < hA->GetNbinsX() + 1; ++i) {
        A = hA->GetBinContent(i);
        N1 = h1->GetBinContent(i);
        N2 = h2->GetBinContent(i);
        N = N1 + N2;
        if (N > 0) dA = sqrt((1.0 - A * A) / N);
        else dA = 0;
        hA->SetBinError(i, dA);
    }
}


void Analysis::MakeHistograms()
{
    double binWidth = 0.05;
    double Emin = 0.0;
    double Emax = 13.0;
    double nbins = (Emax - Emin) / binWidth;
    cout << "Plotting range:   " << Emin << " - " << Emax << " [TeV]\n";
    
    m_mass_ttbar_truth = new vector<double>;
    m_dR_lb_truth = new vector<double>;

    h_pT_l1 = new TH1D("pT_l1", "p^{\\ell^{+}}_{\\mathrm{T}}", 200, 0.0, 2000.0);
    h_pT_l1->Sumw2();
    h_pT_l2 = new TH1D("pT_l2", "p^{\\ell^{-}}_{\\mathrm{T}}", 200, 0.0, 2000.0);
    h_pT_l2->Sumw2();
    h_pT_jets = new TH1D("pT_jet", "p_{\\mathrm{T}}^{jet}", 200, 0.0, 2000.0);
    h_pT_jets->Sumw2();
    h_pT_bjets = new TH1D("pT_bjet", "p_{\\mathrm{T}}^{b-jet}", 200, 0.0, 2000.0);
    h_pT_bjets->Sumw2();
    h_pT_qjets = new TH1D("pT_qjet", "p_{\\mathrm{T}}^{q-jet}", 200, 0.0, 2000.0);
    h_pT_qjets->Sumw2();
    h_pT_top = new TH1D("pT_top", "p_{\\mathrm{T}}^{t}", 200, 0.0, 2000.0);
    h_pT_top->Sumw2();
    h_pT_tbar = new TH1D("pT_tbar", "p_{\\mathrm{T}}^{\\bar{t}}", 200, 0.0, 2000.0);
    h_pT_tbar->Sumw2();
    h_pT_ttbar = new TH1D("pT_ttbar", "p_{\\mathrm{T}}^{t\\bar{t}}", 200, 0.0, 2000.0);
    h_pT_ttbar->Sumw2();
    h_pT_top_truth = new TH1D("pT_top_truth", "p^{truth}_{\\mathrm{T}}^{t}", 200, 0.0, 2000.0);
    h_pT_top_truth->Sumw2();
    h_pT_tbar_truth = new TH1D("pT_tbar_truth", "p^{truth}_{\\mathrm{T}}^{\\bar{t}}", 200, 0.0, 2000.0);
    h_pT_tbar_truth->Sumw2();
    h_pT_ttbar_truth = new TH1D("pT_ttbar_truth", "p^{truth}_{\\mathrm{T}}^{t\\bar{t}}", 200, 0.0, 2000.0);
    h_pT_ttbar_truth->Sumw2();

    h_eta_l1 = new TH1D("eta_l1", "\\eta_{\\ell^{+}}", 200, -5.0, 5.0);
    h_eta_l1->Sumw2();
    h_eta_l2 = new TH1D("eta_l2", "\\eta_{\\ell^{-}}", 200, -5.0, 5.0);
    h_eta_l2->Sumw2();
    h_eta_jets = new TH1D("eta_jet", "\\eta_{jet}", 200, -5.0, 5.0);
    h_eta_jets->Sumw2();
    h_eta_bjets = new TH1D("eta_bjet", "\\eta_{b-jet}", 200, -5.0, 5.0);
    h_eta_bjets->Sumw2();
    h_eta_qjets = new TH1D("eta_qjet", "\\eta_{q-jet}", 200, -5.0, 5.0);
    h_eta_qjets->Sumw2();
    h_eta_top = new TH1D("eta_top", "\\eta_{t}", 200, -5.0, 5.0);
    h_eta_top->Sumw2();
    h_eta_top_truth = new TH1D("eta_top_truth", "\\eta_{t}^{truth}", 200, -5.0, 5.0);
    h_eta_top_truth->Sumw2();
    h_eta_tbar = new TH1D("eta_tbar", "\\eta_{\\bar{t}}", 200, -5.0, 5.0);
    h_eta_tbar->Sumw2();
    h_eta_tbar_truth = new TH1D("eta_tbar_truth", "\\eta^{truth}_{\\bar{t}}", 200, -5.0, 5.0);
    h_eta_tbar_truth->Sumw2();
    h_eta_ttbar = new TH1D("eta_ttbar", "\\eta_{t\\bar{t}}", 200, -5.0, 5.0);
    h_eta_ttbar->Sumw2();
    h_eta_ttbar_truth = new TH1D("eta_ttbar_truth", "\\eta^{truth}_{t\\bar{t}}", 200, -5.0, 5.0);
    h_eta_ttbar_truth->Sumw2();
    
    h_y_top = new TH1D("y_top", "\\y_{t}", 50, -5.0, 5.0);
    h_y_top->Sumw2();
    h_y_top_truth = new TH1D("y_top_truth", "\\y_{t}^{truth}", 50, -5.0, 5.0);
    h_y_top_truth->Sumw2();
    h_y_tbar = new TH1D("y_tbar", "\\y_{\\bar{t}}", 50, -5.0, 5.0);
    h_y_tbar->Sumw2();
    h_y_tbar_truth = new TH1D("y_tbar_truth", "\\y^{truth}_{\\bar{t}}", 50, -5.0, 5.0);
    h_y_tbar_truth->Sumw2();
    h_y_ttbar = new TH1D("y_ttbar", "y_{t\\bar{t}}\\ ", 50, -2.5, 2.5);
    h_y_ttbar->Sumw2();
    h_y_ttbar_truth = new TH1D("y_ttbar_truth", "y^{truth}_{t\\bar{t}}\\ ", 50, -2.5, 2.5);
    h_y_ttbar_truth->Sumw2();

    h_phi_top = new TH1D("phi_top", "\\phi_{t}", 200, -1.0, 1.0);
    h_phi_top->Sumw2();
    h_phi_top_truth = new TH1D("phi_top_truth", "\\phi^{truth}_{t}", 200, -1.0, 1.0);
    h_phi_top_truth->Sumw2();
    h_phi_tbar = new TH1D("phi_tbar", "\\phi_{\\bar{t}}", 200, -1.0, 1.0);
    h_phi_tbar->Sumw2();
    h_phi_tbar_truth = new TH1D("phi_tbar_truth", "\\phi^{truth}_{\\bar{t}}", 200, -1.0, 1.0);
    h_phi_tbar_truth->Sumw2();
    h_phi_ttbar = new TH1D("phi_ttbar", "\\phi_{\\bar{t}}", 200, -1.0, 1.0);
    h_phi_ttbar->Sumw2();
    h_phi_ttbar_truth = new TH1D("phi_ttbar_truth", "\\phi^{truth}_{\\bar{t}}", 200, -1.0, 1.0);
    h_phi_ttbar_truth->Sumw2();

    h_mass_W1 = new TH1D("mass_W1", "m_{W^{+}}\\ ", 150, 0.0, 150.0);
    h_mass_W1->Sumw2();
    h_mass_W2 = new TH1D("mass_W2", "m_{W^{-}}\\ ", 150, 0.0, 150.0);
    h_mass_W2->Sumw2();
    h_mass_top = new TH1D("mass_top", "m_{t}\\ ", 40, 100.0, 300.0);
    h_mass_top->Sumw2();
    h_mass_top_truth = new TH1D("mass_top_truth", "m^{truth}_{t}\\ ", 40, 100.0, 300.0);
    h_mass_top_truth->Sumw2();
    h_mass_tbar = new TH1D("mass_tbar", "m_{\\bar{t}}\\ ", 40, 100.0, 300.0);
    h_mass_tbar->Sumw2();
    h_mass_tbar_truth = new TH1D("mass_tbar_truth", "m^{truth}_{\\bar{t}}\\ ", 40, 100.0, 300.0);
    h_mass_tbar_truth->Sumw2();
    h_mass_ttbar = new TH1D("mass_ttbar", "m_{t\\bar{t}}\\ ", nbins, Emin, Emax);
    h_mass_ttbar->Sumw2();
    h_mass_ttbar_truth = new TH1D("mass_ttbar_truth", "m^{truth}_{t\\bar{t}}\\ ", nbins, Emin, Emax);
    h_mass_ttbar_truth->Sumw2();
    
    h_ETmiss_truth = new TH1D("ET_miss_truth", "E^{\\mathrm{miss,truth}}_{\\mathrm{T}}\\ ", 200, 0, 2000.0);
    h_ETmiss_truth->Sumw2();
    h_ETmiss = new TH1D("ETmiss", "E^{\\mathrm{miss}}_{\\mathrm{T}}\\ ", 200, 0, 2000.0);
    h_ETmiss->Sumw2();

    h_E_top = new TH1D("E_top", "E_{t}", 100, 0.0, 5000.0);
    h_E_top->Sumw2();
    h_E_tbar = new TH1D("E_tbar", "E_{\\bar{t}}", 100, 0.0, 5000.0);
    h_E_tbar->Sumw2();

    // AtFB
    h_mtt_tF = new TH1D("mtt_tF", "m_{t\\bar{t}}^{tF}", nbins, Emin, Emax);
    h_mtt_tF->Sumw2();
    h_mtt_tB = new TH1D("mtt_tB", "m_{t\\bar{t}}^{tB}", nbins, Emin, Emax);
    h_mtt_tB->Sumw2();

    // AlFB
    h_mtt_lF = new TH1D("mtt_lF", "m_{t\\bar{t}}^{\\ell F}", nbins, Emin, Emax);
    h_mtt_lF->Sumw2();
    h_mtt_lB = new TH1D("mtt_lB", "m_{t\\bar{t}}^{\\ell B}", nbins, Emin, Emax);
    h_mtt_lB->Sumw2();

    // AlFB HT
    h_HT_lF = new TH1D("HT_lF", "H_{\\mathrm{T}}^{\\ell F}", nbins, Emin, Emax);
    h_HT_lF->Sumw2();
    h_HT_lB = new TH1D("HT_lB", "H_{\\mathrm{T}}^{\\ell B}", nbins, Emin, Emax);
    h_HT_lB->Sumw2();

    // AlFB KT
    h_KT_lF = new TH1D("KT_lF", "K_{\\mathrm{T}}^{\\ell F}", nbins, Emin, Emax);
    h_KT_lF->Sumw2();
    h_KT_lB = new TH1D("KT_lB", "K_{\\mathrm{T}}^{\\ell B}", nbins, Emin, Emax);
    h_KT_lB->Sumw2();

    // AtC
    h_mtt_tCF = new TH1D("mtt_tCF", "m_{t\\bar{t}}^{tCF}", nbins, Emin, Emax);
    h_mtt_tCF->Sumw2();
    h_mtt_tCB = new TH1D("mtt_tCB", "m_{t\\bar{t}}^{tCB}", nbins, Emin, Emax);
    h_mtt_tCB->Sumw2();

    // AlC
    h_mtt_lCF = new TH1D("mtt_lCF", "m_{t\\bar{t}}^{\\ell CF}", nbins, Emin, Emax);
    h_mtt_lCF->Sumw2();
    h_mtt_lCB = new TH1D("mtt_lCB", "m_{t\\bar{t}}^{\\ell CB}", nbins, Emin, Emax);
    h_mtt_lCB->Sumw2();

    // AlC (HT)
    h_HT_lCF = new TH1D("HT_lCF", "H_{\\mathrm{T}}^{\\ell CF}", nbins, Emin, Emax);
    h_HT_lCF->Sumw2();
    h_HT_lCB = new TH1D("HT_lCB", "H_{\\mathrm{T}}^{\\ell CB}", nbins, Emin, Emax);
    h_HT_lCB->Sumw2();

    // AlC (KT)
    h_KT_lCF = new TH1D("KT_lCF", "K_{\\mathrm{T}}^{\\ell CF}", nbins, Emin, Emax);
    h_KT_lCF->Sumw2();
    h_KT_lCB = new TH1D("KT_lCB", "K_{\\mathrm{T}}^{\\ell CB}", nbins, Emin, Emax);
    h_KT_lCB->Sumw2();

    // Ac1
    h_mtt_c1F = new TH1D("mtt_c1F", "m_{t\\bar{t}}^{c_{1}F}", nbins, Emin, Emax);
    h_mtt_c1F->Sumw2();
    h_mtt_c1B = new TH1D("mtt_c1B", "m_{t\\bar{t}}^{c_{1}B}", nbins, Emin, Emax);
    h_mtt_c1B->Sumw2();

    // Ac2
    h_mtt_c2F = new TH1D("mtt_c2F", "m_{t\\bar{t}}^{c_{2}F}", nbins, Emin, Emax);
    h_mtt_c2F->Sumw2();
    h_mtt_c2B = new TH1D("mtt_c2B", "m_{t\\bar{t}}^{c_{2}B}", nbins, Emin, Emax);
    h_mtt_c2B->Sumw2();

    // Ac1c2
    h_mtt_c1c2F = new TH1D("mtt_c1c2F", "m_{t\\bar{t}}^{c_{1}c_{2}F}", nbins, Emin, Emax);
    h_mtt_c1c2F->Sumw2();
    h_mtt_c1c2B = new TH1D("mtt_c1c2B", "m_{t\\bar{t}}^{c_{1}c_{2}B}", nbins, Emin, Emax);
    h_mtt_c1c2B->Sumw2();

    // cphi
    h_mtt_cPhiF = new TH1D("mtt_cPhiF", "m_{t\\bar{t}}^{\\cos\\phi F}", nbins, Emin, Emax);
    h_mtt_cPhiF->Sumw2();
    h_mtt_cPhiB = new TH1D("mtt_cPhiB", "m_{t\\bar{t}}^{\\cos\\phi B}", nbins, Emin, Emax);
    h_mtt_cPhiB->Sumw2();

    // ADphi
    h_mtt_DphiF = new TH1D("mtt_DphiF", "m_{t\\bar{t}}^{\\Delta\\phi F}", nbins, Emin, Emax);
    h_mtt_DphiF->Sumw2();
    h_mtt_DphiB = new TH1D("mtt_DphiB", "m_{t\\bar{t}}^{\\Delta\\phi B}", nbins, Emin, Emax);
    h_mtt_DphiB->Sumw2();

    // ADphi HT
    h_HT_DphiF = new TH1D("HT_DphiF", "H_{\\mathrm{T}}^{\\Delta\\phi F}", nbins, Emin, Emax);
    h_HT_DphiF->Sumw2();
    h_HT_DphiB = new TH1D("HT_DphiB", "H_{\\mathrm{T}}^{\\Delta\\phi B}", nbins, Emin, Emax);
    h_HT_DphiB->Sumw2();

    // ADphi KT
    h_KT_DphiF = new TH1D("KT_DphiF", "K_{\\mathrm{T}}^{\\Delta\\phi F}", nbins, Emin, Emax);
    h_KT_DphiF->Sumw2();
    h_KT_DphiB = new TH1D("KT_DphiB", "K_{\\mathrm{T}}^{\\Delta\\phi B}", nbins, Emin, Emax);
    h_KT_DphiB->Sumw2();

    h2_mtt_cosThetaStar = new TH2D("mtt_costhetastar", "m_{t\\bar{t}}\\ \\times \\cos\\theta^{*}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosThetaStar->GetXaxis()->SetTitle("m_{t\\bar{t}}");
    h2_mtt_cosThetaStar->GetYaxis()->SetTitle("\\cos\\theta^{*}");
    h2_mtt_cosThetaStar->Sumw2();
    
    h2_mtt_truth_ETmiss_truth = new TH2D("mtt_truth_ETmiss_truth", "m^\\mathrm{truth}_{t\\bar{t}}\\ \\times E^{\\mathrm{miss,truth}}_{\\mathrm{T}}", nbins, Emin, Emax, 200, 0, 2000.0);
    h2_mtt_truth_ETmiss_truth->Sumw2();

    h_cosTheta_ttbar = new TH1D("costheta_tt", "\\cos\\theta_{t,\\bar{t}}", 20, -1.0, 1.0);
    h_cosTheta_ttbar->Sumw2();

    // AtlFB
    // h_mtt_tlF = new TH1D("mtt_tlF", "m_{t\\bar{t}}^{t\\ell F}\\;", nbins, Emin, Emax);
    // h_mtt_tlF->Sumw2();
    // h_mtt_tlB = new TH1D("mtt_tlB", "m_{t\\bar{t}}^{t\\ell B}\\;", nbins, Emin, Emax);
    // h_mtt_tlB->Sumw2();

    // Aphil
    // h_mtt_philF = new TH1D("mtt_philF", "m_{t\\bar{t}}^{\\phi_{\\ell}F}\\;", nbins, Emin, Emax);
    // h_mtt_philF->Sumw2();
    // h_mtt_philB = new TH1D("mtt_philB", "m_{t\\bar{t}}^{\\phi_{\\ell}B}\\;", nbins, Emin, Emax);
    // h_mtt_philB->Sumw2();

    // AlEl
    // h_mtt_ElF = new TH1D("mtt_ElF", "m_{t\\bar{t}}^{E_{\\ell} F}\\;", nbins, Emin, Emax);
    // h_mtt_ElF->Sumw2();
    // h_mtt_ElB = new TH1D("mtt_ElB", "m_{t\\bar{t}}^{E_{\\ell} B}\\;", nbins, Emin, Emax);
    // h_mtt_ElB->Sumw2();

    h_cosTheta = new TH1D("costheta", "\\cos\\theta", 200, -1.0, 1.0);
    h_cosTheta->Sumw2();
    h_cosThetaStar = new TH1D("costheta_star", "\\cos\\theta^{*}", 200, -1.0, 1.0);
    h_cosThetaStar->Sumw2();

    h_cosTheta_l = new TH1D("costheta_l", "\\cos\\theta_{\\ell}", 200, -1.0, 1.0);
    h_cosTheta_l->Sumw2();
    h_cosThetaStar_l = new TH1D("costheta_star_l", "\\cos\\theta^{*}_{\\ell}", 200, -1.0, 1.0);
    h_cosThetaStar_l->Sumw2();

    h_deltaEta_l = new TH1D("deltaEta_l", "\\Delta |\\eta_{\\ell}|", 200, -5.0, 5.0);
    h_deltaEta_l->Sumw2();
    
    h_deltaY_top_truth = new TH1D("DeltaY_top_truth", "\\Delta |y_{t}| (truth)", 200, -5.0, 5.0);
    h_deltaY_top_truth->Sumw2();

    h_deltaY_top = new TH1D("DeltaY_top", "\\Delta |y_{t}|", 200, -5.0, 5.0);
    h_deltaY_top->Sumw2();

    h_HT_all = new TH1D("HT_all", "H^{\\mathrm{all}}_{\\mathrm{T}}", 60, 0.0, 6.0);
    h_HT_all->Sumw2();
    h_HT = new TH1D("HT", "H_{\\mathrm{T}}", 60, 0.0, 6.0);
    h_HT->Sumw2();
    h_HTmet = new TH1D("HTmet", "H_{\\mathrm{T}} + E^{\\mathrm{miss}}_{\\mathrm{T}}", 60, 0.0, 6.0);
    h_HTmet->Sumw2();
    h_HTjMET = new TH1D("HTjMET", "H^{j}_{\\mathrm{T}} + E^{\\mathrm{miss}}_{\\mathrm{T}}}", 60, 0.0, 6.0);
    h_HTjMET->Sumw2();
    h_HT_truth = new TH1D("HT_truth", "H_{\\mathrm{T}}", 60, 0.0, 6.0);
    h_HT_truth->Sumw2();
    h_HTmet_truth = new TH1D("HTmet_truth", "H_{\\mathrm{T}} + E^{\\mathrm{miss}}_{\\mathrm{T}}\\,(\\mathrm{truth})", 60, 0.0, 6.0);
    h_HTmet_truth->Sumw2();

    h_KT_all = new TH1D("KT_all", "K^{\\mathrm{all}}_{\\mathrm{T}}", 50, 0.0, 6.0);
    h_KT_all->Sumw2();
    h_KT = new TH1D("KT", "K_{\\mathrm{T}}", 60, 0.0, 6.0);
    h_KT->Sumw2();
    h_KTj = new TH1D("KTj", "K^{j}_{\\mathrm{T}}", 60, 0.0, 6.0);
    h_KTj->Sumw2();
    h_KT_truth = new TH1D("KT_truth", "K_{\\mathrm{T}}\\,(\\mathrm{truth})", 60, 0.0, 6.0);
    h_KT_truth->Sumw2();

    h_mass_vis_all = new TH1D("mvis_all", "m^{\\mathrm{\\mathrm{all}}}_{\\mathrm{vis}}", 60, 0.0, 6.0);
    h_mass_vis_all->Sumw2();
    h_mass_vis = new TH1D("mass_vis", "m_{\\mathrm{vis}}", 60, 0.0, 6.0);
    h_mass_vis->Sumw2();
    h_mass_bbll = new TH1D("mass_bbll", "m_{bbll}\\ ", 60, 0.0, 6.0);
    h_mass_bbll->Sumw2();

    h_pT_alljets = new TH1D("pT_alljets", "p_{\\mathrm{T}}^{jet}", 1000, 0.0, 5000.0);
    h_pT_alljets->Sumw2();
    h_pT_allel = new TH1D("pT_allel", "p_{\\mathrm{T}}^{e}", 1000, 0.0, 5000.0);
    h_pT_allel->Sumw2();
    h_pT_allmu = new TH1D("pT_allmu", "p_{\\mathrm{T}}^{\\mu}", 1000, 0.0, 5000.0);
    h_pT_allmu->Sumw2();
    h_eta_alljets = new TH1D("eta_alljets", "\\eta^{jet}", 200, -5.0, 5.0);
    h_eta_alljets->Sumw2();
    h_eta_allel = new TH1D("eta_allel", "\\eta^{e}", 200, -5.0, 5.0);
    h_eta_allel->Sumw2();
    h_eta_allmu = new TH1D("eta_allmu", "\\eta^{\\mu}", 200, -5.0, 5.0);
    h_eta_allmu->Sumw2();

    h_cosPhi = new TH1D("cosPhi", "\\cos\\varphi", 10, -1.0, 1.0);
    h_cosPhi->Sumw2();

    h_deltaPhi_ll = new TH1D("delta_phi", "\\Delta\\phi_{\\ell^{+}\\ell^{-}}", 10, 0.0, 1.0);
    h_deltaPhi_ll->Sumw2();
    
    h_deltaPhi_tt_truth = new TH1D("deltaPhi_tt_truth", "\\Delta\\phi_{t\\bar{t}} (truth)", 10, -m_pi, m_pi);
    h_deltaPhi_tt_truth->Sumw2();
    
    h_deltaPhi_tt = new TH1D("deltaPhi_tt", "\\Delta\\phi_{t\\bar{t}}", 10, -m_pi, m_pi);
    h_deltaPhi_tt->Sumw2();
    
    h_cosTheta1_truth = new TH1D("cosTheta_tl1_truth", "\\cos\\theta^{t}_{\\ell^{+}} (truth)", 50, -1.0, 1.0);
    h_cosTheta1_truth->Sumw2();
    h_cosTheta2_truth = new TH1D("cosTheta_tl2_truth", "\\cos\\theta^{\\bar{t}}_{\\ell^{-}} (truth)", 50, -1.0, 1.0);
    h_cosTheta2_truth->Sumw2();
    // h_cos1cos2_truth = new TH1D("cos1cos2_truth", "\\cos\\theta^{t}_{\\ell^{+}}\\cos\\theta^{\\bar{t}}_{\\ell^{-}} (truth)", 10, -1.0, 1.0);
    // h_cos1cos2_truth->Sumw2();

    h_cosTheta1 = new TH1D("cosTheta_tl1", "\\cos\\theta^{t}_{\\ell^{+}}", 10, -1.0, 1.0);
    h_cosTheta1->Sumw2();
    h_cosTheta2 = new TH1D("cosTheta_tl2", "\\cos\\theta^{\\bar{t}}_{\\ell^{-}}", 10, -1.0, 1.0);
    h_cosTheta2->Sumw2();
    h_cos1cos2 = new TH1D("cos1cos2", "\\cos\\theta^{t}_{\\ell^{+}}\\cos\\theta^{\\bar{t}}_{\\ell^{-}}", 10, -1.0, 1.0);
    h_cos1cos2->Sumw2();

    h_n_truthElectrons = new TH1D("n_truth_electrons", "n_{e}\\ ", 10, 0.0, 10.0);
    h_n_truthElectrons->Sumw2();
    h_n_truthMuons = new TH1D("n_truth_muons", "n_{\\mu}\\ ", 10, 0.0, 10.0);
    h_n_truthMuons->Sumw2();
    h_n_truthBquarks = new TH1D("n_truth_bQuarks", "n_{b}\\ ", 10, 0.0, 10.0);
    h_n_truthBquarks->Sumw2();

    h_n_electrons = new TH1D("n_electrons", "n_{e}\\ ", 10, 0.0, 10.0);
    h_n_electrons_truth_tagged = new TH1D("n_electrons_truth_tagged", "n_{e}\\,\\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    h_n_electrons->Sumw2();
    h_n_muons = new TH1D("n_muons", "n_{\\mu}\\ ", 10, 0.0, 10.0);
    h_n_muons->Sumw2();
    h_n_muons_truth_tagged = new TH1D("n_muons_truth_tagged", "n_{\\mu}\\, \\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    h_n_jets = new TH1D("n_jets", "n_{jet}\\ ", 10, 0.0, 10.0);
    h_n_jets->Sumw2();
    h_n_bJets = new TH1D("n_bJets", "n_{b-jet}\\ ", 10, 0.0, 10.0);
    h_n_bJets->Sumw2();
    h_n_jets_truth_tagged = new TH1D("n_jets_truth_tagged", "n_{jet}\\,\\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    h_ntracks_truth_tagged_jets = new TH1D("ntracks_truth_tagged_jets", "n_{tracks}\\,\\mathrm{(truth\\, tagged\\, jets)}\\ ", 50, 0.0, 50.0);
    
    h_n_selElectrons = new TH1D("n_sel_electrons", "n_{e}\\ ", 10, 0.0, 10.0);
    h_n_selElectrons->Sumw2();
    h_n_selElectrons_truth_tagged = new TH1D("n_sel_electrons_truth_tagged", "n_{e}\\,\\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    h_n_selElectrons_truth_tagged->SetNdivisions(4);
    h_n_selMuons = new TH1D("n_sel_muons", "n_{\\mu}\\ ", 10, 0.0, 10.0);
    h_n_selMuons->Sumw2();
    h_n_selMuons_truth_tagged = new TH1D("n_sel_muons_truth_tagged", "n_{\\mu}\\, \\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    h_n_selJets = new TH1D("n_sel_jets", "n_{jet}\\ ", 10, 0.0, 10.0);
    h_n_selJets->Sumw2();
    h_n_selJets_truth_tagged = new TH1D("n_sel_jets_truth_tagged", "n_{jet}\\,\\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    
    h_n_uniqueElectrons = new TH1D("n_unique_electrons", "n_{e}\\ ", 10, 0.0, 10.0);
    h_n_uniqueElectrons->Sumw2();
    h_n_uniqueElectrons_truth_tagged = new TH1D("n_unique_electrons_truth_tagged", "n_{e}\\,\\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    h_n_uniqueMuons = new TH1D("n_unique_muons", "n_{\\mu}\\ ", 10, 0.0, 10.0);
    h_n_uniqueMuons->Sumw2();
    h_n_uniqueMuons_truth_tagged = new TH1D("n_unique_muons_truth_tagged", "n_{\\mu}\\, \\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);
    h_n_uniqueJets = new TH1D("n_unique_jets", "n_{jet}\\ ", 10, 0.0, 10.0);
    h_n_uniqueJets->Sumw2();
    h_n_uniqueJets_truth_tagged = new TH1D("n_unique_jets_truth_tagged", "n_{jet}\\,\\mathrm{(truth\\, tagged)}\\ ", 4, 0.0, 4.0);

    h_n_passElectrons = new TH1D("n_pass_electrons", "n_{e}\\ ", 10, 0.0, 10.0);
    h_n_passElectrons->Sumw2();
    h_n_passMuons = new TH1D("n_pass_muons", "n_{\\mu}\\ ", 10, 0.0, 10.0);
    h_n_passMuons->Sumw2();
    h_n_passJets = new TH1D("n_pass_jets", "n_{jet}\\ ", 10, 0.0, 10.0);
    h_n_passJets->Sumw2();
    h_n_passBjets = new TH1D("n_pass_bJets", "n_{b-jet}\\ ", 10, 0.0, 10.0);
    h_n_passBjets->Sumw2();

    h2_mtt_deltaPhi = new TH2D("mtt_deltaphi", "m_{t\\bar{t}}\\ \\times \\Delta\\phi_{\\ell^{+}\\ell^{-}}", nbins, Emin, Emax, 10, 0.0, 1.0);
    h2_mtt_deltaPhi->GetXaxis()->SetTitle("m_{t\\bar{t}}\\;");
    h2_mtt_deltaPhi->GetYaxis()->SetTitle("\\Delta\\phi_{\\ell^{+}\\ell^{-}}\\;");
    h2_mtt_deltaPhi->Sumw2();

    h2_mtt_cosTheta1 = new TH2D("mtt_costheta_tl1", "m_{t\\bar{t}}\\ \\times \\cos\\theta^{t}_{\\ell^{+}}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosTheta1->GetXaxis()->SetTitle("m_{t\\bar{t}}\\;[TeV]");
    h2_mtt_cosTheta1->GetYaxis()->SetTitle("\\cos\\theta_{\\ell^{+}}");
    h2_mtt_cosTheta1->Sumw2();

    h2_mtt_cosTheta2 = new TH2D("mtt_costheta_tl2", "m_{t\\bar{t}}\\ \\times \\cos\\theta^{\\bar{t}}_{\\ell^{-}}", nbins, Emin, Emax, 10, -1.0, 1.0);
    h2_mtt_cosTheta2->GetXaxis()->SetTitle("m_{t\\bar{t}}\\;[TeV]");
    h2_mtt_cosTheta2->GetYaxis()->SetTitle("\\cos\\theta_{\\ell^{-}}");
    h2_mtt_cosTheta2->Sumw2();

    h2_mtt_cos1cos2 = new TH2D("mtt_cos1cos2", "m_{t\\bar{t}}\\ \\times \\cos\\theta^{t}_{\\ell^{+}}\\cos\\theta^{\\bar{t}}_{\\ell^{-}}", nbins, Emin, Emax, 20, -1.0, 1.0);
    h2_mtt_cos1cos2->GetXaxis()->SetTitle("m_{t\\bar{t}}");
    h2_mtt_cos1cos2->GetYaxis()->SetTitle("\\cos\\theta_{\\ell^{+}}\\cos\\theta_{\\ell^{-}}");
    h2_mtt_cos1cos2->Sumw2();
    
    h2_HT_deltaPhi = new TH2D("HT_deltaphi", "H_{\\mathrm{T}} \\times \\Delta\\phi_{\\ell^{+}\\ell^{-}}", 40, 0.0, 4.0, 10, 0.0, 1.0);
    h2_HT_deltaPhi->GetXaxis()->SetTitle("H_{\\mathrm{T}}\\;[TeV]");
    h2_HT_deltaPhi->GetYaxis()->SetTitle("\\Delta\\phi_{\\ell^{+}\\ell^{-}}\\; [rad/\\pi]");
    h2_HT_deltaPhi->Sumw2();
    
    h2_HTmet_deltaPhi = new TH2D("HTmet_deltaphi", "(H_{\\mathrm{T}} + E^{\\mathrm{miss}}_{\\mathrm{T}}) \\times \\Delta\\phi_{\\ell^{+}\\ell^{-}}", 40, 0.0, 4.0, 10, 0.0, 1.0);
    h2_HTmet_deltaPhi->GetXaxis()->SetTitle("H_{\\mathrm{T}} + E^{\\mathrm{miss}}_{\\mathrm{T}}\\;[TeV]");
    h2_HTmet_deltaPhi->GetYaxis()->SetTitle("\\Delta\\phi_{\\ell^{+}\\ell^{-}}\\; [rad/\\pi]");
    h2_HTmet_deltaPhi->Sumw2();

    h2_mvis_deltaPhi = new TH2D("mvis_deltaphi", "m_{\\mathrm{vis}} \\times \\Delta\\phi_{\\ell^{+}\\ell^{-}}", 40, 0.0, 4.0, 10, 0.0, 1.0);
    h2_mvis_deltaPhi->GetXaxis()->SetTitle("m_{\\mathrm{vis}}\\;[TeV]");
    h2_mvis_deltaPhi->GetYaxis()->SetTitle("\\Delta\\phi_{\\ell^{+}\\ell^{-}}\\;[rad/\\pi]");
    h2_mvis_deltaPhi->Sumw2();

    h2_KT_deltaPhi = new TH2D("KT_deltaphi", "K_{\\mathrm{T}} \\times \\Delta\\phi_{\\ell^{+}\\ell^{-}}\\;", 40, 0.0, 4.0, 10, 0.0, 1.0);
    h2_KT_deltaPhi->GetXaxis()->SetTitle("K_{\\mathrm{T}}\\;[TeV]\\;");
    h2_KT_deltaPhi->GetYaxis()->SetTitle("\\Delta\\phi_{\\ell^{+}\\ell^{-}}\\;[rad/\\pi]");
    h2_KT_deltaPhi->Sumw2();

    h_dR_top = new TH1D("dR_top", "\\Delta R\\left(t_{\\mathrm{truth}},t_{\\mathrm{reco}}\\right)", 100, 0.0, 10.0);
    h_dR_tbar = new TH1D("dR_tbar", "\\Delta R\\left(\\bar{t}_{\\mathrm{truth}},\\bar{t}_{\\mathrm{reco}}\\right)", 100, 0.0, 10.0);
    h_dR_ttbar = new TH1D("dR_ttbar", "\\Delta R\\left(t\\bar{t}_{\\mathrm{truth}}, t\\bar{t}_{\\mathrm{reco}}\\right)", 100, 0.0, 10.0);
    
    h_dR_l = new TH1D("dR_l", "\\Delta R\\left(\\ell_{\\mathrm{reco}}, \\ell_{\\mathrm{truth}}\\right)", 20, 0.0, 0.2);
    
    h_dR_lb_truth = new TH1D("dR_lb_truth", "\\Delta R\\left(\\ell_{\\mathrm{truth}},b_{\\mathrm{truth}}\\right)", 700, 0.0, 7.0);
    // h_dR_l1b1_truth = new TH1D("dR_l1b1_truth", "\\Delta R\\left(\\ell_{\\mathrm{truth}},b_{\\mathrm{truth}}\\right)", 700, 0.0, 7.0);
    // h_dR_l2b2_truth = new TH1D("dR_l2b2_truth", "\\Delta R\\left(\\ell^{-}_{\\mathrm{truth}},\\bar{b}_{\\mathrm{truth}}\\right)", 700, 0.0, 7.0);
    // h_dR_t1l1_truth = new TH1D("dR_t1l1_truth", "\\Delta R\\left(t_{\\mathrm{truth}},\\ell_{\\mathrm{truth}}\\right)", 100, 0.0, 7.0);
    h_dR_t1l1_truth = new TH1D("dR_t1l1_truth", "\\Delta R\\left(t_{\\mathrm{truth}},\\ell^{+}_{\\mathrm{truth}}\\right)", 100, 0.0, 7.0);
    h_dR_t2l2_truth = new TH1D("dR_t2l2_truth", "\\Delta R\\left(\\bar{t}_{\\mathrm{truth}},\\ell^{-}_{\\mathrm{truth}}\\right)", 100, 0.0, 7.0);
    // h_dR_t1b1_truth = new TH1D("dR_t1b1_truth", "\\Delta R\\left(t_{\\mathrm{truth}}, b_{\\mathrm{truth}}\\right)", 100, 0.0, 7.0);
    h_dR_t1b1_truth = new TH1D("dR_t1b1_truth", "\\Delta R\\left(t_{\\mathrm{truth}}, b_{\\mathrm{truth}}\\right)", 100, 0.0, 7.0);
    h_dR_t2b2_truth = new TH1D("dR_t2b2_truth", "\\Delta R\\left(\\bar{t}_{\\mathrm{truth}}, \\bar{b}_{\\mathrm{truth}}\\right)", 100, 0.0, 7.0);
    h_dR_t1t2_truth = new TH1D("dR_t1t2_truth", "\\Delta R\\left(t_{\\mathrm{truth}}, \\bar{t}_{\\mathrm{truth}}\\right)", 100, 0.0, 7.0);
    
    h2_dR_lb_mtt_truth = new TH2D("dR_lb_mtt_truth", " m_{tt}^{truth} \\times \\Delta R(\\ell, b)", nbins, Emin, Emax, 7000, 0.0, 7.0);
    // h2_dR_l1b1_mtt_truth = new TH2D("dR_l1b1_mtt_truth", "dR_l1b1_mtt_truth", nbins, Emin, Emax, 100, 0.0, 7.0);
    // h2_dR_l2b2_mtt_truth = new TH2D("dR_l2b2_mtt_truth", "dR_l2b2_mtt_truth", nbins, Emin, Emax, 100, 0.0, 7.0);
    h2_dR_t1l1_mtt_truth = new TH2D("dR_t1l1_mtt_truth", "dR_t1l1_mtt_truth", nbins, Emin, Emax, 700, 0.0, 7.0);
    h2_dR_t2l2_mtt_truth = new TH2D("dR_t2l2_mtt_truth", "dR_t2l2_mtt_truth", nbins, Emin, Emax, 700, 0.0, 7.0);
    h2_dR_t1b1_mtt_truth = new TH2D("dR_t1b1_mtt_truth", "dR_t1b1_mtt_truth", nbins, Emin, Emax, 700, 0.0, 7.0);
    h2_dR_t2b2_mtt_truth = new TH2D("dR_t2b2_mtt_truth", "dR_t2b2_mtt_truth", nbins, Emin, Emax, 700, 0.0, 7.0);
    h2_dR_l1b1_pTl_truth = new TH2D("dR_l1b1_pTl_truth", "dR_l1b1_pTl_truth", 200, 0, 1000, 700, 0.0, 7.0);
    h2_dR_l2b2_pTl_truth = new TH2D("dR_l2b2_pTl_truth", "dR_l2b2_pTl_truth", 200, 0, 1000, 700, 0.0, 7.0);

    h_perf_mass_top = new TH1D("perf_mass_top", "\\left(m^{\\mathrm{truth}}_{t} - m^{\\mathrm{reco}}_{t}\\right) / m^{\\mathrm{truth}}_{t}\\ ", 100, -3, 3);
    h_perf_pT_top = new TH1D("perf_pT_top", "\\left(p^{t,\\mathrm{truth}}_{T} - p^{t,\\mathrm{reco}}_{T}\\right) / p^{t,\\mathrm{truth}}_{T}\\ ", 100, -3, 3);
    h_perf_eta_top = new TH1D("perf_eta_top", "\\left(\\eta^{\\mathrm{truth}}_{t} - \\eta^{\\mathrm{reco}}_{t}\\right) / \\eta^{\\mathrm{truth}}_{t}}\\ ", 100, -3, 3);
    h_perf_phi_top = new TH1D("perf_phi_top", "\\left(\\phi^{\\mathrm{truth}}_{t} - \\phi^{\\mathrm{reco}}_{t}\\right) / \\phi^{\\mathrm{truth}}_{t}\\ ", 100, -3, 3);
    
    h_perf_mass_tbar = new TH1D("perf_mass_tbar", "\\left(m^{\\mathrm{truth}}_{\\bar{t}} - m^{\\mathrm{reco}}_{\\bar{t}}\\right) / m^{\\mathrm{truth}}_{\\bar{t}}\\ ", 100, -3, 3);
    h_perf_pT_tbar = new TH1D("perf_pT_tbar", "\\left(p^{\\bar{t},\\mathrm{truth}}_{T} - p^{\\bar{t},\\mathrm{reco}}_{T}\\right) / p^{\\bar{t},\\mathrm{truth}}_{T}\\ ", 100, -3, 3);
    h_perf_eta_tbar = new TH1D("perf_eta_tbar", "\\left(\\eta^{\\mathrm{truth}}_{\\bar{t}} - \\eta^{\\mathrm{reco}}_{\\bar{t}}\\right) / \\eta^{\\mathrm{truth}}_{\\bar{t}}}\\ ", 100, -3, 3);
    h_perf_phi_tbar = new TH1D("perf_phi_tbar", "\\left(\\phi^{\\mathrm{truth}}_{\\bar{t}} - \\phi^{\\mathrm{reco}}_{\\bar{t}}\\right) / \\phi^{\\mathrm{truth}}_{\\bar{t}}\\ ", 100, -3, 3);
    
    h_perf_mass_ttbar = new TH1D("perf_mass_ttbar", "\\left(m^{\\mathrm{truth}}_{t\\bar{t}} - m^{\\mathrm{reco}}_{t\\bar{t}}\\right) / m^{\\mathrm{truth}}_{t\\bar{t}}\\ ", 100, -3, 3);
    h_perf_pT_ttbar = new TH1D("perf_pT_ttbar", "\\left(p^{t\\bar{t},\\mathrm{truth}}_{T} - p^{t\\bar{t},\\mathrm{reco}}_{T}\\right) / p^{t\\bar{t},\\mathrm{truth}}_{T}\\ ", 100, -3, 3);
    h_perf_eta_ttbar = new TH1D("perf_eta_ttbar", "\\left(\\eta^{\\mathrm{truth}}_{t\\bar{t}} - \\eta^{\\mathrm{reco}}_{t\\bar{t}}\\right) / \\eta^{\\mathrm{truth}}_{t\\bar{t}}}\\ ", 100, -3, 3);
    h_perf_phi_ttbar = new TH1D("perf_phi_ttbar", "\\left(\\phi^{\\mathrm{truth}}_{t\\bar{t}} - \\phi^{\\mathrm{reco}}_{t\\bar{t}}\\right) / \\phi^{\\mathrm{truth}}_{t\\bar{t}}\\ ", 100, -3, 3);
    
    h_res_pT_bjets = new TH1D("res_bjet", "p^{b,\\mathrm{truth}}_{\\mathrm{T}} - p^{b,\\mathrm{reco}}_{\\mathrm{T}} (high pT)", 100, -3.0, 3.0);
    
    double pT_bins[26] = {0, 40, 50, 60, 70, 80, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000}; 
    h2_resPtBjets_pT = new TH2D("res_bjet_pT", "p^{b,\\mathrm{truth}}_{\\mathrm{T}} - p^{b,\\mathrm{reco}}_{\\mathrm{T}} \\times p_{\\mathrm{T}}", 25, pT_bins, 100, -3, 3);
    // h2_resPtBjets_pT = new TH2D("res_bjet_pT", "p^{b,\\mathrm{truth}}_{\\mathrm{T}} - p^{b,\\mathrm{reco}}_{\\mathrm{T}} \\times p_{\\mathrm{T}}", 10, 0.0, 2000.0, 100, -3.0, 3.0);
    
    h_res_pT_tops = new TH1D("res_pT_tops", "p_{\\mathrm{T}}^{t}\\ resolution", 1000, -2000.0, 2000.0);
    h_res_y_tops = new TH1D("res_y_tops", "y_{t}\\ resolution", 1000, -5.0, 5.0);
    h_res_phi_tops = new TH1D("res_phi_tops", "\\phi_{t}\\ resolution", 1000, -10.0, 10.0);
    h_res_eta_tops = new TH1D("res_eta_tops", "\\eta_{t}\\ resolution", 1000, -5.0, 5.0);
    
    h_res_pT_ttbar = new TH1D("res_pT_ttbar", "p_{\\mathrm{T}}^{t\\bar{t}}\\ resolution", 1000, -2000.0, 2000.0);
    h_res_y_ttbar = new TH1D("res_y_ttbar", "y_{t\\bar{t}}\\ resolution", 1000, -5.0, 5.0);
    h_res_mass_ttbar = new TH1D("res_mass_ttbar", "\\m_{t\\bar{t}}\\ resolution", 1000, -5000.0, 5000.0);
    h_res_eta_ttbar = new TH1D("res_eta_ttbar", "\\eta_{t\\bar{t}}\\ resolution", 1000, -5.0, 5.0);
    
    h_reco_quality = new TProfile("reco_quality", "Reconstruction Quality", 1, 0, 1);

    h_eff_cut_2l_mass_ttbar_truth = new TProfile("eff_cut_2l_mass_ttbar_truth", "eff_cut_2l_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cut_oc_mass_ttbar_truth  = new TProfile("eff_cut_oc_mass_ttbar_truth", "eff_cut_oc_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cut_mll_mass_ttbar_truth  = new TProfile("eff_cut_mll_mass_ttbar_truth", "eff_cut_mll_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cut_mZ_mass_ttbar_truth  = new TProfile("eff_cut_mZ_mass_ttbar_truth", "eff_cut_mZ_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cut_ETmiss_mass_ttbar_truth = new TProfile("eff_cut_ETmiss_mass_ttbar_truth", "eff_cut_ETmiss_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cut_HT_mass_ttbar_truth = new TProfile("eff_cut_HT_mass_ttbar_truth", "eff_cut_HT_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cut_2j_mass_ttbar_truth  = new TProfile("eff_cut_2j_mass_ttbar_truth", "eff_cut_2j_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cut_2b_mass_ttbar_truth = new TProfile("eff_cut_2b_mass_ttbar_truth", "eff_cut_2b_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cuts_mass_ttbar_truth = new TProfile("eff_cuts_mass_ttbar_truth", "eff_cuts_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_cuts_pT_top_truth = new TProfile("eff_cuts_pT_top_truth", "eff_cuts_pT_top_truth", 200, 0, 2000);
    h_eff_cuts_pT_tbar_truth = new TProfile("eff_cuts_pT_tbar_truth", "eff_cuts_pT_tbar_truth", 200, 0, 2000);
    h_eff_reco_mass_ttbar_truth = new TProfile("eff_reco_mass_ttbar_truth", "eff_reco_mass_ttbar_truth", nbins, Emin, Emax);
    h_eff_reco_pT_top_truth = new TProfile("eff_reco_pT_top_truth", "eff_reco_pT_top_truth", 200, 0, 2000);
    h_eff_reco_pT_tbar_truth = new TProfile("eff_reco_pT_tbar_truth", "eff_reco_pT_tbar_truth", 200, 0, 2000);
    
    vector<string> purity_titles = {"cuts","isolation", "overlap removal"};
    int n_pure = purity_titles.size();
    h_lepton_purity = new TProfile("lepton_purity ", "Lepton purity", n_pure, 0.0, n_pure);
    h_jet_purity = new TProfile("jet_purity", "Jet purity", n_pure, 0.0, n_pure);
    for (int i = 0; i < n_pure; ++i) {
        h_lepton_purity->GetXaxis()->SetBinLabel(i + 1, purity_titles[i].data());
        h_jet_purity->GetXaxis()->SetBinLabel(i + 1, purity_titles[i].data());
    }

    h2_perf_mass_ttbar = new TH2D("perf2_mass_ttbar", "perf2_mass_ttbar", nbins, Emin, Emax, 100, -3, 3);
    h2_perf_mass_ttbar_pTtop = new TH2D("perf2_mass_ttbar_pTtop", "perf2_mass_ttbar_pTtop", 200, 0, 2000, 100, -3, 3);
    h2_perf_mass_ttbar_pTtbar = new TH2D("perf2_mass_ttbar_pTtbar", "perf2_mass_ttbar_pTtbar", 200, 0, 2000, 100, -3, 3);

    h2_mass_top_TvR = new TH2D("mass_top_TvR", "mass_top_TvR", 40, 100.0, 300.0, 40, 100.0, 300.0);
    h2_pT_top_TvR = new TH2D("pT_top_TvR", "pT_top_TvR", 200, 0.0, 2000.0, 200, 0.0, 2000.0);
    h2_eta_top_TvR = new TH2D("eta_top_TvR", "eta_top_TvR", 200, -5.0, 5.0, 200, -5.0, 5.0);
    h2_phi_top_TvR = new TH2D("phi_top_TvR", "phi_top_TvR", 200, -m_pi, m_pi, 200, -m_pi, m_pi);

    h2_mass_tbar_TvR = new TH2D("mass_tbar_TvR", "mass_tbar_TvR", 40, 100.0, 300.0, 40, 100.0, 300.0);
    h2_pT_tbar_TvR = new TH2D("pT_tbar_TvR", "pT_tbar_TvR", 200, 0.0, 2000.0, 200, 0.0, 2000.0);
    h2_eta_tbar_TvR = new TH2D("eta_tbar_TvR", "eta_tbar_TvR", 200, -5.0, 5.0, 200, -5.0, 5.0);
    h2_phi_tbar_TvR = new TH2D("phi_tbar_TvR", "phi_tbar_TvR", 200, -m_pi, m_pi, 200, -m_pi, m_pi);

    h2_mass_ttbar_TvR = new TH2D("mass_ttbar_TvR", "mass_ttbar_TvR", nbins, Emin, Emax, nbins, Emin, Emax);
    h2_pT_ttbar_TvR = new TH2D("pT_ttbar_TvR", "pT_ttbar_TvR", nbins, Emin, Emax, nbins, Emin, Emax);
    h2_eta_ttbar_TvR = new TH2D("eta_ttbar_TvR", "eta_ttbar_TvR", 200, -5.0, 5.0, 200, -5.0, 5.0);
    h2_phi_ttbar_TvR = new TH2D("phi_ttbar_TvR", "phi_ttbar_TvR", 200, -m_pi, m_pi, 200, -m_pi, m_pi);
    
    if (m_debug) cout << "Created histograms\n";
}


void Analysis::MakeDistributions()
{
    if (m_debug) cout << "Making distributions ...\n";
    
    h_lepton_purity->Write();
    h_jet_purity->Write();
    
    // int n = m_mass_ttbar_truth->size();
    // double x[n];
    // double y[n];
    // for (int i = 0; i < n; ++i)
    // {
    //     x[i] = m_mass_ttbar_truth->at(i);
    //     y[i] = m_dR_lb_truth->at(i);
    // }
    // 
    // TGraph* g = new TGraph(n, x, y);
    // g->Write();
    // delete m_mass_ttbar_truth;
    // delete m_dR_lb_truth;

    this->WriteEfficiency(h_eff_cut_2l_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Two leptons}");
    this->WriteEfficiency(h_eff_cut_oc_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Opposite Charge}");
    this->WriteEfficiency(h_eff_cut_mll_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Sufficient}\\; m_{\\ell\\ell}");
    this->WriteEfficiency(h_eff_cut_mZ_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Outside}\\; m_{Z}\\; window");
    this->WriteEfficiency(h_eff_cut_ETmiss_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Sufficient}\\; E^{\\mathrm{miss}}_{\\mathrm{T}}");
    this->WriteEfficiency(h_eff_cut_HT_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]","\\mathrm{Sufficient}\\; H_{\\mathrm{T}}");
    this->WriteEfficiency(h_eff_cut_2j_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Sufficient jets}");
    this->WriteEfficiency(h_eff_cut_2b_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Sufficient}\\; b\\mathrm{-tags}");
    this->WriteEfficiency(h_eff_cuts_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{All cuts}");
    this->WriteEfficiency(h_eff_reco_mass_ttbar_truth, "m_{t\\bar{t}}\\, [\\mathrm{TeV}]", "\\mathrm{Top reco}");
    
    this->WriteEfficiency(h_eff_cuts_pT_top_truth, "P_{T}^{t}\\, [\\mathrm{GeV}]", "\\mathrm{All cuts}");
    this->WriteEfficiency(h_eff_cuts_pT_tbar_truth, "P_{T}^{\\bar{t}}\\, [\\mathrm{GeV}]", "\\mathrm{All cuts}");
    this->WriteEfficiency(h_eff_reco_pT_top_truth, "P_{T}^{t}\\, [\\mathrm{GeV}]", "\\mathrm{Top reco}");
    this->WriteEfficiency(h_eff_reco_pT_tbar_truth, "P_{T}^{\\bar{t}}\\, [\\mathrm{GeV}]", "\\mathrm{Top reco}");
    

    // count
    this->MakeDistribution1D(h_n_truthElectrons, "");
    this->MakeDistribution1D(h_n_truthMuons, "");
    this->MakeDistribution1D(h_n_truthBquarks, "");
    this->MakeDistribution1D(h_n_electrons, "");
    this->MakeDistribution1D(h_n_electrons_truth_tagged, "");
    this->MakeDistribution1D(h_n_muons, "");
    this->MakeDistribution1D(h_n_muons_truth_tagged, "");
    this->MakeDistribution1D(h_n_jets, "");
    this->MakeDistribution1D(h_n_jets_truth_tagged, "");
    this->MakeDistribution1D(h_ntracks_truth_tagged_jets, "");
    this->MakeDistribution1D(h_n_bJets, "");
    this->MakeDistribution1D(h_n_selElectrons, "");
    this->MakeDistribution1D(h_n_selElectrons_truth_tagged, "");
    this->MakeDistribution1D(h_n_selMuons, "");
    this->MakeDistribution1D(h_n_selMuons_truth_tagged, "");
    this->MakeDistribution1D(h_n_selJets, "");
    this->MakeDistribution1D(h_n_selJets_truth_tagged, "");
    this->MakeDistribution1D(h_n_uniqueElectrons, "");
    this->MakeDistribution1D(h_n_uniqueElectrons_truth_tagged, "");
    this->MakeDistribution1D(h_n_uniqueMuons, "");
    this->MakeDistribution1D(h_n_uniqueMuons_truth_tagged, "");
    this->MakeDistribution1D(h_n_uniqueJets, "");
    this->MakeDistribution1D(h_n_uniqueJets_truth_tagged, "");
    this->MakeDistribution1D(h_n_passElectrons, "");
    this->MakeDistribution1D(h_n_passMuons, "");
    this->MakeDistribution1D(h_n_passJets, "");
    this->MakeDistribution1D(h_n_passBjets, "");

    // pT
    this->MakeDistribution1D(h_pT_allel, "GeV");
    this->MakeDistribution1D(h_pT_allmu, "GeV");
    this->MakeDistribution1D(h_pT_l1, "TeV");
    this->MakeDistribution1D(h_pT_l2, "TeV");
    this->MakeDistribution1D(h_pT_alljets, "GeV");
    this->MakeDistribution1D(h_pT_jets, "TeV");
    this->MakeDistribution1D(h_pT_bjets, "TeV");
    this->MakeDistribution1D(h_pT_qjets, "TeV");
    this->MakeDistribution1D(h_ETmiss_truth, "TeV");
    this->MakeDistribution1D(h_ETmiss, "GeV");

    // pseudorapidity
    this->MakeDistribution1D(h_eta_allel, "");
    this->MakeDistribution1D(h_eta_allmu, "");
    this->MakeDistribution1D(h_eta_l1, "");
    this->MakeDistribution1D(h_eta_l2, "");
    this->MakeDistribution1D(h_eta_alljets, "");
    this->MakeDistribution1D(h_eta_jets, "");
    this->MakeDistribution1D(h_eta_bjets, "");
    this->MakeDistribution1D(h_eta_qjets, "");
    this->MakeDistribution1D(h_deltaEta_l, "");

    // azimuthal angle
    this->MakeDistribution1D(h_deltaPhi_ll, "rad / \\pi");
    this->MakeDistribution1D(h_deltaPhi_tt, "rad");
    this->MakeDistribution1D(h_deltaPhi_tt_truth, "rad");
    
    this->MakeDistribution1D(h_pT_top_truth, "GeV");
    this->MakeDistribution1D(h_pT_tbar_truth, "GeV");
    this->MakeDistribution1D(h_pT_ttbar_truth, "GeV");
    this->MakeDistribution1D(h_y_top_truth, "");
    this->MakeDistribution1D(h_y_tbar_truth, "");
    this->MakeDistribution1D(h_y_ttbar_truth, "");
    this->MakeDistribution1D(h_eta_top_truth, "");
    this->MakeDistribution1D(h_eta_tbar_truth, "");
    this->MakeDistribution1D(h_eta_ttbar_truth, "");
    this->MakeDistribution1D(h_phi_top_truth, "");
    this->MakeDistribution1D(h_phi_tbar_truth, "");
    this->MakeDistribution1D(h_phi_ttbar_truth, "");
    this->MakeDistribution1D(h_mass_top_truth, "GeV");
    this->MakeDistribution1D(h_mass_tbar_truth, "GeV");
    this->MakeDistribution1D(h_mass_ttbar_truth, "TeV");
    this->MakeDistribution1D(h_deltaY_top_truth, "");
    this->MakeDistribution2D(h2_mtt_truth_ETmiss_truth, "m_{t\\bar{t}}", "TeV", "E^{\\mathrm{miss,truth}}_{\\mathrm{T}}", "GeV");
    
    this->MakeDistribution1D(h_dR_lb_truth, "");
    // this->MakeDistribution1D(h_dR_l1b1_truth, "");
    // this->MakeDistribution1D(h_dR_l2b2_truth, "");
    this->MakeDistribution1D(h_dR_t1l1_truth, "");
    this->MakeDistribution1D(h_dR_t2l2_truth, "");
    this->MakeDistribution1D(h_dR_t1b1_truth, "");
    this->MakeDistribution1D(h_dR_t2b2_truth, "");
    this->MakeDistribution1D(h_dR_t1t2_truth, "");


    // transverse
    this->MakeDistribution1D(h_HT_all, "TeV");
    this->MakeDistribution1D(h_HT, "TeV");
    this->MakeDistribution1D(h_HTmet, "TeV");
    this->MakeDistribution1D(h_HTjMET, "TeV");
    this->MakeDistribution1D(h_HT_truth, "TeV");
    this->MakeDistribution1D(h_HTmet_truth, "TeV");
    this->MakeDistribution1D(h_KT_all, "TeV");
    this->MakeDistribution1D(h_KT, "TeV");
    this->MakeDistribution1D(h_KT_truth, "TeV");
    this->MakeDistribution1D(h_mass_vis_all, "TeV");
    this->MakeDistribution1D(h_mass_vis, "TeV");
    this->MakeDistribution1D(h_mass_bbll, "TeV");

    this->MakeDistribution1D(h_cosTheta_l, "");
    this->MakeDistribution1D(h_cosThetaStar_l, "");
    
    this->MakeDistribution1D(h_cosTheta1_truth, "");
    this->MakeDistribution1D(h_cosTheta2_truth, "");
    // this->MakeDistribution1D(h_cos1cos2_truth, "");

    this->MakeDistribution1D(h_HT_lF, "TeV");
    this->MakeDistribution1D(h_HT_lB, "TeV");
    this->Asymmetry("AlFB_HT", "A^{\\ell}_{FB^{*}}", "m_{t\\bar{t}}", h_HT_lF, h_HT_lB);

    this->MakeDistribution1D(h_KT_lF, "TeV");
    this->MakeDistribution1D(h_KT_lB, "TeV");
    this->Asymmetry("AlFB_KT", "A^{\\ell}_{FB^{*}}", "m_{t\\bar{t}}", h_KT_lF, h_KT_lB);

    this->MakeDistribution1D(h_HT_lCF, "TeV");
    this->MakeDistribution1D(h_HT_lCB, "TeV");
    this->Asymmetry("AlC_HT", "A^{\\ell}_{C}\\;(H_{\\mathrm{T}})", "H_{\\mathrm{T}}", h_HT_lCF, h_HT_lCB);

    this->MakeDistribution1D(h_KT_lCF, "TeV");
    this->MakeDistribution1D(h_KT_lCB, "TeV");
    this->Asymmetry("AlC_KT", "A^{\\ell}_{C}\\;(K_{\\mathrm{T}})", "K_{\\mathrm{T}}", h_KT_lCF, h_KT_lCB);

    this->MakeDistribution1D(h_HT_DphiF, "TeV");
    this->MakeDistribution1D(h_HT_DphiB, "TeV");
    this->Asymmetry("ADphi_HT", "A_{\\Delta\\phi}", "H_{\\mathrm{T}}", h_HT_DphiF, h_HT_DphiB);

    this->MakeDistribution1D(h_KT_DphiF, "TeV");
    this->MakeDistribution1D(h_KT_DphiB, "TeV");
    this->Asymmetry("ADphi_KT", "A_{\\Delta\\phi}", "K_{\\mathrm{T}}", h_KT_DphiF, h_KT_DphiB);
    
    this->MakeDistribution2D(h2_HT_deltaPhi, "H_{\\mathrm{T}}", "GeV", "\\Delta\\phi_{\\ell^{+}\\ell^{-}}", "");
    this->MakeDistribution2D(h2_HTmet_deltaPhi, "H_{\\mathrm{T}} + E^{\\mathrm{miss}}_{\\mathrm{T}}", "GeV", "\\Delta\\phi_{\\ell^{+}\\ell^{-}}", "");
    this->MakeDistribution2D(h2_mvis_deltaPhi, "m_{\\mathrm{vis}}", "GeV", "\\Delta\\phi_{\\ell^{+}\\ell^{-}}", "");
    this->MakeDistribution2D(h2_KT_deltaPhi, "K_{\\mathrm{T}}", "GeV", "\\Delta\\phi_{\\ell^{+}\\ell^{-}}", "");

    if (m_reconstruction == "KIN" or m_reconstruction == "NuW") {
        
        // energy
        this->MakeDistribution1D(h_E_top, "GeV");
        this->MakeDistribution1D(h_E_tbar, "GeV");
        
        this->MakeDistribution1D(h_pT_top, "GeV");
        this->MakeDistribution1D(h_pT_tbar, "GeV");
        this->MakeDistribution1D(h_pT_ttbar, "GeV");
        
        this->MakeDistribution1D(h_eta_top, "");
        this->MakeDistribution1D(h_eta_tbar, "");
        this->MakeDistribution1D(h_eta_ttbar, "");

        this->MakeDistribution1D(h_phi_top, "");
        this->MakeDistribution1D(h_phi_tbar, "");
        this->MakeDistribution1D(h_phi_ttbar, "");
        this->MakeDistribution1D(h_cosPhi, "");

        // invariant mass
        this->MakeDistribution1D(h_mass_W1, "GeV");
        this->MakeDistribution1D(h_mass_W2, "GeV");
        this->MakeDistribution1D(h_mass_top, "GeV");
        this->MakeDistribution1D(h_mass_tbar, "GeV");
        this->MakeDistribution1D(h_mass_ttbar, "TeV");
        this->MakeDistribution1D(h_dR_top, "");
        this->MakeDistribution1D(h_dR_tbar, "");
        this->MakeDistribution1D(h_dR_ttbar, "");
        this->MakeDistribution1D(h_dR_l, "");
        
        this->MakeDistribution1D(h_perf_mass_top, "");
        this->MakeDistribution1D(h_perf_pT_top, "");
        this->MakeDistribution1D(h_perf_eta_top, "");
        this->MakeDistribution1D(h_perf_phi_top, "");
        this->MakeDistribution1D(h_perf_mass_tbar, "");
        this->MakeDistribution1D(h_perf_pT_tbar, "");
        this->MakeDistribution1D(h_perf_eta_tbar, "");
        this->MakeDistribution1D(h_perf_phi_tbar, "");
        this->MakeDistribution1D(h_perf_mass_ttbar, "");
        this->MakeDistribution1D(h_perf_pT_ttbar, "");
        this->MakeDistribution1D(h_perf_eta_ttbar, "");
        this->MakeDistribution1D(h_perf_phi_ttbar, "");
        this->MakeDistribution1D(h_res_pT_bjets, "");
        this->MakeDistribution2D(h2_resPtBjets_pT, "p^{b}_{\\mathrm{T}}\\, ", "GeV",  "(p^{b,\\mathrm{truth}}_{\\mathrm{T}} - p^{b,\\mathrm{reco}}_{\\mathrm{T}})/p^{b,\\mathrm{truth}}_{\\mathrm{T}}", "");
        this->MakeDistribution1D(h_res_pT_tops, "");
        this->MakeDistribution1D(h_res_y_tops, "");
        this->MakeDistribution1D(h_res_phi_tops, "");
        this->MakeDistribution1D(h_res_eta_tops, "");
        this->MakeDistribution1D(h_res_pT_ttbar, "");
        this->MakeDistribution1D(h_res_y_ttbar, "");
        this->MakeDistribution1D(h_res_mass_ttbar, "");
        this->MakeDistribution1D(h_res_eta_ttbar, "");
        this->MakeDistribution1D(h_reco_quality, "");
        this->MakeDistribution2D(h2_perf_mass_ttbar, "m^{\\mathrm{truth}}_{t\\bar{t}}\\ ", "TeV", "\\left(m^{\\mathrm{truth}}_{t\\bar{t}} - m^{\\mathrm{reco}}_{t\\bar{t}}\\right) / m^{\\mathrm{truth}}_{t\\bar{t}}\\ ", "");
        this->MakeDistribution2D(h2_perf_mass_ttbar_pTtop, "p^{t,\\mathrm{truth}}_{T}\\ ", "GeV", "\\left(m^{\\mathrm{truth}}_{t\\bar{t}} - m^{\\mathrm{reco}}_{t\\bar{t}}\\right) / m^{\\mathrm{truth}}_{t\\bar{t}}\\ ", "");
        this->MakeDistribution2D(h2_perf_mass_ttbar_pTtbar, "p^{\\bar{t},\\mathrm{truth}}_{T}\\ ", "GeV", "\\left(m^{\\mathrm{truth}}_{t\\bar{t}} - m^{\\mathrm{reco}}_{t\\bar{t}}\\right) / m^{\\mathrm{truth}}_{t\\bar{t}}\\ ", "");
        
        this->MakeDistribution2D(h2_dR_lb_mtt_truth, "m^{\\mathrm{truth}}_{t\\bar{t}}\\, ", "TeV", "\\Delta R\\left(b_{\\mathrm{truth}},\\ell_{\\mathrm{truth}}\\right)", "");
        this->AverageEachXbin(h2_dR_lb_mtt_truth);
        // this->MakeDistribution2D(h2_dR_l1b1_mtt_truth, "m^{\\mathrm{truth}}_{t\\bar{t}}\\ ", "TeV", "\\Delta R\\left(b_{\\mathrm{truth}},\\ell^{+}_{\\mathrm{truth}}\\right)", "");
        // this->MakeDistribution2D(h2_dR_l2b2_mtt_truth, "m^{\\mathrm{truth}}_{t\\bar{t}}\\ ", "TeV", "\\Delta R\\left(\\bar{b}_{\\mathrm{truth}},\\ell^{-}_{\\mathrm{truth}}\\right)", "");
        this->MakeDistribution2D(h2_dR_t1l1_mtt_truth, "m^{\\mathrm{truth}}_{t\\bar{t}}\\, ", "TeV", "\\Delta R\\left(t_{\\mathrm{truth}},\\ell^{+}_{\\mathrm{truth}}\\right)", "");
        this->MakeDistribution2D(h2_dR_t2l2_mtt_truth, "m^{\\mathrm{truth}}_{t\\bar{t}}\\, ", "TeV", "\\Delta R\\left(\\bar{t}_{\\mathrm{truth}},\\ell^{-}_{\\mathrm{truth}}\\right)", "");
        this->MakeDistribution2D(h2_dR_t1b1_mtt_truth, "m^{\\mathrm{truth}}_{t\\bar{t}}\\, ", "TeV", "\\Delta R\\left(t_{\\mathrm{truth}}, b_{\\mathrm{truth}}\\right)", "");
        this->MakeDistribution2D(h2_dR_t2b2_mtt_truth, "m^{\\mathrm{truth}}_{t\\bar{t}}\\, ", "TeV", "\\Delta R\\left(\\bar{t}_{\\mathrm{truth}}, \\bar{b}_{\\mathrm{truth}}\\right)", "");
        this->MakeDistribution2D(h2_dR_l1b1_pTl_truth, "p^{\\ell^{+},\\mathrm{truth}}_{\\mathrm{T}}\\, ", "TeV", "\\Delta R\\left(\\ell^{+}_{\\mathrm{truth}}, b_{\\mathrm{truth}}\\right)", "");
        this->MakeDistribution2D(h2_dR_l2b2_pTl_truth, "p^{\\ell^{-},\\mathrm{truth}}_{\\mathrm{T}}\\, ", "TeV", "\\Delta R\\left(\\ell^{-}_{\\mathrm{truth}}, \\bar{b}_{\\mathrm{truth}}\\right)", "");

        // rapidity
        this->MakeDistribution1D(h_y_top, "");
        this->MakeDistribution1D(h_y_tbar, "");
        this->MakeDistribution1D(h_y_ttbar, "");
        this->MakeDistribution1D(h_deltaY_top, "");

        // polar angle
        this->MakeDistribution1D(h_cosTheta, "");
        this->MakeDistribution1D(h_cosThetaStar, "");
        this->MakeDistribution1D(h_cosTheta_ttbar, "");
        this->MakeDistribution1D(h_cosTheta1, "");
        this->MakeDistribution1D(h_cosTheta2, "");
        this->MakeDistribution1D(h_cos1cos2, "");

        // asymmetries
        this->MakeDistributionAL(h2_mtt_cosTheta1, "AL1", "A^{\\ell^{+}}_{L}");
        this->MakeDistributionAL(h2_mtt_cosTheta2, "AL2", "A^{\\ell^{-}}_{L}");

        // charge asymmetry
        this->MakeDistribution1D(h_mtt_tF, "TeV");
        this->MakeDistribution1D(h_mtt_tB, "TeV");
        this->Asymmetry("AtFB", "A^{t}_{FB^{*}}\\", "m_{t\\bar{t}}", h_mtt_tF, h_mtt_tB);

        this->MakeDistribution1D(h_mtt_lF, "TeV");
        this->MakeDistribution1D(h_mtt_lB, "TeV");
        this->Asymmetry("AlFB", "A^{\\ell}_{FB^{*}}", "m_{t\\bar{t}}", h_mtt_lF, h_mtt_lB);

        this->MakeDistribution1D(h_mtt_tCF, "TeV");
        this->MakeDistribution1D(h_mtt_tCB, "TeV");
        this->Asymmetry("AtC", "A^{t}_{C}\\", "m_{t\\bar{t}}", h_mtt_tCF, h_mtt_tCB);

        this->MakeDistribution1D(h_mtt_lCF, "TeV");
        this->MakeDistribution1D(h_mtt_lCB, "TeV");
        this->Asymmetry("AlC", "A^{\\ell}_{C}", "m_{t\\bar{t}}", h_mtt_lCF, h_mtt_lCB);

        // this->MakeDistribution1D(h_mtt_tlF, "TeV");
        // this->MakeDistribution1D(h_mtt_tlB, "TeV");
        // this->Asymmetry("AtlFB", "A^{tl}_{FB^{*}}\\;", "m_{t\\bar{t}}", h_mtt_tlF, h_mtt_tlB);

        this->MakeDistribution1D(h_mtt_c1F, "TeV");
        this->MakeDistribution1D(h_mtt_c1B, "TeV");
        this->Asymmetry("Ac1", "A_{c_{1}}", "m_{t\\bar{t}}", h_mtt_c1F, h_mtt_c1B);

        this->MakeDistribution1D(h_mtt_c2F, "TeV");
        this->MakeDistribution1D(h_mtt_c2B, "TeV");
        this->Asymmetry("Ac2", "A_{c_{2}}", "m_{t\\bar{t}}", h_mtt_c2F, h_mtt_c2B);

        this->MakeDistribution1D(h_mtt_c1c2F, "TeV");
        this->MakeDistribution1D(h_mtt_c1c2B, "TeV");
        this->Asymmetry("Ac1c2", "A_{c_{1}c_{2}}", "m_{t\\bar{t}}", h_mtt_c1c2F, h_mtt_c1c2B);

        this->MakeDistribution1D(h_mtt_cPhiF, "TeV");
        this->MakeDistribution1D(h_mtt_cPhiB, "TeV");
        this->Asymmetry("AcPhi", "A_{\\cos\\varphi}", "m_{t\\bar{t}}", h_mtt_cPhiF, h_mtt_cPhiB);

        this->MakeDistribution1D(h_mtt_DphiF, "TeV");
        this->MakeDistribution1D(h_mtt_DphiB, "TeV");
        this->Asymmetry("ADPhi", "A_{\\Delta\\phi}", "m_{t\\bar{t}}", h_mtt_DphiF, h_mtt_DphiB);

        // 2D distributions
        this->MakeDistribution2D(h2_mtt_deltaPhi, "m_{t\\bar{t}}\\ ", "GeV", "\\Delta\\phi_{\\ell^{+}\\ell^{-}}", "");
        this->MakeDistribution2D(h2_mtt_cosThetaStar, "m_{t\\bar{t}}\\ ", "GeV", "\\cos\\theta^{*}", "");
        this->MakeDistribution2D(h2_mtt_cosTheta1, "m_{t\\bar{t}}\\ ", "GeV", "\\cos\\theta_{\\ell^{+}}", "");
        this->MakeDistribution2D(h2_mtt_cosTheta2, "m_{t\\bar{t}}\\ ", "GeV", "\\cos\\theta_{\\ell^{-}}", "");
        this->MakeDistribution2D(h2_mtt_cos1cos2, "m_{t\\bar{t}}\\ ", "GeV", "\\cos\\theta_{\\ell^{+}}\\cos\\theta_{\\ell^{-}}", "");

        this->MakeDistribution2D(h2_mass_top_TvR, "m_{t}\\ ", "TeV", "m^{truth}_{t}\\ ", "TeV");
        this->MakeDistribution2D(h2_pT_top_TvR, "p^{t}_{\\mathrm{T}}\\ ", "GeV", "p^{truth,t}_{\\mathrm{T}}\\ ", "GeV");
        this->MakeDistribution2D(h2_eta_top_TvR, "\\eta_{t}\\ ", "", "\\eta^{truth}_{t}\\ ", "");
        this->MakeDistribution2D(h2_phi_top_TvR, "\\phi_{t}\\ ", "", "\\phi^{truth}_{t}\\ ", "");
        this->MakeDistribution2D(h2_mass_tbar_TvR, "m_{\\bar{t}}\\ ", "TeV", "m^{truth}_{\\bar{t}}\\ ", "TeV");
        this->MakeDistribution2D(h2_pT_tbar_TvR, "p^{\\bar{t}}_{\\mathrm{T}}\\ ", "GeV", "p^{truth,\\bar{t}}_{\\mathrm{T}}\\ ", "GeV");
        this->MakeDistribution2D(h2_eta_tbar_TvR, "\\eta_{\\bar{t}}\\ ", "", "\\eta^{truth}_{\\bar{t}}\\ ", "");
        this->MakeDistribution2D(h2_phi_tbar_TvR, "\\phi_{\\bar{t}}\\ ", "", "\\phi^{truth}_{\\bar{t}}\\ ", "");
        this->MakeDistribution2D(h2_mass_ttbar_TvR, "m_{t\\bar{t}}\\ ", "TeV", "m^{truth}_{t\\bar{t}}\\ ", "TeV");
        this->MakeDistribution2D(h2_pT_ttbar_TvR, "p^{t\\bar{t}}_{\\mathrm{T}}\\ ", "GeV", "p^{truth,t\\bar{t}}_{\\mathrm{T}}\\ ", "GeV");
        this->MakeDistribution2D(h2_eta_ttbar_TvR, "\\eta_{t\\bar{t}}\\ ", "", "\\eta^{truth}_{t\\bar{t}}\\ ", "");
        this->MakeDistribution2D(h2_phi_ttbar_TvR, "\\phi_{t\\bar{t}}\\ ", "", "\\phi^{truth}_{t\\bar{t}}\\ ", "");
    }
}


void Analysis::MakeDistribution1D(TH1D* h, const string& units, bool normalise)
{
    string ytitle, yunits, xunits;
    if (m_xSec) {
        if (m_luminosity > 0) {
            for (int i = 1; i < h->GetNbinsX() + 1; ++i) {
                h->SetBinError(i, sqrt(h->GetBinContent(i)));
                // cout << "N  = " << "" << h->GetBinContent(i) << "\n";
                // cout << "dN = " << "" << h->GetBinError(i) << "\n";
            }
            ytitle = "Expected events";
        }
        else {
            ytitle = "d\\sigma / d" + (string) h->GetTitle();
            if (units != "") {
                yunits = "\\;[\\mathrm{fb/" + units + "}]";
                xunits = "\\;[\\mathrm{" + units + "}]";
            }
            else {
                yunits = "";
                xunits = "";
            }
        }
    }
    else {
        ytitle = "Generated events";
    }
    if (units != "") xunits = "\\;[\\mathrm{" + units + "}]";
    else xunits = "";
    if (normalise) {
        h->Scale(1.0 / h->Integral());
        ytitle = "Normalised";
    }
    h->GetYaxis()->SetTitle((ytitle + yunits).data());
    h->GetXaxis()->SetTitle((h->GetTitle() + xunits).data());
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleOffset(0.95);
    m_outputFile->cd();
    m_outputFile->cd("/");
    h->Write();
}


void Analysis::WriteEfficiency(TH1D* h, const string& xtitle, const string& ytitle)
{
    const string yTitle = "\\epsilon(" + ytitle + ")";
    h->GetYaxis()->SetTitle(yTitle.data());
    h->GetXaxis()->SetTitle(xtitle.data());
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleOffset(0.95);
    m_outputFile->cd();
    m_outputFile->cd("/");
    h->Write();
}


void Analysis::MakeDistribution2D(TH2D* h, string xtitle, string xunits, string ytitle, string yunits)
{
    string ztitle, zunits;
    if (m_xSec) {
        if (m_luminosity > 0) {
            for (int i = 1; i < h->GetNbinsX() + 1; ++i) {
                for (int j = 1; j < h->GetNbinsY() + 1; ++j) {
                    h->SetBinError(i, j, sqrt(h->GetBinContent(i, j)));
                }
            }
            ztitle = "Expected events";
        }
        ztitle = "d\\sigma / d";
        if (xunits != "" and yunits != "") {
            zunits = "\\;[\\mathrm{fb/" + xunits + "/" + yunits + "}]";
            xunits = "\\;[\\mathrm{" + xunits + "}]";
            yunits = "\\;[\\mathrm{" + yunits + "}]";
        }
        else if (xunits != "" and yunits == "") {
            zunits = "\\;[\\mathrm{fb/" + xunits + "}]";
            xunits = "\\;[\\mathrm{" + xunits + "}]";
        }
        else if (xunits == "" and yunits != "") {
            zunits = "\\;[\\mathrm{fb/" + yunits + "}]";
            yunits = "\\;[\\mathrm{" + yunits + "}]";
        }
        else {
            zunits = "";
        }
        h->GetZaxis()->SetTitle((ztitle + h->GetTitle() + zunits).data());
    }
    else {
        if (m_luminosity > 0) {
            ztitle = "Expected events";
        }
        else {
            ztitle = "Generated events";
        }
        if (xunits != "") xunits = "\\;[\\mathrm{" + xunits + "}]";
        if (yunits != "") yunits = "\\;[\\mathrm{" + yunits + "}]";
        h->GetZaxis()->SetTitle((ztitle).data());
    }
    h->GetYaxis()->SetTitle((ytitle + yunits).data());
    h->GetXaxis()->SetTitle((xtitle + xunits).data());
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleOffset(0.95);
    m_outputFile->cd();
    m_outputFile->cd("/");
    h->Write();
}

void Analysis::MakeDistributionAL(TH2D* h, const string& name, const string& title)
{
    TH2D* hist = (TH2D*) h->Clone();
    TF1* func = new TF1("func1", "[0]*x + [1]", -1, 1);
    TObjArray slices;
    this->NormalizeSliceY(hist);
    hist->FitSlicesY(func, 0, -1, 0, "QRN", &slices);
    if (m_debug) for (auto slice : slices) slice->Write();
    TH1D* h_AL = (TH1D*) slices[0]->Clone(name.data());
    slices.Clear();
    h_AL->Scale(2.0 / hist->GetYaxis()->GetBinWidth(1));
    h_AL->SetTitle(title.data());
    h_AL->GetXaxis()->SetTitle("m_{t\\bar{t}}\\;[TeV]");
    h_AL->GetYaxis()->SetTitle(title.data());
    m_outputFile->cd();
    m_outputFile->cd("/");
    h_AL->Write();
}

void Analysis::NormalizeSliceY(TH2D* h)
{
    double integral = 1;
    int k;
    for (int i = 1; i < h->GetNbinsX() + 1; ++i) {
        integral = h->Integral(i, i, 1, h->GetNbinsY());
        for (int j = 1; j < h->GetNbinsY() + 1; ++j) {
            k = h->GetBin(i, j);
            h->SetBinContent(k, h->GetBinContent(k) / integral);
            h->SetBinError(k, h->GetBinError(k) / integral);
        }
    }
}

void Analysis::AverageEachXbin(TH2D* h)
{
    int n = h->GetNbinsX();
    double x[n];
    double y[n];
    for (int i = 0; i < n; ++i) {
        x[i] = h->GetXaxis()->GetBinCenter(i + 1);
        double sum = 0;
        double points = 0;
        for (int j = 0; j < h->GetNbinsY(); ++j) {
            int k = h->GetBin(i + 1, j + 1);
            sum += h->GetBinContent(k) * h->GetYaxis()->GetBinCenter(j + 1);
            points += h->GetBinContent(k);
        }
        double average;
        if (points == 0) average = 0;
        else average = sum / points;
        y[i] = average;
    }
    
    TGraph* g = new TGraph(n, x, y);
    g->Write();
}

void Analysis::GetHardParticles()
{
    // gets particles from the hard process (after initial state radiation)
    // finds tops (after soft gluon radiation) and daughter particles (after soft radiation)
    
    int iTop = -999, iTbar = -999, iWp = -999, iWm = -999;
    for (int i = 0; i < b_Particle->GetEntriesFast(); ++i) {
        GenParticle* particle = (GenParticle*) b_Particle->At(i);
        int pid = particle->PID;
        int m1 = particle->M1;

        // get top (overwrite with radiating top daughter)
        if (pid == 6) {
            // if (particle->M1 == 0 and particle->M2 == 1)
            m_hardTop = particle;
            iTop = i;
        }

        // get tbar (overwrite with radiating tbar daughter)
        if (pid == -6) {
            // if (particle->M1 == 0 and particle->M2 == 1)
            m_hardTbar = particle;
            iTbar = i;
        }

        // get W index (overwrite with radiating daughter)
        if (pid == 24 and (m1 == iTop or m1 == iWp)) iWp = i;

        if (pid == 5 and m1 == iTop) m_hardB = particle;
        
        if ((pid == -11 or pid == -13) and m1 == iWp) m_hardLepP = particle;
        
        if ((pid == 12 or pid == 14) and m1 == iWp) m_hardNu = particle;
        
        // get W index (overwrite with radiating daughter)
        if (pid == -24 and (m1 == iTbar or m1 == iWm)) iWm = i;

        if (pid == -5 and m1 == iTbar) m_hardBbar = particle;

        if ((pid == 11 or pid == 13) and m1 == iWm) m_hardLepM = particle;

        if ((pid == -12 or pid == -14) and m1 == iWm) m_hardNuBar = particle;
    }

    if (m_debug) cout << "Fetched truth particles." << "\n";
}

void Analysis::GetTruthParticles()
{
    m_truthElectrons = new vector<GenParticle*>;
    m_truthMuons = new vector<GenParticle*>;
    m_truthBquarks = new vector<GenParticle*>;
    for (int i = 0; i < b_Particle->GetEntriesFast(); ++i) {
        GenParticle* particle = (GenParticle*) b_Particle->At(i);
        
        // skip if unstable
        if (particle->Status != 1) continue;
        
        int pid = abs(particle->PID);
        int pT = particle->PT;
        int eta = abs(particle->Eta);

        if (pid == 11 and pT > 25.0 and eta < 2.5)
            m_truthElectrons->push_back(particle);

        if (pid == 13 and pT > 25.0 and eta < 2.5)
            m_truthMuons->push_back(particle);

        if (pid == 5 and pT > 25.0 and eta < 2.5)
            m_truthBquarks->push_back(particle);
    }
    if (m_debug) cout << "Fetched truth particles." << "\n";
}

void Analysis::SelectElectrons()
{
    m_electrons = new vector<Electron*>;
    for (int i = 0; i < b_Electron->GetEntries(); ++i) {
        Electron* electron = (Electron*) b_Electron->At(i);
        if (electron->PT < 27.0) continue;
        double eta = abs(electron->Eta);
        if (eta > 2.47) continue;
        if (eta > 1.37 and eta < 1.52) continue;
        m_electrons->push_back(electron);
    }
    
    if (m_debug) cout << "selected electrons\n";
}

void Analysis::SelectMuons()
{
    m_muons = new vector<Muon*>;
    for (int i = 0; i < b_Muon->GetEntries(); ++i) {
        Muon* muon = (Muon*) b_Muon->At(i);
        if (muon->PT < 27.0) continue;
        if (abs(muon->Eta) > 2.5) continue;
        m_muons->push_back(muon);
    }
    
    if (m_debug) cout << "selected muons\n";
}

void Analysis::IsolateElectrons()
{
    if (m_debug) cout << "isolating electrons\n";
    
    // do not include primary track
    double dRmin = 10e-15;
    
    for (int i = 0; i < m_electrons->size();) {
        Electron* electron = (Electron*) m_electrons->at(i);
        double pT = electron->PT;
        double dRmax = min(0.2, 10.0 / pT); 
        double pTcone = 0.0;
        for (int j = 0; j < b_Track->GetEntries(); ++j) {
            Track* track = (Track*) b_Track->At(j);
            double dR = electron->P4().DeltaR(track->P4());
            if (dR < dRmin) continue;
            if (dR < dRmax) pTcone += track->PT;
        }
        if (pTcone / pT > 0.1) m_electrons->erase(m_electrons->begin() + i);
        else ++i;
    }
    if (m_debug) cout << "isolated electrons\n";
}

void Analysis::IsolateMuons()
{
    if (m_debug) cout << "isolating muons\n"; 
    
    // do not include primary track
    double dRmin = 10e-15;
    
    for (int i = 0; i < m_muons->size();) {
        Muon* muon = (Muon*) m_muons->at(i);
        double pT = muon->PT;
        double dRmax = 10.0 / pT;
        double pTcone = 0.0;
        for (int j = 0; j < b_Track->GetEntries(); ++j) {
            Track* track = (Track*) b_Track->At(j);
            double dR = muon->P4().DeltaR(track->P4());
            if (dR < dRmin) continue;
            if (dR < dRmax) pTcone += track->PT;
        }
        if (pTcone / pT > 0.05) m_muons->erase(m_muons->begin() + i);
        else ++i;
    }
    
    if (m_debug) cout << "isolated muons\n"; 
}

void Analysis::SelectJets()
{
    m_jets = new vector<Jet*>;
    for (int i = 0; i < b_Jet->GetEntries(); ++i) {
        Jet *jet = (Jet*) b_Jet->At(i);
        if (jet->PT < 25.0) return;
        if (abs(jet->Eta) > 2.5) return;
        m_jets->push_back(jet);
    }
}

void Analysis::OverlapRemoval()
{    
    // Objects can satisfy both the jet and lepton selection criteria 
    // a procedure called "overlap removal" is applied to
    // associate objects with a unique hypothesis
    this->RemoveJetsCloseToElectrons();
    // this->RemoveJetsCloseToMuons();
    this->RemoveElectronsInsideJets();
    // this->RemoveMuonsInsideJets();
    
    if (m_debug) cout << "removed overlapping objects\n";
}
    
void Analysis::RemoveJetsCloseToElectrons()
{
    // To prevent double-counting of electron energy deposits as jets
    // the closest small-R jet lying dR < 0.2 from a reconstructed electron is discarded.
    for (int i = 0; i < m_electrons->size(); ++i) {
        Electron* electron = (Electron*) m_electrons->at(i);
        double pT = electron->PT;
        double dRmax = 0.2;
        // double dRmax = min(0.2, 10.0 / pT);
        int jkill = -999;
        for (int j = 0; j < m_jets->size(); ++j) {
            Jet *jet = (Jet*) m_jets->at(j);
            double dRmin = DBL_MAX;
            double dR = electron->P4().DeltaR(jet->P4());
            if (dR < dRmax and dR < dRmin) {
                jkill = j;
                dRmin = dR;
            }
        }
        if (jkill != -999) m_jets->erase(m_jets->begin() + jkill);
    }
    if (m_debug) cout << "removed jets close to electrons\n";
}
    
void Analysis::RemoveJetsCloseToMuons()
{
    // if a jet has fewer than three tracks and is dR < 0.4 from a muon
    // the jet is not considered as a top daughter b-quark
    
    double dRmax = 0.4;
    for (int i = 0; i < m_muons->size(); ++i) {
        Muon* muon = (Muon*) m_muons->at(i);
        double pT = muon->PT;
        for (int j = 0; j < m_jets->size();) {
            Jet *jet = (Jet*) m_jets->at(j);
            double dR = muon->P4().DeltaR(jet->P4());
            if (dR < dRmax) {
                int nTracks = 0;
                for (int k = 0; k < jet->Constituents.GetEntriesFast(); ++k) {
                    TObject* object = jet->Constituents.At(k);
                    if (object == 0) continue;
                    if (object->IsA() == Track::Class()) ++nTracks;
                }
                if (nTracks < 3) m_jets->erase(m_jets->begin() + j);
                else ++j;
            }
            else ++j;
        }
    }
    if (m_debug) cout << "removed jets close to muons\n";
}

void Analysis::RemoveElectronsInsideJets()
{
    // electron is removed if it is too close to a small-R jet
    // reduce the impact of non-prompt leptons
    double dRmax = 0.4;
    for (int i = 0; i < m_jets->size(); ++i) {
        Jet *jet = (Jet*) m_jets->at(i);
        for (int j = 0; j < m_electrons->size();) {
            Electron *electron = (Electron*) m_electrons->at(j);
            double dR = jet->P4().DeltaR(electron->P4());
            if (dR < dRmax) {
                m_electrons->erase(m_electrons->begin() + j);
            }
            else ++j;
        }
    }
    if (m_debug) cout << "removed electrons from jets\n";
}

void Analysis::RemoveMuonsInsideJets()
{
    // muon is removed if it is dR < 0.4 from a small-R jet which has at least three tracks
    // tagged unsuitable for being a top daughter
    
    double dRmax = 0.4;
    for (int i = 0; i < m_jets->size(); ++i) {
        Jet *jet = (Jet*) m_jets->at(i);
        for (int j = 0; j < m_muons->size();) {
            Muon *muon = (Muon*) m_muons->at(j);
            double dR = jet->P4().DeltaR(muon->P4());
            if (dR < dRmax) {
                int nTracks = 0;
                for (int k = 0; k < jet->Constituents.GetEntriesFast(); ++k) {
                    TObject* object = jet->Constituents.At(k);
                    if (object == 0) continue;
                    if (object->IsA() == Track::Class()) ++nTracks;
                }
                if (nTracks > 2) m_muons->erase(m_muons->begin() + j);
                else ++j;
            }
            else ++j;
        }
    }
    if (m_debug) cout << "removed muons from jets\n";
}

void Analysis::TruthTagLeptons()
{  
    double dRmax = 0.01;
    
    // delete m_lepton_truth_tags;
    m_electron_truth_tags = new vector<bool>;
    m_muon_truth_tags = new vector<bool>;
    m_lepton_truth_tags = new vector<bool>;
    
    TLorentzVector p_lepP = m_hardLepP->P4();
    TLorentzVector p_lepM = m_hardLepM->P4();
    
    for (int i = 0; i < m_electrons->size(); ++i) {
        Electron* electron = (Electron*) m_electrons->at(i);
        TLorentzVector p_e = electron->P4();
        double dR_lp = p_e.DeltaR(p_lepP);
        double dR_lm = p_e.DeltaR(p_lepM);
        // cout << "dR_lp = " << dR_lp << ", dR_lm =" << dR_lm << "\n";
        if (dR_lp < dRmax or dR_lm < dRmax)  {
            m_electron_truth_tags->push_back(true);
            m_lepton_truth_tags->push_back(true);
        }
        else  {
            m_electron_truth_tags->push_back(false);
            m_electron_truth_tags->push_back(false);
        }
    }
    
    for (int i = 0; i < m_muons->size(); ++i) {
        Muon* muon = (Muon*) m_muons->at(i);
        TLorentzVector p_mu = muon->P4();
        double dR_lp = p_mu.DeltaR(p_lepP);
        double dR_lm = p_mu.DeltaR(p_lepM);
        // cout << "dR_lp = " << dR_lp << ", dR_lm =" << dR_lm << "\n";
        if (dR_lp < dRmax or dR_lm < dRmax)  {
            m_lepton_truth_tags->push_back(true);
            m_muon_truth_tags->push_back(true);
        }
        else {
            m_lepton_truth_tags->push_back(false);
            m_muon_truth_tags->push_back(false);
        }
    }
    
    // for (auto tag : *m_lepton_truth_tags)
    // {
    //     if (!tag) return false;
    // }
    // 
    // return true;
    
    // for (auto tag : *electron_truth_tags) if (!tag) return true;
    // for (bool tag : *muon_truth_tags) if (!tag) return true;
    
    // return false;
}

void Analysis::FillPurities(int bin)
{   
    double midbin = (double) bin - 0.5;
    this->TruthTagLeptons();    
    // bool event_truth_matched = true;
    // for (auto tag : *m_lepton_truth_tags) if (!tag) event_truth_matched = false;
    // h_lepton_purity->Fill(midbin, (int) event_truth_matched);
    int n_tagged = 0;
    int n_untagged = 0;
    for (auto tag : *m_lepton_truth_tags) h_lepton_purity->Fill(midbin, (int) tag);

    this->TruthTagJets();
    // event_truth_matched = true;
    // for (auto tag : *m_jet_truth_tags) if (!tag) event_truth_matched = false;
    // h_jet_purity->Fill(midbin, (int) event_truth_matched);
    for (auto tag : *m_jet_truth_tags) h_jet_purity->Fill(midbin, (int) tag);
}

void Analysis::TruthTagJets()
{  
    const double dRmax = 0.4;
    
    m_jet_truth_tags = new vector<bool>;
    
    TLorentzVector p_b = m_hardB->P4();
    TLorentzVector p_bbar = m_hardBbar->P4();
    
    for (int i = 0; i < m_jets->size(); ++i) {
        Jet* jet = (Jet*) m_jets->at(i);
        TLorentzVector p_j = jet->P4();
        double dR_b = p_j.DeltaR(p_b);
        double dR_bbar = p_j.DeltaR(p_bbar);
        // cout << "dR_lp = " << dR_lp << ", dR_lm =" << dR_lm << "\n";
        if (dR_b < dRmax) {
            h2_resPtBjets_pT->Fill(p_b.Pt(), (p_b.Pt() - jet->PT) /  p_b.Pt());
            if (p_b.Pt() < 200.0) h_res_pT_bjets->Fill((p_b.Pt() - jet->PT) / p_b.Pt());
            m_jet_truth_tags->push_back(true);
        }
        else if (dR_bbar < dRmax) {
            h2_resPtBjets_pT->Fill(p_bbar.Pt(), (p_bbar.Pt() - jet->PT) / p_bbar.Pt());
            if (p_bbar.Pt() < 200.0) h_res_pT_bjets->Fill((p_bbar.Pt() - jet->PT) /  p_bbar.Pt());
            m_jet_truth_tags->push_back(true);
        }
        else m_jet_truth_tags->push_back(false);
    }
}

bool Analysis::ExactlyTwoLeptons()
{
    bool twoLeptons;
    if (m_electrons->size() + m_muons->size() == 2) twoLeptons = true;
    else twoLeptons = false;
    this->UpdateCutflow(c_twoLeptons, twoLeptons);
    return twoLeptons;
}


void Analysis::AssignChannel()
{
    if (m_electrons->size() == 2) m_channel = "ee";
    else if (m_muons->size() == 2) m_channel = "mumu";
    else if (m_electrons->size() == 1 and m_muons->size() == 1) m_channel = "emu";
    else cout << "error: can't assign channel\n";
    if (m_debug) cout << "Channel assigned: " << m_channel << "\n";
}


bool Analysis::OppositeCharge()
{
    double charge1, charge2;
    if (m_channel == "ee") {
        Electron *electron1 = (Electron*) m_electrons->at(0);
        Electron *electron2 = (Electron*) m_electrons->at(1);

        charge1 = electron1->Charge;
        charge2 = electron2->Charge;
    }
    else if (m_channel == "mumu") {
        Muon *muon1 = (Muon*) m_muons->at(0);
        Muon *muon2 = (Muon*) m_muons->at(1);

        charge1 = muon1->Charge;
        charge2 = muon2->Charge;
    }
    else if (m_channel == "emu") {
        Muon *muon = (Muon*) m_muons->at(0);
        Electron *electron = (Electron*) m_electrons->at(0);

        charge1 = muon->Charge;
        charge2 = electron->Charge;
    }
    else cout << "error: invalid channel\n";

    bool oppositeCharge;

    if (charge1 == charge2) oppositeCharge = false;
    else oppositeCharge = true;
    this->UpdateCutflow(c_oppositeCharge, oppositeCharge);

    return oppositeCharge;
}


void Analysis::FillCutsEfficiencies(const vector<double>& eff_values, const int pass) {
    h_eff_cuts_mass_ttbar_truth->Fill(eff_values.at(0), pass);
    h_eff_cuts_pT_top_truth->Fill(eff_values.at(1), pass);
    h_eff_cuts_pT_tbar_truth->Fill(eff_values.at(2), pass);
}


void Analysis::FillRecoEfficiencies(const vector<double>& eff_values, const int pass) {
    h_eff_reco_mass_ttbar_truth->Fill(eff_values.at(0), pass);
    h_eff_reco_pT_top_truth->Fill(eff_values.at(1), pass);
    h_eff_reco_pT_tbar_truth->Fill(eff_values.at(2), pass);
}


pair<TLorentzVector, TLorentzVector> Analysis::GetLeptonMomenta()
{

    pair<TLorentzVector, TLorentzVector> p_l;

    if (m_channel == "ee") {
        Electron *electron1 = m_electrons->at(0);
        Electron *electron2 = m_electrons->at(1);

        double charge = electron1->Charge;

        if (charge > 0) {
            p_l.first.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(electron1->PT, electron1->Eta, electron1->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(electron2->PT, electron2->Eta, electron2->Phi, 0.0);
        }
    }
    else if (m_channel == "mumu") {
        Muon *muon1 = (Muon*) m_muons->at(0);
        Muon *muon2 = (Muon*) m_muons->at(1);

        double charge = muon1->Charge;

        if (charge > 0) {
            p_l.first.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(muon2->PT, muon2->Eta, muon2->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(muon1->PT, muon1->Eta, muon1->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(muon2->PT, muon2->Eta, muon2->Phi, 0.0);
        }
    }
    else if (m_channel == "emu") {
        Electron *electron = (Electron*) m_electrons->at(0);
        Muon *muon = (Muon*) m_muons->at(0);

        double charge = electron->Charge;

        if (charge > 0) {
            p_l.first.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, 0.0);
            p_l.second.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0.0);
        }
        else {
            p_l.second.SetPtEtaPhiM(electron->PT, electron->Eta, electron->Phi, 0.0);
            p_l.first.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0.0);
        }
    }
    else {
        cout << "error: invalid channel\n";
    }

    if (m_debug) p_l.first.Print();
    if (m_debug) p_l.second.Print();

    return p_l;
}

bool Analysis::SufficientMll(const pair<TLorentzVector, TLorentzVector>& p_l)
{
    // suppress hadronic background, e.g. $J/Psi$
    if (m_debug) cout << "cutting on mll...\n";
    bool sufficientMll;
    double mll = (p_l.first + p_l.second).M();
    if (mll > 15.0) sufficientMll = true;
    else sufficientMll = false;
    this->UpdateCutflow (c_sufficientMll, sufficientMll);
    if (m_debug) cout << "cut on mll\n";
    return sufficientMll;
}


bool Analysis::OutsideZmassWindow(const pair<TLorentzVector, TLorentzVector>& p_l)
{
    // suppress DY+jets background
    if (m_debug) cout << "cutting on |mll - mZ| > 10 ...\n";
    bool outsideZmassWindow;
    double mll = (p_l.first + p_l.second).M();
    if ((abs(mll - m_mass_Z)) > 10.0) outsideZmassWindow = true;
    else outsideZmassWindow = false;
    if (m_channel == "emu") outsideZmassWindow = true;
    this->UpdateCutflow (c_outsideZmassWindow, outsideZmassWindow);
    if (m_debug) cout << "cut on |mll - mZ| > 10 ...\n";
    return outsideZmassWindow;
}


bool Analysis::SufficientBtags()
{
    bool sufficientBtags;
    m_bTags = 0;
    for (int i = 0; i < m_jets->size(); ++i) {
        Jet *jet = (Jet*) m_jets->at(i);
        if (jet->BTag > 0) m_bTags++;
    }
    // h_nBjets->Fill(m_bTags, 1.0);
    if (m_bTags >= m_minBtags) sufficientBtags = true;
    else sufficientBtags = false;
    this->UpdateCutflow(c_sufficientBtags, sufficientBtags);
    return sufficientBtags;
}


bool Analysis::SufficientMET(double ET_miss)
{
    // account for the neutrinos
    // further reduce DY background (no cut for e-mu)
    bool sufficientMET = false;
    if (m_channel == "emu") sufficientMET = true;
    else if (m_channel == "ee" or m_channel == "mumu") {
        if (ET_miss > 60.0) sufficientMET = true;
        else sufficientMET = false;
    }
    else cout << "ERROR: invalid channel\n";
    this->UpdateCutflow(c_sufficientMET, sufficientMET);
    return sufficientMET;
}


bool Analysis::SufficientHT(double HT)
{
    // suppress Z/gamma* ( \tau^+ \tau^-) + jets$
    bool sufficientHT = false;
    if (m_channel == "ee" or m_channel == "mumu") {
        sufficientHT = true;
    }
    else if (m_channel == "emu") {
        // ScalarHT *scalarHT = (ScalarHT*) b_ScalarHT->At(0);
        // double HT = scalarHT->HT;
        if (HT > 130.0) sufficientHT = true;
        else sufficientHT = false;
        // sufficientHT = true;
    }
    else cout << "ERROR: invalid channel\n";
    this->UpdateCutflow(c_sufficientHT, sufficientHT);
    return sufficientHT;
}


bool Analysis::SufficientJets()
{
    bool sufficientJets;
    if (m_jets->size() >= 2) sufficientJets = true;
    else sufficientJets = false;
    this->UpdateCutflow(c_sufficientJets, sufficientJets);
    return sufficientJets;
}


void Analysis::PreLoop()
{
    this->SetDataDirectory();
    this->SetupInputFiles();
    this->SetupOutputFile();
    this->InitialiseCutflow();
    this->MakeHistograms();
    cout << "\n";
}

void Analysis::PreLoopSingle()
{
    this->SetDataDirectory();
    this->SetupInputFile();
    this->SetupOutputFile();
    this->InitialiseCutflow();
    this->MakeHistograms();
    cout << "\n";
}


void Analysis::SetDataDirectory()
{
    // sets directory based on hostname

    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    string Hostname(hostname);

    if (Hostname == "Lorkhan")
        m_dataDirectory = "/Users/declan/Data/zprime/";
    else if ((Hostname.find("lxplus") != string::npos) or (Hostname.find("cern") != string::npos))
        m_dataDirectory = "/afs/cern.ch/work/d/demillar/zprime/";
    else if ((Hostname.find("cyan") != string::npos) or (Hostname.find("blue") != string::npos) or (Hostname.find("green") != string::npos))
        m_dataDirectory = "/scratch/dam1g09/zprime/";
    else
        cout << "Hostname " << Hostname << " not recognised.\n";
}


void Analysis::GetBranches()
{
    b_Particle = m_tree->UseBranch("Particle");
    b_EFlowTrack = m_tree->UseBranch("EFlowTrack");
    b_EFlowPhoton = m_tree->UseBranch("EFlowPhoton");
    b_EFlowNeutralHadron = m_tree->UseBranch("EFlowNeutralHadron");
    b_Jet = m_tree->UseBranch("Jet");
    b_Electron = m_tree->UseBranch("Electron");
    b_Muon = m_tree->UseBranch("Muon");
    b_MissingET = m_tree->UseBranch("MissingET");
    b_Track = m_tree->UseBranch("Track");
}




void Analysis::GetGenerationCrossSection(int proc_id)
{
    if (m_debug) cout << "Getting generation cross section ...\n";
    string proc_filename = get<0>(m_processes->at(proc_id));
    if (m_debug) cout << "Process ID = " << proc_id << "\n";
    if (m_debug) cout << "Process File = " << proc_filename << "\n";

    ifstream proc_file;
    proc_file.open(proc_filename);
    if (!proc_file.is_open()) cout << "error: Unable to open " << proc_filename << "\n";
    if (m_debug) cout << "m_processes size = " << m_processes->size() << "\n";
    get<3>(m_processes->at(proc_id)) = get_parameter(&proc_file);
    get<4>(m_processes->at(proc_id)) = get_parameter(&proc_file);
    proc_file.close();
    if (m_debug) cout << "Generation cross section = " << get<3>(m_processes->at(proc_id)) << " +/- " << get<4>(m_processes->at(proc_id)) << " [fb]\n";
}


void Analysis::GetProcessWeight(int proc_id)
{
    int nproc = get<2>(m_processes->at(proc_id));
    if (m_debug) cout << "nproc = " << nproc << "\n";
    int nevents = m_nevents;
    if (nevents != 10000) cout << nevents << " events\n";
    if (m_xSec) {
        get<5>(m_processes->at(proc_id)) = 1.0 * get<3>(m_processes->at(proc_id)) / (nevents * nproc);
        if (m_luminosity > 0) get<5>(m_processes->at(proc_id)) *= m_luminosity;
    }
    else get<5>(m_processes->at(proc_id)) = 1.0;
    if (m_debug) cout << "Event weight = "<< get<5>(m_processes->at(proc_id)) << "\n";
}


void Analysis::Loop()
{    
    cout << "PROGRESS\n";
    int nfiles = m_input->size();
    for (itr_s it = m_input->begin(); it != m_input->end(); ++it) {
        // boost::timer::auto_cpu_timer auto_timer;
        boost::timer::cpu_timer timer;
        int i = it - m_input->begin();
        if (i + 1 < 10) cout << "   ";
        else if (i + 1 < 100) cout << "  ";
        else if (i + 1 < 1000) cout << " ";
        cout << i + 1 << "/" << nfiles << ": ";
        this->EachFile(get<0>(*it));
        std::cout << get<0>(*it) << "\n";
        if (m_nevents_max > 0) m_nevents = m_nevents_max;
        else m_nevents = this->TotalEvents();
        int proc_id = get<1>(m_input->at(i));
        if (m_xSec)  {
            this->GetGenerationCrossSection(proc_id);
            this->GetProcessWeight(proc_id);
        }
        for (Long64_t jevent = 0; jevent < m_nevents; ++jevent) {
            Long64_t ievent = this->IncrementEvent(jevent);
            if (ievent < 0) break;
            if (m_debug) cout << "\n------------------------------\n";
            if (m_debug) cout << "EVENT " << jevent + 1 << "\n";
            this->EveryEvent(get<5>(m_processes->at(proc_id)));
            this->EachEvent(get<5>(m_processes->at(proc_id)));
            if (m_debug) cout << "------------------------------\n";
            ProgressPercentage(jevent, m_nevents - 1, 50);
        }
        this->CleanUp();
        boost::timer::cpu_times elapsed_time(timer.elapsed());
        auto user_time = elapsed_time.user;
        auto system_time = elapsed_time.system;
        auto loop_time = (user_time + system_time) / 1e9;
        m_event_time = loop_time / m_nevents;
        // cout << "elapsed time = " << loop_time << "\n";
        // cout << "event time = " << m_event_time << "\n";
    }
}


void Analysis::EachFile (const string& filename)
{
    if (m_debug) cout << "starting each file ...\n"; 
    this->SetupTreesForNewFile(filename);
    this->GetBranches();
    if (m_debug) cout << "finished each file ...\n"; 
}


Analysis::~Analysis()
{
    delete m_input;
    delete m_outputFile;
}


Long64_t Analysis::TotalEvents()
{
    if (m_tree != 0) return m_tree->GetEntries();
    return -999;
}


Long64_t Analysis::IncrementEvent(Long64_t i)
{
    Long64_t ev(-1);
    if (m_tree != 0) ev = m_tree->ReadEntry(i);
    return ev;
}


void Analysis::SetupTreesForNewFile(const string& filename)
{

    string tree = "Delphes";

    m_chain = new TChain(tree.data(), "");
    string path = filename + "/" + tree;
    m_chain->Add(path.data(), 0);

    m_tree = new ExRootTreeReader(m_chain);
    // Long64_t numberOfEntries = m_tree->GetEntries();
    // cout << "Number of entries = " << numberOfEntries << "\n";
}


void Analysis::CleanUp()
{
    delete m_chain;
    delete m_tree;
}


double Analysis::TotalAsymmetry(TH1D* h_A, TH1D* h_B)
{
    double A = h_A->Integral("width");
    double B = h_B->Integral("width");
    double Atot = (A - B) / (A + B);
    return Atot;
}


void Analysis::InitialiseCutflow()
{

    if (m_debug) cout << "Initialising cutflow\n";
    m_cutflow = vector<int>(m_cuts, -999);
    m_cutNames = vector<string>(
    m_cuts,                            "no name               ");
    m_cutNames[c_events]             = "Events                ";
    m_cutNames[c_twoLeptons]         = "Two leptons           ";
    m_cutNames[c_oppositeCharge]     = "Opposite Charge       ";
    m_cutNames[c_sufficientMll]      = "Sufficient mll        ";
    m_cutNames[c_outsideZmassWindow] = "Outside Z mass window ";
    m_cutNames[c_sufficientMET]      = "Sufficient MET        ";
    m_cutNames[c_sufficientJets]     = "Sufficient jets       ";
    m_cutNames[c_sufficientBtags]    = "Sufficient b-tags     ";
    m_cutNames[c_sufficientHT]       = "Sufficient HT         ";
    m_cutNames[c_validSolution]      = "Top reconstruction    ";

    m_cutTitles = vector<string>(
    m_cuts,                            "no name");
    m_cutTitles[c_events]             = "Events";
    m_cutTitles[c_twoLeptons]         = "Two leptons";
    m_cutTitles[c_oppositeCharge]     = "Opposite Charge";
    m_cutTitles[c_sufficientMll]      = "\\mathrm{Sufficient}\\; m_{\\ell\\ell}";
    m_cutTitles[c_outsideZmassWindow] = "\\mathrm{Outside}\\; m_{Z}\\; window";
    m_cutTitles[c_sufficientMET]      = "\\mathrm{Sufficient}\\; E^{\\mathrm{miss}}_{\\mathrm{T}}";
    m_cutTitles[c_sufficientJets]     = "Sufficient jets";
    m_cutTitles[c_sufficientBtags]    = "\\mathrm{Sufficient}\\; b\\mathrm{-tags}";
    m_cutTitles[c_sufficientHT]       = "\\mathrm{Sufficient}\\; H_{\\mathrm{T}}";
    m_cutTitles[c_validSolution]      = "Top reconstruction";

    h_cutflow = new TH1D("cutflow", "cutflow", m_cuts, 0.0, m_cuts);

    if (m_debug) cout << "Initialised cutflow\n";
}


void Analysis::UpdateCutflow(const int cut, const bool passed)
{
    if (m_cutflow[cut] == -999) m_cutflow[cut] = 0;
    if (passed) m_cutflow[cut] += 1;
}


void Analysis::PrintCutflow()
{
    cout << "\n\nCUTFLOW\n";
    for (int cut = 0; cut < m_cuts; ++cut) {
        if (m_cutflow[cut] == -999) continue;

        h_cutflow->SetBinContent(cut + 1, m_cutflow[cut]);
        h_cutflow->GetXaxis()->SetBinLabel(cut + 1, m_cutTitles[cut].c_str());

        cout << m_cutNames[cut] << " " << m_cutflow[cut] << "\n";
    }
    cout << "\nRECONSTRUCTION PERFORMANCE\n";
    cout << "Reconstruction efficiency [%]   " << double(double(m_cutflow[m_cuts - 1]) / double(m_cutflow[m_cuts - 2]) * 100.0) << "\n";
    cout << "Reconstruction quality [%]      " << 100.0 * h_reco_quality->GetBinContent(1) << "\n";
    cout << "Average CPU time per event [s]  " << m_event_time << "\n";
    cout << "Resolution (RMS) - top and tbar\n";
    cout << "pT [GeV]                        " << h_res_pT_tops->GetRMS() << "\n";
    cout << "y                               " << h_res_y_tops->GetRMS() << "\n";
    cout << "phi                             " << h_res_phi_tops->GetRMS() << " \n";
    cout << "eta                             " << h_res_eta_tops->GetRMS() << "\n";
    cout << "Resolution (RMS) - ttbar system\n";
    cout << "pT [GeV]                        " << h_res_pT_ttbar->GetRMS() << "\n";
    cout << "y                               " << h_res_y_ttbar->GetRMS() << "\n";
    cout << "mass [GeV]                      " << h_res_mass_ttbar->GetRMS() << "\n";
    cout << "eta                             " << h_res_eta_ttbar->GetRMS() << "\n";
    cout << "Resolution (RMS) - b jets\n";
    cout << "pT [GeV]                        " << h_res_pT_bjets->GetRMS() << "\n";
    // cout << "y                               " << h_res_y_ttbar->GetRMS() << "\n";
    // cout << "mass [GeV]                      " << h_res_mass_ttbar->GetRMS() << "\n";
    // cout << "eta                             " << h_res_eta_ttbar->GetRMS() << "\n";
    // cout << "pT [GeV]                        " << h2_resPtBjets_pT->GetRMS(2) << "\n";
    cout << "Resolution (RMS) - b jets in PT ranges \n";
    for (int i = 0; i < h2_resPtBjets_pT->GetNbinsX(); ++i) {
        h2_resPtBjets_pT->GetXaxis()->SetRange(i, i + 1);
        cout << i + 1 << "  pT [GeV]                        " << h2_resPtBjets_pT->GetRMS(2) << "\n";
    }
    if (m_debug) cout << "Writing cut flow ...\n";
    h_cutflow->Write();
    if (m_debug) cout << "Written cut flow\n";
}

