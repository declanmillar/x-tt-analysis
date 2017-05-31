#include "semilepton-reconstructer.hpp"


SemileptonReconstructer::SemileptonReconstructer(double WMass, double TopMass):
    m_WMass(WMass),
    m_TopMass(TopMass)
{
    this->Reset();
}


void SemileptonReconstructer::Reset(){
  m_top   = TLorentzVector();
  m_tbar  = TLorentzVector();
  m_ttbar = TLorentzVector();
  m_nu    = TLorentzVector();
  m_nubar = TLorentzVector();
  m_b     = TLorentzVector();
  m_bbar  = TLorentzVector();
}


bool SemileptonReconstructer::Reconstruct(const std::pair<TLorentzVector, TLorentzVector>& p_l, const std::vector<TLorentzVector>& p_b, const TLorentzVector& p_miss)
{
    // Returns a std::vector of 4-momenta for all 6 particles in the final state with matching of b-quarks to each top
    // and matching of
    // Takes a std::vector of true final-state particle momenta as the argument and the charge of the final
    // state lepton: if +, t decayed leptonically; if -, t~ decayed leptonically.
    // As going from bbllnn->bblnqq/bbqqln requires only a simple re-weighting for parton truth,
    // it saves on storage space and processing time to store all events as bbllnn. However,
    // when we reconstruct the neutrino, we must account for the fact either the top, or the anti-top
    // may decay hadronically. This means there are two distinguishable final states:
    // Q_l = +1 : pp -> b b~ l+ nu q q'
    // Q_l = -1 : pp -> b b~ q q' l- nu
    // Note that the order here is important, as the order of indices in the std::vector of final state momenta
    // relates to the parent particle t=(0,2,3), t~=(1,4,5) and is fixed at the generator level.
    // If we want the results combining each final state, we must add these together.
    // Note: Experimentally p^{x,y}_nu is equated to the MET, of course.

    // this->UpdateCutflow(c_events, true);

    m_nReco++;

    TLorentzVector p_l, p_nu;
    double a, b, c, k, dh, dl, mblv, mjjb, chi2, chi2min = DBL_MAX;
    unsigned int imin, jmin;

    double px_l = p_l.Px(), py_l = p_l.Py(), pz_l = p_l.Pz(), E_l;
    double px_nu = p_nu.Px(), py_nu = p_nu.Py();
    E_l = sqrt(px_l * px_l + py_l * py_l + pz_l * pz_l);
    if (std::abs(E_l - p_l.E()) > 0.00001) printf("ERROR: Lepton energy doesn't match.\n");

    k = m_WMass * m_WMass / 2 + px_l * px_nu + py_l * py_nu; // m_Wmass * m_Wmass / 2 = 3218.42645
    a = px_l * px_l + py_l * py_l;
    b = -2 * k * (pz_l);
    c = (px_nu * px_nu + py_nu * py_nu) * E_l * E_l - k * k;

    double roots[2];
    int nRealRoots = SolveP2(roots, b / a, c / a);
    p_nu_R.clear();
    if (nRealRoots == 2) {
        // two real solutions; pick best match
        // this->UpdateCutflow(c_realSolutions, true);
        for (auto& root : roots) {
            double pz = root;
            TLorentzVector p(px_nu, py_nu, pz, sqrt(px_nu * px_nu + py_nu * py_nu + pz * pz));
            p_nu_R.push_back(p);
        }
    }
    else {
        // no real solutions; take the real part
        double pz = roots[0];
        TLorentzVector p(px_nu, py_nu, pz, sqrt(px_nu * px_nu + py_nu * py_nu + pz * pz));
        p_nu_R.push_back(p);
    }

    std::vector<TLorentzVector> p_q(4);

    p_q[0] = p[0];
    p_q[1] = p[1];
    if (Q_l == +1) {
        p_q[2] = p[4];
        p_q[3] = p[5];
    }
    else if (Q_l == -1) {
        p_q[2] = p[2];
        p_q[3] = p[3];
    }
    std::vector<std::vector<int> > q_perms;
    if (m_btags == 2) q_perms = { {0, 1, 2, 3}, {1, 0, 2, 3} };
    if (m_btags == 1) q_perms = { {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                  {1, 0, 2, 3}, {2, 0, 1, 3}, {3, 0, 1, 2} };
    if (m_btags == 0) q_perms = { {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                  {1, 0, 2, 3}, {2, 0, 1, 3}, {3, 0, 1, 2},
                                  {2, 0, 1, 3}, {2, 1, 0, 3}, {2, 3, 0, 1},
                                  {3, 0, 1, 2}, {3, 1, 0, 2}, {3, 2, 0, 1} };

    imin = 0;
    jmin = 0;
    for (unsigned int i = 0; i < p_nu_R.size(); i++) {
        for (unsigned int j = 0; j < q_perms.size(); j++) {
            mblv = (p_q[q_perms[j][0]] + p_l + p_nu_R[i]).M();
            mjjb = (p_q[q_perms[j][1]] + p_q[q_perms[j][2]] + p_q[q_perms[j][3]]).M();
            dh = mjjb - m_TopMass;
            dl = mblv - m_TopMass;
            chi2 = dh * dh + dl * dl;
            if (chi2 < chi2min) {
                chi2min = chi2;
                imin = i;
                jmin = j;
            }
        }
    }

    if (Q_l == +1) {
        p_R[0] = p_q[q_perms[jmin][0]];
        p_R[1] = p_q[q_perms[jmin][1]];
        p_R[2] = p_l;
        p_R[3] = p_nu_R[imin];
        p_R[4] = p_q[q_perms[jmin][2]];
        p_R[5] = p_q[q_perms[jmin][3]];
    }
    else if (Q_l == -1) {
        p_R[0] = p_q[q_perms[jmin][1]];
        p_R[1] = p_q[q_perms[jmin][0]];
        p_R[2] = p_q[q_perms[jmin][2]];
        p_R[3] = p_q[q_perms[jmin][3]];
        p_R[4] = p_l;
        p_R[5] = p_nu_R[imin];
    }

    // Assess b-matching performance
    unsigned int b_lep;
    if (Q_l == 1) b_lep = 0;
    if (Q_l == -1) b_lep = 1;
    if (b_lep == jmin) m_nQuarksMatched++;

    // Assess neutrino reconstruction performance.
    double pz_nu_truth = p_nu.Pz();
    double Root0MinusTruth;
    double Root1MinusTruth;
    if (nRealRoots == 2) {
        Root0MinusTruth = std::abs(roots[0] - pz_nu_truth);
        Root1MinusTruth = std::abs(roots[1] - pz_nu_truth);
    }
    else if (nRealRoots == 0) {
        Root0MinusTruth = std::abs(roots[0] - pz_nu_truth);
        Root1MinusTruth = std::abs(roots[0] - pz_nu_truth);
    }
    unsigned int bestRoot;
    if (Root0MinusTruth < Root1MinusTruth) bestRoot = 0;
    else if (Root1MinusTruth < Root0MinusTruth) bestRoot = 1;
    else bestRoot = 0;
    if (imin == bestRoot) m_nNeutrinoMatched++;

    if (m_debug) {
        // Print reconstruction performance.
        printf("True pz_nu = %f\n", p_nu.Pz());
        printf("Possible neutrino solutions:\n");
        if (nRealRoots == 2) {
            printf("                             %f\n", roots[0]);
            printf("                             %f\n", roots[1]);
            printf("Chosen solution:             %f\n", roots[imin]);
        }
        else if (nRealRoots == 0) {
            printf("                             %f + %fi\n", roots[0], roots[1]);
            printf("                             %f - %fi\n", roots[0], roots[1]);
            printf("Chosen solution:             %f\n", roots[imin]);
        }
        if (imin == bestRoot) printf("Neutrino solution: correct. \n");
        else printf("Neutrino solution: incorrect. \n");
        if (b_lep == jmin) printf("b-assignment: correct. \n");
        else printf("b-assignment: incorrect. \n");
        printf("--\n");
    }
    return p_R;
}
