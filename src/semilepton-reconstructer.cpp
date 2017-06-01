#include "semilepton-reconstructer.hpp"


SemileptonReconstructer::SemileptonReconstructer(int bTags, double WMass, double TopMass):
    m_bTags(bTags),
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
    m_b     = TLorentzVector();
    m_bbar  = TLorentzVector();
    m_q     = TLorentzVector();
    m_qbar  = TLorentzVector();
}


bool SemileptonReconstructer::Reconstruct(const TLorentzVector& p_l, double charge, const std::vector<TLorentzVector>& p_b, const std::vector<TLorentzVector>& p_q, const TLorentzVector& p_miss)
{
    bool hasRealRoots;
    double px_l = p_l.Px(), py_l = p_l.Py(), pz_l = p_l.Pz(), E_l;
    double px_miss = p_miss.Px(), py_miss = p_miss.Py();
    E_l = sqrt(px_l * px_l + py_l * py_l + pz_l * pz_l);
    if (std::abs(E_l - p_l.E()) > 1E-5) std::cout << "ERROR: Lepton energy doesn't match.\n";

    double k = m_WMass * m_WMass / 2 + px_l * px_miss + py_l * py_miss; // m_Wmass * m_Wmass / 2 = 3218.42645
    double a = px_l * px_l + py_l * py_l;
    double b = -2 * k * (pz_l);
    double c = (px_miss * px_miss + py_miss * py_miss) * E_l * E_l - k * k;

    double roots[2];
    int nRealRoots = SolveP2(roots, b / a, c / a);
    std::vector<TLorentzVector> p_v;
    if (nRealRoots == 2) {
        hasRealRoots = true;
        // two real solutions; pick best match later
        for (auto& root : roots) {
            double pz = root;
            TLorentzVector p(px_miss, py_miss, pz, sqrt(px_miss * px_miss + py_miss * py_miss + pz * pz));
            p_v.push_back(p);
        }
    }
    else {
        // no real solutions; take the real part
        hasRealRoots = false;
        double pz = roots[0];
        TLorentzVector p(px_miss, py_miss, pz, sqrt(px_miss * px_miss + py_miss * py_miss + pz * pz));
        p_v.push_back(p);
    }

    // FIXME Below only works with two b-tags
    double chi2min = DBL_MAX;
    int B1 = 0, B2 = 0, Q1 = 0, Q2 = 0, V = 0;
    for (unsigned int b1 = 0; b1 < p_b.size(); b1++) {
        for (unsigned int b2 = 0; b2 < p_b.size(); b2++) {
            if (b1 == b2) continue;
            for (unsigned int q1 = 0; q1 < p_q.size(); q1++) {
                for (unsigned int q2 = 0; q2 < p_q.size(); q2++) {
                    if (q1 == q2) continue;
                    for (unsigned int v = 0; v < p_q.size(); v++) {
                        double mbjj = (p_b.at(b1) + p_q.at(q1) + p_q.at(q2)).M();
                        double mblv = (p_b.at(b2) + p_l + p_v.at(v)).M();

                        double dhad = mbjj - m_TopMass;
                        double dlep = mblv - m_TopMass;
                        double chi2 = dhad * dhad + dlep * dlep;
                        if (chi2 < chi2min) {
                            chi2min = chi2;
                            B1 = b1;
                            B2 = b2;
                            Q1 = q1;
                            Q2 = q2;
                            V = v;
                        }
                    }
                }
            }
        }
    }
    m_q = p_q.at(Q1);
    m_qbar = p_q.at(Q2);
    m_nu = p_v.at(V);

    if (charge > 0) {
        m_b = p_b.at(B1);
        m_bbar = p_b.at(B2);
        m_top = m_b + m_nu + p_l;
        m_tbar = m_tbar + m_q + m_qbar;
    }
    else {
        m_b = p_b.at(B2);
        m_bbar = p_b.at(B1);
        m_top = m_tbar + m_q + m_qbar;
        m_tbar = m_b + m_nu + p_l;
    }
    m_ttbar = m_top + m_tbar;

    return hasRealRoots;
}

SemileptonReconstructer::~SemileptonReconstructer(){}
