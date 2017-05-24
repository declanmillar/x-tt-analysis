#ifndef _KINEMATIC_RECONSTRUCTION_H_
#define _KINEMATIC_RECONSTRUCTION_H_

#include <fstream>
#include <iostream>

#include "TLorentzVector.h"

class Reconstruction
{
    public:


        int getNumSolutions() const {
            return m_solutions;
        }

        TLorentzVector getv1() const {
            return m_v1;
        }

        TLorentzVector getv2() const {
            return m_v2;
        }

    private:
        int m_solutions;
        TLorentzVector m_v1;
        TLorentzVector m_v2;
};

class KinematicReconstruction
{
    public:
        KinematicReconstruction();
        virtual ~KinematicReconstruction();
        bool apply(const TLorentzVector& l1, const TLorentzVector& b1, const TLorentzVector& l2, const TLorentzVector& b2, const TLorentzVector& ETmiss) const;
        std::string name() const { return "KinematicReconstruction"; }

    private:
        NWSolution solveForNeutrinoEta(const TLorentzVector& lepton, const TLorentzVector& bJet, double topMass, int index) const;
        double neutrino_weight(const TLorentzVector&, const TLorentzVector&, double, double) const;
        double neutrinos[2000][3];
        int etaSize;
        double topMass;
        double m_bmass;
        double m_wmass;
        double sigmax;
        double sigmay;
};

#endif
