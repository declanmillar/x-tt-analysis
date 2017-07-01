#include "kinematic-reconstructer.hpp"


KinematicReconstructer::KinematicReconstructer( double bMass, double WMass, double TopMass ):
    m_bMass( bMass ),
    m_WMass( WMass ),
    m_TopMass( TopMass )
{
    this->Reset();
}


void KinematicReconstructer::Reset(){
  m_top = TLorentzVector();
  m_tbar = TLorentzVector();
  m_ttbar = TLorentzVector();
  m_nu = TLorentzVector();
  m_nubar = TLorentzVector();
  m_b = TLorentzVector();
  m_bbar = TLorentzVector();
}


bool KinematicReconstructer::Reconstruct(
    const std::pair< TLorentzVector, TLorentzVector >& p_l,
    const std::vector< TLorentzVector >& p_b,
    const std::vector< TLorentzVector >& p_q,
    const TLorentzVector& p_miss )
{
    // Uses the Sonnenschein method (http://arxiv.org/abs/hep-ph/0510100)
    // to algebraically solve tt dilepton equations,
    // then selects solution that minimises mtt.

    if ( m_debug ) std::cout << "beginning dilepton reconstruction ...\n";
    std::vector<TLorentzVector> p_bv;

    int n_b = p_b.size();

    if ( m_debug ) std::cout << "number of b-tagged jets = " << n_b << "\n";

    std::vector<std::vector<std::pair<TLorentzVector, TLorentzVector> > > p_vsol;

    // solve quartic for each possible b-jet combination
    for ( int i = 0; i < n_b; i++ )
        for ( int j = 0; j < n_b; j++ )
            if ( i != j ) p_vsol.push_back( this->GetSolutions( p_l.first, p_l.second, p_b[i], p_b[j], p_miss ) );

    if ( m_debug ) std::cout << "selecting best solution ...\n";
    if ( m_debug ) std::cout << "number b-combinations: " << p_vsol.size() << "\n";

    int nsols = 0;
    for ( int i = 0; i < p_vsol.size(); i++ )
    {
        if ( m_debug ) std::cout << "number of real solutions for b-combo " << i << ": " << p_vsol.at(i).size() << "\n";
        nsols += p_vsol.at(i).size();
    }
    if ( m_debug ) std::cout << "total number of solutions = " << nsols << "\n";
    if ( nsols == 0 ) return false;

    // allow no events with null b-combinations
    // bool no_sols = false;
    // for ( int i = 0; i < p_vsol.size(); i++ )
    // {
    //     if ( p_vsol.at(i).size() == 0 ) no_sols = true;
    // }
    // if ( no_sols ) return p_bv;

    int I = 0, J = 0, K = 0, L = 0, l = 0;
    double mtt_min = DBL_MAX;
    for ( int i = 0; i < n_b; i++ )
    {
        for ( int j = 0; j < n_b; j++ )
        {
            if ( i != j )
            {
                for ( int k = 0; k < p_vsol.at(l).size(); k++ )
                {
                    double mtt = ( p_l.first + p_l.second + p_b[i] + p_b[j] + p_vsol.at(l).at(k).first + p_vsol.at(l).at(k).second ).M();
                    // double MT2 =  asymm_mt2_lester_bisect::get_mT2(
                    //            0, p_l.first, pyA,
                    //            0, pxB, pyB,
                    //            pxMiss, pyMiss,
                    //            0, 0,
                    //            desiredPrecisionOnMt2 );
                    if ( mtt < mtt_min )
                    {
                        mtt_min = mtt;
                        I = i; J = j; K = k; L = l;
                    }
                }
                l++;
            }
        }
    }

    double mW1 = ( p_l.first + p_vsol.at(L).at(K).first ).M();
    double mW2 = ( p_l.second + p_vsol.at(L).at(K).second ).M();
    double mt1 = ( p_l.first + p_b[I] + p_vsol.at(L).at(K).first ).M();
    double mt2 = ( p_l.second + p_b[J] + p_vsol.at(L).at(K).second ).M();
    double mtt = ( p_l.first + p_l.second + p_b[I] + p_b[J] + p_vsol.at(L).at(K).first + p_vsol.at(L).at(K).second ).M();

    if ( false )
    {
        // std::cout << "selected a solution\n";
        if ( std::abs( mW1 - m_WMass ) > 10.0 )
        {
            std::cout << "mW+   = " << mW1 << "\n";
            p_l.first.Print();
            p_vsol.at(L).at(K).first.Print();
            std::cout << "number b-combinations: " << p_vsol.size() << "\n";
            for ( int i = 0; i < p_vsol.size(); i++ )
            {
                std::cout << "number of real solutions for b-combo " << i << ": " << p_vsol.at(i).size() << "\n";
            }
            std::cout << "total number of solutions = " << nsols << "\n";
            std::cout << "I = " << I << " J = " << J << " K = " << K << " L = "<<  L << "\n";
        }
        // if ( std::abs( mW2 - m_WMass ) > 1.0 ) std::cout << "mW-   = " << mW2 << "\n";
        // if ( std::abs( mt1 - m_TopMass ) > 0.1 ) std::cout << "mt    = " << mt1 << "\n";
        // if ( std::abs( mt2 - m_TopMass ) > 0.1 ) std::cout << "mtbar = " << mt2 << "\n";
        // std::cout << "mtt   = " << mtt << "\n";
    }

    m_nu = p_vsol.at(L).at(K).first;
    m_nubar = p_vsol.at(L).at(K).second;
    m_b = p_b.at(I);
    m_bbar = p_b.at(J);
    m_top = m_b + p_l.first + m_nu;
    m_tbar = m_bbar + p_l.second + m_nubar;
    m_ttbar = m_top + m_tbar;

    if ( m_debug ) std::cout << "completed dilepton reconstruction.\n";

    return true;
}

std::vector< std::pair< TLorentzVector, TLorentzVector > > KinematicReconstructer::GetSolutions(
    const TLorentzVector& p_l1, const TLorentzVector& p_l2,
    const TLorentzVector& p_b1, const TLorentzVector& p_b2,
    const TLorentzVector& p_miss )
{
    std::vector< std::pair< TLorentzVector, TLorentzVector > > p_v;

    if ( m_debug ) std::cout << "beginning kinematic reconstruction ...\n";

    // store vector components
    double pl1x = p_l1.Px(), pl1y = p_l1.Py(), pl1z = p_l1.Pz(), El1 = p_l1.E();
    double pl2x = p_l2.Px(), pl2y = p_l2.Py(), pl2z = p_l2.Pz(), El2 = p_l2.E();
    double pb1x = p_b1.Px(), pb1y = p_b1.Py(), pb1z = p_b1.Pz(), Eb1 = p_b1.E();
    double pb2x = p_b2.Px(), pb2y = p_b2.Py(), pb2z = p_b2.Pz(), Eb2 = p_b2.E();
    double Emissx = p_miss.Px(), Emissy = p_miss.Py();

    // store on-shell pole masses
    double mt1 = m_TopMass, mt2 = m_TopMass;
    double mw1 = m_WMass, mw2 = m_WMass;
    double mb1 = m_bMass, mb2 = m_bMass;
    double ml1 = 0.0, ml2 = 0.0;

    double a1 = ( Eb1 + El1 ) * ( mw1 * mw1 - ml1 * ml1 )
        - El1 * ( mt1 * mt1 - mb1 * mb1 - ml1 * ml1 )
        + 2 * Eb1 * El1 * El1
        - 2 * El1 * ( pb1x * pl1x + pb1y * pl1y + pb1z * pl1z );

    double a2 = 2 * ( Eb1 * pl1x - El1 * pb1x );
    double a3 = 2 * ( Eb1 * pl1y - El1 * pb1y );
    double a4 = 2 * ( Eb1 * pl1z - El1 * pb1z );

    double b1 = ( Eb2 + El2 ) * ( mw2 * mw2 - ml2 * ml2 )
        - El2 * ( mt2 * mt2 - mb2 * mb2 - ml2 * ml2 )
        + 2 * Eb2 * El2 * El2
        - 2 * El2 * ( pb2x * pl2x + pb2y * pl2y + pb2z * pl2z );

    double b2 = 2 * ( Eb2 * pl2x - El2 * pb2x );
    double b3 = 2 * ( Eb2 * pl2y - El2 * pb2y );
    double b4 = 2 * ( Eb2 * pl2z - El2 * pb2z );

    double c22 = pow(mw1 * mw1 - ml1 * ml1, 2)
        -4 * ( El1 * El1 - pl1z * pl1z ) * ( a1 / a4 ) * ( a1 / a4 )
        -4 * ( mw1 * mw1 - ml1 * ml1 ) * pl1z * a1 / a4;

    double c21 = 4 * ( mw1 * mw1 - ml1 * ml1 ) * ( pl1x - pl1z * a2 / a4 )
        -8 * ( El1 * El1 - pl1z * pl1z ) * a1 * a2 / ( a4 * a4 )
        -8 * pl1x * pl1z * a1 / a4;

    double c20 = -4 * ( El1 * El1 - pl1x * pl1x )
        -4 * ( El1 * El1 - pl1z * pl1z ) * ( a2 / a4 ) * ( a2 / a4 )
        -8 * pl1x * pl1z * a2 / a4;

    double c11 = 4 * ( mw1 * mw1 - ml1 * ml1 ) * ( pl1y - pl1z * a3 / a4 )
        -8 * ( El1 * El1 - pl1z * pl1z ) * a1 * a3 / ( a4 * a4 )
        -8 * pl1y * pl1z * a1 /  a4;

    double c10 = -8 * ( El1 * El1 - pl1z * pl1z) * a2 * a3 / ( a4 * a4 )
        +8 * pl1x * pl1y
        -8 * pl1x * pl1z * a3 / a4
        -8 * pl1y * pl1z * a2 / a4;

    double c00 = -4 * ( El1 * El1 - pl1y * pl1y )
        -4 * ( El1 * El1 - pl1z * pl1z ) * ( a3 / a4 ) * ( a3 / a4 )
        -8 * pl1y * pl1z * a3 / a4;

    c22 = c22 * a4 * a4;
    c21 = c21 * a4 * a4;
    c20 = c20 * a4 * a4;
    c11 = c11 * a4 * a4;
    c10 = c10 * a4 * a4;
    c00 = c00 * a4 * a4;

    double dd22 = pow( mw2 * mw2 - ml2 * ml2, 2 )
        -4 * ( El2 * El2 - pl2z * pl2z ) * ( b1 / b4 ) * ( b1 / b4 )
        -4 * ( mw2 * mw2 - ml2 * ml2 ) * pl2z * b1 / b4;

    double dd21 = 4 * ( mw2 * mw2 - ml2 * ml2 ) * ( pl2x - pl2z * b2 / b4 )
        -8 * ( El2 * El2 - pl2z * pl2z ) * b1 * b2 / ( b4 * b4 )
        -8 * pl2x * pl2z * b1 / b4;

    double dd20 = -4 * ( El2 * El2 - pl2x * pl2x )
        -4 * ( El2 * El2 - pl2z * pl2z ) * ( b2 / b4 ) * ( b2 / b4 )
        -8 * pl2x * pl2z * b2 / b4;

    double dd11 = 4 * ( mw2 * mw2 - ml2 * ml2 ) * ( pl2y - pl2z * b3 / b4 )
        -8 * ( El2 * El2 -pl2z * pl2z ) * b1 * b3 / ( b4 * b4 )
        -8 * pl2y * pl2z * b1 / b4;

    double dd10 = -8 * ( El2 * El2 - pl2z * pl2z ) * b2 * b3 / ( b4 * b4 )
        +8 * pl2x * pl2y
        -8 * pl2x * pl2z * b3 / b4
        -8 * pl2y * pl2z * b2 / b4;

    double dd00 = -4 * ( El2 * El2 - pl2y * pl2y )
        -4 * ( El2 * El2 - pl2z * pl2z ) * ( b3 / b4 ) * ( b3 / b4 )
        -8 * pl2y * pl2z * b3 / b4;

    dd22 = dd22 * b4 * b4;
    dd21 = dd21 * b4 * b4;
    dd20 = dd20 * b4 * b4;
    dd11 = dd11 * b4 * b4;
    dd10 = dd10 * b4 * b4;
    dd00 = dd00 * b4 * b4;

    double d22 = dd22
        + Emissx * Emissx * dd20
        + Emissy * Emissy * dd00
        + Emissx * Emissy * dd10
        + Emissx * dd21
        + Emissy * dd11;

    double d21 = -dd21
        -2 * Emissx * dd20
        - Emissy * dd10;

    double d20 = dd20;

    double d11 = -dd11
        -2 * Emissy * dd00
        - Emissx * dd10;

    double d10 = dd10;
    double d00 = dd00;

    double h4 = c00 * c00 * d22 * d22
        + c11 * d22 * ( c11 * d00 - c00 * d11 )
        + c00 * c22 * ( d11 * d11 - 2 * d00 * d22 )
        + c22 * d00 * ( c22 * d00 - c11 * d11 );

    double h3 = c00 * d21 * ( 2 * c00 * d22 - c11 * d11 )
        + c00 * d11 * ( 2 * c22 * d10 + c21 * d11 )
        + c22 * d00 * ( 2 * c21 * d00 - c11 * d10 )
        - c00 * d22 * ( c11 * d10 + c10 * d11 )
        -2 * c00 * d00 * ( c22 * d21 + c21 * d22 )
        - d00 * d11 * ( c11 * c21 + c10 * c22 )
        + c11 * d00 * ( c11 * d21 + 2 * c10 * d22 );

    double h2 = c00 * c00 * ( 2 * d22 * d20 + d21 * d21 )
        - c00 * d21 * ( c11 * d10 + c10 * d11 )
        + c11 * d20 * ( c11 * d00 - c00 * d11 )
        + c00 * d10 * ( c22 * d10 - c10 * d22 )
        + c00 * d11 * ( 2 * c21 * d10 + c20 * d11 )
        + ( 2 * c22 * c20 + c21 * c21 ) * d00 * d00
        - 2 * c00 * d00 * ( c22 * d20 + c21 * d21 + c20 * d22 )
        + c10 * d00 * ( 2 * c11 * d21 + c10 * d22 )
        - d00 * d10 * ( c11 * c21 + c10 * c22 )
        - d00 * d11 * ( c11 * c20 + c10 * c21 );

    double h1 = c00 * d21 * ( 2 * c00 * d20 - c10 * d10 )
        - c00 * d20 * ( c11 * d10 + c10 * d11 )
        + c00 * d10 * ( c21 * d10 + 2 * c20 * d11 )
        - 2 * c00 * d00 * ( c21 * d20 + c20 * d21 )
        + c10 * d00 * ( 2 * c11 * d20 + c10 * d21 )
        + c20 * d00 * ( 2 * c21 * d00 - c10 * d11 )
        - d00 * d10 * ( c11 * c20 + c10 * c21 );

    double h0 = c00 * c00 * d20 * d20
        + c10 * d20 * ( c10 * d00 - c00 * d10 )
        + c20 * d10 * ( c00 * d10 - c10 * d00 )
        + c20 * d00 * ( c20 * d00 - 2 * c00 * d20 );

    float a[4] = { static_cast< float >( h1 / h0 ), static_cast< float >( h2 / h0 ), static_cast< float >( h3 / h0 ), static_cast< float >( h4 / h0 ) };

    if ( m_debug ) std::cout << "solving quartic ...\n";
    double x[4];
    int nReal = SolveP4( x, a[0], a[1], a[2], a[3] );
    // std::vector< double > x;
    // x = solveQuartic( 1.0, a[0], a[1], a[2], a[3] );
    // int nReal = x.size();

    if ( m_debug ) std::cout << "found " << nReal << " real solutions.\n";
    if ( nReal == 0 ) return p_v;

    double zero_check;
    for ( int i = 0; i < nReal; i++ )
    {
        zero_check = pow( x[i], 4 ) + a[0] * pow( x[i], 3 ) + a[1] * pow ( x[i], 2 ) + a[2] * x[i] + a[3];
        if ( m_debug and zero_check > 10e-3 ) std::cout << "Warning: solution " << i << " gives: " << zero_check << ", should be zero.\n";
    }

    if ( x[0] != x[0] or x[1] != x[1] or x[2] != x[2] or x[3] != x[3] ) std::cout << "Warning: NaN(s) in quartic solutions.\n";

    bool keepComplex = false;

    std::vector<double> pv1xs;
    int nSolutions;

    if ( keepComplex )
    {
        if ( nReal == 4 )
        {
            // they live in x[0], x[1], x[2], x[3].
            nSolutions = 4;
            for ( int i = 0; i < nSolutions; i++ ) pv1xs.push_back( x[i] );
        }
        else if ( nReal == 2 )
        {
            // x[0], x[1] are the real roots and x[2] ± i * x[3] are the complex
            nSolutions = 3;
            for ( int i = 0; i < nSolutions; i++ ) pv1xs.push_back( x[i] );
        }
        else if ( nReal == 0 )
        {
            // the equation has two pairs of pairs of complex conjugate roots in x[0] ± i * x[1] and x[2] ± i * x[3].
            nSolutions = 2;
            pv1xs.push_back( x[0] );
            pv1xs.push_back( x[2] );
        }
    }
    else
    {
        nSolutions = nReal;
        for ( int i = 0; i < nReal; i++ ) pv1xs.push_back( x[i] );
    }

    // Create pairs of neutrino momenta for each real root
    if ( m_debug ) std::cout << "creating neutrino 4-momenta for each real root ...\n";
    for ( int i = 0; i < nSolutions; i++ )
    {
        double pv1x = x[i];
        double pv2x = Emissx - pv1x;

        double c2 = c22 + c21 * pv1x + c20 * pv1x * pv1x;
        double c1 = c11 + c10 * pv1x;
        double c0 = c00;
        double d2 = d22 + d21 * pv1x + d20 * pv1x * pv1x;
        double d1 = d11 + d10 * pv1x;
        double d0 = d00;

        double pv1y = ( c0 * d2 - c2 * d0 ) / ( c1 * d0 - c0 * d1 );
        double pv2y = Emissy - pv1y;

        double pv1z = -( a1 + a2 * pv1x + a3 * pv1y ) / a4;
        double pv2z = -( b1 + b2 * pv2x + b3 * pv2y ) / b4;

        double Ev1 = sqrt( pv1x * pv1x + pv1y * pv1y + pv1z * pv1z );
        double Ev2 = sqrt( pv2x * pv2x + pv2y * pv2y + pv2z * pv2z );

        TLorentzVector pv1( pv1x, pv1y, pv1z, Ev1 );
        TLorentzVector pv2( pv2x, pv2y, pv2z, Ev2 );

        std::pair<TLorentzVector, TLorentzVector> pv( pv1, pv2 );

        p_v.push_back( pv );
    }

    if ( m_debug ) std::cout << "found " << p_v.size() << " pairs of neutrinos.\n";

    return p_v;
}

KinematicReconstructer::~KinematicReconstructer(){}
