// TRANSVERSE VARIABLES

    TVector2 ETmiss = -1*pTcol[0] - pTcol[1] - pTcol[2] - pTcol[4];
    MET = ETmiss.Mod();

    mll = (pcol[2] + pcol[4]).M();

    // calculate invariant mass of visible decay products
    Mbbll = (pcol[0] + pcol[1] + pcol[2] + pcol[4]).M();

    // calculate total scalar sum of transverse energy
    HT = ETcol[0] + ETcol[1] + ETcol[2] + ETcol[4] + MET;

    // ET of visible decay products
    TVector2 pTbbll = pTcol[0] + pTcol[1] + pTcol[2] + pTcol[4];
    ETbbll = std::sqrt(Mbbll*Mbbll + pTbbll.Mod2());

    // scalar sum of visible decay products and MET
    KTbbll = ETbbll + MET;

    ET5 = std::sqrt(mass[2]*mass[2] + pTcol[2].Mod2());
    ET7 = std::sqrt(mass[4]*mass[4] + pTcol[4].Mod2());

    MTll = std::sqrt((ET5 + ET7)*(ET5 + ET7) - (pTcol[2] + pTcol[4]).Mod2());
    MCTll = std::sqrt((ET5 + ET7)*(ET5 + ET7) - (pTcol[2] - pTcol[4]).Mod2());

    m35 = (pcol[0] + pcol[2]).M();
    m47 = (pcol[1] + pcol[3]).M();
    TVector2 pT35 = pTcol[0] + pTcol[2];
    TVector2 pT47 = pTcol[1] + pTcol[4];
    ET35 = std::sqrt(m35*m35 - (pTcol[0] + pTcol[2]).Mod2());
    ET47 = std::sqrt(m47*m47 - (pTcol[1] + pTcol[4]).Mod2());

    MTblbl = std::sqrt((ET35 + ET47)*(ET35 + ET47) - (pTcol[0] + pTcol[2] + pTcol[1] + pTcol[4]).Mod2());
    MCTblbl = std::sqrt((ET35 + ET47)*(ET35 + ET47) - (pTcol[0] + pTcol[2] - pTcol[1] - pTcol[4]).Mod2());

    magnitiude of l+ in collider frame
    double plepmag = std::sqrt(pcol[2].Px()*pcol[2].Px() + pcol[2].Py()*pcol[2].Py() + pcol[2].Pz()*pcol[2].Pz());

    // Create basis vectors for the GRRS based on the top direction
    TVector3 x2, y2, z2;
    TLorentzVector plepGRRS = pcol[2];
    z2[1] = 0;
    z2[2] = 0;
    z2[3] = 1;
    double p5xp, p5yp, p5zp;

    // top direction is x'
    for (int i = 1; i < 3; i++) {
      x2[i] = ptcol[i]/sqrt(ptcol[1]*ptcol[1] + ptcol[2]*ptcol[2]);
    }
    x2[3] = 0;

    // y' obtained using cross product
    y2[1] = z2[2]*x2[3] - z2[3]*x2[2];
    y2[2] = z2[3]*x2[1] - z2[1]*x2[3];
    y2[3] = z2[1]*x2[2] - z2[2]*x2[1];

    double plepGRRSx, plepGRRSy, plepGRRSz;
    for (int i = 1; i < 4; i++) {
      p5xp = p5xp + pcol[3][i]*x2[i];
      p5yp = p5yp + pcol[3][i]*y2[i];
      p5zp = p5zp + pcol[3][i]*z2[i];
    }

    double plepmagGRRS = std::sqrt(p5xp*p5xp + p5yp*p5yp + p5zp*p5zp);

    if (std::abs(plepmag - plepmagGRRS) >= 1e-11) printf("Error: coordinate transform mismatch.\n");

    phi_l = atan2(p5yp,p5xp)

    if (phi_l < 0)phi_l = phi_l + 2*pi
    call rootadddouble(phi_l, "fl")

    cosfl = cos(phi_l)
    call rootadddouble(cosfl, "cosphil")
            end if

/ negative velocity of t/tb in collider frame
  TVector3 vtcol = -1*ptcol.BoostVector();
  TVector3 vtbcol = -1*ptbcol.BoostVector();

  TLorentzVector plepptop;
  TLorentzVector plepmtop;



  if (m_channel == "bbllnn") {
    TLorentzVector plepptop = pcol[2];
    TLorentzVector plepmtop = pcol[4];
    plepptop.Boost(vtcol);
    plepmtop.Boost(vtbcol);    
  }