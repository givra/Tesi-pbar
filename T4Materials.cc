#include "T4Materials.hh"

T4Materials* T4Materials::materials = NULL;

T4Materials* T4Materials::getInstance(void)
{
    if (materials == NULL) {
        materials = new T4Materials();
    }
    return materials;
}

pair<G4Material*,G4Material*> T4Materials::BuildDYtarget(double packing_factor1, double packing_factor2, int year)
{
    string name;
    double density;
    double A;
    double fractionmass;

    // Cell 1
    name = "dy_ptcell1_" + to_string(year);
    density = ammonia->GetDensity() * packing_factor1 + lhe->GetDensity() * (1-packing_factor1);
    A = lhe->GetDensity()/ammonia->GetDensity() * (1-packing_factor1) / packing_factor1;
    fractionmass = 1/(1+A);

    G4Material *cell1 = new G4Material(name, density, 2);
    cell1->AddMaterial(this->ammonia,fractionmass);
    cell1->AddMaterial(this->lhe,1-fractionmass);

    // Cell 2
    name = "dy_ptcell2_" + to_string(year);
    density = ammonia->GetDensity() * packing_factor2 + lhe->GetDensity() * (1-packing_factor2);
    A = lhe->GetDensity()/ammonia->GetDensity() * (1-packing_factor2) / packing_factor2;
    fractionmass = 1/(1+A);

    G4Material *cell2 = new G4Material(name, density, 2);
    cell2->AddMaterial(this->ammonia,fractionmass);
    cell2->AddMaterial(this->lhe,1-fractionmass);

    return make_pair(cell1,cell2);
}

T4Materials::T4Materials(void)
{
    useOptical = T4SettingsFile::getInstance()->isOpticalPhysicsActivated();

    DefineIndex();

    manager = G4NistManager::Instance();
    elH = manager->FindOrBuildElement("H");
    elC = manager->FindOrBuildElement("C");
    elN = manager->FindOrBuildElement("N");
    elO = manager->FindOrBuildElement("O");
    elF = manager->FindOrBuildElement("F");
    elFe = manager->FindOrBuildElement("Fe");
    elRb = manager->FindOrBuildElement("Rb");
    elMo = manager->FindOrBuildElement("Mo");
    elSb = manager->FindOrBuildElement("Sb");
    elCs = manager->FindOrBuildElement("Cs");
    elNi = manager->FindOrBuildElement("Ni");
    elCr = manager->FindOrBuildElement("Cr");
    elPb = manager->FindOrBuildElement("Pb");
    elCe = manager->FindOrBuildElement("Ce");
    elAs = manager->FindOrBuildElement("As");
    elK = manager->FindOrBuildElement("K");
    elSi = manager->FindOrBuildElement("Si");
    elLi = manager->FindOrBuildElement("Li");
    elCl = manager->FindOrBuildElement("Cl");
    elHe = manager->FindOrBuildElement("He");
    elBa = manager->FindOrBuildElement("Ba");
    elTi = manager->FindOrBuildElement("Ti");

    BuildAir();
    BuildBc408();
    BuildBcf10();
    BuildPMMA();
    BuildPlexiglass();
    BuildGlass();
    BuildAluminium();
    BuildBialkali();
    BuildC4F10();
    BuildN2();

    c2h6 = new G4Material("c2h6", 1.290 * CLHEP::mg / CLHEP::cm3, 2);
    c2h6->AddElement(elC, 0.25);
    c2h6->AddElement(elH, 0.75);

    neon = manager->FindOrBuildMaterial("G4_Ne");
    neon->SetName("neon");

    borosilicate = manager->FindOrBuildMaterial("G4_Pyrex_Glass");
    borosilicate->SetName("Borosilicate");

    plasticFoam = new G4Material("PlasticFoam", 0.15 * CLHEP::g / CLHEP::cm3, 3);
    plasticFoam->AddElement(elC, 2);
    plasticFoam->AddElement(elH, 2);
    plasticFoam->AddElement(elO, 1);

    carbonFibre = new G4Material("CarbonFibre", 1.8 * CLHEP::g / CLHEP::cm3, 2);
    carbonFibre->AddElement(elC, 10);
    carbonFibre->AddElement(elN, 1);

    carbonLH2 = new G4Material("carbonLH2", 2.073 * CLHEP::g / CLHEP::cm3, 2);
    carbonLH2->AddElement(elC, 10);
    carbonLH2->AddElement(elN, 1);

    nomex = new G4Material("nomex", 0.03 * CLHEP::g / CLHEP::cm3, 1);
    nomex->AddElement(elC, 1.0);

    epoxy = new G4Material("epoxy", 1.2 * CLHEP::g / CLHEP::cm3, 3);
    epoxy->AddElement(elC, 21);
    epoxy->AddElement(elH, 25);
    epoxy->AddElement(elO, 5);

    epoxy_09 = new G4Material("epoxy_09", 0.9 * CLHEP::g / CLHEP::cm3, 1);
    epoxy_09->AddMaterial(epoxy, 1.0);

    concrete = manager->FindOrBuildMaterial("G4_CONCRETE");
    concrete->SetName("Concrete");

    lh2 = manager->FindOrBuildMaterial("G4_lH2");
    lh2->SetName("LH2");

    mylar = manager->FindOrBuildMaterial("G4_MYLAR");
    mylar->SetName("Mylar");

    nylon = manager->FindOrBuildMaterial("G4_NYLON-6/6");
    nylon->SetName("Nylon");

    stainlessSteel = new G4Material(
        "StainlessSteel",
        7.9 * CLHEP::g / CLHEP::cm3, 3, kStateSolid, 50.0 * CLHEP::kelvin,
        1.0 * CLHEP::atmosphere
        );
    stainlessSteel->AddElement(elC, 1);
    stainlessSteel->AddElement(elFe, 10);
    stainlessSteel->AddElement(elMo, 1);

    polyethylene = manager->FindOrBuildMaterial("G4_POLYETHYLENE");
    polyethylene->SetName("Polyethylene");

    polypropylene = manager->FindOrBuildMaterial("G4_POLYPROPYLENE");
    polypropylene->SetName("Polypropylene");

    polystyrene = manager->FindOrBuildMaterial("G4_POLYSTYRENE");
    polystyrene->SetName("Polystyrene");

    pvc = manager->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
    pvc->SetName("pvc");

    airexC70 = new G4Material("airexC70", 40 * CLHEP::kg / CLHEP::m3, 1);
    airexC70->AddMaterial(pvc, 1.0);

    iron = manager->FindOrBuildMaterial("G4_Fe");
    iron->SetName("Iron");

    silver = manager->FindOrBuildMaterial("G4_Ag");
    silver->SetName("silver");

    helium = manager->FindOrBuildMaterial("G4_He");
    helium->SetName("Helium");

    G4Isotope* elHe3  = new G4Isotope("He3", 2., 3, 3.016*CLHEP::g/CLHEP::mole);
    G4Isotope* elHe4  = new G4Isotope("He4", 2., 4, 4.003*CLHEP::g/CLHEP::mole);
    G4Element* elDilutedLHe = new G4Element("DilutedLHe", "He", 2);

    double ratioHe34 = 0.10;
    elDilutedLHe->AddIsotope(elHe4,1-ratioHe34);
    elDilutedLHe->AddIsotope(elHe3,ratioHe34);

    double diluted_density = 0.08243*ratioHe34 + 0.1451*(1-ratioHe34);
    lhe = new G4Material("DilutedLHeCryo", diluted_density * CLHEP::g/CLHEP::cm3, 1, kStateLiquid, 0.5965*CLHEP::kelvin, 0.68*CLHEP::bar);
    lhe->AddElement(elDilutedLHe,1);
    //Previous implementation
    //lhe = new G4Material("LHeCryo", 0.1451 * CLHEP::g/CLHEP::cm3, 1, kStateLiquid, 0.5965*CLHEP::kelvin, 0.68*CLHEP::bar);
    //lhe->AddElement(elHe,1);

    ammonia = new G4Material("NH3CryoCrystal", 0.853 * CLHEP::g/CLHEP::cm3, 2, kStateSolid, 0.5965*CLHEP::kelvin, 0.68*CLHEP::bar);
    ammonia->AddElement(elN, 1);
    ammonia->AddElement(elH, 3);

    G4Isotope* elH1   = new G4Isotope("H1" , 1, 1, 1.000*CLHEP::g/CLHEP::mole);
    G4Isotope* elH2   = new G4Isotope("H2" , 1, 2, 2.012*CLHEP::g/CLHEP::mole);
    G4Element* elD    = new G4Element("D","H",2);
    double ratioH = 1.00;
    elD->AddIsotope(elH1,1-ratioH);
    elD->AddIsotope(elH2,  ratioH);
    elD->SetName("Deuterium");

    // Create material using the element
    ld2 = new G4Material("LD2", 0.162 * CLHEP::g/CLHEP::cm3, 1, kStateLiquid, 24.04*CLHEP::kelvin, 1.13*CLHEP::bar);
    ld2->AddElement(elD, 1.0);

    LiDeuterated = new G4Material("LiDCryoCrystal", 0.82 * CLHEP::g/CLHEP::cm3, 2, kStateSolid, 0.50*CLHEP::kelvin, 0.68*CLHEP::bar);
    LiDeuterated->AddElement(elLi, 1);
    LiDeuterated->AddElement(elD , 2);

    packing_factor_1st_2014 = 0.5212;
    packing_factor_2nd_2014 = 0.4558;
    pair<G4Material*,G4Material*> dy_PTtarget_2014 = BuildDYtarget(packing_factor_1st_2014,packing_factor_2nd_2014, 2014);
    nh3_he_1st_2014 = dy_PTtarget_2014.first;
    nh3_he_2nd_2014 = dy_PTtarget_2014.second;

    packing_factor_1st_2015 = 0.5657;
    packing_factor_2nd_2015 = 0.4797;
    pair<G4Material*,G4Material*> dy_PTtarget_2015 = BuildDYtarget(packing_factor_1st_2015, packing_factor_2nd_2015, 2015);
    nh3_he_1st_2015 = dy_PTtarget_2015.first;
    nh3_he_2nd_2015 = dy_PTtarget_2015.second;

    packing_factor_1st_2018 = 0.558;
    packing_factor_2nd_2018 = 0.526;
    pair<G4Material*,G4Material*> dy_PTtarget_2018 = BuildDYtarget(packing_factor_1st_2018, packing_factor_2nd_2018, 2018);
    nh3_he_1st_2018 = dy_PTtarget_2018.first;
    nh3_he_2nd_2018 = dy_PTtarget_2018.second;

    nh3_he_fluka = new G4Material("nh3_he_fluka", 0.425 * CLHEP::g / CLHEP::cm3, 2);
    nh3_he_fluka->AddMaterial(helium, 0.5);
    nh3_he_fluka->AddMaterial(ammonia, 0.5);

    nh3_he       = new G4Material("nh3_he", 0.425 * CLHEP::g / CLHEP::cm3, 2);
    nh3_he      ->AddMaterial(helium, 0.5);
    nh3_he      ->AddMaterial(ammonia, 0.5);

    dli_he_2021  = new G4Material("dli_he_2021", 0.425 * CLHEP::g / CLHEP::cm3, 2);
    dli_he_2021 ->AddMaterial(helium, 0.5);
    dli_he_2021 ->AddMaterial(LiDeuterated, 0.5);

    //Target holder for 2015 (diff from 2014)
    C2F3Cl = new G4Material("C2F3Cl",2.13 * CLHEP::g / CLHEP::cm3,3);
    C2F3Cl->AddElement(elC,2);
    C2F3Cl->AddElement(elF,3);
    C2F3Cl->AddElement(elCl,1);

    //Define the equivalent of "aluminum low dens"
    AlMylar = new G4Material("Al_mylar", 0.729 * CLHEP::g / CLHEP::cm3, 1);
    AlMylar->AddMaterial(aluminium_noOptical, 1.0);

    argon = manager->FindOrBuildMaterial("G4_Ar");
    argon->SetName("Argon");

    co2 = manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    co2->SetName("CO2");

    kapton = new G4Material("Kapton", 1.42 * CLHEP::g / CLHEP::cm3, 4);
    kapton->AddElement(elH, 0.0273);
    kapton->AddElement(elC, 0.7213);
    kapton->AddElement(elN, 0.0765);
    kapton->AddElement(elO, 0.1749);

    Cu = manager->FindOrBuildMaterial("G4_Cu");
    Cu->SetName("Copper");

    Nb = manager->FindOrBuildMaterial("G4_Nb");
    Nb->SetName("Niobium");

    Ti = manager->FindOrBuildMaterial("G4_Ti");
    Ti->SetName("Titanium");

    Ba = manager->FindOrBuildMaterial("G4_Ba");
    Ba->SetName("Barium");

    // For capacitors
    BaTiO3 = new G4Material("BaTiO3",6.02 * CLHEP::g / CLHEP::cm3, 3);
    BaTiO3->AddElement(elBa,1);
    BaTiO3->AddElement(elTi,1);
    BaTiO3->AddElement(elO,3);

    NbTi = new G4Material("NbTi", 5.749 * CLHEP::g / CLHEP::cm3, 2);
    NbTi->AddMaterial(Nb, 0.3);
    NbTi->AddMaterial(Ti, 0.7);

    NbTi_Cu = new G4Material("NbTi_Cu", 7.5937 * CLHEP::g / CLHEP::cm3, 2);
    NbTi_Cu->AddMaterial(Cu, 0.5745);
    NbTi_Cu->AddMaterial(NbTi, 0.4255);

    superconducter = new G4Material("superconducter", 4.95239 * CLHEP::g / CLHEP::cm3, 4);
    superconducter->AddMaterial(Cu, 0.45);
    superconducter->AddMaterial(Nb, 0.045);
    superconducter->AddMaterial(Ti, 0.105);
    superconducter->AddMaterial(lhe, 0.4);

    kevlar = manager->FindOrBuildMaterial("G4_KEVLAR");
    kevlar->SetName("Kevlar");

    Si = manager->FindOrBuildMaterial("G4_Si");
    Si->SetName("Silicon");

    Cu_10 = new G4Material("Cu_10", 0.896 * CLHEP::g / CLHEP::cm3, 1);
    Cu_10->AddMaterial(Cu, 1.0);

    Cu_Kap = new G4Material("Cu_Kap", 1.78 * CLHEP::g / CLHEP::cm3, 2);
    Cu_Kap->AddMaterial(Cu, 0.66);
    Cu_Kap->AddMaterial(kapton, 0.34);

    Be = manager->FindOrBuildMaterial("G4_Be");
    Be->SetName("Beryllium");

    tungsten = manager->FindOrBuildMaterial("G4_W");
    tungsten->SetName("tungsten");

    nickel = manager->FindOrBuildMaterial("G4_Ni");
    nickel->SetName("nickel");

    gold = manager->FindOrBuildMaterial("G4_Au");
    gold->SetName("Gold");

    lead = manager->FindOrBuildMaterial("G4_Pb");
    lead->SetName("Lead");

    ethan = manager->FindOrBuildMaterial("G4_ETHANE");
    ethan->SetName("Ethan");

    graphit = manager->FindOrBuildMaterial("G4_C");
    graphit->SetName("graphit");

    Li2CO3 = new G4Material("Li2C03", 1.36 * CLHEP::g / CLHEP::cm3, 3);
    Li2CO3->AddElement(elLi, 2);
    Li2CO3->AddElement(elC,1);
    Li2CO3->AddElement(elO,3);
    //comment

    cf4 = new G4Material("CF4", 3.72 * CLHEP::mg / CLHEP::cm3, 2);
    cf4->AddElement(elC, 1);
    cf4->AddElement(elF, 4);

    mmSupport = new G4Material("mmSupport", 0.125 * CLHEP::g / CLHEP::cm3, 1);
    mmSupport->AddMaterial(neon, 1);

    g10 = new G4Material("G10", 1.7 * CLHEP::g / CLHEP::cm3, 3);
    g10->AddElement(elC, 5);
    g10->AddElement(elH, 8);
    g10->AddElement(elO, 2);

    g10_10 = new G4Material("G10_10", 0.17 * CLHEP::g / CLHEP::cm3, 1);
    g10_10->AddMaterial(g10, 1.0);

    g11 = new G4Material("G11", 1.9 * CLHEP::g / CLHEP::cm3, 3);
    g11->AddElement(elC, 5);
    g11->AddElement(elH, 8);
    g11->AddElement(elO, 2);

    ArCo2 = new G4Material("ArCo2", 1.78 * CLHEP::mg / CLHEP::cm3, 2);
    ArCo2->AddMaterial(argon, 0.70);
    ArCo2->AddMaterial(co2, 0.30);

    dc4gas = new G4Material("dc4gas", 1.65 * CLHEP::mg / CLHEP::cm3, 3);
    dc4gas->AddMaterial(argon, 0.45);
    dc4gas->AddMaterial(ethan, 0.45);
    dc4gas->AddMaterial(cf4, 0.1);

    dcEBox = new G4Material("dcEBox", 4.323 * CLHEP::g / CLHEP::cm3, 4);
    dcEBox->AddMaterial(aluminium_noOptical, 0.2 / 1.1);
    dcEBox->AddMaterial(stainlessSteel, 0.3 / 1.1);
    dcEBox->AddMaterial(g11, 0.5 / 1.1);
    dcEBox->AddMaterial(Cu, 0.1 / 1.1);

    strawGas = new G4Material("strawGas", 2.08 * CLHEP::mg / CLHEP::cm3, 3);
    strawGas->AddMaterial(argon, 0.74);
    strawGas->AddMaterial(co2, 0.06);
    strawGas->AddMaterial(cf4, 0.20);

    w45Gas = new G4Material("w45Gas", 1.88 * CLHEP::mg / CLHEP::cm3, 3);
    w45Gas->AddMaterial(argon, 0.85);
    w45Gas->AddMaterial(co2, 0.05);
    w45Gas->AddMaterial(cf4, 0.10);

    // wire: diameter = 100 micrometer, pitch = 2 mm
    // area full: 2 mm * 100 micrometer
    // area wire: PI * (50 micrometer)^2
    // ==> 3.927%
    // density = 0.96073*1.88 + 0.03927*1848 = 74.377g/cm3
    w45CathodeGas = new G4Material("w45CathodeGas", 74.377 * CLHEP::mg / CLHEP::cm3, 2);
    w45CathodeGas->AddMaterial(w45Gas, 0.96073);
    w45CathodeGas->AddMaterial(Be, 0.03927);

    mdtGas = new G4Material("mdtGas", 1.716 * CLHEP::mg / CLHEP::cm3, 2);
    mdtGas->AddMaterial(argon, 0.7);
    mdtGas->AddMaterial(co2, 0.3);

    mw2Gas = new G4Material("mw2Gas", 1.41 * CLHEP::mg / CLHEP::cm3, 2);     // density = 1.17 * density_air
    mw2Gas->AddMaterial(argon, 0.75);
    mw2Gas->AddMaterial(ethan, 0.25);

    plastic = manager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    plastic->SetName("plastic");

    rohacell = new G4Material("rohacell", 51.3 * CLHEP::kg / CLHEP::m3, 4);     // C9 H13 N O2  0.0513 g/cc
    rohacell->AddElement(elC, 9);
    rohacell->AddElement(elH, 13);
    rohacell->AddElement(elN, 1);
    rohacell->AddElement(elO, 2);

    rohacell_30 = new G4Material("rohacell_30", 30.0 * CLHEP::kg / CLHEP::m3, 1);
    rohacell_30->AddMaterial(rohacell, 1.0);

    inox = new G4Material("innox", 7.90 * CLHEP::g / CLHEP::cm3, 3);
    inox->AddElement(elC, 0.002737583);
    inox->AddElement(elCr, 0.664841611);
    inox->AddElement(elNi, 0.332420806);

    siliconMaterial = manager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    siliconMaterial->SetName("silicon-material");

    gemMaterial = new G4Material("gemMaterial", 1.7 * CLHEP::g / CLHEP::cm3, 3);
    gemMaterial->AddElement(elC, 0.625);
    gemMaterial->AddElement(elH, 0.04166667);
    gemMaterial->AddElement(elO, 0.33333333);

    mwpcGas = new G4Material("mwpcGas", 2.08 * CLHEP::mg / CLHEP::cm3, 3);
    mwpcGas->AddMaterial(argon, 0.74);
    mwpcGas->AddMaterial(co2, 0.06);
    mwpcGas->AddMaterial(cf4, 0.20);

    pixelGEMMaterial = new G4Material("pgem-material",
                                      4.49 * CLHEP::g / CLHEP::cm3, 2);
    pixelGEMMaterial->AddMaterial(kapton, 0.91);
    pixelGEMMaterial->AddMaterial(Cu, 0.09);

    // lead glass modules
    PbO = new G4Material("PbO", 9.53 * CLHEP::g / CLHEP::cm3, 2);
    PbO->AddElement(elPb, 1);
    PbO->AddElement(elO, 1);

    SiO2 = new G4Material("SiO2", 2.65 * CLHEP::g / CLHEP::cm3, 2);
    SiO2->AddElement(elSi, 1);
    SiO2->AddElement(elO, 2);

    K2O = new G4Material("K2O", 2.32 * CLHEP::g / CLHEP::cm3, 2);
    K2O->AddElement(elK, 2);
    K2O->AddElement(elO, 1);

    As2O3 = new G4Material("As2O3", 3.74 * CLHEP::g / CLHEP::cm3, 2);
    As2O3->AddElement(elAs, 2);
    As2O3->AddElement(elO, 3);

    CeO2 = new G4Material("CeO2", 7.215 * CLHEP::g / CLHEP::cm3, 2);
    CeO2->AddElement(elCe, 1);
    CeO2->AddElement(elO, 2);

    tf1 = new G4Material("tf1", 3.86 * CLHEP::g / CLHEP::cm3, 4);
    tf1->AddMaterial(PbO, 0.512);
    tf1->AddMaterial(SiO2, 0.413);
    tf1->AddMaterial(K2O, 0.07);
    tf1->AddMaterial(As2O3, 0.005);

    tf101 = new G4Material("tf101", 3.86 * CLHEP::g / CLHEP::cm3, 4);
    tf101->AddMaterial(PbO, 0.512);
    tf101->AddMaterial(SiO2, 0.416);
    tf101->AddMaterial(K2O, 0.07);
    tf101->AddMaterial(CeO2, 0.002);

    sf57 = new G4Material("sf57", 5.51 * CLHEP::g / CLHEP::cm3, 3);
    sf57->AddMaterial(PbO, 0.748);
    sf57->AddMaterial(SiO2, 0.239);
    sf57->AddMaterial(K2O, 0.013);

    sf5 = new G4Material("sf5", 4.08 * CLHEP::g / CLHEP::cm3, 2);
    sf5->AddMaterial(PbO, 0.55);
    sf5->AddMaterial(SiO2, 0.45);

    Al2O3 = manager->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    Al2O3->SetName("Al2O3");

    pctfe = manager->FindOrBuildMaterial("G4_POLYTRIFLUOROCHLOROETHYLENE");
    pctfe->SetName("pctfe");

    Klegecell = new G4Material("Klegecell", 0.055 * CLHEP::g / CLHEP::cm3, 3);
    Klegecell->AddElement(elCl, 0.57);
    Klegecell->AddElement(elH, 0.05);
    Klegecell->AddElement(elC, 0.38);

    RICHWindowEff = new G4Material("RICHWindowEff", 0.125 * CLHEP::g / CLHEP::cm3, 2);
    RICHWindowEff->AddMaterial(Klegecell, 0.5);
    RICHWindowEff->AddMaterial(aluminium_noOptical, 0.5);

        paper = new G4Material("mcote_paper", 0.85*CLHEP::g / CLHEP::cm3, 3);
        paper->AddElement(elH , 10);
        paper->AddElement(elO , 5);
        paper->AddElement(elC , 6);

    if (useOptical) {
        BuildSurfaceBc408Air();
        BuildSurfaceAirBc408();
        BuildSurfaceBc408Plexiglass();
        BuildSurfacePlexiglassGlass();
        BuildSurfacePlexiglassAir();
        BuildSurfaceAirPlexiglass();
        BuildSurfaceGlassVacuum();
        BuildSurfaceVacuumBialkali();
        BuildSurfaceAluminium();
        BuildSurfaceMirror();
        BuildSurfaceBcf10();
    }

    int verboseLevel =
                T4SettingsFile::getInstance()->getStructManager()->getGeneral()->verboseLevel;
    // Print material table and element table
    if(verboseLevel > 0)
        G4cout << *(G4Material::GetMaterialTable()) << endl;
    if(verboseLevel > 2)
        G4cout << *(G4Element::GetElementTable()) << endl;
}

void T4Materials::DefineIndex(void)
{
    bc408Index = 1.58;
//  bc408AbsorptionLength = 3.80 * CLHEP::m; // bulk light attenuation length
    bc408AbsorptionLength = 1.80 * CLHEP::m;     // bulk light attenuation length
    bc408Yield = 0.64 * 2 / (100.0 * CLHEP::eV);     // 64% anthracene
    plexiglassIndex = 1.49;
    plexiglassAbsLength = 10.0 * CLHEP::m;
    glassIndex = 1.474;
    aluminiumReflectivity = 0.97;
    bcf10Index = 1.60;
    pmmaIndex = 1.49;

    specularLobe[0] = 1.0;
    specularLobe[1] = 1.0;
    specularSpike[0] = 0.0;
    specularSpike[1] = 0.0;
    backScatter[0] = 0.0;
    backScatter[1] = 0.0;

    photonEnergyShort[0] = 2.0 * CLHEP::eV;     // 619.92 nm
    photonEnergyShort[1] = 7.14 * CLHEP::eV;     // 173.65 nm
    photonEnergyLong[0] = 2.08 * CLHEP::eV;     // 596.08 nm
    photonEnergyLong[1] = 2.38 * CLHEP::eV;     // 520.94 nm
    photonEnergyLong[2] = 2.58 * CLHEP::eV;     // 480.56 nm
    photonEnergyLong[3] = 2.70 * CLHEP::eV;     // 459.20 nm
    photonEnergyLong[4] = 2.76 * CLHEP::eV;     // 449.22 nm
    photonEnergyLong[5] = 2.82 * CLHEP::eV;     // 439.66 nm
    photonEnergyLong[6] = 2.92 * CLHEP::eV;     // 424.60 nm
    photonEnergyLong[7] = 2.95 * CLHEP::eV;     // 420.29 nm
    photonEnergyLong[8] = 3.02 * CLHEP::eV;     // 410.54 nm
    photonEnergyLong[9] = 3.10 * CLHEP::eV;     // 399.95 nm
    photonEnergyLong[10] = 3.26 * CLHEP::eV;     // 380.32 nm
    photonEnergyLong[11] = 3.44 * CLHEP::eV;     // 360.42 nm

    airRefractiveIndex[0] = 1.0;
    airRefractiveIndex[1] = 1.0;
    airAbsorpiton[0] = 10.0 * CLHEP::m;
    airAbsorpiton[1] = 10.0 * CLHEP::m;
    plexiglassRefractiveIndex[0] = plexiglassIndex;
    plexiglassRefractiveIndex[1] = plexiglassIndex;
    plexiglassAbsorpiton[0] = plexiglassAbsLength;
    plexiglassAbsorpiton[1] = plexiglassAbsLength;
    glassRefractiveIndex[0] = glassIndex;
    glassRefractiveIndex[1] = glassIndex;
    bc408RefractiveIndex[0] = bc408Index;
    bc408RefractiveIndex[1] = bc408Index;
    bcf10RefractiveIndex[0] = bcf10Index;
    bcf10RefractiveIndex[1] = bcf10Index;
    pmmaRefractiveIndex[0] = pmmaIndex;
    pmmaRefractiveIndex[1] = pmmaIndex;
        
}

void T4Materials::BuildAir(void)
{
    air_noOptical = manager->FindOrBuildMaterial("G4_AIR");
    air_noOptical->SetName("air_noOptical");
/*    G4Material* N2 = new G4Material("N2", 0.808 * CLHEP::g / CLHEP::mL, 1);
    N2->AddElement(elN, 2);
    G4Material* O2= new G4Material("O2", 1.429 * CLHEP::g / CLHEP::cm3, 1);
    O2->AddElement(elO, 2);
    double density = 1.290*CLHEP::mg/CLHEP::cm3;
    G4Material* air = new G4Material("Air  ", density, 2);
    air->AddMaterial(O2, 0.3);
    air->AddMaterial(N2, 0.7);*/

    vacuum_noOptical = manager->FindOrBuildMaterial("G4_Galactic");
    vacuum_noOptical->SetName("vacuum_noOptical");

    if (useOptical) {
        air_optical = new G4Material("air_optical",
                                     0.00120479 * CLHEP::g / CLHEP::cm3, 1);
        air_optical->AddMaterial(air_noOptical, 1.0);

        vacuum_optical = new G4Material("vacuum_optical",
                                        1e-25 * CLHEP::g / CLHEP::cm3, 1);
        vacuum_optical->AddMaterial(vacuum_noOptical, 1.0);

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                    airRefractiveIndex, nEntryShort);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", photonEnergyShort,
                                                    airAbsorpiton, nEntryShort);
        air_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
        vacuum_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        air_optical = air_noOptical;
    }
}

void T4Materials::BuildBc408(void)
{
    bc408_noOptical = new G4Material("bc408_noOptical",
                                     1.032 * CLHEP::g / CLHEP::cm3, 2);
    bc408_noOptical->AddElement(elC, 9);
    bc408_noOptical->AddElement(elH, 10);

    if (useOptical) {
        bc408_optical = new G4Material("bc408_optical", 1.032 * CLHEP::g / CLHEP::cm3,
                                       1);
        bc408_optical->AddMaterial(bc408_noOptical, 1.0);

        G4double _bc408RefractiveIndex[nEntryLong];
        G4double bc408Absorpiton[nEntryLong];
        for (G4int i = 0; i < nEntryLong; i++) {
            _bc408RefractiveIndex[i] = bc408Index;
            bc408Absorpiton[i] = bc408AbsorptionLength;
        }
        G4double bc408Scintillation[nEntryLong] = { 0.01, 0.03, 0.17, 0.40, 0.55,
                                                    0.83, 1.00, 0.84, 0.49, 0.20, 0.07, 0.04 };

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyLong,
                                                    _bc408RefractiveIndex, nEntryLong);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", photonEnergyLong,
                                                    bc408Absorpiton, nEntryLong);
        materialPropertiesTable.back()->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergyLong,
                                                    bc408Scintillation, nEntryLong);
        materialPropertiesTable.back()->AddConstProperty("SCINTILLATIONYIELD",
                                                         bc408Yield);
        materialPropertiesTable.back()->AddConstProperty("RESOLUTIONSCALE", 1.0);
        materialPropertiesTable.back()->AddConstProperty("SCINTILLATIONTIMECONSTANT1",
                                                         0.9 * CLHEP::ns);         // rise time
        materialPropertiesTable.back()->AddConstProperty("SCINTILLATIONTIMECONSTANT2",
                                                         2.1 * CLHEP::ns);         // decay time
        materialPropertiesTable.back()->AddConstProperty("SCINTILLATIONYIELD1", 0.27);

        bc408_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
        bc408_optical->GetIonisation()->SetBirksConstant(
            0.126 * CLHEP::mm / CLHEP::MeV);
    } else {
        bc408_optical = bc408_noOptical;
    }
}

void T4Materials::BuildBcf10(void)
{
    bcf10_noOptical = new G4Material("bcf10_noOptical", 1.05 * CLHEP::g / CLHEP::cm3,
                                     2);
    bcf10_noOptical->AddElement(elC, 8);
    bcf10_noOptical->AddElement(elH, 8);

    if (useOptical) {
        bcf10_optical = new G4Material("bcf10_optical", 1.05 * CLHEP::g / CLHEP::cm3,
                                       1);
        bcf10_optical->AddMaterial(bcf10_noOptical, 1.0);

        G4double _bcf10RefractiveIndex[nEntryLong];
        G4double bcf10Absorpiton[nEntryLong];
        for (G4int i = 0; i < nEntryLong; i++) {
            _bcf10RefractiveIndex[i] = bcf10Index;
            bcf10Absorpiton[i] = 190.0 * CLHEP::cm;
        }
        G4double bcf10Scintillation[nEntryLong] = { 0.001, 0.13, 0.42, 0.76, 0.85,
                                                    1.00, 0.85, 0.76, 0.38, 0.09, 0.001, 0.0001 };

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyLong,
                                                    _bcf10RefractiveIndex, nEntryLong);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", photonEnergyLong,
                                                    bcf10Absorpiton, nEntryLong);
        materialPropertiesTable.back()->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergyLong,
                                                    bcf10Scintillation, nEntryLong);
        G4double bcf10Yield = 8000.0 / CLHEP::MeV;
        materialPropertiesTable.back()->AddConstProperty("SCINTILLATIONYIELD",
                                                         bcf10Yield);
        materialPropertiesTable.back()->AddConstProperty("RESOLUTIONSCALE", 1.0);
        materialPropertiesTable.back()->AddConstProperty("SCINTILLATIONTIMECONSTANT1",
                                                         0.9 * CLHEP::ns);         // rise time
        materialPropertiesTable.back()->AddConstProperty("SCINTILLATIONTIMECONSTANT2",
                                                         2.7 * CLHEP::ns);         // decay time

        bcf10_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        bcf10_optical = bcf10_noOptical;
    }
}

void T4Materials::BuildPMMA(void)
{
    pmma_noOptical = new G4Material("pmma_noOptical", 1.18 * CLHEP::g / CLHEP::cm3,
                                    3);
    pmma_noOptical->AddElement(elC, 5);
    pmma_noOptical->AddElement(elO, 2);
    pmma_noOptical->AddElement(elH, 8);

    if (useOptical) {
        pmma_optical = new G4Material("pmma_optical", 1.18 * CLHEP::g / CLHEP::cm3,
                                      1);
        pmma_optical->AddMaterial(pmma_noOptical, 1.0);

        G4double pmmaAbsorpiton[nEntryShort];
        for (G4int i = 0; i < nEntryShort; i++) {
            pmmaAbsorpiton[i] = 190.0 * CLHEP::cm;             // ?
        }

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                    pmmaRefractiveIndex, nEntryShort);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", photonEnergyShort,
                                                    pmmaAbsorpiton, nEntryShort);
        pmma_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        pmma_optical = pmma_noOptical;
    }
}

void T4Materials::BuildPlexiglass(void)
{
    plexiglass_noOptical = manager->FindOrBuildMaterial("G4_PLEXIGLASS");
    plexiglass_noOptical->SetName("plexiglass_noOptical");

    if (useOptical) {
        plexiglass_optical = new G4Material("plexiglass_optical",
                                            1.19 * CLHEP::g / CLHEP::cm3, 1);
        plexiglass_optical->AddMaterial(plexiglass_noOptical, 1.0);

        G4double plexiglassEnergy[nEntryShort] = {2.0 * CLHEP::eV, 7.14 * CLHEP::eV};
        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", plexiglassEnergy,
                                                    plexiglassRefractiveIndex, nEntryShort);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", plexiglassEnergy,
                                                    plexiglassAbsorpiton, nEntryShort);
        plexiglass_optical->SetMaterialPropertiesTable(
            materialPropertiesTable.back());
    } else {
        plexiglass_optical = plexiglass_noOptical;
    }
}

void T4Materials::BuildGlass(void)
{
    glass_noOptical = manager->FindOrBuildMaterial("G4_GLASS_PLATE");
    glass_noOptical->SetName("glass_noOptical");

    if (useOptical) {
        glass_optical = new G4Material("glass_optical", 2.4 * CLHEP::g / CLHEP::cm3,
                                       1);
        glass_optical->AddMaterial(glass_noOptical, 1.0);

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                    glassRefractiveIndex, nEntryShort);
        glass_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        glass_optical = glass_noOptical;
    }
}

void T4Materials::BuildAluminium(void)
{
    aluminium_noOptical = manager->FindOrBuildMaterial("G4_Al");
    aluminium_noOptical->SetName("aluminium_noOptical");

    if (useOptical) {
        aluminium_optical = new G4Material("aluminium_optical",
                                           2.699 * CLHEP::g / CLHEP::cm3, 1);
        aluminium_optical->AddMaterial(aluminium_noOptical, 1.0);

        G4double aluminiumEnergy[nEntryShort] = { 2.0 * CLHEP::eV, 3.5 * CLHEP::eV };
        G4double aluminiumRefractiveIndex[nEntryShort] = { 1.51, 1.61 };
        G4double aluminiumAbsorpiton[nEntryShort] = { 1.0e-20 * CLHEP::m, 1.0e-20
                                                      * CLHEP::m };
        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", aluminiumEnergy,
                                                    aluminiumRefractiveIndex, nEntryShort);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", aluminiumEnergy,
                                                    aluminiumAbsorpiton, nEntryShort);
        aluminium_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        aluminium_optical = aluminium_noOptical;
    }
}

void T4Materials::BuildBialkali(void)
{
    bialkali_noOptical = new G4Material("bialkali_noOptical",
                                        3.0 * CLHEP::g / CLHEP::cm3, 3, kStateSolid);
    bialkali_noOptical->AddElement(elSb, 1);
    bialkali_noOptical->AddElement(elRb, 1);
    bialkali_noOptical->AddElement(elCs, 1);

    if (useOptical) {
        bialkali_optical = new G4Material("bialkali_optical",
                                          3.0 * CLHEP::g / CLHEP::cm3, 1);
        bialkali_optical->AddMaterial(bialkali_noOptical, 1.0);

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        G4double bialkaliAbsorption[nEntryShort] = { 1.0e-6 * CLHEP::mm, 1.0e-6
                                                     * CLHEP::mm };
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", photonEnergyShort,
                                                    bialkaliAbsorption, nEntryShort);
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                    glassRefractiveIndex, nEntryShort);
        bialkali_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        bialkali_optical = bialkali_noOptical;
    }
}

void T4Materials::BuildC4F10(void)
{
    c4f10_noOptical = new G4Material("c4f10_noOptical",
                                     9.912 * CLHEP::kg / CLHEP::m3, 2, kStateGas, 298 * CLHEP::kelvin,
                                     1.0 * CLHEP::atmosphere);
    c4f10_noOptical->AddElement(elC, 4);
    c4f10_noOptical->AddElement(elF, 10);

    if (useOptical) {
        c4f10_optical = new G4Material("c4f10_optical", 9.912 * CLHEP::kg / CLHEP::m3,
                                       1);
        c4f10_optical->AddMaterial(c4f10_noOptical, 1.0);

        // see COMPASS note 2014-2
        // Studies on RICH MC, Federica Sozzi
        // .x evaluate_CGcards_RICHindex.C(1.001309, 1.001451)
        G4double photonEnergy[10] = {1.00, 2.00, 3.00, 4.00, 4.50, 5.00, 5.50, 5.99, 6.00, 8.00};
        G4double c4f10RefractiveIndex[10] = {1.001171, 1.001203, 1.001225, 1.001257, 1.001277, 1.001301, 1.001328, 1.001358, 1.001359, 1.001533};
        G4double c4f10Absorpiton[10];

        for (G4int i = 0; i < 10; i++) {
            photonEnergy[i] *= CLHEP::eV;
            c4f10Absorpiton[i] = 10.00 * CLHEP::m;
        }

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergy,
                                                    c4f10RefractiveIndex, 10);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", photonEnergy,
                                                    c4f10Absorpiton, 10);
        c4f10_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        c4f10_optical = c4f10_noOptical;
    }
}

void T4Materials::BuildN2(void)
{
    n2_noOptical = manager->FindOrBuildMaterial("G4_N");
    n2_noOptical->SetName("n2_noOptical");

    if (useOptical) {
        n2_optical = new G4Material("n2_optical", 0.0011652 * CLHEP::g / CLHEP::cm3, 1);
        n2_optical->AddMaterial(n2_noOptical, 1.0);

        // http://refractiveindex.info/?shelf=main&book=N2&page=Griesmann
        G4double photonEnergy[6] = {1.0072, 1.9774, 5.0095, 6.0099, 6.985, 7.999};
        G4double n2RefractiveIndex[6] = {1.0002794, 1.0002823, 1.00032265, 1.00033786, 1.00035797, 1.0003868};
        G4double n2Absorpiton[6];

        for (G4int i = 0; i < 6; i++) {
            photonEnergy[i] *= CLHEP::eV;
            n2Absorpiton[i] = 10.00 * CLHEP::m;
        }

        materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
        materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergy,
                                                    n2RefractiveIndex, 6);
        materialPropertiesTable.back()->AddProperty("ABSLENGTH", photonEnergy,
                                                    n2Absorpiton, 6);
        n2_optical->SetMaterialPropertiesTable(materialPropertiesTable.back());
    } else {
        n2_optical = n2_noOptical;
    }
}

void T4Materials::BuildSurfaceBc408Air(void)
{
    surfaceBc408Air = new G4OpticalSurface("bc408AirSurface");
    surfaceBc408Air->SetType(dielectric_dielectric);
    surfaceBc408Air->SetModel(unified);
    surfaceBc408Air->SetFinish(ground);
    surfaceBc408Air->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                airRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfaceBc408Air->SetMaterialPropertiesTable(materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceAirBc408(void)
{
    surfaceAirBc408 = new G4OpticalSurface("airBc408Surface");
    surfaceAirBc408->SetType(dielectric_dielectric);
    surfaceAirBc408->SetModel(unified);
    surfaceAirBc408->SetFinish(ground);
    surfaceAirBc408->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                bc408RefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfaceAirBc408->SetMaterialPropertiesTable(materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceBc408Plexiglass(void)
{
    surfaceBc408Plexiglass = new G4OpticalSurface("bc408PlexiglassSurface");
    surfaceBc408Plexiglass->SetType(dielectric_dielectric);
    surfaceBc408Plexiglass->SetModel(unified);
    surfaceBc408Plexiglass->SetFinish(ground);
    surfaceBc408Plexiglass->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                plexiglassRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfaceBc408Plexiglass->SetMaterialPropertiesTable(
        materialPropertiesTable.back());
}

void T4Materials::BuildSurfacePlexiglassGlass(void)
{
    surfacePlexiglassGlass = new G4OpticalSurface("plexiglassGlassSurface");
    surfacePlexiglassGlass->SetType(dielectric_dielectric);
    surfacePlexiglassGlass->SetModel(unified);
    surfacePlexiglassGlass->SetFinish(ground);
    surfacePlexiglassGlass->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                glassRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfacePlexiglassGlass->SetMaterialPropertiesTable(
        materialPropertiesTable.back());
}

void T4Materials::BuildSurfacePlexiglassAir(void)
{
    surfacePlexiglassAir = new G4OpticalSurface("plexiglassAirSurface");
    surfacePlexiglassAir->SetType(dielectric_dielectric);
    surfacePlexiglassAir->SetModel(unified);
    surfacePlexiglassAir->SetFinish(ground);
    surfacePlexiglassAir->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                airRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfacePlexiglassAir->SetMaterialPropertiesTable(
        materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceAirPlexiglass(void)
{
    surfaceAirPlexiglass = new G4OpticalSurface("airPlexiglassSurface");
    surfaceAirPlexiglass->SetType(dielectric_dielectric);
    surfaceAirPlexiglass->SetModel(unified);
    surfaceAirPlexiglass->SetFinish(ground);
    surfaceAirPlexiglass->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                plexiglassRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfaceAirPlexiglass->SetMaterialPropertiesTable(
        materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceGlassVacuum(void)
{
    surfaceGlassVacuum = new G4OpticalSurface("glassVacuumSurface");
    surfaceGlassVacuum->SetType(dielectric_dielectric);
    surfaceGlassVacuum->SetModel(unified);
    surfaceGlassVacuum->SetFinish(ground);     //groundair
    surfaceGlassVacuum->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                airRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfaceGlassVacuum->SetMaterialPropertiesTable(
        materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceVacuumBialkali(void)
{
    surfaceVacuumBialkali = new G4OpticalSurface("vacuumBialkaliSurface");
    surfaceVacuumBialkali->SetType(dielectric_dielectric);
    surfaceVacuumBialkali->SetModel(unified);
    surfaceVacuumBialkali->SetFinish(ground);
    surfaceVacuumBialkali->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                glassRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfaceVacuumBialkali->SetMaterialPropertiesTable(
        materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceAluminium(void)
{
    surfaceAluminium = new G4OpticalSurface("aluminiumSurface");
    surfaceAluminium->SetType(dielectric_metal);
    surfaceAluminium->SetModel(unified);
    surfaceAluminium->SetFinish(ground);
    surfaceAluminium->SetPolish(0.0);

    G4double aluminiumReflectivityShort[nEntryShort] = { aluminiumReflectivity,
                                                         aluminiumReflectivity };
    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("REFLECTIVITY", photonEnergyShort,
                                                aluminiumReflectivityShort, nEntryShort);
    surfaceAluminium->SetMaterialPropertiesTable(materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceMirror(void)
{
    surfaceMirror = new G4OpticalSurface("mirrorSurface");
    surfaceMirror->SetModel(unified);
    surfaceMirror->SetType(dielectric_dielectric);
    surfaceMirror->SetFinish(polishedfrontpainted);

    G4double surfaceReflectivityShort[nEntryShort] = { 0.86, 0.86 };
    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("REFLECTIVITY", photonEnergyShort,
                                                surfaceReflectivityShort, nEntryShort);
    surfaceMirror->SetMaterialPropertiesTable(materialPropertiesTable.back());
}

void T4Materials::BuildSurfaceBcf10(void)
{
    surfaceBcf10 = new G4OpticalSurface("surfaceBcf10");
    surfaceBcf10->SetType(dielectric_dielectric);
    surfaceBcf10->SetModel(unified);
    surfaceBcf10->SetFinish(ground);
    surfaceBcf10->SetSigmaAlpha(0.035);

    materialPropertiesTable.push_back(new G4MaterialPropertiesTable());
    materialPropertiesTable.back()->AddProperty("RINDEX", photonEnergyShort,
                                                pmmaRefractiveIndex, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARLOBECONSTANT",
                                                photonEnergyShort, specularLobe, nEntryShort);
    materialPropertiesTable.back()->AddProperty("SPECULARSPIKECONSTANT",
                                                photonEnergyShort, specularSpike, nEntryShort);
    materialPropertiesTable.back()->AddProperty("BACKSCATTERCONSTANT",
                                                photonEnergyShort, backScatter, nEntryShort);
    surfaceBcf10->SetMaterialPropertiesTable(materialPropertiesTable.back());
}

T4Materials::~T4Materials(void)
{
    if (useOptical) {
        delete air_optical;
        delete vacuum_optical;
        delete bc408_optical;
        delete bcf10_optical;
        delete pmma_optical;
        delete plexiglass_optical;
        delete glass_optical;
        delete aluminium_optical;
        delete bialkali_optical;
        delete c4f10_optical;
    }

    delete air_noOptical;
    delete vacuum_noOptical;
    delete bc408_noOptical;
    delete bcf10_noOptical;
    delete pmma_noOptical;
    delete plexiglass_noOptical;
    delete glass_noOptical;
    delete aluminium_noOptical;
    delete bialkali_noOptical;
    delete c4f10_noOptical;

    delete c2h6;
    delete plasticFoam;
    delete carbonFibre;
    delete carbonLH2;
    delete nomex;
    delete epoxy;
    delete stainlessSteel;
    delete airexC70;
    delete lhe;
    delete ammonia;
    if(nh3_he) delete nh3_he;
    delete nh3_he_1st_2014;
    delete nh3_he_2nd_2014;
    delete nh3_he_1st_2015;
    delete nh3_he_2nd_2015;
    delete nh3_he_1st_2018;
    delete nh3_he_2nd_2018;
    delete C2F3Cl;
    delete elCl;
    delete nh3_he_fluka;
    delete dli_he_2021;
    delete elLi;
    delete Li2CO3;
    delete AlMylar;
    delete cf4;
    delete g10;
    delete g10_10;
    delete g11;
    delete ArCo2;
    delete dc4gas;
    delete dcEBox;
    delete strawGas;
    delete w45Gas;
    delete mdtGas;
    delete mw2Gas;
    delete kapton;
    delete Cu_10;
    delete Cu_Kap;
    delete rohacell;
    delete rohacell_30;
    delete inox;
    delete gemMaterial;
    delete mwpcGas;
    delete mmSupport;
    delete pixelGEMMaterial;
    delete PbO;
    delete SiO2;
    delete K2O;
    delete As2O3;
    delete CeO2;
    delete tf1;
    delete tf101;
    delete sf57;
    delete sf5;
    delete paper;

    if (useOptical) {
        delete surfaceBc408Air;
        delete surfaceAirBc408;
        delete surfaceBc408Plexiglass;
        delete surfacePlexiglassGlass;
        delete surfacePlexiglassAir;
        delete surfaceAirPlexiglass;
        delete surfaceGlassVacuum;
        delete surfaceVacuumBialkali;
        delete surfaceAluminium;
        delete surfaceMirror;
        delete surfaceBcf10;

        for (unsigned int i = 0; i < materialPropertiesTable.size(); i++)
            delete materialPropertiesTable.at(i);
        materialPropertiesTable.clear();
    }
}
