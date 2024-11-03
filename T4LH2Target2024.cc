#include "T4LH2Target2024.hh"
#include <iostream>
#include "TMath.h"


T4LH2Target2024::T4LH2Target2024(T4SDetector* lh2)
{
  yearSetup = strToDouble(T4EventManager::getInstance()->getTriggerPlugin()->getPluginName().substr(0,4));

  setPosition(lh2->position);
  useMechanicalStructure = lh2->useMechanicalStructure;

 
  // all length measurements must be half of the actual length

  cavityRadius1 = 318 / 2 * CLHEP::mm;
  cavityThickness1 = 3 * CLHEP::mm;
  cavityLength1 =  1068 / 2 * CLHEP::mm;

  cavityRadius2 = 164.3 / 2 * CLHEP::mm;
  cavityThickness2 = 2 * CLHEP::mm;
  cavityLength2 = 338 / 2 * CLHEP::mm;

  cavityRadius3 =  163.0 / 2 * CLHEP::mm;
  cavityThickness3 = 2.6 * CLHEP::mm;
  cavityLength3 =  335 / 2 * CLHEP::mm;

  cavityRadius4 = 163.0 / 2 * CLHEP::mm;
  cavityThickness4 = 2.6 * CLHEP::mm;
  cavityLength4 = 286 / 2 * CLHEP::mm;
  
  cavityRadius4_1 = 163.0 / 2 * CLHEP::mm;
  cavityThickness4_1 = 2.6 * CLHEP::mm;
  cavityLength4_1 = 60 / 2 * CLHEP::mm;

  EndLength1 = 17. / 2 * CLHEP::mm;
  EndRadiusMin1 = 318. / 2 * CLHEP::mm;
  EndRadiusMax1 = 370. / 2 * CLHEP::mm;

  EndLength_O = 5.6 / 2 * CLHEP::mm;
  EndRadiusMin_O =  341 / 2 * CLHEP::mm;
  EndRadiusMax_O = 358 / 2 * CLHEP::mm;

  SideLength1 = 11.5 / 2 * CLHEP::mm;
  SideRadiusMin1 = 80. / 2 * CLHEP::mm;
  SideRadiusMax1 = 139.5 / 2 * CLHEP::mm;

  SideLength2 = 15.5 / 2 * CLHEP::mm;
  SideRadiusMin2 = 80. / 2 * CLHEP::mm;
  SideRadiusMax2 = 219.5 / 2 * CLHEP::mm;

  JunctionLength12_1 = 17. / 2 * CLHEP::mm;
  JunctionRadiusMin12_1 = 168. / 2 * CLHEP::mm;
  JunctionRadiusMax12_1 = 370. / 2 * CLHEP::mm;

  JunctionLength12_O = 5.6 / 2 * CLHEP::mm;   
  JunctionRadiusMin12_O = 313.6 / 2 * CLHEP::mm;
  JunctionRadiusMax12_O = 358 / 2 * CLHEP::mm;

  JunctionLength12_2 = 17. / 2 * CLHEP::mm;
  JunctionRadiusMin12_2 = 318 / 2 * CLHEP::mm;
  JunctionRadiusMax12_2 = 370 / 2 * CLHEP::mm;

  JunctionLength23_1 =  15.0 / 2 * CLHEP::mm;
  JunctionRadiusMin23_1 =  110 / 2 * CLHEP::mm;
  JunctionRadiusMax23_1 =  239.5 / 2 * CLHEP::mm;

  JunctionLength23_2 = 14.0 / 2 * CLHEP::mm;      // junction O-ring parameters between C2 and C3
  JunctionRadiusMin23_2 = 164.3 / 2 * CLHEP::mm;
  JunctionRadiusMax23_2 = 239 / 2 * CLHEP::mm;

  shift1 = 29. * CLHEP::mm;
  shift2 = 103. * CLHEP::mm;                              
  dz1 = 480. * CLHEP::mm; 
  dz2 = 204. * CLHEP::mm; 

  targetRadius =  39.8 / 2 * CLHEP::mm;                   
  targetLength = 1400 / 2 * CLHEP::mm;                   // already contains in its definition the 2 radii
  mylarRadius =  40 / 2 * CLHEP::mm;                     // for mylar semi-spheres end caps
  mylarThickness =  0.125 * CLHEP::mm;                   // 125 um
  mylarCylLength = 6.2 / 2 * CLHEP::mm;

  trapThickness = 20. / 2 * CLHEP::mm;                     // target holder thickness
  trapDx = 12.4 * CLHEP::mm;
  trapSide = 269.6 * CLHEP::mm;
  trapHeight = 244.4 * CLHEP::mm;
  HoleD = 51.5 * CLHEP::mm;
  yBox =  14.5 / 2 * CLHEP::mm;
  xBox = 30.4 / 2 * CLHEP::mm;

  mylarWindowRadius = 370. / 2 * CLHEP::mm;
  mylarWindowThick = 0.35 * CLHEP::mm;
  
  kaptonLength =  1265 / 2 * CLHEP::mm;
  kaptonRadius =  40 / 2 * CLHEP::mm;
  kaptonThickness = 0.125 * CLHEP::mm;

  steelRadius = 39.8 / 2 * CLHEP::mm;         
  steelThickness = 0.4 * CLHEP::mm;
  steelLength = 95 / 2 * CLHEP::mm;

  // caseThickness1 = 8.5 * CLHEP::mm;
  // caseThickness1_1 = 3.5 * CLHEP::mm;
  // caseThickness2 = 7.5 / 2 * CLHEP::mm;
  caseThickness = 30 * 11 * CLHEP::um;
  caseLength1 = 1347.5 / 2 * CLHEP::mm;
  //caseLength1_1 = 16. / 2 * CLHEP::mm;
  //caseCone = 7. / 2 * CLHEP::mm;
  //caseLength1 = 120. / 2 * CLHEP::mm;
  //caseLength2 = 468. / 2 * CLHEP::mm;
  //caseLength3 = 468. / 2 * CLHEP::mm;
  //caseLength4 = 194. / 2 * CLHEP::mm;
  caseLength2 = 20. / 2 * CLHEP::mm;

  rotate1 = new G4RotationMatrix();
  rotate1->rotateX(90 * CLHEP::deg);
  rotate2 = new G4RotationMatrix();
  rotate2->rotateX(-90 * CLHEP::deg);
  
  rotate120 = new G4RotationMatrix();
	rotate120->rotateZ(120 * CLHEP::degree);

  rotate60 = new G4RotationMatrix();
	rotate60->rotateZ(60 * CLHEP::degree);

  setDimension(positionVector.z() - targetLength, positionVector.z() + targetLength);
  setELossParams(0., 0.);
}

void T4LH2Target2024::construct(G4LogicalVolume* world_log)
{

  // _____________________________________ biggest stainless-steel cylinder #1 _____________________________________

    double posZ_C1 = shift1 + shift2 + JunctionLength12_1 + 2*JunctionLength12_O;
    tubs.push_back(new G4Tubs("C1", cavityRadius1, cavityRadius1 + cavityThickness1, cavityLength1, 0., 2.*M_PI));
    // logical volume
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0.,posZ_C1), log.back(), "Cavity1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    // junctions between C1 and C2 
    double posZ_C1j2 = posZ_C1 - cavityLength1;
    tubs.push_back(new G4Tubs("C1j2", JunctionRadiusMin12_2, JunctionRadiusMax12_2, JunctionLength12_2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C1j2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C1j2), log.back(), "Cavity1j2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->orange);

    double posZ_O12 = posZ_C1j2 - JunctionLength12_2 - JunctionLength12_O;
    tubs.push_back(new G4Tubs("C1_O", JunctionRadiusMin12_O, JunctionRadiusMax12_O, JunctionLength12_O, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C1_O_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_O12), log.back(), "Cavity1_O", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->gray);
    
    double posZ_C1j1 = posZ_O12 - JunctionLength12_1 - JunctionLength12_O;
    tubs.push_back(new G4Tubs("C1j1", JunctionRadiusMin12_1, JunctionRadiusMax12_1, JunctionLength12_1, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C1j1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C1j1), log.back(), "Cavity1j1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->yellow_50);

  // _____________________________________ small stainless-steel cylinder #2 _____________________________________

    double posZ_C2 = posZ_C1j1 - cavityLength2;
    tubs.push_back(new G4Tubs("C2", cavityRadius2, cavityRadius2 + cavityThickness2, cavityLength2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C2), log.back(), "Cavity2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    // junctions between #2 and #3
    double posZ_C2j2 = posZ_C2 - cavityLength2 - JunctionLength23_2;
    tubs.push_back(new G4Tubs("C2j2", JunctionRadiusMin23_2, JunctionRadiusMax23_2, JunctionLength23_2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C2j2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C2j2 ), log.back(), "Cavity2j2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgreen);

    double posZ_C2j1 = posZ_C2j2 - JunctionLength23_2 - JunctionLength23_1;
    tubs.push_back(new G4Tubs("C2j1", JunctionRadiusMin23_1, JunctionRadiusMax23_1, JunctionLength23_1, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C2j1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C2j1), log.back(), "Cavity2j1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->lightgreen);

  // _____________________________________ #3 stainless-steel cylinder _____________________________________

    double posZ_C3 = posZ_C2j1 - JunctionLength23_1 - cavityLength3;
    tubs.push_back(new G4Tubs("C3", cavityRadius3, cavityRadius3 + cavityThickness3, cavityLength3, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C3_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C3), log.back(), "Cavity3", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

  // _______________________________________ mother log volumes for cylinder 4 _______________________________________
   
    double posZ_motherC4 = posZ_C3 - cavityLength3 - cavityLength4;
    tubs.push_back(new G4Tubs("motherC4_vol", 0., 2*cavityRadius4, 1.5*cavityLength4, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "motherC4_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_motherC4), log.back(), "MotherC4", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    MotherC4_log = log.back();

  // _____________________________________ #4 stainless-steel cylinder _____________________________________
    
    tubs.push_back(new G4Tubs("C4", cavityRadius4, cavityRadius4 + cavityThickness4, cavityLength4, 0., 2.*M_PI));      

    // disc to create a hole in Cylinder 4 along y axis
    G4Tubs* disc4 = new G4Tubs("disc4", 0., cavityRadius4_1 + cavityThickness4_1, cavityLength4_1, 0., 2.*M_PI);
    G4ThreeVector Pos_C4_1(0.,cavityRadius4,0.);

    subtractionC4 = new G4SubtractionSolid("C4-disc4", tubs.back(), disc4, 0, Pos_C4_1);
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "subtractionC4_log"));
    new G4PVPlacement(0, positionVector, log.back(), "SubtractionC4", MotherC4_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    tubs.push_back(new G4Tubs("C4_1", cavityRadius4_1, cavityRadius4_1 + cavityThickness4_1, cavityLength4_1, 0., 2.*M_PI));      
    //log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4_1_log"));
    //new G4PVPlacement(rotate1, positionVector + G4ThreeVector(0.,cavityRadius4,0.), log.back(), "Cavity4_1", MotherC4_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->darkgray_50);

   // unionC4 = new G4UnionSolid("subtractionC4+C4_1",subtractionC4,tubs.back(),rotate1,G4ThreeVector(0.,cavityRadius4,0.));                              
   // log.push_back(new G4LogicalVolume(union1, materials->stainlessSteel, "UnionC4_log"));
   // new G4PVPlacement(0, positionVector, log.back(), "UnionC4", MotherC4_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->darkgray_50);

    // O-rings at the beginning of C4

    tubs.push_back(new G4Tubs("C4_side3", cavityRadius4, SideRadiusMax2, SideLength2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C4Side3_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., -cavityLength4 - SideLength2), log.back(), "Cavity4Side3", MotherC4_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->cyan);

    tubs.push_back(new G4Tubs("C4_side2", SideRadiusMin2, SideRadiusMax2, SideLength2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C4Side2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., -cavityLength4 - 3*SideLength2), log.back(), "Cavity4Side2", MotherC4_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->dodgerblue);
    
    tubs.push_back(new G4Tubs("C4_side1", SideRadiusMin1, SideRadiusMax1, SideLength1, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C4Side1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., -cavityLength4 - 4*SideLength2 - SideLength1), log.back(), "Cavity4Side1", MotherC4_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->orange);

  // ________________________________________ mylar window _________________________________________

    double posZ_MW = targetLength + shift2 + EndLength1 + 2*EndLength_O + mylarWindowThick;
    tubs.push_back(new G4Tubs("mylarWin", 0., mylarWindowRadius, mylarWindowThick, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "MW_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MW), log.back(), "MylarWindow", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->blue);

    // O-rings for mylar window
    double posZ_MWj1 = targetLength + shift2;
    tubs.push_back(new G4Tubs("mylarWin_j1", EndRadiusMin1, EndRadiusMax1, EndLength1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "MWj1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MWj1), log.back(), "MylarWindow_j1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->mediumpurple);

    double posZ_MW_O = posZ_MWj1 + EndLength1 + EndLength_O;
    tubs.push_back(new G4Tubs("mylarWin_O", EndRadiusMin_O, EndRadiusMax_O, EndLength_O, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "MW_O_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MW_O), log.back(), "MylarWindow_O", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->gray);

    double posZ_MWj2 = posZ_MW + mylarWindowThick + EndLength1;
    tubs.push_back(new G4Tubs("mylarWin_j2", EndRadiusMin1, EndRadiusMax1, EndLength1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "MWj2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MWj2), log.back(), "MylarWindow_j2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->magenta);

  // _______________________________________ mother log volumes for target caps and insulation _______________________________________
    
    // #1 referrs to the end of target, closest to mylar window
    // dimensions of mother volume don't matter much, just ensure they're big enough to contain all the solids
    tubs.push_back(new G4Tubs("mother1_vol", 0., 2*targetRadius + 2*caseThickness, 2*targetLength, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "mother1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength - 2*targetRadius), log.back(), "MotherVolLog_1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    CapMother1_log = log.back();

    tubs.push_back(new G4Tubs("mother2_vol", 0., 2*targetRadius + 2*caseThickness, targetLength, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "mother2_log"));
    new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,targetLength - 2*targetRadius), log.back(), "MotherVolLog_2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    CapMother2_log = log.back();

  // ___________________________________________ H2 target ___________________________________________

   // cylindrical part
    tubs.push_back(new G4Tubs("target_H2", 0., targetRadius, targetLength - 2*targetRadius, 0., 2.*M_PI));
    //log.push_back(new G4LogicalVolume(tubs.back(), materials->lh2, "targetH2_log"));
    //addToTarget(new G4PVPlacement(0, positionVector, log.back(), "LH2Hadron", world_log, 0, 0, checkOverlap));
    //log.back()->SetVisAttributes(colour->lightblue);

    // adding first semisphere to cylindrical lh2 target     
    sphere.push_back(new G4Sphere("LH2Sphere_1", 0., targetRadius, 0., M_PI, 0., M_PI));
    //log.push_back(new G4LogicalVolume(sphere.back(), materials->lh2, "LH2BegCap_log"));
    //new G4PVPlacement(rotate2, positionVector , log.back(), "LH2BegCap", CapMother1_log, 0, 0, checkOverlap);   
    //log.back()->SetVisAttributes(colour->lightblue);              

    unionLH2_1 = new G4UnionSolid("LH2sphere_1+target_H2", tubs.back(), sphere.back(),rotate2,- G4ThreeVector(0.,0.,targetLength - 2*targetRadius));
    //log.push_back(new G4LogicalVolume(unionLH2_1, materials->lh2, "targetCap1_log"));
    //new G4PVPlacement(0, positionVector, log.back(), "targetCap_1", world_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->invisible);

    // adding second semisphere to cylindrical lh2 target
    sphere.push_back(new G4Sphere("LH2Sphere_2", 0., targetRadius, 0., M_PI, 0., M_PI));
    //log.push_back(new G4LogicalVolume(sphere.back(), materials->lh2, "LH2EndCap_log"));
    //new G4PVPlacement(rotate1, positionVector, log.back(), "LH2EndCap", CapMother2_log, 0, 0, checkOverlap);   
    //log.back()->SetVisAttributes(colour->lightblue);

    unionLH2_2 = new G4UnionSolid("LH2Sphere_2+unionLH2_1", unionLH2_1, sphere.back(),rotate1,G4ThreeVector(0.,0.,targetLength - 2*targetRadius));
    log.push_back(new G4LogicalVolume(unionLH2_2, materials->lh2, "targetCap2_log"));
    new G4PVPlacement(0, positionVector, log.back(), "targetCap_2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->lightblue);

    // ___________________________________________ mylar caps ___________________________________________

    // position of mylar caps depends on definition of targetLength
    // mylar cap at the beginning of target
    sphere.push_back(new G4Sphere("BegCapSphere", mylarRadius, mylarRadius + mylarThickness, 0., M_PI, 0., M_PI));    
    log.push_back(new G4LogicalVolume(sphere.back(), materials->mylar, "BegCap_log"));
    new G4PVPlacement(rotate1, positionVector, log.back(), "BegCap", CapMother2_log, 0, 0, checkOverlap);      
    log.back()->SetVisAttributes(colour->white);

    tubs.push_back(new G4Tubs("cylinder_mylar1", mylarRadius, mylarRadius + mylarThickness, mylarCylLength, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "mylarCyl1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,mylarCylLength), log.back(), "myCyl1", CapMother2_log, 0, 0, checkOverlap);      
    log.back()->SetVisAttributes(colour->white);

    //union1 = new G4UnionSolid("semisphere1+cylinder1", sphere.back(), tubs.back(),rotate1,G4ThreeVector(0.,0.,mylarCylLength));                              
    //log.push_back(new G4LogicalVolume(union1, materials->mylar, "mylarU1_log"));
    //new G4PVPlacement(0, positionVector, log.back(), "myU1", CapMother2_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->white);

    // mylar cap at the end of the target
    sphere.push_back(new G4Sphere("EndCapSphere", mylarRadius, mylarRadius + mylarThickness, 0., M_PI, 0., M_PI));    
    log.push_back(new G4LogicalVolume(sphere.back(), materials->mylar, "EndCap_log"));
    new G4PVPlacement(rotate2, positionVector, log.back(), "EndCap", CapMother1_log, 0, 0, checkOverlap);        
    log.back()->SetVisAttributes(colour->white);

    tubs.push_back(new G4Tubs("cylinder_mylar2", mylarRadius, mylarRadius + mylarThickness, mylarCylLength, 0., 2.* M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "mylarCyl2_log"));
    new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,mylarCylLength), log.back(), "myCyl2", CapMother1_log, 0, 0, checkOverlap);     
    log.back()->SetVisAttributes(colour->white);

    //union2 = new G4UnionSolid("semisphere2+cylinder2", sphere.back(), tubs.back());                            // second mylar cap at end of target
    //log.push_back(new G4LogicalVolume(union2, materials->mylar, "mylarU2_log"));
    //new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength - mylarCylLength/2), log.back(), "myU2", CapMother1_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->white);

    // ___________________________________________ kapton cylinder around target __________________________________________

    tubs.push_back(new G4Tubs("kapton_tube", kaptonRadius, kaptonRadius + kaptonThickness, kaptonLength, 0., 2 * M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->kapton, "kapton_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., + steelLength + kaptonLength + mylarCylLength), log.back(), "KaptonTube", CapMother2_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->orange_50);

    // __________________________________________ stainless steel cylinder around target ________________________________

    tubs.push_back(new G4Tubs("steelCylinder", steelRadius, steelRadius + steelThickness, steelLength, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "steelCyl_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., +steelLength), log.back(), "steelCyl", CapMother2_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray);


  // ________________________________________ target insulation case __________________________________________

    // first cap of case at beginning of target
    // cylinder 2 
    G4Tubs* cylinder = new G4Tubs("case_BegCap1", targetRadius + mylarThickness, targetRadius + mylarThickness + caseThickness, caseLength1, 0., 2.*M_PI);
    //tubs.push_back(new G4Tubs("case_BegCap1", targetRadius + mylarThickness + steelThickness, targetRadius + mylarThickness + steelThickness + caseThickness, caseLength2, 0., 2.*M_PI));
    //log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapCyl_log1"));      
    //new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0., targetRadius + mylarThickness - caseLength2), log.back(), "targetCaseCyl_1", CapMother2_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->silver);
    // disc
    tubs.push_back(new G4Tubs("case_BegCap2", 0., targetRadius + mylarThickness + caseThickness, caseThickness, 0., 2.*M_PI));
    //log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapDisc_log1"));      
    //new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0., targetRadius + mylarThickness + caseThickness), log.back(), "targetCaseDisc_1", CapMother2_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->silver);

   // unionCaseCap = new G4UnionSolid("disc+cylinder", tubs.back(), cylinder);
   // log.push_back(new G4LogicalVolume(unionCaseCap, materials->AlMylar, "unionCaseCap_log"));        
   // new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,targetLength + mylarThickness - caseThickness), log.back(), "caseCap", CapMother2_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->green);

   // ************************************* first cap of case at beginning of target ******************************
   // G4Tubs* cylinder = new G4Tubs("case_BegCap1", targetRadius + mylarThickness, targetRadius + mylarThickness + caseThickness, caseLength1, 0., 2.*M_PI);
   // log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapCyl_log1"));      // material to be checked
   // new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0., targetRadius + mylarThickness - caseLength2), log.back(), "targetCaseCyl1", CapMother2_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->yellow);

    // rest of case
    // cylinder 1
    tubs.push_back(new G4Tubs("case_EndCap1", targetRadius + mylarThickness, targetRadius + mylarThickness + caseThickness, caseLength1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseEndCapCyl_log1"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetRadius  + mylarThickness - caseLength1), log.back(), "targetCaseEndCyl_1", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // disc
    tubs.push_back(new G4Tubs("case_EndCap2", 0., targetRadius + mylarThickness + caseThickness, caseThickness, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseEndCapDisc_log2"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetRadius + mylarThickness + caseThickness), log.back(), "targetCaseEndDisc_2", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);

   // G4UnionSolid* unionCaseCap1 = new G4UnionSolid("disc+cylinder", cylinder1, tubs.back());
   // log.push_back(new G4LogicalVolume(unionCaseCap1, materials->AlMylar, "unionCaseCap1_log"));        
   // new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength + mylarThickness + caseThickness), log.back(), "caseCap1", CapMother1_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->silver);


  // ************************************* to be deleted probably ***************************************
    // 1st cone
/*    double posZ_cone1 = targetRadius + mylarThickness - 2*caseLength1 - caseCone;
    cons.push_back(new G4Cons("cone1",targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, targetRadius + kaptonThickness,targetRadius + kaptonThickness + caseThickness1, caseCone, 0., 2*M_PI));
    log.push_back(new G4LogicalVolume(cons.back(), materials->AlMylar, "caseCone1_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone1), log.back(), "targetCaseCone_1", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // thin cylinder under target holder
    tubs.push_back(new G4Tubs("casePiece1_1", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, caseLength1_1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapCyl1_1_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., posZ_cone1 - caseLength1_1 - caseCone), log.back(), "targetCaseCyl1_1", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // 2nd cone
    double posZ_cone2 = posZ_cone1 - 2*caseLength1_1 - 2*caseCone;
    cons.push_back(new G4Cons("cone2", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1, targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, caseCone, 0., 2*M_PI));
    log.push_back(new G4LogicalVolume(cons.back(), materials->AlMylar, "caseCone2_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone2), log.back(), "targetCaseCone_2", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // cylinder 2
    tubs.push_back(new G4Tubs("case_C2", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1, caseLength2, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseC2_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone2 - caseLength2 - caseCone), log.back(), "targetCaseC2", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // 3rd cone
    double posZ_cone3 = posZ_cone2 - 2*caseCone - 2*caseLength2;
    cons.push_back(new G4Cons("cone3",targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, targetRadius + kaptonThickness,targetRadius + kaptonThickness + caseThickness1, caseCone, 0., 2*M_PI));
    log.push_back(new G4LogicalVolume(cons.back(), materials->AlMylar, "caseCone3_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone3), log.back(), "targetCaseCone_3", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // thin cylinder under target holder
    tubs.push_back(new G4Tubs("casePiece2_1", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, caseLength1_1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapCyl2_1_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., posZ_cone3 - caseLength1_1 - caseCone), log.back(), "targetCaseCyl2_1", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // 4th cone
    double posZ_cone4 = posZ_cone3 - 2*caseLength1_1 - 2*caseCone;
    cons.push_back(new G4Cons("cone4", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1, targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, caseCone, 0., 2*M_PI));
    log.push_back(new G4LogicalVolume(cons.back(), materials->AlMylar, "caseCone4_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone4), log.back(), "targetCaseCone_4", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // cylinder 3
    tubs.push_back(new G4Tubs("case_C3", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1, caseLength3, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseC3_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone4 - caseLength3 - caseCone), log.back(), "targetCaseC3", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // 5th cone
    double posZ_cone5 = posZ_cone4 - 2*caseCone - 2*caseLength3;
    cons.push_back(new G4Cons("cone5",targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, targetRadius + kaptonThickness,targetRadius + kaptonThickness + caseThickness1, caseCone, 0., 2*M_PI));
    log.push_back(new G4LogicalVolume(cons.back(), materials->AlMylar, "caseCone5_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone5), log.back(), "targetCaseCone_5", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // thin cylinder under target holder
    tubs.push_back(new G4Tubs("casePiece3_1", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, caseLength1_1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapCyl3_1_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., posZ_cone5 - caseLength1_1 - caseCone), log.back(), "targetCaseCyl3_1", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // 6th cone
    double posZ_cone6 = posZ_cone5 - 2*caseLength1_1 - 2*caseCone;
    cons.push_back(new G4Cons("cone6", targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1, targetRadius + kaptonThickness, targetRadius + kaptonThickness + caseThickness1_1, caseCone, 0., 2*M_PI));
    log.push_back(new G4LogicalVolume(cons.back(), materials->AlMylar, "caseCone6_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone6), log.back(), "targetCaseCone_6", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
    // cylinder 4
    tubs.push_back(new G4Tubs("case_C4", targetRadius + kaptonThickness + steelThickness, targetRadius + kaptonThickness + steelThickness + caseThickness1, caseLength4, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseC4_log"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_cone6 - caseLength4 - caseCone), log.back(), "targetCaseC4", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->silver);
*/
  // __________________________________________ mother vol for target holder1 (TH) ______________________________________

    double posZ1 = posZ_C2j2 - JunctionLength12_2 + dz2 + 2*dz1 + 5*trapThickness;       
    box.push_back(new G4Box("motherTH1_vol", 2*cavityRadius1, 2*cavityRadius1, trapThickness));
    log.push_back(new G4LogicalVolume(box.back(), materials->vacuum_noOptical, "motherTH1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ1), log.back(), "MotherTH1", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    MotherTH1_log = log.back();

  // ____________________________ triangular-ish target holder 1 ____________________________

  double trapDh = TMath::Tan(60 * CLHEP::degree)*trapDx;
  double trapDL = trapDx / TMath::Sin(30 * CLHEP::degree);
  double trapL = trapSide + 2*trapDL;

  // the center of the triangle doesn't coincide with trapHeight/2
  double trap_h1 = trapL / (2*TMath::Cos(30 * CLHEP::degree));       // longest segment
  double trap_h2 = trapHeight - trap_h1;                             // shortest segment
  double trap_l1 = trap_h1 / (2*TMath::Cos(30 * CLHEP::degree));
  double trap_l2 = trapL - trap_l1;

  // vertices coordinates must be given in anti-clockwise order, from bottom to top along z axis.

    double X1 = 0.;
    double Y1 = -trap_h1;
    double X2 = (trap_l1 + trap_l2)* TMath::Sin(30 * CLHEP::degree);
    double Y2 = trap_h2;

  // G4Trapezoid must have 8 vertices so the 1st and 4th coincide
  vector<G4TwoVector> vertices = {{X1,Y1},{-X2,Y2},{X2,Y2},{X1,Y1},{X1,Y1},{-X2,Y2},{X2,Y2},{X1,Y1}};

  trap.push_back(new G4GenericTrap("targetHolder1", trapThickness, vertices));
  
  // ring to remove the edges
  // subtraction volume cannot be the same thickness as the main volume otherwise the result is undefined
  tubs.push_back(new G4Tubs("targetHolderEdge", cavityRadius1, cavityRadius1 + trapDh, 2*trapThickness, 0., 2*M_PI));
  subtractionEdge = new G4SubtractionSolid("targetHolder1-targetHolderEdge",trap.back(),tubs.back(),0,G4ThreeVector(0.,0.,0.));         
  // target holder hole
  tubs.push_back(new G4Tubs("targetHolderHole", 0., kaptonRadius + kaptonThickness + caseThickness, 2*trapThickness, 0., 2*M_PI));
  subtractionHole = new G4SubtractionSolid("subtractionEdge-targetHolderHole",subtractionEdge,tubs.back(),0,G4ThreeVector(0.,0.,0.));         
 
  // first rectangular hole
  double dR = HoleD - 2*(kaptonRadius + kaptonThickness + caseThickness) - yBox;
  G4ThreeVector posBox1(0., kaptonRadius + kaptonThickness + caseThickness + dR,0.);
  box.push_back(new G4Box("box1", xBox, yBox, 2*trapThickness));
  subtractionBox1 = new G4SubtractionSolid("subtractionHole-box1", subtractionHole, box.back(),0,-posBox1);

  // second  and third rectangular hole 
  double posX_box = (kaptonRadius + kaptonThickness + caseThickness + dR) * TMath::Cos(30 * CLHEP::degree);
  double posY_box = (kaptonRadius + kaptonThickness + caseThickness + dR) * TMath::Sin(30 * CLHEP::degree);
  G4ThreeVector posBox2(-posX_box,posY_box,0.);
  G4ThreeVector posBox3(posX_box,posY_box,0.);

  box.push_back(new G4Box("box2", xBox, yBox, 2*trapThickness));
  subtractionBox2 = new G4SubtractionSolid("subtractionBox1-box2", subtractionBox1, box.back(),rotate120,posBox2);
  box.push_back(new G4Box("box3", xBox, yBox, 2*trapThickness));
  subtractionBox3 = new G4SubtractionSolid("subtractionBox2-box3", subtractionBox2, box.back(),rotate60,posBox3);
  log.push_back(new G4LogicalVolume(subtractionBox3, materials->rohacell, "subtractionBox3_log"));
  new G4PVPlacement(0, positionVector, log.back(), "SubtractionBox3", MotherTH1_log, 0, 0, checkOverlap);
  log.back()->SetVisAttributes(colour->magenta);

  // __________________________________________ mother vol for target holder2 (TH) ______________________________________

    double posZ2 = posZ_C2j2 - JunctionLength12_2 + dz2 + dz1 + 3*trapThickness;       
    box.push_back(new G4Box("motherTH2_vol", 2*cavityRadius1, 2*cavityRadius1, trapThickness));
    log.push_back(new G4LogicalVolume(box.back(), materials->vacuum_noOptical, "motherTH2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ), log.back(), "MotherTH2", CapMother1_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    MotherTH2_log = log.back();

  // ____________________________ triangular-ish target holder 2 ____________________________

  trap.push_back(new G4GenericTrap("targetHolder1", trapThickness, vertices));
  
  // ring to remove the edges
  // subtraction volume cannot be the same thickness as the main volume otherwise the result is undefined
  tubs.push_back(new G4Tubs("targetHolderEdge", cavityRadius1, cavityRadius1 + trapDh, 2*trapThickness, 0., 2*M_PI));
  subtractionEdge = new G4SubtractionSolid("targetHolder1-targetHolderEdge",trap.back(),tubs.back(),0,G4ThreeVector(0.,0.,0.));         
  // target holder hole
  tubs.push_back(new G4Tubs("targetHolderHole", 0., kaptonRadius + kaptonThickness + caseThickness, 2*trapThickness, 0., 2*M_PI));
  subtractionHole = new G4SubtractionSolid("subtractionEdge-targetHolderHole",subtractionEdge,tubs.back(),0,G4ThreeVector(0.,0.,0.));         
 
  // first rectangular hole
  double dR = HoleD - 2*(kaptonRadius + kaptonThickness + caseThickness) - yBox;
  G4ThreeVector posBox1(0., kaptonRadius + kaptonThickness + caseThickness + dR,0.);
  box.push_back(new G4Box("box1", xBox, yBox, 2*trapThickness));
  subtractionBox1 = new G4SubtractionSolid("subtractionHole-box1", subtractionHole, box.back(),0,-posBox1);

  // second  and third rectangular hole 
  double posX_box = (kaptonRadius + kaptonThickness + caseThickness + dR) * TMath::Cos(30 * CLHEP::degree);
  double posY_box = (kaptonRadius + kaptonThickness + caseThickness + dR) * TMath::Sin(30 * CLHEP::degree);
  G4ThreeVector posBox2(-posX_box,posY_box,0.);
  G4ThreeVector posBox3(posX_box,posY_box,0.);

  box.push_back(new G4Box("box2", xBox, yBox, 2*trapThickness));
  subtractionBox2 = new G4SubtractionSolid("subtractionBox1-box2", subtractionBox1, box.back(),rotate120,posBox2);
  box.push_back(new G4Box("box3", xBox, yBox, 2*trapThickness));
  subtractionBox3 = new G4SubtractionSolid("subtractionBox2-box3", subtractionBox2, box.back(),rotate60,posBox3);
  log.push_back(new G4LogicalVolume(subtractionBox3, materials->rohacell, "subtractionBox3_log"));
  new G4PVPlacement(0, positionVector, log.back(), "SubtractionBox3", MotherTH1_log, 0, 0, checkOverlap);
  log.back()->SetVisAttributes(colour->magenta);

}



void T4LH2Target2024::getTargetDetDat(std::vector<T4STargetInformation> &targetInformation) {
  T4STargetInformation target; 
  target.name = "targetH2_2024";
  target.rotMatrix = TGEANT::ROT_XtoZ;
  target.sh = 5;
  target.xSize = targetRadius;
  target.ySize = targetRadius;
  target.zSize = targetLength;
  target.xCen = positionVector.getX();
  target.yCen = positionVector.getY();
  target.zCen = positionVector.getZ();

  targetInformation.push_back(target);
}

T4LH2Target2024::~T4LH2Target2024(void)
{
 for (unsigned int i = 0; i < box.size();i++)
    delete box.at(i);
  box.clear();

  for (unsigned int i = 0; i < tubs.size();i++)
    delete tubs.at(i);
  tubs.clear();

  for (unsigned int i = 0; i < log.size();i++)
    delete log.at(i);
  log.clear();

  for (unsigned int i = 0; i < sphere.size();i++)
    delete sphere.at(i);
  sphere.clear();

  for (unsigned int i = 0; i < trap.size();i++)
    delete trap.at(i);
  trap.clear();

}