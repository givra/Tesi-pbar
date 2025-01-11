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

  JunctionLength12_1 = 16. / 2 * CLHEP::mm;
  JunctionRadiusMin12_1 = 168. / 2 * CLHEP::mm;
  JunctionRadiusMax12_1 = 370. / 2 * CLHEP::mm;

  JunctionLength12_O = 5.6 / 2 * CLHEP::mm;   
  JunctionRadiusMin12_O = 313.6 / 2 * CLHEP::mm;
  JunctionRadiusMax12_O = 358 / 2 * CLHEP::mm;

  JunctionLength12_2 = 16. / 2 * CLHEP::mm;
  JunctionRadiusMin12_2 = 318 / 2 * CLHEP::mm;
  JunctionRadiusMax12_2 = 370 / 2 * CLHEP::mm;

  JunctionLength23_1 =  15.0 / 2 * CLHEP::mm;
  JunctionRadiusMin23_1 =  110 / 2 * CLHEP::mm;
  JunctionRadiusMax23_1 =  239.5 / 2 * CLHEP::mm;

  JunctionLength23_2 = 14.0 / 2 * CLHEP::mm;      // junction O-ring parameters between C2 and C3
  JunctionRadiusMin23_2 = 164.3 / 2 * CLHEP::mm;
  JunctionRadiusMax23_2 = 239 / 2 * CLHEP::mm;

  shift = 103. * CLHEP::mm;                              
  dz1 = 480. * CLHEP::mm; 
  dz2 = 204. * CLHEP::mm; 

  targetRadius =  39.8 / 2 * CLHEP::mm;                   
  targetLength = 1400 / 2 * CLHEP::mm;                   // already contains in its definition the 2 radii
  mylarRadius =  40 / 2 * CLHEP::mm;                     // for mylar semi-spheres end caps
  mylarThickness =  0.125 * CLHEP::mm;                   // 125 um
  mylarCylLength = 6.2 / 2 * CLHEP::mm;

  trapThickness = 20. / 2 * CLHEP::mm;                     // target holder thickness
  trapDx = 12.4 /2 * CLHEP::mm;
  trapSide = 269.6 /2 * CLHEP::mm;
  HoleD = 51.5 * CLHEP::mm;
  yBox =  14.5 / 2 * CLHEP::mm;
  xBox = 30.4 / 2 * CLHEP::mm;

  trapDx_3 = 14.9 / 2 * CLHEP::mm;    
  trapSide_3 = 134.8 / 2 * CLHEP::mm;  

  mylarWindowRadius = 370. / 2 * CLHEP::mm;
  mylarWindowThick = 0.35 * CLHEP::mm;
  mylarWindowRadius_beg = 139.5 / 2 * CLHEP::mm;
  mylarWindowThick_beg = 0.125 * CLHEP::mm;
  
  kaptonLength =  1265 / 2 * CLHEP::mm;
  kaptonRadius =  40 / 2 * CLHEP::mm;
  kaptonThickness = 0.125 * CLHEP::mm;

  steelRadius = 39.8 / 2 * CLHEP::mm;         
  steelThickness = 0.4 * CLHEP::mm;
  steelLength = 95 / 2 * CLHEP::mm;

  caseThickness = 30 * 11 * CLHEP::um;
  caseLength1 = 1347.5 / 2 * CLHEP::mm;
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

    double posZ_C1 = targetLength + shift - cavityLength1; 
    tubs.push_back(new G4Tubs("C1", cavityRadius1, cavityRadius1 + cavityThickness1, cavityLength1, 0., 2.*M_PI));
  //  // logical volume
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0.,posZ_C1), log.back(), "Cavity1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    // junctions between C1 and C2 
    double posZ_C1j2 = posZ_C1 - cavityLength1;
    tubs.push_back(new G4Tubs("C1j2", JunctionRadiusMin12_2, JunctionRadiusMax12_2, JunctionLength12_2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C1j2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C1j2), log.back(), "Cavity1j2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    double posZ_O12 = posZ_C1j2 - JunctionLength12_2 - JunctionLength12_O;
    tubs.push_back(new G4Tubs("C1_O", JunctionRadiusMin12_O, JunctionRadiusMax12_O, JunctionLength12_O, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "C1_O_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_O12), log.back(), "Cavity1_O", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->gray);
    
    double posZ_C1j1 = posZ_O12 - JunctionLength12_1 - JunctionLength12_O;
    tubs.push_back(new G4Tubs("C1j1", JunctionRadiusMin12_1, JunctionRadiusMax12_1, JunctionLength12_1, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C1j1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C1j1), log.back(), "Cavity1j1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

  // _____________________________________ small stainless-steel cylinder #2 _____________________________________

    double posZ_C2 = posZ_C1j1 - cavityLength2;
    tubs.push_back(new G4Tubs("C2", cavityRadius2, cavityRadius2 + cavityThickness2, cavityLength2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C2), log.back(), "Cavity2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    // junctions between #2 and #3
    double posZ_C2j2 = posZ_C2 - cavityLength2 - JunctionLength23_2;
    tubs.push_back(new G4Tubs("C2j2", JunctionRadiusMin23_2, JunctionRadiusMax23_2, JunctionLength23_2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C2j2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C2j2 ), log.back(), "Cavity2j2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    double posZ_C2j1 = posZ_C2j2 - JunctionLength23_2 - JunctionLength23_1;
    tubs.push_back(new G4Tubs("C2j1", JunctionRadiusMin23_1, JunctionRadiusMax23_1, JunctionLength23_1, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C2j1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C2j1), log.back(), "Cavity2j1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

  // _____________________________________ #3 stainless-steel cylinder _____________________________________

    double posZ_C3 = posZ_C2j1 - JunctionLength23_1 - cavityLength3;
    tubs.push_back(new G4Tubs("C3", cavityRadius3, cavityRadius3 + cavityThickness3, cavityLength3, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C3_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_C3), log.back(), "Cavity3", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

  // _______________________________________ mother log volumes for cylinder 4 _______________________________________
   
    double posZ_motherC4 = posZ_C3 - cavityLength3 - cavityLength4;
/*    tubs.push_back(new G4Tubs("motherC4_vol", 0., 2*cavityRadius4, 1.5*cavityLength4, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "motherC4_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_motherC4), log.back(), "MotherC4", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->green);
    MotherC4_log = log.back();
*/
  // _____________________________________ #4 stainless-steel cylinder _____________________________________
    
    tubs.push_back(new G4Tubs("C4", cavityRadius4, cavityRadius4 + cavityThickness4, cavityLength4, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_motherC4), log.back(), "Cavity4", world_log, 0, 0, checkOverlap);    // MotherC4_log
    log.back()->SetVisAttributes(colour->darkgray_50);      

    // disc to create a hole in Cylinder 4 along y axis
    //G4Tubs* disc4 = new G4Tubs("disc4", 0., cavityRadius4_1 + cavityThickness4_1, cavityLength4_1, 0., 2.*M_PI);
    //G4ThreeVector Pos_C4_1(0.,cavityRadius4,0.);

    //subtractionC4 = new G4SubtractionSolid("C4-disc4", tubs.back(), disc4, 0, Pos_C4_1);
    

    // tubs.push_back(new G4Tubs("C4_1", cavityRadius4_1, cavityRadius4_1 + cavityThickness4_1, cavityLength4_1, 0., 2.*M_PI));      
    // log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4_1_log"));
    // new G4PVPlacement(rotate1, positionVector + G4ThreeVector(0.,cavityRadius4,0.), log.back(), "Cavity4_1", MotherC4_log, 0, 0, checkOverlap);
    // log.back()->SetVisAttributes(colour->darkgray_50);

   // unionC4 = new G4UnionSolid("subtractionC4+C4_1",subtractionC4,tubs.back(),rotate1,G4ThreeVector(0.,cavityRadius4,0.));                              
   // log.push_back(new G4LogicalVolume(union1, materials->stainlessSteel, "UnionC4_log"));
   // new G4PVPlacement(0, positionVector, log.back(), "UnionC4", MotherC4_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->darkgray_50);

    // O-rings at the beginning of C4

    tubs.push_back(new G4Tubs("C4_side3", cavityRadius4, SideRadiusMax2, SideLength2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4Side3_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_motherC4 -cavityLength4 - SideLength2), log.back(), "Cavity4Side3", world_log, 0, 0, checkOverlap);    // MotherC4_log   
    log.back()->SetVisAttributes(colour->darkgray_50);

    tubs.push_back(new G4Tubs("C4_side2", SideRadiusMin2, SideRadiusMax2, SideLength2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4Side2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_motherC4 -cavityLength4 - 3*SideLength2), log.back(), "Cavity4Side2", world_log, 0, 0, checkOverlap);     // MotherC4_log
    log.back()->SetVisAttributes(colour->darkgray_50);
    
    tubs.push_back(new G4Tubs("C4_MW_beg", 0., mylarWindowRadius_beg, mylarWindowThick_beg, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "C4_MW_beg_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_motherC4 -cavityLength4 - 3*SideLength2 - mylarWindowThick_beg), log.back(), "Cavity4_MW_beg", world_log, 0, 0, checkOverlap);      // MotherC4_log
    log.back()->SetVisAttributes(colour->blue);

    tubs.push_back(new G4Tubs("C4_side1", SideRadiusMin1, SideRadiusMax1, SideLength1, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4Side1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_motherC4 -cavityLength4 - 4*SideLength2 -2* mylarWindowThick_beg - SideLength1), log.back(), "Cavity4Side1", world_log, 0, 0, checkOverlap);    // MotherC4_log
    log.back()->SetVisAttributes(colour->darkgray_50);

  // ________________________________________ mylar window _________________________________________

    double posZ_MW = targetLength + shift + EndLength1 + 2*EndLength_O + mylarWindowThick;
    tubs.push_back(new G4Tubs("mylarWin", 0., mylarWindowRadius, mylarWindowThick, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "MW_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MW), log.back(), "MylarWindow", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->blue);

    // O-rings for mylar window
    double posZ_MWj1 = targetLength + shift;
    tubs.push_back(new G4Tubs("mylarWin_j1", EndRadiusMin1, EndRadiusMax1, EndLength1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "MWj1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MWj1), log.back(), "MylarWindow_j1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    double posZ_MW_O = posZ_MWj1 + EndLength1 + EndLength_O;
    tubs.push_back(new G4Tubs("mylarWin_O", EndRadiusMin_O, EndRadiusMax_O, EndLength_O, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->Al2O3, "MW_O_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MW_O), log.back(), "MylarWindow_O", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->gray);

    double posZ_MWj2 = posZ_MW + mylarWindowThick + EndLength1;
    tubs.push_back(new G4Tubs("mylarWin_j2", EndRadiusMin1, EndRadiusMax1, EndLength1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "MWj2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_MWj2), log.back(), "MylarWindow_j2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

  // _______________________________________ mother log volumes for target and mylar caps _______________________________________
    
    // #1 referrs to the end of target, closest to mylar window
    // dimensions of mother volume: ensure they're big enough to contain all the solids but not too big to overlap with other mother volumes
    double posZ_mother1 = targetLength - targetRadius;
    tubs.push_back(new G4Tubs("mother1_vol", 0., targetRadius + mylarThickness, 2*targetRadius, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "mother1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ_mother1), log.back(), "MotherVolLog_1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    CapMother1_log = log.back();

    tubs.push_back(new G4Tubs("mother2_vol", 0., targetRadius + mylarThickness, 2*targetRadius, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "mother2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,- posZ_mother1), log.back(), "MotherVolLog_2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    CapMother2_log = log.back();

  // ___________________________________________ H2 target ___________________________________________

   string YearSetup = T4EventManager::getInstance()->getTriggerPlugin()->getPluginName().substr(0,6);
   G4Material* tarMaterial = materials->lh2;
   if( YearSetup == "2024-D" ) tarMaterial = materials->ld2;      

   // cylindrical part
    tubs.push_back(new G4Tubs("target_H2", 0., targetRadius, targetLength - targetRadius, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), tarMaterial, "targetH2_log"));
    addToTarget(new G4PVPlacement(0, positionVector, log.back(), "LH2Hadron", world_log, 0, 0, checkOverlap));
    log.back()->SetVisAttributes(colour->lightblue);
    if( YearSetup == "2024-D" )  log.back()->SetVisAttributes(colour->green);

    // adding first semisphere to cylindrical lh2 target     
    sphere.push_back(new G4Sphere("LH2Sphere_1", 0., targetRadius, 0., M_PI, 0., M_PI));
    log.push_back(new G4LogicalVolume(sphere.back(), tarMaterial, "LH2BegCap_log"));
    addToTarget(new G4PVPlacement(rotate1, positionVector, log.back(), "LH2BegCap", CapMother2_log, 0, 0, checkOverlap));
    //new G4PVPlacement(rotate1, positionVector , log.back(), "LH2BegCap", CapMother2_log, 0, 0, checkOverlap);   
    log.back()->SetVisAttributes(colour->lightblue);              

    //unionLH2_1 = new G4UnionSolid("LH2sphere_1+target_H2", tubs.back(), sphere.back(),rotate2,- G4ThreeVector(0.,0.,targetLength - 2*targetRadius));
    //log.push_back(new G4LogicalVolume(unionLH2_1, materials->lh2, "targetCap1_log"));
    //new G4PVPlacement(0, positionVector, log.back(), "targetCap_1", world_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->invisible);

    // adding second semisphere to cylindrical lh2 target
    sphere.push_back(new G4Sphere("LH2Sphere_2", 0., targetRadius, 0., M_PI, 0., M_PI));
    log.push_back(new G4LogicalVolume(sphere.back(), tarMaterial, "LH2EndCap_log"));
    addToTarget(new G4PVPlacement(rotate2, positionVector, log.back(), "LH2EndCap", CapMother1_log, 0, 0, checkOverlap));
    //new G4PVPlacement(rotate2, positionVector, log.back(), "LH2EndCap", CapMother1_log, 0, 0, checkOverlap);   
    log.back()->SetVisAttributes(colour->lightblue);

    //unionLH2_2 = new G4UnionSolid("LH2Sphere_2+unionLH2_1", unionLH2_1, sphere.back(),rotate1,G4ThreeVector(0.,0.,targetLength - 2*targetRadius));
    //log.push_back(new G4LogicalVolume(unionLH2_2, materials->lh2, "targetCap2_log"));
    //new G4PVPlacement(0, positionVector, log.back(), "targetCap_2", world_log, 0, 0, checkOverlap);
    //log.back()->SetVisAttributes(colour->lightblue);

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

    double posZ_kapton = - posZ_mother1 + 2*steelLength + kaptonLength;
    tubs.push_back(new G4Tubs("kapton_tube", kaptonRadius, kaptonRadius + kaptonThickness, kaptonLength, 0., 2 * M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->kapton, "kapton_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., posZ_kapton), log.back(), "KaptonTube", world_log, 0, 0, checkOverlap);     // CapMother2_log
    log.back()->SetVisAttributes(colour->orange_50);

    // __________________________________________ stainless steel cylinder around target ________________________________
    
    double posZ_steel = - posZ_mother1 + steelLength;
    tubs.push_back(new G4Tubs("steelCylinder", steelRadius, steelRadius + steelThickness, steelLength, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "steelCyl_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., posZ_steel), log.back(), "steelCyl", world_log, 0, 0, checkOverlap);     // CapMother2_log
    log.back()->SetVisAttributes(colour->darkgray);


  // ________________________________________ target insulation case __________________________________________

    // first cap of case at beginning of target
    // cylinder 2 
    //  G4Tubs* cylinder = new G4Tubs("case_BegCap1", targetRadius + mylarThickness, targetRadius + mylarThickness + caseThickness, caseLength1, 0., 2.*M_PI);
    tubs.push_back(new G4Tubs("case_BegCap1", targetRadius + mylarThickness + steelThickness, targetRadius + mylarThickness + steelThickness + caseThickness, caseLength2, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapCyl_log1"));      
    new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0., posZ_mother1 + targetRadius + mylarThickness - caseLength2), log.back(), "targetCaseCyl_1", world_log, 0, 0, checkOverlap);     // CapMother2_log
    log.back()->SetVisAttributes(colour->silver);
    // disc
    tubs.push_back(new G4Tubs("case_BegCap2", 0., targetRadius + mylarThickness + caseThickness, caseThickness, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseCapDisc_log1"));      
    new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0., posZ_mother1 + targetRadius + mylarThickness + caseThickness), log.back(), "targetCaseDisc_1", world_log, 0, 0, checkOverlap);       // CapMother2_log
    log.back()->SetVisAttributes(colour->silver);

   // unionCaseCap = new G4UnionSolid("disc+cylinder", tubs.back(), cylinder);
   // log.push_back(new G4LogicalVolume(unionCaseCap, materials->AlMylar, "unionCaseCap_log"));        
   // new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,targetLength + mylarThickness - caseThickness), log.back(), "caseCap", CapMother2_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->green);

    // rest of case
    // cylinder 1
    tubs.push_back(new G4Tubs("case_EndCap1", targetRadius + mylarThickness, targetRadius + mylarThickness + caseThickness, caseLength1, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseEndCapCyl_log1"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength + mylarThickness - caseLength1), log.back(), "targetCaseEndCyl_1", world_log, 0, 0, checkOverlap);      // CapMother1_log
    log.back()->SetVisAttributes(colour->silver);
    // disc
    tubs.push_back(new G4Tubs("case_EndCap2", 0., targetRadius + mylarThickness + caseThickness, caseThickness, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->AlMylar, "caseEndCapDisc_log2"));      
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength + mylarThickness + caseThickness), log.back(), "targetCaseEndDisc_2", world_log, 0, 0, checkOverlap);   // CapMother1_log
    log.back()->SetVisAttributes(colour->silver);

   // G4UnionSolid* unionCaseCap1 = new G4UnionSolid("disc+cylinder", cylinder1, tubs.back());
   // log.push_back(new G4LogicalVolume(unionCaseCap1, materials->AlMylar, "unionCaseCap1_log"));        
   // new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength + mylarThickness + caseThickness), log.back(), "caseCap1", CapMother1_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->silver);

  // ____________________________ triangular-ish target holder 1 ____________________________

  double trapDh = TMath::Tan(60 * CLHEP::degree)*trapDx;
  double trapDL = trapDx / TMath::Sin(60 * CLHEP::degree);
  double trapL = trapSide + trapDL;
  double trapHeight = 2*trapL*TMath::Cos(30 * CLHEP::degree);

  // the center of the triangle doesn't coincide with trapHeight/2
  double trap_h1 = trapL / TMath::Cos(30 * CLHEP::degree);           // longest segment
  double trap_h2 = trapHeight - trap_h1;                             // shortest segment
  double trap_l1 = trap_h1 / TMath::Cos(30 * CLHEP::degree);
  double trap_l2 = 2*trapL - trap_l1;

  // __________________________________________ mother vol for target holder1 (TH) ______________________________________

    double posZ1 = posZ_C2j2 - JunctionLength23_2 + dz2 + 2*dz1 + 5*trapThickness;       
    tubs.push_back(new G4Tubs("motherTH1_vol", 0, trap_h1, trapThickness, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "motherTH1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ1), log.back(), "MotherTH1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    MotherTH1_log = log.back();


  trap.push_back(new G4Trap("targetHolder1", trapThickness,    // half lenght along z
                          0.0,                                 // theta polar of the line joining the centres of the bases at -/+pDz
                          0.0,                                 // phi azimuthal angle of the line joining the centre of the base at -pDz to the centre of the base at +pDz
                          (trapHeight - trapDh)/2,             // dy1 (half y length at -z)
                          trapDx,                              // dx1 (half x length at -z, smaller y)
                          trapL,                               // dx2 (half x length at -z, bigger y)
                          0.0,                                 // alpha1 (angle along x at -z)
                          (trapHeight - trapDh)/2,             // dy2 (half y length at +z)
                          trapDx,                              // dx3 (half x length at +z, smaller y)
                          trapL,                               // dx4 (half x length at +z, bigger y)
                          0.0));

  
  double shift_TH = (trapHeight - trapDh)/2 - trap_h2;

  // ring to remove the edges
  // subtraction volume cannot be the same thickness as the main volume otherwise the result is undefined
  tubs.push_back(new G4Tubs("targetHolderEdge1", cavityRadius1, cavityRadius1 + 2*trapDh, 1.5*trapThickness, 0., 2*M_PI));
  subtractionEdge1 = new G4SubtractionSolid("targetHolder1-targetHolderEdge1",trap.back(),tubs.back(),0,G4ThreeVector(0.,shift_TH,0.));
  // target holder hole
  tubs.push_back(new G4Tubs("targetHolderHole1", 0., kaptonRadius + kaptonThickness + caseThickness, 1.75*trapThickness, 0., 2*M_PI));
  subtractionHole1 = new G4SubtractionSolid("subtractionEdge1-targetHolderHole1",subtractionEdge1,tubs.back(),0,G4ThreeVector(0.,shift_TH,0.));         
 
  // first rectangular hole
  G4ThreeVector posBox1(0., -shift_TH + kaptonRadius + kaptonThickness + caseThickness,0.);
  box.push_back(new G4Box("box1_1", xBox, yBox, 2*trapThickness));
  subtractionBox1_1 = new G4SubtractionSolid("subtractionHole1-box1_1", subtractionHole1, box.back(),0,-posBox1);

  // second  and third rectangular hole 
  double posX_box = (kaptonRadius + kaptonThickness + caseThickness) * TMath::Cos(30 * CLHEP::degree);
  double posY_box = (kaptonRadius + kaptonThickness + caseThickness) * TMath::Sin(30 * CLHEP::degree);
  G4ThreeVector posBox2(-posX_box,posY_box + shift_TH,0.);
  G4ThreeVector posBox3(posX_box,posY_box + shift_TH,0.);

  box.push_back(new G4Box("box2_1", xBox, yBox, 2.25*trapThickness));
  subtractionBox2_1 = new G4SubtractionSolid("subtractionBox1_1-box2_1", subtractionBox1_1, box.back(),rotate120,posBox2);
  box.push_back(new G4Box("box3_1", xBox, yBox, 2.5*trapThickness));
  subtractionBox3_1 = new G4SubtractionSolid("subtractionBox2_1-box3_1", subtractionBox2_1, box.back(),rotate60,posBox3);
  log.push_back(new G4LogicalVolume(subtractionBox3_1, materials->rohacell, "subtractionBox3_1_log"));
  new G4PVPlacement(0, positionVector - G4ThreeVector(0.,shift_TH,0.), log.back(), "SubtractionBox3_1", MotherTH1_log, 0, 0, checkOverlap);   
  log.back()->SetVisAttributes(colour->white);

  // __________________________________________ mother vol for target holder2 (TH) ______________________________________

    double posZ2 = posZ_C2j2 - JunctionLength23_2 + dz2 + dz1 + 3*trapThickness;       
    tubs.push_back(new G4Tubs("motherTH2_vol", 0, trap_h1, trapThickness, 0., 2.*M_PI));     // cavityRadius1 + 2*trapDh
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "motherTH2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ2), log.back(), "MotherTH2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    MotherTH2_log = log.back();

  // ____________________________ triangular-ish target holder 2 ____________________________

  trap.push_back(new G4Trap("targetHolder2", trapThickness,    // half lenght along z
                          0.0,                                 // theta (polar angle)
                          0.0,                                 // phi (azimuthal angle)
                          (trapHeight - trapDh)/2,             // dy1 (half y length at -z)
                          trapDx,                              // dx1 (half x length at -z, smaller y)
                          trapL,                               // dx2 (half x length at -z, bigger y)
                          0.0,                                 // alpha1 (angle along x at -z)
                          (trapHeight - trapDh)/2,             // dy2 (half y length at +z)
                          trapDx,                              // dx3 (half x length at +z, smaller y)
                          trapL,                               // dx4 (half x length at +z, bigger y)
                          0.0));
  
  // ring to remove the edges
  // subtraction volume cannot be the same thickness as the main volume otherwise the result is undefined
  tubs.push_back(new G4Tubs("targetHolderEdge2", cavityRadius1, cavityRadius1 + 2*trapDh, 1.5*trapThickness, 0., 2*M_PI));
  subtractionEdge2 = new G4SubtractionSolid("targetHolder2-targetHolderEdge2",trap.back(),tubs.back(),0,G4ThreeVector(0.,shift_TH,0.));         
  // target holder hole
  tubs.push_back(new G4Tubs("targetHolderHole2", 0., kaptonRadius + kaptonThickness + caseThickness, 1.75*trapThickness, 0., 2*M_PI));
  subtractionHole2 = new G4SubtractionSolid("subtractionEdge2-targetHolderHole2",subtractionEdge2,tubs.back(),0,G4ThreeVector(0.,shift_TH,0.));

  // rectangular holes
  box.push_back(new G4Box("box1_2", xBox, yBox, 2*trapThickness));
  subtractionBox1_2 = new G4SubtractionSolid("subtractionHole2_2-box1_2", subtractionHole2, box.back(),0,-posBox1);
  box.push_back(new G4Box("box2_2", xBox, yBox, 2.25*trapThickness));
  subtractionBox2_2 = new G4SubtractionSolid("subtractionBox1_2-box2_2", subtractionBox1_2, box.back(),rotate120,posBox2);
  box.push_back(new G4Box("box3_2", xBox, yBox, 2.5*trapThickness));
  subtractionBox3_2 = new G4SubtractionSolid("subtractionBox2_2-box3_2", subtractionBox2_2, box.back(),rotate60,posBox3);
  log.push_back(new G4LogicalVolume(subtractionBox3_2, materials->rohacell, "subtractionBox3_2_log"));
  new G4PVPlacement(0, positionVector - G4ThreeVector(0.,shift_TH,0.), log.back(), "SubtractionBox3_2", MotherTH2_log, 0, 0, checkOverlap);
  log.back()->SetVisAttributes(colour->white);

  // ____________________________ triangular-ish target holder 3 ____________________________

  double trapDh_3 = TMath::Tan(60 * CLHEP::degree)*trapDx_3;
  double trapDL_3 = trapDx_3 / TMath::Sin(30 * CLHEP::degree);
  double trapL_3 = trapSide_3 + trapDL_3;
  double trapHeight_3 = 2*trapL_3*TMath::Cos(30 * CLHEP::degree);

  // the center of the triangle doesn't coincide with trapHeight/2
  double trap_h1_3 = trapL_3 / TMath::Cos(30 * CLHEP::degree);         // longest segment
  double trap_h2_3 = trapHeight_3 - trap_h1_3;                             // shortest segment
  double trap_l1_3 = trap_h1_3 / TMath::Cos(30 * CLHEP::degree);
  double trap_l2_3 = 2*trapL_3 - trap_l1_3;

// __________________________________________ mother vol for target holder3 (TH) ______________________________________

    double posZ3 = posZ_C2j2 - JunctionLength23_2 + dz2 + trapThickness;       
    tubs.push_back(new G4Tubs("motherTH3_vol", 0, trap_h1_3, trapThickness, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "motherTH3_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,posZ3), log.back(), "MotherTH3", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    MotherTH3_log = log.back();

trap.push_back(new G4Trap("targetHolder3", trapThickness,    // half lenght along z
                          0.0,                                 // theta (polar angle)
                          0.0,                                 // phi (azimuthal angle)
                          (trapHeight_3 - trapDh_3)/2,             // dy1 (half y length at -z)
                          trapDx_3,                              // dx1 (half x length at -z, smaller y)
                          trapL_3,                               // dx2 (half x length at -z, bigger y)
                          0.0,                                 // alpha1 (angle along x at -z)
                          (trapHeight_3 - trapDh_3)/2,             // dy2 (half y length at +z)
                          trapDx_3,                              // dx3 (half x length at +z, smaller y)
                          trapL_3,                               // dx4 (half x length at +z, bigger y)
                          0.0));

  double shift_TH_3 = (trapHeight_3 - trapDh_3)/2 - trap_h2_3;
  G4ThreeVector posBox1_3(0., -shift_TH_3 + kaptonRadius + kaptonThickness + caseThickness,0.);
  G4ThreeVector posBox2_3(-posX_box,posY_box + shift_TH_3,0.);
  G4ThreeVector posBox3_3(posX_box,posY_box + shift_TH_3,0.);
  
  // ring to remove the edges
  // subtraction volume cannot be the same thickness as the main volume otherwise the result is undefined
  tubs.push_back(new G4Tubs("targetHolderEdge3", cavityRadius2, trap_h1_3, 1.5*trapThickness, 0., 2*M_PI));
  subtractionEdge3 = new G4SubtractionSolid("targetHolder3-targetHolderEdge3",trap.back(),tubs.back(),0,G4ThreeVector(0.,shift_TH_3,0.));         
  // target holder hole
  tubs.push_back(new G4Tubs("targetHolderHole3", 0., kaptonRadius + kaptonThickness + caseThickness, 1.75*trapThickness, 0., 2*M_PI));
  subtractionHole3 = new G4SubtractionSolid("subtractionEdge3-targetHolderHole3",subtractionEdge3,tubs.back(),0,G4ThreeVector(0.,shift_TH_3,0.));         
  box.push_back(new G4Box("box1_3", xBox, yBox, 2*trapThickness));
  subtractionBox1_3 = new G4SubtractionSolid("subtractionHole3-box1_3", subtractionHole3, box.back(),0,-posBox1_3);
  box.push_back(new G4Box("box2_3", xBox, yBox, 2.25*trapThickness));
  subtractionBox2_3 = new G4SubtractionSolid("subtractionBox1_3-box2_3", subtractionBox1_3, box.back(),rotate120,posBox2_3);
  box.push_back(new G4Box("box3_3", xBox, yBox, 2.5*trapThickness));
  subtractionBox3_3 = new G4SubtractionSolid("subtractionBox2_3-box3_3", subtractionBox2_3, box.back(),rotate60,posBox3_3);
  log.push_back(new G4LogicalVolume(subtractionBox3_3, materials->rohacell, "subtractionBox3_3_log"));
  new G4PVPlacement(0, positionVector - G4ThreeVector(0.,shift_TH_3,0.), log.back(), "SubtractionBox3_3", MotherTH3_log, 0, 0, checkOverlap);
  log.back()->SetVisAttributes(colour->white);

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