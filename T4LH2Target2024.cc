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

  shift = 103. * CLHEP::mm;                              // shift along target axis of outer cavity center wrt target center

  targetRadius =  39.8 / 2 * CLHEP::mm;                    // o 40mm????? lh2 volume
  targetLength = 1400 / 2 * CLHEP::mm;                   // already contains in its definition the 2 radii
  mylarRadius =  40 / 2 * CLHEP::mm;                     // for mylar semi-spheres end caps
  mylarThickness =  0.125 * CLHEP::mm;                   // 125 um
  mylarCylLength = 6.2 / 2 * CLHEP::mm;

  trapThickness = 20.1 / 2 * CLHEP::mm;                     // target holder thickness
  trapDx = 12.4 * CLHEP::mm;
  trapSide = 269.6 / 2 * CLHEP::mm;

  mylarWindowRadius = 370. / 2 * CLHEP::mm;
  mylarWindowThick = 0.35 * CLHEP::mm;
  xBox =  6.1 * CLHEP::mm;                          // half length
  yBox = 15.2 * CLHEP::mm;

  kaptonLength =  1265 / 2 * CLHEP::mm;
  kaptonRadius =  40 / 2 * CLHEP::mm;
  kaptonThickness = 0.125 * CLHEP::mm;

  steelRadius = 39.8 / 2 * CLHEP::mm;         // inner or outer diameter?
  steelThickness = 0.4 * CLHEP::mm;
  steelLength = 95 / 2 * CLHEP::mm;

  caseThickness1 = 8.5 * CLHEP::mm;
  caseThickness1_1 = 3.5 * CLHEP::mm;
  caseThickness2 = 7.5 / 2 * CLHEP::mm;

  caseLength1_1 = 16. / 2 * CLHEP::mm;
  caseCone = 7 / 2 * CLHEP::mm;
  caseLength1 = 123.5 / 2 * CLHEP::mm;
  caseLength2 = 467 / 2 * CLHEP::mm;
  caseLength3 = 467 / 2 * CLHEP::mm;
  caseLength4 = 194.5 * CLHEP::mm;
  caseLength5 = 25 / 2 * CLHEP::mm;

  rotate1 = new G4RotationMatrix();
  rotate1->rotateX(90 * CLHEP::deg);
  rotate2 = new G4RotationMatrix();
  rotate2->rotateX(-90 * CLHEP::deg);
 

  setDimension(positionVector.z() - targetLength, positionVector.z() + targetLength);
  setELossParams(0., 0.);
}

void T4LH2Target2024::construct(G4LogicalVolume* world_log)
{

  // _____________________________________ biggest stainless-steel cylinder _____________________________________
    tubs.push_back(new G4Tubs("C1", cavityRadius1, cavityRadius1 + cavityThickness1, cavityLength1, 0., 2.*M_PI));
    // logical volume
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., cavityLength2 + shift), log.back(), "Cavity1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    // ***************** corona circolare tra 1 e 2 ***************************
    
  // _____________________________________ small stainless-steel cylinder _____________________________________
    tubs.push_back(new G4Tubs("C2", cavityRadius2, cavityRadius2 + cavityThickness2, cavityLength2, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C2_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., -cavityLength1 + shift), log.back(), "Cavity2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

  // _____________________________________ #3 stainless-steel cylinder _____________________________________
    tubs.push_back(new G4Tubs("C3", cavityRadius3, cavityRadius3 + cavityThickness3, cavityLength3, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C3_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., -cavityLength1 + shift - cavityLength2 - cavityLength3), log.back(), "Cavity3", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    // ***************** mother volume just for cylinder 4??? ************************

  // _____________________________________ #4 stainless-steel cylinder _____________________________________
    
    G4Tubs* disc4 = new G4Tubs("disc4", 0., cavityRadius4_1 + cavityThickness4_1, cavityLength4_1, 0., 2.*M_PI);
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "disc4_log"));
    new G4PVPlacement(rotate1, positionVector + G4ThreeVector(0.,cavityRadius4, -cavityLength1 + shift - cavityLength2 - 2*cavityLength3 - cavityLength4), log.back(), "Disc4", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);

    tubs.push_back(new G4Tubs("C4", cavityRadius4, cavityRadius4 + cavityThickness4, cavityLength4, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., -cavityLength1 + shift - cavityLength2 - 2*cavityLength3 - cavityLength4), log.back(), "Cavity4", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    subtraction2 = new G4SubtractionSolid("cylinder4 - disc4", tubs.back(), disc4);
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "intersection1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., -cavityLength1 + shift - cavityLength2 - 2*cavityLength3 - cavityLength4), log.back(), "intersection1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    tubs.push_back(new G4Tubs("C4_1", cavityRadius4_1, cavityRadius4_1 + cavityThickness4_1, cavityLength4_1, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4_1_log"));
    new G4PVPlacement(rotate1, positionVector + G4ThreeVector(0.,cavityRadius4, -cavityLength1 + shift - cavityLength2 - 2*cavityLength3 - cavityLength4), log.back(), "Cavity4_1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);

    tubs.push_back(new G4Tubs("C4_side", 0., cavityRadius4 + cavityThickness4, cavityThickness4, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "C4Side_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., -cavityLength1 + shift - cavityLength2 - 2*cavityLength3 - 2*cavityLength4 - cavityThickness4), log.back(), "Cavity4Side", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray_50);



  // ________________________________________ mylar window _________________________________________

    tubs.push_back(new G4Tubs("mylarWin", 0., mylarWindowRadius, mylarWindowThick, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "MW_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., targetLength + shift - mylarWindowThick), log.back(), "MylarWindow", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->blue);

  // _______________________________________ mother log volumes for target caps _______________________________________
    
    // #1 referrs to the end of target, closest to mylar window
    // dimensions of mother volume don't matter much, just ensure they're big enough to contain all the solids
    tubs.push_back(new G4Tubs("mother1_vol", 0., 2*targetRadius + 2*caseThickness1, targetLength, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "mother1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength - 2*targetRadius), log.back(), "MotherVolLog_1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    CapMother1_log = log.back();

    tubs.push_back(new G4Tubs("mother2_vol", 0., 2*targetRadius + 2*caseThickness1, targetLength, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->vacuum_noOptical, "mother2_log"));
    new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,targetLength - 2*targetRadius), log.back(), "MotherVolLog_2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);
    CapMother2_log = log.back();

  // ___________________________________________ H2 target ___________________________________________

    // cylindrical part
    tubs.push_back(new G4Tubs("target_H2", 0., targetRadius, targetLength - 2*targetRadius, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->lh2, "targetH2_log"));
    addToTarget(new G4PVPlacement(0, positionVector, log.back(), "LH2Hadron", world_log, 0, 0, checkOverlap));
    log.back()->SetVisAttributes(colour->lightblue);

    // adding first semisphere to cylindrical lh2 target     
    sphere.push_back(new G4Sphere("LH2Sphere_1", 0., targetRadius, 0., M_PI, 0., M_PI));
    log.push_back(new G4LogicalVolume(sphere.back(), materials->lh2, "LH2BegCap_log"));

    new G4PVPlacement(rotate2, positionVector , log.back(), "LH2BegCap", CapMother1_log, 0, 0, checkOverlap);   
    log.back()->SetVisAttributes(colour->lightblue);              // - G4ThreeVector(0.,0.,targetLength - 2*targetRadius)

    unionLH2_1 = new G4UnionSolid("LH2sphere_1 + target_H2", sphere.back(), tubs.back());
    log.push_back(new G4LogicalVolume(unionLH2_1, materials->lh2, "targetCap1_log"));
    new G4PVPlacement(0, positionVector, log.back(), "targetCap_1", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);

    // adding second semisphere to cylindrical lh2 target
    sphere.push_back(new G4Sphere("LH2Sphere_2", 0., targetRadius, 0., M_PI, 0., M_PI));
    log.push_back(new G4LogicalVolume(sphere.back(), materials->lh2, "LH2EndCap_log"));

    new G4PVPlacement(rotate1, positionVector, log.back(), "LH2EndCap", CapMother2_log, 0, 0, checkOverlap);   // + G4ThreeVector(0.,0., targetLength - 2*targetRadius)
    log.back()->SetVisAttributes(colour->lightblue);

    unionLH2_2 = new G4UnionSolid("LH2Sphere_2 + unionLH2_1", sphere.back(), unionLH2_1);   //unionLH2_1
    log.push_back(new G4LogicalVolume(unionLH2_2, materials->lh2, "targetCap2_log"));
    new G4PVPlacement(0, positionVector, log.back(), "targetCap_2", world_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->invisible);

    // ___________________________________________ mylar caps ___________________________________________

    // position of mylar caps depends on definition of targetLength
    // mylar cap at the beginning of target
    sphere.push_back(new G4Sphere("BegCapSphere", mylarRadius, mylarRadius + mylarThickness, 0., M_PI, 0., M_PI));    
    log.push_back(new G4LogicalVolume(sphere.back(), materials->mylar, "BegCap_log"));
    new G4PVPlacement(rotate1, positionVector, log.back(), "BegCap", CapMother2_log, 0, 0, checkOverlap);      // - G4ThreeVector(0.,0.,targetLength - 2*targetRadius)
    log.back()->SetVisAttributes(colour->white);

    tubs.push_back(new G4Tubs("cylinder_mylar1", mylarRadius, mylarRadius + mylarThickness, mylarCylLength, 0., 2.*M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "mylarCyl1_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,mylarCylLength), log.back(), "myCyl1", CapMother2_log, 0, 0, checkOverlap);      // - G4ThreeVector(0.,0.,targetLength - 2*targetRadius - mylarCylLength)
    log.back()->SetVisAttributes(colour->white);

  //  union1 = new G4UnionSolid("semisphere1 + cylinder1", sphere.back(), tubs.back());                              
  //  log.push_back(new G4LogicalVolume(union1, materials->mylar, "mylarU1_log"));
  //  new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,targetLength - mylarCylLength), log.back(), "myU1", targetMother_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->red_50);

    // mylar cap at the end of the target
    sphere.push_back(new G4Sphere("EndCapSphere", mylarRadius, mylarRadius + mylarThickness, 0., M_PI, 0., M_PI));    
    log.push_back(new G4LogicalVolume(sphere.back(), materials->mylar, "EndCap_log"));
    new G4PVPlacement(rotate2, positionVector, log.back(), "EndCap", CapMother1_log, 0, 0, checkOverlap);        // + G4ThreeVector(0.,0., targetLength - 2*targetRadius)
    log.back()->SetVisAttributes(colour->white);

    tubs.push_back(new G4Tubs("cylinder_mylar2", mylarRadius, mylarRadius + mylarThickness, mylarCylLength, 0., 2.* M_PI));
    log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "mylarCyl2_log"));
    new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,mylarCylLength), log.back(), "myCyl2", CapMother1_log, 0, 0, checkOverlap);     // 0.,0.,targetLength - 2*targetRadius - mylarCylLength
    log.back()->SetVisAttributes(colour->white);

  //  union2 = new G4UnionSolid("semisphere2 + cylinder2", sphere.back(), tubs.back());                            // second mylar cap at end of target
  //  log.push_back(new G4LogicalVolume(union2, materials->mylar, "mylarU2_log"));
  //  new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength - mylarCylLength/2), log.back(), "myU2", targetMother_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->white);

    // ___________________________________________ kapton cylinder around target __________________________________________

    tubs.push_back(new G4Tubs("kapton_tube", kaptonRadius, kaptonRadius + kaptonThickness, kaptonLength, 0., 2 * M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->kapton, "kapton_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., + steelLength + kaptonLength + mylarCylLength), log.back(), "KaptonTube", CapMother2_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->orange_50);

    // __________________________________________ stainless steel cylinder around target ________________________________

    tubs.push_back(new G4Tubs("steelCylinder", steelRadius, steelRadius + steelThickness, steelLength, 0., 2.*M_PI));      
    log.push_back(new G4LogicalVolume(tubs.back(), materials->stainlessSteel, "steelCyl_log"));
    new G4PVPlacement(0, positionVector + G4ThreeVector(0., 0., +steelLength + mylarCylLength), log.back(), "steelCyl", CapMother2_log, 0, 0, checkOverlap);
    log.back()->SetVisAttributes(colour->darkgray);


  // ________________________________________ target case __________________________________________

    // first cap of case at beginning of target
  //  G4Tubs* cylinder = new G4Tubs("case_BegCap1", targetRadius, targetRadius + mylarThickness + caseThickness1, caseLength5, 0., 2.*M_PI);
  //  log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "caseCapCyl_log1"));      // material to be checked
  //  new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0., targetRadius + mylarThickness - caseLength5), log.back(), "targetCaseCyl1", CapMother2_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->violet);

  //  tubs.push_back(new G4Tubs("case_BegCap2", 0., targetRadius + mylarThickness + caseThickness1, caseThickness2, 0., 2.*M_PI));
  //  log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "caseCapDisc_log1"));      // material to be checked
  //  new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0., targetRadius + mylarThickness + caseThickness2), log.back(), "targetCaseDisc1", CapMother2_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->orange);

  //  unionCaseCap = new G4UnionSolid("disc + cylinder", tubs.back(), cylinder);
  //  log.push_back(new G4LogicalVolume(unionCaseCap, materials->mylar, "unionCaseCap_log"));        // material to be checked
  //  new G4PVPlacement(0, positionVector - G4ThreeVector(0.,0.,targetLength + mylarThickness - caseThickness2), log.back(), "caseCap", caseMother_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->green);

    // rest of case
  //  G4Tubs* cylinder1 = new G4Tubs("case_EndCap1", targetRadius, targetRadius + mylarThickness + caseThickness1, caseLength1, 0., 2.*M_PI);
  //  log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "caseCapCyl_log2"));      // material to be checked
  //  new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., targetRadius + mylarThickness - caseLength1), log.back(), "targetCaseCyl2", CapMother1_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->blue);

  //  tubs.push_back(new G4Tubs("case_EndCap2", 0., targetRadius + mylarThickness + caseThickness1, caseThickness2, 0., 2.*M_PI));
  //  log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "caseCapDisc_log2"));      // material to be checked
  //  new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetRadius + mylarThickness + caseThickness2), log.back(), "targetCaseDisc2", CapMother1_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->magenta);

   // G4UnionSolid* unionCaseCap1 = new G4UnionSolid("disc + cylinder", cylinder1, tubs.back());
   // log.push_back(new G4LogicalVolume(unionCaseCap1, materials->mylar, "unionCaseCap1_log"));        // material to be checked
   // new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0.,targetLength + mylarThickness + caseThickness2), log.back(), "caseCap1", caseMother_log, 0, 0, checkOverlap);
   // log.back()->SetVisAttributes(colour->magenta);

  //  G4Tubs* cylinder1_1 = new G4Tubs("casePiece1_1", targetRadius, targetRadius + kaptonThickness + caseThickness1_1, caseLength1_1, 0., 2.*M_PI);
  //  log.push_back(new G4LogicalVolume(tubs.back(), materials->mylar, "caseCapCyl1_1_log"));      // material to be checked
  //  new G4PVPlacement(0, positionVector + G4ThreeVector(0.,0., - 2*caseLength1 - caseCone - caseLength1_1), log.back(), "targetCaseCyl1_1", CapMother1_log, 0, 0, checkOverlap);
  //  log.back()->SetVisAttributes(colour->yellow);

  // G4trapezoid can have at most 8 vertices in total, this has 12

  //____________________________ triangular-ish target holders ____________________________

/*  double trapHeight = trapSide * TMath::Sin(M_PI/3) + trapDx * TMath::Cos(M_PI/6);          // height of target holder
  
  // vertices coordinates must be given in anti-clockwise order, from bottom to top.
  double X1 = trapDx;
  double Y1 = trapHeight/2;
  double X2 = trapSide * TMath::Cos(M_PI/3) + trapDx;
  double Y2 = trapSide/2 * TMath::Sin(M_PI/3);
  double X3 = trapDx + trapSide*TMath::Cos(M_PI/3) - trapDx*TMath::Sin(M_PI/6);

   vector<G4TwoVector> vertices = {{-X1,-Y1},{-X2,-Y2},{-X3,Y1},{X3,Y1},{X2,Y2},{X1,-Y1}};
   trap.push_back(new G4GenericTrap("targetHolder1", trapThickness, vertices));
*/ // new G4PVPlacement(0, positionVector, log.back(), "Holder1", targetMother_log, 0, 0, checkOverlap);     // da sistemare la posizione

  // try intersecating 3 triangles
  // mock triangle just to see if it works
/*  vector<G4TwoVector> vertices = {{0.,-20},{20,20},{-20,20},{-20,20},{0.,-20},{20,20},{-20,20},{-20,20}};
  trap.push_back(new G4GenericTrap("targetHolder1", trapThickness, vertices));
   log.push_back(new G4LogicalVolume(trap.back(), materials->stainlessSteel, "holder1_log"));
   new G4PVPlacement(0, positionVector - G4ThreeVector(0., 0., cavityLength1 + 2* cavityLength2), log.back(), "Holder1", targetMother_log, 0, 0, checkOverlap);
   
   // target holder hole
  tubs.push_back(new G4Tubs("targetHolderHole", 0., mylarRadius, trapThickness, 0., 2*M_PI));
  new G4PVPlacement(0, positionVector  - G4ThreeVector(0., 0., cavityLength1 + 2* cavityLength2), log.back(), "Hole1", targetMother_log, 0, 0, checkOverlap);        // stessa posizione del trap
  
  subtraction1 = new G4SubtractionSolid("trap - cylinder",trap.back(),tubs.back());         // last two arguments are rotationMatrix and Translation
  log.push_back(new G4LogicalVolume(subtraction1, materials->stainlessSteel, "holder1_log"));
  new G4PVPlacement(0, positionVector - G4ThreeVector(0., 0., cavityLength1 + 2* cavityLength2), log.back(), "Hole1", targetMother_log, 0, 0, checkOverlap);
//  log.back()->SetVisAttributes(colour->green);
*/
  // box.push_back(new G4Box("box1", xBox, yBox, trapThickness));
  // box.push_back(new G4Box("box2", xBox, yBox, trapThickness));
  // box.push_back(new G4Box("box3", xBox, yBox, trapThickness));

}



void T4LH2Target2024::getTargetDetDat(std::vector<T4STargetInformation> &targetInformation) {
  T4STargetInformation target;
  target.name = "targetH2_2024";
  target.rotMatrix = TGEANT::ROT_XtoZ;
  target.sh = 5;
  // target.xSize = targetRadius;
  // target.ySize = targetLength * 2.;
  target.zSize = 0;
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
}