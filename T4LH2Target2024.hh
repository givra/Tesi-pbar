#ifndef T4LH2TARGET2024_HH_
#define T4LH2TARGET2024_HH_

#include "T4BaseDetector.hh"
#include "T4TargetBackend.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4GenericTrap.hh"     // generic trapezoid
#include "G4SubtractionSolid.hh"
#include <G4UnionSolid.hh>
#include "G4IntersectionSolid.hh"
#include "G4Cons.hh"

class T4LH2Target2024 : public T4BaseDetector, public T4TargetBackend
{
  public:
    T4LH2Target2024(T4SDetector*);
    virtual ~T4LH2Target2024(void);

    void construct(G4LogicalVolume*);

    void getTargetDetDat(std::vector<T4STargetInformation>&);
    
    void setCavityRadius1(double r1) {cavityRadius1 = r1;}        // referring to the biggest cavity
    void setCavityRadius2(double r2) {cavityRadius2 = r2;}        // at the left of #1
    //void setCavityRadius2(double r3) {cavityRadius2 = r3;}        // at the left of #2
    //void setCavityRadius2(double r4) {cavityRadius2 = r4;}        // at the left of #3
    void setCavityThickness1(double d1) {cavityThickness1 = d1;}
    void setCavityThickness2(double d2) {cavityThickness2 = d2;}
    //void setCavityThickness2(double d3) {cavityThickness2 = d3;}
    //void setCavityThickness2(double d4) {cavityThickness2 = d4;}
    //void setTargetRadius(double r) {targetRadius = r;}
    
    //G4LogicalVolume* getCapMother1_Log(void) {return CapMother1_log;}
    //G4LogicalVolume* getCapMother2_Log(void) {return CapMother2_log;}
    //G4LogicalVolume* getMotherTH1_log(void) {return MotherTH1_log;}


  private:
    G4double yearSetup;
    vector<G4Box*> box;
    vector<G4Tubs*> tubs;
    vector<G4Sphere*> sphere;
    vector<G4GenericTrap*> trap;
    vector<G4Cons*> cons;
    vector<G4LogicalVolume*> log;

    G4SubtractionSolid* subtractionEdge;        // used for target holders
    G4SubtractionSolid* subtractionHole;
    G4SubtractionSolid* subtractionC4;
    G4SubtractionSolid* subtractionTriangle1;
    G4SubtractionSolid* subtractionTriangle2;
    G4SubtractionSolid* subtractionTriangle3;
    G4SubtractionSolid* subtractionBox1;
    G4SubtractionSolid* subtractionBox2;
    G4SubtractionSolid* subtractionBox3;
    G4UnionSolid* union1;
    G4UnionSolid* union2;
    G4UnionSolid* unionC4;
    G4UnionSolid* unionLH2_1;
    G4UnionSolid* unionLH2_2;
    G4UnionSolid* unionCaseCap;

    G4LogicalVolume* CapMother1_log;
    G4LogicalVolume* CapMother2_log;
    G4LogicalVolume* MotherC4_log;
    G4LogicalVolume* MotherTH1_log;

    G4RotationMatrix* rotate1;
    G4RotationMatrix* rotate2;
    G4RotationMatrix* rotate120;    // used for target holders
    G4RotationMatrix* rotate60;    // used for target holders

    G4double cavityRadius1;        // bigger vacuum volume radius
    G4double cavityRadius2;        // smaller vacuum volume radius to the left of #1
    G4double cavityRadius3;        // to the left of #2
    G4double cavityRadius4;        // to the left of #3
    G4double cavityRadius4_1;      // upwards tube in cavity #4
    G4double cavityThickness1;
    G4double cavityThickness2;
    G4double cavityThickness3;
    G4double cavityThickness4;
    G4double cavityThickness4_1;
    G4double cavityLength1;
    G4double cavityLength2;
    G4double cavityLength3;
    G4double cavityLength4;
    G4double cavityLength4_1;

    G4double EndLength1;      // O-ring parameters around mylar window at the end of C1
    G4double EndRadiusMin1;
    G4double EndRadiusMax1;

    G4double EndLength_O;        // junction O-ring parameters between "End_1" and "End_2"
    G4double EndRadiusMin_O;
    G4double EndRadiusMax_O;

    G4double SideLength1;      // small O-ring parameters for cylinder 4
    G4double SideRadiusMin1;
    G4double SideRadiusMax1;

    G4double SideLength2;      // big O-ring parameters for cylinder 4
    G4double SideRadiusMin2;
    G4double SideRadiusMax2;

    G4double JunctionLength12_1;      // junction O-ring parameters between C1 and C2
    G4double JunctionRadiusMin12_1;
    G4double JunctionRadiusMax12_1;

    G4double JunctionLength12_O;        // junction O-ring parameters between "Junction12_1" and "Junction12_2"
    G4double JunctionRadiusMin12_O;
    G4double JunctionRadiusMax12_O;

    G4double JunctionLength12_2;     
    G4double JunctionRadiusMin12_2;
    G4double JunctionRadiusMax12_2;

    G4double JunctionLength23_1;      // junction O-ring parameters between C2 and C3
    G4double JunctionRadiusMin23_1;
    G4double JunctionRadiusMax23_1;

    G4double JunctionLength23_2;      // junction O-ring parameters between C2 and C3
    G4double JunctionRadiusMin23_2;
    G4double JunctionRadiusMax23_2;


    G4double shift1;               // distance from tip of target to Junction23_2
    G4double shift2;               // distance from tip of target to mylar window
    G4double dz1;                  // distance between big target holders
    G4double dz2;                  // distance up to small target holder

    G4double targetRadius;        // volume containing H2
    G4double targetLength;

      /* the length is measured from tip to tip of target which looks something like this:
         ______________________________
        /                              \
        \______________________________/
        the whole volume is described as a union of semisphere + cylinder + semisphere
        therefore if you just need the cylinder length always remember to -2*targetRadius
      */

    G4double kaptonLength;        // kapton cylinder around H2 volume
    G4double kaptonRadius;
    G4double kaptonThickness;
      
    G4double mylarRadius;         // mylar layer containing H2
    G4double mylarThickness;
    G4double mylarCylLength;      // length of mylar cylinder that will be unified with mylarCap_sphere

    G4double steelRadius;         // dimensions of stainless steel cylinder next to kapton
    G4double steelThickness;
    G4double steelLength;
    
    G4double mylarWindowRadius;           
    G4double mylarWindowThick;           
    G4double caseThickness;    // aluminium case thickness
    G4double caseLength1;
    G4double caseLength2;
    //G4double caseThickness1_1;    // aluminium case thickness in corrispondence of triangular-ish holder
    //G4double caseThickness2;    // thickness of cylindrical end
        
    //G4double caseCone;           // height of cone   
    //G4double caseLength1_1;     // small segment for Target Holder
    //G4double caseLength1;
    //G4double caseLength2;
    //G4double caseLength3;
    //G4double caseLength4;
    //G4double caseLength5;
    
    G4double trapThickness;       // thickness of target holder
    G4double trapDx;              // half of small edge of target holder
    G4double trapSide;            // side length of target holder
    G4double trapHeight;
    G4double HoleD;               // referring to hole in target holder; distance external side of rectangular hole - circular hole
    G4double xBox;                // sides for rectangular holes in target holder
    G4double yBox;
    

};

#endif