#ifndef T4MATERIALS_HH_
#define T4MATERIALS_HH_

#include "CLHEP/Units/SystemOfUnits.h"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4LogicalBorderSurface.hh"

#include "T4SettingsFile.hh"

/*! \class T4Materials
 *  \brief Singleton class to manage all materials.
 *
 *  This class is used in each T4BaseDetector class to give all
 *  logical volumes a material.
 */
class T4Materials
{
public:
        /*! \brief Get the instance of this class.*/
        static T4Materials* getInstance(void);
        /*! \brief Default destructor.*/
        ~T4Materials(void);

        /*! \brief Some materials with optical properties.*/
        G4Material* air_optical;
        G4Material* vacuum_optical;
        G4Material* bc408_optical;
        G4Material* bcf10_optical;
        G4Material* pmma_optical;
        G4Material* plexiglass_optical;
        G4Material* glass_optical;
        G4Material* aluminium_optical;
        G4Material* bialkali_optical;
        G4Material* c4f10_optical;
        G4Material* n2_optical;

        /*! \brief The same materials like above w/o optical properties.*/
        G4Material* air_noOptical;
        G4Material* vacuum_noOptical;
        G4Material* bc408_noOptical;
        G4Material* bcf10_noOptical;
        G4Material* pmma_noOptical;
        G4Material* plexiglass_noOptical;
        G4Material* glass_noOptical;
        G4Material* aluminium_noOptical;
        G4Material* bialkali_noOptical;
        G4Material* c4f10_noOptical;
        G4Material* n2_noOptical;

        /*! \brief Some materials.*/
        G4Material* neon;
        G4Material* borosilicate;
        G4Material* plasticFoam;
        G4Material* carbonFibre;
        G4Material* carbonLH2;
        G4Material* nomex;
        G4Material* mwpcGas;
        G4Material* epoxy;
        G4Material* epoxy_09;
        G4Material* concrete;
        G4Material* lh2;
        G4Material* ld2;
        G4Material* mylar;
        G4Material* nylon;
        G4Material* stainlessSteel;
        G4Material* polyethylene;
        G4Material* polypropylene;
        G4Material* polystyrene;
        G4Material* pvc;
        G4Material* airexC70;
        G4Material* iron;
        G4Material* silver;
        G4Material* helium;
        G4Material* lhe;
        G4Material* ammonia;
        G4Material* LiDeuterated;
        G4Material* nh3_dy;
        G4Material* nh3_he = NULL;
        G4Material* nh3_he_1st_2014; // Composition of first cell 2014
        G4Material* nh3_he_1st_2015; // Composition of first cell 2015
        G4Material* nh3_he_1st_2018; // Composition of first cell 2018
        G4Material* nh3_he_2nd_2014; // Composition of first cell 2014
        G4Material* nh3_he_2nd_2015; // Composition of first cell 2015
        G4Material* nh3_he_2nd_2018; // Composition of first cell 2018
        G4Material* nh3_he_fluka;
        G4Material* dli_he_2021;
        G4Material* C2F3Cl;
        G4Material* Li2CO3;
        G4Material* AlMylar;
        G4Material* argon;
        G4Material* ArCo2;
        G4Material* co2;
        G4Material* graphit;
        G4Material* Si;
        G4Material* Cu;
        G4Material* Nb;
        G4Material* Ti;
        G4Material* Ba;
        G4Material* BaTiO3;
        G4Material* NbTi;
        G4Material* NbTi_Cu;
        G4Material* superconducter;
        G4Material* kevlar;
        G4Material* Cu_10;
        G4Material* Cu_Kap;
        G4Material* Be;
        G4Material* tungsten;
        G4Material* nickel;
        G4Material* gold;
        G4Material* lead;
        G4Material* ethan;
        G4Material* cf4;
        G4Material* g10;
        G4Material* g10_10;
        G4Material* g11;
        G4Material* dc4gas;
        G4Material* dcEBox;
        G4Material* strawGas;
        G4Material* w45Gas;
        G4Material* w45CathodeGas;
        G4Material* mmSupport;
        G4Material* mdtGas;
        G4Material* mw2Gas;
        G4Material* kapton;
        G4Material* plastic;
        G4Material* rohacell;
        G4Material* rohacell_30;
        G4Material* inox;
        G4Material* siliconMaterial;
        G4Material* gemMaterial;
        G4Material* pixelGEMMaterial;
        G4Material* c2h6;
        G4Material* PbO;
        G4Material* SiO2;
        G4Material* K2O;
        G4Material* As2O3;
        G4Material* CeO2;
        G4Material* tf1; // gams
        G4Material* tf101; // gams rh
        G4Material* sf57; // mainz
        G4Material* sf5; // olga
        G4Material* Al2O3;
        G4Material* pctfe;
        G4Material* paper;
        G4Material* RICHWindowEff;
        G4Material* Klegecell;

        /*! \brief Some G4OpticalSurface used for surfaces between volumes with optical properties.*/
        G4OpticalSurface* surfaceBc408Air;
        G4OpticalSurface* surfaceAirBc408;
        G4OpticalSurface* surfaceBc408Plexiglass;
        G4OpticalSurface* surfacePlexiglassGlass;
        G4OpticalSurface* surfacePlexiglassAir;
        G4OpticalSurface* surfaceAirPlexiglass;
        G4OpticalSurface* surfaceGlassVacuum;
        G4OpticalSurface* surfaceVacuumBialkali;
        G4OpticalSurface* surfaceAluminium;
        G4OpticalSurface* surfaceMirror;
        G4OpticalSurface* surfaceBcf10;

        /*! \brief Reset and destruct the instance.*/
        static void resetInstance(void) {
                delete materials; materials = NULL;
        }

private:
        /*! \brief Private constructor.*/
        T4Materials(void);
        /*! \brief Private class pointer.*/
        static T4Materials* materials;
        /*! \brief True if FULL physics activated.*/
        bool useOptical;

        /*! \brief Define some indices.*/
        void DefineIndex(void);
        /*! \brief Build air and vacuum w/ and w/o optical properties.*/
        void BuildAir(void);
        /*! \brief Build Bc408 w/ and w/o optical properties.*/
        void BuildBc408(void);
        /*! \brief Build Bcf10 w/ and w/o optical properties.*/
        void BuildBcf10(void);
        /*! \brief Build PMMA w/ and w/o optical properties.*/
        void BuildPMMA(void);
        /*! \brief Build Plexiglass w/ and w/o optical properties.*/
        void BuildPlexiglass(void);
        /*! \brief Build Glass w/ and w/o optical properties.*/
        void BuildGlass(void);
        /*! \brief Build Aluminium w/ and w/o optical properties.*/
        void BuildAluminium(void);
        /*! \brief Build Bialkali w/ and w/o optical properties.*/
        void BuildBialkali(void);
        /*! \brief Build C4F10 w/ and w/o optical properties.*/
        void BuildC4F10(void);
        /*! \brief Build N2 w/ and w/o optical properties.*/
        void BuildN2(void);

        /*! \brief Build the optical surface between Bc408 and Air.*/
        void BuildSurfaceBc408Air(void);
        /*! \brief Build the optical surface between Air and Bc408.*/
        void BuildSurfaceAirBc408(void);
        /*! \brief Build the optical surface between Bc408 and Plexiglass.*/
        void BuildSurfaceBc408Plexiglass(void);
        /*! \brief Build the optical surface between Plexiglass and Glass.*/
        void BuildSurfacePlexiglassGlass(void);
        /*! \brief Build the optical surface between Plexiglass and Air.*/
        void BuildSurfacePlexiglassAir(void);
        /*! \brief Build the optical surface between Air and Plexiglass.*/
        void BuildSurfaceAirPlexiglass(void);
        /*! \brief Build the optical surface between Glass and Vacuum.*/
        void BuildSurfaceGlassVacuum(void);
        /*! \brief Build the optical surface between Vacuum and Bialkali.*/
        void BuildSurfaceVacuumBialkali(void);
        /*! \brief Build the optical surface for Aluminium.*/
        void BuildSurfaceAluminium(void);
        /*! \brief Build the optical surface for the RICH mirror.*/
        void BuildSurfaceMirror(void);
        /*! \brief Build the optical surface for Bcf10.*/
        void BuildSurfaceBcf10(void);

        pair<G4Material*,G4Material*> BuildDYtarget(double, double, int);

        /*! \brief Some elements to build up materials.*/
        G4Element* elH;
        G4Element* elC;
        G4Element* elN;
        G4Element* elO;
        G4Element* elF;
        G4Element* elFe;
        G4Element* elRb;
        G4Element* elMo;
        G4Element* elSb;
        G4Element* elCs;
        G4Element* elCr;
        G4Element* elNi;
        G4Element* elPb;
        G4Element* elCe;
        G4Element* elAs;
        G4Element* elK;
        G4Element* elSi;
        G4Element* elLi;
        G4Element* elCl;
        G4Element* elHe;
        G4Element* elBa;
        G4Element* elTi;

        /*! \brief G4NistManager pointer.*/
        G4NistManager* manager;

        /*! \brief Some indices.*/
        static const G4int nEntryLong = 12;
        static const G4int nEntryShort = 2;
        G4double bc408Index;
        G4double bc408AbsorptionLength;
        G4double bc408Yield;
        G4double plexiglassIndex;
        G4double plexiglassAbsLength;
        G4double glassIndex;
        G4double aluminiumReflectivity;
        G4double bcf10Index;
        G4double pmmaIndex;

        G4double specularLobe[nEntryShort];
        G4double specularSpike[nEntryShort];
        G4double backScatter[nEntryShort];

        G4double photonEnergyShort[nEntryShort];
        G4double photonEnergyLong[nEntryLong];

        G4double airRefractiveIndex[nEntryShort];
        G4double airAbsorpiton[nEntryShort];
        G4double plexiglassRefractiveIndex[nEntryShort];
        G4double plexiglassAbsorpiton[nEntryShort];
        G4double glassRefractiveIndex[nEntryShort];
        G4double bc408RefractiveIndex[nEntryShort];
        G4double bcf10RefractiveIndex[nEntryShort];
        G4double pmmaRefractiveIndex[nEntryShort];

        //DY target parameters 2014
        G4double packing_factor_1st_2014; //NH3 packing factor of first cell
        G4double packing_factor_2nd_2014; //NH3 packing factor of second cell
        //DY target parameters 2015
        G4double packing_factor_1st_2015; //NH3 packing factor of second cell
        G4double packing_factor_2nd_2015; //NH3 packing factor of second cell
        //DY target parameters 2018
        G4double packing_factor_1st_2018; //NH3 packing factor of first cell
        G4double packing_factor_2nd_2018; //NH3 packing factor of second cell

        /*! \brief Container for all G4MaterialPropertiesTable.*/
        std::vector<G4MaterialPropertiesTable*> materialPropertiesTable;
};

#endif /* T4MATERIALS_HH_ */