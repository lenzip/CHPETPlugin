#ifndef CaloCubeCHPETDetectorConstruction_h
#define CaloCubeCHPETDetectorConstruction_h 1

#include "geometry/GGSVGeometryConstruction.h"
#include "globals.hh"

class G4VPhysicalVolume;
class G4OpticalSurface;
class G4GenericMessenger;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.
class CHPETDetectorConstruction : public GGSVGeometryConstruction 
{
  public:
  enum {XPLUSSIDE, 
        XMINUSSIDE,
        YPLUSSIDE,
        YMINUSSIDE,
        ZPLUSSIDE,
        ZMINUSSIDE};

  public:
    CHPETDetectorConstruction();
    virtual ~CHPETDetectorConstruction();
    void SetReflectivity(G4double, G4OpticalSurface*, G4double&);

    G4VPhysicalVolume* Construct();
    G4VPhysicalVolume* GetVolume(){return _Physical_World;};

  private:
    G4LogicalVolume* buildPD(); 

    G4VPhysicalVolume * _Physical_World;
    G4double _sizeX, _sizeY, _sizeZ;
    G4int _detectorside;
    G4double _scintillationYield;
    G4double _Diode_Size_X;
    G4double _Diode_Size_Y;
    G4double _Diode_Size_Z;
    G4double _Resin_Size_X;
    G4double _Resin_Size_Y;
    G4double _Resin_Size_Z;
    G4double _Case_Size_X;
    G4double _Case_Size_Y;
    G4double _Case_Size_Z;

    G4PVPlacement * _Physical_Diode ;
    G4PVPlacement * _Physical_Resin ;
    G4LogicalVolume * _Logic_Resin;

    static G4double PhotonEnergy[];
    static G4double PhotonEnergy_2[];

    static G4double RefractiveIndex_CsI[];
    static G4double Absorption_CsI[];

    static G4double ScintilFast[];
    static G4double RefractiveIndex_Air[];

    G4GenericMessenger* _messenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

