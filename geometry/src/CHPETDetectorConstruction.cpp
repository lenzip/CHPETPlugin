#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Transform3D.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4GenericMessenger.hh"
#include "G4UImanager.hh"

#include "../include/CHPETDetectorConstruction.h"
#include "geometry/pluginmanagers/GGSGeoPluginMacros.h"

G4double* wavelenghToEnergy(const G4double* wavelenghts, G4int n){
  G4double* energies = new G4double[n];
  G4cout << m << G4endl;
  const G4double hc = 1239.80*MeV*1e-15*m;
  for (G4int i = 0; i < n; ++i){
    *(energies+i) = (hc / *(wavelenghts+i));
    G4cout << "i: " << i <<" lambda: " << wavelenghts[i]/nm << " with energy: " << *(energies+i)/eV << G4endl;
  }
  return energies;
}


GeometryPlugin(CHPETDetectorConstruction);
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CHPETDetectorConstruction::CHPETDetectorConstruction()
{
 _messenger = new G4GenericMessenger(this, "/GGS/geometry/CHPETDetectorConstruction/");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CHPETDetectorConstruction::~CHPETDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* CHPETDetectorConstruction::Construct()
{

  G4NistManager* man = G4NistManager::Instance();

  if (_geoDataCard != "") {
      G4UImanager::GetUIpointer()->ApplyCommand(G4String("/control/execute " + _geoDataCard));
   }

//	------------- Materials -------------

  G4double a, z, density;
  G4int nelements;
  G4double fraction;

// Air
//
  G4Material * Air = man->FindOrBuildMaterial("G4_AIR");;

// siuv 
//
  //numbers from https://www.escooptics.com/material-data/s1-uv-ultraviolet-grade-fused-silica.html
  G4Material * lead = man->FindOrBuildMaterial("G4_Pb");
  G4Material * SiO2 = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  density = 2.203*g/cm3;

// Silicon
// 
  G4Material * SiliconMaterial = new G4Material("SiliconMaterial", z=14, 
					       a=28.09*g/mole, 
					       density=2.329*g/cm3);

// Optical Resin
//
  G4Material * ResinMaterial = new G4Material("ResinMaterial", z=6, a=12*g/mole, 
					     density=1*g/cm3);

// Ceramic
//
  G4Material * CeramicMaterial = new G4Material("ResinMaterial", z=6, a=12*g/mole, 
					     density=2*g/cm3);
//lead glass
  G4Material * PbO  = man->FindOrBuildMaterial("G4_LEAD_OXIDE");
  G4Material * K2O  = man->FindOrBuildMaterial("G4_POTASSIUM_IODIDE");
  G4Material * Na2O = man->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");
  density = 4.06*g/cm3;
  /*
  G4Material * Lead_Glass = new G4Material("LEAD_GLASS_CEREN25", density, nelements=4);
  Lead_Glass->AddMaterial(SiO2,fraction=0.39);
  Lead_Glass->AddMaterial(PbO,fraction=0.55);
  Lead_Glass->AddMaterial(K2O,fraction=0.03);
  Lead_Glass->AddMaterial(Na2O,fraction=0.03); 
  */
  G4Material * Lead_Glass = new G4Material("LEAD_GLASS_CEREN25", density, nelements=2);
  Lead_Glass->AddMaterial(SiO2,fraction=0.70);
  Lead_Glass->AddMaterial(PbO,fraction=0.30);

//
// ------------ Generate & Add Material Properties Table ------------
//
  const G4int nEntries_SiO2 = 12;
  G4double wavelenghs_SiO2[nEntries_SiO2] =
    { 180.*nm, 214.*nm, 230.*nm, 240*nm, 265*nm, 280*nm, 334*nm, 320*nm, 347*nm, 365*nm, 405*nm, 436*nm };

  //G4double refractiveIndex_Si02[nEntries_SiO2] =
  //  {1.53,     1.53,    1.52,    1.51,   1.50,    1.49,  1.49,   1.48,   1.48,    1.47,  1.47,   1.47   };
  G4double refractiveIndex_Si02[nEntries_SiO2] =
  { 1.64,    1.63,    1.61,     1.61,   1.59,   1.58, 1.57,    1.57,   1.57,   1.56,   1.56,   1.55};

  G4double absorption_SiO2[nEntries_SiO2] =
    { 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 250.*cm, 
      250.*cm, 250.*cm, 250.*cm, 250.*cm };

  G4double scintillation[nEntries_SiO2] = { 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

  G4MaterialPropertiesTable * SiO2MPT = new G4MaterialPropertiesTable();
  
  G4double* energies_SiO2 = wavelenghToEnergy(wavelenghs_SiO2, nEntries_SiO2);

  SiO2MPT->AddProperty("RINDEX", energies_SiO2, refractiveIndex_Si02, nEntries_SiO2);
  SiO2MPT->AddProperty("ABSLENGTH", energies_SiO2, absorption_SiO2, nEntries_SiO2);
  //SiO2MPT->AddProperty("FASTCOMPONENT", energies_SiO2, scintillation, nEntries_SiO2);
  //SiO2MPT->AddConstProperty ("SCINTILLATIONYIELD", 1000./MeV);
  //SiO2MPT->AddConstProperty("RESOLUTIONSCALE",1.0);
  //SiO2MPT->AddConstProperty("FASTTIMECONSTANT",2.*ns);
  //SiO2MPT->AddConstProperty("YIELDRATIO",1.0);

  SiO2->SetMaterialPropertiesTable(SiO2MPT);
  Lead_Glass->SetMaterialPropertiesTable(SiO2MPT);

  // Set the Birks Constant for the scintillator
  //Lead_Glass->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

//
// Air
//
  G4double RefractiveIndex_Air[nEntries_SiO2] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* airMPT = new G4MaterialPropertiesTable();
  airMPT->AddProperty("RINDEX", energies_SiO2, RefractiveIndex_Air, nEntries_SiO2);
  
  Air->SetMaterialPropertiesTable(airMPT);

//
// Optical Resin
//
  G4double RefractiveIndex_Resin[nEntries_SiO2] =
   //{ 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55 };
  {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  G4double Absorption_Resin[nEntries_SiO2] =
    {100*cm,100*cm,100.*cm,100.*cm,100.*cm,100.*cm,100.*cm,100*cm,100*cm,100*cm,100*cm,100*cm};

  G4MaterialPropertiesTable * resinMPT = new G4MaterialPropertiesTable();
  resinMPT->AddProperty("RINDEX", energies_SiO2, RefractiveIndex_Resin, nEntries_SiO2);
  resinMPT->AddProperty("ABSLENGTH", energies_SiO2, Absorption_Resin, nEntries_SiO2);
  
  ResinMaterial->SetMaterialPropertiesTable(resinMPT);

//
//	------------- Volumes --------------

// The experimental Hall
//
  G4double _World_Size_X    = 30*cm;
  G4double _World_Size_Y    = 30*cm;
  G4double _World_Size_Z    = 60*cm;

  G4double _Cube_Size_X     = 2*cm;
  G4double _Cube_Size_Y     = 3*cm;
  G4double _Cube_Size_Z     = 3*cm;


  G4double _Absorber_Size_Z = 0.*cm;
  G4double _Absorber_Size_X = _Cube_Size_X;
  G4double _Absorber_Size_Y = _Cube_Size_Y;

  G4double _Resin_Size_X = 1.4*cm;
  G4double _Resin_Size_Y = 1.4*cm;
  G4double _Resin_Size_Z = 2*mm;

  G4double abs_cub_gap = 0*cm;


  G4Box * _Solid_World = new G4Box("World", _World_Size_X / 2., 
			   _World_Size_Y / 2., _World_Size_Z / 2.);

  G4LogicalVolume * _Logic_World = new G4LogicalVolume(_Solid_World, Air, "World");

  G4PVPlacement * _Physical_World
    = new G4PVPlacement(0, G4ThreeVector(), _Logic_World ,"World", 0, false, 0);

  // the detector box
  G4Box * _detectorBox = new  G4Box ("DetectorBox", (_Cube_Size_X)/2., (_Cube_Size_Y)/2., (_Cube_Size_Z+_Absorber_Size_Z+_Resin_Size_Z)/2.);
  G4LogicalVolume * _logicDetector = new G4LogicalVolume(_detectorBox, Air, "DetectorBox");
  G4Box * _detectorBoxDetectorSide = new  G4Box ("DetectorBoxDetectorSide", (_Cube_Size_X)/2., (_Cube_Size_Y)/2., (_Cube_Size_Z+_Absorber_Size_Z)/2.);
  //G4LogicalVolume * _logicBoxDetectorSide = new G4LogicalVolume(_detectorBoxDetectorSide, Air, "DetectorBoxDetectorSide");
  //G4Box * _detectorBoxSensorSide = new  G4Box ("DetectorBoxSensorSide", _Cube_Size_X/2., _Cube_Size_Y/2., _Resin_Size_Z/2.);
  //G4LogicalVolume * _logicBoxSensorSide = new G4LogicalVolume(_detectorBoxSensorSide, Air, "SensorBoxSensorSide");

  //G4PVPlacement * placementBDS = new G4PVPlacement(0, G4ThreeVector(0, 0, -(_Resin_Size_Z)/2.), _logicBoxDetectorSide, "DetectorBoxDetectorSide", _logicDetector, false, 0);
  //G4PVPlacement * placementBSS = new G4PVPlacement(0, G4ThreeVector(0, 0, (_Cube_Size_Z+_Absorber_Size_Z)/2.), _logicBoxSensorSide, "SensorBoxSensorSide", _logicDetector, false, 0);



// The SiO2 cube
//      
  G4Box * _Solid_Cube = new G4Box("SiO2_Cube", _Cube_Size_X / 2.,
                                 _Cube_Size_Y / 2., _Cube_Size_Z / 2.);

  G4LogicalVolume * _Logic_Cube  = new G4LogicalVolume(_Solid_Cube, SiO2, 
                                                      "SiO2_Cube");

  G4double cube_Z = -(_Resin_Size_Z)/2.;

  G4PVPlacement * _Physical_Cube =
    new G4PVPlacement(0,G4ThreeVector(
                                      0.,
                                      0.,
                                      cube_Z
                                      ),
                      _Logic_Cube , "SiO2_Cube", _logicDetector, false, 0);
/*
// The Absorber 
//      
  G4Box * _Solid_Absorber = new G4Box("Pb_Absorber", _Absorber_Size_X / 2.,
                                 _Absorber_Size_Y / 2., _Absorber_Size_Z / 2.);

  G4LogicalVolume * _Logic_Absorber  = new G4LogicalVolume(_Solid_Absorber, Lead_Glass,
                                                      "Pb_Absorber");

  G4PVPlacement * _Physical_Absorber =
    new G4PVPlacement(0,G4ThreeVector(
                                      0.,
                                      0.,
                                      -_Cube_Size_Z/2.
                                      ),
                      _Logic_Absorber , "Pb_Absorber", _logicBoxDetectorSide, false, 0);  
*/
  
// Optical Resin

  //G4Box * _Solid_Resin = new G4Box("Optical_Resin", _Resin_Size_X / 2., 
			   //_Resin_Size_Y / 2., _Resin_Size_Z / 2.);
  G4Tubs * _Solid_Resin = new  G4Tubs("Optical_Resin", 0.*mm, 4*mm, _Resin_Size_Z / 2., 0*degree, 360*degree);                          

  G4LogicalVolume * _Logic_Resin  = new G4LogicalVolume(_Solid_Resin, 
				      ResinMaterial, 
				      "Optical_Resin");

  G4PVPlacement * _Physical_Resin = 
    new G4PVPlacement(0, G4ThreeVector(
				       0., 
				       0., 
				       (_Cube_Size_Z+_Absorber_Size_Z)/2.
				       ), 
		      _Logic_Resin , "Optical_Resin", _logicDetector, 
		      false, 0);

//Surfaces

  //SiO2_Air surface
  G4OpticalSurface * Op_SiO2_Surface = new G4OpticalSurface("SiO2_Surface");
  Op_SiO2_Surface->SetType(dielectric_dielectric);
  //Op_SiO2_Surface->SetFinish(groundbackpainted);
  Op_SiO2_Surface->SetFinish(groundfrontpainted);
  Op_SiO2_Surface->SetModel(unified);
  Op_SiO2_Surface->SetSigmaAlpha(0.3);
  Op_SiO2_Surface->SetPolish(1.);

  G4LogicalBorderSurface * SiO2_Air_Surf =
    new G4LogicalBorderSurface("SiO2_to_Air", _Physical_Cube,
                               _Physical_World, Op_SiO2_Surface);

  G4LogicalBorderSurface * Air_SiO2_Surf =
     new G4LogicalBorderSurface("Air_to_SiO2", _Physical_World,
                               _Physical_Cube, Op_SiO2_Surface);

  const G4int num = 2;

  G4double wavelenghts_resin[num] = {180.*nm, 436*nm};
  G4double* Ephoton = wavelenghToEnergy(wavelenghts_resin, num);
  G4double Reflectivity_SiO2[num] = {1., 1.};
  G4double Efficiency_SiO2[num]   = {0., 0.};

  G4double specularlobe_SiO2[num] = {1., 1.};
  G4double specularspike_SiO2[num] = {0., 0.};
  G4double backscatter_SiO2[num] = {0., 0.};

  G4MaterialPropertiesTable *CubeSPT = new G4MaterialPropertiesTable();

  CubeSPT->AddProperty("REFLECTIVITY", Ephoton, Reflectivity_SiO2, num);
  CubeSPT->AddProperty("EFFICIENCY",   Ephoton, Efficiency_SiO2,   num);
  CubeSPT->AddProperty("SPECULARLOBECONSTANT",
                        Ephoton, specularlobe_SiO2, num);
  CubeSPT->AddProperty("SPECULARSPIKECONSTANT",
                        Ephoton, specularspike_SiO2, num);
  CubeSPT->AddProperty("BACKSCATTERCONSTANT",
                        Ephoton, backscatter_SiO2, num);

  Op_SiO2_Surface->SetMaterialPropertiesTable(CubeSPT);

  // SiO2 - Resin Surface
  G4OpticalSurface * Op_Resin_Surface = new G4OpticalSurface("Resin_Surface");
  Op_Resin_Surface->SetType(dielectric_dielectric);
  Op_Resin_Surface->SetFinish(polished);
  Op_Resin_Surface->SetModel(unified);
  Op_Resin_Surface->SetSigmaAlpha(0.);
  Op_Resin_Surface->SetPolish(1.);

  G4LogicalBorderSurface * SiO2_Res_Surf =
    new G4LogicalBorderSurface("SiO2_to_Resin", _Physical_Cube,
                               _Physical_Resin, Op_Resin_Surface);

  G4LogicalBorderSurface * Res_SiO2_Surf =
     new G4LogicalBorderSurface("Resin_to_SiO2", _Physical_Resin,
                               _Physical_Cube, Op_Resin_Surface);

  //const G4int num = 2;

  //G4double wavelenghts_resin[num] = {180.*nm, 436*nm};
  //G4double* Ephoton = wavelenghToEnergy(wavelenghts_resin, num);
  G4double Reflectivity_Resin[num] = {1., 1.};
  G4double Efficiency_Resin[num]   = {0., 0.};

  G4double specularlobe_Resin[num] = {1., 1.};
  G4double specularspike_Resin[num] = {0., 0.};
  G4double backscatter_Resin[num] = {0., 0.};

  G4MaterialPropertiesTable *resinSPT = new G4MaterialPropertiesTable();

  resinSPT->AddProperty("REFLECTIVITY", Ephoton, Reflectivity_Resin, num);
  resinSPT->AddProperty("EFFICIENCY",   Ephoton, Efficiency_Resin,   num);
  resinSPT->AddProperty("SPECULARLOBECONSTANT",
                        Ephoton, specularlobe_Resin, num);
  resinSPT->AddProperty("SPECULARSPIKECONSTANT",
                        Ephoton, specularspike_Resin, num);
  resinSPT->AddProperty("BACKSCATTERCONSTANT",
                        Ephoton, backscatter_Resin, num);

  Op_Resin_Surface->SetMaterialPropertiesTable(resinSPT);

/*
  //Resin Air surace
  G4OpticalSurface * Op_Air_Surface = new G4OpticalSurface("Air_Surface");
  Op_Resin_Surface->SetType(dielectric_dielectric);
  Op_Resin_Surface->SetFinish(groundbackpainted);
  Op_Resin_Surface->SetModel(unified);
  Op_Resin_Surface->SetSigmaAlpha(0.);
  Op_Resin_Surface->SetPolish(1.);

  G4LogicalBorderSurface * Res_Air_Surf =
    new G4LogicalBorderSurface("Air_to_Resin", _Physical_Resin,
                               placementBSS, Op_Air_Surface);
  G4LogicalBorderSurface * Air_Res_Surf =
      new G4LogicalBorderSurface("Resin_to_Air", _Physical_World,
                                     _Physical_Resin, Op_Air_Surface);
  

  G4double Reflectivity_Air[num] = {0., 0.};
  G4double Efficiency_Air[num]   = {1., 1.};

  G4double specularlobe_Air[num] = {1., 1.};
  G4double specularspike_Air[num] = {0., 0.};
  G4double backscatter_Air[num] = {0., 0.};

  G4MaterialPropertiesTable *airSPT = new G4MaterialPropertiesTable();

  airSPT->AddProperty("REFLECTIVITY", Ephoton, Reflectivity_Air, num);
  airSPT->AddProperty("EFFICIENCY",   Ephoton, Efficiency_Air,   num);
  airSPT->AddProperty("SPECULARLOBECONSTANT",
                        Ephoton, specularlobe_Air, num);
  airSPT->AddProperty("SPECULARSPIKECONSTANT",
                        Ephoton, specularspike_Air, num);
  airSPT->AddProperty("BACKSCATTERCONSTANT",
                        Ephoton, backscatter_Air, num);

  Op_Air_Surface->SetMaterialPropertiesTable(airSPT); 
*/                      


  G4PVPlacement * _Physical_Detector = new G4PVPlacement(0, G4ThreeVector(0, 0, 5*cm), _logicDetector, "DetectorBox1", _Logic_World, false, 0);  
  G4RotationMatrix * rotm = new G4RotationMatrix();
  rotm->rotateX(180*deg);
  G4PVPlacement * _Physical_Detector2 = new G4PVPlacement(rotm, G4ThreeVector(0, 0, -5*cm), _logicDetector, "DetectorBox2", _Logic_World, false, 1);

//always return the physical World
  return _Physical_World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CHPETDetectorConstruction::SetReflectivity(G4double ref, 
					       G4OpticalSurface *surface, G4double& cubeReflectivity ){
  G4double ephot[2] = {1.65*eV, 4.28*eV};
  G4double reflect[2];
  reflect[0]=ref;
  reflect[1]=ref;
  G4MaterialPropertiesTable *mpt = surface->GetMaterialPropertiesTable();
  mpt->RemoveProperty("REFLECTIVITY");
  mpt->AddProperty("REFLECTIVITY", ephot, reflect, 2);
  surface->SetMaterialPropertiesTable(mpt);

  cubeReflectivity = ref;
}

