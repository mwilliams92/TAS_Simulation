//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4/B4c/src/DetectorConstruction.cc
/// \brief Implementation of the B4c::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace B4c
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Sodium Iodide
  G4Element* Na = nistManager->FindOrBuildElement("Na");
  G4Element* I  = nistManager->FindOrBuildElement("I");
  G4Element* Tl = nistManager->FindOrBuildElement("Tl");

  //G4double NaI_density = 3.67 * g/cm3;

  G4Material* NaITl = new G4Material("NaI_Tl", density= 3.67 * g/cm3, 3);
  NaITl->AddElement(Na, 0.4999);
  NaITl->AddElement(I,  0.4999);
  NaITl->AddElement(Tl, 0.0002); // example doping fraction for NaI

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers = 8;

  auto worldSizeXY = 1.*m; // * calorSizeXY;
  auto worldSizeZ  = 1*m; // * calorThickness;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");

  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operation
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps


  // Detector Module
  G4double ringThicknessZ = 10.*cm;

  G4double halfLength = ringThicknessZ/2.;
  G4double innerRadius = 4.5*cm;
  G4double outerRadius = 20.*cm;
  G4double startPhi = 0.*deg;
  G4double deltaPhi = 180.*deg;

  auto detectorS 
	= new G4Tubs(
    	"Detector", innerRadius, outerRadius, halfLength, startPhi, deltaPhi);

  auto detectorLV 
	= new G4LogicalVolume(detectorS, absorberMaterial, "DetectorLV");

  G4RotationMatrix* rotBottom = new G4RotationMatrix();
  rotBottom->rotateZ(180.*deg);

  for (G4int i = 0; i < 4; i++)
  {
    G4double zPos = (-1.5 + i) * ringThicknessZ;

    // Upper half
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0,0,zPos),
        detectorLV,
        "UpperDetector",
        worldLV,
        false,
        2*i,
        fCheckOverlaps);

    // Lower half
    new G4PVPlacement(
        rotBottom,
        G4ThreeVector(0,0,zPos),
        detectorLV,
        "LowerDetector",
        worldLV,
        false,
        2*i+1,
        fCheckOverlaps);
  }
  
  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  detectorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //
  auto detSD
    = new CalorimeterSD("DetectorSD", "DetectorHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(detSD);
  SetSensitiveDetector("DetectorLV",detSD);

  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
