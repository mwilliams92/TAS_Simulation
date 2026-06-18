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
/// \file B4/B4c/src/EventAction.cc
/// \brief Implementation of the B4c::EventAction class

#include "EventAction.hh"
#include "CalorimeterSD.hh"
#include "CalorHit.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>

namespace B4c
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorHitsCollection*
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection
    = static_cast<CalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(
                              G4double absoEdep, G4double absoTrackLength) const
{
  // print event statistics
  G4cout
     << "   Absorber: total energy: "
     << std::setw(7) << G4BestUnit(absoEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SmearNaI(G4double E)
{
    // resolution ~ 7% at 662 keV
    const G4double a = 0.07;

    // sigma scaling (convert FWHM → sigma)
    G4double sigma = a * std::sqrt(E/0.662) * E / 2.355;

    G4double smearedE = G4RandGauss::shoot(E, sigma);

    if (smearedE < 0.) return 0.;
    return smearedE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{

  G4double theta[8] =
  {
    132.0*deg, 132.0*deg,   // ring 1 (segments)
    106.7*deg,106.7*deg, 	// ring 2
    42.0*deg, 42.0*deg, 	// ring 3
    16.7*deg, 16.7*deg  	// ring 4
  };

  const G4double beta = 0.3;
  const G4double gamma = 1.0 / std::sqrt(1 - beta*beta);

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Get hits collections IDs (only once)
  if ( fDetHCID == -1 ) {
    fDetHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("DetectorHitsCollection");
  }

  // Get hits collections
  auto detHC = GetHitsCollection(fDetHCID, event);

  // Get hit with total values
  auto detHit = (*detHC)[detHC->entries()-1];

  G4double etot = 0.;
  G4double etot_dop = 0.;

  for(G4int i=0; i< detHC->entries()-1; i++)
  {
	auto hit = (*detHC)[i];
	if (!hit) continue;

    G4double edep = hit->GetEdep();
    G4int id = hit->GetID();

	if (edep > 0.)
    {
        G4double smeared = SmearNaI(edep);
		G4double cosTheta = std::cos(theta[id]);
		G4double E_dop = smeared * gamma * (1.0 - beta * cosTheta);

		analysisManager->FillH1(1,smeared);
		analysisManager->FillNtupleDColumn(1, smeared);

		analysisManager->FillH1(3,E_dop);
		analysisManager->FillNtupleDColumn(3, E_dop);

		etot += smeared;
		etot_dop += E_dop;
	}

  }

  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;

    PrintEventStatistics(
      detHit->GetEdep(), detHit->GetTrackLength());
  }

  // Fill histograms, ntuple
  //

  // fill histograms
  analysisManager->FillH1(0, etot);
  analysisManager->FillH1(2, etot_dop);

  // fill ntuple
  analysisManager->FillNtupleDColumn(0, etot);
  analysisManager->FillNtupleDColumn(2, etot_dop);
  analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
