/*
 * SDScoringExample.cpp
 *
 *  Created on: 30 Aug 2013
 *      Author: Nicola Mori
 */

#include "../include/SDScoringExample.h"

#include "montecarlo/pluginmanagers/GGSMCPluginMacros.h"

//#include "G4GenericMessenger.hh"
#include "G4SDManager.hh"

RegisterSD(SDScoringExample);

SDScoringExample::SDScoringExample(G4String name) :
    G4VSensitiveDetector(name), _messenger(NULL), _onlyNeutral(false), _hitsColl(NULL) {

  _messenger = new G4GenericMessenger(this, G4String("/GGS/examples/scoring/SDScoringExample/").append(name).append("/"));
  _messenger->DeclareProperty("onlyNeutral", _onlyNeutral, "Save only energy deposited by neutral particles").SetDefaultValue(
      "true");
  collectionName.insert(name + ".SDScoringExample");
}

SDScoringExample::~SDScoringExample() {
  delete _messenger;
}

void SDScoringExample::Initialize(G4HCofThisEvent* hitCollection) {
  _hitsColl = new HitsCollection(SensitiveDetectorName, collectionName[0]);
  hitCollection->AddHitsCollection(0, _hitsColl);

}

G4bool SDScoringExample::ProcessHits(G4Step* aStep, G4TouchableHistory* ROHist) {

  if (_onlyNeutral && aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge() != 0.)
    return true;
  else {
    for (int iHit = 0; iHit < _hitsColl->entries(); iHit++) {
      HitScoringExample *hit = (*(_hitsColl))[iHit]; // The :: scope resolution operator is needed to avoid confusion with the method G4VSDScoringExampleHitScoringExample.
      if (hit->GetVolName() == aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()) { // This is an old hit
        hit->AddRelease(aStep->GetTotalEnergyDeposit());
        return true;
      }
    }

    // This is a new hit
    HitScoringExample* hit = new HitScoringExample(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());
    hit->AddRelease(aStep->GetTotalEnergyDeposit());
    _hitsColl->insert(hit);
    return true;
  }
}

void SDScoringExample::EndOfEvent(G4HCofThisEvent* hitCollection) {

  G4int IHID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hitCollection->AddHitsCollection(IHID, _hitsColl);
}
