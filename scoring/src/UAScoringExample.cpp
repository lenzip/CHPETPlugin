/*
 * UAScoringExample.cpp
 *
 *  Created on: 30 Aug 2013
 *      Author: Nicola mori
 */

#include "../include/UAScoringExample.h"
#include "../include/HitScoringExample.h"

#include "montecarlo/pluginmanagers/GGSMCPluginMacros.h"

//#include "G4GenericMessenger.hh"
#include "G4Event.hh"

RegisterUA(UAScoringExample);

UAScoringExample::UAScoringExample() :
    _outFileName("UAScoringExample.txt"), _outFile(NULL) {

  G4cout << "[UAScoringExample::UAScoringExample]" << G4endl;
  _messenger = new G4GenericMessenger(this, "/GGS/examples/scoring/UAScoringExample/");
  _messenger->DeclareProperty("outFileName", _outFileName, "Name of the output file").SetStates(G4State_Idle, G4State_PreInit);
}

UAScoringExample::~UAScoringExample() {
  G4cout << "[UAScoringExample::~UAScoringExample]" << G4endl;
  delete _messenger;
}

void UAScoringExample::BeginOfRunAction(const G4Run *run) {
  G4cout << "[UAScoringExample::BeginOfRunAction]" << G4endl;

  _outFile = new std::ofstream(_outFileName);
  (*_outFile) << "Total energy release for each event (GeV)\n";
  (*_outFile) << "-----------------------------------------\n";
}

void UAScoringExample::BeginOfEventAction(const G4Event *event) {

}

void UAScoringExample::EndOfEventAction(const G4Event *event) {
  (*_outFile) << "Event " << event->GetEventID() << ":\n";

  G4HCofThisEvent *hcEvent = event->GetHCofThisEvent();
  for (G4int iColl = 0; iColl < hcEvent->GetNumberOfCollections(); iColl++) {
    HitsCollection *hitColl = dynamic_cast<HitsCollection*>(hcEvent->GetHC(iColl));
    if (hitColl && hitColl->GetName().contains(".SDScoringExample")) {
      for (size_t iHitScoringExample = 0; iHitScoringExample < hitColl->GetSize(); iHitScoringExample++) {
        HitScoringExample *hit = (HitScoringExample*) (hitColl->GetHit(iHitScoringExample));
        (*_outFile) << "  " << hit->GetVolName() << " " << hit->GetRelease() / GeV << "\n";
      }
    }
  }
}

void UAScoringExample::EndOfRunAction(const G4Run *run) {
  G4cout << "[UAScoringExample::EndOfRunAction]" << G4endl;
  _outFile->close();
  delete _outFile;
}
