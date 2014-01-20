/*
 * UserActionExample.cpp
 *
 *  Created on: 29 Aug 2013
 *      Author: Nicola mori
 */

#include "../include/UserActionExample.h"

#include "montecarlo/pluginmanagers/GGSMCPluginMacros.h"

//#include "G4GenericMessenger.hh"
#include "G4Event.hh"

RegisterUA(UserActionExample);

UserActionExample::UserActionExample() :
    _outFileName("UserActionExample.txt"), _outFile(NULL) {

  G4cout << "[UserActionExample::UserActionExample]" << G4endl;
  _messenger = new G4GenericMessenger(this, "/GGS/examples/useraction/UserActionExample/");
  _messenger->DeclareProperty("outFileName", _outFileName, "Name of the output file").SetStates(G4State_Idle, G4State_PreInit);
}

UserActionExample::~UserActionExample() {
  G4cout << "[UserActionExample::~UserActionExample]" << G4endl;
  delete _messenger;
}

void UserActionExample::BeginOfRunAction(const G4Run *run) {
  G4cout << "[UserActionExample::BeginOfRunAction]" << G4endl;

  _outFile = new std::ofstream(_outFileName);
  (*_outFile) << "Number of particles in each event\n";
  (*_outFile) << "---------------------------------\n";
}

void UserActionExample::BeginOfEventAction(const G4Event *event) {

}

void UserActionExample::EndOfEventAction(const G4Event *event) {
  (*_outFile) << "Event " << event->GetEventID() << ": " <<  event->GetTrajectoryContainer()->size() << " particles \n";
}

void UserActionExample::EndOfRunAction(const G4Run *run) {
  G4cout << "[UserActionExample::EndOfRunAction]" << G4endl;
  _outFile->close();
  delete _outFile;
}
