/*
 * SDScoringExample.h
 *
 *  Created on: 30 Aug 2013
 *      Author: Nicola Mori
 */

#ifndef SDSCORINGEXAMPLE_H_
#define SDSCORINGEXAMPLE_H_

#include "G4VSensitiveDetector.hh"

#include "../include/HitScoringExample.h"

class G4GenericMessenger;

class SDScoringExample: public G4VSensitiveDetector {

public:
  SDScoringExample(G4String name);
  ~SDScoringExample();

  void Initialize(G4HCofThisEvent* hitCollection);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROHist);
  void EndOfEvent(G4HCofThisEvent* hitCollection);

private:

  G4GenericMessenger *_messenger;
  bool _onlyNeutral;
  HitsCollection *_hitsColl;
};

#endif /* SDSCORINGEXAMPLE_H_ */
