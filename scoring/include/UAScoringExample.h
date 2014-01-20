/*
 * UAScoringExample.h
 *
 *  Created on: 30 Aug 2013
 *      Author: Nicola Mori
 */

#ifndef UASCORINGEXAMPLE_H_
#define UASCORINGEXAMPLE_H_

#include "montecarlo/useractions/GGSUserAction.h"

#include "G4String.hh"

#include <ostream>

class G4GenericMessenger;

class UAScoringExample: public GGSUserAction {

public:

  UAScoringExample();
  ~UAScoringExample();

  void BeginOfEventAction(const G4Event *event);
  void EndOfEventAction(const G4Event *event);

  void BeginOfRunAction(const G4Run *run);
  void EndOfRunAction(const G4Run *run);

private:

  G4GenericMessenger *_messenger;
  G4String _outFileName;
  std::ofstream *_outFile;
};

#endif /* UASCORINGEXAMPLE_H_ */
