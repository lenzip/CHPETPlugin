/*
 * UserActionExample.h
 *
 *  Created on: 29 Aug 2013
 *      Author: Nicola Mori
 */

#ifndef USERACTIONEXAMPLE_H_
#define USERACTIONEXAMPLE_H_

#include "montecarlo/useractions/GGSUserAction.h"

#include "G4String.hh"

#include <ostream>

class G4GenericMessenger;

class UserActionExample: public GGSUserAction {

public:

  UserActionExample();
  ~UserActionExample();

  void BeginOfEventAction(const G4Event *event);
  void EndOfEventAction(const G4Event *event);

  void BeginOfRunAction(const G4Run *run);
  void EndOfRunAction(const G4Run *run);

private:

  G4GenericMessenger *_messenger;
  G4String _outFileName;
  std::ofstream *_outFile;
};

#endif /* USERACTIONEXAMPLE_H_ */
