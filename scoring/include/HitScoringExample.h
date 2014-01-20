/*
 * HitScoringExample.h
 *
 *  Created on: 30/ago/2013
 *      Author: Nicola Mori
 */

#ifndef HITSCORINGEXAMPLE_H_
#define HITSCORINGEXAMPLE_H_

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4String.hh"

class HitScoringExample: public G4VHit {

public:
  HitScoringExample(const G4String &volName) :
      _energyRelease(0), _volName(volName) {
  }
  ~HitScoringExample() {
  }

  void AddRelease(double release);
  double GetRelease();

  const G4String& GetVolName();

  void* operator new(size_t);
  void operator delete(void* aHitScoringExample);

private:

  double _energyRelease;
  const G4String &_volName;
};

inline void HitScoringExample::AddRelease(double release) {
  _energyRelease += release;
}

inline double HitScoringExample::GetRelease() {
  return _energyRelease;
}

inline const G4String& HitScoringExample::GetVolName() {
  return _volName;
}

typedef G4THitsCollection<HitScoringExample> HitsCollection;
extern G4Allocator<HitScoringExample> HitScoringExampleAllocator;

inline void* HitScoringExample::operator new(size_t) {
  void* aHitScoringExample;
  aHitScoringExample = (void*) HitScoringExampleAllocator.MallocSingle();
  return aHitScoringExample;
}

inline void HitScoringExample::operator delete(void* aHitScoringExample) {
  HitScoringExampleAllocator.FreeSingle((HitScoringExample*) aHitScoringExample);
}
#endif /* HITSCORINGEXAMPLE_H_ */
