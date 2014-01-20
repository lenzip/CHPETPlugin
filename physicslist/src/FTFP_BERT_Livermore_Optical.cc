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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   
//
// Author: 2007 Gunter Folger
//
//   created from FTFP
//
// Modified:
// 19.06.2008 G.Folger: don't use chips quasielastic in FTF
// 27.11.2009 G.Folger: Remobe experimental status
// 04.06.2010 G.Folger: Use new ctor for builders
// 16.08.2010 H.Kurashige: Remove inclusion of G4ParticleWithCuts 
// 26.06.2012 A.Ribon:  Use FTF/Preco and BERT for nuclear capture at rest.
// 27.07.2012 A.Ribon:  Use the new class G4BertiniAndFritiofStoppingPhysics
// 16.10.2012 A.Ribon:  Renamed the physics classes used
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4DataQuestionaire.hh"
#include "HadronPhysicsFTFP_BERT.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4OpticalPhysics.hh"
#include "../include/FTFP_BERT_Livermore_Optical.hh"

#include "montecarlo/pluginmanagers/GGSMCPluginMacros.h"

PhysicsListPlugin(FTFP_BERT_Livermore_Optical);


FTFP_BERT_Livermore_Optical::FTFP_BERT_Livermore_Optical(G4int ver) 
{
  G4cout << "GIULIO!!!!!!!!! Calling: FTFP_BERT_Livermore_Optical::FTFP_BERT_Livermore_Optical" << G4endl;
  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*CLHEP::mm;
  G4DataQuestionaire it(photon);
  G4cout <<G4endl;
  this->defaultCutValue = 0.7*CLHEP::mm;  
  this->SetVerboseLevel(ver);

 // EM Physics
  this->RegisterPhysics( new G4EmLivermorePhysics(ver));

  // Synchroton Radiation & GN Physics
  this->RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays 
  this->RegisterPhysics( new G4DecayPhysics(ver) );

   // Hadron Elastic scattering
  this->RegisterPhysics( new G4HadronElasticPhysics(ver) );

   // Hadron Physics
  this->RegisterPhysics(  new HadronPhysicsFTFP_BERT(ver));

  // Stopping Physics
  this->RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  this->RegisterPhysics( new G4IonPhysics(ver));
  
  // Neutron tracking cut
  this->RegisterPhysics( new G4NeutronTrackingCut(ver));

  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();

  opticalPhysics->SetWLSTimeProfile("delta");

  opticalPhysics->SetScintillationYieldFactor(1.);
  opticalPhysics->SetScintillationExcitationRatio(1.0);

  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  this->RegisterPhysics( opticalPhysics );
  this->SetVerboseLevel(500);

  G4int iphys =0;
  while (this->GetPhysics(iphys)){
    G4cout << "Registered " << this->GetPhysics(iphys)->GetPhysicsName() << G4endl;
    iphys++;
  }
}

FTFP_BERT_Livermore_Optical::~FTFP_BERT_Livermore_Optical()
{
}

void FTFP_BERT_Livermore_Optical::SetCuts()
{
  if (this->verboseLevel >1){
    G4cout << "FTFP_BERT::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 

  this->SetCutsWithDefault();   
 
//  if (this->verboseLevel > 0)
//    G4VUserPhysicsList::DumpCutValuesTable();  
}
