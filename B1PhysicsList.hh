// src/B1PhysicsList.hh
#ifndef B1PhysicsList_h
#define B1PhysicsList_h 1

#include "G4VUserPhysicsList.hh"

// Incluye los constructores de física necesarios
#include "G4EmStandardPhysics_option4.hh"      // Física EM estándar, opción 4 (recomendada para baja energía)
#include "G4HadronPhysicsQGSP_BERT_HP.hh"   // Física hadrónica: QGSP_BERT con neutrones de alta precisión (HP)
#include "G4NeutronTrackingCut.hh"          // Para cortar el tracking de neutrones de baja energía
#include "G4DecayPhysics.hh"                // Para la desintegración de partículas
#include "G4RadioactiveDecayPhysics.hh"     // Para la desintegración radiactiva de núcleos

class B1PhysicsList : public G4VUserPhysicsList
{
public:
  B1PhysicsList();
  virtual ~B1PhysicsList();

  // Métodos virtuales de G4VUserPhysicsList
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();

private:
  // Punteros a los constructores de física
  G4EmStandardPhysics_option4* fEmPhysicsList;
  G4HadronPhysicsQGSP_BERT_HP* fHadronPhysicsList;
  G4NeutronTrackingCut* fNeutronTrackingCut;
  G4DecayPhysics* fDecayPhysicsList;
  G4RadioactiveDecayPhysics* fRadioactiveDecayPhysicsList;
};

#endif
