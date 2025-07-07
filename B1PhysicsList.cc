
//
#include "B1PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh" // Para G4BestUnit


// Para registrar constructores de partículas
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

B1PhysicsList::B1PhysicsList()
: G4VUserPhysicsList(),
  fEmPhysicsList(0),
  fHadronPhysicsList(0),
  fNeutronTrackingCut(0),
  fDecayPhysicsList(0),
  fRadioactiveDecayPhysicsList(0)
{
  G4int verb = 1; // Nivel de verbosidad (0=silencioso, 1=estándar, 2=detallado)
  SetVerboseLevel(verb);

  // Inicializar los constructores de física
  fEmPhysicsList = new G4EmStandardPhysics_option4(verb);
  fHadronPhysicsList = new G4HadronPhysicsQGSP_BERT_HP(verb); // ¡Crucial para fisión y neutrones térmicos!
  fNeutronTrackingCut = new G4NeutronTrackingCut(verb);
  fDecayPhysicsList = new G4DecayPhysics(verb);
  fRadioactiveDecayPhysicsList = new G4RadioactiveDecayPhysics(verb); // Para desintegración de productos de fisión
}

B1PhysicsList::~B1PhysicsList()
{
  delete fEmPhysicsList;
  delete fHadronPhysicsList;
  delete fNeutronTrackingCut;
  delete fDecayPhysicsList;
  delete fRadioactiveDecayPhysicsList;
}

void B1PhysicsList::ConstructParticle()
{
  // Registra todas las partículas predefinidas por Geant4
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor; // Importante para el neutrón
  pShortLivedConstructor.ConstructParticle();
}

void B1PhysicsList::ConstructProcess()
{
  // En B1PhysicsList::ConstructProcess()

auto particleIterator = GetParticleIterator();
particleIterator->reset();

while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    // ¡Añade el transporte a TODAS las partículas!
    pmanager->AddProcess(new G4Transportation(), -1, 0, 0);

    // ... aquí añades los otros procesos (ionización, compton, etc.)
}

  // Añadir la física electromagnética
  fEmPhysicsList->ConstructProcess();

  // Añadir la física hadrónica (incluye dispersión elástica, inelástica, fisión, captura)
  fHadronPhysicsList->ConstructProcess();

  // Añadir el corte de tracking de neutrones (optimización)
  fNeutronTrackingCut->ConstructProcess();

  // Añadir procesos de desintegración de partículas
  fDecayPhysicsList->ConstructProcess();

  // Añadir procesos de desintegración radiactiva de núcleos (para productos de fisión)
  fRadioactiveDecayPhysicsList->ConstructProcess();
}

void B1PhysicsList::SetCuts()
{
  // Define los cortes de producción para partículas secundarias
  // Es importante ajustar esto según la precisión y eficiencia deseadas.
  // Un corte pequeño significa más precisión, pero más tiempo de cálculo.
  G4double defaultCut = 0.1 * mm; // Puedes ajustar este valor
  SetDefaultCutValue(defaultCut);

  // Puedes establecer cortes específicos para diferentes partículas si es necesario:
  // SetCutValue(1.0*um, "gamma");
  // SetCutValue(0.1*mm, "e-");
  // SetCutValue(0.1*mm, "e+");

  //G4cout << "Cuts for particles: " << G4endl;
  //DumpCutValuesTable(); // Imprime los valores de corte actuales
}


