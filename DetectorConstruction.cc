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

//

/// \file B1/src/DetectorConstruction.cc

/// \brief Implementation of the B1::DetectorConstruction class



#include "DetectorConstruction.hh"



#include "G4RunManager.hh"

#include "G4NistManager.hh"

#include "G4Box.hh"

#include "G4Cons.hh"

#include "G4Orb.hh"

#include "G4Tubs.hh"


#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"


#include "G4LogicalVolume.hh"

#include "G4PVPlacement.hh"

#include "G4SystemOfUnits.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

namespace B1

{



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4VPhysicalVolume* DetectorConstruction::Construct()

{

  // Get nist material manager

  G4NistManager* nist = G4NistManager::Instance();


  // Option to switch on/off checking of volumes overlaps

  //

  G4bool checkOverlaps = true;


  //

  // Mundo

  //

  G4double rmin = 0*cm, rmax = 202*cm;

  G4double z = 1120*cm;
  
  G4double ai = 0*deg, af = 360*deg;

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");



  auto solidWorld = new G4Tubs("World",                           // its name

    rmin, rmax, 0.5*z, ai, af);  // its size



  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid

    world_mat,                                       // its material

    "World");                                        // its name



  auto physWorld = new G4PVPlacement(nullptr,  // no rotation

    G4ThreeVector(),                           // at (0,0,0)

    logicWorld,                                // its logical volume

    "World",                                   // its name

    nullptr,                                   // its mother  volume

    false,                                     // no boolean operation

    0,                                         // copy number

    checkOverlaps);                            // overlaps checking



  //

  // Piscina

  //
  
  G4double r1max = 200*cm;

  G4double z1 = 1100*cm;

  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  auto solidEnv = new G4Tubs("Envelope",                    // its name

    rmin, r1max, 0.5*z1, ai, af);  // its size



  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid

    env_mat,                                     // its material

    "Envelope");                                 // its name



  new G4PVPlacement(nullptr,  // no rotation

    G4ThreeVector(),          // at (0,0,0)

    logicEnv,                 // its logical volume

    "Envelope",               // its name

    logicWorld,               // its mother  volume

    false,                    // no boolean operation

    0,                        // copy number

    checkOverlaps);           // overlaps checking



  // NUCLEO

  //

  // Forma 

  // 

  G4ThreeVector pos1 = G4ThreeVector(0, 0, 0);


  G4double dim_x = 20*cm;

  G4double dim_y = 20*cm;

  G4double dim_z = 200*cm;

  auto solidShape1 = new G4Box("Detector", dim_x / 2, dim_y / 2, dim_z / 2);


  //

  // material del nucleo

  //

  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_Al");

  

  // Rejilla

   G4double cut_sizeXY = 9.5 * cm;

   G4double cut_sizeZ = 200.1 * cm;

  auto solidCut = new G4Box("Cut", cut_sizeXY / 2, cut_sizeXY / 2, cut_sizeZ / 2);

  

   G4ThreeVector positions[4] = {

//        G4ThreeVector(-45 * cm, 45 * cm, 0 * cm), //1

//        G4ThreeVector(-35 * cm, 45 * cm, 0 * cm), //2

//      G4ThreeVector(-25 * cm, 45 * cm, 0 * cm), //3

//        G4ThreeVector(-15 * cm, 45 * cm, 0 * cm), //4

//        G4ThreeVector(-5 * cm, 45 * cm, 0 * cm), //5

//        G4ThreeVector(5* cm, 45 * cm, 0 * cm), //6

//        G4ThreeVector(15 * cm, 45 * cm, 0 * cm), //7

//        G4ThreeVector(25 * cm, 45 * cm, 0 * cm), //8

//        G4ThreeVector(35 * cm, 45 * cm, 0 * cm), //9

//        G4ThreeVector(45 * cm, 45* cm, 0 * cm), //10

//        G4ThreeVector(-45 * cm, 35 * cm, 0 * cm), //11

//        G4ThreeVector(-35 * cm, 35* cm, 0 * cm), //12

//        G4ThreeVector(-25 * cm, 35 * cm, 0 * cm), //13

//        G4ThreeVector(-15 * cm, 35 * cm, 0 * cm), //14

//        G4ThreeVector(-5 * cm, 35 * cm, 0 * cm), //15
        
//        G4ThreeVector(5 * cm, 35 * cm, 0 * cm), //16

//        G4ThreeVector(15 * cm, 35 * cm, 0 * cm), //17

//        G4ThreeVector(25 * cm, 35 * cm, 0 * cm), //18

//        G4ThreeVector(35 * cm, 35 * cm, 0 * cm), //19

 //       G4ThreeVector(45 * cm, 35 * cm, 0 * cm) //20 observacion con la coma

//        G4ThreeVector(-45 * cm, 25 * cm, 0 * cm), //21

//        G4ThreeVector(-35 * cm, 25 * cm, 0 * cm), //22

//        G4ThreeVector(-25 * cm, 25 * cm, 0 * cm), //23

//        G4ThreeVector(-15 * cm, 25* cm, 0 * cm), //24

//        G4ThreeVector(-5 * cm, 25 * cm, 0 * cm), //25

//        G4ThreeVector(5 * cm, 25 * cm, 0 * cm), //26

//        G4ThreeVector(15 * cm, 25 * cm, 0 * cm), //27

//        G4ThreeVector(25 * cm, 25 * cm, 0 * cm), //28

//        G4ThreeVector(35 * cm, 25 * cm, 0 * cm), //29

//        G4ThreeVector(45 * cm, 25 * cm, 0 * cm), //30
        
//        G4ThreeVector(-45 * cm, 15 * cm, 0 * cm), //31

//        G4ThreeVector(-35 * cm, 15 * cm, 0 * cm), //32

//        G4ThreeVector(-25 * cm, 15 * cm, 0 * cm), //33

//        G4ThreeVector(-15 * cm, 15 * cm, 0 * cm), //34

//        G4ThreeVector(-5 * cm, 15 * cm, 0 * cm), //35

//        G4ThreeVector(5 * cm, 15 * cm, 0 * cm), //36

//        G4ThreeVector(15 * cm, 15 * cm, 0 * cm), //37

//        G4ThreeVector(25 * cm, 15 * cm, 0 * cm), //38
  
//        G4ThreeVector(35 * cm, 15 * cm, 0 * cm), //39

//        G4ThreeVector(45 * cm, 15 * cm, 0 * cm), //40

//        G4ThreeVector(-45 * cm, 5 * cm, 0 * cm), //41

//        G4ThreeVector(-35 * cm, 5 * cm, 0 * cm), //42

//        G4ThreeVector(-25 * cm, 5 * cm, 0 * cm), //43

//        G4ThreeVector(-15 * cm, 5 * cm, 0 * cm), //44

        G4ThreeVector(-5 * cm, 5 * cm, 0 * cm), //45
        
        G4ThreeVector(5 * cm, 5 * cm, 0 * cm), //46

//        G4ThreeVector(15 * cm, 5 * cm, 0 * cm), //47

//        G4ThreeVector(25 * cm, 5 * cm, 0 * cm), //48

//        G4ThreeVector(35 * cm, 5 * cm, 0 * cm), //49

//        G4ThreeVector(45 * cm, 5 * cm, 0 * cm), //50*****
        
//        G4ThreeVector(-45 * cm, 5 * cm, 0 * cm), //1

//        G4ThreeVector(-35 * cm, 5 * cm, 0 * cm), //2

//        G4ThreeVector(-25 * cm, 5 * cm, 0 * cm), //3

//        G4ThreeVector(-15 * cm, 5 * cm, 0 * cm), //4

        G4ThreeVector(-5 * cm, -5 * cm, 0 * cm), //5

        G4ThreeVector(5* cm, -5 * cm, 0 * cm) //6

//        G4ThreeVector(15 * cm, 5 * cm, 0 * cm), //7

//        G4ThreeVector(25 * cm, 5 * cm, 0 * cm), //8

//        G4ThreeVector(35 * cm, 5 * cm, 0 * cm), //9

//        G4ThreeVector(45 * cm, 5* cm, 0 * cm), //10

//        G4ThreeVector(-45 * cm, 15 * cm, 0 * cm), //11

//        G4ThreeVector(-35 * cm, 15* cm, 0 * cm), //12

//        G4ThreeVector(-25 * cm, 15 * cm, 0 * cm), //13

//        G4ThreeVector(-15 * cm, 15 * cm, 0 * cm), //14

  //      G4ThreeVector(-5 * cm, 15 * cm, 0 * cm), //15
        
    //    G4ThreeVector(5 * cm, 15 * cm, 0 * cm), //16

   //     G4ThreeVector(15 * cm, 15 * cm, 0 * cm) //17

//        G4ThreeVector(25 * cm, 15 * cm, 0 * cm), //18

//        G4ThreeVector(35 * cm, 15 * cm, 0 * cm), //19

//        G4ThreeVector(45 * cm, 15 * cm, 0 * cm), //20

//        G4ThreeVector(-45 * cm, 25 * cm, 0 * cm), //21

//        G4ThreeVector(-35 * cm, 25 * cm, 0 * cm), //22

//        G4ThreeVector(-25 * cm, 25 * cm, 0 * cm), //23

//        G4ThreeVector(-15 * cm, 25* cm, 0 * cm), //24

//        G4ThreeVector(-5 * cm, 25 * cm, 0 * cm), //25

//        G4ThreeVector(5 * cm, 25 * cm, 0 * cm), //26

//        G4ThreeVector(15 * cm, 25 * cm, 0 * cm), //27

//        G4ThreeVector(25 * cm, 25 * cm, 0 * cm), //28

//        G4ThreeVector(35 * cm, 25 * cm, 0 * cm), //29

//        G4ThreeVector(45 * cm, 25 * cm, 0 * cm), //30
        
//        G4ThreeVector(-45 * cm, 35 * cm, 0 * cm), //31

//        G4ThreeVector(-35 * cm, 35 * cm, 0 * cm), //32

//        G4ThreeVector(-25 * cm, 35 * cm, 0 * cm), //33

//        G4ThreeVector(-15 * cm, 35 * cm, 0 * cm), //34

//        G4ThreeVector(-5 * cm, 35 * cm, 0 * cm), //35

//        G4ThreeVector(5 * cm, 35 * cm, 0 * cm), //36

//        G4ThreeVector(15 * cm, 35 * cm, 0 * cm), //37

//        G4ThreeVector(25 * cm, 35 * cm, 0 * cm), //38
  
//        G4ThreeVector(35 * cm, 35 * cm, 0 * cm), //39

//        G4ThreeVector(45 * cm, 35 * cm, 0 * cm), //40

//        G4ThreeVector(-45 * cm, 45 * cm, 0 * cm), //41

//        G4ThreeVector(-35 * cm, 45 * cm, 0 * cm), //42

//        G4ThreeVector(-25 * cm, 45 * cm, 0 * cm), //43

//        G4ThreeVector(-15 * cm, 45 * cm, 0 * cm), //44

//        G4ThreeVector(-5 * cm, 45 * cm, 0 * cm), //45
        
//        G4ThreeVector(5 * cm, 45 * cm, 0 * cm), //46

//        G4ThreeVector(15 * cm, 45 * cm, 0 * cm), //47

//        G4ThreeVector(25 * cm, 45 * cm, 0 * cm), //48

//        G4ThreeVector(35 * cm, 45 * cm, 0 * cm), //49

//        G4ThreeVector(45 * cm, 45 * cm, 0 * cm) //50*****
        
    };

    

    G4VSolid* subtractedDetector = solidShape1;



    for (int i = 0; i < 4; i++) {

        subtractedDetector = new G4SubtractionSolid("DetectorCut", subtractedDetector, solidCut, nullptr, positions[i]);

    }

  

  auto logicShape2 = new G4LogicalVolume(subtractedDetector,  // its solid

    shape2_mat,                                        // its material

    "Volumen2");                                         // its name



  new G4PVPlacement(nullptr,  // no rotation

    G4ThreeVector(),                     // at position

    logicShape2,              // its logical volume

    "Volumen2",                 // its name

    logicEnv,                 // its mother  volume

    false,                    // no boolean operation

    0,                        // copy number

    checkOverlaps);           // overlaps checking

    

    

    // ELEMENTO COMBUSTIBLE NORMAL ECN

    // Barra 
    
    G4Material* barra_mat = nist->FindOrBuildMaterial("G4_Al");
    
    G4double bx1 = 9*cm, by1 = 8*cm, bz1 = 77*cm;
    G4double bx2 = 8.75*cm, by2 = 7.75*cm, bz2 = 77.1*cm;
    
    auto solidBarra1 = new G4Box("barra1", bx1/2, by1/2, bz1/2);
    
    auto solidBarra2 = new G4Box("barra2", bx2/2, by2/2, bz2/2);
    
    //G4VSolid* subtractionbarra = solidBarra;
    
    G4SubtractionSolid* subtractionbarra = new G4SubtractionSolid("BarraCorte", solidBarra1, solidBarra2);
    
    auto logicBarraCorte = new G4LogicalVolume(subtractionbarra, barra_mat, "BarraCorte");
    
    
    new G4PVPlacement(nullptr,
    	positions[0],
    	logicBarraCorte,
    	"BarraCorte",
    	logicEnv,
    	false,
    	0,
    	checkOverlaps);
    	
    	
    // placas de elemento combustible normal
    
    //G4Material* placa_mat = nist->FindOrBuildMaterial("G4_U");
    
    // parametros de los materiales
    G4String name, symbol;
    G4double a, density;
    G4int iz, n;
    G4int ncomponents, natoms;
    G4double abundance;
    
    // creacion del elemento combustible
    G4Element* elSi = nist->FindOrBuildElement("Si");
    // relacion isotopica del uranio
    G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
    G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);
    G4Element* elU  = new G4Element(name="Uranio enrequecido", symbol="U", ncomponents=2);
    elU->AddIsotope(U5, abundance=80.*perCent);
    elU->AddIsotope(U8, abundance=20.*perCent);
    // molecula de siliciuro
    density = 4.8*g/cm3;
    G4Material* U3Si2 = new G4Material(name="siliciuro de uranio", density, ncomponents=2);
    U3Si2->AddElement(elU, natoms=3);
    U3Si2->AddElement(elSi, natoms=2);
    
    G4double px1 = 0.2*cm, py1 = 7.7*cm, pz1 = 68*cm;
    
    
    G4ThreeVector positions2[16] = {
    	G4ThreeVector(-0.25 * cm, 0 * cm, 0 * cm), //1

        G4ThreeVector(-0.75 * cm, 0 * cm, 0 * cm), //2

        G4ThreeVector(-1.25 * cm, 0 * cm, 0 * cm), //3

        G4ThreeVector(-1.75 * cm, 0 * cm, 0 * cm), //4

        G4ThreeVector(-2.25 * cm, 0 * cm, 0 * cm), //5

        G4ThreeVector(-2.75 * cm, 0 * cm, 0 * cm), //6

        G4ThreeVector(-3.25 * cm, 0 * cm, 0 * cm), //7

        G4ThreeVector(-3.75 * cm, 0 * cm, 0 * cm), //8*****

        G4ThreeVector(0.25 * cm, 0 * cm, 0 * cm), //1

        G4ThreeVector(0.75 * cm, 0 * cm, 0 * cm), //2

        G4ThreeVector(1.25 * cm, 0 * cm, 0 * cm), //3

        G4ThreeVector(1.75 * cm, 0 * cm, 0 * cm), //4

        G4ThreeVector(2.25 * cm, 0 * cm, 0 * cm), //5

        G4ThreeVector(2.75 * cm, 0 * cm, 0 * cm), //6

        G4ThreeVector(3.25 * cm, 0 * cm, 0 * cm), //7

        G4ThreeVector(3.75 * cm, 0 * cm, 0 * cm) //8
    
    };
    
    
    auto solidPlacas = new G4Box("placas", px1/2, py1/2, pz1/2);
    
    auto logicPlaca = new G4LogicalVolume(solidPlacas, U3Si2, "placas");
    
    for (int i = 0; i<16; i++){
    	new G4PVPlacement(nullptr,
	    	positions2[i],
	    	logicPlaca,
	    	"placas",
	    	logicBarraCorte,
	    	false,
	    	i,
	    	checkOverlaps);
    	}



  // ELEMENTO COMBUSTIBLE DE CONTROL ECC
  
  // Barra
  
    G4SubtractionSolid* subtractionbarra2 = new G4SubtractionSolid("BarraCorte2", solidBarra1, solidBarra2);
    
    auto logicBarraCorte2 = new G4LogicalVolume(subtractionbarra2, barra_mat, "BarraCorte2");
    
    new G4PVPlacement(nullptr,
    	positions[1],
    	logicBarraCorte2,
    	"BarraCorte2",
    	logicEnv,
    	false,
    	0,
    	checkOverlaps);
  
  // placas de elemento combustible de control
  
    G4ThreeVector positions3[12] = {
    	G4ThreeVector(-0.35 * cm, 0 * cm, 0 * cm), //1

        G4ThreeVector(-1.05 * cm, 0 * cm, 0 * cm), //2

        G4ThreeVector(-1.75 * cm, 0 * cm, 0 * cm), //3

        G4ThreeVector(-2.45 * cm, 0 * cm, 0 * cm), //4

        G4ThreeVector(-3.15 * cm, 0 * cm, 0 * cm), //5

        G4ThreeVector(-3.85 * cm, 0 * cm, 0 * cm), //6

        G4ThreeVector(0.35 * cm, 0 * cm, 0 * cm), //7

        G4ThreeVector(1.05 * cm, 0 * cm, 0 * cm), //8*****

        G4ThreeVector(1.75 * cm, 0 * cm, 0 * cm), //1

        G4ThreeVector(2.45 * cm, 0 * cm, 0 * cm), //2

        G4ThreeVector(3.15 * cm, 0 * cm, 0 * cm), //3

        G4ThreeVector(3.85 * cm, 0 * cm, 0 * cm)
    
    };
   
    auto solidPlacas2 = new G4Box("placas2", px1/2, py1/2, pz1/2);
    
    auto logicPlaca2 = new G4LogicalVolume(solidPlacas2, U3Si2, "placas2");
    
    for (int i = 0; i<12; i++){
    	new G4PVPlacement(nullptr,
	    	positions3[i],
	    	logicPlaca2,
	    	"placas2",
	    	logicBarraCorte2,
	    	false,
	    	i,
	    	checkOverlaps);
    	}
  
  
  // REFLECTOR DE BERILIO
  
  G4Material* berilio_mat = nist->FindOrBuildMaterial("G4_Be");
  
  G4double rminb = 0*cm, rmaxb = 4*cm;
  G4double zb = 68*cm;
  G4double aib = 0*deg, afb = 360*deg;
   
  auto solidBerilio = new G4Tubs("berilio", rminb, rmaxb, zb, aib, afb);
  
  auto logicBerilio = new G4LogicalVolume(solidBerilio, berilio_mat, "berilio");
  
  new G4PVPlacement(nullptr,
  	positions[2],
  	logicBerilio,
  	"berilio",
  	logicEnv,
  	false,
  	0,
  	checkOverlaps);
  	
  	
  // DETECTOR DE HE/3
  
  G4Material* detector_mat = nist->FindOrBuildMaterial("G4_Al");
  
  G4double rmind = 0*cm, rmaxd = 3*cm;
  G4double zd = 120*cm;
  G4double aid = 0*deg, afd = 360*deg;
   
  auto solidDetec = new G4Tubs("detector", rmind, rmaxd, zd, aid, afd);
  
  auto logicDetec = new G4LogicalVolume(solidDetec, detector_mat, "detector");
  
  new G4PVPlacement(nullptr,
  	positions[3],
  	logicDetec,
  	"detector",
  	logicEnv,
  	false,
  	0,
  	checkOverlaps);
  
  

  // Set Shape as scoring volume

  //

  fScoringVolume = logicDetec;



  //

  //always return the physical World

  //

  return physWorld;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



}

