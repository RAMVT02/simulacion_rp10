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

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

namespace B1

{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()

{

  // Get nist material manager

  G4NistManager* nist = G4NistManager::Instance();
  
  //
  
  // materials
  
  //
  
  G4Material* aire = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Material* agua = nist->FindOrBuildMaterial("G4_WATER");
  
  G4Material* mat_Al = nist->FindOrBuildMaterial("G4_Al");
  
  G4Element* Al = nist->FindOrBuildElement("Al");
  
  G4Element* Mg = nist->FindOrBuildElement("Mg");
 
  G4Element* B = nist->FindOrBuildElement("B");

  G4Element* Ni = nist->FindOrBuildElement("Ni");
  
  G4Element* Be = nist->FindOrBuildElement("Be");

  G4Element* C = nist->FindOrBuildElement("C");
  
  G4Element* Si = nist->FindOrBuildElement("Si");

  G4Element* Ag = nist->FindOrBuildElement("Ag");

  G4Element* In = nist->FindOrBuildElement("In");

  G4Element* Cd = nist->FindOrBuildElement("Cd");
  
  G4Material* acero_inoxidable = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    
    // cronstruccion de materiales
    G4double abundance;
    G4int natoms;
    G4double densityU3Si2 = 4.8*g/cm3, densityAlMg = 2.67*g/cm3;
    G4double densityC = 2.1*g/cm3, densityBe = 1.8*g/cm3;
    
    // combustible
    // relacion isotopica del uranio
    G4Isotope* U5 = new G4Isotope("U235", 92, 235, 235.01*g/mole);
    G4Isotope* U8 = new G4Isotope("U238", 92, 238, 238.03*g/mole);
    G4Element* elU  = new G4Element("Uranio enrequecido", "U", 2);
    elU->AddIsotope(U5, abundance=20.*perCent);
    elU->AddIsotope(U8, abundance=80.*perCent);
    // molecula de siliciuro
    G4Material* U3Si2 = new G4Material("siliciuro de uranio", densityU3Si2, 2);
    U3Si2->AddElement(elU, natoms=3);
    U3Si2->AddElement(Si, natoms=2);
    
    // barra contenedora del combustible
    G4Material* AlMg = new G4Material("siliciuro de uranio", densityAlMg, 2);
    AlMg->AddElement(Al, natoms=1);
    AlMg->AddElement(Mg, natoms=1);

    // berilio grado nuclear
    G4double niquelMassFraction = 15.0e-6;
    G4double berilioMassFraction = 1.0 - niquelMassFraction;
    G4Material* berilio_nuclear = new G4Material("berilio_alta_puerza", densityBe, 2);
    berilio_nuclear->AddElement(Be, berilioMassFraction);
    berilio_nuclear->AddElement(Ni, niquelMassFraction);

    // grafito grado nuclear
    G4double boronMassFraction = 15.0e-6;
    G4double carbonMassFraction = 1.0 - boronMassFraction;
    G4Material* grafito_nuclear = new G4Material("grafito_alta_puerza", densityC, 2);
    grafito_nuclear->AddElement(C, carbonMassFraction);
    grafito_nuclear->AddElement(B, boronMassFraction);

  // aleacion AgInCd
  G4double densityAgInCd = 10.17 * g/cm3;
  G4double Ag_fraction = 0.80; 
  G4double In_fraction = 0.15;
  G4double Cd_fraction = 0.05; 
  G4Material* AgInCd = new G4Material("AgInCd_aleacion", densityAgInCd, 3);
    AgInCd->AddElement(Ag, Ag_fraction);
    AgInCd->AddElement(In, In_fraction);
    AgInCd->AddElement(Cd, Cd_fraction);

 
  // Option to switch on/off checking of volumes overlaps

  G4bool checkOverlaps = true;


  //

  // Mundo

  //

  G4double rmin = 0*cm, rmax = 202*cm;

  G4double z = 1120*cm;
  
  G4double ai = 0*deg, af = 360*deg;

  auto solidWorld = new G4Tubs("World",                           // its name

    rmin, rmax, 0.5*z, ai, af);  // its size



  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid

    aire,                                       // its material

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

  // dimensiones reales - altura 1100cm y radio 200cm
  
  G4double r1max = 50*cm;

  G4double z1 = 200*cm;

  auto solidEnv = new G4Tubs("Envelope",                    // its name

    rmin, r1max, 0.5*z1, ai, af);  // its size



  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid

    agua,                                     // its material

    "Envelope");                                 // its name



  new G4PVPlacement(nullptr,  // no rotation

    G4ThreeVector(),          // at (0,0,0)

    logicEnv,                 // its logical volume

    "Envelope",               // its name

    logicWorld,               // its mother  volume

    false,                    // no boolean operation

    0,                        // copy number

    checkOverlaps);           // overlaps checking


  
  //

  // NUCLEO

  //


  // Forma 

  G4ThreeVector pos = G4ThreeVector(0, 0, 0);

  G4double dim_x = 20*cm;

  G4double dim_y = 20*cm;

  G4double dim_z = 200*cm;

  auto solidShape1 = new G4Box("Detector", dim_x / 2, dim_y / 2, dim_z / 2);

  // Rejilla

   G4double cut_sizeXY = 9.75*cm;

   G4double cut_sizeZ = 200.01*cm;

  auto solidCut = new G4Box("Cut", cut_sizeXY / 2, cut_sizeXY / 2, cut_sizeZ / 2);

   G4ThreeVector pos0[4] = {

        G4ThreeVector(-5 * cm, 5 * cm, 0 * cm), //1 ECN
        
        G4ThreeVector(5 * cm, 5 * cm, 0 * cm), //2 ECC

        G4ThreeVector(-5 * cm, -5 * cm, 0 * cm), //3 BERILIO

        G4ThreeVector(5* cm, -5 * cm, 0 * cm) //4 GRAFITO

    };

    G4VSolid* subtractedDetector = solidShape1;

    for (int i = 0; i < 4; i++) {

        subtractedDetector = new G4SubtractionSolid("DetectorCut", subtractedDetector, solidCut, nullptr, pos0[i]);

    }

  auto logicShape2 = new G4LogicalVolume(subtractedDetector,  // its solid

    mat_Al,                                        // its material

    "Volumen2");                                         // its name



  new G4PVPlacement(nullptr,  // no rotation

    G4ThreeVector(),                     // at position

    logicShape2,              // its logical volume

    "Volumen2",                 // its name

    logicEnv,                 // its mother  volume

    false,                    // no boolean operation

    0,                        // copy number

    checkOverlaps);           // overlaps checking
    


    //

    // ELEMENTO COMBUSTIBLE NORMAL ECN

    //
    

    // barra contenedora 
    
    auto box1 = new G4Box("box1", 0.5*7.7*cm, 0.5*8.2*cm, 0.5*67.75*cm);

    auto box2 = new G4Box("box2", 0.5*6.72*cm, 0.5*8.124*cm, 0.5*67.75*cm);
    
   G4SubtractionSolid* subtractionbox1y2 = new G4SubtractionSolid("subtractionbox1y2", box1, box2);
    
    auto logicbox1y2 = new G4LogicalVolume(subtractionbox1y2, mat_Al, "subtractionbox1y2");
    
    G4VisAttributes* logicbox1y2VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // verde
    logicbox1y2VisAtt->SetVisibility(true);
    logicbox1y2VisAtt->SetForceWireframe(true); 
    logicbox1y2->SetVisAttributes(logicbox1y2VisAtt);
    
    new G4PVPlacement(nullptr,
    	pos0[0],
    	logicbox1y2,
    	"subtractionbox1y2",
    	logicEnv,
    	false,
    	0,
    	checkOverlaps);
    	
    // placas de combustible
    
    auto bar1 = new G4Box("bar1", 0.5*6.72*cm, 0.5*0.176*cm, 0.5*67.75*cm); // 14 ECN
    
    auto bar2 = new G4Box("bar2", 0.5*6.72*cm, 0.5*0.19*cm, 0.5*67.75*cm); // 2 ECN
    	
    auto cort1y2 = new G4Box("barcort1", 0.5*6.275*cm, 0.5*0.1*cm, 0.5*67.75*cm);
    
    G4SubtractionSolid* subtractionbar1cort1y2 = new G4SubtractionSolid("subtractionbar1cort1y2", bar1, cort1y2);
    
    G4SubtractionSolid* subtractionbar2cort1y2 = new G4SubtractionSolid("subtractionbar2cort1y2", bar2, cort1y2);
    
    auto logicbar1 = new G4LogicalVolume( subtractionbar1cort1y2, AlMg, "subtractionbar1cort1y2");
    
    auto logicbar2 = new G4LogicalVolume( subtractionbar2cort1y2, AlMg, "subtractionbar1cort1y2");
    
    G4ThreeVector pos1[14] = {
    
      G4ThreeVector(0*cm, 3.289*cm, 0*cm),
      
      G4ThreeVector(0*cm, 2.783*cm, 0*cm),
      
      G4ThreeVector(0*cm, 2.277*cm, 0*cm),
      
      G4ThreeVector(0*cm, 1.771*cm, 0*cm),
      
      G4ThreeVector(0*cm, 1.265*cm, 0*cm),
      
      G4ThreeVector(0*cm, 0.759*cm, 0*cm),
      
      G4ThreeVector(0*cm, 0.253*cm, 0*cm),
      
      G4ThreeVector(0*cm, -0.253*cm, 0*cm),
      
      G4ThreeVector(0*cm, -0.759*cm, 0*cm),
      
      G4ThreeVector(0*cm, -1.265*cm, 0*cm),
      
      G4ThreeVector(0*cm, -1.771*cm, 0*cm),
      
      G4ThreeVector(0*cm, -2.277*cm, 0*cm),
      
      G4ThreeVector(0*cm, -2.783*cm, 0*cm),
      
      G4ThreeVector(0*cm, -3.289*cm, 0*cm)
    };
    
    for (int i = 0; i<14; i++){
      new G4PVPlacement(nullptr,
      pos1[i],
      logicbar1,
      "bar1cort1y2",
      logicbox1y2,
      false,
      i,
      checkOverlaps
      );
    }
    
    G4ThreeVector pos2[2] = {
      
      G4ThreeVector(0*cm, 3.802*cm, 0*cm),
      
      G4ThreeVector(0*cm, -3.802*cm, 0*cm)
    };
    
    for (int i = 0; i<2; i++){
      new G4PVPlacement(nullptr,
      pos2[i],
      logicbar2,
      "bar1cort1y2",
      logicbox1y2,
      false,
      i,
      checkOverlaps
      );
    }
    
    // combustible 
    
    auto comb = new G4Box("combustible", 0.5*6.274*cm, 0.5*0.99*cm, 0.5*67.75*cm);
    
    auto logiccomb = new G4LogicalVolume(comb, U3Si2, "combustible");
    
    for (int i = 0; i<14; i++){
      new G4PVPlacement(nullptr,
      	pos1[i],
      	logiccomb,
      	"combustible1",
      	logicbox1y2,
      	false,
      	i,
      	checkOverlaps);
    }
    
    for (int i = 0; i<2; i++){
      new G4PVPlacement(nullptr,
      	pos2[i],
      	logiccomb,
      	"combustible2",
      	logicbox1y2,
      	false,
      	i,
      	checkOverlaps);
    }
    
    //
    
    // ELEMENTO COMBUSTIBLE DE CONTROL ECC
    
    //
    
    
    // barra contenedora 

    auto box3 = new G4Box("box3", 0.5*7.13*cm, 0.5*8.2*cm, 0.5*150*cm);

    auto box4 = new G4Box("box4", 0.5*6.15*cm, 0.5*8.201*cm, 0.5*150*cm);
    
   G4SubtractionSolid* subtractionbox3y4 = new G4SubtractionSolid("subtractionbox3y4", box3, box4);
    
    auto logicbox3y4 = new G4LogicalVolume(subtractionbox3y4, mat_Al, "subtractionbox3y4");
    
    G4VisAttributes* logicbox3y4VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // rojo
    logicbox3y4VisAtt->SetVisibility(true);
    logicbox3y4VisAtt->SetForceWireframe(true); 
    logicbox3y4->SetVisAttributes(logicbox3y4VisAtt);
    
    new G4PVPlacement(nullptr,
    	pos0[1],
    	logicbox3y4,
    	"subtractionbox3y4",
    	logicEnv,
    	false,
    	0,
    	checkOverlaps);
    
    
    // placas de combustible 

    auto bar3 = new G4Box("bar3", 0.5*6.15*cm, 0.5*0.127*cm, 0.5*150*cm); // 4 ECN

    auto bar4 = new G4Box("bar2", 0.5*6.575*cm, 0.5*0.44*cm, 0.5*150*cm); // 2 ECN
    	
    auto cort4 = new G4Box("barcort4", 0.5*6.275*cm, 0.5*0.29*cm, 0.5*150*cm);
    
    G4SubtractionSolid* subtractionbar4cort4 = new G4SubtractionSolid("subtractionbar4cort4", bar4, cort4);
    
    auto logicbar3 = new G4LogicalVolume( bar3, AlMg, "placa_combustible");
    
    auto logicbar4cort4 = new G4LogicalVolume( subtractionbar4cort4, AlMg, "subtractionbar4cort4");

    G4ThreeVector pos3[12] = {
      
      G4ThreeVector(0*cm, 2.783*cm, 0*cm),
      
      G4ThreeVector(0*cm, 2.277*cm, 0*cm),
      
      G4ThreeVector(0*cm, 1.771*cm, 0*cm),
      
      G4ThreeVector(0*cm, 1.265*cm, 0*cm),
      
      G4ThreeVector(0*cm, 0.759*cm, 0*cm),
      
      G4ThreeVector(0*cm, 0.253*cm, 0*cm),
      
      G4ThreeVector(0*cm, -0.253*cm, 0*cm),
      
      G4ThreeVector(0*cm, -0.759*cm, 0*cm),
      
      G4ThreeVector(0*cm, -1.265*cm, 0*cm),
      
      G4ThreeVector(0*cm, -1.771*cm, 0*cm),
      
      G4ThreeVector(0*cm, -2.277*cm, 0*cm),
      
      G4ThreeVector(0*cm, -2.783*cm, 0*cm)
    };

    for (int i = 0; i<12; i++){
      new G4PVPlacement(nullptr,
      pos3[i],
      logicbar1,
      "bar1ecc",
      logicbox3y4,
      false,
      i,
      checkOverlaps
      );
    }

    G4ThreeVector pos4[4] = {
      
      G4ThreeVector(0*cm, 3.9115*cm, 0*cm),

      G4ThreeVector(0*cm, 3.1745*cm, 0*cm),

      G4ThreeVector(0*cm, -3.1745*cm, 0*cm),
      
      G4ThreeVector(0*cm, -3.9115*cm, 0*cm)
    };

    for (int i = 0; i<4; i++){
      new G4PVPlacement(nullptr,
      pos4[i],
      logicbar3,
      "bar3",
      logicbox3y4,
      false,
      i,
      checkOverlaps
      );
    }

    G4ThreeVector pos5[2] = {
      
      G4ThreeVector(0*cm, 3.543*cm, 0*cm),
      
      G4ThreeVector(0*cm, -3.543*cm, 0*cm)
    };

    for (int i = 0; i<2; i++){
      new G4PVPlacement(nullptr,
      pos5[i],
      logicbar4cort4,
      "bar4cort4",
      logicbox3y4,
      false,
      i,
      checkOverlaps
      );
    }
   
   // aleacion 

   auto control = new G4Box("aleacion", 0.5*6.274*cm, 0.5*0.28*cm, 0.5*150*cm);

   auto logiccontrol = new G4LogicalVolume(control, AgInCd, "aleacion");

  for (int i = 0; i<2; i++){
	  new G4PVPlacement(nullptr,
		pos5[i],
		logiccontrol,
		"aleacion",
		logicbox3y4,
		false,
		i,
		checkOverlaps);
  }

   // combustible

   for (int i = 0; i<12; i++){
      new G4PVPlacement(nullptr,
      	pos3[i],
      	logiccomb,
      	"combustible1ecc",
      	logicbox3y4,
      	false,
      	i,
      	checkOverlaps);
    }



 //

  // REFLECTOR DE BERILIO

 //


  // barra contenedora
  
  auto barcort1 = new G4Box("barcort1", 0.5*9.024*cm, 0.5*7.6*cm, 0.5*67.4*cm);

  auto barcort2 = new G4Box("barcort2", 0.5*8.124*cm, 0.5*6.7*cm, 0.5*67.4*cm);

  G4SubtractionSolid* subtractionbarcort1y2 = new G4SubtractionSolid("subtractionbarcort1y2", barcort1, barcort2);

  auto logicbarcortbe = new G4LogicalVolume(subtractionbarcort1y2, AlMg, "barcort1y2");
  
    G4VisAttributes* logicbarcortbeVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)); // magenta
    logicbarcortbeVisAtt->SetVisibility(true);
    logicbarcortbeVisAtt->SetForceWireframe(true); 
    logicbarcortbe->SetVisAttributes(logicbarcortbeVisAtt);

  new G4PVPlacement(nullptr,
	pos0[2],
	logicbarcortbe,
	"barcort1y2",
	logicEnv,
	false,
	0,
	checkOverlaps);

  // berilio

  auto solidbe = new G4Box("berilio", 0.5*8.123*cm, 0.5*6.69*cm, 0.5*67.4*cm);

  auto logicbe = new G4LogicalVolume(solidbe, berilio_nuclear, "berilio");

  new G4PVPlacement(nullptr,
	pos0[2],
	logicbe,
	"berilio",
	logicEnv,
	false,
	0,
	checkOverlaps);



   //

  // REFLECTOR DE GRAFITO

  //


  // barra contenedora
  
  auto barcort3 = new G4Box("barcort3", 0.5*8.474*cm, 0.5*8.2*cm, 0.5*70.5*cm);

  auto barcort4 = new G4Box("barcort4", 0.5*7.874*cm, 0.5*7.6*cm, 0.5*70.5*cm);

  G4SubtractionSolid* subtractionbarcort3y4 = new G4SubtractionSolid("subtractionbarcort3y4", barcort3, barcort4);

  auto logicbarcortgrafito = new G4LogicalVolume(subtractionbarcort3y4, AlMg, "barcort3y4");
  
    G4VisAttributes* logicbarcortgrafitoVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0)); // amarillo
    logicbarcortgrafitoVisAtt->SetVisibility(true);
    logicbarcortgrafitoVisAtt->SetForceWireframe(true); 
    logicbarcortgrafito->SetVisAttributes(logicbarcortgrafitoVisAtt);

  new G4PVPlacement(nullptr,
	pos0[3],
	logicbarcortgrafito,
	"barcort3y4",
	logicEnv,
	false,
	0,
	checkOverlaps);


  // grafito

  auto solidgrafito = new G4Box("grafito", 0.5*7.873*cm, 0.5*7.59*cm, 0.5*70.5*cm);

  auto logicgrafito = new G4LogicalVolume(solidgrafito, grafito_nuclear, "grafito");

  new G4PVPlacement(nullptr,
	pos0[3],
	logicgrafito,
	"grafito",
	logicEnv,
	false,
	0,
	checkOverlaps);



  // Set Shape as scoring volume

  //

  fScoringVolume = logicgrafito;



  //

  //always return the physical World

  //

  return physWorld;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



}


