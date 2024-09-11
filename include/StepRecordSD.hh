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
// StepRecordSD.hh
// Developer : Chao Gu
// History:
//   May 2018, C. Gu, Add for beam energy loss study.
//
// WARNING: do not use it for large density materials!
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef StepRecordSD_h
#define StepRecordSD_h 1

#include "G4VSensitiveDetector.hh"

#include "StandardHit.hh"

#include "G4String.hh"

#include <vector>

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;
class TTree;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StepRecordSD : public G4VSensitiveDetector
{
public:
    StepRecordSD(G4String name, G4String abbrev);
    virtual ~StepRecordSD();

    virtual void Initialize(G4HCofThisEvent *);
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *);

protected:
    virtual void Register(TTree *);
    virtual void Clear();

    G4int fID;
    G4String fAbbrev;

    bool fRegistered;

    int fN;
    std::vector<int> fPID; // Particle ID
    std::vector<int> fTID; // Track ID
    std::vector<int> fPTID; // Parent Track ID
    std::vector<double> fX;
    std::vector<double> fY;
    std::vector<double> fZ;
    std::vector<double> fMomentum;
    std::vector<double> fTheta;
    std::vector<double> fPhi;
    std::vector<int> fProcessID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
