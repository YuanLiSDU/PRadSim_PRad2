Using tool package CADMesh, it can read STL files and import in Geant4.
1. Donwlod CADMesh, git clone https://github.com/christopherpoole/CADMesh.git
2. copy ~/CADMesh/CADMesh.hh to ~/PRadSim/include/
3. In ~/PRadSim/src/DetectorConstruction.cc, write #include "CADMesh.hh"

Get Model.stl file from CAD step file.
1. Using FreeCAD to open the step file. 
2. Unfold the toolbox, choose "Mesh Design".
3. Click the model you want in the small window on the left.
4. Unfold "Meshes" on the top, choose "creat mesh from shape".
5. Meshing options: standard, write surface deviation and angular deviation, then click "OK".
6. Right-click the Model(meshed), choose "export mesh",and "ASCLL STL(.stl)", then we get the file Model(mesh).stl
8. Copy the file to ~/PRadSim/

Code to import in Geant4 DetectorConstruction.cc:
	// CADMesh :: STL //
	auto Model = CADMesh::TessellatedMesh::FromSTL("./Model(mesh).stl");
	Model->SetScale(1);
	Model->SetOffset(0,0,0);
        G4LogicalVolume *logicalModel = new G4LogicalVolume(Model->GetSolid(), Matirial,"ModelNameLV");
        new G4PVPlacement(0, G4ThreeVector(0,0,0), logicalModel, "ModelName", logicWorld, false, 0);

make 
./PRadSim
        
