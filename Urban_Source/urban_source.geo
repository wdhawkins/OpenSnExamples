SetFactory("OpenCASCADE");
Merge "city.step";

Coherence;
BooleanDifference(100) = { Volume{3}; Delete; } { Volume{1, 2}; };
Physical Volume("0") = {1};
Physical Volume("1") = {2};
Physical Volume("2") = {100};
Mesh.CompoundClassify = 1;
Mesh.CharacteristicLengthFactor = 0.2;  // Adjust for finer mesh
