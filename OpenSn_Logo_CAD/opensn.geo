SetFactory("OpenCASCADE");
Merge "opensn.step";
Coherence;
BooleanDifference(100) = { Volume{5}; Delete; } { Volume{9}; };
BooleanDifference(101) = { Volume{3}; Delete; } { Volume{10}; };
BooleanDifference(102) = { Volume{8}; Delete; } { Volume{100, 101, 2, 4, 6, 7}; };
Physical Volume("0") = {2, 4, 6, 7, 100, 101};
Physical Volume("1") = {102, 9, 10};
Mesh.CompoundClassify = 1;
Mesh.CharacteristicLengthFactor = 0.2;  // Adjust for finer mesh
