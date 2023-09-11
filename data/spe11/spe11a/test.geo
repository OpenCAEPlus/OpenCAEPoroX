// SPDX-FileCopyrightText: 2023 Dennis Gläser <dennis.glaeser@iws.uni-stuttgart.de>
//
// SPDX-License-Identifier: MIT
DefineConstant[ refinement_factor = 1.0 ];

cl__1 = 0.1*refinement_factor;
cl__2 = 0.1*refinement_factor;
cl__3 = 0.1*refinement_factor;
cl__4 = 0.1*refinement_factor;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {1, 0, 0, cl__1};
Point(3) = {1, 1, 0, cl__1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

Physical Curve("Bottom_Boundary", 5) = {1};
Physical Curve("Right_Boundary", 6) = {2};
Physical Curve("Top_Boundary", 7) = {3};


// Facies 1
Curve Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};
Physical Surface("Facies 1", 1) = {1};