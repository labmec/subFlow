// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
L = 1;
h = 1;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {L, 0, 0, 1.0};
//+
Point(3) = {L, h, 0, 1.0};
//+
Point(4) = {0, h, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {1,2,3,4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom", 2) = {1};
//+
Physical Curve("right", 3) = {2};
//+
Physical Curve("top", 4) = {3};
//+
Physical Curve("left", 5) = {4};
//+
Physical Surface("dom", 1) = {1};
//+
Transfinite Surface {1} = {1,2,3,4};
//+
Transfinite Curve {1,3} = 21 Using Progression 1;
//+
Transfinite Curve {2,4} = 21 Using Progression 1;
//+
Recombine Surface {1};
