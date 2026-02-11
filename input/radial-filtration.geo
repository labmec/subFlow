// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0.0325, 0, 0, 1.0};
//+
Point(2) = {0.335, 0, 0, 1.0};
//+
Point(3) = {0.335, 0.87, 0, 1.0};
//+
Point(4) = {0.0325, 0.87, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom", 2) = {1};
//+
Physical Curve("diffuser", 3) = {2};
//+
Physical Curve("lid", 4) = {3};
//+
Physical Curve("sandscreen", 5) = {4};
//+
Physical Surface("dom", 1) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {4, 2} = 11 Using Progression 1;
//+
Transfinite Curve {3, 1} = 21 Using Progression 1;
//+
Recombine Surface {1};
