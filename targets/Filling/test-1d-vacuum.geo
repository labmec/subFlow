// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
rin = 0; //sandscreen tube radius
rout = 1; //external radius
h = 1; //height
//+
Point(1) = {rin, 0, 0, 1.0};
//+
Point(2) = {rout, 0, 0, 1.0};
//+
Point(3) = {rout, h, 0, 1.0};
//+
Point(4) = {rin, h, 0, 1.0};
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
Physical Curve("diffusor", 3) = {2};
//+
Physical Curve("lid", 4) = {3};
//+
Physical Curve("sandscreen", 5) = {4};
//+
Physical Surface("dom", 1) = {1};
//+
Transfinite Surface {1} = {1,2,3,4};
//+
Transfinite Curve {1,3} = 2 Using Progression 1;
//+
Transfinite Curve {2,4} = 3 Using Progression 1;
//+
Recombine Surface {1};
