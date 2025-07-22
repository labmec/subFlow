// Gmsh project created on Fri Jul  5 14:35:05 2024
SetFactory("OpenCASCADE");
rin = 30.15; //sandscreen tube radius
rout = 30.15+345; //external radius
h = 1000; //height
//+
Point(1) = {rin, 0, 0, 1.0};
//+
Point(2) = {rin+310.5, 0, 0, 1.0};
//+
Point(3) = {rout, 0, 0, 1.0};
//+
Point(4) = {rout, h, 0, 1.0};
//+
Point(5) = {rin+310.5, h, 0, 1.0};
//+
Point(6) = {rin, h, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {1,2,3,4,5,6};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bottom", 2) = {1,2};
//+
Physical Curve("diffusor", 3) = {3};
//+
Physical Curve("gap", 4) = {4};
//+
Physical Curve("lid", 5) = {5};
//+
Physical Curve("sandscreen", 6) = {6};
//+
Physical Surface("dom", 1) = {1};
//+
Transfinite Surface {1} = {1,3,4,6};
//+
Transfinite Curve {1,5} = 19 Using Progression 1;
//+
Transfinite Curve {2,4} = 3 Using Progression 1;
//+
Transfinite Curve {3,6} = 21 Using Progression 1;
//+
Recombine Surface {1};
