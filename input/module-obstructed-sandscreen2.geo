// Axisymmetric module, 0.8 m high, split into three surfaces (0.1 m, 0.52 m, 0.18 m)
SetFactory("OpenCASCADE");
//+
Point(1) = {0.0325, 0, 0, 1.0};
//+
Point(2) = {0.335, 0, 0, 1.0};
//+
Point(3) = {0.335, 0.1, 0, 1.0};
//+
Point(4) = {0.0325, 0.1, 0, 1.0};
//+
Point(5) = {0.335, 0.62, 0, 1.0};
//+
Point(6) = {0.0325, 0.62, 0, 1.0};
//+
Point(7) = {0.335, 0.8, 0, 1.0};
//+
Point(8) = {0.0325, 0.8, 0, 1.0};
//+
// Horizontal lines
Line(1) = {1, 2};
//+
Line(2) = {4, 3};
//+
Line(3) = {6, 5};
//+
Line(4) = {8, 7};
//+
// Right (diffuser) vertical lines
Line(5) = {2, 3};
//+
Line(6) = {3, 5};
//+
Line(7) = {5, 7};
//+
// Left (sandscreen side) vertical lines
Line(8) = {4, 1};
//+
Line(9) = {6, 4};
//+
Line(10) = {8, 6};
//+
Curve Loop(1) = {1, 5, -2, 8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 6, -3, 9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 7, -4, 10};
//+
Plane Surface(3) = {3};
//+
Physical Curve("bottom", 2) = {1};
//+
Physical Curve("diffuser", 3) = {5, 6, 7};
//+
Physical Curve("lid", 4) = {4};
//+
Physical Curve("sandscreen", 5) = {9};
//+
Physical Curve("noflow", 6) = {8, 10};
//+
Physical Surface("dom0", 1) = {1};
//+
Physical Surface("dom1", 7) = {2};
//+
Physical Surface("dom2", 8) = {3};
//+
Transfinite Surface {1, 2, 3};
//+
// Vertical lines: 0.01 m elements
Transfinite Curve {5, 8} = 11 Using Progression 1;
//+
Transfinite Curve {6, 9} = 53 Using Progression 1;
//+
Transfinite Curve {7, 10} = 19 Using Progression 1;
//+
Transfinite Curve {1, 2, 3, 4} = 41 Using Progression 1;
//+
Recombine Surface {1, 2, 3};
