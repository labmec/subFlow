// Gmsh project: 3D Hollowed Cylinder with structured hex mesh in cylindrical coordinates
SetFactory("OpenCASCADE");

// Parameters
r_inner = 0.0325;    // Inner radius (m)
r_outer = 0.335;     // Outer radius (m)
height = 0.87;       // Height (m)

// Number of elements along each direction
n_radial = 10;       // Radial direction elements
n_tangential = 4;    // Tangential direction elements  
n_height = 10;       // Height direction elements

// Create center points for bottom and top
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, height, 0, 1.0};

// Inner circle (bottom) - create 4 points at every 90 degrees
Point(3) = {r_inner, 0, 0, 1.0};
Point(4) = {0, 0, r_inner, 1.0};
Point(5) = {-r_inner, 0, 0, 1.0};
Point(6) = {0, 0, -r_inner, 1.0};

// Outer circle (bottom) - create 4 points at every 90 degrees
Point(7) = {r_outer, 0, 0, 1.0};
Point(8) = {0, 0, r_outer, 1.0};
Point(9) = {-r_outer, 0, 0, 1.0};
Point(10) = {0, 0, -r_outer, 1.0};

// Inner circle (top) - create 4 points at every 90 degrees
Point(11) = {r_inner, height, 0, 1.0};
Point(12) = {0, height, r_inner, 1.0};
Point(13) = {-r_inner, height, 0, 1.0};
Point(14) = {0, height, -r_inner, 1.0};

// Outer circle (top) - create 4 points at every 90 degrees
Point(15) = {r_outer, height, 0, 1.0};
Point(16) = {0, height, r_outer, 1.0};
Point(17) = {-r_outer, height, 0, 1.0};
Point(18) = {0, height, -r_outer, 1.0};

// Inner circle (bottom)
Circle(1) = {3, 1, 4};
Circle(2) = {4, 1, 5};
Circle(3) = {5, 1, 6};
Circle(4) = {6, 1, 3};

// Outer circle (bottom)
Circle(5) = {7, 1, 8};
Circle(6) = {8, 1, 9};
Circle(7) = {9, 1, 10};
Circle(8) = {10, 1, 7};

// Inner circle (top)
Circle(9) = {11, 2, 12};
Circle(10) = {12, 2, 13};
Circle(11) = {13, 2, 14};
Circle(12) = {14, 2, 11};

// Outer circle (top)
Circle(13) = {15, 2, 16};
Circle(14) = {16, 2, 17};
Circle(15) = {17, 2, 18};
Circle(16) = {18, 2, 15};

// Inner cylinder vertical lines
Line(17) = {3, 11};
Line(18) = {4, 12};
Line(19) = {5, 13};
Line(20) = {6, 14};

// Outer cylinder vertical lines
Line(21) = {7, 15};
Line(22) = {8, 16};
Line(23) = {9, 17};
Line(24) = {10, 18};

// Create radial lines connecting the inner to the outer circles
// Radial lines (bottom)
Line(25) = {3, 7};
Line(26) = {4, 8};
Line(27) = {5, 9};
Line(28) = {6, 10};

// Radial lines (top)
Line(29) = {11, 15};
Line(30) = {12, 16};
Line(31) = {13, 17};
Line(32) = {14, 18};

//Create the circular surfaces
// Bottom surfaces
Curve Loop(1) = {1, 26, -5, -25};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 27, -6, -26};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 28, -7, -27};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 25, -8, -28};
Plane Surface(4) = {4};

// Top surfaces
Curve Loop(5) = {9, 30, -13, -29};
Plane Surface(5) = {5};
Curve Loop(6) = {10, 31, -14, -30};
Plane Surface(6) = {6};
Curve Loop(7) = {11, 32, -15, -31};
Plane Surface(7) = {7};
Curve Loop(8) = {12, 29, -16, -32};
Plane Surface(8) = {8};

// Inner cylindrical surfaces
Curve Loop(9) = {1, 18, -9, -17};
Surface(9) = {9};
Curve Loop(11) = {4, 17, -12, -20};
Surface(10) = {11};
Curve Loop(13) = {3, 20, -11, -19};
Surface(11) = {13};
Curve Loop(15) = {2, 19, -10, -18};
Surface(12) = {15};

// Outer cylyndrical surfaces
Curve Loop(17) = {6, 23, -14, -22};
Surface(13) = {17};
Curve Loop(19) = {5, 22, -13, -21};
Surface(14) = {19};
Curve Loop(21) = {8, 21, -16, -24};
Surface(15) = {21};
Curve Loop(23) = {7, 24, -15, -23};
Surface(16) = {23};

// Axial surfaces
Curve Loop(25) = {17, 29, -21, -25};
Surface(17) = {25};
Curve Loop(27) = {28, 24, -32, -20};
Surface(18) = {27};
Curve Loop(29) = {27, 23, -31, -19};
Surface(19) = {29};
Curve Loop(31) = {26, 22, -30, -18};
Surface(20) = {31};

// Volumes
Surface Loop(1) = {13, 2, 12, 6, 20, 19};
Volume(1) = {1};
Surface Loop(2) = {14, 5, 9, 1, 17, 20};
Volume(2) = {2};
Surface Loop(3) = {15, 8, 10, 4, 18, 17};
Volume(3) = {3};
Surface Loop(4) = {16, 3, 11, 7, 19, 18};
Volume(4) = {4};

// Physical surfaces
Physical Surface("bottom", 2) = {1, 2, 3, 4};
Physical Surface("diffuser", 3) = {13, 14, 15, 16};
Physical Surface("lid", 4) = {5, 6, 7, 8};
Physical Surface("sandscreen", 5) = {9, 10, 11, 12};

// Physical volume
Physical Volume("dom", 1) = {1, 2, 3, 4};

// Mesh sizes
Transfinite Curve {24, 22, 23, 21, 17, 18, 19, 20} = n_height+1 Using Progression 1;
Transfinite Curve {32, 31, 30, 29, 27, 28, 26, 25} = n_radial+1 Using Progression 1;
Transfinite Curve {11, 9, 10, 12, 14, 13, 16, 15, 7, 6, 5, 8, 3, 4, 1, 2} = n_tangential+1 Using Progression 1;

// Set all surfaces to transfinite
Transfinite Surface {1:20};

// Recombine surfaces to generate hexahedra
Recombine Surface {1:20};

// Set all volumes to transfinite
Transfinite Volume {1:4};