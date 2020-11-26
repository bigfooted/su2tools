// Gmsh project created on Tue Nov 24 08:58:51 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {-1.50, 0, 0, 1.0};
//+
Point(3) = {-1.50, -0.015, 0, 1.0};
//+
Point(4) = {0, -0.015, 0, 1.0};
//+
Point(5) = {0.315, 0, 0, 1.0};
//+
Point(6) = {0.315,-0.0705, 0, 1.0};
//+
Point(7) = {1.065, 0, 0, 1.0};
//+
Point(8) = {1.065, -0.0705, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 4};
//+
Line(4) = {3, 4};
//+
Line(5) = {1, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {4, 6};
//+
Line(8) = {5, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {6, 8};
//+
Curve Loop(1) = {8, 9, -10, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, -7, -3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 3, -4, -1};
//+
Plane Surface(3) = {3};
//+
Physical Curve("inlet") = {1};
//+
Physical Curve("outlet") = {9};
//+
Physical Curve("wall_top") = {2, 5, 8};
//+
Physical Curve("wall_bottom") = {4, 7, 10};
//+
Transfinite Surface {3} = {2, 1, 4, 3};
//+
Transfinite Surface {2} = {1, 5, 6, 4};
//+
Transfinite Surface {1} = {5, 7, 8, 6};
//+
Transfinite Curve {1, 3} = 64 Using Bump 0.2;
//+
Transfinite Curve {3, 6} = 64 Using Bump 0.2;
//+
Transfinite Curve {6, 9} = 64 Using Bump 0.2;
//+
Transfinite Curve {2, 4} = 1000 Using Progression 0.998;
//+
Transfinite Curve {5, 7} = 400 Using Bump 0.88;
//+
Transfinite Curve {8, 10} = 300 Using Progression 1.006;
//+
Recombine Surface {3, 2, 1};
