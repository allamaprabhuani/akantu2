LA = 128;
L = 202;
cl1 = 8;
Point(1) = {0, 0, -L, cl1};
Point(2) = {LA, 0, -L, cl1};
Point(3) = {0, LA, -L, cl1};
Point(4) = {LA, LA, -L, cl1};
Point(5) = {0, 0, L, cl1};
Point(6) = {LA, 0, L, cl1};
Point(7) = {0, LA, L, cl1};
Point(8) = {LA, LA, L, cl1};
Line(1) = {7, 8};
Line(2) = {8, 4};
Line(3) = {4, 3};
Line(4) = {3, 7};
Line(5) = {1, 5};
Line(6) = {5, 6};
Line(7) = {6, 2};
Line(8) = {2, 1};
Line(9) = {3, 1};
Line(10) = {7, 5};
Line(11) = {8, 6};
Line(12) = {4, 2};
Line Loop(14) = {1, 11, -6, -10};
Plane Surface(14) = {14};
Transfinite Surface {14};
Line Loop(16) = {3, 4, 1, 2};
Plane Surface(16) = {16};
Transfinite Surface {16};
Line Loop(18) = {6, 7, 8, 5};
Plane Surface(18) = {18};
Transfinite Surface {18};
Line Loop(20) = {3, 9, -8, -12};
Plane Surface(20) = {20};
Transfinite Surface {20};
Line Loop(22) = {4, 10, -5, -9};
Plane Surface(22) = {22};
Transfinite Surface {22};
Line Loop(24) = {11, 7, -12, -2};
Plane Surface(24) = {24};
Transfinite Surface {24};
Surface Loop(26) = {24, 14, 16, 20, 22, 18};
Volume(26) = {26};
Transfinite Volume {26};

Physical Volume(30) = {26};