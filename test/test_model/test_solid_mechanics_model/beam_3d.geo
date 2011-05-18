cl1 = 0.3;
Point(1) = {0, 0, 0, cl1};
Point(2) = {10, 0, 0, cl1};
Point(3) = {10, 1, 0, cl1};
Point(4) = {5, 1, 0, cl1};
Point(5) = {0, 1, 0, cl1};
Point(6) = {0, 0, 1, cl1};
Point(7) = {10, 0, 1, cl1};
Point(11) = {10, 1, 1, cl1};
Point(15) = {5, 1, 1, cl1};
Point(19) = {0, 1, 1, cl1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Line(10) = {6, 7};
Line(11) = {7, 11};
Line(12) = {11, 15};
Line(13) = {15, 19};
Line(14) = {19, 6};
Line(16) = {1, 6};
Line(17) = {2, 7};
Line(21) = {3, 11};
Line(25) = {4, 15};
Line(29) = {5, 19};
Line Loop(7) = {1, 2, 3, 4, 5};
Plane Surface(7) = {7};
Line Loop(18) = {1, 17, -10, -16};
Ruled Surface(18) = {18};
Line Loop(22) = {2, 21, -11, -17};
Ruled Surface(22) = {22};
Line Loop(26) = {3, 25, -12, -21};
Ruled Surface(26) = {26};
Line Loop(30) = {4, 29, -13, -25};
Ruled Surface(30) = {30};
Line Loop(34) = {5, 16, -14, -29};
Ruled Surface(34) = {34};
Line Loop(35) = {10, 11, 12, 13, 14};
Plane Surface(35) = {35};
Surface Loop(37) = {7, 18, 22, 26, 30, 34, 35};
Volume(37) = {37};
Physical Volume(1) = {37};
