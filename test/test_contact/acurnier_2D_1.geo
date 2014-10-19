cl1 = 2;
cl2 = -0.15;
cl3 = -0.5;
Point(1) = {-1, cl3, 0, cl1};
Point(2) = {-1, 0.5, 0, cl1};
Point(3) = {1, 0.5, 0, cl1};
Point(4) = {1, cl3, 0, cl1};
Point(5) = {0, 0.7+cl2, 0, cl1};
Point(6) = {0.7, 1.4+cl2, 0, cl1};
Point(7) = {0, 1.4+cl2, 0, cl1};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 5};
Line Loop(9) = {1, 2, 3, 4};
Plane Surface(9) = {9};
Line Loop(11) = {5, 6, 7};
Plane Surface(11) = {11};
Physical Line("Bottom") = {1};
Physical Line("Top") = {6};
Physical Line("Contact") = {3};
Physical Line("Dummy") = {2, 4, 5, 7};
Physical Surface("Top_body") = {11};
Physical Surface("Bottom_body") = {9};
