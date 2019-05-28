a = 1.0;
r = 2*a;
angle = 15 * Pi/180.;
height = a/4.0 * Sin(angle);
gap = 0;
nbpoints = 10.0;
cl1 = a/(4.0*nbpoints);
cl2 = r/(20.0*nbpoints);

Point(1) = {-a/4, height + gap, 0, 2*cl1};
Point(2) = {a/4, height + gap, 0, 2*cl1};
Point(3) = {a/4, a/4., 0, cl1};
Point(4) = {-a/4, a/4., 0, cl1};
Point(5) = {r/2, 0, 0, cl2};
Point(6) = {-r/2, 0, 0, cl2};
Point(7) = {-r/2, -r/2, 0, 5*cl2};
Point(8) = { r/2, -r/2, 0, 5*cl2};
Point(9) = {0, 0 + gap, 0, 2*cl1};

Line(1) = {1, 9};
Line(2) = {9, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 1};
Line(6) = {7, 8};
Line(7) = {8, 5};
Line(8) = {5, 6};
Line(9) = {6, 7};

Line Loop(1) = {1, 2, 3, 4, 5};
Line Loop(2) = {6, 7, 8, 9};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line("contact_top") = {1, 2};
Physical Line("contact_bottom") = {8};
Physical Line("top") = {4};
Physical Line("bottom") = {6};

Physical Surface("bot_body") = {2};
Physical Surface("top_body") = {1};