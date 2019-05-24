cl1 = 0.001;
cl2 = 0.005;
cl3 = 0.005;
Dy = 0.0;
radius = 0.1;
y = 0.1;
epsilon = 1e-10;
Point(1) = {0, y + epsilon, 0, cl1};
Point(2) = {radius, radius + y + epsilon, 0, cl2};
Point(3) = {-radius, radius + y + epsilon, 0, cl2};

Point(11) = {0, y, 0, cl2};
Point(12) = {radius, -radius + y, 0, cl2};
Point(13) = {-radius,- radius + y, 0, cl2};

Point(8) = {0, radius + y, 0, cl2};
Point(18) = {0, -radius + y, 0, cl2};

Circle(1) = {3, 8, 1};
Circle(2) = {1, 8, 2};

Circle(11) = {13, 18, 11};
Circle(12) = {11, 18, 12};

Line(3) = {2, 8};
Line(13) = {8, 3};

Line(30) = {12, 18};
Line(31) = {18, 13};


Line Loop(9) = {2, 3, 13, 1};
Plane Surface(9) = {9};

Line Loop(19) = {-12, -30, -31, -11};
Plane Surface(19) = {19};

Physical Line("contact_bottom") = {11, 12};
Physical Line("contact_top") = {1, 2};
Physical Line("bottom") = {30, 31};
Physical Line("top") = {3, 13};
Physical Surface("bot_body") = {19};
Physical Surface("top_body") = {9};