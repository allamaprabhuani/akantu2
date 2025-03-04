cl1 = 0.0025;
cl2 = 0.0025;
cl3 = 0.005;
Dy = 0.0;
radius = 0.1;
y = 0.1;
epsilon = 0.0;
Point(1) = {0, y, 0, cl1};
Point(2) = {radius, radius + y, 0, cl2};
Point(3) = {-radius, radius + y, 0, cl2};
Point(4) = {0.25, 0.1 -epsilon, 0, cl3};
Point(5) = {-0.25, 0.1-epsilon, 0, cl3};
Point(6) = {-0.25, -0.2, 0, 10*cl2};
Point(7) = {0.25, -0.2, 0, 10*cl2};
Point(8) = {0, radius + y, 0, cl2};
Circle(1) = {3, 8, 1};
Circle(2) = {1, 8, 2};
Line(3) = {2, 8};
Line(13) = {8, 3};
Line(4) = {6, 7};
Line(5) = {7, 4};
Line(6) = {4, 5};
Line(7) = {5, 6};
Line Loop(9) = {2, 3, 13, 1};
Plane Surface(9) = {9};
Line Loop(11) = {6, 7, 4, 5};
Plane Surface(11) = {11};
Physical Line("contact_bottom") = {6};
Physical Line("contact_top") = {1, 2};
Physical Line("top") = {3, 13};
Physical Line("bottom") = {4};
Physical Surface("top_body") = {9};
Physical Surface("bot_body") = {11};