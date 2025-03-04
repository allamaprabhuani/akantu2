el_size = 0.1;
inclusion_radius = 0.25;
box_size = 1.0; // this is the edge length of the box

Point(1) = {-box_size/2, -box_size/2, 0, el_size};
Point(2) = {-box_size/2 + inclusion_radius, -box_size/2, 0, el_size};
Point(3) = {-box_size/2, -box_size/2 + inclusion_radius, 0, el_size};
Point(4) = {-box_size/2 - inclusion_radius, -box_size/2, 0, el_size};
Point(5) = {-box_size/2, -box_size/2 - inclusion_radius, 0, el_size};

Point(6) = {0, 0, 0, el_size};
Point(7) = {-box_size, 0, 0, el_size};
Point(8) = {-box_size, -box_size, 0, el_size};
Point(9) = {0, -box_size, 0, el_size};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line(5) = {9, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line Loop(9) = {4, 1, 2, 3};
Plane Surface(10) = {9};
Line Loop(11) = {5, 6, 7, 8};
Plane Surface(12) = {9, -11};
Physical Surface("aggregates") = {10};
Physical Surface("paste") = {12};
Physical Line("bottom") = {8};
Physical Line("right") = {5};
Physical Line("top") = {6};
Physical Line("left") = {7};
