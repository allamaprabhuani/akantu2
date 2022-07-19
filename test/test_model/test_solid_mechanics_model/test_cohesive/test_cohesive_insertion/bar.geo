// Lenght of the bar (m)
L = 50e-3;
x0 = -L/2;
xf = L/2;

// Number of triangular elements (n_el)
n_el = 2;

// Length of each linear element (h)
h = L/(n_el/2);


Point(1) = { x0, 0, 0, h };
Point(2) = { xf, 0, 0, h };
Point(3) = { xf, h, 0, h };
Point(4) = { x0, h, 0, h };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};
Physical Line("YBlocked") = {1,3};
Physical Line("right") = {2};
Physical Line("left") = {4};

Transfinite Surface {1} Left;
