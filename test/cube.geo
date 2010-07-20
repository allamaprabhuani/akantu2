// CUBE ON CUBE (2 BODIES) 3D:

// REMARKS:
//   Physical Surfaces are defined so that load_mesh_msh can load
// directly the surface element structure. When creating Plane Surface,
// the surface normal has to point to the inside of the body. Otherwise the
// solvecontact3d algorithm won't work. The Physical Surface number 
// corresponds to the face number. The face numbering has to start with 1!
//   Physical Volume defines the material of each body. This can be read
// by the load_mesh_msh function in adlib.

// Mesh size
h = 0.5;    // Top cube

// Dimensions of top cube
Lx = 1;
Ly = 1;
Lz = 1;

// ------------------------------------------
// Geometry
// ------------------------------------------

// Base Cube
Point(101) = { 0.0, 0.0, 0.0, h}; // Bottom Face
Point(102) = { Lx,  0.0, 0.0, h}; // Bottom Face
Point(103) = { Lx,  Ly, 0.0,  h}; // Bottom Face
Point(104) = { 0.0, Ly, 0.0,  h}; // Bottom Face

Point(111) = { 0.0, 0.0, Lz,  h}; // Top Face
Point(112) = { Lx,  0.0, Lz,  h}; // Top Face
Point(113) = { Lx,  Ly,  Lz,  h}; // Top Face
Point(114) = { 0.0, Ly,  Lz,  h}; // Top Face

// Base Cube
Line(101) = {101,102}; // Bottom Face
Line(102) = {102,103}; // Bottom Face
Line(103) = {103,104}; // Bottom Face
Line(104) = {104,101}; // Bottom Face

Line(111) = {111,112}; // Top Face
Line(112) = {112,113}; // Top Face
Line(113) = {113,114}; // Top Face
Line(114) = {114,111}; // Top Face

Line(121) = {101,111}; // Side Faces
Line(122) = {102,112}; // Side Faces
Line(123) = {103,113}; // Side Faces
Line(124) = {104,114}; // Side Faces

// Base Cube
Line Loop(101) = {101:104}; // Base
Line Loop(111) = {111:114}; // Top

Line Loop(121) = {101, 122, -111, -121}; // Side
Line Loop(122) = {102, 123, -112, -122}; // Side
Line Loop(123) = {103, 124, -113, -123}; // Side
Line Loop(124) = {104, 121, -114, -124}; // Side

// Base Cube
Plane Surface(101) = {101};  // Base
Plane Surface(111) = {-111}; // Top

Plane Surface(121)  = {-121}; // Side
Plane Surface(122)  = {-122}; // Side
Plane Surface(123)  = {-123}; // Side
Plane Surface(124)  = {-124}; // Side

// Base Cube
Physical Surface(4) = {111}; // Top Surface of Base Cube --> Contact Surface
Physical Surface(5) = {121, 122, 123, 124};
Physical Surface(6) = {101}; // Base Surface of Base Cube --> Boundary Condition


// Base Cube
Surface Loop(101) = {-101, -111, -121, -122, -123, -124};

// Base Cube
Volume(101) = {101};

// Base Cube
Physical Volume(2) = {101};
