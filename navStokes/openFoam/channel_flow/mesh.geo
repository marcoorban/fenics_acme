lc = 1;
// delta is defined to be the distance between the midstream and the plate.
delta = 1;
// Streamwise length defined as X = (4 * pi * delta)
X = 12.56;
// Spanwise length defined as Z = (2 * pi * delta)
Z = 6.28;
Point(1) = {0, -delta, 0, lc};
Point(2) = {X, -delta, 0, lc};
Point(3) = {X, delta, 0, lc};
Point(4) = {0, delta, 0, lc};
Line(11) = {1, 2}; //bottom
Line(12) = {2, 3}; //right
Line(13) = {3, 4}; //top
Line(14) = {4, 1}; //left
Curve Loop(15) = {11, 12, 13, 14};
Plane Surface(21) = {15};

// Set up the transfinite curves
// Transfinite Curve{c}=N sets N *points* (N-1 divisions).
Transfinite Curve{11} = 192;
Transfinite Curve{13} = 192;
// Wall-normal grading: 129 pts, Bump 0.01054 puts the first cell CENTER
// at y+ = 0.05 for Re_tau=180 (KMM 1987), i.e. first node off the wall
// at 5.556e-4 m (= 2 * 2.778e-4 m, since the cell center sits halfway
// between the wall and the first interior node). Growth ratio ~1.10
// near the wall, tuned empirically via an isolated single-line test
// (see log.md Learnings).
Transfinite Curve{12} = 129 Using Bump 0.011566;
Transfinite Curve{14} = 129 Using Bump 0.011566;
Transfinite Surface{21} = {1, 2, 3, 4};

Recombine Surface{21};

// Physical zones
// Layers{N} sets N divisions directly (verified empirically -- unlike
// Transfinite Curve, no off-by-one here, and no separate Transfinite
// Curve is needed on the extruded spanwise edges).
out[] = Extrude {0, 0, Z} { Surface{21}; Layers{159}; Recombine; };
// Extrude returns a list with the created surfaces
// in the following order:
// out[0] = outer cap (z=1), out[1] = volume,
// out[2..5] = side surfaces extruded from Line(11..14):
//   out[2] <- Line(11), y=-1 (bottom wall)
//   out[3] <- Line(12), x=10 (right / outlet)
//   out[4] <- Line(13), y=1  (top wall)
//   out[5] <- Line(14), x=0  (left / inlet)

Physical Surface("wall")  = {out[2], out[4]};
Physical Surface("left")  = {out[5]};
Physical Surface("right") = {out[3]};
Physical Surface("inner") = {21};
Physical Surface("outer") = {out[0]};
Physical Volume("domain") = {out[1]};

Mesh 3;
// gmshToFoam (openfoam.org) mis-parses gmsh's default MSH 4.1 output --
// force the older 2.2 format so gmshToFoam reads it correctly.
Mesh.MshFileVersion = 2.2;
Save "mesh.msh";

