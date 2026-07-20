// Double branched cover invariant for planar knotoids.
//
// C++ port of double_branched_cover/two_fold_py3.py (D. Goundaroulis et al.).
// Given the extended oriented Gauss code of a planar knotoid, this builds the
// double branched cover (a knot in a solid torus) and evaluates its Jones
// polynomial (in variables A and v).
//
// Input Gauss code (all three parts required):
//   gauss  : crossing sequence, positive = overpass, negative = underpass
//            (e.g. {1,-2,-1,2}). Arc labels are the positions 0..gauss.size().
//   signs  : sign per crossing (+1/-1), indexed by |label|-1 (e.g. {-1,-1}).
//   outer  : arc labels touching the outer region (e.g. {0,2,3}).

#ifndef DOUBLEBRANCHEDCOVER_HH
#define DOUBLEBRANCHEDCOVER_HH

#include <ctime>
#include <vector>
#include <Polynomial.hh>
#include <PlanarDiagram.hh>

// The double branched cover of a planar knotoid: a knot in a solid torus,
// described by its Gauss code plus the arcs bounding the two distinguished
// (outer / inner) regions of the annulus.
struct DBCLift
{
  std::vector<int> gauss;    // Gauss code of the lift
  std::vector<int> signs;    // crossing signs of the lift
  std::vector<int> outside;  // arc labels on the outer region boundary
  std::vector<int> inside;   // arc labels on the inner region boundary
};

// Build the double branched cover ("gc_in_st" in the reference).
DBCLift double_branched_cover_lift(const std::vector<int> &gauss,
                                   const std::vector<int> &signs,
                                   const std::vector<int> &outer);

// Jones polynomial (in A and v) of a knot in a solid torus ("jones_poly_st").
// timeout (in seconds, 0 = no timeout) throws std::runtime_error if exceeded.
Polynomial solid_torus_jones_polynomial(const DBCLift &lift, time_t timeout = 0);

// Full invariant from a planar knotoid's extended oriented Gauss code.
Polynomial double_branched_cover_polynomial(const std::vector<int> &gauss,
                                            const std::vector<int> &signs,
                                            const std::vector<int> &outer,
                                            time_t timeout = 0);

// Full invariant from a KnotoID planar PlanarDiagram (flag_planar == true).
Polynomial double_branched_cover_polynomial(PlanarDiagram &diagram, time_t timeout = 0);

#endif
