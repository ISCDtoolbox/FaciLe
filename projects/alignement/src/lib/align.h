#ifndef ALIGNDEF
#define ALIGNDEF

#include "mesh.h"

#include "super4pcs/algorithms/4pcs.h"
#include "super4pcs/algorithms/super4pcs.h"
#include "super4pcs/io/io.h"
#include "super4pcs/utils/geometry.h"

#include "icpPointToPlane.h"
#include "icpPointToPoint.h"

const float* super4PCS(Mesh &sourceMesh, Mesh &targetMesh, float overlap, float delta, int n_points);
const float* icp(Mesh &sourceMesh, Mesh &targetMesh, int maxIt, float inlierDist);

#endif
