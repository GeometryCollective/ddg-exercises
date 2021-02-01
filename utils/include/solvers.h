#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <complex>

using namespace geometrycentral;
using namespace geometrycentral::surface;

SparseMatrix<double> sparseInverseDiagonal(SparseMatrix<double>& M);

double residual(const SparseMatrix<std::complex<double>>& A, const Vector<std::complex<double>>& x);

Vector<std::complex<double>> solveInversePowerMethod(const SparseMatrix<std::complex<double>>& A);