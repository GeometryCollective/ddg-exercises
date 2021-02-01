#pragma once

#include <algorithm>
#include <cmath>

#include "geometrycentral/surface/vertex_position_geometry.h"
using namespace geometrycentral;

std::array<double, 3> mapToColor(double x, double min, double max, const std::string& colormap_name);

std::array<double, 3> hsv(double h, double s, double v);