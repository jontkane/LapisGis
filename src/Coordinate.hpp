#pragma once
#ifndef LP_COORDINATE_H
#define LP_COORDINATE_H

#include"LapisGisTypeDefs.hpp"

namespace lapis {
	struct CoordXY {
		coord_t x, y;
		CoordXY() : x(0), y(0) {}
		CoordXY(coord_t x, coord_t y) : x(x), y(y) {}
		bool operator==(const CoordXY& other) const = default;
	};
	struct CoordXYZ {
		coord_t x, y, z;
		CoordXYZ() : x(0), y(0), z(0) {}
		CoordXYZ(coord_t x, coord_t y, coord_t z) : x(x), y(y), z(z) {}
		bool operator==(const CoordXYZ& other) const = default;
	};
}

#endif