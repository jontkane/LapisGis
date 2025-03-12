#pragma once
#ifndef lp_lapisgistypedefs_h
#define lp_lapisgistypedefs_h

namespace lapis {

	using coord_t = double;
	using cell_t = int64_t;
	using rowcol_t = int32_t;
	constexpr coord_t LAPIS_EPSILON = 0.0001;
	using intensity_t = uint32_t;
	using metric_t = float;
}

#endif