#pragma once
#ifndef LAPIS_GIS_HPP
#define LAPIS_GIS_HPP

#include"Alignment.hpp"
#include"Coordinate.hpp"
#include"CoordRef.hpp"
#include"CoordTransform.hpp"
#include"CoordVector.hpp"
#include"CropView.hpp"
#include"CurrentLasPoint.hpp"
#include"Extent.hpp"
#include"GDALWrappers.hpp"
#include"Geometry.hpp"
#include"GisExceptions.hpp"
#include"LapisGisTypeDefs.hpp"
#include"LasExtent.hpp"
#include"LasFilter.hpp"
#include"LasIO.hpp"
#include"LasReader.hpp"
#include"MultiBandRaster.hpp"
#include"ProjWrappers.hpp"
#include"QuadExtent.hpp"
#include"Raster.hpp"
#include"RasterAlgos.hpp"
#include"Unit.hpp"
#include"Vector.hpp"

namespace lapis {

    //Clears certain problematic environment variables which might be set unknowingly by the user (e.g. in a conda environment with an additional install of gdal or proj)
    //Also tells proj where to find proj.db. If you pass the empty string, it will try to find proj.db in the same folder as the executable. If you pass a path, it will look for proj.db in that path.
    bool lapisGisInit(const std::string& projDbPath = "");
}

#endif