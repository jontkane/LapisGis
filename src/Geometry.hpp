#pragma once
#ifndef lapis_geometry_h
#define lapis_geometry_h

#include"gis_pch.hpp"
#include"CoordRef.hpp"
#include"Extent.hpp"
#include"QuadExtent.hpp"

namespace lapis {

	class WrongGeometryTypeException : public std::runtime_error {
	public:
		WrongGeometryTypeException(const std::string& error);
	};

	class Geometry {
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbUnknown;
		using GdalEquivalent = OGRGeometry;

		virtual OGRwkbGeometryType gdalGeometryType() const;

		virtual std::unique_ptr<OGRGeometry> gdalGeometryGeneric() const = 0;

		virtual Extent boundingBox() const = 0;

		const CoordRef& crs() const;
		void setCrs(const CoordRef& crs);

		void projectInPlace(const CoordRef& newCrs);
		virtual void projectInPlace(const CoordTransform& transform) = 0;

		virtual ~Geometry() = default;

	protected:
		CoordRef _crs;
		Geometry() = default;
		Geometry(const Geometry&) = default;
	};
	class Point : public Geometry {
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbPoint;
		using GdalEquivalent = OGRPoint;

		Point() = default;
		Point(const OGRGeometry& geom);
		Point(const OGRGeometry& geom, const CoordRef& crs);
		Point(coord_t x, coord_t y);
		Point(coord_t x, coord_t y, const CoordRef& crs);
		Point(CoordXY xy);
		Point(CoordXY xy, const CoordRef& crs);

		OGRwkbGeometryType gdalGeometryType() const override;
		std::unique_ptr<OGRPoint> gdalGeometry() const;
		std::unique_ptr<OGRGeometry> gdalGeometryGeneric() const override;

		coord_t x() const;
		coord_t y() const;

		Extent boundingBox() const override;

        void projectInPlace(const CoordTransform& transform) override;
	private:
		CoordXY _point;
		void _sharedConstructorFromGdal(const OGRGeometry& geom);
	};
	class Polygon : public Geometry {
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbPolygon;
		using GdalEquivalent = OGRPolygon;

		Polygon() = default;
		Polygon(const OGRGeometry& geom);
		Polygon(const OGRGeometry& geom, const CoordRef& crs);
		Polygon(const std::vector<CoordXY>& outerRing);
		Polygon(const std::vector<CoordXY>& outerRing, const CoordRef& crs);
		Polygon(const Extent& e);
		Polygon(const QuadExtent& q);

		OGRwkbGeometryType gdalGeometryType() const override;
		std::unique_ptr<OGRPolygon> gdalGeometry() const;
		std::unique_ptr<OGRGeometry> gdalGeometryGeneric() const override;

		void addInnerRing(const std::vector<CoordXY>& innerRing);

		const std::vector<CoordXY>& getOuterRing() const;
		int nInnerRings() const;
		const std::vector<CoordXY>& getInnerRing(int index) const;
		const std::vector<CoordXY>& getInnerRingUnsafe(int index) const;
		

		Extent boundingBox() const override;
		bool containsPoint(coord_t x, coord_t y) const;
		bool containsPoint(CoordXY xy) const;
		bool containsPoint(Point p) const;

		coord_t area() const;

        void projectInPlace(const CoordTransform& transform) override;
	private:
		std::vector<CoordXY> _outerRing;
        std::vector<std::vector<CoordXY>> _innerRings;

		static coord_t _areaFromRing(const std::vector<CoordXY>& ring);
		void _sharedConstructorFromGdal(const OGRGeometry& geom);
	};
	class MultiPolygon : public Geometry {
	private:
		class iterator;
		class const_iterator;
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbMultiPolygon;
		using GdalEquivalent = OGRMultiPolygon;

		MultiPolygon() = default;
		MultiPolygon(const OGRGeometry& geom);
		MultiPolygon(const OGRGeometry& geom, const CoordRef& crs);

		size_t nPolygon() const;

		OGRwkbGeometryType gdalGeometryType() const override;
		std::unique_ptr<OGRMultiPolygon> gdalGeometry() const;
		std::unique_ptr<OGRGeometry> gdalGeometryGeneric() const override;

		std::vector<Polygon>::iterator begin();
		std::vector<Polygon>::iterator end();
		std::vector<Polygon>::const_iterator begin() const;
		std::vector<Polygon>::const_iterator end() const;

		void addPolygon(const Polygon& polygon);

		Extent boundingBox() const override;
		bool containsPoint(coord_t x, coord_t y) const;
		bool containsPoint(CoordXY xy) const;
		bool containsPoint(Point p) const;

		coord_t area() const;

        void projectInPlace(const CoordTransform& transform) override;
	private:
        std::vector<Polygon> _polygons;
		void _sharedConstructorFromGdal(const OGRGeometry& geom, const CoordRef& crs);
	};
}

#endif