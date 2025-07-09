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

		virtual const OGRGeometry& gdalGeometryGeneric() const = 0;

		virtual Extent boundingBox() const = 0;

		const CoordRef& crs() const;
		void setCrs(const CoordRef& crs);

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
		const OGRPoint& gdalGeometry() const;
		const OGRGeometry& gdalGeometryGeneric() const override;

		coord_t x() const;
		coord_t y() const;

		Extent boundingBox() const override;
	private:
		OGRPoint _point;
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
		const OGRPolygon& gdalGeometry() const;
		const OGRGeometry& gdalGeometryGeneric() const override;

		void addInnerRing(const std::vector<CoordXY>& innerRing);

		std::vector<CoordXY> getOuterRing() const;
		int nInnerRings() const;
		std::vector<CoordXY> getInnerRing(int index) const;

		Extent boundingBox() const override;
		bool containsPoint(coord_t x, coord_t y) const;
		bool containsPoint(CoordXY xy) const;
		bool containsPoint(Point p) const;

		coord_t area() const;
	private:
		static std::vector<CoordXY> _coordsFromRing(const OGRLinearRing* ring);
		OGRPolygon _polygon;
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

		OGRwkbGeometryType gdalGeometryType() const override;
		const OGRMultiPolygon& gdalGeometry() const;
		const OGRGeometry& gdalGeometryGeneric() const override;

		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;

		void addPolygon(const Polygon& polygon);

		Extent boundingBox() const override;
		bool containsPoint(coord_t x, coord_t y) const;
		bool containsPoint(CoordXY xy) const;
		bool containsPoint(Point p) const;

		coord_t area() const;
	private:
		OGRMultiPolygon _multiPolygon;

		class iterator {
		public:
			iterator(OGRMultiPolygon* multiPolygon, size_t index);
			iterator& operator++();
			bool operator==(const iterator& other) const = default;
			Polygon operator*();
		private:
			OGRMultiPolygon* _multiPolygon;
			size_t _index;
		};
		class const_iterator {
		public:
			const_iterator(const OGRMultiPolygon* multiPolygon, size_t index);
			const_iterator& operator++();
			bool operator==(const const_iterator& other) const = default;
			const Polygon operator*();
		private:
			const OGRMultiPolygon* _multiPolygon;
			size_t _index;
		};
	};
}

#endif