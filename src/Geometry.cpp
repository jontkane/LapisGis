#include"Geometry.hpp"
#include"GisExceptions.hpp"

namespace lapis {

    WrongGeometryTypeException::WrongGeometryTypeException(const std::string& error) : std::runtime_error(error) {}

    OGRwkbGeometryType Geometry::gdalGeometryType() const
    {
        return OGRwkbGeometryType::wkbUnknown;
    }
    const CoordRef& Geometry::crs() const
    {
        return _crs;
    }
    void Geometry::setCrs(const CoordRef& crs)
    {
        _crs = crs;
    }

    void Geometry::projectInPlace(const CoordRef& newCrs)
    {
        const CoordTransform& transform = CoordTransformFactory::getTransform(_crs, newCrs);
        projectInPlace(transform);
    }
    Point::Point(const OGRGeometry& geom)
    {
        _sharedConstructorFromGdal(geom);
        setCrs(CoordRef(geom.getSpatialReference()));
    }
    Point::Point(const OGRGeometry& geom, const CoordRef& crs)
    {
        _sharedConstructorFromGdal(geom);
        setCrs(crs);
    }
    Point::Point(coord_t x, coord_t y)
        : _point(x, y)
    {
    }
    Point::Point(coord_t x, coord_t y, const CoordRef& crs)
        : _point(x,y)
    {
        _crs = crs;
    }
    Point::Point(CoordXY xy)
        : _point(xy)
    {
    }
    Point::Point(CoordXY xy, const CoordRef& crs)
        : _point(xy)
    {
        _crs = crs;
    }
    OGRwkbGeometryType Point::gdalGeometryType() const
    {
        return OGRwkbGeometryType::wkbPoint;
    }
    std::unique_ptr<OGRPoint> Point::gdalGeometry() const
    {
        std::unique_ptr<OGRPoint> gdalPoint = std::make_unique<OGRPoint>(_point.x, _point.y);
        gdalPoint->assignSpatialReference(_crs.gdalSpatialRef());
        return gdalPoint;
    }
    std::unique_ptr<OGRGeometry> Point::gdalGeometryGeneric() const
    {
        std::unique_ptr<OGRPoint> gdalPoint = gdalGeometry();
        return std::unique_ptr<OGRGeometry>(dynamic_cast<OGRGeometry*>(gdalPoint.release()));
    }
    coord_t Point::x() const
    {
        return _point.x;
    }
    coord_t Point::y() const
    {
        return _point.y;
    }

    Extent Point::boundingBox() const
    {
        return Extent(x(), x(), y(), y(), _crs);
    }

    void Point::projectInPlace(const CoordTransform& transform)
    {
        _point = transform.transformSingleXY(_point.x, _point.y);
        _crs = transform.dst();
    }

    bool Point::overlapsExtent(const Extent& e) const
    {
        if (_crs.isConsistentHoriz(e.crs())) {
            return overlapsExtentSameCrs(e);
        }
        else {
            const CoordTransform& transform = CoordTransformFactory::getTransform(_crs, e.crs());
            CoordXY transformedPoint = transform.transformSingleXY(_point.x, _point.y);
            return e.contains(transformedPoint.x, transformedPoint.y);
        }
    }

    bool Point::overlapsExtentSameCrs(const Extent& e) const
    {
        return e.contains(x(), y());
    }

    void Point::_sharedConstructorFromGdal(const OGRGeometry& geom)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbPoint) {
            throw WrongGeometryTypeException("Wrong geometry; expected Point");
        }
        const OGRPoint* gdalPoint = geom.toPoint();
        _point = CoordXY{ gdalPoint->getX(), gdalPoint->getY() };
    }


    Polygon::Polygon(const OGRGeometry& geom)
    {
        _sharedConstructorFromGdal(geom);
        setCrs(CoordRef(geom.getSpatialReference()));
    }
    Polygon::Polygon(const OGRGeometry& geom, const CoordRef& crs)
    {
        _sharedConstructorFromGdal(geom);
        setCrs(crs);
    }
    Polygon::Polygon(const std::vector<CoordXY>& outerRing)
    {
        _outerRing = outerRing;
        if (_outerRing.back() != _outerRing.front()) {
            _outerRing.push_back(_outerRing.front());
        }
        if (_outerRing.size() < 4) {
            throw std::runtime_error("Rings must have at least 3 points");
        }
    }
    Polygon::Polygon(const std::vector<CoordXY>& outerRing, const CoordRef& crs)
        : Polygon(outerRing)
    {
        setCrs(crs);
    }
    Polygon::Polygon(const Extent& e)
    {
        setCrs(e.crs());

        _outerRing.reserve(5);
        _outerRing.emplace_back(e.xmin(), e.ymax());
        _outerRing.emplace_back(e.xmin(), e.ymin());
        _outerRing.emplace_back(e.xmax(), e.ymin());
        _outerRing.emplace_back(e.xmax(), e.ymax());
        _outerRing.push_back(_outerRing.front());
    }
    Polygon::Polygon(const QuadExtent& q)
    {
        setCrs(q.crs());
        const CoordXYVector& coords = q.coords();
        _outerRing.reserve(coords.size() + 1);
        for (const CoordXY& xy : coords) {
            _outerRing.push_back(xy);
        }
        _outerRing.push_back(_outerRing.front());
    }
    OGRwkbGeometryType Polygon::gdalGeometryType() const
    {
        return OGRwkbGeometryType::wkbPolygon;
    }
    std::unique_ptr<OGRPolygon> Polygon::gdalGeometry() const
    {
        std::unique_ptr<OGRPolygon> gdalPolygon = std::make_unique<OGRPolygon>();
        OGRLinearRing* gdalOuterRing = new OGRLinearRing();
        for (const CoordXY& xy : _outerRing) {
            gdalOuterRing->addPoint(xy.x, xy.y);
        }
        gdalOuterRing->closeRings(); //shouldn't be necessary but just in case
        gdalPolygon->addRingDirectly(gdalOuterRing);

        for (const std::vector<CoordXY>& innerRing : _innerRings) {
            OGRLinearRing* gdalInnerRing = new OGRLinearRing();
            for (const CoordXY& xy : innerRing) {
                gdalInnerRing->addPoint(xy.x, xy.y);
            }
            gdalInnerRing->closeRings();
            gdalPolygon->addRingDirectly(gdalInnerRing);
        }
        gdalPolygon->assignSpatialReference(_crs.gdalSpatialRef());
        return gdalPolygon;
    }
    std::unique_ptr<OGRGeometry> Polygon::gdalGeometryGeneric() const
    {
        return std::unique_ptr<OGRGeometry>(dynamic_cast<OGRGeometry*>(gdalGeometry().release()));
    }
    void Polygon::addInnerRing(const std::vector<CoordXY>& innerRing)
    {
        std::vector<CoordXY> copy = innerRing;
        if (copy.back() != copy.front()) {
            copy.push_back(copy.front());
        }
        if (copy.size() < 4) {
            throw std::runtime_error("Rings must have at least 3 points");
        }
        _innerRings.push_back(std::move(copy));
    }
    const std::vector<CoordXY>& Polygon::getOuterRing() const {
        return _outerRing;
    }
    int Polygon::nInnerRings() const
    {
        return (int)_innerRings.size();
    }
    const std::vector<CoordXY>& Polygon::getInnerRing(int index) const
    {
        return _innerRings.at(index);
    }
    const std::vector<CoordXY>& Polygon::getInnerRingUnsafe(int index) const
    {
        return _innerRings[index];
    }
    Extent Polygon::boundingBox() const
    {
        coord_t xmin = std::numeric_limits<coord_t>::max();
        coord_t xmax = std::numeric_limits<coord_t>::lowest();
        coord_t ymin = std::numeric_limits<coord_t>::max();
        coord_t ymax = std::numeric_limits<coord_t>::lowest();
        for (const CoordXY& xy : _outerRing) {
            if (xy.x < xmin) xmin = xy.x;
            if (xy.x > xmax) xmax = xy.x;
            if (xy.y < ymin) ymin = xy.y;
            if (xy.y > ymax) ymax = xy.y;
        }
        return Extent{ xmin, xmax, ymin, ymax, _crs };
    }
    bool Polygon::containsPoint(coord_t x, coord_t y) const
    {
        auto pointInRing = [](const std::vector<CoordXY>& ring, coord_t x, coord_t y) {
            int nvert = (int)(ring.size() - 1);
            bool within = false;
            for (int i = 0, j = nvert - 1; i < nvert; j = i++) {
                coord_t ix = ring[i].x;
                coord_t iy = ring[i].y;
                coord_t jx = ring[j].x;
                coord_t jy = ring[j].y;
                if (((iy > y) != (jy > y)) &&
                    (x < (jx - ix) * (y - iy) / (jy - iy) + ix)) {
                    within = !within;
                }
            }
            return within;
            };

        if (!pointInRing(_outerRing, x, y)) {
            return false;
        }
        for (const std::vector<CoordXY>& innerRing : _innerRings) {
            if (pointInRing(innerRing, x, y)) {
                return false;
            }
        }
        return true;
    }
    bool Polygon::containsPoint(CoordXY xy) const
    {
        return containsPoint(xy.x, xy.y);
    }
    bool Polygon::containsPoint(Point p) const
    {
        return containsPoint(p.x(), p.y());
    }
    coord_t Polygon::area() const
    {
        coord_t totalArea = _areaFromRing(_outerRing);
        for (const std::vector<CoordXY>& innerRing : _innerRings) {
            totalArea -= _areaFromRing(innerRing);
        }
        return totalArea;
    }
    void Polygon::_sharedConstructorFromGdal(const OGRGeometry& geom) {
        if (wkbFlatten(geom.getGeometryType()) != wkbPolygon) {
            throw WrongGeometryTypeException("Wrong geometry; expected Polygon");
        }
        const OGRPolygon* gdalPolygon = geom.toPolygon();

        const OGRLinearRing* exteriorRing = gdalPolygon->getExteriorRing();
        for (const OGRPoint& point : *exteriorRing) {
            _outerRing.emplace_back(point.getX(), point.getY());
        }
        for (int i = 0; i < gdalPolygon->getNumInteriorRings(); i++) {
            const OGRLinearRing* innerRing = gdalPolygon->getInteriorRing(i);
            std::vector<CoordXY> innerCoords;
            for (const OGRPoint& point : *innerRing) {
                innerCoords.emplace_back(point.getX(), point.getY());
            }
            _innerRings.emplace_back(std::move(innerCoords));
        }
    }

    void Polygon::projectInPlace(const CoordTransform& transform)
    {
        transform.transformXY(_outerRing);
        for (std::vector<CoordXY>& innerRing : _innerRings) {
            transform.transformXY(innerRing);
        }
        setCrs(transform.dst());
    }

    bool Polygon::overlapsExtent(const Extent& e) const
    {
        if (_crs.isConsistentHoriz(e.crs())) {
            return overlapsExtentSameCrs(e);
        }

        Polygon otherPoly{ QuadExtent(e,CoordTransformFactory::getTransform(e.crs(), _crs)) };

        //they overlap if either contains a vertex of the other, or if any pair of edges intersect.
        //first check if bounding box overlaps
        if (!boundingBox().overlapsUnsafe(otherPoly.boundingBox())) {
            return false;
        }

        //check if other vertices are inside this
        //note that otherpoly is guaranteed to not have inner rings
        for (const CoordXY& xy : otherPoly._outerRing) {
            if (containsPoint(xy)) {
                return true;
            }
        }

        //check if this polygon's vertices are inside the other
        for (const CoordXY& xy : _outerRing) {
            if (otherPoly.containsPoint(xy)) {
                return true;
            }
        }
        for (const std::vector<CoordXY>& innerRing : _innerRings) {
            for (const CoordXY& xy : innerRing) {
                if (otherPoly.containsPoint(xy)) {
                    return true;
                }
            }
        }

        auto lineSegmentsIntersect = [](coord_t x1, coord_t y1, coord_t x2, coord_t y2, coord_t x3, coord_t y3, coord_t x4, coord_t y4) {
            // Check if line segments (x1,y1)-(x2,y2) and (x3,y3)-(x4,y4) intersect
            auto orientation = [](coord_t ax, coord_t ay, coord_t bx, coord_t by, coord_t cx, coord_t cy) {
                coord_t val = (by - ay) * (cx - bx) - (bx - ax) * (cy - by);
                if (val == 0) return 0; // collinear
                return (val > 0) ? 1 : 2; // clock or counterclock wise
            };
            int o1 = orientation(x1, y1, x2, y2, x3, y3);
            int o2 = orientation(x1, y1, x2, y2, x4, y4);
            int o3 = orientation(x3, y3, x4, y4, x1, y1);
            int o4 = orientation(x3, y3, x4, y4, x2, y2);
            if (o1 != o2 && o3 != o4)
                return true;
            return false; // Doesn't fall in any of the above cases
            };

        //check if any edge of this polygon intersects any edge of the other
        for (size_t i = 0; i < _outerRing.size() - 1; i++) {
            for (size_t j = 0; j < otherPoly._outerRing.size() - 1; j++) {
                if (lineSegmentsIntersect(_outerRing[i].x, _outerRing[i].y, _outerRing[i + 1].x, _outerRing[i + 1].y,
                    otherPoly._outerRing[j].x, otherPoly._outerRing[j].y, otherPoly._outerRing[j + 1].x, otherPoly._outerRing[j + 1].y)) {
                    return true;
                }
            }
        }
        for (const std::vector<CoordXY>& innerRing : _innerRings) {
            for (size_t i = 0; i < innerRing.size() - 1; i++) {
                for (size_t j = 0; j < otherPoly._outerRing.size() - 1; j++) {
                    if (lineSegmentsIntersect(innerRing[i].x, innerRing[i].y, innerRing[i + 1].x, innerRing[i + 1].y,
                        otherPoly._outerRing[j].x, otherPoly._outerRing[j].y, otherPoly._outerRing[j + 1].x, otherPoly._outerRing[j + 1].y)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    bool Polygon::overlapsExtentSameCrs(const Extent& e) const
    {
        //they overlap if either contains a vertex of the other, or if any pair of edges intersect.

        //first check if bounding box overlaps
        if (!boundingBox().overlapsUnsafe(e)) {
            return false;
        }

        //check if polygon vertices are inside the extent
        for (const CoordXY& xy : _outerRing) {
            if (e.contains(xy.x, xy.y)) {
                return true;
            }
        }

        //check if extent corners are inside the polygon
        std::vector<CoordXY> extentCorners = {
            {e.xmin(), e.ymin()},
            {e.xmin(), e.ymax()},
            {e.xmax(), e.ymin()},
            {e.xmax(), e.ymax()}
        };
        for (const CoordXY& xy : extentCorners) {
            if (containsPoint(xy)) {
                return true;
            }
        }

        //check if a line segment intersects a horizontal line segment
        auto segmentIntersectsHorizontal = [](coord_t x1, coord_t y1, coord_t x2, coord_t y2, coord_t y, coord_t xMin, coord_t xMax) {
            if ((y1 > y) != (y2 > y)) {
                coord_t intersectX = (x2 - x1) * (y - y1) / (y2 - y1) + x1;
                return (intersectX >= xMin) && (intersectX <= xMax);
            }
            return false;
        };
        //same but for vertical line segments
        auto segmentIntersectsVertical = [](coord_t x1, coord_t y1, coord_t x2, coord_t y2, coord_t x, coord_t yMin, coord_t yMax) {
            if ((x1 > x) != (x2 > x)) {
                coord_t intersectY = (y2 - y1) * (x - x1) / (x2 - x1) + y1;
                return (intersectY >= yMin) && (intersectY <= yMax);
            }
            return false;
            };

        auto segmentIntersectsExtent = [&](coord_t x1, coord_t y1, coord_t x2, coord_t y2) {
            //check bounding box first
            coord_t segXMin = std::min(x1, x2);
            coord_t segXMax = std::max(x1, x2);
            coord_t segYMin = std::min(y1, y2);
            coord_t segYMax = std::max(y1, y2);

            Extent segmentExtent(segXMin, segXMax, segYMin, segYMax);
            if (!segmentExtent.overlapsUnsafe(e)) {
                return false;
            }

            //check if it intersects any of the four edges of the extent
            if (segmentIntersectsHorizontal(x1, y1, x2, y2, e.ymin(), e.xmin(), e.xmax())) {
                return true;
            };
            if (segmentIntersectsHorizontal(x1, y1, x2, y2, e.ymax(), e.xmin(), e.xmax())) {
                return true;
            };
            if (segmentIntersectsVertical(x1, y1, x2, y2, e.xmin(), e.ymin(), e.ymax())) {
                return true;
            };
            if (segmentIntersectsVertical(x1, y1, x2, y2, e.xmax(), e.ymin(), e.ymax())) {
                return true;
            };
            return false;
            };

        //check if any edge of the polygon intersects the extent
        for (size_t i = 0; i < _outerRing.size() - 1; i++) {
            if (segmentIntersectsExtent(_outerRing[i].x, _outerRing[i].y, _outerRing[i + 1].x, _outerRing[i + 1].y)) {
                return true;
            }
        }
        for (const std::vector<CoordXY>& innerRing : _innerRings) {
            for (size_t i = 0; i < innerRing.size() - 1; i++) {
                if (segmentIntersectsExtent(innerRing[i].x, innerRing[i].y, innerRing[i + 1].x, innerRing[i + 1].y)) {
                    return true;
                }
            }
        }
        return false;
    }

    coord_t Polygon::_areaFromRing(const std::vector<CoordXY>& ring)
    {
        int nPoints = (int)ring.size();
        if (nPoints < 4) { //3 for a triangle, plus the duplicated point
            return 0; 
        }
        CoordXY pointOne, pointTwo;
        coord_t area = 0;
        pointOne = ring[0];

        //this is the shoelace formula
        for (int i = 1; i < nPoints; ++i) {
            pointTwo = ring[i];

            coord_t x1 = pointOne.x;
            coord_t y1 = pointOne.y;
            coord_t x2 = pointTwo.x;
            coord_t y2 = pointTwo.y;
            area += (x1 * y2) - (x2 * y1);

            std::swap(pointOne, pointTwo);
        }

        return std::abs(area) / 2.0;
    }

    MultiPolygon::MultiPolygon(const OGRGeometry& geom)
    {
        CoordRef crs{ geom.getSpatialReference() };
        _sharedConstructorFromGdal(geom, crs);
        setCrs(crs);
    }
    MultiPolygon::MultiPolygon(const OGRGeometry& geom, const CoordRef& crs)
    {
        _sharedConstructorFromGdal(geom, crs);
        setCrs(crs);
    }
    OGRwkbGeometryType MultiPolygon::gdalGeometryType() const
    {
        return OGRwkbGeometryType::wkbMultiPolygon;
    }
    std::unique_ptr<OGRMultiPolygon> MultiPolygon::gdalGeometry() const
    {
        std::unique_ptr<OGRMultiPolygon> gdalMultiPolygon = std::make_unique<OGRMultiPolygon>();
        for (const Polygon& poly : _polygons) {
            gdalMultiPolygon->addGeometry(poly.gdalGeometry().get());
        }
        gdalMultiPolygon->assignSpatialReference(_crs.gdalSpatialRef());
        return gdalMultiPolygon;
    }
    std::unique_ptr<OGRGeometry> MultiPolygon::gdalGeometryGeneric() const
    {
        return std::unique_ptr<OGRGeometry>(dynamic_cast<OGRGeometry*>(gdalGeometry().release()));
    }
    void MultiPolygon::addPolygon(const Polygon& polygon)
    {
        if (!polygon.crs().isConsistent(_crs)) {
            throw CRSMismatchException("Polygon CRS does not match MultiPolygon CRS");
        }
        _polygons.push_back(polygon);
    }
    std::vector<Polygon>::iterator MultiPolygon::begin() {
        return _polygons.begin();
    }
    std::vector<Polygon>::iterator MultiPolygon::end() {
        return _polygons.end();
    }
    std::vector<Polygon>::const_iterator MultiPolygon::begin() const
    {
        return _polygons.begin();
    }
    std::vector<Polygon>::const_iterator MultiPolygon::end() const
    {
        return _polygons.end();
    }
    Extent MultiPolygon::boundingBox() const
    {
        if (!_polygons.size()) {
            return Extent{ 0, 0, 0, 0, _crs };
        }
        Extent out = _polygons[0].boundingBox();
        for (size_t i = 1; i < _polygons.size(); i++) {
            out = extendExtent(out, _polygons[i].boundingBox());
        }
        return out;
    }
    bool MultiPolygon::containsPoint(coord_t x, coord_t y) const
    {
        for (const Polygon& poly : _polygons) {
            if (poly.containsPoint(x, y)) {
                return true;
            }
        }
        return false;
    }
    bool MultiPolygon::containsPoint(CoordXY xy) const
    {
        return containsPoint(xy.x, xy.y);
    }
    bool MultiPolygon::containsPoint(Point p) const
    {
        return containsPoint(p.x(), p.y());
    }
    coord_t MultiPolygon::area() const
    {
        coord_t totalArea = 0;
        for (Polygon poly : *this) {
            totalArea += poly.area();
        }
        return totalArea;
    }
    void MultiPolygon::projectInPlace(const CoordTransform& transform)
    {
        for (Polygon& poly : _polygons) {
            poly.projectInPlace(transform);
        }
        setCrs(transform.dst());
    }
    bool MultiPolygon::overlapsExtent(const Extent& e) const
    {
        for (const Polygon& poly : _polygons) {
            if (poly.overlapsExtent(e)) {
                return true;
            }
        }
        return false;
    }
    bool MultiPolygon::overlapsExtentSameCrs(const Extent& e) const
    {
        for (const Polygon& poly : _polygons) {
            if (poly.overlapsExtentSameCrs(e)) {
                return true;
            }
        }
        return false;
    }
    void MultiPolygon::_sharedConstructorFromGdal(const OGRGeometry& geom, const CoordRef& crs)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbMultiPolygon
            && wkbFlatten(geom.getGeometryType()) != wkbPolygon) {
            throw WrongGeometryTypeException("Wrong geometry; expected MultiPolygon");
        }
        if (wkbFlatten(geom.getGeometryType()) == wkbPolygon) {
            _polygons.emplace_back(geom, crs);
            return;
        }
        const OGRMultiPolygon* gdalMultiPolygon = geom.toMultiPolygon();
        for (int i = 0; i < gdalMultiPolygon->getNumGeometries(); i++) {
            _polygons.emplace_back(*gdalMultiPolygon->getGeometryRef(i), crs);
        }
    }
    size_t MultiPolygon::nPolygon() const
    {
        return _polygons.size();
    }
}