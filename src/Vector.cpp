#include"vector.hpp"

namespace lapis {

    static bool pointInGdalRing(const OGRLinearRing* ring, coord_t x, coord_t y) {

        //hoping that this algorithm is faster than gdal's, which is very very slow
        int nvert = ring->getNumPoints() - 1; //gdal stores the first point twice; but this algorithm cares about the true number of vertices
        OGRPoint iPoint;
        OGRPoint jPoint;
        bool within = false;
        for (int i = 0, j = nvert-1; i < nvert; j = i++) {
            ring->getPoint(i, &iPoint);
            ring->getPoint(j, &jPoint);
            coord_t ix = iPoint.getX();
            coord_t iy = iPoint.getY();
            coord_t jx = jPoint.getX();
            coord_t jy = jPoint.getY();
            if (((iy > y) != (jy > y)) &&
                (x < (jx - ix) * (y - iy) / (jy - iy) + ix)) {
                within = !within;
            }
        }
        return within;
    }
    static bool pointInGdalPolygon(const OGRPolygon* poly, coord_t x, coord_t y) {
        for (int i = 0; i < poly->getNumInteriorRings(); i++) {
            if (pointInGdalRing(poly->getInteriorRing(i), x, y)) {
                return false;
            }
        }
        return pointInGdalRing(poly->getExteriorRing(), x, y);
    }

    WrongGeometryTypeException::WrongGeometryTypeException(const std::string& error) : std::runtime_error(error) {}
    WrongFieldTypeException::WrongFieldTypeException(const std::string& error) : std::runtime_error(error) {}

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
    Point::Point(const OGRGeometry& geom)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbPoint) {
            throw WrongGeometryTypeException("Wrong geometry; expected Point");
        }
        _point = *geom.toPoint();
        setCrs(CoordRef(_point.getSpatialReference()));
    }
    Point::Point(const OGRGeometry& geom, const CoordRef& crs)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbPoint) {
            throw WrongGeometryTypeException("Wrong geometry; expected Point");
        }
        _point = *geom.toPoint();
        setCrs(crs);
    }
    Point::Point(coord_t x, coord_t y)
    {
        _point = OGRPoint(x, y);
    }
    Point::Point(coord_t x, coord_t y, const CoordRef& crs)
    {
        _point = OGRPoint(x, y);
        _crs = crs;
    }
    Point::Point(CoordXY xy)
    {
        _point = OGRPoint(xy.x, xy.y);
    }
    Point::Point(CoordXY xy, const CoordRef& crs)
    {
        _point = OGRPoint(xy.x, xy.y);
        _crs = crs;
    }
    OGRwkbGeometryType Point::gdalGeometryType() const
    {
        return OGRwkbGeometryType::wkbPoint;
    }
    const OGRPoint& Point::gdalGeometry() const
    {
        return _point;
    }
    const OGRGeometry& Point::gdalGeometryGeneric() const
    {
        return _point;
    }
    coord_t Point::x() const
    {
        return _point.getX();
    }
    coord_t Point::y() const
    {
        return _point.getY();
    }

    Extent Point::boundingBox() const
    {
        return Extent(x(), x(), y(), y(), _crs);
    }

    static void addPointToRing(OGRLinearRing& ring, coord_t x, coord_t y) {
        OGRPoint point;
        point.setX(x);
        point.setY(y);
        ring.addPoint(&point);
    }
    Polygon::Polygon(const OGRGeometry& geom)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbPolygon) {
            throw WrongGeometryTypeException("Wrong geometry; expected Polygon");
        }
        _polygon = *geom.toPolygon();
        setCrs(CoordRef(_polygon.getSpatialReference()));
    }
    Polygon::Polygon(const OGRGeometry& geom, const CoordRef& crs)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbPolygon) {
            throw WrongGeometryTypeException("Wrong geometry; expected Polygon");
        }
        _polygon = *geom.toPolygon();
        setCrs(crs);
    }
    Polygon::Polygon(const std::vector<CoordXY>& outerRing)
    {
        OGRLinearRing gdalRing;
        for (const CoordXY& xy : outerRing) {
            addPointToRing(gdalRing, xy.x, xy.y);
        }
        gdalRing.closeRings();
        _polygon.addRing(&gdalRing);
    }
    Polygon::Polygon(const std::vector<CoordXY>& outerRing, const CoordRef& crs)
    {
        OGRLinearRing gdalRing;
        for (const CoordXY& xy : outerRing) {
            addPointToRing(gdalRing, xy.x, xy.y);
        }
        gdalRing.closeRings();
        _polygon.addRing(&gdalRing);
        _crs = crs;
    }
    Polygon::Polygon(const Extent& e)
    {
        _crs = e.crs();
        OGRLinearRing gdalRing;
        addPointToRing(gdalRing, e.xmin(), e.ymax());
        addPointToRing(gdalRing, e.xmin(), e.ymin());
        addPointToRing(gdalRing, e.xmax(), e.ymin());
        addPointToRing(gdalRing, e.xmax(), e.ymax());
        gdalRing.closeRings();
        _polygon.addRing(&gdalRing);
    }
    Polygon::Polygon(const QuadExtent& q)
    {
        _crs = q.crs();
        const CoordXYVector& coords = q.coords();
        OGRLinearRing gdalRing;
        for (size_t i = 0; i < coords.size(); ++i) {
            addPointToRing(gdalRing, coords[i].x, coords[i].y);
        }
        gdalRing.closeRings();
        _polygon.addRing(&gdalRing);
    }
    OGRwkbGeometryType Polygon::gdalGeometryType() const
    {
        return OGRwkbGeometryType::wkbPolygon;
    }
    const OGRPolygon& Polygon::gdalGeometry() const
    {
        return _polygon;
    }
    const OGRGeometry& Polygon::gdalGeometryGeneric() const
    {
        return _polygon;
    }
    void Polygon::addInnerRing(const std::vector<CoordXY>& innerRing)
    {
        OGRLinearRing gdalRing;
        for (const CoordXY& xy : innerRing) {
            addPointToRing(gdalRing, xy.x, xy.y);
        }
        gdalRing.closeRings();
        _polygon.addRing(&gdalRing);
    }
    std::vector<CoordXY> Polygon::getOuterRing() const {
        return _coordsFromRing(_polygon.getExteriorRing());
    }
    int Polygon::nInnerRings() const
    {
        return _polygon.getNumInteriorRings();
    }
    std::vector<CoordXY> Polygon::getInnerRing(int index) const
    {
        return _coordsFromRing(_polygon.getInteriorRing(index));
    }
    Extent Polygon::boundingBox() const
    {
        OGREnvelope g;
        _polygon.getEnvelope(&g);
        return Extent{ g.MinX,g.MaxX,g.MinY,g.MaxY,_crs };
    }
    bool Polygon::containsPoint(coord_t x, coord_t y) const
    {
        /*OGRPoint p;
        p.setX(x);
        p.setY(y);
        return p.Within(&_polygon);*/
        return pointInGdalPolygon(&_polygon, x, y);
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
        return _polygon.get_Area();
    }
    std::vector<CoordXY> Polygon::_coordsFromRing(const OGRLinearRing* ring)
    {
        std::vector<CoordXY> out;
        int nCoords = ring->getNumPoints();
        out.reserve(nCoords);
        for (OGRPoint point : ring) {
            out.emplace_back(point.getX(), point.getY());
        }
        return out;
    }

    MultiPolygon::MultiPolygon(const OGRGeometry& geom)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbMultiPolygon
            && wkbFlatten(geom.getGeometryType()) != wkbPolygon) {
            throw WrongGeometryTypeException("Wrong geometry; expected MultiPolygon");
        }
        _multiPolygon = *geom.toMultiPolygon();
        setCrs(CoordRef(_multiPolygon.getSpatialReference()));
    }
    MultiPolygon::MultiPolygon(const OGRGeometry& geom, const CoordRef& crs)
    {
        if (wkbFlatten(geom.getGeometryType()) != wkbMultiPolygon
            && wkbFlatten(geom.getGeometryType()) != wkbPolygon) {
            throw WrongGeometryTypeException("Wrong geometry; expected MultiPolygon");
        }
        _multiPolygon = *geom.toMultiPolygon();
        setCrs(crs);
    }
    OGRwkbGeometryType MultiPolygon::gdalGeometryType() const
    {
        return OGRwkbGeometryType::wkbMultiPolygon;
    }
    const OGRMultiPolygon& MultiPolygon::gdalGeometry() const
    {
        return _multiPolygon;
    }
    const OGRGeometry& MultiPolygon::gdalGeometryGeneric() const
    {
        return _multiPolygon;
    }
    void MultiPolygon::addPolygon(const Polygon& polygon)
    {
        _multiPolygon.addGeometry(&polygon.gdalGeometry());
    }
    MultiPolygon::iterator MultiPolygon::begin() {
        return iterator(&_multiPolygon, 0);
    }
    MultiPolygon::iterator MultiPolygon::end() {
        return iterator(&_multiPolygon, _multiPolygon.getNumGeometries());
    }
    MultiPolygon::const_iterator MultiPolygon::begin() const
    {
        return const_iterator(&_multiPolygon, 0);
    }
    MultiPolygon::const_iterator MultiPolygon::end() const
    {
        return const_iterator(&_multiPolygon, _multiPolygon.getNumGeometries());
    }
    Extent MultiPolygon::boundingBox() const
    {
        OGREnvelope g;
        _multiPolygon.getEnvelope(&g);
        return Extent{ g.MinX,g.MaxX,g.MinY,g.MaxY,_crs };
    }
    bool MultiPolygon::containsPoint(coord_t x, coord_t y) const
    {
        /*OGRPoint p;
        p.setX(x);
        p.setY(y);
        return p.Within(&_multiPolygon);*/
        int nPoly = _multiPolygon.getNumGeometries();
        for (int i = 0; i < _multiPolygon.getNumGeometries(); i++) {
            const OGRPolygon* poly = _multiPolygon.getGeometryRef(i);
            if (pointInGdalPolygon(poly, x, y)) {
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

    MultiPolygon::iterator::iterator(OGRMultiPolygon* multiPolygon, size_t index) : _multiPolygon(multiPolygon), _index(index) {}
    MultiPolygon::iterator& MultiPolygon::iterator::operator++() {
        _index++;
        return *this;
    }
    Polygon MultiPolygon::iterator::operator*() {
        OGRPolygon* poly = _multiPolygon->getGeometryRef((int)_index)->toPolygon();
        Polygon out{ *poly };
        return Polygon(*poly);
    }

    MultiPolygon::const_iterator::const_iterator(const OGRMultiPolygon* multiPolygon, size_t index)
        : _multiPolygon(multiPolygon), _index(index)
    {
    }
    MultiPolygon::const_iterator& MultiPolygon::const_iterator::operator++() {
        _index++;
        return *this;
    }
    const Polygon MultiPolygon::const_iterator::operator*() {
        const OGRPolygon* poly = _multiPolygon->getGeometryRef((int)_index)->toPolygon();
        Polygon out{ *poly };
        return Polygon(*poly);
    }

    void AttributeTable::setStringField(size_t index, const std::string& name, const std::string& value)
    {
        if (_fields.at(name).type != FieldType::String) {
            throw WrongFieldTypeException("Wrong field type; expected string");
        }
        FixedWidthString& fws = std::get<FixedWidthString>(_fields.at(name).values.at(index));
        fws.set(value, _fields.at(name).width);
    }
    void AttributeTable::setIntegerField(size_t index, const std::string& name, int64_t value)
    {
        if (_fields.at(name).type != FieldType::Integer) {
            throw WrongFieldTypeException("Wrong field type; expected integer");
        }
        _fields.at(name).values.at(index) = value;
    }
    void AttributeTable::setRealField(size_t index, const std::string& name, double value)
    {
        if (_fields.at(name).type != FieldType::Real) {
            throw WrongFieldTypeException("Wrong field type; expected real");
        }
        _fields.at(name).values.at(index) = value;
    }
    void AttributeTable::_reserve(size_t n)
    {
        for (auto& [key, value] : _fields) {
            value.values.reserve(n);
        }
    }
    void AttributeTable::addStringField(const std::string& name, size_t width)
    {
        std::vector<Variant> values = std::vector<Variant>(_nrow, Variant(FixedWidthString()));
        Field newField{ FieldType::String, width, std::move(values) };
        _fields.emplace(name, std::move(newField));
        _fieldNamesInOrder.push_back(name);
    }
    void AttributeTable::addIntegerField(const std::string& name)
    {
        std::vector<Variant> values = std::vector<Variant>(_nrow, Variant((int64_t)0));
        Field newField{ FieldType::Integer, 0, std::move(values) };
        _fields.emplace(name, std::move(newField));
        _fieldNamesInOrder.push_back(name);
    }
    void AttributeTable::addRealField(const std::string& name)
    {
        std::vector<Variant> values = std::vector<Variant>(_nrow, Variant((double)0));
        Field newField{ FieldType::Real, 0, std::move(values) };
        _fields.emplace(name, std::move(newField));
        _fieldNamesInOrder.push_back(name);
    }
    void AttributeTable::resize(size_t nrow)
    {
        _nrow = nrow;
        for (auto& keyValue : _fields) {
            switch (keyValue.second.type) {
            case FieldType::Integer:
                keyValue.second.values.resize(_nrow, Variant((int64_t)0));
                break;
            case FieldType::Real:
                keyValue.second.values.resize(_nrow, Variant((double)0.));
                break;
            case FieldType::String:
                keyValue.second.values.resize(_nrow, Variant(FixedWidthString()));
                break;
            }
        }
    }
    void AttributeTable::addRow()
    {
        _nrow++;
        resize(_nrow);
    }
    size_t AttributeTable::nFeature() const
    {
        return _nrow;
    }
    std::vector<std::string> AttributeTable::getAllFieldNames() const
    {
        return _fieldNamesInOrder;
    }
    FieldType AttributeTable::getFieldType(const std::string& name) const
    {
        return _fields.at(name).type;
    }
    size_t AttributeTable::getStringFieldWidth(const std::string& name) const
    {
        if (_fields.at(name).type != FieldType::String) {
            throw WrongFieldTypeException("Wrong field type; expected string");
        }
        return _fields.at(name).width;
    }
    const std::string& AttributeTable::getStringField(size_t index, const std::string& name) const
    {
        if (_fields.at(name).type != FieldType::String) {
            throw WrongFieldTypeException("Wrong field type; expected string");
        }
        return std::get<FixedWidthString>(_fields.at(name).values.at(index)).get();
    }
    int64_t AttributeTable::getIntegerField(size_t index, const std::string& name) const
    {
        if (_fields.at(name).type != FieldType::Integer) {
            throw WrongFieldTypeException("Wrong field type; expected integer");
        }
        return std::get<int64_t>(_fields.at(name).values.at(index));
    }
    double AttributeTable::getRealField(size_t index, const std::string& name) const
    {
        if (_fields.at(name).type != FieldType::Real) {
            throw WrongFieldTypeException("Wrong field type; expected real");
        }
        return std::get<double>(_fields.at(name).values.at(index));
    }

    const std::string& AttributeTable::FixedWidthString::get() const
    {
        return _data;
    }
    void AttributeTable::FixedWidthString::set(const std::string& value, size_t width)
    {
        if (value.size() > width) {
            _data = value.substr(0, width);
        }
        else {
            _data = value + std::string(width - value.size(), '\0');
        }
    }

}