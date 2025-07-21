#include"vector.hpp"

namespace lapis {

    WrongFieldTypeException::WrongFieldTypeException(const std::string& error) : std::runtime_error(error) {}

    AttributeTable::AttributeTable(const std::string& filename)
    {
        gdalAllRegisterThreadSafe();
        UniqueGdalDataset shp = vectorGDALWrapper(filename);
        if (!shp) {
            return;
        }
        OGRLayer* layer = shp->GetLayer(0);

        bool initFields = false;
        resize(layer->GetFeatureCount());

        for (const OGRFeatureUniquePtr& feature : layer) {
            if (!initFields) {
                for (int i = 0; i < feature->GetFieldCount(); ++i) {
                    OGRFieldDefn* field = feature->GetFieldDefnRef(i);
                    switch (field->GetType()) {
                    case OFTInteger:
                    case OFTInteger64:
                        addIntegerField(field->GetNameRef());
                        break;
                    case OFTReal:
                        addRealField(field->GetNameRef());
                        break;
                    case OFTString:
                        addStringField(field->GetNameRef(), field->GetWidth());
                        break;
                    default:
                        throw std::runtime_error("unimplemented field type when reading shapefile");
                    }
                }
                initFields = true;
            }
            OGRGeometry* gdalGeometry = feature->GetGeometryRef();
            for (int i = 0; i < feature->GetFieldCount(); ++i) {
                OGRFieldDefn* field = feature->GetFieldDefnRef(i);
                switch (field->GetType()) {
                case OFTInteger:
                case OFTInteger64:
                    setIntegerField(i, field->GetNameRef(), feature->GetFieldAsInteger64(field->GetNameRef()));
                    break;
                case OFTReal:
                    setRealField(i, field->GetNameRef(), feature->GetFieldAsDouble(field->GetNameRef()));
                    break;
                case OFTString:
                    setStringField(i, field->GetNameRef(), feature->GetFieldAsString(field->GetNameRef()));
                    break;
                default:
                    throw WrongFieldTypeException("Unimplemented field type when reading shapefile");
                }
            }
        }
    }
    AttributeTable::AttributeTable(const std::filesystem::path& filename)
        : AttributeTable(filename.string())
    {
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
    void AttributeTable::removeField(const std::string& name)
    {
        auto it = _fields.find(name);
        if (it != _fields.end()) {
            _fields.erase(it);
            _fieldNamesInOrder.erase(std::remove(_fieldNamesInOrder.begin(), _fieldNamesInOrder.end(), name), _fieldNamesInOrder.end());
        } else {
            throw std::runtime_error("Field not found: " + name);
        }
    }

    bool AttributeTable::fieldExists(const std::string& name) const
    {
        return _fields.contains(name);
    }

    void AttributeTable::resize(size_t nrow) {
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

    const std::vector<std::string>& AttributeTable::getAllFieldNames() const
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
    void AttributeTable::setStringFieldWidth(const std::string& name, size_t width)
    {
        if (_fields.at(name).type != FieldType::String) {
            throw WrongFieldTypeException("Wrong field type; expected string");
        }
        _fields.at(name).width = width;
        for (auto& value : _fields.at(name).values) {
            std::get<FixedWidthString>(value).set(std::get<FixedWidthString>(value).get(), width);
        }
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

    ConstAttributeRow AttributeTable::getRow(size_t index) const
    {
        return ConstAttributeRow(this, index);
    }
    MutableAttributeRow AttributeTable::getRow(size_t index)
    {
        return MutableAttributeRow(this, index);
    }
    AttributeTable::iterator AttributeTable::begin()
    {
        return iterator(this, 0);
    }
    AttributeTable::iterator AttributeTable::end()
    {
        return iterator(this, _nrow);
    }
    AttributeTable::const_iterator AttributeTable::begin() const
    {
        return const_iterator(this, 0);
    }
    AttributeTable::const_iterator AttributeTable::end() const
    {
        return const_iterator(this, _nrow);
    }
    MutableAttributeRow AttributeTable::front()
    {
        return MutableAttributeRow(this, 0);
    }
    MutableAttributeRow AttributeTable::back()
    {
        return MutableAttributeRow(this, _nrow - 1);
    }
    ConstAttributeRow AttributeTable::front() const
    {
        return ConstAttributeRow(this, 0);
    }
    ConstAttributeRow AttributeTable::back() const
    {
        return ConstAttributeRow(this, _nrow - 1);
    }

    void AttributeTable::_reserve(size_t n)
    {
        for (auto& [key, value] : _fields) {
            value.values.reserve(n);
        }
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
    const std::string& AttributeTable::FixedWidthString::get() const
    {
        return _data;
    }

}