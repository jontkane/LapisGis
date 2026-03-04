#include"gis_pch.hpp"
#include"Unit.hpp"

namespace lapis {
    Unit::Unit()
        : _type(UnitType::unknown), _name("unknown"), _convToSI(1.), _isUnknown(true)
    {
    }

    Unit::Unit(UnitType type, const std::string& name, double convToSI)
        : _type(type), _name(name), _convToSI(convToSI), _isUnknown(false)
    {
    }

    const std::string& Unit::name() const {
        return _name;
    }

    UnitType Unit::type() const {
        return _type;
    }

    bool Unit::isUnknown() const {
        return _isUnknown;
    }

    double Unit::convertOneToSI(double value) const {
        return value * _convToSI;
    }

    double Unit::convertOneFromSI(double value) const {
        return value / _convToSI;
    }

    double Unit::convertOneToThis(double value, const Unit& sourceUnits) const {
        if (_type != sourceUnits._type) {
            throw UnitMismatchException("Cannot convert from " + sourceUnits.name() + " to " + name() + " because they have different types");
        }
        return convertOneFromSI(sourceUnits.convertOneToSI(value));
    }

    double Unit::convertOneFromThis(double value, const Unit& destUnits) const {
        if (_type != destUnits._type) {
            throw UnitMismatchException("Cannot convert from " + name() + " to " + destUnits.name() + " because they have different types");
        }
        return destUnits.convertOneFromSI(convertOneToSI(value));
    }

    double Unit::convertOneToThis(double value, const std::optional<Unit>& sourceUnits) const {
        if (sourceUnits.has_value()) {
            if (_type != sourceUnits.value()._type) {
                throw UnitMismatchException("Cannot convert from " + sourceUnits.value().name() + " to " + name() + " because they have different types");
            }
            return convertOneToThis(value, sourceUnits.value());
        }
        else {
            return value;
        }
    }

    double Unit::convertOneFromThis(double value, const std::optional<Unit>& destUnits) const {
        if (destUnits.has_value()) {
            if (_type != destUnits.value()._type) {
                throw UnitMismatchException("Cannot convert from " + name() + " to " + destUnits.value().name() + " because they have different types");
            }
            return convertOneFromThis(value, destUnits.value());
        }
        else {
            return value;
        }
    }

    bool Unit::isConsistent(const Unit& other) const {
        bool consistentType = _type == other._type || _type == UnitType::unknown || other._type == UnitType::unknown;
        return consistentType && (isUnknown() || other.isUnknown() || std::abs(_convToSI - other._convToSI) < LAPIS_EPSILON);
    }


    UnitConverter::UnitConverter()
        : _conv(1.)
    {
    }

    UnitConverter::UnitConverter(const Unit& source, const Unit& dest)
        : _conv(source.convertOneToSI(1.) / dest.convertOneToSI(1.))
    {
        if (source._type != dest._type) {
            throw UnitMismatchException("Cannot convert from " + source.name() + " to " + dest.name() + " because they have different types");
        }
    }

    UnitConverter::UnitConverter(const std::optional<Unit>& source, const std::optional<Unit>& dest)
        : _conv(1.)
    {
        if (source.has_value() && dest.has_value()) {
            if (source.value()._type != dest.value()._type) {
                throw UnitMismatchException("Cannot convert from " + source.value().name() + " to " + dest.value().name() + " because they have different types");
            }
            _conv = source.value().convertOneToSI(1.) / dest.value().convertOneToSI(1.);
        }
    }

    double UnitConverter::convertOne(double value) const {
        return value * _conv;
    }

    double UnitConverter::operator()(double value) const {
        return convertOne(value);
    }

}