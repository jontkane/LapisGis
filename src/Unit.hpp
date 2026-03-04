#pragma once
#ifndef lp_unit_h
#define lp_unit_h

#include"GisExceptions.hpp"
#include"LapisGisTypeDefs.hpp"

namespace lapis {

	enum class UnitType {
		linear,
		angular,
		area,
		volume,
		mass,
		unknown
	};

	class UnitConverter;
	template<UnitType TYPE>
	class SpecificUnitConverter;

	class Unit {
		friend class UnitConverter;
	public:
		Unit();
        Unit(UnitType type, const std::string& name, double convToSI);

        const std::string& name() const;
        UnitType type() const;
        bool isUnknown() const;

		double convertOneToSI(double value) const;
        double convertOneFromSI(double value) const;
		double convertOneToThis(double value, const Unit& sourceUnits) const;
        double convertOneFromThis(double value, const Unit& destUnits) const;
		double convertOneToThis(double value, const std::optional<Unit>& sourceUnits) const;
        double convertOneFromThis(double value, const std::optional<Unit>& destUnits) const;

		bool isConsistent(const Unit& other) const; //returns true if they have the same type and conversion factor
		bool operator==(const Unit& other) const = default;
	protected:
		UnitType _type;
		std::string _name;
		double _convToSI; //multiply by this to convert to SI
		bool _isUnknown;
	};

	template<UnitType TYPE>
	class SpecificUnit : public Unit {
        friend class SpecificUnitConverter<TYPE>;
	public:
		SpecificUnit();
        SpecificUnit(const std::string& name, double convToSI);

        double convertOneToThis(double value, const SpecificUnit<TYPE>& sourceUnits) const;
		double convertOneFromThis(double value, const SpecificUnit<TYPE>& destUnits) const;
		double convertOneToThis(double value, const std::optional<SpecificUnit<TYPE>>& sourceUnits) const;
        double convertOneFromThis(double value, const std::optional<SpecificUnit<TYPE>>& destUnits) const;

        bool isConsistent(const SpecificUnit<TYPE>& other) const; //returns true if the conversion between these units is the identity, or within a rounding error of identity
        bool operator==(const SpecificUnit<TYPE>& other) const = default;

		template<UnitType OTHER_TYPE>
			requires (OTHER_TYPE != TYPE)
		double convertOneToThis(double value, const SpecificUnit<OTHER_TYPE>& sourceUnits) const = delete
			template<UnitType OTHER_TYPE>
			requires (OTHER_TYPE != TYPE)
		double convertOneFromThis(double value, const SpecificUnit<OTHER_TYPE>&destUnits) const = delete;
		template<UnitType OTHER_TYPE>
			requires (OTHER_TYPE != TYPE)
		double convertOneToThis(double value, const std::optional<SpecificUnit<OTHER_TYPE>>& sourceUnits) const = delete;
		template<UnitType OTHER_TYPE>
			requires (OTHER_TYPE != TYPE)
		double convertOneFromThis(double value, const std::optional<SpecificUnit<OTHER_TYPE>>& destUnits) const = delete;
		using Unit::convertOneToThis;
		using Unit::convertOneFromThis;
	};

    using LinearUnit = SpecificUnit<UnitType::linear>;
    using AngularUnit = SpecificUnit<UnitType::angular>;
    using AreaUnit = SpecificUnit<UnitType::area>;
    using VolumeUnit = SpecificUnit<UnitType::volume>;
    using MassUnit = SpecificUnit<UnitType::mass>;

	class UnitConverter {
	public:
		UnitConverter();
		UnitConverter(const Unit& source, const Unit& dest);
		UnitConverter(const std::optional<Unit>& source, const std::optional<Unit>& dest);

        double convertOne(double value) const;
        double operator()(double value) const;

		//This function converts values in a contiguous container, which is not necessarily an array of coord_t
		//The span variable indicates the distance between successive values to convert
		//For example, if we have: struct XYZ {T x, y, z;};
		//and a std::vector<XYZ> of size 1000 named v, and we want to convert the z values only, the call would be:
		//convertManyInPlace(&v[0].z, 1000, sizeof(XYZ));
		//if span is unspecified, then the assumption is that po is a pointer to an array of T
		template<class T>
		void convertManyInPlace(T* po, size_t count, size_t span = sizeof(T)) const;
	protected:
		double _conv; //the value to multiply the source from to get the dest
	};

	template<UnitType TYPE>
	class SpecificUnitConverter : public UnitConverter {
	public:
		SpecificUnitConverter();
		SpecificUnitConverter(const SpecificUnit<TYPE>& source, const SpecificUnit<TYPE>& dest);
		SpecificUnitConverter(const std::optional<SpecificUnit<TYPE>>& source, const std::optional<SpecificUnit<TYPE>>& dest);
	private:
	};

    using LinearUnitConverter = SpecificUnitConverter<UnitType::linear>;
    using AngularUnitConverter = SpecificUnitConverter<UnitType::angular>;
    using AreaUnitConverter = SpecificUnitConverter<UnitType::area>;
    using VolumeUnitConverter = SpecificUnitConverter<UnitType::volume>;
    using MassUnitConverter = SpecificUnitConverter<UnitType::mass>;

	namespace linearUnitPresets {
		const LinearUnit meter{ "metre", 1. }; //Spelled 'metre' to match the convention in WKT strings
		const LinearUnit internationalFoot{ "foot",0.3048 };
		const LinearUnit usSurveyFoot{ "US survey foot", 0.30480060960122 };
		const LinearUnit internationalInch{ "inch", 0.0254 };
		const LinearUnit usSurveyInch{ "US survey inch", 0.0254000508001016 };
		const LinearUnit centimeter{ "centimetre", 0.01 };
		const LinearUnit unknownLinear{};
	}
	namespace angularUnitPresets {
		const AngularUnit radian{ "radian", 1. };
		const AngularUnit degree{ "degree", 0.0174532925199433 };
		const AngularUnit unknownAngular{};
	}

	namespace areaUnitPresets {
		const AreaUnit squareMeter{ "square metre", 1. };
		const AreaUnit squareFoot{ "square foot", 0.09290304 };
		const AreaUnit squareUSSurveyFoot{ "square US survey foot", 0.092903411613275 };
		const AreaUnit squareInch{ "square inch", 0.00064516 };
		const AreaUnit squareUSSurveyInch{ "square US survey inch", 0.000645161290322581 };
		const AreaUnit squareCentimeter{ "square centimetre", 0.0001 };
		const AreaUnit hectare{ "hectare", 10000. };
		const AreaUnit acre{ "acre", 4046.8564224 };
		const AreaUnit squareKilometer{ "square kilometre", 1000000. };
		const AreaUnit unknownArea{};
	}

	namespace volumeUnitPresets {
		const VolumeUnit cubicMeter{ "cubic metre", 1. };
		const VolumeUnit cubicFoot{ "cubic foot", 0.028316846592 };
		const VolumeUnit cubicUSSurveyFoot{ "cubic US survey foot", 0.0283170189843 };
		const VolumeUnit cubicInch{ "cubic inch", 0.000016387064 };
		const VolumeUnit cubicUSSurveyInch{ "cubic US survey inch", 0.0000163871431186369 };
		const VolumeUnit cubicCentimeter{ "cubic centimetre", 0.000001 };
		const VolumeUnit unknownVolume{};
	}


	template<UnitType TYPE>
	inline SpecificUnit<TYPE>::SpecificUnit()
	{
		_type = TYPE;
		_name = "unknown";
		_convToSI = 1.;
        _isUnknown = true;
	}

    template<UnitType TYPE>
    inline SpecificUnit<TYPE>::SpecificUnit(const std::string& name, double convToSI)
		: Unit(TYPE, name, convToSI)
	{
    }

    template<UnitType TYPE>
	inline double SpecificUnit<TYPE>::convertOneToThis(double value, const SpecificUnit<TYPE>& sourceUnits) const {
        return value * sourceUnits._convToSI / _convToSI;
    }

    template<UnitType TYPE>
	inline double SpecificUnit<TYPE>::convertOneFromThis(double value, const SpecificUnit<TYPE>& destUnits) const {
		return value * _convToSI / destUnits._convToSI;
    }

    template<UnitType TYPE>
	inline double SpecificUnit<TYPE>::convertOneToThis(double value, const std::optional<SpecificUnit<TYPE>>& sourceUnits) const {
		if (sourceUnits.has_value()) {
			return convertOneToThis(value, sourceUnits.value());
		}
		else {
			return value;
		}
    }

    template<UnitType TYPE>
	inline double SpecificUnit<TYPE>::convertOneFromThis(double value, const std::optional<SpecificUnit<TYPE>>& destUnits) const {
		if (destUnits.has_value()) {
			return convertOneFromThis(value, destUnits.value());
		}
		else {
			return value;
		}
	}

    template<UnitType TYPE>
	inline bool SpecificUnit<TYPE>::isConsistent(const SpecificUnit<TYPE>& other) const {
		return isUnknown() || other.isUnknown() || std::abs(_convToSI - other._convToSI) < LAPIS_EPSILON;
	}

	template<class T>
	inline void UnitConverter::convertManyInPlace(T* po, size_t count, size_t span) const
	{
		for (size_t i = 0; i < count; i++) {
			*po = convertOne(*po);
			po = (T*)((char*)po + span);
        }
	}

	template<UnitType TYPE>
	inline SpecificUnitConverter<TYPE>::SpecificUnitConverter()
		: UnitConverter()
	{
    }

    template<UnitType TYPE>
	inline SpecificUnitConverter<TYPE>::SpecificUnitConverter(const SpecificUnit<TYPE>& source, const SpecificUnit<TYPE>& dest)
		: UnitConverter(source, dest)
	{
	}

    template<UnitType TYPE>
	inline SpecificUnitConverter<TYPE>::SpecificUnitConverter(const std::optional<SpecificUnit<TYPE>>& source, const std::optional<SpecificUnit<TYPE>>& dest)
		: UnitConverter(source, dest)
	{
	}
}

#endif