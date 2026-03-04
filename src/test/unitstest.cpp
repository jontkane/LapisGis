#include "test_pch.hpp"
#include "..\Unit.hpp"

namespace lapis {
    TEST(UnitTest, DefaultConstructor) {
        Unit u;
        EXPECT_EQ(u.type(), UnitType::unknown);
        EXPECT_EQ(u.name(), "unknown");
        EXPECT_TRUE(u.isUnknown());
        EXPECT_NEAR(u.convertOneToSI(1.0), 1.0, 0.0001);
    }

    TEST(UnitTest, ParameterizedConstructor) {
        Unit u(UnitType::linear, "meter", 1.0);
        EXPECT_EQ(u.type(), UnitType::linear);
        EXPECT_EQ(u.name(), "meter");
        EXPECT_FALSE(u.isUnknown());

        Unit u2(UnitType::linear, "foot", 0.3048);
        EXPECT_EQ(u2.name(), "foot");
        EXPECT_NEAR(u2.convertOneToSI(1.0), 0.3048, 0.0001);
    }

    TEST(UnitTest, ConversionToSI) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        EXPECT_NEAR(meter.convertOneToSI(5.0), 5.0, 0.0001);
        EXPECT_NEAR(foot.convertOneToSI(1.0), 0.3048, 0.0001);
        EXPECT_NEAR(foot.convertOneToSI(10.0), 3.048, 0.0001);
    }

    TEST(UnitTest, ConversionFromSI) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        EXPECT_NEAR(meter.convertOneFromSI(5.0), 5.0, 0.0001);
        EXPECT_NEAR(foot.convertOneFromSI(0.3048), 1.0, 0.0001);
        EXPECT_NEAR(foot.convertOneFromSI(3.048), 10.0, 0.0001);
    }

    TEST(UnitTest, ConversionBetweenUnits) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        // Convert 10 feet to meters
        EXPECT_NEAR(meter.convertOneToThis(10.0, foot), 3.048, 0.0001);

        // Convert 3.048 meters to feet
        EXPECT_NEAR(foot.convertOneToThis(3.048, meter), 10.0, 0.0001);

        // Convert from this
        EXPECT_NEAR(meter.convertOneFromThis(3.048, foot), 10.0, 0.0001);
        EXPECT_NEAR(foot.convertOneFromThis(10.0, meter), 3.048, 0.0001);
    }

    TEST(UnitTest, ConversionBetweenUnitsWithOptional) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        std::optional<Unit> optFoot = foot;
        std::optional<Unit> emptyOpt;

        // With value
        EXPECT_NEAR(meter.convertOneToThis(10.0, optFoot), 3.048, 0.0001);

        // Without value (should return original)
        EXPECT_NEAR(meter.convertOneToThis(10.0, emptyOpt), 10.0, 0.0001);

        // From this with optional
        EXPECT_NEAR(meter.convertOneFromThis(3.048, optFoot), 10.0, 0.0001);
        EXPECT_NEAR(meter.convertOneFromThis(10.0, emptyOpt), 10.0, 0.0001);
    }

    TEST(UnitTest, ConversionMismatchedTypesThrows) {
        Unit linear(UnitType::linear, "meter", 1.0);
        Unit angular(UnitType::angular, "radian", 1.0);

        EXPECT_THROW(linear.convertOneToThis(5.0, angular), UnitMismatchException);
        EXPECT_THROW(linear.convertOneFromThis(5.0, angular), UnitMismatchException);

        std::optional<Unit> optAngular = angular;
        EXPECT_THROW(linear.convertOneToThis(5.0, optAngular), UnitMismatchException);
        EXPECT_THROW(linear.convertOneFromThis(5.0, optAngular), UnitMismatchException);
    }

    TEST(UnitTest, IsConsistent) {
        Unit meter1(UnitType::linear, "meter", 1.0);
        Unit meter2(UnitType::linear, "metre", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);
        Unit unknown;
        Unit angular(UnitType::angular, "radian", 1.0);

        EXPECT_TRUE(meter1.isConsistent(meter2));
        EXPECT_FALSE(meter1.isConsistent(foot));
        EXPECT_FALSE(meter1.isConsistent(angular)); // Different types

        // Unknown units are consistent with everything of the same type
        EXPECT_TRUE(unknown.isConsistent(meter1));
        EXPECT_TRUE(meter1.isConsistent(unknown));
    }

    TEST(UnitTest, SpecificUnitDefaultConstructor) {
        LinearUnit lu;
        EXPECT_EQ(lu.type(), UnitType::linear);
        EXPECT_EQ(lu.name(), "unknown");
        EXPECT_TRUE(lu.isUnknown());

        AngularUnit au;
        EXPECT_EQ(au.type(), UnitType::angular);

        AreaUnit areu;
        EXPECT_EQ(areu.type(), UnitType::area);

        VolumeUnit vu;
        EXPECT_EQ(vu.type(), UnitType::volume);

        MassUnit mu;
        EXPECT_EQ(mu.type(), UnitType::mass);
    }

    TEST(UnitTest, SpecificUnitParameterizedConstructor) {
        LinearUnit meter("meter", 1.0);
        EXPECT_EQ(meter.type(), UnitType::linear);
        EXPECT_EQ(meter.name(), "meter");
        EXPECT_FALSE(meter.isUnknown());

        LinearUnit foot("foot", 0.3048);
        EXPECT_EQ(foot.name(), "foot");
    }

    TEST(UnitTest, SpecificUnitConversion) {
        LinearUnit meter("meter", 1.0);
        LinearUnit foot("foot", 0.3048);

        EXPECT_NEAR(meter.convertOneToThis(10.0, foot), 3.048, 0.0001);
        EXPECT_NEAR(foot.convertOneToThis(3.048, meter), 10.0, 0.0001);
        EXPECT_NEAR(meter.convertOneFromThis(3.048, foot), 10.0, 0.0001);
        EXPECT_NEAR(foot.convertOneFromThis(10.0, meter), 3.048, 0.0001);
    }

    TEST(UnitTest, SpecificUnitConversionWithOptional) {
        LinearUnit meter("meter", 1.0);
        LinearUnit foot("foot", 0.3048);

        std::optional<LinearUnit> optFoot = foot;
        std::optional<LinearUnit> emptyOpt;

        EXPECT_NEAR(meter.convertOneToThis(10.0, optFoot), 3.048, 0.0001);
        EXPECT_NEAR(meter.convertOneToThis(10.0, emptyOpt), 10.0, 0.0001);
        EXPECT_NEAR(meter.convertOneFromThis(3.048, optFoot), 10.0, 0.0001);
        EXPECT_NEAR(meter.convertOneFromThis(10.0, emptyOpt), 10.0, 0.0001);
    }

    TEST(UnitTest, SpecificUnitIsConsistent) {
        LinearUnit meter1("meter", 1.0);
        LinearUnit meter2("metre", 1.0);
        LinearUnit foot("foot", 0.3048);
        LinearUnit unknown;

        EXPECT_TRUE(meter1.isConsistent(meter2));
        EXPECT_FALSE(meter1.isConsistent(foot));
        EXPECT_TRUE(unknown.isConsistent(meter1));
        EXPECT_TRUE(meter1.isConsistent(unknown));
    }

    TEST(UnitTest, UnitConverterDefaultConstructor) {
        UnitConverter uc;
        EXPECT_NEAR(uc.convertOne(5.0), 5.0, 0.0001);
        EXPECT_NEAR(uc(10.0), 10.0, 0.0001);
    }

    TEST(UnitTest, UnitConverterWithUnits) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        UnitConverter footToMeter(foot, meter);
        EXPECT_NEAR(footToMeter.convertOne(10.0), 3.048, 0.0001);
        EXPECT_NEAR(footToMeter(10.0), 3.048, 0.0001);

        UnitConverter meterToFoot(meter, foot);
        EXPECT_NEAR(meterToFoot.convertOne(3.048), 10.0, 0.0001);
    }

    TEST(UnitTest, UnitConverterWithOptionalUnits) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        std::optional<Unit> optMeter = meter;
        std::optional<Unit> optFoot = foot;
        std::optional<Unit> emptyOpt;

        UnitConverter uc1(optFoot, optMeter);
        EXPECT_NEAR(uc1.convertOne(10.0), 3.048, 0.0001);

        // One or both empty should give identity conversion
        UnitConverter uc2(emptyOpt, optMeter);
        EXPECT_NEAR(uc2.convertOne(10.0), 10.0, 0.0001);

        UnitConverter uc3(optFoot, emptyOpt);
        EXPECT_NEAR(uc3.convertOne(10.0), 10.0, 0.0001);

        UnitConverter uc4(emptyOpt, emptyOpt);
        EXPECT_NEAR(uc4.convertOne(10.0), 10.0, 0.0001);
    }

    TEST(UnitTest, UnitConverterMismatchedTypesThrows) {
        Unit linear(UnitType::linear, "meter", 1.0);
        Unit angular(UnitType::angular, "radian", 1.0);

        EXPECT_THROW(UnitConverter(linear, angular), UnitMismatchException);

        std::optional<Unit> optLinear = linear;
        std::optional<Unit> optAngular = angular;
        EXPECT_THROW(UnitConverter(optLinear, optAngular), UnitMismatchException);
    }

    TEST(UnitTest, UnitConverterConvertManyInPlace) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        UnitConverter footToMeter(foot, meter);

        std::vector<double> values = { 1.0, 2.0, 3.0, 4.0, 5.0 };
        footToMeter.convertManyInPlace(values.data(), values.size());

        EXPECT_NEAR(values[0], 0.3048, 0.0001);
        EXPECT_NEAR(values[1], 0.6096, 0.0001);
        EXPECT_NEAR(values[2], 0.9144, 0.0001);
        EXPECT_NEAR(values[3], 1.2192, 0.0001);
        EXPECT_NEAR(values[4], 1.5240, 0.0001);
    }

    TEST(UnitTest, UnitConverterConvertManyInPlaceWithSpan) {
        Unit meter(UnitType::linear, "meter", 1.0);
        Unit foot(UnitType::linear, "foot", 0.3048);

        UnitConverter footToMeter(foot, meter);

        struct Point {
            double x;
            double y;
            double z;
        };

        std::vector<Point> points = {
            {1.0, 100.0, 200.0},
            {2.0, 100.0, 200.0},
            {3.0, 100.0, 200.0}
        };

        // Convert only x values
        footToMeter.convertManyInPlace(&points[0].x, 3, sizeof(Point));

        EXPECT_NEAR(points[0].x, 0.3048, 0.0001);
        EXPECT_NEAR(points[1].x, 0.6096, 0.0001);
        EXPECT_NEAR(points[2].x, 0.9144, 0.0001);

        // y and z should be unchanged
        EXPECT_NEAR(points[0].y, 100.0, 0.0001);
        EXPECT_NEAR(points[0].z, 200.0, 0.0001);
    }

    TEST(UnitTest, SpecificUnitConverterTest) {
        LinearUnit meter("meter", 1.0);
        LinearUnit foot("foot", 0.3048);

        LinearUnitConverter footToMeter(foot, meter);
        EXPECT_NEAR(footToMeter.convertOne(10.0), 3.048, 0.0001);

        LinearUnitConverter meterToFoot(meter, foot);
        EXPECT_NEAR(meterToFoot.convertOne(3.048), 10.0, 0.0001);
    }

    TEST(UnitTest, SpecificUnitConverterWithOptional) {
        LinearUnit meter("meter", 1.0);
        LinearUnit foot("foot", 0.3048);

        std::optional<LinearUnit> optFoot = foot;
        std::optional<LinearUnit> optMeter = meter;
        std::optional<LinearUnit> emptyOpt;

        LinearUnitConverter uc1(optFoot, optMeter);
        EXPECT_NEAR(uc1.convertOne(10.0), 3.048, 0.0001);

        LinearUnitConverter uc2(emptyOpt, optMeter);
        EXPECT_NEAR(uc2.convertOne(10.0), 10.0, 0.0001);
    }

    TEST(UnitTest, ComplexConversionChain) {
        using namespace linearUnitPresets;

        // Convert from feet to inches to meters
        LinearUnit inch("inch", 0.0254);

        LinearUnitConverter footToInch(internationalFoot, inch);
        double valueInInches = footToInch.convertOne(1.0);
        EXPECT_NEAR(valueInInches, 12.0, 0.0001);

        LinearUnitConverter inchToMeter(inch, meter);
        double valueInMeters = inchToMeter.convertOne(valueInInches);
        EXPECT_NEAR(valueInMeters, 0.3048, 0.0001);
    }
}