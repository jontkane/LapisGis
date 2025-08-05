#include"test_pch.hpp"
#include"..\Vector.hpp"

namespace lapis {
    TEST(VectorTest, pointsConstructor) {
        std::string file = LAPISGISTESTFILES;
        file += "/testpoints.shp";
        VectorDataset<Point> p{ file };

        ASSERT_EQ(p.nFeature(), 3);
        ASSERT_EQ(p.getAllFieldNames().size(), 3);
        ASSERT_EQ(p.getAllFieldNames()[0], "integer");
        ASSERT_EQ(p.getAllFieldNames()[1], "double");
        ASSERT_EQ(p.getAllFieldNames()[2], "string");
        ASSERT_EQ(p.getFieldType(p.getAllFieldNames()[0]), FieldType::Integer);
        ASSERT_EQ(p.getFieldType(p.getAllFieldNames()[1]), FieldType::Real);
        ASSERT_EQ(p.getFieldType(p.getAllFieldNames()[2]), FieldType::String);

        int i = 0;
        std::vector<std::string> names = { "one","two","three" };
        for (auto feature : p) {
            ++i;
            EXPECT_EQ(feature.getGeometry().x(), i);
            EXPECT_EQ(feature.getGeometry().y(), i);
            EXPECT_EQ(feature.getIntegerField("integer"), i);
            EXPECT_EQ(feature.getRealField("double"), i);
            EXPECT_STREQ(feature.getStringField("string").c_str(), names[i - 1].c_str()); //the attribute table pads with \0 which will make std::string::operator== unreliable
        }

        EXPECT_TRUE(p.crs().isConsistent("2285"));
    }

    TEST(VectorTest, PolygonArea) {
        // Square: (0,0), (0,1), (1,1), (1,0)
        std::vector<CoordXY> square = { {0,0}, {0,1}, {1,1}, {1,0} };
        Polygon poly(square);
        EXPECT_NEAR(poly.area(), 1.0, 0.01);

        // Triangle: (0,0), (1,0), (0,1)
        std::vector<CoordXY> triangle = { {0,0}, {1,0}, {0,1} };
        Polygon tri(triangle);
        EXPECT_NEAR(tri.area(), 0.5, 0.01);

        // Rectangle: (0,0), (0,2), (3,2), (3,0)
        std::vector<CoordXY> rect = { {0,0}, {0,2}, {3,2}, {3,0} };
        Polygon rectangle(rect);
        EXPECT_NEAR(rectangle.area(), 6.0, 0.01);

        // Polygon with a hole (outer: square, inner: smaller square)
        std::vector<CoordXY> outer = { {0,0}, {0,4}, {4,4}, {4,0} };
        std::vector<CoordXY> inner = { {1,1}, {1,2}, {2,2}, {2,1} };
        Polygon polyWithHole(outer);
        polyWithHole.addInnerRing(inner);
        EXPECT_NEAR(polyWithHole.area(), 16.0 - 1.0, 0.01);
    }

    TEST(VectorTest, MultiPolygonArea) {
        // Two non-overlapping squares
        std::vector<CoordXY> square1 = { {0,0}, {0,1}, {1,1}, {1,0} };
        std::vector<CoordXY> square2 = { {2,2}, {2,3}, {3,3}, {3,2} };
        Polygon poly1(square1);
        Polygon poly2(square2);

        MultiPolygon mp;
        mp.addPolygon(poly1);
        mp.addPolygon(poly2);

        EXPECT_NEAR(mp.area(), 2.0,0.01);

        // MultiPolygon with a polygon with a hole
        std::vector<CoordXY> outer = { {0,0}, {0,4}, {4,4}, {4,0} };
        std::vector<CoordXY> inner = { {1,1}, {1,2}, {2,2}, {2,1} };
        Polygon polyWithHole(outer);
        polyWithHole.addInnerRing(inner);

        MultiPolygon mp2;
        mp2.addPolygon(polyWithHole);
        mp2.addPolygon(poly1);

        EXPECT_NEAR(mp2.area(), (16.0 - 1.0) + 1.0, 0.01);
    }
}