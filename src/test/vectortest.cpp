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

    TEST(VectorTest, Append) {
        std::string file1 = LAPISGISTESTFILES;
        file1 += "/testpoints.shp";
        VectorDataset<Point> p1{ file1 };
        VectorDataset<Point> append{ file1 };

        std::string file2 = LAPISGISTESTFILES;
        file2 += "/testpoints.shp";
        VectorDataset<Point> p2{ file2 };
        size_t n1 = p1.nFeature();
        size_t n2 = p2.nFeature();
        append.appendFile(file2);
        ASSERT_EQ(append.nFeature(), n1 + n2);

        for (size_t i = 0; i < n1; ++i) {
            auto expected = p1.getFeature(i);
            auto feature1 = append.getFeature(i);
            auto feature2 = append.getFeature(i + n1);

            EXPECT_EQ(feature1.getGeometry().x(), expected.getGeometry().x());
            EXPECT_EQ(feature1.getGeometry().y(), expected.getGeometry().y());
            EXPECT_EQ(feature1.getIntegerField("integer"), expected.getIntegerField("integer"));
            EXPECT_EQ(feature1.getRealField("double"), expected.getRealField("double"));
            EXPECT_EQ(feature1.getStringField("string"), expected.getStringField("string"));
            
            EXPECT_EQ(feature2.getGeometry().x(), expected.getGeometry().x());
            EXPECT_EQ(feature2.getGeometry().y(), expected.getGeometry().y());
            EXPECT_EQ(feature2.getIntegerField("integer"), expected.getIntegerField("integer"));
            EXPECT_EQ(feature2.getRealField("double"), expected.getRealField("double"));
            EXPECT_EQ(feature2.getStringField("string"), expected.getStringField("string"));

        }
    }

    TEST(VectorTest, MultiFileConstructor) {
        std::string file1 = LAPISGISTESTFILES;
        file1 += "/testpoints.shp";
        std::string file2 = LAPISGISTESTFILES;
        file2 += "/testpoints.shp";
        std::vector<std::filesystem::path> files = { file1, file2 };
        VectorDataset<Point> p{ files };
        //expecting this to be equivalent to appending file2 to file1
        VectorDataset<Point> expected{ file1 };
        expected.appendFile(file2);
        ASSERT_EQ(p.nFeature(), expected.nFeature());
        for (size_t i = 0; i < p.nFeature(); ++i) {
            auto featureP = p.getFeature(i);
            auto featureExpected = expected.getFeature(i);
            EXPECT_EQ(featureP.getGeometry().x(), featureExpected.getGeometry().x());
            EXPECT_EQ(featureP.getGeometry().y(), featureExpected.getGeometry().y());
            EXPECT_EQ(featureP.getIntegerField("integer"), featureExpected.getIntegerField("integer"));
            EXPECT_EQ(featureP.getRealField("double"), featureExpected.getRealField("double"));
            EXPECT_EQ(featureP.getStringField("string"), featureExpected.getStringField("string"));
        }
    }

    TEST(VectorTest, AddFeature) {
        std::string file = LAPISGISTESTFILES;
        file += "/testpoints.shp";
        VectorDataset<Point> p{ file };
        size_t nOriginal = p.nFeature();
        auto featureToAdd = p.getFeature(0);
        p.addFeature(featureToAdd);
        ASSERT_EQ(p.nFeature(), nOriginal + 1);
        auto addedFeature = p.getFeature(nOriginal);
        featureToAdd = p.getFeature(0); //we invalidated the iterators by adding an element
        EXPECT_EQ(addedFeature.getGeometry().x(), featureToAdd.getGeometry().x());
        EXPECT_EQ(addedFeature.getGeometry().y(), featureToAdd.getGeometry().y());
        EXPECT_EQ(addedFeature.getIntegerField("integer"), featureToAdd.getIntegerField("integer"));
        EXPECT_EQ(addedFeature.getRealField("double"), featureToAdd.getRealField("double"));
        EXPECT_EQ(addedFeature.getStringField("string"), featureToAdd.getStringField("string"));
    }

    TEST(VectorTest, EmptyDatasetFromTemplate) {
        std::string exampleFile = LAPISGISTESTFILES;
        exampleFile += "/testpoints.shp";
        VectorDataset<Point> p{ exampleFile };

        VectorDataset<Point> fromTemplate = emptyVectorDatasetFromTemplate<Point>(p);
        ASSERT_EQ(fromTemplate.nFeature(), 0);
        ASSERT_TRUE(fromTemplate.crs().isConsistent(p.crs()));
        ASSERT_EQ(fromTemplate.getAllFieldNames().size(), p.getAllFieldNames().size());
        for (const auto& fieldName : p.getAllFieldNames()) {
            ASSERT_EQ(fromTemplate.getFieldType(fieldName), p.getFieldType(fieldName));
            if (p.getFieldType(fieldName) == FieldType::String) {
                ASSERT_EQ(fromTemplate.getStringFieldWidth(fieldName), p.getStringFieldWidth(fieldName));
            }
        }
    }
}