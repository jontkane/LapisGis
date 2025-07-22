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
}