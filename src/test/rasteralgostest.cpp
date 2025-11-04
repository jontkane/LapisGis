#include"test_pch.hpp"
#include"..\rasteralgos.hpp"

namespace lapis {

	class RasterAlgosTest : public ::testing::Test {
	public:

		Raster<double> r;
		Alignment a;

		void SetUp() override {
			a = Alignment(Extent(0, 4, 0, 4), 2, 2);
			r = Raster<double>(Alignment(Extent(0, 4, 0, 4), 4, 4));
			for (cell_t cell = 0; cell < r.ncell(); ++cell) {
				r[cell].value() = (int)cell;
				r[cell].has_value() = true;
			}
			r[0].has_value() = false;
			r[1].has_value() = false;
			r[4].has_value() = false;
			r[5].has_value() = false;
			r[15].has_value() = false;
		}
	};
	
	TEST_F(RasterAlgosTest, AggregateMax) {
		ViewFunc<double, double> f{ &viewMax<double, double> };
		auto out = aggregate(r, a, f);

		EXPECT_EQ((Alignment)out, a);
		
		EXPECT_FALSE(out[0].has_value());
		EXPECT_TRUE(out[1].has_value());
		EXPECT_TRUE(out[2].has_value());
		EXPECT_TRUE(out[3].has_value());

		EXPECT_EQ(out[1].value(), 7);
		EXPECT_EQ(out[2].value(), 13);
		EXPECT_EQ(out[3].value(), 14);
	}

	TEST_F(RasterAlgosTest, AggregateMean) {
		ViewFunc<double, double> f{ &viewMean<double, double> };
		auto out = aggregate(r, a, f);

		EXPECT_EQ((Alignment)out, a);

		EXPECT_FALSE(out[0].has_value());
		EXPECT_TRUE(out[1].has_value());
		EXPECT_TRUE(out[2].has_value());
		EXPECT_TRUE(out[3].has_value());

		EXPECT_NEAR(out[1].value(), 4.5, 0.01);
		EXPECT_NEAR(out[2].value(), 10.5, 0.01);
		EXPECT_NEAR(out[3].value(), 11.7, 0.1);
	}

	TEST_F(RasterAlgosTest, AggregateStdDev) {
		ViewFunc<double, double> f{ &viewStdDev<double, double> };
		auto out = aggregate(r, a, f);

		EXPECT_EQ((Alignment)out, a);

		EXPECT_FALSE(out[0].has_value());
		EXPECT_TRUE(out[1].has_value());
		EXPECT_TRUE(out[2].has_value());
		EXPECT_TRUE(out[3].has_value());

		EXPECT_NEAR(out[1].value(), 2.06, 0.01);
		EXPECT_NEAR(out[2].value(), 2.06, 0.01);
		EXPECT_NEAR(out[3].value(), 1.7, 0.01);
	}

	TEST_F(RasterAlgosTest, Mosaic) {
		// 1. Output alignment is of expected size
		Alignment a1(Extent(0, 4, 0, 4), 2, 2); // 2x2
		Alignment a2(Extent(2, 6, 2, 6), 2, 2); // 2x2, offset
		Raster<double> r1(a1);
		Raster<double> r2(a2);
		for (cell_t cell = 0; cell < r1.ncell(); ++cell) {
			r1[cell].value() = 10 + (double)cell;
			r1[cell].has_value() = true;
		}
		for (cell_t cell = 0; cell < r2.ncell(); ++cell) {
			r2[cell].value() = 20 + (double)cell;
			r2[cell].has_value() = true;
		}
		// Set some cells to missing
		r1[0].has_value() = false;
		r2[3].has_value() = false;
		std::vector<Raster<double>*> rasters{ &r1, &r2 };
		// Use maxCombiner for overlap
		auto out = mosaic(rasters, maxCombiner<double>);
		Alignment expectedAlign = extendAlignment(a1, a2, SnapType::near);
		EXPECT_TRUE(out.isSameAlignment(expectedAlign));
		EXPECT_EQ(out.ncell(), expectedAlign.ncell());
		// 2. Unaligned rasters should throw
		Alignment badAlign(Extent(0, 4, 0, 4), 3, 3);
		Raster<double> badRaster(badAlign);
		badRaster[0].value() = 99;
		badRaster[0].has_value() = true;
		std::vector<Raster<double>*> badRasters{ &r1, &badRaster };
		EXPECT_THROW(mosaic(badRasters, maxCombiner<double>), AlignmentMismatchException);
		// 3. Output has expected values
		for (cell_t cell = 0; cell < out.ncell(); ++cell) {
			coord_t x = out.xFromCell(cell);
			coord_t y = out.yFromCell(cell);
			bool inR1 = a1.contains(x, y);
			bool inR2 = a2.contains(x, y);
			bool r1Has = inR1 ? r1[r1.cellFromXY(x, y)].has_value() : false;
			bool r2Has = inR2 ? r2[r2.cellFromXY(x, y)].has_value() : false;
			if (r1Has && r2Has) {
				cell_t c1 = r1.cellFromXY(x, y);
				cell_t c2 = r2.cellFromXY(x, y);
				EXPECT_EQ(out[cell].value(), std::max(r1[c1].value(), r2[c2].value()));
			}
			else if (r1Has) {
				cell_t c1 = r1.cellFromXY(x, y);
				EXPECT_EQ(out[cell].value(), r1[c1].value());
			}
			else if (r2Has) {
				cell_t c2 = r2.cellFromXY(x, y);
				EXPECT_EQ(out[cell].value(), r2[c2].value());
			}
			else if (inR1 || inR2) {
				EXPECT_FALSE(out[cell].has_value());
			}
			else {
				EXPECT_FALSE(out[cell].has_value());
			}
		}
		// 4. Order-sensitive combiner
		auto firstCombiner = [](xtl::xoptional<double> a, xtl::xoptional<double> b) {
			if (a.has_value()) return a;
			return b;
			};
		auto out1 = mosaic(rasters, firstCombiner);
        std::vector<Raster<double>*> rastersRev{ &r2, &r1 };
		auto out2 = mosaic(rastersRev, firstCombiner);
		for (cell_t cell = 0; cell < out.ncell(); ++cell) {
			coord_t x = out.xFromCell(cell);
			coord_t y = out.yFromCell(cell);
			bool inR1 = a1.contains(x, y);
			bool inR2 = a2.contains(x, y);
			bool r1Has = inR1 ? r1[r1.cellFromXY(x, y)].has_value() : false;
			bool r2Has = inR2 ? r2[r2.cellFromXY(x, y)].has_value() : false;
			if (r1Has) {
				cell_t c1 = r1.cellFromXY(x, y);
				EXPECT_EQ(out1[cell].value(), r1[c1].value());
			}
			else if (r2Has) {
				cell_t c2 = r2.cellFromXY(x, y);
				EXPECT_EQ(out1[cell].value(), r2[c2].value());
			}
			else {
				EXPECT_FALSE(out1[cell].has_value());
			}
			if (r2Has) {
				cell_t c2 = r2.cellFromXY(x, y);
				EXPECT_EQ(out2[cell].value(), r2[c2].value());
			}
			else if (r1Has) {
				cell_t c1 = r1.cellFromXY(x, y);
				EXPECT_EQ(out2[cell].value(), r1[c1].value());
			}
			else {
				EXPECT_FALSE(out2[cell].has_value());
			}
		}
	}
}