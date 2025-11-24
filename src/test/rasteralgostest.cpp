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

	TEST_F(RasterAlgosTest, EuclideanDistanceTransform) {
		//a small 5x5 raster with two disjoint true cells
		//cellsize is non-1 to test conversion from pixel units to real units
		Raster<int> r(Alignment(0, 0, 5, 5, 0.5, 0.5));

		auto predicate = [](const int& v) {return v == 9; }; //an unusual predicate to make sure it's applied properly



		//fill the raster with has_value=true, value=anything but 9
		for (rowcol_t row = 0; row < r.nrow(); ++row) {
			for (rowcol_t col = 0; col < r.ncol(); ++col) {
				r.atRCUnsafe(row, col).has_value() = true;
				int value = row + col;
				if (value == 9) {
					value = 8;
				}
				r.atRCUnsafe(row, col).value() = value;
			}
		}

		//0=false according to predicate
		//1=true according to predicate
		/*
		* 0, 1, 0, 0, 0
		* NA,0, 0, 0, 0
		* 0, 0, 0, 1, 0
		* 0, NA,0, 0, 0
		* 0, 0, 0, 0, 0	
		*/
		//set some specific cells to things to test:
		//nodata which is 9, nodata which is not 9, and valid data which is 9
		r.atRCUnsafe(0, 1).value() = 9; //true cell
		r.atRCUnsafe(1, 0).has_value() = false; //nodata
		r.atRCUnsafe(2, 3).value() = 9; //true cell
		r.atRCUnsafe(3, 1).has_value() = false; //nodata, value = 9
		r.atRCUnsafe(3, 1).value() = 9;

		//expected (in pixel units):
		/*
		* 1,0,1,2,sqrt(5)
		* NA,1,sqrt(2),1,sqrt(2)
		* sqrt(5),2,1,0,1
		* sqrt(10),NA,sqrt(2),1,sqrt(2)
		* sqrt(13),sqrt(8),sqrt(5),2,sqrt(5)
		*/
		auto out = euclideanDistanceTransform(r, predicate);
		ASSERT_EQ(out.nrow(), r.nrow());
		ASSERT_EQ(out.ncol(), r.ncol());

		Raster<coord_t> expected((Alignment)r);
		for (cell_t cell : CellIterator(expected)) {
			expected.atCellUnsafe(cell).has_value() = true;
		}
		expected.atRCUnsafe(0, 0).value() = 1.0;
		expected.atRCUnsafe(0, 1).value() = 0.0;
        expected.atRCUnsafe(0, 2).value() = 1.0;
        expected.atRCUnsafe(0, 3).value() = 2.0;
        expected.atRCUnsafe(0, 4).value() = std::sqrt(5.0);
        
		expected.atRCUnsafe(1, 0).has_value() = false;
		expected.atRCUnsafe(1,1).value() = 1.0;
		expected.atRCUnsafe(1, 2).value() = std::sqrt(2.0);
		expected.atRCUnsafe(1, 3).value() = 1.0;
        expected.atRCUnsafe(1, 4).value() = std::sqrt(2.0);

        expected.atRCUnsafe(2, 0).value() = std::sqrt(5.0);
        expected.atRCUnsafe(2, 1).value() = 2.0;
        expected.atRCUnsafe(2, 2).value() = 1.0;
        expected.atRCUnsafe(2, 3).value() = 0.0;
        expected.atRCUnsafe(2, 4).value() = 1.0;

        expected.atRCUnsafe(3, 0).value() = std::sqrt(10.0);
        expected.atRCUnsafe(3, 1).has_value() = false;
        expected.atRCUnsafe(3, 2).value() = std::sqrt(2.0);
		expected.atRCUnsafe(3, 3).value() = 1.0;
        expected.atRCUnsafe(3, 4).value() = std::sqrt(2.0);

        expected.atRCUnsafe(4, 0).value() = std::sqrt(13.0);
        expected.atRCUnsafe(4, 1).value() = std::sqrt(8.0);
        expected.atRCUnsafe(4, 2).value() = std::sqrt(5.0);
        expected.atRCUnsafe(4, 3).value() = 2.0;
        expected.atRCUnsafe(4, 4).value() = std::sqrt(5.0);

        //convert expected from pixel units to real units (cellsize = 0.5, so multiply by 0.5)
		expected = expected * 0.5;

		for (cell_t cell : CellIterator(out)) {
			if (!expected.atCellUnsafe(cell).has_value()) {
				EXPECT_FALSE(out.atCellUnsafe(cell).has_value()) << " at cell " << cell;
			}
			else {
				ASSERT_TRUE(out.atCellUnsafe(cell).has_value()) << " at cell " << cell;
				EXPECT_NEAR(out.atCellUnsafe(cell).value(), expected.atCellUnsafe(cell).value(), 0.0001) << " at cell " << cell;
			}
        }
	}

	TEST_F(RasterAlgosTest, ConnectedComponentsTest) {
		Raster<int> r{ Alignment(Extent(0,5,0,5),5,5) };
		/*
		*  NA  1  NA  2  2
		*   1  1  NA  NA  2
		*  NA NA  NA  5 NA
		*   3  6   4  4  5
		*   3  3  NA NA  NA
		*/
        r.atRCUnsafe(0, 0).has_value() = false;
        r.atRCUnsafe(0, 1).value() = 1; r.atRCUnsafe(0, 1).has_value() = true;
        r.atRCUnsafe(0, 2).has_value() = false;
        r.atRCUnsafe(0, 3).value() = 2; r.atRCUnsafe(0, 3).has_value() = true;
		r.atRCUnsafe(0, 4).value() = 2; r.atRCUnsafe(0, 4).has_value() = true;
        r.atRCUnsafe(1, 0).value() = 1; r.atRCUnsafe(1, 0).has_value() = true;
        r.atRCUnsafe(1, 1).value() = 1; r.atRCUnsafe(1, 1).has_value() = true;
        r.atRCUnsafe(1, 2).has_value() = false;
        r.atRCUnsafe(1, 3).has_value() = false;
        r.atRCUnsafe(1, 4).value() = 2; r.atRCUnsafe(1, 4).has_value() = true;
        r.atRCUnsafe(2, 0).has_value() = false;
        r.atRCUnsafe(2, 1).has_value() = false;
        r.atRCUnsafe(2, 2).has_value() = false;
        r.atRCUnsafe(2, 3).value() = 5; r.atRCUnsafe(2, 3).has_value() = true;
        r.atRCUnsafe(2, 4).has_value() = false;
        r.atRCUnsafe(3, 0).value() = 3; r.atRCUnsafe(3, 0).has_value() = true;
        r.atRCUnsafe(3, 1).value() = 6; r.atRCUnsafe(3, 1).has_value() = true;
        r.atRCUnsafe(3, 2).value() = 4; r.atRCUnsafe(3, 2).has_value() = true;
        r.atRCUnsafe(3, 3).value() = 4; r.atRCUnsafe(3, 3).has_value() = true;
		r.atRCUnsafe(3, 4).value() = 5; r.atRCUnsafe(3, 4).has_value() = true;
		r.atRCUnsafe(4, 0).value() = 3; r.atRCUnsafe(4, 0).has_value() = true;
		r.atRCUnsafe(4, 1).value() = 3; r.atRCUnsafe(4, 1).has_value() = true;
		r.atRCUnsafe(4, 2).has_value() = false;
		r.atRCUnsafe(4, 3).has_value() = false;
        r.atRCUnsafe(4, 4).has_value() = false;

        Raster<cell_t> withDiagonals = connectedComponents(r, true);
        Raster<cell_t> withoutDiagonals = connectedComponents(r, false);

		ASSERT_EQ((Alignment)r, (Alignment(withDiagonals)));
        ASSERT_EQ((Alignment)r, (Alignment(withoutDiagonals)));

        //check NAs are in the right places, first
		for (cell_t cell : CellIterator(r)) {
            EXPECT_EQ(r.atCellUnsafe(cell).has_value(), withDiagonals.atCellUnsafe(cell).has_value()) << " at cell " << cell;
            EXPECT_EQ(r.atCellUnsafe(cell).has_value(), withoutDiagonals.atCellUnsafe(cell).has_value()) << " at cell " << cell;
        }

		//we can't expect any specific labels, only that each component has a single label
        cell_t label1 = withDiagonals.atRCUnsafe(0, 1).value();
        EXPECT_EQ(label1, withDiagonals.atRCUnsafe(1, 0).value());
        EXPECT_EQ(label1, withDiagonals.atRCUnsafe(1, 1).value());

		cell_t label2 = withDiagonals.atRCUnsafe(0, 3).value();
        EXPECT_EQ(label2, withDiagonals.atRCUnsafe(0, 4).value());
        EXPECT_EQ(label2, withDiagonals.atRCUnsafe(1, 4).value());

		cell_t label3 = withDiagonals.atRCUnsafe(3, 0).value();
        EXPECT_EQ(label3, withDiagonals.atRCUnsafe(4, 0).value());
        EXPECT_EQ(label3, withDiagonals.atRCUnsafe(4, 1).value());

		cell_t label4 = withDiagonals.atRCUnsafe(3, 2).value();
        EXPECT_EQ(label4, withDiagonals.atRCUnsafe(3, 3).value());

        cell_t label5 = withDiagonals.atRCUnsafe(2, 3).value();
        EXPECT_EQ(label5, withDiagonals.atRCUnsafe(3, 4).value());

		cell_t label6 = withDiagonals.atRCUnsafe(3, 1).value();

		std::unordered_set<cell_t > labelsWithDiagonals{ label1, label2, label3, label4, label5, label6 };
        EXPECT_EQ(labelsWithDiagonals.size(), 6);

        //without diagonals, some the 5s should be split into two labels
        label1 = withoutDiagonals.atRCUnsafe(0, 1).value();
        EXPECT_EQ(label1, withoutDiagonals.atRCUnsafe(1, 0).value());
        EXPECT_EQ(label1, withoutDiagonals.atRCUnsafe(1, 1).value());

		label2 = withoutDiagonals.atRCUnsafe(0, 3).value();
		EXPECT_EQ(label2, withoutDiagonals.atRCUnsafe(0, 4).value());
		EXPECT_EQ(label2, withoutDiagonals.atRCUnsafe(1, 4).value());

		label3 = withoutDiagonals.atRCUnsafe(3, 0).value();
		EXPECT_EQ(label3, withoutDiagonals.atRCUnsafe(4, 0).value());
		EXPECT_EQ(label3, withoutDiagonals.atRCUnsafe(4, 1).value());

		label4 = withoutDiagonals.atRCUnsafe(3, 2).value();
		EXPECT_EQ(label4, withoutDiagonals.atRCUnsafe(3, 3).value());

		label5 = withoutDiagonals.atRCUnsafe(2, 3).value();
        cell_t label5b = withoutDiagonals.atRCUnsafe(3, 4).value();

		label6 = withoutDiagonals.atRCUnsafe(3, 1).value();

		std::unordered_set<cell_t > labelsWithoutDiagonals{ label1, label2, label3, label4, label5, label5b, label6 };
        EXPECT_EQ(labelsWithoutDiagonals.size(), 7);

	}

	TEST_F(RasterAlgosTest, Nibble) {
		/* inraster. NA cells should remain NA in the output; other cells should either remain or change based on the mask raster
        * NA  1  2  3  4
		* NA  1  2  3  4
		* NA  1  2  3  4
		* NA  1  2  3  4
		* NA  1  2  3  4
		*/

        /* maskraster. NA = nibble (fill from neighbors); 1 = keep original value
		* NA NA 1  1  1
        * NA NA 1  1  1
		* 1  1  NA NA NA
        * 1  1  NA NA NA
		* 1  1  NA NA NA
		*/
		Raster<int> inraster{ Alignment(Extent(0,5,0,5),5,5) };
		Raster<int> maskraster{ Alignment(Extent(0,5,0,5),5,5) };
		for (rowcol_t row = 0; row < inraster.nrow(); ++row) {
			for (rowcol_t col = 0; col < inraster.ncol(); ++col) {
                auto inv = inraster.atRCUnsafe(row, col);
                auto maskv = maskraster.atRCUnsafe(row, col);
				if (col == 0) {
					inv.has_value() = false;
				}
				else {
					inv.has_value() = true;
					inv.value() = col;
				}

				if ( (col <= 1) && (row <= 1) ) {
					//upper left square of NAs
					maskv.has_value() = false;
				}
				else {
					if (row <= 1) {
						//upper right square of values
                        maskv.has_value() = true;
						maskv.value() = row + col; //semi-random value because it shouldn't matter
					}
					else {
                        //lower left square of values
						if (col <= 1) {
							maskv.has_value() = true;
							maskv.value() = row + col; //semi-random value because it shouldn't matter
						}
						else {
							//lower right square of NAs
							maskv.has_value() = false;
                        }
					}
                }
			}
		}

        Raster<int> outraster = nibble(inraster, maskraster);

		for (rowcol_t row = 0; row < inraster.nrow(); ++row) {
			for (rowcol_t col = 0; col < inraster.ncol(); ++col) {
				//NAs in input should stay NA
				if (!inraster.atRCUnsafe(row, col).has_value()) {
					EXPECT_FALSE(outraster.atRCUnsafe(row, col).has_value()) << " at row " << row << ", col " << col;
					continue;
				}
				//non-NAs in input
				if (maskraster.atRCUnsafe(row, col).has_value()) {
					//mask says keep original
					EXPECT_TRUE(outraster.atRCUnsafe(row, col).has_value()) << " at row " << row << ", col " << col;
					EXPECT_EQ(outraster.atRCUnsafe(row, col).value(), inraster.atRCUnsafe(row, col).value()) << " at row " << row << ", col " << col;
					continue;
				}
				//mask says nibble (fill from neighbors)
				EXPECT_TRUE(outraster.atRCUnsafe(row, col).has_value()) << " at row " << row << ", col " << col;

				//the algorithm used doesn't guarantee which neighbor it fills from
				if (col <= 1) {
					//upper-left corner, valid fills are 1 or 2
                    int val = outraster.atRCUnsafe(row, col).value();
                    EXPECT_TRUE((val == 1) || (val == 2)) << " at row " << row << ", col " << col;
				}
				else {
					//lower-right corner, 1-4 are all valid fills
					int val = outraster.atRCUnsafe(row, col).value();
                    EXPECT_TRUE((val >= 1) && (val <= 4)) << " at row " << row << ", col " << col;
				}
			}
		}
		

	}

	TEST_F(RasterAlgosTest, MosaicInside) {
        //two 3x3 rasters, with two columns of overlap. Values are numbered based on the column number they expect to end up in in the mosaic. Second raster is that value +10 so I can track which raster the value came from.
		/*
		* 0 1 2
		* 0 1 2
		* 0 1 2
		* 
		* 11 12 13
		* 11 12 13
		* 11 12 13
		*/
		//top row of the second raster is NA, to make sure that's handled properly


        Raster<int> r1{ Alignment(Extent(0,3,0,3),3,3) };
        Raster<int> r2{ Alignment(Extent(1,4,0,3),3,3) };
		for (rowcol_t row = 0; row < r1.nrow(); ++row) {
			for (rowcol_t col = 0; col < r1.ncol(); ++col) {
				r1.atRCUnsafe(row, col).has_value() = true;
				r1.atRCUnsafe(row, col).value() = (int)col;
			}
        }
		for (rowcol_t row = 0; row < r1.nrow(); ++row) {
			for (rowcol_t col = 0 ; col < r1.ncol(); ++col) {
				if (row == 0) {
					r2.atRCUnsafe(row, col).has_value() = false;
				}
				else {
					r2.atRCUnsafe(row, col).has_value() = true;
					r2.atRCUnsafe(row, col).value() = 11 + (int)col;
				}
            }
		}

		/* expected output:
		* 0 1 2 NA
		* 0 1 12 13
		* 0 1/11 2/12 13
		*/

		std::vector<Raster<int>*> rasters{ &r1, &r2 };
		auto out = mosaicInside(rasters);
		Alignment expectedAlign = extendAlignment(r1, r2, SnapType::near);
		EXPECT_TRUE(out.isSameAlignment(expectedAlign));
		EXPECT_EQ(out.ncell(), expectedAlign.ncell());
		for (rowcol_t row = 0; row < out.nrow(); ++row) {
			for (rowcol_t col = 0; col < out.ncol(); ++col) {
                auto outv = out.atRCUnsafe(row, col);
				if ( (row == 0) && (col == 3) ) {
					EXPECT_FALSE(outv.has_value()) << " at row " << row << ", col " << col;
				}
				if ((row == 0) && (col == 2)) {
                    EXPECT_TRUE(outv.has_value()) << " at row " << row << ", col " << col;
                    EXPECT_EQ(outv.value(), 2) << " at row " << row << ", col " << col;
				}
				if (col <= 1) {
					EXPECT_TRUE(outv.has_value()) << " at row " << row << ", col " << col;
					EXPECT_EQ(outv.value(), (int)col) << " at row " << row << ", col " << col;
				}
				if (col >= 2 && row == 1) {
					EXPECT_TRUE(outv.has_value()) << " at row " << row << ", col " << col;
                    EXPECT_EQ(outv.value(), 10 + (int)col) << " at row " << row << ", col " << col;
				}
				if (col == 0 && row == 3) {
					EXPECT_TRUE(outv.has_value()) << " at row " << row << ", col " << col;
					EXPECT_EQ(outv.value(), 0) << " at row " << row << ", col " << col;
				}
				if (col == 3 && row == 3) {
					EXPECT_TRUE(outv.has_value()) << " at row " << row << ", col " << col;
                    EXPECT_EQ(outv.value(), 13) << " at row " << row << ", col " << col;
				}
				if (col >= 1 && col <= 2 && row == 3) {
					EXPECT_TRUE(outv.has_value()) << " at row " << row << ", col " << col;
					int val = outv.value();
					EXPECT_TRUE((val == col) || (val == 10 + col)) << " at row " << row << ", col " << col;
				}
			}
		}
	}
}