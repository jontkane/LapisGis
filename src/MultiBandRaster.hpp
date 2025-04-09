#pragma once
#ifndef LP_MULTIBANDRASTER_H
#define LP_MULTIBANDRASTER_H

#include"Raster.hpp"


namespace lapis {

	//returns the number of bands the file has, without reading all the data in it.
	//returns 0 if the file is not a valid raster
	band_t nBandsForFile(const std::string& file);

	template<typename T>
	class MultiBandRaster : public Alignment {
	public:
		MultiBandRaster(const std::string& filename);
		MultiBandRaster(const std::string& filename, const Extent& e, SnapType snap);
		MultiBandRaster(const Alignment& a, band_t nBands);

		band_t nBands() const;

		Raster<T>& bandAt(band_t band);
		Raster<T>& bandAtUnsafe(band_t band);

		const Raster<T>& bandAt(band_t band) const;
		const Raster<T>& bandAtUnsafe(band_t band) const;

		auto atCell(cell_t cell, band_t band);
		auto atCellUnsafe(cell_t cell, band_t band);
		auto atRC(rowcol_t row, rowcol_t col, band_t band);
		auto atRCUnsafe(rowcol_t row, rowcol_t col, band_t band);
		auto atXY(coord_t x, coord_t y, band_t band);
		auto atXYUnsafe(coord_t x, coord_t y, band_t band);

		const auto atCell(cell_t cell, band_t band) const;
		const auto atCellUnsafe(cell_t cell, band_t band) const;
		const auto atRC(rowcol_t row, rowcol_t col, band_t band) const;
		const auto atRCUnsafe(rowcol_t row, rowcol_t col, band_t band) const;
		const auto atXY(coord_t x, coord_t y, band_t band) const;
		const auto atXYUnsafe(coord_t x, coord_t y, band_t band) const;

		void writeRaster(const std::string& fileName, const std::string driver = "GTiff", const T navalue = std::numeric_limits<T>::lowest(), GDALDataType gdt = GDT_Unknown);
	private:
		std::vector<Raster<T>> _bands;

		void _checkBand(band_t band) const;
	};
	template<typename T>
	inline MultiBandRaster<T>::MultiBandRaster(const std::string& filename) {
		using namespace lapis;
		UniqueGdalDataset wgd = rasterGDALWrapper(filename);
		if (!wgd) {
			throw InvalidRasterFileException("Unable to open " + filename + " as a raster");
		}
		alignmentInitFromGDALRaster(wgd, getGeoTrans(wgd, filename));
		checkValidAlignment();

		band_t nBand = wgd->GetRasterCount();
		_bands.resize(nBand);
		for (band_t i = 1; i <= nBand; ++i) {
			_bands[i - 1] = Raster<T>(filename, i);
		}
	}
	template<typename T>
	inline MultiBandRaster<T>::MultiBandRaster(const std::string& filename, const Extent& e, SnapType snap)
	{
		UniqueGdalDataset wgd = rasterGDALWrapper(filename);
		if (!wgd) {
			throw InvalidRasterFileException("Unable to open " + filename + " as a raster");
		}
		Alignment a = Alignment(filename);
		Extent snapE = a.alignExtent(e, snap);
		snapE = cropExtent(snapE, a);
		auto rc = a.rowColExtent(snapE, snap);

		_crs = a.crs();
		_xmin = snapE.xmin();
		_xmax = snapE.xmax();
		_ymin = snapE.ymin();
		_ymax = snapE.ymax();
		_nrow = rc.maxrow - rc.minrow + 1;
		_ncol = rc.maxcol - rc.mincol + 1;
		_xres = a.xres();
		_yres = a.yres();
		checkValidAlignment();

		band_t nBand = wgd->GetRasterCount();
		_bands.resize(nBand);
		for (band_t i = 1; i <= nBand; ++i) {
			_bands[i - 1] = Raster<T>(filename, e, snap, i);
		}
	}
	template<typename T>
	inline MultiBandRaster<T>::MultiBandRaster(const Alignment& a, band_t nBands)
		: Alignment(a)
	{
		for (band_t band = 0; band < nBands; ++band) {
			_bands.emplace_back(a);
		}
	}
	template<typename T>
	inline band_t MultiBandRaster<T>::nBands() const {
		return (band_t)_bands.size();
	}
	template<typename T>
	inline Raster<T>& MultiBandRaster<T>::bandAt(band_t band) {
		_checkBand(band);
		return bandAtUnsafe(band);
	}
	template<typename T>
	inline Raster<T>& MultiBandRaster<T>::bandAtUnsafe(band_t band) {
		return _bands[band - 1];
	}
	template<typename T>
	inline const Raster<T>& MultiBandRaster<T>::bandAt(band_t band) const {
		_checkBand(band);
		return bandAtUnsafe(band);
	}
	template<typename T>
	inline const Raster<T>& MultiBandRaster<T>::bandAtUnsafe(band_t band) const {
		return _bands[band - 1];
	}
	template<typename T>
	inline auto MultiBandRaster<T>::atCell(cell_t cell, band_t band) {
		return bandAt(band).atCell(cell);
	}
	template<typename T>
	inline auto MultiBandRaster<T>::atCellUnsafe(cell_t cell, band_t band) {
		return bandAtUnsafe(band).atCellUnsafe(cell);
	}
	template<typename T>
	inline auto MultiBandRaster<T>::atRC(rowcol_t row, rowcol_t col, band_t band) {
		return bandAt(band).atRC(row, col);
	}
	template<typename T>
	inline auto MultiBandRaster<T>::atRCUnsafe(rowcol_t row, rowcol_t col, band_t band) {
		return bandAtUnsafe(band).atRCUnsafe(row, col);
	}
	template<typename T>
	inline auto MultiBandRaster<T>::atXY(coord_t x, coord_t y, band_t band) {
		return bandAt(band).atXY(x, y);
	}
	template<typename T>
	inline auto MultiBandRaster<T>::atXYUnsafe(coord_t x, coord_t y, band_t band) {
		return bandAtUnsafe(band).atXYUnsafe(x, y);
	}
	template<typename T>
	inline const auto MultiBandRaster<T>::atCell(cell_t cell, band_t band) const {
		return bandAt(band).atCell(cell);
	}
	template<typename T>
	inline const auto MultiBandRaster<T>::atCellUnsafe(cell_t cell, band_t band) const {
		return bandAtUnsafe(band).atCellUnsafe(cell);
	}
	template<typename T>
	inline const auto MultiBandRaster<T>::atRC(rowcol_t row, rowcol_t col, band_t band) const {
		return bandAt(band).atRC(row, col);
	}
	template<typename T>
	inline const auto MultiBandRaster<T>::atRCUnsafe(rowcol_t row, rowcol_t col, band_t band) const {
		return bandAtUnsafe(band).atRCUnsafe(row, col);
	}
	template<typename T>
	inline const auto MultiBandRaster<T>::atXY(coord_t x, coord_t y, band_t band) const {
		return bandAt(band).atXY(x, y);
	}
	template<typename T>
	inline const auto MultiBandRaster<T>::atXYUnsafe(coord_t x, coord_t y, band_t band) const {
		return bandAtUnsafe(band).atXYUnsafe(x, y);
	}
	template<typename T>
	inline void MultiBandRaster<T>::writeRaster(const std::string& fileName, const std::string driver, const T navalue, GDALDataType gdt)
	{
		if (gdt == GDT_Unknown) {
			gdt = _bands.front().GDT();
		}
		gdalAllRegisterThreadSafe();
		GDALDriver* d = GetGDALDriverManager()->GetDriverByName(driver.c_str());
		UniqueGdalDataset ugd = makeUniqueGdalDataset(d->Create(fileName.c_str(), ncol(), nrow(), nBands(), gdt, nullptr));
		if (!ugd) {
			throw InvalidRasterFileException("Unable to open " + fileName + " as a raster");
		}
		std::array<double, 6> gt = { _xmin, _xres,0,_ymax,0,-(_yres) };
		ugd->SetGeoTransform(gt.data());
		ugd->SetProjection(_crs.getCompleteWKT().c_str());
		for (band_t band = 0; band < nBands(); ++band) {
			for (cell_t cell = 0; cell < ncell(); ++cell) {
				auto v = atCellUnsafe(cell, band + 1);
				if (!v.has_value()) {
					v.value() = navalue;
				}
			}
		}
		for (band_t band = 0; band < nBands(); ++band) {
			auto thisBand = ugd->GetRasterBand(band+1);
			thisBand->SetNoDataValue((double)navalue);
			thisBand->RasterIO(GF_Write, 0, 0, _ncol, _nrow, _bands[band]._data.value().data(), _ncol, _nrow, gdt, 0, 0);
		}
	}
	template<typename T>
	inline void MultiBandRaster<T>::_checkBand(band_t band) const {
		if (band < 1 || band > _bands.size()) {
			throw std::out_of_range("Band out of range");
		}
	}

}

#endif