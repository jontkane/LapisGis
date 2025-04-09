#include"MultiBandRaster.hpp"

namespace lapis {
    band_t nBandsForFile(const std::string& file)
    {
		UniqueGdalDataset wgd = rasterGDALWrapper(file);
		if (!wgd) {
			return 0;
		}
		return wgd->GetRasterCount();
    }
}