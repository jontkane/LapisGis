#include "LapisGis.hpp"

namespace lapis {

	static std::string executableFilePath() {
		int length = wai_getExecutablePath(nullptr, 0, nullptr);
		std::vector<char> path;
		path.resize((size_t)length + 1);
		wai_getExecutablePath(path.data(), length, nullptr);
		path[length] = '\0';
		std::string asString{ path.begin(),path.end() };
		std::filesystem::path asPath{ asString };
		return asPath.parent_path().string();
	}

    bool lapisGisInit(const std::string& projDbPath)
    {
        static bool initialized = false;
        static std::mutex initMutex;
        std::scoped_lock lock{ initMutex };
		if (initialized) {
			return true;
        }
        if (projDbPath.empty()) {
            setProjDefaultDirectory(executableFilePath());
		}
		else {
            setProjDefaultDirectory(projDbPath);
		}
		
#ifdef WIN32
        _putenv_s("PROJ_NETWORK", "OFF");
        _putenv_s("GDAL_DRIVER_PATH", "");
        _putenv_s("GDAL_SKIP", "");
#endif
		initialized = true;
		return true;
    }

#ifdef LAPISGIS_AUTO_INIT
	static bool lapisGisAutoInit = lapisGisInit("");
#endif
}
