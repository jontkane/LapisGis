#include"gis_pch.hpp"
#include"projwrappers.hpp"
#include"GDALWrappers.hpp"

namespace lapis {

	std::string executableFilePath() {
		int length = wai_getExecutablePath(nullptr, 0, nullptr);
		std::vector<char> path;
		path.resize((size_t)length + 1);
		wai_getExecutablePath(path.data(), length, nullptr);
		path[length] = '\0';
		std::string asString{ path.begin(),path.end() };
		std::filesystem::path asPath{ asString };
		return asPath.parent_path().string();
	}

	PJ_CONTEXT* ProjContextByThread::get()
	{
		static std::mutex mut;
		std::scoped_lock<std::mutex> lock{ mut };
		std::thread::id thisthread = std::this_thread::get_id();
		if (!_ctxs.count(thisthread)) {
			_ctxs.emplace(thisthread, getNewPJContext());

#ifdef LAPISGIS_PROJDB_IN_EXE_DIR
			setProjDirectory(executableFilePath(), _ctxs.at(thisthread).get());
#endif
		}
		return _ctxs.at(thisthread).get();
	}
	SharedPJ makeSharedPJ(PJ* pj)
	{
		return SharedPJ(pj,
			[](PJ* pj) {
				if (pj) {
					proj_destroy(pj);
				}
			}
		);
	}
	SharedPJCtx getNewPJContext()
	{
		return SharedPJCtx(proj_context_create(),
			[](PJ_CONTEXT* pjc) {
				if (pjc) {
					proj_context_destroy(pjc);
				}
			}
		);
	}
	SharedPJ projCreateWrapper(const std::string& s)
	{
		return makeSharedPJ(proj_create(ProjContextByThread::get(), s.c_str()));
	}
	SharedPJ projCrsToCrsWrapper(SharedPJ from, SharedPJ to)
	{
		//static const char* ballpark = "ALLOW_BALLPARK=YES";
        //static const char* options[2] = { ballpark, nullptr };
		SharedPJ out = makeSharedPJ(proj_create_crs_to_crs_from_pj(
			ProjContextByThread::get(), from.get(), to.get(), nullptr, nullptr)
		);
		if (out) {
			out = makeSharedPJ(
				proj_normalize_for_visualization(ProjContextByThread::get(), out.get())
			);
		}
		return out;
	}
	SharedPJ getSubCrs(const SharedPJ base, int index)
	{
		return makeSharedPJ(
			proj_crs_get_sub_crs(ProjContextByThread::get(), base.get(), index)
		);
	}
	SharedPJ sharedPJFromOSR(const OGRSpatialReference& osr)
	{
		UniqueGdalString wgs = exportToWktWrapper(osr);
		return projCreateWrapper(wgs.get());
	}
	SharedPJ getHorizontalCrs(const SharedPJ& crs)
	{
		if (!crs) {
            return SharedPJ();
		}
        PJ_TYPE t = proj_get_type(crs.get());
		SharedPJ temp;
		switch (t) {
		case PJ_TYPE_UNKNOWN:
		case PJ_TYPE_ELLIPSOID:
		case PJ_TYPE_PRIME_MERIDIAN:
		case PJ_TYPE_GEODETIC_REFERENCE_FRAME:
		case PJ_TYPE_DYNAMIC_GEODETIC_REFERENCE_FRAME:
		case PJ_TYPE_VERTICAL_REFERENCE_FRAME:
        case PJ_TYPE_DYNAMIC_VERTICAL_REFERENCE_FRAME:
		case PJ_TYPE_DATUM_ENSEMBLE:
		case PJ_TYPE_VERTICAL_CRS:
		case PJ_TYPE_CONVERSION:
		case PJ_TYPE_TRANSFORMATION:
		case PJ_TYPE_CONCATENATED_OPERATION:
		case PJ_TYPE_OTHER_COORDINATE_OPERATION:
		case PJ_TYPE_TEMPORAL_DATUM:
		case PJ_TYPE_ENGINEERING_DATUM:
		case PJ_TYPE_PARAMETRIC_DATUM:
		case PJ_TYPE_COORDINATE_METADATA:
			return SharedPJ();
		case PJ_TYPE_CRS:
		case PJ_TYPE_GEOCENTRIC_CRS:
		case PJ_TYPE_GEOGRAPHIC_CRS:
		case PJ_TYPE_GEOGRAPHIC_2D_CRS:
		case PJ_TYPE_PROJECTED_CRS:
		case PJ_TYPE_TEMPORAL_CRS:
		case PJ_TYPE_ENGINEERING_CRS:
		case PJ_TYPE_OTHER_CRS:
			return crs;
		case PJ_TYPE_GEOGRAPHIC_3D_CRS:
            return getHorizontalCrs(makeSharedPJ(proj_crs_demote_to_2D(ProjContextByThread::get(), nullptr, crs.get())));
		case PJ_TYPE_COMPOUND_CRS:
			temp = getHorizontalCrs(getSubCrs(crs, 0));
			if (!temp) {
                return getHorizontalCrs(getSubCrs(crs, 1));
			}
			return temp;
		case PJ_TYPE_BOUND_CRS:
		case PJ_TYPE_DERIVED_PROJECTED_CRS:
            return getHorizontalCrs(makeSharedPJ(proj_get_source_crs(ProjContextByThread::get(), crs.get())));
		default:
			return SharedPJ();
		}
	}
	bool setProjDirectory(const std::string& path, PJ_CONTEXT* context)
	{
		namespace fs = std::filesystem;
		std::string folder;
		if (!fs::is_directory(path)) {
			folder = fs::path(path).parent_path().string();
		}
		else {
			folder = path;
		}
		char* data = folder.data();
		proj_context_set_search_paths(context, 1, &data);

		return true;
	}
	PJIdentifyWrapper::PJIdentifyWrapper(const SharedPJ& p, const std::string& auth) : _obj(nullptr), _confidence(nullptr)
	{
		_obj = proj_identify(ProjContextByThread::get(), p.get(), auth.c_str(), nullptr, &_confidence);
	}
	PJIdentifyWrapper::~PJIdentifyWrapper()
	{
		if (_confidence) {
			proj_int_list_destroy(_confidence);
		}
		if (_obj) {
			proj_list_destroy(_obj);
		}
	}
	const int* PJIdentifyWrapper::getConfidence() const
	{
		return _confidence;
	}
	const PJ_OBJ_LIST* PJIdentifyWrapper::getObjList() const
	{
		return _obj;
	}

#ifdef LAPISGIS_PROJDB_IN_EXE_DIR
	bool ProjContextByThread::set_proj_db_for_null_context = setProjDirectory(executableFilePath(), nullptr);
	bool ProjContextByThread::set_proj_lib = _putenv(("PROJ_LIB="+executableFilePath()).c_str());
	bool ProjContextByThread::set_proj_data = _putenv(("PROJ_DATA=" + executableFilePath()).c_str());
#else
	bool ProjContextByThread::set_proj_db_for_null_context = false;
	bool ProjContextByThread::set_proj_lib = false;
    bool ProjContextByThread::set_proj_data = false;
#endif
}