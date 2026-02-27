#pragma once
#ifndef lapis_projwrappers_h
#define lapis_projwrappers_h

#include"gis_pch.hpp"
#include"LapisGisTypeDefs.hpp"


namespace lapis {

	std::string executableFilePath();

	using SharedPJ = std::shared_ptr<PJ>;
	SharedPJ makeSharedPJ(PJ* pj);

	using SharedPJCtx = std::shared_ptr<PJ_CONTEXT>;
	SharedPJCtx getNewPJContext();

	class ProjContextByThread {
	private:
		inline static std::unordered_map<std::thread::id, SharedPJCtx> _ctxs;
	public:
		static PJ_CONTEXT* get();
	};

	SharedPJ projCreateWrapper(const std::string& s);
	SharedPJ projCrsToCrsWrapper(SharedPJ from, SharedPJ to);
	SharedPJ getSubCrs(const SharedPJ base, int index);

	SharedPJ sharedPJFromOSR(const OGRSpatialReference& osr);

    SharedPJ getHorizontalCrs(const SharedPJ& crs);

	class PJIdentifyWrapper {
	public:
		PJIdentifyWrapper(const SharedPJ& p, const std::string& auth = "EPSG");
		~PJIdentifyWrapper();

		const int* getConfidence() const;
		const PJ_OBJ_LIST* getObjList() const;
	private:
		int* _confidence;
		PJ_OBJ_LIST* _obj;
	};

	bool setProjDefaultDirectory(const std::string& path);
	bool setProjDirectory(const std::string& path, PJ_CONTEXT* context);

	bool isProjDirectorySet();
	std::string getProjDirectory();
	bool projDbExists();
}

#endif