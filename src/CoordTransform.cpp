#include"gis_pch.hpp"
#include"CoordTransform.hpp"


namespace lapis{
	CoordTransform::CoordTransform(const CoordRef& src, const CoordRef& dst) : _src(src), _dst(dst) {
		_conv = LinearUnitConverter(src.getZUnits(), dst.getZUnits());
		_needZConv = !src.isConsistentZUnits(dst);
		_needXYConv = !src.isConsistentHoriz(dst);
		if (_needXYConv) {
			_tr = projCrsToCrsWrapper(src.getSharedPtr(), dst.getSharedPtr());
		}
	}

	PJ* CoordTransform::getPtr() {
		return _tr.get();
	}

	const PJ * CoordTransform::getPtr() const {
		return _tr.get();
	}

	const SharedPJ& CoordTransform::getSharedPtr() const
	{
		return _tr;
	}

	CoordXY CoordTransform::transformSingleXY(coord_t x, coord_t y) const
	{
		proj_trans_generic(_tr.get(), PJ_FWD,
			&x, 0, 1,
			&y, 0, 1,
			nullptr, 0, 0,
			nullptr, 0, 0);
		return { x,y };
	}
	const CoordRef& CoordTransform::src() const
	{
		return _src;
	}
	const CoordRef& CoordTransform::dst() const
	{
		return _dst;
	}
	size_t CoordTransformFactory::CoordRefPairHasher::operator()(const CoordRefPair& p) const
	{
		size_t h1 = CoordRefHasher()(p.first);
		size_t h2 = CoordRefHasher()(p.second);
        return h1 ^ (h2 << 1);
	}
	bool CoordTransformFactory::CoordRefPairEqual::operator()(const CoordRefPair& a, const CoordRefPair& b) const
	{
        return a.first.equalForHash(b.first) && a.second.equalForHash(b.second);
	}
	const CoordTransform& CoordTransformFactory::getTransform(const CoordRef& src, const CoordRef& dst)
	{
		CoordRefPair p = std::make_pair(src, dst);
		{
			std::shared_lock lock{ _mut };
			auto it = _cache.find(p);
			if (it != _cache.end()) {
				return *(it->second);
			}
		}

		{
			std::unique_lock lock{ _mut };
			auto it = _cache.find(p);
			if (it == _cache.end()) {
				_cache[p] = std::make_unique<CoordTransform>(src, dst);
				it = _cache.find(p);
			}
			return *(it->second);
		}
	}
	std::shared_mutex CoordTransformFactory::_mut = std::shared_mutex{};
    std::unordered_map<CoordTransformFactory::CoordRefPair, std::unique_ptr<CoordTransform>, CoordTransformFactory::CoordRefPairHasher, CoordTransformFactory::CoordRefPairEqual> CoordTransformFactory::_cache = std::unordered_map<CoordTransformFactory::CoordRefPair, std::unique_ptr<CoordTransform>, CoordTransformFactory::CoordRefPairHasher, CoordTransformFactory::CoordRefPairEqual>{};
}
