#pragma once
#include"gis_pch.hpp"
#include"CoordRef.hpp"
#include"Raster.hpp"

namespace lapis {

	class WrongGeometryTypeException : public std::runtime_error {
	public:
		WrongGeometryTypeException(const std::string& error);
	};
	class WrongFieldTypeException : public std::runtime_error {
	public:
		WrongFieldTypeException(const std::string& error);
	};

	class Geometry {
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbUnknown;
		using GdalEquivalent = OGRGeometry;

		virtual OGRwkbGeometryType gdalGeometryType() const;

		virtual const OGRGeometry& gdalGeometryGeneric() const = 0;

		const CoordRef& crs() const;
		void setCrs(const CoordRef& crs);

	protected:
		CoordRef _crs;
		Geometry() = default;
		Geometry(const Geometry&) = default;
	};
	class Point : public Geometry {
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbPoint;
		using GdalEquivalent = OGRPoint;

		Point() = delete;
		Point(const OGRGeometry& geom);
		Point(coord_t x, coord_t y);
		Point(coord_t x, coord_t y, const CoordRef& crs);
		Point(CoordXY xy);
		Point(CoordXY xy, const CoordRef& crs);

		OGRwkbGeometryType gdalGeometryType() const override;
		const OGRPoint& gdalGeometry() const;
		const OGRGeometry& gdalGeometryGeneric() const override;

		coord_t x() const;
		coord_t y() const;
	private:
		OGRPoint _point;
	};
	class Polygon : public Geometry {
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbPolygon;
		using GdalEquivalent = OGRPolygon;

		Polygon() = delete;
		Polygon(const OGRGeometry& geom);
		Polygon(const std::vector<CoordXY>& outerRing);
		Polygon(const std::vector<CoordXY>& outerRing, const CoordRef& crs);
		Polygon(const Extent& e);
		Polygon(const QuadExtent& q);

		OGRwkbGeometryType gdalGeometryType() const override;
		const OGRPolygon& gdalGeometry() const;
		const OGRGeometry& gdalGeometryGeneric() const override;
		
		void addInnerRing(const std::vector<CoordXY>& innerRing);
	private:
		OGRPolygon _polygon;
	};
	class MultiPolygon : public Geometry {
	public:
		constexpr static OGRwkbGeometryType gdalGeometryTypeStatic = wkbMultiPolygon;
		using GdalEquivalent = OGRMultiPolygon;

		MultiPolygon() = delete;
		MultiPolygon(const OGRGeometry& geom);

		OGRwkbGeometryType gdalGeometryType() const override;
		const OGRMultiPolygon& gdalGeometry() const;
		const OGRGeometry& gdalGeometryGeneric() const override;

		void addPolygon(const Polygon& polygon);
	private:
		OGRMultiPolygon _multiPolygon;
	};

	enum class FieldType {
		Integer,
		Real,
		String
	};

	class AttributeTable {
	public:
		void addStringField(const std::string& name, size_t width);
		void addIntegerField(const std::string& name);
		void addRealField(const std::string& name);
		template<class T>
		void addNumericField(const std::string& name);
		
		void resize(size_t nrow);
		void addRow();

		size_t nFeature() const;

		std::vector<std::string> getAllFieldNames() const;
		FieldType getFieldType(const std::string& name) const;
		size_t getStringFieldWidth(const std::string& name) const;

		const std::string& getStringField(size_t index, const std::string& name) const;
		int64_t getIntegerField(size_t index, const std::string& name) const;
		double getRealField(size_t index, const std::string& name) const;
		template<class T>
		T getNumericField(size_t index, const std::string& name) const;

		void setStringField(size_t index, const std::string& name, const std::string& value);
		void setIntegerField(size_t index, const std::string& name, int64_t value);
		void setRealField(size_t index, const std::string& name, double value);
		template<class T>
		void setNumericField(size_t index, const std::string& name, T value);
	private:

		class FixedWidthString {
		public:
			const std::string& get() const;
			void set(const std::string& value, size_t width);
		private:
			std::string _data;
		};

		using Variant = std::variant<int64_t, double, FixedWidthString>;
		struct Field {
			FieldType type;
			size_t width; //only used if the type is String
			std::vector<Variant> values;
		};

		size_t _nrow = 0;
		std::unordered_map<std::string, Field> _fields;
		std::vector<std::string> _fieldNamesInOrder;
	};

	template<class GEOMETRY>
	class VectorsAndAttributes : public AttributeTable {

	private:
		class SingleGeometryWithAttributes;
		class iterator;

	public:
		VectorsAndAttributes() = default;
		explicit VectorsAndAttributes(const CoordRef& crs);
		VectorsAndAttributes(const std::string& filename);
		VectorsAndAttributes(const std::filesystem::path& filename);

		void writeShapefile(const std::filesystem::path& filename);

		iterator begin();
		iterator end();

		const CoordRef& crs() const;
		void setCrs(const CoordRef& crs);

		const GEOMETRY& getGeometry(size_t index) const;
		void replaceGeometry(size_t index, const GEOMETRY& geom);
		void addGeometry(const GEOMETRY& g);

		SingleGeometryWithAttributes getFeature(size_t index);
		SingleGeometryWithAttributes front();
		SingleGeometryWithAttributes back();
	private:
		CoordRef _crs;
		std::vector<GEOMETRY> _geometries;

		class SingleGeometryWithAttributes {
		public:
			SingleGeometryWithAttributes(AttributeTable& fullAttributeTable, const GEOMETRY& geometry, size_t attributeIndex);

			const GEOMETRY& getGeometry() const;

			const CoordRef& crs() const;

			std::vector<std::string> getAllFieldNames() const;
			FieldType getFieldType(const std::string& name) const;
			size_t getStringFieldWidth(const std::string& name) const;

			const std::string& getStringField(const std::string& name) const;
			int64_t getIntegerField(const std::string& name) const;
			double getRealField(const std::string& name) const;
			template<class T>
			T getNumericField(const std::string& name) const;

			void setStringField(const std::string& name, const std::string& value);
			void setIntegerField(const std::string& name, int64_t value);
			void setRealField(const std::string& name, double value);
			template<class T>
			void setNumericField(const std::string& name, T value);
		private:
			const GEOMETRY& _geometry;
			AttributeTable& _attributes;
			size_t _attributeIndex;
		};

		class iterator {
		public:
			iterator(AttributeTable* attributes, std::vector<GEOMETRY>::iterator geomIt, size_t attributeIndex);

			iterator& operator++();
			bool operator==(const iterator& other) const = default;
			SingleGeometryWithAttributes operator*();
		private:
			AttributeTable* _attributes;
			size_t _attributeIndex;
			std::vector<GEOMETRY>::iterator _geomIt;
		};

		void _constructFromFilename(const std::string& filename);
	};

	template<class T>
	inline void AttributeTable::addNumericField(const std::string& name)
	{
		if constexpr (std::is_integral<T>()) {
			addIntegerField(name);
		}
		else if constexpr (std::is_floating_point<T>()) {
			addRealField(name);
		}
		else {
			[] <bool flag = false>()
			{
				static_assert(flag, "incorrect type in addTemplateField");
			}();
		}
	}
	template<class T>
	inline T AttributeTable::getNumericField(size_t index, const std::string& name) const
	{
		if constexpr (std::is_integral<T>()) {
			return (T)getIntegerField(index, name);
		}
		else if constexpr (std::is_floating_point<T>()) {
			return (T)getRealField(index, name);
		}
		else {
			[] <bool flag = false>()
			{
				static_assert(flag, "incorrect type in getNumericField");
			}();
			return 0;
		}
	}
	template<class T>
	inline void AttributeTable::setNumericField(size_t index, const std::string& name, T value)
	{
		if constexpr (std::is_integral<T>()) {
			setIntegerField(index, name, value);
		}
		else if constexpr (std::is_floating_point<T>()) {
			setRealField(index, name, value);
		}
		else {
			[] <bool flag = false>()
			{
				static_assert(flag, "incorrect type in setNumericField");
			}();
		}
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::VectorsAndAttributes(const CoordRef& crs)
	{
		_crs = crs;
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::VectorsAndAttributes(const std::string& filename)
	{
		_constructFromFilename(filename);
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::VectorsAndAttributes(const std::filesystem::path& filename)
	{
		_constructFromFilename(filename.string());
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::writeShapefile(const std::filesystem::path& filename)
	{
		gdalAllRegisterThreadSafe();
		UniqueGdalDataset outshp = gdalCreateWrapper("ESRI Shapefile", filename.string().c_str(), 0, 0, GDT_Unknown);
		OGRSpatialReference crs;
		crs.importFromWkt(_crs.getCleanEPSG().getCompleteWKT().c_str());

		OGRLayer* layer = outshp->CreateLayer("layer", &crs, GEOMETRY::gdalGeometryTypeStatic, nullptr);

		for (const auto& fieldName : getAllFieldNames()) {
			OGRFieldDefn newField = OGRFieldDefn(fieldName.c_str(), OFTString);
			switch (getFieldType(fieldName)) {
			case FieldType::String:
				newField.SetType(OFTString);
				layer->CreateField(&newField);
				break;
			case FieldType::Real:
				newField.SetType(OFTReal);
				layer->CreateField(&newField);
				break;
			case FieldType::Integer:
				newField.SetType(OFTInteger64);
				layer->CreateField(&newField);
				break;
			}
		}
		for (SingleGeometryWithAttributes feature : *this) {
			UniqueOGRFeature gdalFeature = createFeatureWrapper(layer);
			for (const auto& fieldName : getAllFieldNames()) {
				switch (getFieldType(fieldName)) {
				case FieldType::String:
					gdalFeature->SetField(fieldName.c_str(), feature.getStringField(fieldName).c_str());
					break;
				case FieldType::Real:
					gdalFeature->SetField(fieldName.c_str(), feature.getRealField(fieldName));
					break;
				case FieldType::Integer:
					gdalFeature->SetField(fieldName.c_str(), feature.getIntegerField(fieldName));
					break;
				}
			}
			const GEOMETRY::GdalEquivalent& geometry = feature.getGeometry().gdalGeometry();
			gdalFeature->SetGeometry(&geometry);
			layer->CreateFeature(gdalFeature.get());
		}
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::iterator VectorsAndAttributes<GEOMETRY>::begin()
	{
		return iterator(this, _geometries.begin(), 0);
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::iterator VectorsAndAttributes<GEOMETRY>::end()
	{
		return iterator(this, _geometries.end(), nFeature() - 1);
	}
	template<class GEOMETRY>
	inline const CoordRef& VectorsAndAttributes<GEOMETRY>::crs() const
	{
		return _crs;
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::setCrs(const CoordRef& crs)
	{
		_crs = crs;
	}
	template<class GEOMETRY>
	inline const GEOMETRY& VectorsAndAttributes<GEOMETRY>::getGeometry(size_t index) const
	{
		_geometries[index];
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::replaceGeometry(size_t index, const GEOMETRY& geom)
	{
		_geometries[index] = geom;
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::addGeometry(const GEOMETRY& g)
	{
		_geometries.push_back(g);
		addRow();
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes VectorsAndAttributes<GEOMETRY>::getFeature(size_t index)
	{
		return SingleGeometryWithAttributes(*this, _geometries[index], index);
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes VectorsAndAttributes<GEOMETRY>::front()
	{
		return getFeature(0);
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes VectorsAndAttributes<GEOMETRY>::back()
	{
		return getFeature(nFeature() - 1);
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::_constructFromFilename(const std::string& filename)
	{
		gdalAllRegisterThreadSafe();
		UniqueGdalDataset shp = vectorGDALWrapper(filename);
		if (!shp) {
			return;
		}
		OGRLayer* layer = shp->GetLayer(0);
		if (wkbFlatten(layer->GetGeomType()) != GEOMETRY::gdalGeometryTypeStatic) {
			throw InvalidVectorFileException(filename + " is not the expected geometry type");
		}
		bool initFields = false;
		for (const OGRFeatureUniquePtr& feature : layer) {
			if (!initFields) {
				for (int i = 0; i < feature->GetFieldCount(); ++i) {
					OGRFieldDefn* field = feature->GetFieldDefnRef(i);
					switch (field->GetType()) {
					case OFTInteger:
					case OFTInteger64:
						addIntegerField(field->GetNameRef());
						break;
					case OFTReal:
						addRealField(field->GetNameRef());
						break;
					case OFTString:
						addStringField(field->GetNameRef(), field->GetWidth());
						break;
					default:
						throw std::runtime_error("unimplemented field type when reading shapefile");
					}
				}
				initFields = true;
			}
			OGRGeometry* gdalGeometry = feature->GetGeometryRef();
			GEOMETRY lapisGeometry{ *gdalGeometry };
			addGeometry(lapisGeometry);
			for (int i = 0; i < feature->GetFieldCount(); ++i) {
				OGRFieldDefn* field = feature->GetFieldDefnRef(i);
				switch (field->GetType()) {
				case OFTInteger:
				case OFTInteger64:
					setIntegerField(_geometries.size() - 1, field->GetNameRef(), feature->GetFieldAsInteger64(field->GetNameRef()));
					break;
				case OFTReal:
					setRealField(_geometries.size() - 1, field->GetNameRef(), feature->GetFieldAsDouble(field->GetNameRef()));
					break;
				case OFTString:
					setStringField(_geometries.size() - 1, field->GetNameRef(), feature->GetFieldAsString(field->GetNameRef()));
					break;
				default:
					throw WrongFieldTypeException("Unimplemented field type when reading shapefile");
				}
			}
		}
		OGRSpatialReference* osr = layer->GetSpatialRef();
		_crs = CoordRef(osr);
	}

	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::SingleGeometryWithAttributes
	(AttributeTable& fullAttributeTable, const GEOMETRY& geometry, size_t attributeIndex)
		: _attributes(fullAttributeTable), _geometry(geometry), _attributeIndex(attributeIndex)
	{}
	template<class GEOMETRY>
	inline const GEOMETRY& VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getGeometry() const
	{
		return _geometry;
	}
	template<class GEOMETRY>
	inline const CoordRef& VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::crs() const
	{
		return _geometry.crs();
	}
	template<class GEOMETRY>
	inline std::vector<std::string> VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getAllFieldNames() const
	{
		return _attributes.getAllFieldNames();
	}
	template<class GEOMETRY>
	inline FieldType VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getFieldType(const std::string& name) const
	{
		return _attributes.getFieldType(name);
	}
	template<class GEOMETRY>
	inline size_t VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getStringFieldWidth(const std::string& name) const
	{
		return _attributes.getStringFieldWidth(name);
	}
	template<class GEOMETRY>
	inline const std::string& VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getStringField(const std::string& name) const
	{
		return _attributes.getStringField(_attributeIndex, name);
	}
	template<class GEOMETRY>
	inline int64_t VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getIntegerField(const std::string& name) const
	{
		return _attributes.getIntegerField(_attributeIndex, name);
	}
	template<class GEOMETRY>
	inline double VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getRealField(const std::string& name) const
	{
		return _attributes.getRealField(_attributeIndex, name);
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::setStringField(const std::string& name, const std::string& value)
	{
		_attributes.setStringField(_attributeIndex, name, value);
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::setIntegerField(const std::string& name, int64_t value)
	{
		_attributes.setIntegerField(_attributeIndex, name, value);
	}
	template<class GEOMETRY>
	inline void VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::setRealField(const std::string& name, double value)
	{
		_attributes.setRealField(_attributeIndex, name, value);
	}
	template<class GEOMETRY>
	template<class T>
	inline T VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::getNumericField(const std::string& name) const
	{
		return _attributes.getNumericField<T>(_attributeIndex, name);
	}
	template<class GEOMETRY>
	template<class T>
	inline void VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes::setNumericField(const std::string& name, T value)
	{
		_attributes.setNumericField<T>(_attributeIndex, name, value);
	}

	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::iterator::iterator(AttributeTable* attributes, std::vector<GEOMETRY>::iterator geomIt, size_t attributeIndex)
		: _attributes(attributes), _geomIt(geomIt), _attributeIndex(attributeIndex)
	{}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::iterator& VectorsAndAttributes<GEOMETRY>::iterator::operator++()
	{
		_geomIt++;
		_attributeIndex++;
		return *this;
	}
	template<class GEOMETRY>
	inline VectorsAndAttributes<GEOMETRY>::SingleGeometryWithAttributes VectorsAndAttributes<GEOMETRY>::iterator::operator*()
	{
		return SingleGeometryWithAttributes(*_attributes, *_geomIt, _attributeIndex);
	}
}