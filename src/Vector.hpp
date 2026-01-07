#pragma once
#ifndef lapis_vector_h
#define lapis_vector_h

#include"gis_pch.hpp"
#include"Geometry.hpp"
#include"GisExceptions.hpp"

namespace lapis {
	class WrongFieldTypeException : public std::runtime_error {
	public:
		WrongFieldTypeException(const std::string& error);
	};

	enum class FieldType {
		Integer,
		Real,
		String
	};

	namespace {
		inline std::string getFieldName(OGRFieldDefn* field) {
			const char* altName = field->GetAlternativeNameRef();
			if (altName && altName[0] != '\0') {
				return std::string(altName);
			}
			return std::string(field->GetNameRef());
		}
	}

	class AttributeTable;
    template<class attribute_pointer>
    class AttributeRow;
	using MutableAttributeRow = AttributeRow<AttributeTable*>;
    using ConstAttributeRow = AttributeRow<const AttributeTable*>;
	template<class GEOM>
	class VectorDataset;

	class AttributeTable {
		template<class GEOM>
        friend class VectorDataset;
	public:
		using FeatureType = MutableAttributeRow;
        using ConstFeatureType = ConstAttributeRow;

		AttributeTable() = default;
        AttributeTable(const std::string& filename);
        AttributeTable(const std::filesystem::path& filename);

        void addStringField(const std::string& name, size_t width);
        void addIntegerField(const std::string& name);
        void addRealField(const std::string& name);
		template<class T>
		void addNumericField(const std::string& name);
		void removeField(const std::string& name);
        bool fieldExists(const std::string& name) const;

        void resize(size_t nrow);
		void addRow();
		size_t nFeature() const;

		const std::vector<std::string>& getAllFieldNames() const;
		FieldType getFieldType(const std::string& name) const;
		size_t getStringFieldWidth(const std::string& name) const;
        void setStringFieldWidth(const std::string& name, size_t width);

		std::string getStringField(size_t index, const std::string& name) const;
		int64_t getIntegerField(size_t index, const std::string& name) const;
        double getRealField(size_t index, const std::string& name) const;
		template<class T>
		T getNumericField(size_t index, const std::string& name) const;

		void setStringField(size_t index, const std::string& name, const std::string& value);
		void setIntegerField(size_t index, const std::string& name, int64_t value);
		void setRealField(size_t index, const std::string& name, double value);
		template<class T>
		void setNumericField(size_t index, const std::string& name, T value);

		template<class attribute_pointer>
		class base_iterator {
		public:
			base_iterator(attribute_pointer attributes, size_t index);
			base_iterator<attribute_pointer>& operator++();
			bool operator==(const base_iterator<attribute_pointer>& other) const = default;
            AttributeRow<attribute_pointer> operator*();
		private:
			attribute_pointer _attributes;
            size_t _index;
		};
        using iterator = base_iterator<AttributeTable*>;
        using const_iterator = base_iterator<const AttributeTable*>;

        ConstAttributeRow getRow(size_t index) const;
        MutableAttributeRow getRow(size_t index);
        iterator begin();
        iterator end();
        const_iterator begin() const;
        const_iterator end() const;
        ConstAttributeRow front() const;
        ConstAttributeRow back() const;
		MutableAttributeRow front();
		MutableAttributeRow back();

	private:
		void _reserve(size_t n);

		class FixedWidthString {
		public:
			std::string get() const;
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

    template<class attribute_pointer>
	class AttributeRow {
	public:

        constexpr static bool is_const = std::is_const_v<std::remove_pointer_t<attribute_pointer>>;

		AttributeRow(attribute_pointer attributes, size_t index);

		std::vector<std::string> getAllFieldNames() const;
		FieldType getFieldType(const std::string& name) const;
		size_t getStringFieldWidth(const std::string& name) const;
		bool fieldExists(const std::string& name) const;

		std::string getStringField(const std::string& name) const;
		int64_t getIntegerField(const std::string& name) const;
		double getRealField(const std::string& name) const;
		template<class T>
		T getNumericField(const std::string& name) const;
		
		void setStringField(const std::string& name, const std::string& value) requires !is_const {
            _attributes->setStringField(_index, name, value);
		}
		void setIntegerField(const std::string& name, int64_t value) requires !is_const{
            _attributes->setIntegerField(_index, name, value);
		}
		void setRealField(const std::string& name, double value) requires !is_const{
            _attributes->setRealField(_index, name, value);
		}
		template<class T>
		void setNumericField(const std::string& name, T value) requires !is_const{
            _attributes->setNumericField<T>(_index, name, value);
		}

        using const_attribute_pointer = std::add_const_t<std::remove_pointer_t<attribute_pointer>>*;
		operator AttributeRow<const_attribute_pointer>() const {
            return AttributeRow<const_attribute_pointer>(_attributes, _index);
		}

	private:
		attribute_pointer _attributes;
        size_t _index;
	};

	template<class GEOM, class attribute_pointer>
	class Feature;
	template<class GEOM>
    using MutableFeature = Feature<GEOM, AttributeTable*>;
    template<class GEOM>
    using ConstFeature = Feature<GEOM, const AttributeTable*>;

	template<class GEOM>
	class VectorDataset {
	public:
		using FeatureType = MutableFeature<GEOM>;
        using ConstFeatureType = ConstFeature<GEOM>;
		using Geometry = GEOM;

        VectorDataset() = default;
		explicit VectorDataset(const CoordRef& crs);
		VectorDataset(const std::string& filename);
		explicit VectorDataset(const std::filesystem::path& filename);
        VectorDataset(const std::vector<std::filesystem::path>& filenames);

		void writeShapefile(const std::filesystem::path& filename) const;
		void writeShapefile(const std::string& filename) const;

		void addStringField(const std::string& name, size_t width);
		void addIntegerField(const std::string& name);
		void addRealField(const std::string& name);
		template<class T>
		void addNumericField(const std::string& name);
		void removeField(const std::string& name);
		bool fieldExists(const std::string& name) const;

		const GEOM& getGeometry(size_t index) const;
		void replaceGeometry(size_t index, const GEOM& geom);
		void replaceGeometryPreciseExtent(size_t index, const GEOM& geom);
		void addGeometry(const GEOM& g);
		void addFeature(const ConstFeature<GEOM>& feature);
		void addFeatureUnsafe(const ConstFeature<GEOM>& feature);

		const CoordRef& crs() const;
		const Extent& extent() const;
		size_t nFeature() const;

		const std::vector<std::string>& getAllFieldNames() const;
		FieldType getFieldType(const std::string& name) const;
		size_t getStringFieldWidth(const std::string& name) const;
		void setStringFieldWidth(const std::string& name, size_t width);

		std::string getStringField(size_t index, const std::string& name) const;
		int64_t getIntegerField(size_t index, const std::string& name) const;
		double getRealField(size_t index, const std::string& name) const;
		template<class T>
		T getNumericField(size_t index, const std::string& name) const;

		void setStringField(size_t index, const std::string& name, const std::string& value);
		void setIntegerField(size_t index, const std::string& name, int64_t value);
		void setRealField(size_t index, const std::string& name, double value);
		template<class T>
		void setNumericField(size_t index, const std::string& name, T value);

        void projectInPlace(const CoordRef& newCrs);
        void projectInPlacePreciseExtent(const CoordRef& newCrs);

		void appendFile(const std::string& filename);
        void appendFile(const std::filesystem::path& filename);

		void reserve(size_t n);


		template<class attribute_pointer>
		class base_iterator {
		public:
			static constexpr bool is_const = std::is_const_v<std::remove_pointer_t<attribute_pointer>>;
			using geom_it = std::conditional_t<is_const, typename std::vector<GEOM>::const_iterator, typename std::vector<GEOM>::iterator>;

			base_iterator(AttributeTable::base_iterator<attribute_pointer> attributesIt, geom_it _geometryIt);
			base_iterator<attribute_pointer>& operator++();
            bool operator==(const base_iterator<attribute_pointer>& other) const = default;
			Feature<GEOM, attribute_pointer> operator*();
		private:
			AttributeTable::base_iterator<attribute_pointer> _attributesIt;
            geom_it _geometryIt;
		};
		using iterator = base_iterator<AttributeTable*>;
		using const_iterator = base_iterator<const AttributeTable*>;

		ConstFeature<GEOM> getFeature(size_t index) const;
		MutableFeature<GEOM> getFeature(size_t index);
		iterator begin();
		iterator end();
		const_iterator begin() const;
		const_iterator end() const;
		ConstFeature<GEOM> front() const;
		ConstFeature<GEOM> back() const;
		MutableFeature<GEOM> front();
		MutableFeature<GEOM> back();

	private:
		std::vector<GEOM> _geometries;
		Extent _extent;
		AttributeTable _attributes;

		void _projectInPlaceShared(const CoordRef& newCrs);
	};

	template<class GEOM>
    VectorDataset<GEOM> emptyVectorDatasetFromTemplate(const VectorDataset<GEOM>& templateDataset);

    template<class GEOM, class attribute_pointer>
	class Feature {
	public:

		static constexpr bool is_const = std::is_const_v<std::remove_pointer_t<attribute_pointer>>;
		using geom_pointer = std::conditional_t<is_const, const GEOM*, GEOM*>;

		Feature(AttributeRow<attribute_pointer> row, geom_pointer geometry);

		const GEOM& getGeometry() const;
		const CoordRef& crs() const;

		std::vector<std::string> getAllFieldNames() const;
		FieldType getFieldType(const std::string& name) const;
		size_t getStringFieldWidth(const std::string& name) const;
		bool fieldExists(const std::string& name) const;

		std::string getStringField(const std::string& name) const;
		int64_t getIntegerField(const std::string& name) const;
		double getRealField(const std::string& name) const;
		template<class T>
		T getNumericField(const std::string& name) const;

		void setStringField(const std::string& name, const std::string& value) requires !is_const{
            _attributeRow.setStringField(name, value);
		}
		void setIntegerField(const std::string& name, int64_t value) requires !is_const{
            _attributeRow.setIntegerField(name, value);
		}
		void setRealField(const std::string& name, double value) requires !is_const{
            _attributeRow.setRealField(name, value);
		}
        template<class T>
		void setNumericField(const std::string& name, T value) requires !is_const {
            _attributeRow.setNumericField<T>(name, value);
		}

        using const_attribute_pointer = std::add_const_t<std::remove_pointer_t<attribute_pointer>>*;
		operator Feature<GEOM, const_attribute_pointer>() const {
            return Feature<GEOM, const_attribute_pointer>(static_cast<ConstAttributeRow>(_attributeRow), _geometry);
		}

		operator AttributeRow<attribute_pointer>() const requires !is_const {
            return _attributeRow;
		}

		operator AttributeRow<const_attribute_pointer>() const {
            return static_cast<AttributeRow<const_attribute_pointer>>(_attributeRow);
        }


	private:
        AttributeRow<attribute_pointer> _attributeRow;
        geom_pointer _geometry;
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
		auto type = getFieldType(name);
		switch (type) {
		case FieldType::Integer:
			return static_cast<T>(getIntegerField(index, name));
		case FieldType::Real:
			return static_cast<T>(getRealField(index, name));
		default:
			throw WrongFieldTypeException("Cannot cast field " + name + " to numeric type");
		}
	}
	template<class T>
	inline void AttributeTable::setNumericField(size_t index, const std::string& name, T value)
	{
		auto type = getFieldType(name);
		switch (type) {
		case FieldType::Integer:
			setIntegerField(index, name, static_cast<int64_t>(value));
			break;
		case FieldType::Real:
			setRealField(index, name, static_cast<double>(value));
			break;
		default:
			throw WrongFieldTypeException("Cannot cast field " + name + " to numeric type");
        }
	}
	template<class attribute_pointer>
	inline AttributeTable::base_iterator<attribute_pointer>::base_iterator(attribute_pointer attributes, size_t index)
        : _attributes(attributes), _index(index)
	{
	}
    template<class attribute_pointer>
	inline typename AttributeTable::base_iterator<attribute_pointer>& AttributeTable::base_iterator<attribute_pointer>::operator++()
	{
		++_index;
		return *this;
	}
    template<class attribute_pointer>
	inline AttributeRow<attribute_pointer> AttributeTable::base_iterator<attribute_pointer>::operator*()
	{
		return AttributeRow<attribute_pointer>(_attributes, _index);
	}

    template<class attribute_pointer>
	inline AttributeRow<attribute_pointer>::AttributeRow(attribute_pointer attributes, size_t index)
		: _attributes(attributes), _index(index)
	{
	}

    template<class attribute_pointer>
	inline std::vector<std::string> AttributeRow<attribute_pointer>::getAllFieldNames() const
	{
        return _attributes->getAllFieldNames();
	}
    template<class attribute_pointer>
	inline FieldType AttributeRow<attribute_pointer>::getFieldType(const std::string& name) const
	{
		return _attributes->getFieldType(name);
	}
    template<class attribute_pointer>
	inline size_t AttributeRow<attribute_pointer>::getStringFieldWidth(const std::string& name) const
	{
		return _attributes->getStringFieldWidth(name);
	}
	template<class attribute_pointer>
	inline bool AttributeRow<attribute_pointer>::fieldExists(const std::string& name) const {
        return _attributes->fieldExists(name);
	}

    template<class attribute_pointer>
	inline std::string AttributeRow<attribute_pointer>::getStringField(const std::string& name) const
	{
		return _attributes->getStringField(_index, name);
	}
    template<class attribute_pointer>
	inline int64_t AttributeRow<attribute_pointer>::getIntegerField(const std::string& name) const
	{
		return _attributes->getIntegerField(_index, name);
	}
    template<class attribute_pointer>
	inline double AttributeRow<attribute_pointer>::getRealField(const std::string& name) const
	{
		return _attributes->getRealField(_index, name);
	}
    template<class attribute_pointer>
    template<class T>
	inline T AttributeRow<attribute_pointer>::getNumericField(const std::string& name) const
	{
		return _attributes->getNumericField<T>(_index, name);
	}

	template<class GEOM>
	inline VectorDataset<GEOM>::VectorDataset(const CoordRef& crs)
	{
		_extent.defineCRS(crs);
	}
    template<class GEOM>
	inline VectorDataset<GEOM>::VectorDataset(const std::string& filename)
	{
		gdalAllRegisterThreadSafe();
		UniqueGdalDataset shp = vectorGDALWrapper(filename);
		if (!shp) {
            throw InvalidVectorFileException("Could not open " + filename);
		}
		OGRLayer* layer = shp->GetLayer(0);
		if (!layer) {
			throw InvalidVectorFileException("Could not read layer from " + filename);
        }
		if (wkbFlatten(layer->GetGeomType()) != GEOM::gdalGeometryTypeStatic) {
			if (wkbFlatten(layer->GetGeomType()) == wkbPolygon && GEOM::gdalGeometryTypeStatic == wkbMultiPolygon) {
				//esri shp files do not have a MultiPolygon type, so this mismatch is common
			}
			else {
				throw WrongGeometryTypeException(filename + " is not the expected geometry type");
			}
		}
		bool initFields = false;
		reserve(layer->GetFeatureCount());
		OGRSpatialReference* osr = layer->GetSpatialRef();
		CoordRef crs(osr);
		OGREnvelope envelope;
		layer->GetExtent(&envelope);
		_extent = Extent(envelope.MinX, envelope.MaxX, envelope.MinY, envelope.MaxY, crs);

		for (const OGRFeatureUniquePtr& feature : layer) {
			if (!initFields) {
				for (int i = 0; i < feature->GetFieldCount(); ++i) {
					OGRFieldDefn* field = feature->GetFieldDefnRef(i);
					switch (field->GetType()) {
					case OFTInteger:
					case OFTInteger64:
						addIntegerField(getFieldName(field));
						break;
					case OFTReal:
						addRealField(getFieldName(field));
						break;
					case OFTString:
						addStringField(getFieldName(field), field->GetWidth());
						break;
					default:
						throw std::runtime_error("unimplemented field type when reading shapefile");
					}
				}
				initFields = true;
			}
			OGRGeometry* gdalGeometry = feature->GetGeometryRef();
			GEOM lapisGeometry{ *gdalGeometry, crs };
			addGeometry(lapisGeometry);
			for (int i = 0; i < feature->GetFieldCount(); ++i) {
				OGRFieldDefn* field = feature->GetFieldDefnRef(i);
				switch (field->GetType()) {
				case OFTInteger:
				case OFTInteger64:
					setIntegerField(_geometries.size() - 1, getFieldName(field), feature->GetFieldAsInteger64(field->GetNameRef()));
					break;
				case OFTReal:
					setRealField(_geometries.size() - 1, getFieldName(field), feature->GetFieldAsDouble(field->GetNameRef()));
					break;
				case OFTString:
					setStringField(_geometries.size() - 1, getFieldName(field), feature->GetFieldAsString(field->GetNameRef()));
					break;
				default:
					throw WrongFieldTypeException("Unimplemented field type when reading shapefile");
				}
			}
		}
	}
    template<class GEOMETRY>
	inline VectorDataset<GEOMETRY>::VectorDataset(const std::filesystem::path& filename)
		: VectorDataset<GEOMETRY>(filename.string())
	{
    }

	template<class GEOM>
	inline VectorDataset<GEOM>::VectorDataset(const std::vector<std::filesystem::path>& filenames)
	{
		namespace fs = std::filesystem;
		//first do a pass through the files to ensure they have the correct geometry type
		//and all have matching CRSes and attributes
		//keep a running count of features so we can reserve the sum
        size_t totalFeatures = 0;
        CoordRef commonCrs;
		std::vector<std::string> fieldNames;
		std::vector<FieldType> fieldTypes;
        std::vector<size_t> fieldWidths;
        bool firstFile = true;
		for (const fs::path& filename : filenames) {
			gdalAllRegisterThreadSafe();
			UniqueGdalDataset shp = vectorGDALWrapper(filename.string());
			if (!shp) {
				throw InvalidVectorFileException("Could not open " + filename.string());
			}
			OGRLayer* layer = shp->GetLayer(0);
			if (!layer) {
				throw InvalidVectorFileException("Could not read layer from " + filename.string());
			}
			if (wkbFlatten(layer->GetGeomType()) != GEOM::gdalGeometryTypeStatic) {
				if (wkbFlatten(layer->GetGeomType()) == wkbPolygon && GEOM::gdalGeometryTypeStatic == wkbMultiPolygon) {
					//esri shp files do not have a MultiPolygon type, so this mismatch is common
				}
				else {
					throw WrongGeometryTypeException(filename.string() + " is not the expected geometry type");
				}
			}

			if (firstFile) {
                commonCrs = CoordRef(layer->GetSpatialRef());
				for (int i = 0; i < layer->GetLayerDefn()->GetFieldCount(); ++i) {
					OGRFieldDefn* field = layer->GetLayerDefn()->GetFieldDefn(i);
					fieldNames.push_back(getFieldName(field));
					switch (field->GetType()) {
					case OFTInteger:
					case OFTInteger64:
						fieldTypes.push_back(FieldType::Integer);
                        fieldWidths.push_back(0);
						break;
					case OFTReal:
						fieldTypes.push_back(FieldType::Real);
                        fieldWidths.push_back(0);
						break;
					case OFTString:
						fieldTypes.push_back(FieldType::String);
                        fieldWidths.push_back(static_cast<size_t>(field->GetWidth()));
						break;
					default:
						throw std::runtime_error("unimplemented field type when reading shapefile");
                    }
				}
				totalFeatures += layer->GetFeatureCount();
				firstFile = false;
				continue;
			}

            CoordRef fileCrs(layer->GetSpatialRef());
			if (!commonCrs.isConsistent(fileCrs)) {
				throw CRSMismatchException("CRS of " + filename.string() + " does not match other files in the set");
            }
			if (layer->GetLayerDefn()->GetFieldCount() != static_cast<int>(fieldNames.size())) {
				throw std::runtime_error("Attribute table of " + filename.string() + " does not match other files in the set");
			}
			for (int i = 0; i < layer->GetLayerDefn()->GetFieldCount(); ++i) {
				OGRFieldDefn* field = layer->GetLayerDefn()->GetFieldDefn(i);
				if (getFieldName(field) != fieldNames[i]) {
					throw std::runtime_error("Attribute table of " + filename.string() + " does not match other files in the set");
				}
				FieldType expectedType = fieldTypes[i];
				FieldType fieldType;
				switch (field->GetType()) {
				case OFTInteger:
				case OFTInteger64:
					fieldType = FieldType::Integer;
					break;
				case OFTReal:
					fieldType = FieldType::Real;
					break;
				case OFTString:
					fieldType = FieldType::String;
					break;
				default:
					throw std::runtime_error("unimplemented field type when reading shapefile");
				}
				if (fieldType != expectedType) {
					throw std::runtime_error("Attribute table of " + filename.string() + " does not match other files in the set");
				}
				if (fieldType == FieldType::String) {
					size_t width = static_cast<size_t>(field->GetWidth());
                    fieldWidths[i] = std::max(fieldWidths[i], width);
                }
            }

            totalFeatures += layer->GetFeatureCount();
		}

        reserve(totalFeatures);
        _extent.defineCRS(commonCrs);
		for (size_t i = 0; i < fieldNames.size(); ++i) {
			switch (fieldTypes[i]) {
			case FieldType::Integer:
				addIntegerField(fieldNames[i]);
				break;
			case FieldType::Real:
				addRealField(fieldNames[i]);
				break;
			case FieldType::String:
				addStringField(fieldNames[i], fieldWidths[i]);
				break;
			}
        }
		for (const fs::path& filename : filenames) {
			appendFile(filename);
        }
	}

    template<class GEOM>
	inline void VectorDataset<GEOM>::writeShapefile(const std::string& filename) const
	{
		gdalAllRegisterThreadSafe();
		UniqueGdalDataset outshp = gdalCreateWrapperVector(filename.c_str());
		if (!outshp) {
			throw InvalidVectorFileException("Could not create shapefile " + filename);
		}
		OGRSpatialReference crs;
		crs.importFromWkt(_extent.crs().getCleanEPSG().getCompleteWKT().c_str());

		OGRLayer* layer = outshp->CreateLayer("layer", &crs, GEOM::gdalGeometryTypeStatic, nullptr);
		if (layer == nullptr) {
            std::string gdalError = CPLGetLastErrorMsg();
            throw InvalidVectorFileException("Error when writing " + filename + ": " + gdalError);
		}

		for (const auto& fieldName : getAllFieldNames()) {
			OGRFieldDefn newField = OGRFieldDefn(fieldName.c_str(), OFTString);
			newField.SetAlternativeName(fieldName.c_str());
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

        const std::vector<std::string>& fieldNames = getAllFieldNames();

		for (auto&& feature : *this) {
			UniqueOGRFeature gdalFeature = createFeatureWrapper(layer);
			for (size_t i = 0; i < fieldNames.size(); ++i) {
                const std::string& fieldName = fieldNames[i];
				switch (getFieldType(fieldName)) {
				case FieldType::String:
					gdalFeature->SetField(static_cast<int>(i), feature.getStringField(fieldName).c_str());
					break;
				case FieldType::Real:
					gdalFeature->SetField(static_cast<int>(i), feature.getRealField(fieldName));
					break;
				case FieldType::Integer:
					gdalFeature->SetField(static_cast<int>(i), feature.getIntegerField(fieldName));
					break;
				}
			}
			std::unique_ptr<OGRGeometry> geometry = feature.getGeometry().gdalGeometryGeneric();
			gdalFeature->SetGeometry(geometry.get());
			layer->CreateFeature(gdalFeature.get());
		}
	}
    template<class GEOM>
	inline void VectorDataset<GEOM>::writeShapefile(const std::filesystem::path& filename) const
	{
		writeShapefile(filename.string());
	}

    template<class GEOM>
	inline void VectorDataset<GEOM>::addStringField(const std::string& name, size_t width)
	{
		_attributes.addStringField(name, width);
	}
    template<class GEOM>
	inline void VectorDataset<GEOM>::addIntegerField(const std::string& name)
	{
		_attributes.addIntegerField(name);
	}
    template<class GEOM>
	inline void VectorDataset<GEOM>::addRealField(const std::string& name)
	{
        _attributes.addRealField(name);
	}
    template<class GEOM>
	template<class T>
	inline void VectorDataset<GEOM>::addNumericField(const std::string& name)
	{
		_attributes.addNumericField<T>(name);
	}
    template<class GEOM>
	inline void VectorDataset<GEOM>::removeField(const std::string& name)
	{
		_attributes.removeField(name);
	}

    template<class GEOM>
	inline const GEOM& VectorDataset<GEOM>::getGeometry(size_t index) const
	{
		return _geometries[index];
	}
	template<class GEOM>
	inline void VectorDataset<GEOM>::replaceGeometry(size_t index, const GEOM& geom)
	{
		_geometries[index] = geom;
        _extent = extendExtent(_extent, geom.boundingBox());
	}
	template<class GEOM>
	inline void VectorDataset<GEOM>::replaceGeometryPreciseExtent(size_t index, const GEOM& geom)
	{
		_geometries[index] = geom;
		bool init = false;
		for (const Feature<GEOM, AttributeTable*>& f : *this) {
			if (!init) {
				_extent = f.getGeometry().boundingBox();
				init = true;
			}
			else {
				_extent = extendExtent(_extent, f.getGeometry().boundingBox());
			}
        }
	}
    template<class GEOM>
	inline void VectorDataset<GEOM>::addGeometry(const GEOM& g)
	{
		_geometries.push_back(g);
		_attributes.addRow();
		if (nFeature() == 1) {
			_extent = g.boundingBox();
		}
		else {
			_extent = extendExtent(_extent, g.boundingBox());
        }
	}

	template<class GEOM>
	inline void VectorDataset<GEOM>::addFeature(const ConstFeature<GEOM>& feature)
	{
		if (!feature.getGeometry().crs().isConsistent(_extent.crs())) {
			throw CRSMismatchException("Feature CRS does not match dataset CRS");
        }
		for (const auto& fieldName : feature.getAllFieldNames()) {
			if (!fieldExists(fieldName)) {
				throw std::runtime_error("Feature attribute table does not match dataset attribute table");
			}
			if (getFieldType(fieldName) != feature.getFieldType(fieldName)) {
				throw WrongFieldTypeException("Feature attribute table does not match dataset attribute table");
			}
        }
		addFeatureUnsafe(feature);
	}

	template<class GEOM>
	inline void VectorDataset<GEOM>::addFeatureUnsafe(const ConstFeature<GEOM>& feature)
	{
		//making a deep copy first to avoid invalidating the feature if it's from this vectordataset

        GEOM geometryCopy = feature.getGeometry();
        std::unordered_map<std::string, std::variant<int64_t, double, std::string>> attributeValues;
		for (const auto& fieldName : getAllFieldNames()) {
			switch (getFieldType(fieldName)) {
			case FieldType::String:
				attributeValues[fieldName] = feature.getStringField(fieldName);
				break;
			case FieldType::Real:
				attributeValues[fieldName] = feature.getRealField(fieldName);
				break;
			case FieldType::Integer:
				attributeValues[fieldName] = feature.getIntegerField(fieldName);
				break;
            }
		}
        addGeometry(geometryCopy);
		for (const auto& fieldName : getAllFieldNames()) {
			switch (getFieldType(fieldName)) {
			case FieldType::String:
                back().setStringField(fieldName, std::get<std::string>(attributeValues[fieldName]));
                break;
            case FieldType::Real:
				back().setRealField(fieldName,  std::get<double>(attributeValues[fieldName]));
                break;
            case FieldType::Integer:
				back().setIntegerField(fieldName,  std::get<int64_t>(attributeValues[fieldName]));
                break;
			}
		}


		if (nFeature() == 1) {
			_extent = geometryCopy.boundingBox();
		}
		else {
			_extent = extendExtent(_extent, geometryCopy.boundingBox());
        }
	}

    template<class GEOM>
	inline const CoordRef& VectorDataset<GEOM>::crs() const
	{
		return _extent.crs();
	}
    template<class GEOM>
	inline const Extent& VectorDataset<GEOM>::extent() const
	{
		return _extent;
	}
    template<class GEOM>
	inline size_t VectorDataset<GEOM>::nFeature() const
	{
		return _geometries.size();
	}

    template<class GEOM>
	inline const std::vector<std::string>& VectorDataset<GEOM>::getAllFieldNames() const
	{
		return _attributes.getAllFieldNames();
	}
    template<class GEOM>
	inline FieldType VectorDataset<GEOM>::getFieldType(const std::string& name) const
	{
		return _attributes.getFieldType(name);
	}
    template<class GEOM>
	inline size_t VectorDataset<GEOM>::getStringFieldWidth(const std::string& name) const
	{
		return _attributes.getStringFieldWidth(name);
	}
    template<class GEOM>
	inline void VectorDataset<GEOM>::setStringFieldWidth(const std::string& name, size_t width)
	{
		_attributes.setStringFieldWidth(name, width);
	}
	template<class GEOM>
	inline bool VectorDataset<GEOM>::fieldExists(const std::string& name) const {
        return _attributes.fieldExists(name);
	}

    template<class GEOM>
	inline std::string VectorDataset<GEOM>::getStringField(size_t index, const std::string& name) const
	{
		return _attributes.getStringField(index, name);
	}
    template<class GEOM>
	inline int64_t VectorDataset<GEOM>::getIntegerField(size_t index, const std::string& name) const
	{
		return _attributes.getIntegerField(index, name);
	}
    template<class GEOM>
	inline double VectorDataset<GEOM>::getRealField(size_t index, const std::string& name) const
	{
		return _attributes.getRealField(index, name);
	}
    template<class GEOM>
    template<class T>
	inline T VectorDataset<GEOM>::getNumericField(size_t index, const std::string& name) const
	{
		return _attributes.getNumericField<T>(index, name);
	}

    template<class GEOM>
	inline void VectorDataset<GEOM>::setStringField(size_t index, const std::string& name, const std::string& value)
	{
		_attributes.setStringField(index, name, value);
	}
	template<class GEOM>
	inline void VectorDataset<GEOM>::setIntegerField(size_t index, const std::string& name, int64_t value)
	{
		_attributes.setIntegerField(index, name, value);
    }
    template<class GEOM>
	inline void VectorDataset<GEOM>::setRealField(size_t index, const std::string& name, double value)
	{
        _attributes.setRealField(index, name, value);
	}
	template<class GEOM>
    template<class T>
	inline void VectorDataset<GEOM>::setNumericField(size_t index, const std::string& name, T value)
	{
        _attributes.setNumericField<T>(index, name, value);
	}

	template<class GEOM>
	inline void VectorDataset<GEOM>::projectInPlace(const CoordRef& newCrs)
	{
		if (_extent.crs().isConsistent(newCrs)) {
			return;
		}
		_projectInPlaceShared(newCrs);
		
		_extent = QuadExtent(_extent, newCrs).outerExtent();
    }
	template<class GEOM>
	inline void VectorDataset<GEOM>::projectInPlacePreciseExtent(const CoordRef& newCrs)
	{
		if (_extent.crs().isConsistent(newCrs)) {
			return;
		}
		_projectInPlaceShared(newCrs);

		bool init = false;
		Extent newExtent;
		for (const GEOM& geometry : _geometries) {
			if (!init) {
				newExtent = geometry.boundingBox();
				init = true;
			}
			else {
				newExtent = extendExtent(newExtent, geometry.boundingBox());
			}
		}
		_extent = newExtent;
	}

	template<class GEOM>
	inline void VectorDataset<GEOM>::appendFile(const std::string& filename)
	{
		gdalAllRegisterThreadSafe();
		UniqueGdalDataset otherShp = vectorGDALWrapper(filename);
		if (!otherShp) {
            throw InvalidVectorFileException("Could not open vector file " + filename + " to append");
		}
		OGRLayer* layer = otherShp->GetLayer(0);
		if (!layer) {
			throw InvalidVectorFileException("Could not open vector file " + filename + " to append");
        }
		if (wkbFlatten(layer->GetGeomType()) != GEOM::gdalGeometryTypeStatic) {
			if (wkbFlatten(layer->GetGeomType()) == wkbPolygon && GEOM::gdalGeometryTypeStatic == wkbMultiPolygon) {
				//esri shp files do not have a MultiPolygon type, so this mismatch is common
			}
			else {
				throw WrongGeometryTypeException(filename + " is not the expected geometry type");
			}
		}
		if (layer->GetFeatureCount()==0) {
			return;
        }

		OGRSpatialReference* osr = layer->GetSpatialRef();
		CoordRef otherCrs(osr);
		if (!otherCrs.isConsistentHoriz(_extent.crs())) {
            throw CRSMismatchException("Cannot append vector file " + filename + " with different CRS");
        }

		OGRFeature* exampleFeature = layer->GetFeature(0);
		for (const std::string& fieldName : getAllFieldNames()) {
            int idx = exampleFeature->GetFieldIndex(fieldName.c_str());
			if (idx == -1) {
                throw std::runtime_error("Field " + fieldName + " not found in " + filename + " when appending vector dataset");
			}
			OGRFieldDefn* field = exampleFeature->GetFieldDefnRef(idx);
			switch (field->GetType()) {
			case OFTInteger:
			case OFTInteger64:
				if (getFieldType(fieldName) != FieldType::Integer) {
					throw WrongFieldTypeException("Field " + fieldName + " has different type in " + filename + " when appending vector dataset");
                }
				break;
			case OFTReal:
				if (getFieldType(fieldName) != FieldType::Real) {
					throw WrongFieldTypeException("Field " + fieldName + " has different type in " + filename + " when appending vector dataset");
                }
				break;
			case OFTString:
				if (getFieldType(fieldName) != FieldType::String) {
					throw WrongFieldTypeException("Field " + fieldName + " has different type in " + filename + " when appending vector dataset");
                }
                setStringFieldWidth(fieldName, std::max(getStringFieldWidth(fieldName), static_cast<size_t>(field->GetWidth())));
				break;
			default:
				throw std::runtime_error("unimplemented field type when reading shapefile");
			}
		}

        reserve(nFeature() + layer->GetFeatureCount());
		for (const OGRFeatureUniquePtr& feature : layer) {
			OGRGeometry* gdalGeometry = feature->GetGeometryRef();
			GEOM lapisGeometry{ *gdalGeometry, otherCrs };
			addGeometry(lapisGeometry);
			for (const std::string& fieldName : getAllFieldNames()) {
				int idx = feature->GetFieldIndex(fieldName.c_str());
				OGRFieldDefn* field = feature->GetFieldDefnRef(idx);
				switch (getFieldType(fieldName)) {
                case FieldType::Integer:
					back().setIntegerField(fieldName, feature->GetFieldAsInteger64(fieldName.c_str()));
					break;
                case FieldType::Real:
					back().setRealField(fieldName, feature->GetFieldAsDouble(fieldName.c_str()));
					break;
                case FieldType::String:
					back().setStringField(fieldName, feature->GetFieldAsString(fieldName.c_str()));
					break;
				default:
					throw WrongFieldTypeException("Unimplemented field type when reading shapefile");
				}
            }
		}
	}

	template<class GEOM>
	inline void VectorDataset<GEOM>::appendFile(const std::filesystem::path& filename)
	{
        appendFile(filename.string());
	}

	template<class GEOM>
	inline void VectorDataset<GEOM>::_projectInPlaceShared(const CoordRef& newCrs) {
		if (_extent.crs().isConsistent(newCrs)) {
			return;
        }
		const CoordTransform& transform = CoordTransformFactory::getTransform(_extent.crs(), newCrs);

		for (GEOM& geometry : _geometries) {
			geometry.projectInPlace(transform);
        }
	}

    template<class GEOM>
    template<class attribute_pointer>
    inline VectorDataset<GEOM>::base_iterator<attribute_pointer>::base_iterator(
		AttributeTable::base_iterator<attribute_pointer> attributesIt, typename base_iterator<attribute_pointer>::geom_it geometryIt)
        : _attributesIt(attributesIt), _geometryIt(geometryIt)
	{ }
    template<class GEOM>
    template<class attribute_pointer>
	inline typename VectorDataset<GEOM>::base_iterator<attribute_pointer>& VectorDataset<GEOM>::base_iterator<attribute_pointer>::operator++()
	{
		++_geometryIt;
		++_attributesIt;
		return *this;
	}
	template<class GEOM>
    template<class attribute_pointer>
	inline Feature<GEOM, attribute_pointer> VectorDataset<GEOM>::base_iterator<attribute_pointer>::operator*()
	{
		return Feature<GEOM, attribute_pointer>(*_attributesIt, &*_geometryIt);
	}

	template<class GEOM>
	ConstFeature<GEOM> VectorDataset<GEOM>::getFeature(size_t index) const
	{
		return ConstFeature<GEOM>(_attributes.getRow(index), &_geometries[index]);
    }
    template<class GEOM>
	MutableFeature<GEOM> VectorDataset<GEOM>::getFeature(size_t index)
	{
		return MutableFeature<GEOM>(_attributes.getRow(index), &_geometries[index]);
	}
    template<class GEOM>
	inline VectorDataset<GEOM>::iterator VectorDataset<GEOM>::begin()
	{
		return iterator(_attributes.begin(), _geometries.begin());
	}
	template<class GEOM>
	inline VectorDataset<GEOM>::iterator VectorDataset<GEOM>::end()
	{
		return iterator(_attributes.end(), _geometries.end());
    }
    template<class GEOM>
	inline VectorDataset<GEOM>::const_iterator VectorDataset<GEOM>::begin() const
	{
        return const_iterator(_attributes.begin(), _geometries.begin());
	}
    template<class GEOM>
	inline VectorDataset<GEOM>::const_iterator VectorDataset<GEOM>::end() const
	{
		return const_iterator(_attributes.end(), _geometries.end());
	}
    template<class GEOM>
	inline ConstFeature<GEOM> VectorDataset<GEOM>::front() const
	{
		return getFeature(0);
	}
    template<class GEOM>
	inline ConstFeature<GEOM> VectorDataset<GEOM>::back() const
	{
		return getFeature(nFeature() - 1);
	}
    template<class GEOM>
	inline MutableFeature<GEOM> VectorDataset<GEOM>::front()
	{
		return getFeature(0);
	}
    template<class GEOM>
	inline MutableFeature<GEOM> VectorDataset<GEOM>::back()
	{
		return getFeature(nFeature() - 1);
	}

	template<class GEOM>
	inline void VectorDataset<GEOM>::reserve(size_t n) {
		_attributes._reserve(n);
        _geometries.reserve(n);
	}

    template<class GEOM, class attribute_pointer>
	inline Feature<GEOM, attribute_pointer>::Feature(AttributeRow<attribute_pointer> row, geom_pointer geometry)
		: _attributeRow(row), _geometry(geometry)
	{
	}
    template<class GEOM, class attribute_pointer>
	inline const GEOM& Feature<GEOM, attribute_pointer>::getGeometry() const
	{
        return *_geometry;
	}
    template<class GEOM, class attribute_pointer>
	inline const CoordRef& Feature<GEOM, attribute_pointer>::crs() const
	{
        return _geometry->crs();
	}

    template<class GEOM, class attribute_pointer>
	inline std::vector<std::string> Feature<GEOM, attribute_pointer>::getAllFieldNames() const
	{
		return _attributeRow.getAllFieldNames();
	}
    template<class GEOM, class attribute_pointer>
	inline FieldType Feature<GEOM, attribute_pointer>::getFieldType(const std::string& name) const
	{
		return _attributeRow.getFieldType(name);
	}
    template<class GEOM, class attribute_pointer>
	inline size_t Feature<GEOM, attribute_pointer>::getStringFieldWidth(const std::string& name) const
	{
		return _attributeRow.getStringFieldWidth(name);
	}
	template<class GEOM, class attribute_pointer>
	inline bool Feature<GEOM, attribute_pointer>::fieldExists(const std::string& name) const {
        return _attributeRow.fieldExists(name);
	}

    template<class GEOM, class attribute_pointer>
	inline std::string Feature<GEOM, attribute_pointer>::getStringField(const std::string& name) const
	{
		return _attributeRow.getStringField(name);
	}
    template<class GEOM, class attribute_pointer>
	inline int64_t Feature<GEOM, attribute_pointer>::getIntegerField(const std::string& name) const
	{
        return _attributeRow.getIntegerField(name);
	}
    template<class GEOM, class attribute_pointer>
	inline double Feature<GEOM, attribute_pointer>::getRealField(const std::string& name) const
	{
		return _attributeRow.getRealField(name);
	}
	template<class GEOM, class attribute_pointer>
    template<class T>
	inline T Feature<GEOM, attribute_pointer>::getNumericField(const std::string& name) const
	{
        return _attributeRow.getNumericField<T>(name);
	}
	template<class GEOM>
	VectorDataset<GEOM> emptyVectorDatasetFromTemplate(const VectorDataset<GEOM>& templateDataset)
	{
        VectorDataset<GEOM> newDataset(templateDataset.crs());
		for (const auto& fieldName : templateDataset.getAllFieldNames()) {
			switch (templateDataset.getFieldType(fieldName)) {
			case FieldType::String:
				newDataset.addStringField(fieldName, templateDataset.getStringFieldWidth(fieldName));
				break;
			case FieldType::Real:
				newDataset.addRealField(fieldName);
				break;
			case FieldType::Integer:
				newDataset.addIntegerField(fieldName);
				break;
            }
		}
		return newDataset;
	}
}

#endif