#include "Geographies.h"

double Geographies::calFraction(OGRFeature* fea, const int& footprintIdx)
{

	OGRGeometry* geo = fea->GetGeometryRef();
	OGRwkbGeometryType gtype = geo->getGeometryType();
	if (gtype == wkbLineString || gtype == wkbLineString25D || gtype == wkbMultiLineString)
	{
		double len1 = Utils::calPolylineLength(fea->GetGeometryRef());// *FOOTPRINT_SCALE_FACTOR;
		double len2 = fea->GetFieldAsDouble(footprintIdx);
		return len1 / len2;
	}
	else if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		double area1 = Utils::calPolygonArea(fea->GetGeometryRef());
		double area2 = fea->GetFieldAsDouble(footprintIdx);
		return area1 / area2;
	}

	return 1;
}

void Geographies::setup(std::string srcShpFile, std::vector<OGRFieldData> fields, std::string destShpFile)
{
	m_AttributeData.clear();
	m_destFileName = destShpFile;
	std::vector<int> srcIndices;
	std::vector<int> destIndices;
	ShapeFile srcShapes(srcShpFile);
	ShapeFile destShapes(destShpFile, 1);
	int numfeas = destShapes.poLayer->GetFeatureCount();
	for (int i = 0; i < srcShapes.poLayer->GetLayerDefn()->GetFieldCount(); i++)
	{
		OGRFieldDefn* field = srcShapes.poLayer->GetLayerDefn()->GetFieldDefn(i);
		std::string fieldname = field->GetNameRef();
		for (int j = 0; j < fields.size(); j++)
		{
			if (fieldname == fields[j].srcName)
			{
				OGRFieldData fdata = fields[j];
				fdata.tp = field->GetType();
				fdata.srcIdx = i;
				OGRFieldDefn destField(field);
				destField.SetName(fdata.destName.data());
				fdata.destIdx = destShapes.getOrCreateField(&destField);
				fdata.destName = destShapes.poLayer->GetLayerDefn()->GetFieldDefn(fdata.destIdx)->GetNameRef();
				m_AttributeData.push_back(fdata);
				m_AttributeData[m_AttributeData.size() - 1].data.resize(numfeas);
			}
		}

	}

}
void Geographies::gather(std::string srcShpFile, std::string fidName, OperatorType tp)
{

	    ShapeFile input(srcShpFile);
		int fidIndex = input.getField(fidName.data());
		OGRFeature *poFeature;
		input.poLayer->ResetReading();

		OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
		int footprintIndex = -1;
		if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
		{
			footprintIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("length");
		}
		if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
		{
			footprintIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("area");
		}
		int idx = 0;
		for (int ifield = 0; ifield < m_AttributeData.size(); ifield++)
		{
			m_AttributeData[ifield].begin();
		}
		while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
		{

			double fraction = calFraction(poFeature, footprintIndex);
			int fid = poFeature->GetFieldAsInteger(fidIndex);
			for (int ifield = 0; ifield < m_AttributeData.size(); ifield++)
			{
				OGRFieldData& fdata = m_AttributeData[ifield];
				double val = poFeature->GetFieldAsDouble(fdata.srcIdx);
				if (fdata.operatorType == OperatorType::SUM)
				{
					val = val * fraction;
				}
				fdata.set(val, fid);
			}
			OGRFeature::DestroyFeature(poFeature);
			idx++;
		}
		for (int ifield = 0; ifield < m_AttributeData.size(); ifield++)
		{
			m_AttributeData[ifield].end();
		}
	
}

void Geographies::update()
{
	ShapeFile input(m_destFileName,1);
	OGRFeature *poFeature;
	input.poLayer->ResetReading();
	int idx = 0;
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{

		for (int ifield = 0; ifield < m_AttributeData.size(); ifield++)
		{
			m_AttributeData[ifield].setField(poFeature,idx);
		}
		input.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
		idx++;
	}
	input.close();
}
