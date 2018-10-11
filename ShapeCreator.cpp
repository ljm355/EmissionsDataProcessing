#include "ShapeCreator.h"
#include <limits>
#include "qdir.h"
#include <fstream>
#include <sstream>
#include "Utils.h"
ShapeCreator::ShapeCreator()
{
}


ShapeCreator::~ShapeCreator()
{
}

std::vector<std::string> ShapeCreator::split(const char& delimiter, const std::string& line)
{
	std::string str = line;
	std::vector<std::string> splits;
	if (str.length() < 1)
		return splits;
	std::string::iterator iter = str.begin();
	int index = 0;
	int lastIndex = 0;
	while (iter != str.end() && index < str.length())
	{
		char& c = *iter;
		if (c == delimiter || index == str.length() - 1)
		{
			if (index == str.length() - 1)
				index = str.length();
			if (index - lastIndex>0)
				splits.push_back(str.substr(lastIndex, index - lastIndex));
			lastIndex = index + 1;
		}
		index++;
		iter++;
	}
	return splits;
}
void ShapeCreator::split(const char& delimiter, const std::string& line, std::string* psplits)
{
	std::string str = line;
	if (str.length() < 1)
		return;
	std::string::iterator iter = str.begin();
	int index = 0;
	int lastIndex = 0;
	while (iter != str.end() && index < str.length())
	{
		char& c = *iter;
		if (c == delimiter || index == str.length() - 1)
		{
			if (index == str.length() - 1)
				index = str.length();
			if (index - lastIndex > 0)
			{
				*psplits = str.substr(lastIndex, index - lastIndex);
				psplits++;
			}
			lastIndex = index + 1;
		}
		index++;
		iter++;
	}
}



//std::vector
union DataValue
{
	std::string s;
	double d;
	int i;
};


//struct ShapePointRecord
//{
//	/*2 IND
//	3 COM
//	6 NonRoad
//	7 Railroad
//	9 Airport*/
//	int SectorID;
//	double longitude;
//	double latitude;
//	int FIPS;
//	std::string facility_name;
//	std::string scc;
//	double New_LAT;
//	double New_LON;
//	double COMNG;
//	double COMpetrol;
//	double INDNG;
//	double INDpetrol;
//	double ELE;
//	double AIR;
//	double RAIL;
//	double Nonroad;
//
//	//facility_name scc New_LAT New_LON COM(tC / yr)	COM NG	COM petrol	IND(tC / yr)	IND NG	IND petrol	ELE(tC / yr)	ELE NG	ELE petrol	AIR(tC / yr)	RAIL(tC / yr)	Nonroad(tC / yr)
//	ShapePointRecord(std::string line)
//	{
//		//name	lat	lon	scc	facilityid	source	sector	ca11
//		scc = double2string(fea->GetFieldAsDouble("scc"));
//
//		OGRPoint*  p = (OGRPoint*)fea->GetGeometryRef();
//		lon = p->getX();
//		lat = p->getY();
//		ca = fea->GetFieldAsDouble("SR_CO2_tC");
//
//		facilityid = double2string(fea->GetFieldAsDouble("facility_i"));
//
//		sector = (int)fea->GetFieldAsDouble("Sector_ID");
//		FIPS = (int)fea->GetFieldAsDouble("FIPS");
//
//		name = fea->GetFieldAsString("facility_n");
//
//	}
//
//	void setFields(OGRFeature* poFeature)
//	{
//		poFeature->SetField("name", name.data());
//		poFeature->SetField("scc", scc.data());
//		poFeature->SetField("ca11", ca);
//		poFeature->SetField("facilityid", facilityid.data());
//
//		OGRPoint pt;
//		pt.setX(lon);
//		pt.setY(lat);
//		poFeature->SetGeometry(&pt);
//	}
//
//};
void ShapeCreator::create(std::string csvfielname, std::string outfilename, std::vector<ValueField> fields, ValueField lonfield, ValueField latfield, ValueField filterfield)
{
	QFileInfo fileInfo(csvfielname.data());
	
	ShapeFile shp;
	shp.create(outfilename.data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);


	std::ifstream ifs;
	ifs.open(csvfielname.data());
	std::string headerline;
	std::getline(ifs, headerline);
	std::vector<std::string> splits = Utils::splitCSV(',', headerline);
	for (size_t i = 0; i < fields.size(); i++)
	{
		fields[i].findIndex(splits);
		fields[i].createField(shp.poLayer);
	}
	filterfield.findIndex(splits);

	ValueField reffilterfield(filterfield.tp, filterfield.destname, filterfield.srcname, filterfield.srcname2);
	reffilterfield.findIndex(splits);
	latfield.tp = OGRFieldType::OFTReal;
	lonfield.tp = OGRFieldType::OFTReal;
	latfield.findIndex(splits);
	lonfield.findIndex(splits);
	int numrecords = 0;
	while (ifs.peek() != -1)
	{
		std::string line;
		std::getline(ifs, line);
		splits = Utils::splitCSV(',', line);
		latfield.set(splits);
		lonfield.set(splits);
		reffilterfield.set(splits);
		if (!reffilterfield.equal(&filterfield))
			continue;
		for (size_t i = 0; i < fields.size(); i++)
		{
			fields[i].set(splits);
		}

		OGRFeature* poFeature = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());

		OGRPoint pt;
		pt.setX(lonfield.d);
		pt.setY(latfield.d);
		poFeature->SetGeometry(&pt);
		for (size_t i = 0; i < fields.size(); i++)
		{
			fields[i].get(poFeature);
		}
		shp.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
		numrecords++;
	}
	if (numrecords == 0)
	{
		const char *pszDriverName = "ESRI Shapefile";
		OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName)->Delete(outfilename.data());
	}

}
void ShapeCreator::createBySector(std::string filename, std::string dir)
{
	QFileInfo fileInfo(filename.data());
	std::vector<ShapeFileBySector> shpBySectors;


	//std::string dir = "Z:/Hestia_Baltimore/gridPrep_SHP_master/";
	ShapeFile IndPoint;
	IndPoint.create((dir + "IndPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&IndPoint, "IND", 2));

	ShapeFile ComPoint;
	ComPoint.create((dir + "ComPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&ComPoint, "COM", 3));

	ShapeFile NonRoadPoint;
	NonRoadPoint.create((dir + "NonRoadPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&NonRoadPoint, "NonRoad", 6));

	ShapeFile createRailroadPoint;
	createRailroadPoint.create((dir + "RailroadPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&createRailroadPoint, "RAIL", 7));

	ShapeFile Airport;
	Airport.create((dir + "Airport.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&Airport, "AIR", 9));



	std::ifstream ifs;
	ifs.open(filename.data());
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);

	std::getline(ifs, line);
	std::vector<std::string> fieldtypeline = Utils::splitCSV(',', line);

	std::vector<bool> isColumnEmpty;
	std::vector<OGRFieldType> fieldType;
	int latidx = -1;
	int lonidx = -1;
	int sectIdx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		OGRFieldType tp = (OGRFieldType)atoi(fieldtypeline[i].data());
		fieldType.push_back(tp);
		OGRFieldDefn field(headerline[i].data(), tp);
		isColumnEmpty.push_back(false);

		for (size_t j = 0; j < shpBySectors.size(); j++)
		{
			ShapeFileBySector& shpBySector = shpBySectors[j];
			shpBySector.SHP->poLayer->CreateField(&field);
		}
		if (headerline[i] == "lat")
			latidx = i;
		if (headerline[i] == "lon")
			lonidx = i;
		if (headerline[i] == "sector")
			sectIdx = i;
	}




	std::vector<std::string> splits = headerline;
	int numfields = splits.size();
	std::string* psplits = &splits[0];
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		splits = Utils::splitCSV(',', line);
		idx++;
		for (size_t i = 0; i < shpBySectors.size(); i++)
		{
			ShapeFileBySector& shpBySector = shpBySectors[i];
			if (shpBySector.BC != psplits[sectIdx])
				continue;
			OGRFeature* poFeature = OGRFeature::CreateFeature(shpBySector.SHP->poLayer->GetLayerDefn());
			std::string latstr = psplits[latidx];
			std::string lonstr = psplits[lonidx];
			OGRPoint pt;
			pt.setX(atof(lonstr.data()));
			pt.setY(atof(latstr.data()));
			poFeature->SetGeometry(&pt);


			for (int n = 0; n < numfields; n++)
			{

				//typedef enum
				//{
				//	/** Simple 32bit integer */                   OFTInteger = 0,
				//	/** List of 32bit integers */                 OFTIntegerList = 1,
				//	/** Double Precision floating point */        OFTReal = 2,
				//	/** List of doubles */                        OFTRealList = 3,
				//	/** String of ASCII chars */                  OFTString = 4,
				//	/** Array of strings */                       OFTStringList = 5,
				//	/** deprecated */                             OFTWideString = 6,
				//	/** deprecated */                             OFTWideStringList = 7,
				//	/** Raw Binary data */                        OFTBinary = 8,
				//	/** Date */                                   OFTDate = 9,
				//	/** Time */                                   OFTTime = 10,
				//	/** Date and Time */                          OFTDateTime = 11,
				//	OFTMaxType = 11
				//} OGRFieldType;
				if (psplits[n] == "")
					continue;
				if (fieldType[n] == OFTReal)
				{
					poFeature->SetField(n, atof(psplits[n].data()));
				}
				else if (fieldType[n] == OFTInteger)
				{
					poFeature->SetField(n, atoi(psplits[n].data()));

				}
				else
				{
					poFeature->SetField(n, psplits[n].data());
				}
			}
			shpBySector.SHP->poLayer->CreateFeature(poFeature);
			OGRFeature::DestroyFeature(poFeature);
		}

	}

}
void ShapeCreator::create(std::string configFile, std::string outfile)
{

	ShapeFile ElecProd;
	ElecProd.create(outfile.data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	std::ifstream ifs;
	ifs.open(configFile.data());
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = split('\t', line);

	std::getline(ifs, line);
	std::vector<std::string> fieldtypeline = split('\t', line);

	std::vector<bool> isColumnEmpty;
	std::vector<OGRFieldType> fieldType;
	int latidx = -1;
	int lonidx = -1;
	int sectIdx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		OGRFieldType tp = (OGRFieldType)atoi(fieldtypeline[i].data());
		fieldType.push_back(tp);
		OGRFieldDefn field(headerline[i].data(), tp);
		isColumnEmpty.push_back(false);


		ElecProd.poLayer->CreateField(&field);

		if (headerline[i] == "lat")
			latidx = i;
		if (headerline[i] == "lon")
			lonidx = i;
		if (headerline[i] == "sector")
			sectIdx = i;
	}

	std::vector<std::string> splits = headerline;
	int numfields = splits.size();
	std::string* psplits = &splits[0];
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		split('\t', line, psplits);
		idx++;


		OGRFeature* poFeature = OGRFeature::CreateFeature(ElecProd.poLayer->GetLayerDefn());
		std::string latstr = psplits[latidx];
		std::string lonstr = psplits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		poFeature->SetGeometry(&pt);


		for (int n = 0; n < numfields; n++)
		{

			//typedef enum
			//{
			//	/** Simple 32bit integer */                   OFTInteger = 0,
			//	/** List of 32bit integers */                 OFTIntegerList = 1,
			//	/** Double Precision floating point */        OFTReal = 2,
			//	/** List of doubles */                        OFTRealList = 3,
			//	/** String of ASCII chars */                  OFTString = 4,
			//	/** Array of strings */                       OFTStringList = 5,
			//	/** deprecated */                             OFTWideString = 6,
			//	/** deprecated */                             OFTWideStringList = 7,
			//	/** Raw Binary data */                        OFTBinary = 8,
			//	/** Date */                                   OFTDate = 9,
			//	/** Time */                                   OFTTime = 10,
			//	/** Date and Time */                          OFTDateTime = 11,
			//	OFTMaxType = 11
			//} OGRFieldType;
			if (psplits[n] == "")
				continue;
			if (fieldType[n] == OFTReal)
			{
				poFeature->SetField(n, atof(psplits[n].data()));
			}
			else if (fieldType[n] == OFTInteger)
			{
				poFeature->SetField(n, atoi(psplits[n].data()));

			}
			else
			{
				poFeature->SetField(n, psplits[n].data());
			}
		}
		ElecProd.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);


	}

}