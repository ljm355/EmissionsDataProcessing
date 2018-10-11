#pragma once
#include "ShapeFile.h"
#include <string>
#include <qstring.h>
#include <vector>
class ScaleByFuel
{
public:
	ScaleByFuel();
	~ScaleByFuel();
	void spatialAllocationTotal(ShapeFile * shp, std::vector<std::string> fuelfieldnames, std::string outputfieldname);
	void spatialAllocationForNonPoint(ShapeFile * shp, int baseyear, std::string weightField, std::string fuel, double total, std::string outputfieldname);
	void spatialAllocation(ShapeFile * shp, std::string weightField, double total, std::string outputfieldname);
	void spatialAllocationElecProd(ShapeFile * shp, std::string weightField, double total, std::string outputfieldname);
	void scale(std::string configfile, std::string dir);
	void scaleNonRoad(std::string rootdir, std::string cfgfile, std::string weightfield = "ca");
	void scaleNonRoad(std::vector<std::string> shapefiles, std::vector<std::string>  fields, std::vector<double> totals, std::string weightField);
	double sumField(std::string shapefile, std::string fieldname,std::vector<double>& values);
};

