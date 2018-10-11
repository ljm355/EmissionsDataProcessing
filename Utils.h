#pragma once
#include "ShapeFile.h"
#include <vector>
class Utils
{
public:
	Utils();
	~Utils();
public:
	static double calPolygonArea(OGRGeometry* geom);
	static double calPolylineLength(OGRGeometry* geom);
	static void updateArea(ShapeFile* file, double scale = 1.0);
	static void updateLengh(ShapeFile* file, double scale = 1.0);
	static void updateFootprint(std::string filename,double scale = 1.0);
	static double getDistanceFromLatLonInMeter(double lat1, double lon1, double lat2, double lon2);
	static double getDegreeToMeter(OGRLayer* layer);
	static void updateFootprintForDir(std::string indir, bool force = false, double scale = 1.0);
	static void dropFieldsForDir(std::string indir, std::vector<std::string> fields2drop);
	static void renameFieldsForDir(std::string indir, std::vector<std::string> srcFieldNames, std::vector<std::string> destFieldNames);
	static void keepFieldsForDir(std::string indir, std::vector<std::string> fields2keep);
	static void copyFolder(std::string drcdir, std::string destdir, std::vector<std::string> fields2keep);
	static std::vector<std::string> split(const char& delimiter, const std::string& line);
	static std::vector<std::string> splitCSV(const char& delimiter, const std::string& line);
	static std::vector<std::string> findSubdirectories(std::string indir);
	static std::vector<std::string> findFiles(std::string indir, std::string match);
	static std::vector<std::string> buildVector(std::string dir, std::string * filearr, int numoffiles);
	static std::vector<std::string> buildVector(std::string * filearr, int numoffiles);
	static void dbf2csv(std::string shpfile,std::string csvfile);
	static std::string getProjFromShapefile(std::string filename);
	static std::string int2string(int num);
};

class SubdirManager
{
public:
	std::string m_dirpath;
	SubdirManager(std::string dir);
	std::vector<std::string> findFilesMatch(std::string match);
	std::vector<std::string> findFilesMatch(std::vector<std::string> matches);
};