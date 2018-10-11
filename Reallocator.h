#pragma once
#include "Utils.h"
#include "ShapeFile.h"


struct IntersectedFeature
{
	int id;
	double* timestructure;
	double total;
};

class CASector
{
public:
	std::string sectorname;
	std::string cafieldname;
	std::vector<std::string> idfieldnames;
	std::vector<std::string> shapefiles;
	std::map<std::string, double*> timestructMap;
	std::vector<IntersectedFeature> features;
	virtual void loadAttribute(std::string _cafieldname, std::vector<std::string> _idfieldnames);
	std::vector<double> aggregrate(std::string basemapfile,std::string outfile);
};

class Reallocator : public CASector
{
public:
	Reallocator();
	~Reallocator();
	std::string intersectionShapefile;
	std::string idFieldName;
	std::vector<CASector*> sectors;
	char* idbuf;
	double* fractionbuf;
//private:
	void loadAttribute(std::string _cafieldname, std::vector<std::string> _idfieldnames);
	void addSector(std::vector<std::string> shapefiles, std::string sectorname);
	std::vector<double> aggregrate(std::string basemapfile, std::string outfile);
	void loadtimestruct(std::string txtdir, std::string binaryfile,int numofhours);
	void clearSectors();

};

