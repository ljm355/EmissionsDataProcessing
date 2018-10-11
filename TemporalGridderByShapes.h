#pragma once
#include "TemporalGridding.h"

class FFCO2Sector
{
public:
	std::string sectorname;
	double* cells;
	std::map<std::string, double*> timestructMap;
	int numhours;
	int numdays;
	int numcells;
	std::string fieldname;
	std::vector<FeatureCell> features;
	std::vector<std::string> shapefiles;
	virtual void loadAttribute(std::string _fieldname);
	virtual void getTimeSlice(double* gridcells, int houridx);
	virtual void getTotal(double* gridcells);
	virtual void getTimeSlice(int houridx);
	virtual void getTotal();
	virtual void makeHourlyTotal(std::string hourlyFilename);
	virtual void makeDailyTotal(std::string dailyFilename);
	virtual void makeAnnualTotal(std::string totalFilename);
	virtual void addCells(double* grid1, double* grid2, int numcells);
	FFCO2Sector(int _numhours)
	{
		numhours = _numhours;
		numdays = numhours / 24;
	}
};

class TemporalGridderByShapes : public FFCO2Sector
{
public:
	char* idbuf;
	double* fractionbuf;
	std::vector<TimeStructContainer*> timestructs;
	TemporalGridderByShapes(int _numhours)
		:FFCO2Sector(_numhours)
	{

		idbuf = NULL;
		fractionbuf = NULL;
	}
	~TemporalGridderByShapes();
	std::vector<FFCO2Sector*> sectors;
	virtual void loadAttribute(std::string fieldname);
	virtual void normalizeFractions(double*& fractions);
	virtual void clearSectors();
	virtual void loadtimestruct(std::string txtdir, std::string binaryfile);
	virtual void loadtimestruct(std::string binaryfile);
	virtual void clearTimeStruct();
	virtual void addSectorGrid(std::vector<std::string> shapefiles, std::string sectorname);
	virtual void getTotal();
	virtual void makeHourlyTotal(std::string hourlyFilename);
	virtual void makeDailyTotal(std::string dailyFilename);
	virtual void makeAnnualTotal(std::string totalFilename);
};

