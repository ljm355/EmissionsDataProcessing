#pragma once
#include "Utils.h"
#include "Preprocessor.h"
#include "GDAL_DS.h"


class Vulcan2014
{
public:
	std::string m_inputdir;
	std::string m_outputdir;

	Vulcan2014();
	~Vulcan2014();
	void createAirport();
	void createRailroadPoint();
	void createNonroadPoint();
	void createElecProd();
	void createIndPoint();
	void createComPoint();
	void createCementPoint();
	//void createNonPoints();
	void sumBySectorRaster(std::string dir, std::vector<std::vector<std::string>> sectorfiles, std::vector<std::string> sectornames, std::string outcsvfile);
	void sumBySectorRaster(std::string dir, std::string outcsvfile);
	void sumBySectorShape(std::string dir, std::vector<std::vector<std::string>> sectorfiles, std::vector<std::string> sectornames, std::string outcsvfile);
	void sumBySectorShape(std::string dir, std::string outcsvfile);
};

