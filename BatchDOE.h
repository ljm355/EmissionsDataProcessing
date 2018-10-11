#pragma once
#include <string>
#include "ShapeFile.h"
class BatchDOE
{
public:
	struct WeatherStation
	{
		std::string State;
		std::string Site_Name;
		std::string Elev;
		int Filter;
		std::string DataFileName;
		double Lon;
		double Lat;
		int ID;
	};
	BatchDOE();
	~BatchDOE();
	void init(std::string bdindir, std::string doedir, std::string stationShapeFile);
	void run(std::string bdPrototype, std::string year, WeatherStation& sta);
	void runAll(std::string year);
	void output2shapes();
	void gapfill();
	void findNearestStations(std::string bkFile);
private:
	std::string _bdindir;
	std::string _doedir;
	std::string _stationShapeFile;
	std::vector<WeatherStation> _weatherStations;
	std::vector<std::string> _bdPrototypes;
	std::vector<std::string> readTextFile(std::string filename);
	void writeTextFile(std::vector<std::string>& lines,std::string filename);
	void formatName(std::string& name);
	void copyPrototype(std::string bdPrototype, std::string year , WeatherStation& weatherStation);
	bool parseSimResults(std::string filename, std::vector<std::string>& results);
	void cleanup(std::string dir);
};

