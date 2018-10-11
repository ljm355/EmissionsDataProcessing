#pragma once
#include "ANNSearchEngine.h"
#include "PeMSStation.h"
#include "ShapeFile.h"
struct TimeProfileCoords{
	double X;
	double Y;
};
class RoadTimeIDW
{
public:
	RoadTimeIDW();
	~RoadTimeIDW();
	double m_maxDist;
	int m_maxnumpoints;
	int m_minnumpoints;
	std::string m_outdir;
	std::map<int,int> m_stationProfileIDs;
	std::map<int, TimeProfileCoords> m_TimeProfileCoords;
	ANNSearchEngine m_searchEngine;
	std::vector<PeMSStation*> m_stations;
	static void getPoints(OGRGeometry * geom, std::vector<OGRPoint>& points);
	void loadStations(std::string stationdir,std::string outputdir);
	void loadStations(std::string shapefile, std::string timestructField, std::string timestructDir);
	void loadCCSStations(std::string shapefile, std::string temporalFile);
	void createIDW();
	void stations2shapes(std::string shapefile);
	void coordinatesFromShapefile(std::string shapefile, std::string idField = "");
	void idw(std::string roadshapefile, std::string fieldname, std::string profilePrefixSector, int idoffset);
	void idw(std::string roadshapefile, std::string fieldnameold, std::string fieldnamenew,std::string profilePrefix, std::string profilePrefixFip, int idoffset = 0);
	void nearest(std::string roadshapefile, std::string fieldname, std::string profilePrefixSector, std::string profilePrefixFip, int idoffset = 0);
	void idw_ccs(std::string roadshapefile, OGREnvelope bound, double resol);

	void fillgap_ccsp(std::string roadshapefile, std::string samplefile);
	
};

