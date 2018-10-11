#pragma once
#include "Grid.h"
struct GeomInfo
{
	OGRGeometry* geom;
	int id;
	int fishnetid;
	OGREnvelope bb;
	GeomInfo(OGRGeometry* _geom, int _id, OGREnvelope _bb)
	{
		geom = _geom;
		id = _id;
		bb = _bb;
	}
	GeomInfo(OGRGeometry* _geom, int _id, int _fishnetid, OGREnvelope _bb)
	{
		geom = _geom;
		id = _id;
		fishnetid = _fishnetid;
		bb = _bb;
	}
	GeomInfo()
	{
		geom = NULL;
	}

};


typedef std::vector<GeomInfo> GEOMGROUPS;


class MultitheadIntersection
{
public:
	MultitheadIntersection();
	~MultitheadIntersection();
	Grid* grid;
	OGRPolygon* createBB(OGREnvelope bb);
	void setupGrid(OGREnvelope bb,double cellsize);
	void setupGrid(std::string tiff);
	void intersect(std::string inshapefile, std::string outshapefile);
	void intersectParallel(std::string inshapefile, std::string outshapefile, int numThreads);
};

