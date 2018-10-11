#pragma once
#include "ShapeFile.h"
#include <vector>
#include "gdal_priv.h"
#include "Utils.h"
class Grid
{
public:
	Grid();
	Grid(const OGREnvelope& bb, const double& gridcellsize, int expandByCells = 1);
	Grid(double * adfGeoTransform, int xsize, int yxize);
	void fromFishnetShape(std::string fishnetfile);
	void fromFishnetRaster(std::string fishnetfile, bool load = false);
	~Grid();
	int getNumberOfValidCells();
	void reset(std::string attributeName = "ca11");
	double sum();
	void resetValue(std::string attributeName, double initValue=1);
	void reset(std::vector<std::string> dimensions);
	static Grid* createFishnet(const std::string& shapefile, double gridcellsize);
	OGRPolygon* toOGRPolygon(OGRLayer* layer, const OGREnvelope& bb);
	void normalizedByArea();
	void toShape(std::string wkt, const std::string& outputfile, bool writeAttribute = false);
	void toShape(std::string wkt, const std::string& outputfile, int cellID);
	void Grid::toShape(std::string wkt, const std::string& outputfile, OGREnvelope destBound, bool writeAttribute = false);
	void toBoundaryShape(const std::string & outputfile, const std::string & georeffile);
	void toRaster(std::string filename, std::string wkt = "");
	Grid* upscale(int scale);
	void gatherCells(ShapeFile* input, const char* attributeID = "ca11", double scale = 1.0);
	void gatherCells(std::string shapefile, const char* attributeID = "ca11", double scale = 1.0);
	void intersect(std::string inshapefile, std::string outshapefile);
	//void gatherCells(const std::string& inputfile);
private:
	double calFraction(OGRFeature* fea, const int& footprintIdx);
public:
	double* cells;
	OGREnvelope bound;
	int nrows;
	int ncols;
	int slice;
	double _adfGeoTransform[6];
	std::vector<std::string> dims;
	std::string proj;
};

