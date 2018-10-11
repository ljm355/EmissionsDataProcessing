#pragma once
#include "Utils.h"
#include "Grid.h"
#include "ShapeFile.h"
#include <netcdf.h>
class Preprocessor
{
public:
	//EXTERNL int
	//nc_get_att_text(int ncid, int varid, const char *name, char *ip);
	struct NetCDFAttribute
	{
		char name[255];
		char text[255];
		size_t len;
		nc_type vtype;
	};

	struct NetCDFVar
	{
		char name[255];
		int id;
		int ndims;
		int dimids[4];
		int varnatt;
		nc_type vtype;
		std::vector<NetCDFAttribute> attributes;
	};
	struct NetCDFDim
	{
		char name[255];
		size_t len;
		int id;
	};
	Preprocessor();
	~Preprocessor();
	static void updateFieldAfterIntersection(std::string filename, double scale = 1.0);
	static void updateFieldAfterIntersection(std::string filename, std::vector<std::string> fieldnames, double scale = 1.0);
	static void updateFieldAfterIntersectionForDir(std::string indir, double scale = 1.0);
	static void gridFolderByRaster(std::string indir, std::string outdir, std::string fishnetrasterfile, double scale = 1.0);
	static void gridFolderByShape(std::string indir, std::string outdir, std::string fishneshapefile, double scale = 1.0);
	static void reprojectWithArcGIS(std::string inputfile, std::string outputfile, std::string pythonScriptFile);
	static void reprojectDir(std::string indir, std::string outdir, std::string pythonFile);
	static void intersectWithArcGIS(std::string inputfile, std::string fishnet, std::string outputfile);
	static void clipWithArcGIS(std::string inputfile1, std::string inputfile2, std::string outputfile);
	static void intersectWithFishnet(Grid* grid, const std::string& fishnetfile, const std::string& inputfile, const std::string& outputfile);
	static bool convertPolyline(std::string infile, std::string outfile);
	static bool convertPolylineGeometry(OGRGeometry* geom, OGRGeometry*& newgeom);
	static OGRMultiLineString* convertOGRMultiLineStringFrom25D(OGRGeometry* geom);
	static OGRLineString* convertOGRLineStringFrom25D(OGRGeometry* geom);
	static bool has25D(OGRGeometry* geom);
	static std::string getFootprintFieldName(OGRwkbGeometryType gtype);
	static void gridShapeFile(std::string infile, std::string outfile, std::string fishnetfile, const char* fieldname = "ca11");
	//static void readFromZipFile(const char * zipFile, std::vector<float*> dataArray, int xsize);
	static void readHourlyFFDASNC(const char * filename, std::vector<float*> dataArrays);
	static void extractFFDASHourly(std::string dir, std::string outfile, double* adftransform, int startrow, int nrows, int startcol, int ncols);
	static void extractEDGASHourly(std::string infile, std::string outfile, double* adftransform, int startrow, int nrows, int startcol, int ncols);
	static void extractGlobalRaster(std::string infile, std::string outfile, double * adftransform, int startrow, int nrows, int startcol, int ncols);
	static std::vector<float*> readData(const char * filename, int & xsize, int & ysize);
	static float *readData(const char* filename);
};


