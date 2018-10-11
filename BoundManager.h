#pragma once
#include "ShapeFile.h"
#include <sstream>
#include <fstream>
#include <iomanip>
class BoundManager
{
public:
	BoundManager();
	~BoundManager();
	static OGREnvelope readBoundFromShape(std::string shapefile)
	{
		ShapeFile input(shapefile);
		OGREnvelope bound;
		input.poLayer->GetExtent(&bound, true);
		return bound;
	}
	static OGREnvelope readBoundFromRaster(std::string rasterfile)
	{
		OGREnvelope bound;
		double adfGeoTransform[6];
		GDALDataset* dataset = (GDALDataset*)GDALOpen(rasterfile.data(), GA_ReadOnly);
		if (!dataset)
			return bound;
		dataset->GetGeoTransform(adfGeoTransform);
		int ncols = dataset->GetRasterXSize();
		int nrows = dataset->GetRasterYSize();

		bound.MinX = adfGeoTransform[0];
		bound.MaxY = adfGeoTransform[3];
		bound.MaxX = adfGeoTransform[0] + adfGeoTransform[1] * ncols;
		bound.MinY = adfGeoTransform[3] + adfGeoTransform[5] * nrows;
		if (dataset)
			GDALClose(dataset);
		dataset = NULL;
		return bound;
	}
	static OGREnvelope readBoundFromShapes(std::vector<std::string> shapefiles)
	{
		bool init = false;
		OGREnvelope masterbb;
		for (size_t i = 0; i < shapefiles.size(); i++)
		{
			ShapeFile input(shapefiles[i]);
			
			OGREnvelope bound;
			input.poLayer->GetExtent(&bound, true);
			double sum = bound.MaxX + bound.MinX + bound.MaxY + bound.MinY;
			//if (sum != sum || isinf(sum) || abs(sum) > 1000000000)
			//{
			//	printf("%s\n", shapefiles[i].data());
			//	continue;
			//}
			/*if (bound.MinX < 0 || bound.MinY < 0 || bound.MaxX < 0 || bound.MaxY < 0)
			{
				printf("%s\n", shapefiles[i].data());
				continue;
			}*/
			if (!init)
			{
				masterbb = bound;
				init = true;
			}
			else
			{
				masterbb.Merge(bound);
			}
		}

		return masterbb;
	}
	static OGREnvelope readBound(std::string filename)
	{
		std::ifstream ifs(filename.data());
		OGREnvelope bound;
		ifs >> bound.MinX >> bound.MaxX >> bound.MinY >> bound.MaxY;
		ifs.close();
		return bound;
	}
	static void writeBound(OGREnvelope bound, std::string outfilename)
	{
		std::ofstream ofs(outfilename.data());
		ofs << std::setprecision(16) << bound.MinX << " " << bound.MaxX << " " << bound.MinY << " " << bound.MaxY;
		ofs.close();
	}
	static void writeBound(std::string shapefile, std::string outfilename)
	{
		OGREnvelope bound = readBoundFromShape(shapefile);
		writeBound(bound, outfilename);
	}

};

