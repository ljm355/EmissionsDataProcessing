#pragma once
#include "Utils.h"
#include "ShapeFile.h"
#include "ogrsf_frmts.h"
#include <QFileinfo>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <gdal_priv.h>
#include "qdir.h"
#include "qpainter.h"
#include "ExtractByPolygon.h"
#include "Grid.h"

struct FeatureFootprintCount
{
	int fid;
	int cellid;
	double cellvalue;
	double area;
	double totalarea;
	FeatureFootprintCount()
	{

	}
	FeatureFootprintCount(int _fid, int _cellid, double _area)
	{
		fid = _fid; cellid = _cellid; area = _area; totalarea = _area;
	}

};

class TileRenderer
{
public:
	OGREnvelope envelope;
	double topleftx;
	double toplefty;
	int tilesizex;
	int tilesizey;
	int tilecols;
	int tilerows;
	double cellsize;
	int  ncols;
	int  nrows;
	int numtiles;
	void init(const double& _cellsize, OGREnvelope _envelope)
	{
		envelope = _envelope;
		cellsize = _cellsize;
		ncols = (int)(ceil((envelope.MaxX - envelope.MinX) / cellsize));
		nrows = (int)(ceil((envelope.MaxY - envelope.MinY) / cellsize));
		int minsize = 10;
		int maxsize = 8192;

		while (ncols < 100 || nrows < 100)
		{
			cellsize = cellsize / 2;
			ncols = (int)((envelope.MaxX - envelope.MinX) / cellsize);
			nrows = (int)((envelope.MaxY - envelope.MinY) / cellsize);
		}
		tilesizex = maxsize;
		tilesizey = maxsize;
		if (ncols <= maxsize && nrows <= maxsize)
		{
			tilesizex = ncols;
			tilesizey = nrows;
		}
		tilecols = (int)(ceil((double)ncols / (double)tilesizex));
		tilerows = (int)(ceil((double)nrows / (double)tilesizey));
		numtiles = tilecols * tilerows;
		topleftx = envelope.MinX;
		toplefty = envelope.MaxY;

	}
	void drawPolygon(QPainter* painter, const QBrush& brush, const int& xsize, const int& ysize, const OGRPolygon *polygon)
	{
		QPainterPath outerpath;
		QPainterPath innerpath;
		const OGRLinearRing* exteriorRing = polygon->getExteriorRing();
		QPolygonF outerPolygon;
		for (size_t i = 0; i < exteriorRing->getNumPoints(); i++)
		{
			OGRPoint pt;
			exteriorRing->getPoint(i, &pt);
			outerPolygon.push_back(QPointF(pt.getX(), pt.getY()));
			//printf("%f,%f\n", pt.getX(), pt.getY());
		}
		outerpath.addPolygon(outerPolygon);
		for (size_t n = 0; n < polygon->getNumInteriorRings(); n++)
		{
			QPolygonF innerPolygon;
			const OGRLinearRing* interiorRing = polygon->getInteriorRing(n);
			for (size_t i = 0; i < interiorRing->getNumPoints(); i++)
			{
				OGRPoint pt;
				interiorRing->getPoint(i, &pt);
				innerPolygon.push_back(QPointF(pt.getX(), pt.getY()));
			}
			innerpath.addPolygon(innerPolygon);
			//path.addPolygon(outerPolygon);
		}
		if (polygon->getNumInteriorRings())
			outerpath = outerpath.subtracted(innerpath);
		painter->fillPath(outerpath, brush);

	}
	void drawPolyline(QPainter* painter, const QBrush& brush, const int& xsize, const int& ysize, const OGRLineString *polyline)
	{
		QVector<QPointF> qpoints;
		painter->setPen(Qt::black);
		for (size_t i = 0; i < polyline->getNumPoints() - 1; i++)
		{
			OGRPoint pt;
			polyline->getPoint(i, &pt);
			qpoints.push_back(QPointF(pt.getX(), pt.getY()));
			polyline->getPoint(i + 1, &pt);
			qpoints.push_back(QPointF(pt.getX(), pt.getY()));
		}

		painter->drawLines(qpoints);

	}
	RasterMask* drawTile(int tileid, OGRGeometry* geom, bool ispolygon = false)
	{

		int tiley = tileid / tilecols;
		int tilex = tileid % tilecols;

		if (tileid > numtiles - 1)
			return NULL;
		OGREnvelope tilebound;
		tilebound.MinX = topleftx + tilex * tilesizex * cellsize;
		tilebound.MaxY = toplefty - tiley * tilesizey * cellsize;
		tilebound.MaxX = tilebound.MinX + tilesizex * cellsize;
		tilebound.MinY = tilebound.MaxY - tilesizey * cellsize;
		QImage* qImage = new QImage(tilesizex, tilesizey, QImage::Format_ARGB32);
		qImage->fill(qRgba(0, 0, 0, 0));
		QBrush br(Qt::white);
		QPainter painter(qImage);
		painter.begin(qImage);
		//painter.setRenderHint(QPainter::Antialiasing); //make it look nicer
		painter.setRenderHint(QPainter::Antialiasing); //make it look nicer
													   //Antialiasing = 0x01,
													   //TextAntialiasing = 0x02,
													   //SmoothPixmapTransform = 0x04,
													   //HighQualityAntialiasing = 0x08,
													   //NonCosmeticDefaultPen = 0x10

		QTransform transform = QTransform().translate(-tilebound.MinX, -tilebound.MaxY) * QTransform().scale((float)tilesizex / (tilebound.MaxX - tilebound.MinX), -(float)tilesizey / (tilebound.MaxY - tilebound.MinY));
		painter.setTransform(transform);
		if (ispolygon)
		{
			painter.setBrush(br);
			drawPolygon(&painter, br, tilesizex, tilesizey, (OGRPolygon*)geom);
			painter.end();
		}
		else
		{
			drawPolyline(&painter, br, tilesizex, tilesizey, (OGRLineString*)geom);
			painter.end();
		}
		RasterMask* mask = new RasterMask(qImage, tilebound.MinX, tilebound.MinY, tilebound.MaxX, tilebound.MaxY);
		return mask;
	}
};


class RasterizationGridder
{
public:
	RasterizationGridder();
	~RasterizationGridder();
	void gatherFromPolygonMap(RasterMask* mask, Grid* grid, int fid, std::map<int, FeatureFootprintCount>& featureCountMap)
	{


		size_t stride = mask->_stride;
		double* geoTransform = mask->_adfGeoTransform;
		double maskCellSize = abs(geoTransform[5]);

		size_t num = 0;
		size_t numpixels = mask->_ysize*mask->_xsize;
		//uchar *pData = mask->_image->bits();
		unsigned int* pData = (unsigned int*)mask->_image->bits();
		//int whiteidx = 0;
		/*for (size_t i = 0; i < numpixels; i++)
		{
		if (pData[0] == 255 && pData[3] == 255)
		{
		num++;
		}
		pData += stride;
		}*/
		//4294967295
		for (size_t i = 0; i < numpixels; i++)
		{
			if (*pData == 4294967295)
				num++;
			pData++;
		}

		pData = (unsigned int*)mask->_image->bits();
		double fraction = (double)1 / (double)num;
		double total2 = 0;
		double cellsize = grid->_adfGeoTransform[1];
		int maxcellidx = grid->nrows * grid->ncols;
		for (size_t i = 0; i < mask->_ysize; i++)
		{
			float y = geoTransform[3] + geoTransform[5] * 0.5 + geoTransform[5] * i;
			unsigned int yindex = (unsigned int)((grid->bound.MaxY - y) / cellsize);

			for (size_t j = 0; j < mask->_xsize; j++)
			{
				float x = geoTransform[0] + geoTransform[1] * 0.5 + geoTransform[1] * j;
				unsigned int xindex = (unsigned int)((x - grid->bound.MinX) / cellsize);
				if (*pData == 4294967295)
				{
					//total2 += fraction * total;
					int cellidx = xindex + yindex * grid->ncols;
					if (cellidx > 0 && cellidx < maxcellidx)
					{
						std::map<int, FeatureFootprintCount>::iterator findIter = featureCountMap.find(cellidx);

						if (findIter != featureCountMap.end())
						{
							findIter->second.area += maskCellSize;
						}
						else
						{
							FeatureFootprintCount fea;
							fea.cellid = cellidx;
							fea.area = maskCellSize;
							fea.fid = fid;
							featureCountMap[cellidx] = fea;
						}
						//grid->cells[xindex + yindex * grid->ncols] += (fraction * total);
						//grid->cells[xindex + yindex * grid->ncols] ++;
						//grid->cells[xindex + yindex * grid->ncols] += fraction;

					}

				}
				pData++;
			}
		}
		//printf("%f,%f\n", total, total2);
		//if (total2 / total > 0.95)
		//	return true;
		//return false;
	}
	void gatherCellsFromPolygonMap(RasterMask* mask, Grid* grid, std::vector<unsigned long>& cells)
	{

		size_t stride = mask->_stride;
		double* geoTransform = mask->_adfGeoTransform;
		double maskCellSize = abs(geoTransform[5]);

		size_t num = 0;
		size_t numpixels = mask->_ysize*mask->_xsize;
		//uchar *pData = mask->_image->bits();
		unsigned int* pData = (unsigned int*)mask->_image->bits();

		for (size_t i = 0; i < numpixels; i++)
		{
			if (*pData == 4294967295)
				num++;
			pData++;
		}

		pData = (unsigned int*)mask->_image->bits();
		double fraction = (double)1 / (double)num;
		double total2 = 0;
		double cellsize = grid->_adfGeoTransform[1];
		int maxcellidx = grid->nrows * grid->ncols;
		for (size_t i = 0; i < mask->_ysize; i++)
		{
			float y = geoTransform[3] + geoTransform[5] * 0.5 + geoTransform[5] * i;
			unsigned long yindex = (unsigned long)((grid->bound.MaxY - y) / cellsize);

			for (size_t j = 0; j < mask->_xsize; j++)
			{
				float x = geoTransform[0] + geoTransform[1] * 0.5 + geoTransform[1] * j;
				unsigned long xindex = (unsigned long)((x - grid->bound.MinX) / cellsize);
				if (*pData == 4294967295)
				{
					unsigned long cellidx = xindex + yindex * grid->ncols;
					if (cellidx >= 0 && cellidx < maxcellidx)
					{
						cells.push_back(cellidx);
					}

				}
				pData++;
			}
		}
	}
	void gatherFromPolyline(const OGRLineString  *polyline, double& sampleDist, Grid* grid, int fid, std::map<int, FeatureFootprintCount>& featureCountMap)
	{
		unsigned int maxcellidx = grid->nrows * grid->ncols;
		double cellsize = grid->_adfGeoTransform[1];

		//std::vector<OGRPoint> points;
		//for (int i = 0; i < polyline->getNumPoints() - 1; i++)
		//{
		//	OGRPoint pt;
		//	PT->getNextPoint(&pt);
		//}
		for (int i = 0; i < polyline->getNumPoints() - 1; i++)
		{
			OGRPoint pt1;
			OGRPoint pt2;
			polyline->getPoint(i, &pt1);
			polyline->getPoint(i + 1, &pt2);
			//OGRPoint pt1 = points[i];
			//OGRPoint pt2 = points[i + 1];
			double xdif = pt2.getX() - pt1.getX();
			double ydif = pt2.getY() - pt1.getY();
			double dist = sqrt(xdif*xdif + ydif*ydif);
			int ndivisions = ceil(dist / sampleDist);
			double adjustedSampleDist = dist / ndivisions;
			double xstep = xdif / ndivisions;
			double ystep = ydif / ndivisions;

			double p1x = pt1.getX();
			double p1y = pt1.getY();
			for (int n = 0; n < ndivisions + 1; n++)
			{
				double px = p1x + xstep * n;
				double py = p1y + ystep * n;
				unsigned int xindex = (unsigned int)((px - grid->bound.MinX) / cellsize);
				unsigned int yindex = (unsigned int)((grid->bound.MaxY - py) / cellsize);
				int cellidx = xindex + yindex * grid->ncols;
				if (cellidx > 0 && cellidx < maxcellidx)
				{
					std::map<int, FeatureFootprintCount>::iterator findIter = featureCountMap.find(cellidx);

					if (findIter != featureCountMap.end())
					{
						findIter->second.area += adjustedSampleDist;
					}
					else
					{
						FeatureFootprintCount fea;
						fea.cellid = cellidx;
						fea.area = adjustedSampleDist;
						fea.fid = fid;
						featureCountMap[cellidx] = fea;
					}

				}
			}

		}

	}
	void intersect(Grid* grid, std::string inputfile)
	{
		ShapeFile input(inputfile.data());
		/*char wkt[512];
		char* pwkt = wkt;
		if (input.poLayer->GetSpatialRef())
		input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);*/
		int fieldIdx = input.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
		if (fieldIdx < 0)
			return;
		double cellsize = grid->_adfGeoTransform[1];
		OGRFeature* poFeature;
		size_t fid = -1;
		TileRenderer renderer;
		OGREnvelope bound = grid->bound;

		double RASTERIZATION_RESOLUTION = cellsize / 10;
		double LINE_SAMPLE_DIST = cellsize / 100;
		while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
		{
			fid++;
			double value = poFeature->GetFieldAsDouble(fieldIdx);
			OGRGeometry* poGeometry = poFeature->GetGeometryRef();
			if (!poGeometry)
			{
				OGRFeature::DestroyFeature(poFeature);
				continue;
			}
			OGRwkbGeometryType gtype = poGeometry->getGeometryType();
			if (gtype == OGRwkbGeometryType::wkbPoint || gtype == OGRwkbGeometryType::wkbPoint25D)
			{
				OGRPoint* pt = (OGRPoint*)poGeometry;
				unsigned int row = (unsigned int)((bound.MaxY - pt->getY()) / cellsize);
				unsigned int col = (unsigned int)((pt->getX() - bound.MinX) / cellsize);
				unsigned int cellidx = col + row*grid->ncols;
				if (cellidx < grid->slice)
					grid->cells[cellidx] += value;
			}
			else
			{
				bool ispolygon = false;
				if (gtype == wkbPolygon || gtype == wkbPolygon25D)
					ispolygon = true;
				std::vector<OGRGeometry*> geomlist;
				if (dynamic_cast<OGRGeometryCollection*>(poGeometry))
				{
					OGRGeometryCollection* geomcollect = dynamic_cast<OGRGeometryCollection*>(poGeometry);
					for (int igeom = 0; igeom < geomcollect->getNumGeometries(); igeom++)
					{
						geomlist.push_back(geomcollect->getGeometryRef(igeom));
					}
				}
				else
				{
					geomlist.push_back(poGeometry);
				}
				std::map<int, FeatureFootprintCount> feamap;

				for (int igeom = 0; igeom < geomlist.size(); igeom++)
				{
					OGRGeometry* poGeometry = geomlist[igeom];
					OGREnvelope geobb;
					poGeometry->getEnvelope(&geobb);
					if (!bound.Intersects(geobb))
						continue;
					if (ispolygon)
					{
						renderer.init(RASTERIZATION_RESOLUTION, geobb);
						for (int itile = 0; itile < renderer.numtiles; itile++)
						{
							RasterMask* mask = renderer.drawTile(itile, poGeometry, ispolygon);
							gatherFromPolygonMap(mask, grid, fid, feamap);
							////uchar* data = mask->_image->bits();
							////for (int ipixel = 0; ipixel < mask->_xsize * mask->_ysize; ipixel++)
							////{
							////	if (data[1] == 255)
							////	{
							////		data[1] = rand() % 255;
							////		data[2] = rand() % 255;
							////		data[3] = rand() % 255;
							////	}
							////	data += 4;
							////}
							////std::stringstream ss;
							////ss << outdir << fid << "_" << igeom <<  "_" << itile <<  ".tif";
							////mask->save(ss.str());
							delete mask;
						}
					}
					else
					{
						if (dynamic_cast<OGRLineString*>(poGeometry))
						{
							OGRLineString* linestr = (OGRLineString *)poGeometry;
							gatherFromPolyline(linestr, LINE_SAMPLE_DIST, grid, fid, feamap);
						}
						else if (dynamic_cast<OGRCircularString*>(poGeometry))
						{
							OGRCircularString* circularstr = (OGRCircularString *)poGeometry;
							OGRLineString* linestr = dynamic_cast<OGRLineString*>(circularstr->getLinearGeometry());
							gatherFromPolyline(linestr, LINE_SAMPLE_DIST, grid, fid, feamap);
						}

					}


				}
				std::map<int, FeatureFootprintCount>::iterator feamapIter = feamap.begin();
				double totalarea = 0;

				while (feamapIter != feamap.end())
				{
					totalarea += feamapIter->second.area;
					feamapIter++;
				}

				feamapIter = feamap.begin();

				while (feamapIter != feamap.end())
				{
					FeatureFootprintCount& feacount = feamapIter->second;
					feacount.totalarea = totalarea;
					grid->cells[feacount.cellid] += (feacount.area / totalarea * value);
					//intersectedFeaturesCount.push_back(feacount);
					feamapIter++;
				}
			}

			OGRFeature::DestroyFeature(poFeature);
			if (fid % 100 == 0)
				printf("%d\n", fid);
		}
	}
	void intersect(Grid* grid, std::string inputfile, std::string outputfile)
	{
		ShapeFile input(inputfile.data());
		char wkt[512];
		char* pwkt = wkt;
		if (input.poLayer->GetSpatialRef())
			input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		intersect(grid, inputfile);
		grid->toRaster(outputfile + ".tif", pwkt);
		grid->toShape(pwkt, outputfile + ".shp", true);
	}
	void intersectFolder(std::string indir, std::string fishnetfile, std::string outdir)
	{

		QDir qoutdir(outdir.data());
		if (!qoutdir.exists())
			qoutdir.mkpath(".");
		outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();

		std::vector<std::string> files;
		QDir input_dir(indir.data());
		input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
		input_dir.setSorting(QDir::Name);
		indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
		QFileInfoList list = input_dir.entryInfoList();

		char wkt[512];
		char* pwkt = NULL;
		Grid* sumGrid = NULL;
		for (int i = 0; i < list.size(); ++i) {
			QFileInfo fileInfo = list.at(i);
			std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
			if (!fileInfo.fileName().endsWith(".shp", Qt::CaseInsensitive))
				continue;
			std::string inputfile = fileInfo.absoluteFilePath().toLocal8Bit().data();
			if (!pwkt)
			{
				ShapeFile input(inputfile.data());
				pwkt = wkt;
				if (input.poLayer->GetSpatialRef())
					input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
				input.close();
			}

			Grid grid;

			std::string outputfile = outdir + fileInfo.completeBaseName().toLocal8Bit().data() + ".tif";
			if (QFileInfo(outputfile.data()).exists())
			{
				grid.fromFishnetRaster(outputfile, true);
			}
			else
			{
				grid.fromFishnetRaster(fishnetfile);
				grid.reset();
				intersect(&grid, inputfile);
			}

			if (sumGrid == NULL)
			{
				sumGrid = new Grid;
				sumGrid->fromFishnetRaster(fishnetfile);
				sumGrid->reset();
			}

			for (int icell = 0; icell < grid.slice; icell++)
			{
				double cellval = grid.cells[icell];
				if (cellval <= 0)
					continue;
				sumGrid->cells[icell] += cellval;
			}
			//
			//grid->toShape(pwkt, outputfile + ".shp", true);

			//grid.toRaster(outputfile, pwkt);
			printf("%s\n", outputfile.data());
		}
		for (int icell = 0; icell < sumGrid->slice; icell++)
		{
			double cellval = sumGrid->cells[icell];
			//if (cellval <= 0)
			//	continue;
			sumGrid->cells[icell] = log10(cellval);
		}
		sumGrid->toRaster(outdir + "total.tif", pwkt);
		delete sumGrid;
	}
	void intersect(OGREnvelope bound, double gridcellsize, std::string inputfile, std::string outputfile)
	{
		Grid grid(bound, gridcellsize, 0);
		grid.reset();
		intersect(&grid, inputfile, outputfile);
	}
	void intersect(std::string fishnet_rasterfile, std::string inputfile, std::string outputfile)
	{

		Grid grid;
		grid.fromFishnetRaster(fishnet_rasterfile);
		grid.reset();
		intersect(&grid, inputfile, outputfile);
	}

	std::vector<double> readAttribute(std::string shapefile, std::string fieldname)
	{
		std::vector<double> values;
		ShapeFile input(shapefile.data());
		int fieldIdx = input.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data());
		OGRFeature* poFeature;
		while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
		{

			OGRFeature::DestroyFeature(poFeature);
			values.push_back(poFeature->GetFieldAsDouble(fieldIdx));
		}
		return values;
	}
	void compareFiles(std::string file1, std::string file2, std::string fieldname)
	{

		std::vector<double> values1 = readAttribute(file1, fieldname);
		std::vector<double> values2 = readAttribute(file2, fieldname);
		double max = -999999999;
		double min = 999999999;
		double sum = 0;
		int count = 0;
		for (size_t i = 0; i < values1.size(); i++)
		{
			if (isnan(values1[i]) || isnan(values2[i]))
				continue;
			if (values1[i] == 0)
				continue;
			double dif = ((abs)(values2[i] - values1[i])) / values1[i];
			if (isnan(dif) || dif != dif)
				continue;
			if (max < dif)
				max = dif;
			if (min > dif)
				min = dif;
			if (dif > 0.1)
			{
				printf("%d\n", i);
			}

			sum += dif;
			count++;
		}

		double mean = sum / count;
		printf("mean=%f,min=%f,max=%f", mean, min, max);
	}
	void createPolygon(const std::string& outputfile)
	{

		ShapeFile polygonShapefile;
		polygonShapefile.create(outputfile.data());

		OGRFeatureDefn *poFDefn = polygonShapefile.poLayer->GetLayerDefn();
		OGRFieldDefn def("ca11", OGRFieldType::OFTReal);
		polygonShapefile.poLayer->CreateField(&def);
		OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(polygonShapefile.poLayer->GetLayerDefn());

		OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
		OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);

		double radius = 0.5;
		for (int i = 0; i < 360; i++)
		{
			float degreeInRad = i*0.0174533;
			double x = cos(degreeInRad)*radius;
			double y = sin(degreeInRad)*radius;
			linearRing->addPoint(x, y);
		}
		poPolygon->addRing(linearRing);
		poFeaPolygon->SetGeometry(poPolygon);
		poFeaPolygon->SetField("ca11", 1);
		polygonShapefile.poLayer->CreateFeature(poFeaPolygon);

		OGRFeature::DestroyFeature(poFeaPolygon);

	}
};

