#include "Grid.h"

Grid::Grid(const OGREnvelope& bb, const double& gridcellsize, int expandByCells)
{
	//ShapeFile input(inputfile.data());
	ncols = (bb.MaxX - bb.MinX) / gridcellsize + expandByCells;
	if (gridcellsize * ncols < (bb.MaxX - bb.MinX))
	{
		ncols = ncols + 1;
	}

	nrows = (bb.MaxY - bb.MinY) / gridcellsize + expandByCells;

	//printf("%f,%f\n", bb.MaxX - bb.MinX, bb.MaxY - bb.MinY)
	if (gridcellsize * nrows < (bb.MaxY - bb.MinY))
	{
		nrows = nrows + 1;
	}

	OGREnvelope newbound;
	newbound.MinX = bb.MinX - ((ncols * gridcellsize) - (bb.MaxX - bb.MinX)) * 0.5;
	newbound.MaxX = bb.MaxX + ((ncols * gridcellsize) - (bb.MaxX - bb.MinX)) * 0.5;
	newbound.MinY = bb.MinY - ((nrows * gridcellsize) - (bb.MaxY - bb.MinY)) * 0.5;
	newbound.MaxY = bb.MaxY + ((nrows * gridcellsize) - (bb.MaxY - bb.MinY)) * 0.5;
	bound = newbound;
	_adfGeoTransform[0] = newbound.MinX;
	_adfGeoTransform[1] = gridcellsize;
	_adfGeoTransform[2] = 0;
	_adfGeoTransform[3] = newbound.MaxY;
	_adfGeoTransform[4] = 0;
	_adfGeoTransform[5] = -gridcellsize;
	cells = NULL;
	//reset();
}
Grid::Grid(double* adfGeoTransform,int xsize,int yxize)
{
	//ShapeFile input(inputfile.data());
	ncols = xsize;
	nrows = yxize;
	for (size_t i = 0; i < 6; i++)
	{
		_adfGeoTransform[i] = adfGeoTransform[i];
	}
	bound.MinX = _adfGeoTransform[0];
	bound.MaxX = bound.MinX + xsize*_adfGeoTransform[1];
	bound.MaxY = _adfGeoTransform[3];
	bound.MinY = bound.MaxY + yxize*_adfGeoTransform[5];
	cells = NULL;
	//reset();
}


void Grid::fromFishnetShape(std::string fishnetfile)
{
	ShapeFile input(fishnetfile.data());
	input.poLayer->GetExtent(&bound);
	OGRFeature *poFeature;
	input.poLayer->ResetReading();
	OGREnvelope cellBB;
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{

		((OGRPolygon*)poFeature->GetGeometryRef())->getEnvelope(&cellBB);
		OGRFeature::DestroyFeature(poFeature);
		break;
	}

	char wkt[512];
	char* pwkt = wkt;
	if (input.poLayer->GetSpatialRef())
	{
		input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		proj = pwkt;
	}


	double gridcellsize = (cellBB.MaxX - cellBB.MinX);
	ncols = (int)((bound.MaxX - bound.MinX + gridcellsize * 0.5) / gridcellsize);
	nrows = (int)((bound.MaxY - bound.MinY + gridcellsize * 0.5) / gridcellsize);

	_adfGeoTransform[0] = bound.MinX;
	_adfGeoTransform[1] = gridcellsize;
	_adfGeoTransform[2] = 0;
	_adfGeoTransform[3] = bound.MaxY;
	_adfGeoTransform[4] = 0;
	_adfGeoTransform[5] = -gridcellsize;
	cells = NULL;
}


void Grid::fromFishnetRaster(std::string fishnetfile,bool load)
{

	cells = NULL;
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDataset* pDataset = (GDALDataset*)GDALOpen(fishnetfile.data(), GA_ReadOnly);
	pDataset->GetGeoTransform(_adfGeoTransform);
	ncols = pDataset->GetRasterXSize();
	nrows = pDataset->GetRasterYSize();
	proj = pDataset->GetProjectionRef();
	if (load)
	{
		reset();
		pDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, ncols, nrows, cells, ncols, nrows, GDT_Float64, 0, 0);
	}
	GDALClose((GDALDatasetH)pDataset);
	bound.MinX = _adfGeoTransform[0];
	bound.MaxY = _adfGeoTransform[3];
	bound.MaxX = _adfGeoTransform[0] + _adfGeoTransform[1] * ncols;
	bound.MinY = _adfGeoTransform[3] + _adfGeoTransform[5] * nrows;

	
	//reset();
}

Grid::Grid()
{
	cells = NULL;
}


Grid::~Grid()
{
	if (cells)
		delete[] cells;
}
double Grid::calFraction(OGRFeature* fea, const int& footprintIdx)
{

	OGRGeometry* geo = fea->GetGeometryRef();
	OGRwkbGeometryType gtype = geo->getGeometryType();
	if (gtype == wkbLineString || gtype == wkbLineString25D || gtype == wkbMultiLineString)
	{
		double len1 = Utils::calPolylineLength(fea->GetGeometryRef());// *FOOTPRINT_SCALE_FACTOR;
		double len2 = fea->GetFieldAsDouble(footprintIdx);
		return len1 / len2;
	}
	else if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		double area1 = Utils::calPolygonArea(fea->GetGeometryRef());
		double area2 = fea->GetFieldAsDouble(footprintIdx);
		return area1 / area2;
	}

	return 1;
}

int Grid::getNumberOfValidCells()
{
	int count = 0;
	for (size_t i = 0; i < slice; i++)
	{
		double val = cells[i];
		if (val <= 0 || isinf(val) || val != val)
			continue;
		count++;
	}
	return count;
}
void Grid::reset(std::string attributeName)
{
	if (cells)
		delete[] cells;
	cells = new double[ncols*nrows];
	dims.clear();
	dims.push_back(attributeName);
	memset(cells, 0, sizeof(double)*ncols*nrows);
	slice = ncols*nrows;
}
double Grid::sum()
{
	double sumval = 0;
	if(cells == NULL)
		return 0.0;
	for (size_t i = 0; i < slice; i++)
	{
		double val = cells[i];
		if (val < 0 || val != val || isinf(val) > 99999999999)
			continue;
		sumval += val;
	}
	return sumval;
}
void Grid::resetValue(std::string attributeName,double initValue)
{
	if (cells)
		delete[] cells;
	cells = new double[ncols*nrows];
	dims.clear();
	dims.push_back(attributeName);
	memset(cells, 0, sizeof(double)*ncols*nrows);
	slice = ncols*nrows;

	for (size_t i = 0; i < slice; i++)
	{
		cells[i] = initValue;
	}
}
void Grid::reset(std::vector<std::string> dimensions)
{
	dims = dimensions;
	if (cells)
		delete[] cells;
	cells = new double[ncols*nrows*dims.size()];
	memset(cells, 0, sizeof(double)*ncols*nrows*dims.size());
	slice = ncols*nrows;
}

Grid* Grid::createFishnet(const std::string& shapefile, double gridcellsize)
{
	ShapeFile input(shapefile.data());
	OGREnvelope bound;
	input.poLayer->GetExtent(&bound, true);
	Grid* grid = new Grid(bound, gridcellsize);

	return grid;
}

OGRPolygon* Grid::toOGRPolygon(OGRLayer* layer, const OGREnvelope& bb)
{

	OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
	OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);
	linearRing->addPoint(bb.MinX, bb.MinY);
	linearRing->addPoint(bb.MinX, bb.MaxY);
	linearRing->addPoint(bb.MaxX, bb.MaxY);
	linearRing->addPoint(bb.MaxX, bb.MinY);
	linearRing->addPoint(bb.MinX, bb.MinY);

	poPolygon->addRing(linearRing);//also crashed
	return poPolygon;
}
void Grid::normalizedByArea()
{
	double area = _adfGeoTransform[1] * _adfGeoTransform[1];
	for (size_t i = 0; i < slice; i++)
	{
		cells[i] = cells[i] / area;
	}
}

void Grid::toShape(std::string wkt, const std::string& outputfile, bool writeAttribute)
{
	
	ShapeFile fishnet;
	if (wkt != "")
	{
		OGRSpatialReference spatialref;
		char* proj = (char*)wkt.data();
		spatialref.importFromWkt(&proj);
		fishnet.create(outputfile.data(), &spatialref);
	}
	else
	{
		fishnet.create(outputfile.data());
	}

	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();
	double gridcellsize =  _adfGeoTransform[1];
	int idIdx =  fishnet.getOrCreateField("gridId", OGRFieldType::OFTInteger);
	//int colIdx = fishnet.getOrCreateField("col", OGRFieldType::OFTInteger);
	//int rowIdx = fishnet.getOrCreateField("row", OGRFieldType::OFTInteger);

	//int caIndexOutput = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
	std::vector<int> fields;
	if (writeAttribute)
	{
		//if (caIndexOutput < 0)
		//{

		//OGRFieldDefn def("ca11", OGRFieldType::OFTReal);
		//fishnet.poLayer->CreateField(&def);
		//caIndexOutput = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
		//}

		for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
		{
			OGRFieldDefn def(dims[iattribute].data(), OGRFieldType::OFTReal);
			fishnet.poLayer->CreateField(&def);
			fields.push_back(fishnet.poLayer->GetLayerDefn()->GetFieldIndex(dims[iattribute].data()));
		}
	}

	int id = 0;
	int slice = nrows * ncols;
	for (size_t i = 0; i < nrows; i++)
	{
		//std::vector<GridCell>& newrow = gridcells[i];
		OGREnvelope tmpBound;
		tmpBound.MaxY = bound.MaxY - gridcellsize * i;
		tmpBound.MinY = bound.MaxY - gridcellsize * (i + 1);
		for (size_t j = 0; j < ncols; j++)
		{

			OGREnvelope bb = tmpBound;
			bb.MinX = bound.MinX + gridcellsize * j;
			bb.MaxX = bound.MinX + gridcellsize * (j + 1);

		/*	if (bb.MinX >(-133.177071456802 - 0.1) && bb.MaxX < (-104.855640115368 + 0.1) && bb.MinY >(25.1208667507254 - 0.1) && bb.MaxY < (46.2050395140432 + 0.1))
			{*/
				OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());

				poFeaPolygon->SetField(idIdx, id);
	/*			poFeaPolygon->SetField(rowIdx, (int)i);
				poFeaPolygon->SetField(colIdx, (int)j);*/
				bool hasnonzero = false;
				if (writeAttribute)
				{
					//poFeaPolygon->SetField(caIndexOutput, cells[id]);
					for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
					{
						double val = cells[slice*iattribute + id];
						if (val != 0)
							hasnonzero = true;
						poFeaPolygon->SetField(fields[iattribute], cells[slice*iattribute + id]);
					}
				}
				OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bb);
				poFeaPolygon->SetGeometry(poPolygon);

				if (!writeAttribute || (writeAttribute && hasnonzero))
					fishnet.poLayer->CreateFeature(poFeaPolygon);

				OGRFeature::DestroyFeature(poFeaPolygon);
			//}
			
			id++;
		}

	}
	fishnet.close();
}
void Grid::toShape(std::string wkt, const std::string& outputfile, OGREnvelope destBound, bool writeAttribute)
{

	ShapeFile fishnet;
	if (wkt != "")
	{
		OGRSpatialReference spatialref;
		char* proj = (char*)wkt.data();
		spatialref.importFromWkt(&proj);
		fishnet.create(outputfile.data(), &spatialref);
	}
	else
	{
		fishnet.create(outputfile.data());
	}

	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();
	double gridcellsize = _adfGeoTransform[1];
	//int idIdx = fishnet.getOrCreateField("gridId", OGRFieldType::OFTInteger);
	std::vector<int> fields;
	if (writeAttribute)
	{
		for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
		{
			OGRFieldDefn def(dims[iattribute].data(), OGRFieldType::OFTReal);
			fishnet.poLayer->CreateField(&def);
			fields.push_back(fishnet.poLayer->GetLayerDefn()->GetFieldIndex(dims[iattribute].data()));
		}
	}

	int id = 0;
	int slice = nrows * ncols;
	for (size_t i = 0; i < nrows; i++)
	{
		OGREnvelope tmpBound;
		tmpBound.MaxY = bound.MaxY - gridcellsize * i;
		tmpBound.MinY = bound.MaxY - gridcellsize * (i + 1);
		for (size_t j = 0; j < ncols; j++)
		{

			OGREnvelope bb = tmpBound;
			bb.MinX = bound.MinX + gridcellsize * j;
			bb.MaxX = bound.MinX + gridcellsize * (j + 1);
			if (destBound.Intersects(bb) || destBound.Contains(bb))
			{
				OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());
				//poFeaPolygon->SetField(idIdx, id);
				bool hasnonzero = false;
				if (writeAttribute)
				{
					for (size_t iattribute = 0; iattribute < dims.size(); iattribute++)
					{
						double val = cells[slice*iattribute + id];
						poFeaPolygon->SetField(fields[iattribute], cells[slice*iattribute + id]);
					}
				}
				OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bb);
				poFeaPolygon->SetGeometry(poPolygon);

				fishnet.poLayer->CreateFeature(poFeaPolygon);
				OGRFeature::DestroyFeature(poFeaPolygon);
			}
			id++;
		}

	}
	fishnet.close();
}
void Grid::toShape(std::string wkt, const std::string& outputfile, int cellID)
{

	ShapeFile fishnet;
	if (wkt != "")
	{
		OGRSpatialReference spatialref;
		char* proj = (char*)wkt.data();
		spatialref.importFromWkt(&proj);
		fishnet.create(outputfile.data(), &spatialref);
	}
	else
	{
		fishnet.create(outputfile.data());
	}

	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();
	double gridcellsize = _adfGeoTransform[1];
	int idIdx = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("gridId");

	OGRFieldDefn def("gridId", OGRFieldType::OFTInteger);
	fishnet.poLayer->CreateField(&def);
	idIdx = fishnet.poLayer->GetLayerDefn()->GetFieldIndex("gridId");

	int nrow = cellID / ncols;
	int ncol = cellID % ncols;
	int id = 0;

	OGREnvelope bb;
	bb.MaxY = bound.MaxY - gridcellsize * nrow;
	bb.MinY = bound.MaxY - gridcellsize * (nrow + 1);
	bb.MinX = bound.MinX + gridcellsize * ncol;
	bb.MaxX = bound.MinX + gridcellsize * (ncol + 1);

	OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());

	poFeaPolygon->SetField(idIdx, id);
	OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bb);
	poFeaPolygon->SetGeometry(poPolygon);

	fishnet.poLayer->CreateFeature(poFeaPolygon);

	OGRFeature::DestroyFeature(poFeaPolygon);
	fishnet.close();
}
void Grid::toBoundaryShape(const std::string& outputfile,const std::string& georeffile)
{

	ShapeFile fishnet;
	if (georeffile != "")
	{
		ShapeFile copyfrom(georeffile.data());
		fishnet.create(outputfile.data(), copyfrom.poLayer->GetSpatialRef());
		copyfrom.close();
	}
	else
	{
		fishnet.create(outputfile.data());
	}
	OGRFeatureDefn *poFDefn = fishnet.poLayer->GetLayerDefn();

 
	OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());

	OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, bound);
	poFeaPolygon->SetGeometry(poPolygon);
	fishnet.poLayer->CreateFeature(poFeaPolygon);
	OGRFeature::DestroyFeature(poFeaPolygon);
	fishnet.close();


}

void Grid::toRaster(std::string filename, std::string wkt)
{
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset* pDataset = poDriver->Create(filename.data(), ncols, nrows, 1, GDT_Float64, papszOptions);
	pDataset->SetGeoTransform(_adfGeoTransform);
	double* pValue = cells;
	double _nodata = 0;
	GDALRasterBand *pBand = pDataset->GetRasterBand(1);
	if (cells != NULL)
	{
		for (size_t i = 0; i < ncols * nrows; i++)
		{
			double val = cells[i];
			if (val <= 0)
				cells[i] = _nodata;
		}
		pBand->RasterIO(GF_Write, 0, 0, ncols, nrows, cells, ncols, nrows, GDT_Float64, 0, 0);
	}
	pBand->SetNoDataValue(-9999);
	if (wkt == "")
		wkt == proj;
	//if (wkt != NULL)
	//{
		pDataset->SetProjection(wkt.data());
	//}
	GDALClose((GDALDatasetH)pDataset);
}

Grid* Grid::upscale(int scale)
{
	Grid* upscaledGrid = new Grid(bound, _adfGeoTransform[1] * scale, 0);
	upscaledGrid->reset();
	int idx = 0;
	for (int nrow = 0; nrow < nrows; nrow++)
	{
		int nrow2 = nrow / scale;
		for (int ncol = 0; ncol < ncols; ncol++)
		{
			int ncol2 = ncol / scale;
			upscaledGrid->cells[ncol2 + nrow2*upscaledGrid->ncols] += cells[idx++];
		}
	}
	return upscaledGrid;
}
void Grid::gatherCells(ShapeFile* input,const char* attributeID, double scale)
{

	OGRFeature *poFeature;
	input->poLayer->ResetReading();
	int caIndexOld = input->poLayer->GetLayerDefn()->GetFieldIndex(attributeID);
	int idIndex = input->poLayer->GetLayerDefn()->GetFieldIndex("gridId");
	if (caIndexOld == -1)
		return;
	OGRwkbGeometryType gtype = input->poLayer->GetGeomType();
	int footprintIndexOld = -1;
	int footprintIndexNew = -1;
	if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
	{
		footprintIndexOld = input->poLayer->GetLayerDefn()->GetFieldIndex("length");
	}
	if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		footprintIndexOld = input->poLayer->GetLayerDefn()->GetFieldIndex("area");
	}
	while ((poFeature = input->poLayer->GetNextFeature()) != NULL)
	{
		

		double fraction = calFraction(poFeature, footprintIndexOld) * scale;
		double ca = poFeature->GetFieldAsDouble(caIndexOld);

		if (ca > 0)
		{
			int id = poFeature->GetFieldAsInteger(idIndex);
			cells[id] += ca;// *fraction;

		}
		OGRFeature::DestroyFeature(poFeature);
	}

}
void Grid::gatherCells(std::string shapefile, const char* attributeID, double scale)
{
	ShapeFile input(shapefile);
	OGRFeature *poFeature;
	input.poLayer->ResetReading();
	int caIndexOld = input.poLayer->GetLayerDefn()->GetFieldIndex(attributeID);
	int idIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("gridId");
	if (caIndexOld == -1)
		return;
	OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
	int footprintIndexOld = -1;
	int footprintIndexNew = -1;
	if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
	{
		footprintIndexOld = input.poLayer->GetLayerDefn()->GetFieldIndex("length");
	}
	if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		footprintIndexOld = input.poLayer->GetLayerDefn()->GetFieldIndex("area");
	}
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{


		//double fraction = calFraction(poFeature, footprintIndexOld) * scale;
		double ca = poFeature->GetFieldAsDouble(caIndexOld);

		if (ca > 0)
		{
			int id = poFeature->GetFieldAsInteger(idIndex);
			cells[id] += ca;// *fraction;

		}
		OGRFeature::DestroyFeature(poFeature);
	}

}
//void Grid::gatherCells(const std::string& inputfile)
//{
//	ShapeFile input(inputfile.data());
//	gatherCells(&input);
//
//}