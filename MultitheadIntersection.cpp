#include "MultitheadIntersection.h"



MultitheadIntersection::MultitheadIntersection()
{
	grid = NULL;
}


MultitheadIntersection::~MultitheadIntersection()
{

}

OGRPolygon* MultitheadIntersection::createBB(OGREnvelope bb)
{
	OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
	OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);
	linearRing->addPoint(bb.MinX, bb.MinY);
	linearRing->addPoint(bb.MinX, bb.MaxY);
	linearRing->addPoint(bb.MaxX, bb.MaxY);
	linearRing->addPoint(bb.MaxX, bb.MinY);
	OGRPoint firstp;
	linearRing->getPoint(0, &firstp);
	linearRing->addPoint(&firstp);
	poPolygon->addRing(linearRing);//also crashed
	return poPolygon;
}
void MultitheadIntersection::setupGrid(OGREnvelope bb, double cellsize)
{
	if (grid)
		delete grid;
	grid = new Grid(bb, cellsize);
}

void MultitheadIntersection::setupGrid(std::string tiff)
{
	if (grid)
		delete grid;
	grid = new Grid;
	grid->fromFishnetRaster(tiff);
}

void MultitheadIntersection::intersect(std::string inshapefile, std::string outshapefile)
{
	int numcells = grid->ncols * grid->nrows;
	double gridcellsize = grid->_adfGeoTransform[1];
	GeomInfo* fishnet = new GeomInfo[numcells];
	int id = 0;
	for (int i = 0; i < grid->nrows; i++)
	{
		OGREnvelope tmpBound;
		tmpBound.MaxY = grid->bound.MaxY - gridcellsize * i;
		tmpBound.MinY = grid->bound.MaxY - gridcellsize * (i + 1);
		for (size_t j = 0; j < grid->ncols; j++)
		{
			OGREnvelope bb = tmpBound;
			bb.MinX = grid->bound.MinX + gridcellsize * j;
			bb.MaxX = grid->bound.MinX + gridcellsize * (j + 1);
			OGRPolygon* cell = createBB(bb);
			OGREnvelope geombb;
			cell->getEnvelope(&geombb);

			GeomInfo geominfo(cell, id, geombb);
			fishnet[id] = geominfo;
			id++;
		}

	}
	std::string idname = "FID_";
	std::string filename;
	bool found = false;
	for (int i = inshapefile.length() - 1; i >= 0; i--)
	{
		if (inshapefile[i] == '/' || inshapefile[i] == '\\')
		{
			filename = inshapefile.substr(i + 1, inshapefile.length() - (i + 1));
			filename = filename.substr(0, filename.length() - 4);
			if (filename.size() > 6)
				filename = filename.substr(0, 6);
			found = true;
			break;
		}
	}
	if (!found)
	{
		filename = inshapefile.substr(0, inshapefile.length() - 4);
		if (filename.size() > 6)
			filename = filename.substr(0, 6);
		found = true;
	}
	idname = idname + filename;
	GEOMGROUPS geoms;
	ShapeFile inshp(inshapefile.data());
	ShapeFile outshp;
	outshp.create(outshapefile, inshp.poLayer->GetSpatialRef(), 0, inshp.poLayer->GetGeomType());
	int idfield = outshp.getOrCreateField(idname.data(), OGRFieldType::OFTInteger);
	int fishfield = outshp.getOrCreateField("Id", OGRFieldType::OFTInteger);
	inshp.poLayer->ResetReading();
	OGRFeature *poFeature;
	id = 0;

	while ((poFeature = inshp.poLayer->GetNextFeature()) != NULL)
	{
		std::vector<OGRGeometry*> childrenGeoms;
		if (!dynamic_cast<OGRGeometryCollection*>(poFeature->GetGeometryRef()))
		{
			childrenGeoms.push_back(poFeature->GetGeometryRef());
		}
		else
		{
			OGRGeometryCollection* geomcollect = dynamic_cast<OGRGeometryCollection*>(poFeature->GetGeometryRef());
			for (int ichild = 0; ichild < geomcollect->getNumGeometries(); ichild++)
			{
				childrenGeoms.push_back(geomcollect->getGeometryRef(ichild));
			}
		}
		//if (geom->getGeometryType() == wkbMultiPolygon25D || geom->getGeometryType() == wkbMultiPolygon ||
		//	geom->getGeometryType() == wkbMultiLineString25D || geom->getGeometryType() == wkbMultiLineString)
		//{
		//	
		//}
		for (int igeom = 0; igeom <childrenGeoms.size(); igeom++)
		{
			OGREnvelope geombb;
			childrenGeoms[igeom]->getEnvelope(&geombb);
			GeomInfo geominfo(childrenGeoms[igeom], id, geombb);
			geoms.push_back(geominfo);

			int startCol = (int)((geombb.MinX - grid->bound.MinX) / gridcellsize);
			int startRow = (int)((grid->bound.MaxY - geombb.MaxY) / gridcellsize);
			int endCol = (int)((geombb.MaxX - grid->bound.MinX) / gridcellsize);
			int endRow = (int)((grid->bound.MaxY - geombb.MinY) / gridcellsize);

			for (size_t row = startRow; row <= endRow; row++)
			{
				for (size_t col = startCol; col <= endCol; col++)
				{
					int cellid = col + row * grid->ncols;
					OGRGeometry* p2 = fishnet[cellid].geom;
					OGRGeometry* pIntersection = childrenGeoms[igeom]->Intersection(p2);
					OGRFeature *poFeatureNew = OGRFeature::CreateFeature(outshp.poLayer->GetLayerDefn());
					poFeatureNew->SetGeometry(pIntersection);
					poFeatureNew->SetField(fishfield, cellid);
					poFeatureNew->SetField(idfield, id);
					outshp.poLayer->CreateFeature(poFeatureNew);
					OGRFeature::DestroyFeature(poFeatureNew);
				}

			}
		}
		if (id % 100 == 0)
			printf("%d\n", id);
		id++;
		OGRFeature::DestroyFeature(poFeature);
	}
	for (size_t i = 0; i < numcells; i++)
	{
		OGRGeometryFactory::destroyGeometry(fishnet[i].geom);
	}
	delete[] fishnet;

	inshp.close();
	outshp.close();
}

int _numcells;
int _ncols;
int _nrows;
OGREnvelope _bound;
double _gridcellsize;
GEOMGROUPS _geoms;
GeomInfo* _fishnet;
GEOMGROUPS* _results;
#include <thread> 
#include <mutex>  
void intersectFunc(int taskID,int startGeom,int numGeoms)
{
	for (int igeom = startGeom; igeom < startGeom + numGeoms; igeom++)
	{
		GeomInfo& geominfo = _geoms[igeom];
		int startCol = (int)((geominfo.bb.MinX - _bound.MinX) / _gridcellsize);
		int startRow = (int)((_bound.MaxY - geominfo.bb.MaxY) / _gridcellsize);
		int endCol = (int)((geominfo.bb.MaxX - _bound.MinX) / _gridcellsize);
		int endRow = (int)((_bound.MaxY - geominfo.bb.MinY) / _gridcellsize);
		GEOMGROUPS& result = _results[taskID];

		for (size_t row = startRow; row <= endRow; row++)
		{
			for (size_t col = startCol; col <= endCol; col++)
			{
				int cellid = col + row * _ncols;
				OGRGeometry* p2 = _fishnet[cellid].geom;
				OGRGeometry* pIntersection = geominfo.geom->Intersection(p2);
				if (pIntersection && pIntersection->IsValid())
				{
					OGREnvelope geombb;
					pIntersection->getEnvelope(&geombb);
					result.push_back(GeomInfo(pIntersection, geominfo.id, cellid, geombb));
				}
			}

		}
		//printf("%d\n", igeom);
	}
}

void MultitheadIntersection::intersectParallel(std::string inshapefile, std::string outshapefile, int numThreads)
{
	_results = new GEOMGROUPS[numThreads];
	_ncols = grid->ncols;
	_nrows = grid->nrows;
	_numcells = _ncols * _nrows;
	_bound = grid->bound;
	_gridcellsize = grid->_adfGeoTransform[1];
	_fishnet = new GeomInfo[_numcells];

	int id = 0;
	for (int i = 0; i < grid->nrows; i++)
	{
		OGREnvelope tmpBound;
		tmpBound.MaxY = _bound.MaxY - _gridcellsize * i;
		tmpBound.MinY = _bound.MaxY - _gridcellsize * (i + 1);
		for (size_t j = 0; j < grid->ncols; j++)
		{
			OGREnvelope bb = tmpBound;
			bb.MinX = _bound.MinX + _gridcellsize * j;
			bb.MaxX = _bound.MinX + _gridcellsize * (j + 1);
			OGRPolygon* cell = createBB(bb);
			OGREnvelope geombb;
			cell->getEnvelope(&geombb);

			GeomInfo geominfo(cell, id, geombb);
			_fishnet[id] = geominfo;
			id++;
		}

	}
	std::string idname = "FID_";
	std::string filename;
	bool found = false;
	for (int i = inshapefile.length() - 1; i >= 0; i--)
	{
		if (inshapefile[i] == '/' || inshapefile[i] == '\\')
		{
			filename = inshapefile.substr(i + 1, inshapefile.length() - (i + 1));
			filename = filename.substr(0, filename.length() - 4);
			if (filename.size() > 6)
				filename = filename.substr(0, 6);
			found = true;
			break;
		}
	}
	if (!found)
	{
		filename = inshapefile.substr(0, inshapefile.length() - 4);
		if (filename.size() > 6)
			filename = filename.substr(0, 6);
		found = true;
	}
	idname = idname + filename;

	ShapeFile inshp(inshapefile.data());
	ShapeFile outshp;
	outshp.create(outshapefile, inshp.poLayer->GetSpatialRef(), 0, inshp.poLayer->GetGeomType());
	int idfield = outshp.getOrCreateField(idname.data(), OGRFieldType::OFTInteger);
	int fishfield = outshp.getOrCreateField("Id", OGRFieldType::OFTInteger);
	inshp.poLayer->ResetReading();
	OGRFeature *poFeature;
	id = 0;

	while ((poFeature = inshp.poLayer->GetNextFeature()) != NULL)
	{
		std::vector<OGRGeometry*> childrenGeoms;
		if (!dynamic_cast<OGRGeometryCollection*>(poFeature->GetGeometryRef()))
		{
			childrenGeoms.push_back(poFeature->GetGeometryRef()->clone());
		}
		else
		{
			OGRGeometryCollection* geomcollect = dynamic_cast<OGRGeometryCollection*>(poFeature->GetGeometryRef());
			for (int ichild = 0; ichild < geomcollect->getNumGeometries(); ichild++)
			{
				childrenGeoms.push_back(geomcollect->getGeometryRef(ichild)->clone());
			}
		}
		for (int igeom = 0; igeom <childrenGeoms.size(); igeom++)
		{
			OGREnvelope geombb;
			childrenGeoms[igeom]->getEnvelope(&geombb);
			GeomInfo geominfo(childrenGeoms[igeom], id, geombb);
			_geoms.push_back(geominfo);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	inshp.close();

	std::thread tasks[1000];
	int numgeoms_per_task = _geoms.size() / numThreads;
	int startGeoms = 0;
	for (int itask = 0; itask < numThreads; itask++)
	{
		int numgeoms = numgeoms_per_task;
		if (startGeoms + numgeoms > _geoms.size())
			numgeoms = _geoms.size() - startGeoms;
		tasks[itask] =  std::thread(intersectFunc, itask, startGeoms, numgeoms);
		startGeoms += numgeoms;
	}
	clock_t t;
	int f;
	t = clock();
	for (int itask = 0; itask < numThreads; itask++)
	{
		tasks[itask].join();
	}
	t = clock() - t;
	printf("It took me %d clicks (%f seconds).\n", t, ((float)t) / CLOCKS_PER_SEC);
	for (int itask = 0; itask < numThreads; itask++)
	{
		GEOMGROUPS result = _results[itask];
		for (int igeom = 0; igeom < result.size(); igeom++)
		{
			if (result[igeom].geom && result[igeom].geom->IsValid())
			{
				OGRFeature *poFeatureNew = OGRFeature::CreateFeature(outshp.poLayer->GetLayerDefn());
				poFeatureNew->SetGeometry(result[igeom].geom);
				poFeatureNew->SetField(fishfield, result[igeom].fishnetid);
				poFeatureNew->SetField(idfield, result[igeom].id);
				outshp.poLayer->CreateFeature(poFeatureNew);
				OGRFeature::DestroyFeature(poFeatureNew);
			}

		}

	}

	outshp.close();

	for (size_t i = 0; i < _numcells; i++)
	{
		OGRGeometryFactory::destroyGeometry(_fishnet[i].geom);
	}

	for (size_t i = 0; i < _geoms.size(); i++)
	{
		OGRGeometryFactory::destroyGeometry(_geoms[i].geom);
	}

	delete[] _fishnet;


}
