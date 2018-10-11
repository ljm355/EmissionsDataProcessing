#pragma once
#include "Qt/qimage.h"
#include "Qt/qpainter.h"
#include "ogrsf_frmts.h"
#include <QFileinfo>
#include "GDAL_DS.h"
#include "ShapeFile.h"
class RasterMask
{

public:
	QImage* _image;
	double _adfGeoTransform[6];
	int _xsize;
	int _ysize;
	uchar _nodata;
	int _stride;
	RasterMask(QImage* img, double xmin, double ymin, double xmax, double ymax)
	{
		_image = img;
		_adfGeoTransform[0] = xmin;
		_adfGeoTransform[1] = (xmax - xmin) / (img->width());
		_adfGeoTransform[2] = 0;
		_adfGeoTransform[3] = ymax;
		_adfGeoTransform[4] = 0;
		_adfGeoTransform[5] = -(ymax - ymin) / (img->height());
		_nodata = 0;
		_xsize = img->width();
		_ysize = img->height();
		_stride = 4;
	}
	~RasterMask()
	{
		delete _image;
	}

	//adfGeoTransform[0] /* top left x */
	//adfGeoTransform[1] /* w-e pixel resolution */
	//adfGeoTransform[2] /* 0 */
	//adfGeoTransform[3] /* top left y */
	//adfGeoTransform[4] /* 0 */
	//adfGeoTransform[5] /* n-s pixel resolution (negative value) */
};

class PolygonFiller
{
public :
	static void fillPolygon(QPainter * painter, const QBrush & brush, const int & xsize, const int & ysize, const OGRPolygon * polygon);
	static RasterMask * fillPolygon(const QBrush & br, const float & cellsize, OGRGeometry * poGeometry);
};

//class PolygonGeom :public OGREnvelope
//{
//public:
//
//	PolygonGeom* children[4];
//	OGREnvelope bound;
//	std::vector<OGRGeometry*> geometryList;
//	PolygonGeom()
//	{
//		for (size_t i = 0; i < 4; i++)
//		{
//			children[i] = NULL;
//		}
//	}
//
//	void insertChild(OGRGeometry* geom)
//	{
//
//	}
//
//};



class PolygonGeom :public OGREnvelope
{
public:
	std::vector<QPolygonF> polygons;
	std::vector<PolygonGeom*> children;
	PolygonGeom(OGRGeometry* geom)
	{
		OGREnvelope bound;
		geom->getEnvelope(&bound);
		this->MinX = bound.MinX;
		this->MaxX = bound.MaxX;
		this->MinY = bound.MinY;
		this->MaxY = bound.MaxY;
		std::vector<OGRLinearRing*> toRings;
		OGRwkbGeometryType totype = geom->getGeometryType();
		if (totype == wkbPolygon || totype == wkbPolygon25D)
		{
			toRings.push_back(((OGRPolygon*)geom)->getExteriorRing());
		}
		else if (totype == wkbMultiPolygon || totype == wkbMultiPolygon25D)
		{
			OGRMultiPolygon* multiPoly = (OGRMultiPolygon*)geom;
			for (size_t ipolygon = 0; ipolygon < multiPoly->getNumGeometries(); ipolygon++)
			{
				toRings.push_back(((OGRPolygon*)multiPoly->getGeometryRef(ipolygon))->getExteriorRing());
			}
		}
		for (size_t itopolygon = 0; itopolygon < toRings.size(); itopolygon++)
		{
			OGRLinearRing* toouterring = toRings[itopolygon];
			QPolygonF qpolygon;
			for (size_t itopoint = 0; itopoint < toouterring->getNumPoints(); itopoint++)
			{
				OGRPoint topt;
				toouterring->getPoint(itopoint, &topt);
				qpolygon.push_back(QPointF(topt.getX(), topt.getY()));
			}
			qpolygon.push_back(qpolygon.first());
			polygons.push_back(qpolygon);
		}


	}
	~PolygonGeom()
	{
	}
};

template <class T>
class ExtractByPolygon
{
public:
	GDAL_DS<T>* g_pDS;
	T* g_pData;
	ExtractByPolygon(GDAL_DS<T>* ds)
	{
		g_pDS = ds;
		g_pData = ds->readData(1);
	}

	ExtractByPolygon()
	{
		g_pDS = NULL;
		g_pData = NULL;
	}
	double sumVolumn(RasterMask* mask)
	{
	
		size_t   xsizeVolumn = g_pDS->ncols;
		size_t   ysizeVolumn = g_pDS->nrows;
		size_t   ncellsVolumn = xsizeVolumn*ysizeVolumn;
		T*   dataVolumn = g_pData;
	
		double* geoTransformVolumn = g_pDS->adfGeoTransform;
		T nodataVolumn = g_pDS->getNoData(1);

		size_t   xsizeMask = mask->_xsize;
		size_t   ysizeMask = mask->_ysize;
		size_t   ncellsMask = xsizeMask*ysizeMask;
		const uchar *dataMask = mask->_image->bits();
		uchar nodataMask = mask->_nodata;
		size_t stride = mask->_stride;
		double* geoTransformMask = mask->_adfGeoTransform;
	
		double sum = 0;
		size_t validCount = 0;
		const uchar* pdataMask = dataMask;
		/*float* outdata = NULL;
		if (outputfile != "")
		{
			outdata = new float[xsizeMask*ysizeMask];
			memset(outdata, 0, xsizeMask*ysizeMask*sizeof(float));
		}*/
		size_t curIdx = 0;
	
		for (size_t yindexMask = 0; yindexMask < ysizeMask; yindexMask++)
		{
			float yMask = geoTransformMask[3] + geoTransformMask[5] * 0.5 + geoTransformMask[5] * yindexMask;
			size_t yindexVolumn = (int)((yMask - geoTransformVolumn[3]) / geoTransformVolumn[5]);
			float yVolumn = geoTransformVolumn[3] + geoTransformVolumn[5] * 0.5 + geoTransformVolumn[5] * yindexVolumn;
			for (size_t xindexMask = 0; xindexMask < xsizeMask; xindexMask++)
			{
				//float xMask = geoTransformVolumn[0] + geoTransformVolumn[1] * 0.5 + geoTransformVolumn[1] * xindexMask
				curIdx++;
				if (*pdataMask == nodataMask)
				{
					pdataMask += stride;
					continue;
				}
	
				float xMask = geoTransformMask[0] + geoTransformMask[1] * 0.5 + geoTransformMask[1] * xindexMask;
				size_t xindexVolumn = (int)((xMask - geoTransformVolumn[0]) / geoTransformVolumn[1]);
				float xVolumn = geoTransformVolumn[0] + geoTransformVolumn[1] * 0.5 + geoTransformVolumn[1] * xindexVolumn;
				int idx = yindexVolumn*xsizeVolumn + xindexVolumn;
	
				double valVolumn = dataVolumn[idx];
				if (valVolumn == nodataVolumn)
				{
					pdataMask += stride;
					continue;
				}
				sum += valVolumn;
				validCount++;
				pdataMask += stride;
	/*			if (outdata)
					outdata[curIdx - 1] = result;*/
			}
		}
		//if (outdata)
		//{
		//	char **papszOptions = NULL;
		//	const char *pszFormat = "GTiff";
		//	GDALDriver *poDriver;
		//	char **papszMetadata;
		//	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
		//	GDALDataset *destDataset = poDriver->Create(outputfile.data(), xsizeMask, ysizeMask, 1, GDT_Float32, papszOptions);
		//	GDALRasterBand *poBandDest = destDataset->GetRasterBand(1);
		//	poBandDest->RasterIO(GF_Write, 0, 0, xsizeMask, ysizeMask, outdata, xsizeMask, ysizeMask, GDT_Float32, 0, 0);
		//	destDataset->SetGeoTransform(mask->_adfGeoTransform);
		//	poBandDest->SetNoDataValue(0);
		//	destDataset->FlushCache();
		//	GDALClose((GDALDatasetH)destDataset);
	
		//	delete[] outdata;
		//}

		return sum;
	}
	void extract(std::string polyonshapefile,std::string outAttibute)
	{
		ShapeFile shp(polyonshapefile,1);
		int ioutfield = shp.getOrCreateField(outAttibute.data(), OGRFieldType::OFTReal);
		OGRFeature *poFeature;
		double sum = 0;
		shp.poLayer->ResetReading();
		int idx = -1;
		QBrush br(Qt::white);
		double cellsize = g_pDS->adfGeoTransform[1] / 10;
		while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
		{
			idx++;
			OGRGeometry* poGeometry = poFeature->GetGeometryRef();
			RasterMask* mask = PolygonFiller::fillPolygon(br, cellsize, poGeometry);
			double sum = sumVolumn(mask);
			poFeature->SetField(ioutfield, sum);
			shp.poLayer->SetFeature(poFeature);
			OGRFeature::DestroyFeature(poFeature);
		}
		shp.close();

	}
	//static bool pointInPolygon(OGRPoint* pt, OGRLinearRing* polygon)
	//{

	//	float x[8];
	//	float y[8];
	//	float x1, x2;

	//	/* The coordinates of the point */
	//	float px, py;

	//	/* How many times the ray crosses a line-segment */
	//	int crossings = 0;

	//	/* Coordinates of the points */
	//	x[0] = 100;	y[0] = 100;
	//	x[1] = 200;	y[1] = 200;
	//	x[2] = 300;	y[2] = 200;
	//	x[3] = 300;	y[3] = 170;
	//	x[4] = 240;	y[4] = 170;
	//	x[5] = 240;	y[5] = 90;
	//	x[6] = 330;	y[6] = 140;
	//	x[7] = 270;	y[7] = 30;

	//	/* Iterate through each line */
	//	for (int i = 0; i < 8; i++) {

	//		/* This is done to ensure that we get the same result when
	//		the line goes from left to right and right to left */
	//		if (x[i] < x[(i + 1) % 8]) {
	//			x1 = x[i];
	//			x2 = x[(i + 1) % 8];
	//		}
	//		else {
	//			x1 = x[(i + 1) % 8];
	//			x2 = x[i];
	//		}

	//		/* First check if the ray is possible to cross the line */
	//		if (px > x1 && px <= x2 && (py < y[i] || py <= y[(i + 1) % 8])) {
	//			static const float eps = 0.000001;

	//			/* Calculate the equation of the line */
	//			float dx = x[(i + 1) % 8] - x[i];
	//			float dy = y[(i + 1) % 8] - y[i];
	//			float k;

	//			if (fabs(dx) < eps) {
	//				k = INFINITY;	// math.h
	//			}
	//			else {
	//				k = dy / dx;
	//			}

	//			float m = y[i] - k * x[i];

	//			/* Find if the ray crosses the line */
	//			float y2 = k * px + m;
	//			if (py <= y2) {
	//				crossings++;
	//			}
	//		}
	//	}

	//	printf("The point is crossing %d lines", crossings);
	//	if (crossings % 2 == 1) {
	//		printf(" thus it is inside the polygon");
	//	}
	//	printf("\n");
	//}

	static void link(std::string fromshapefile, std::string toshapefile, std::string linkid, std::vector<std::string> outAttibutes)
	{
		ShapeFile fromshp(fromshapefile, 1);
		ShapeFile toshp(toshapefile, 0);
		std::vector<int> tofields;
		std::vector<int> fromfields;
		for (size_t i = 0; i < outAttibutes.size(); i++)
		{
			 int ifield = toshp.getField(outAttibutes[i].data());
			 tofields.push_back(ifield);
			 if (fromshp.poLayer->GetLayerDefn()->GetFieldIndex(outAttibutes[i].data()) < 0)
			 {
			    fromshp.poLayer->CreateField(toshp.poLayer->GetLayerDefn()->GetFieldDefn(ifield));
		
			 }
			 fromfields.push_back(fromshp.poLayer->GetLayerDefn()->GetFieldIndex(outAttibutes[i].data()));
		}

		
		//int idfield = fromshp.getOrCreateField("toID", OGRFieldType::OFTInteger);
		int idfield = fromshp.getOrCreateField(linkid.data(), OGRFieldType::OFTInteger);
		OGRFeature *fromFeature;
		OGRFeature *toFeature;

		int idx = -1;		
		//OGRwkbGeometryType fromtype = fromshp.poLayer->GetGeomType();
		//memcpy(buf, ds->projection.data(), ds->projection.size());
		//oTargetSRS.importFromWkt(&buf);
		//OGRSpatialReference oSourceSRS;
		//oSourceSRS.SetWellKnownGeogCS("WGS84");
		toshp.poLayer->ResetReading();
		std::vector<PolygonGeom> polygonGeoms;
		while ((toFeature = toshp.poLayer->GetNextFeature()) != NULL)
		{
			polygonGeoms.push_back(toFeature->GetGeometryRef());
			OGRFeature::DestroyFeature(toFeature);
		}

		//OGREnvelope fromShpBound;
		//fromshp.poLayer->GetExtent(&fromShpBound);

		//OGREnvelope toShpBound;
		//toshp.poLayer->GetExtent(&toShpBound);

		//double gridcellsize1 = (fromShpBound.MaxX - fromShpBound.MinX) *  (fromShpBound.MaxY - fromShpBound.MinY) / fromshp.poLayer->GetFeatureCount();
		//double gridcellsize2 = (toShpBound.MaxX - toShpBound.MinX) *  (toShpBound.MaxY - toShpBound.MinY) / toShpBound.poLayer->GetFeatureCount();
		//double gridcellsize = 

		fromshp.poLayer->ResetReading();
		OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(fromshp.poLayer->GetSpatialRef(), toshp.poLayer->GetSpatialRef());
		while ((fromFeature = fromshp.poLayer->GetNextFeature()) != NULL)
		{
			idx++;
			OGRGeometry* fromGeometry = fromFeature->GetGeometryRef();

			//toshp.poLayer->SetSpatialFilter(fromGeometry);
			OGREnvelope fromRect;
			fromGeometry->getEnvelope(&fromRect);
			OGRPoint centroid;
			fromGeometry->Centroid(&centroid);

			OGRwkbGeometryType fromtype = fromGeometry->getGeometryType();
			//poCT->Transform(1, &fromRect.MinX, &fromRect.MinY);
			//poCT->Transform(1, &fromRect.MaxX, &fromRect.MaxY);

			toshp.poLayer->SetSpatialFilterRect(fromRect.MinX, fromRect.MinY, fromRect.MaxX, fromRect.MaxY);
			toshp.poLayer->ResetReading();
		
		/*	std::vector<OGRLinearRing*> fromRings;
			if (fromtype == wkbPolygon || fromtype == wkbPolygon25D)
			{
				fromRings.push_back(((OGRPolygon*)fromGeometry)->getExteriorRing());
			}
			else if (fromtype == wkbMultiPolygon || fromtype == wkbMultiPolygon25D)
			{
				OGRMultiPolygon* multiPoly = (OGRMultiPolygon*)fromGeometry;
				for (size_t ipolygon = 0; ipolygon < multiPoly->getNumGeometries(); ipolygon++)
				{
					fromRings.push_back(((OGRPolygon*)multiPoly->getGeometryRef(ipolygon))->getExteriorRing());
				}
			}*/

			//while ((toFeature = toshp.poLayer->GetNextFeature()) != NULL)
			//{
			//	OGRGeometry* toGeometry = toFeature->GetGeometryRef();
			//	OGRwkbGeometryType totype = toGeometry->getGeometryType();
			//	std::vector<OGRLinearRing*> toRings;
			//	if (totype == wkbPolygon || totype == wkbPolygon25D)
			//	{
			//		toRings.push_back(((OGRPolygon*)toGeometry)->getExteriorRing());
			//	}
			//	else if (totype == wkbMultiPolygon || totype == wkbMultiPolygon25D)
			//	{
			//		OGRMultiPolygon* multiPoly = (OGRMultiPolygon*)toGeometry;
			//		for (size_t ipolygon = 0; ipolygon < multiPoly->getNumGeometries(); ipolygon++)
			//		{
			//			toRings.push_back(((OGRPolygon*)multiPoly->getGeometryRef(ipolygon))->getExteriorRing());
			//		}
			//	
			//	}
			//	int numpoints = 0;
			//	int numpointsinpolygon = 0;
			//	for (size_t itopolygon = 0; itopolygon < toRings.size(); itopolygon++)
			//	{
			//		OGRLinearRing* toouterring = toRings[itopolygon];
			//		QPolygonF qpolygon;
			//		for (size_t itopoint = 0; itopoint < toouterring->getNumPoints(); itopoint++)
			//		{
			//			OGRPoint topt;
			//			toouterring->getPoint(itopoint, &topt);
			//			qpolygon.push_back(QPointF(topt.getX(), topt.getY()));
			//		}
			//		qpolygon.push_back(qpolygon.first());
			//		//printf("%d\n", toFeature->GetFID());


			//		
			//		QPointF qpt(centroid.getX(), centroid.getY());
			//		numpoints++;
			//		//if (qpolygon.containsPoint(qpt, Qt::OddEvenFill))
			//		//{
			//		//	toID = toFeature->GetFID();
			//		//	numpointsinpolygon++;
			//		//	break;
			//		//}
		
			//		//for (int ifrompolygon = 0; ifrompolygon < fromRings.size(); ifrompolygon++)
			//		//{
			//		//	OGRLinearRing* fromouterring = fromRings[ifrompolygon];
			//		//	int factor = 1;
			//		//	if (fromouterring->getNumPoints() > 10)
			//		//		factor = fromouterring->getNumPoints() / 10;
			//		//	for (int ifrompoint = 0; ifrompoint < fromouterring->getNumPoints(); ifrompoint++)
			//		//	{
			//		//		if (ifrompoint % factor == 0)
			//		//		{
			//		//			OGRPoint frompt;
			//		//			fromouterring->getPoint(ifrompoint, &frompt);
			//		//			double ptx = frompt.getX();
			//		//			double pty = frompt.getY();
			//		//			//poCT->Transform(1, &ptx, &pty);
			//		//			QPointF qpt(ptx, pty);
			//		//			numpoints++;
			//		//			if (qpolygon.containsPoint(qpt, Qt::OddEvenFill))
			//		//				numpointsinpolygon++;
			//		//		}

			//		//	}
			//		//}
			//		
			//	}
			//	OGRFeature::DestroyFeature(toFeature);
			//	if (toID > -1)
			//		break;
			//}
			int toID = -1;
			QPointF qpt(centroid.getX(), centroid.getY());
			for (size_t i = 0; i < polygonGeoms.size(); i++)
			{
				PolygonGeom& polygeom = polygonGeoms[i];
				if (!polygeom.Intersects(fromRect))
					continue;
				for (size_t j = 0; j < polygeom.polygons.size(); j++)
				{
					if (polygeom.polygons[j].containsPoint(qpt, Qt::OddEvenFill))
					{
						toID = i;
						break;
					}
				}
			}


			fromFeature->SetField(idfield, toID);
			//if (toID > -1)
			//{
			//	toFeature = toshp.poLayer->GetFeature(toID);
			//	for (size_t iattribute = 0; iattribute < tofields.size(); iattribute++)
			//	{
			//		fromFeature->SetField(fromfields[iattribute], toFeature->GetRawFieldRef(tofields[iattribute]));
			//	}
			//	//OGRFeature::DestroyFeature(toFeature);
			//}

			fromshp.poLayer->SetFeature(fromFeature);
			OGRFeature::DestroyFeature(fromFeature);
			if (idx % 100 == 0)
				printf("%d\n", idx);
		}
		fromshp.close();

	}

	static void link2(std::string fromshapefile, std::string toshapefile, std::string linkid, std::vector<std::string> outAttibutes)
	{
		ShapeFile fromshp(fromshapefile, 0);
		ShapeFile toshp(toshapefile, 1);
		std::vector<int> tofields;
		std::vector<int> fromfields;
		for (size_t i = 0; i < outAttibutes.size(); i++)
		{
			int ifield = fromshp.getField(outAttibutes[i].data());
			tofields.push_back(ifield);
			if (toshp.poLayer->GetLayerDefn()->GetFieldIndex(outAttibutes[i].data()) < 0)
			{
				toshp.poLayer->CreateField(fromshp.poLayer->GetLayerDefn()->GetFieldDefn(ifield));

			}
			fromfields.push_back(toshp.poLayer->GetLayerDefn()->GetFieldIndex(outAttibutes[i].data()));
		}


		//int idfield = fromshp.getOrCreateField("toID", OGRFieldType::OFTInteger);
		int idfield = toshp.getOrCreateField(linkid.data(), OGRFieldType::OFTInteger);
		OGRFeature *fromFeature;
		fromshp.poLayer->ResetReading();
		int idx = -1;
		//OGRwkbGeometryType fromtype = fromshp.poLayer->GetGeomType();
		//memcpy(buf, ds->projection.data(), ds->projection.size());
		//oTargetSRS.importFromWkt(&buf);
		//OGRSpatialReference oSourceSRS;
		//oSourceSRS.SetWellKnownGeogCS("WGS84");
		OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(toshp.poLayer->GetSpatialRef(), fromshp.poLayer->GetSpatialRef());

		while ((fromFeature = fromshp.poLayer->GetNextFeature()) != NULL)
		{
			idx++;
			OGRGeometry* fromGeometry = fromFeature->GetGeometryRef();
			//toshp.poLayer->SetSpatialFilter(fromGeometry);
			OGREnvelope fromRect;
			fromGeometry->getEnvelope(&fromRect);
			OGRwkbGeometryType fromtype = fromGeometry->getGeometryType();
			//poCT->Transform(1, &fromRect.MinX, &fromRect.MinY);
			//poCT->Transform(1, &fromRect.MaxX, &fromRect.MaxY);

			toshp.poLayer->SetSpatialFilterRect(fromRect.MinX, fromRect.MinY, fromRect.MaxX, fromRect.MaxY);
			toshp.poLayer->ResetReading();
			int toID = -1;
			std::vector<OGRLinearRing*> fromRings;
			if (fromtype == wkbPolygon || fromtype == wkbPolygon25D)
			{
				fromRings.push_back(((OGRPolygon*)fromGeometry)->getExteriorRing());
			}
			else if (fromtype == wkbMultiPolygon || fromtype == wkbMultiPolygon25D)
			{
				OGRMultiPolygon* multiPoly = (OGRMultiPolygon*)fromGeometry;
				for (size_t ipolygon = 0; ipolygon < multiPoly->getNumGeometries(); ipolygon++)
				{
					fromRings.push_back(((OGRPolygon*)multiPoly->getGeometryRef(ipolygon))->getExteriorRing());
				}
			}
			OGRFeature *toFeature;
			while ((toFeature = toshp.poLayer->GetNextFeature()) != NULL)
			{
				OGRGeometry* toGeometry = toFeature->GetGeometryRef();
				OGRwkbGeometryType totype = toGeometry->getGeometryType();
				std::vector<OGRLinearRing*> toRings;
				if (totype == wkbPolygon || totype == wkbPolygon25D)
				{
					toRings.push_back(((OGRPolygon*)toGeometry)->getExteriorRing());
				}
				else if (totype == wkbMultiPolygon || totype == wkbMultiPolygon25D)
				{
					OGRMultiPolygon* multiPoly = (OGRMultiPolygon*)toGeometry;
					for (size_t ipolygon = 0; ipolygon < multiPoly->getNumGeometries(); ipolygon++)
					{
						toRings.push_back(((OGRPolygon*)multiPoly->getGeometryRef(ipolygon))->getExteriorRing());
					}

				}
				int numpoints = 0;
				int numpointsinpolygon = 0;
				for (int ifrompolygon = 0; ifrompolygon < fromRings.size(); ifrompolygon++)
				{
					OGRLinearRing* fromouterring = fromRings[ifrompolygon];
					QPolygonF qpolygon;
					for (int ifrompoint = 0; ifrompoint < fromouterring->getNumPoints(); ifrompoint++)
					{
						OGRPoint frompt;
						fromouterring->getPoint(ifrompoint, &frompt);
						double ptx = frompt.getX();
						double pty = frompt.getY();
						//poCT->Transform(1, &ptx, &pty);
						qpolygon.push_back(QPointF(ptx, pty));
					}
					qpolygon.push_back(qpolygon.first());
					for (size_t itopolygon = 0; itopolygon < toRings.size(); itopolygon++)
					{
						OGRLinearRing* toouterring = toRings[itopolygon];
						int factor = 1;
						if (toouterring->getNumPoints() > 10)
						{
							factor = toouterring->getNumPoints() / 10;
							//factor = toouterring->getNumPoints() / factor;
						}

						for (size_t itopoint = 0; itopoint < toouterring->getNumPoints(); itopoint++)
						{
							if (itopoint % factor == 0)
							{
								numpoints++;
								OGRPoint topt;
								toouterring->getPoint(itopoint, &topt);
								QPointF qpt(topt.getX(), topt.getY());
								if (qpolygon.containsPoint(qpt, Qt::OddEvenFill))
									numpointsinpolygon++;
							}
			
						}
						if ((double)numpointsinpolygon / (double)numpoints > 0.5)
						{
							toID = fromFeature->GetFID();
							toFeature->SetField(idfield, toID);
							break;
						}
					}
				}


				toshp.poLayer->SetFeature(toFeature);
				OGRFeature::DestroyFeature(toFeature);
			}

			//if (toID > -1)
			//{
			//	toFeature = toshp.poLayer->GetFeature(toID);
			//	for (size_t iattribute = 0; iattribute < tofields.size(); iattribute++)
			//	{
			//		fromFeature->SetField(fromfields[iattribute], toFeature->GetRawFieldRef(tofields[iattribute]));
			//	}
			//	OGRFeature::DestroyFeature(toFeature);
			//}
			OGRFeature::DestroyFeature(fromFeature);
			//if (idx % 100 == 0)
			printf("%d\n", idx);
		}
		toshp.close();
		fromshp.close();
	}
	~ExtractByPolygon() {};
};

