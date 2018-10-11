#pragma once
#include "gdal_priv.h"
#include "gdal_rat.h"
#include "ogrsf_frmts.h"

#include "qt4/QtGui/qapplication.h"
#include "qt4/QtGui/qpainter.h"
#include "qt4/QtGui/qimage.h"
#include <osg/Image>
#include <osgDB/WriteFile>
#include <osg/Vec3>
#include <osg/Vec2>
#include <osg/Vec2d>
#include <osg/Vec3d>
//adfGeoTransform[0] /* top left x */
//adfGeoTransform[1] /* w-e pixel resolution */
//adfGeoTransform[2] /* 0 */
//adfGeoTransform[3] /* top left y */
//adfGeoTransform[4] /* 0 */
//adfGeoTransform[5] /* n-s pixel resolution (negative value) */
class GDAL_DSInfo
{
public:
	OGREnvelope bound;
	size_t nrows;
	size_t ncols;
	size_t slice;
	int numbands;
	double adfGeoTransform[6];
	std::string filename;
	std::string pszFormat;
	std::string projection;
};
template <class T>
class GDAL_DS : public GDAL_DSInfo
{
public:
	GDAL_DS()
	{
		dataset = NULL;
		pszFormat = "GTiff";
	}
	~GDAL_DS() {
		close();
	}
	bool isInside(double x, double y)
	{
		if (x >= bound.MinX && x <= bound.MaxX && y >= bound.MinY && y <= bound.MaxY)
			return true;
		return false;
	}
	bool open(std::string fname, GDALAccess mode = GA_ReadOnly)
	{
		close();
		filename = fname;
		dataset = (GDALDataset*)GDALOpen(filename.data(), mode);
		if (!dataset)
			return false;
		dataset->GetGeoTransform(adfGeoTransform);
		ncols = dataset->GetRasterXSize();
		nrows = dataset->GetRasterYSize();
		slice = nrows*ncols;
		numbands = dataset->GetRasterCount();
		projection = dataset->GetProjectionRef();
		setGeoTransform(adfGeoTransform);
	}
	void create(std::string fname)
	{
		close();
		filename = fname;
		char **papszOptions = NULL;
		GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat.data());
		poDriver->Delete(filename.data());
		papszOptions = CSLSetNameValue(papszOptions, "TILED", "YES");
		papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");
		dataset = poDriver->Create(filename.data(), ncols, nrows, numbands, getType(), papszOptions);
		dataset->SetGeoTransform(adfGeoTransform);
		dataset->SetProjection(projection.data());
	}
	void setGeoTransform(double* transform)
	{
		memcpy(adfGeoTransform, transform, sizeof(double) * 6);
		bound.MinX = adfGeoTransform[0];
		bound.MaxY = adfGeoTransform[3];
		bound.MaxX = adfGeoTransform[0] + adfGeoTransform[1] * ncols;
		bound.MinY = adfGeoTransform[3] + adfGeoTransform[5] * nrows;
	}
	T* readData(int bandidx)
	{
		if (!dataset)
			return NULL;
		T* data = new T[slice];
		GDALDataType gdt = getType();
		dataset->GetRasterBand(bandidx)->RasterIO(GF_Read, 0, 0, ncols, nrows, data, ncols, nrows, gdt, 0, 0);
		return data;
	}
	void scaleTotal(int bandidx, double total)
	{
		T* data = readData(bandidx);
		T nodata = (T)getNoData(bandidx);
		T* pdata = data;
		double sum = 0;
		for (size_t i = 0; i < slice; i++)
		{
			if (*pdata != nodata)
			{
				sum += *pdata;
			}
			pdata++;
		}

		pdata = data;
		for (size_t i = 0; i < slice; i++)
		{
			if (*pdata != nodata)
			{
				double fraction = *pdata / sum;
				*pdata = (T)(total * fraction);
			}
			pdata++;
		}
		writeData(bandidx, data, nodata);
		delete[] data;
	}
	void readData(int bandidx, T* data)
	{
		if (!dataset)
			return;
		GDALDataType gdt = getType();
		dataset->GetRasterBand(bandidx)->RasterIO(GF_Read, 0, 0, ncols, nrows, data, ncols, nrows, gdt, 0, 0);
	}

	void writeData(int bandidx, T* data, double nodata)
	{
		GDALRasterBand *pBand = dataset->GetRasterBand(bandidx);
		pBand->RasterIO(GF_Write, 0, 0, ncols, nrows, data, ncols, nrows, getType(), 0, 0);
		pBand->SetNoDataValue(nodata);
	}
	void writeData(int bandidx, T* data)
	{
		GDALRasterBand *pBand = dataset->GetRasterBand(bandidx);
		pBand->RasterIO(GF_Write, 0, 0, ncols, nrows, data, ncols, nrows, getType(), 0, 0);
	}
	double sum(int bandindex)
	{
		T* data = readData(bandindex);
		T nodata = (T)getNoData(bandindex);
		T* pdata = data;
		double sum = 0;
		for (size_t i = 0; i < slice; i++)
		{
			if (*pdata != nodata)
			{
				sum += *pdata;
			}
			pdata++;
		}
		delete[] data;
		return sum;
	}
	void scale(int bandindex, T factor)
	{
		T* data = readData(bandindex);
		T nodata = (T)getNoData(bandindex);
		T* pdata = data;
		double sum = 0;
		for (size_t i = 0; i < slice; i++)
		{
			if (*pdata != nodata)
			{
				*pdata = *pdata * factor;
			}
			pdata++;
		}
		writeData(bandindex, data, nodata);
		delete[] data;
	}
	double getNoData(int bandidx)
	{
		return dataset->GetRasterBand(bandidx)->GetNoDataValue();
	}
	void setNoData(int bandidx,double nodata)
	{
		return dataset->GetRasterBand(bandidx)->SetNoDataValue(nodata);
	}
	void close(){
		if (dataset)
			GDALClose(dataset);
		dataset = NULL;
	}
	GDALDataType getType()
	{
		if (typeid(T) == typeid(unsigned char) || typeid(T) == typeid(char))
			return GDT_Byte;
		if (typeid(T) == typeid(unsigned short))
			return GDT_UInt16;
		if (typeid(T) == typeid(short))
			return GDT_Int16;
		if (typeid(T) == typeid(unsigned int))
			return GDT_UInt32;
		if (typeid(T) == typeid(int))
			return GDT_Int32;
		if (typeid(T) == typeid(float))
			return GDT_Float32;
		if (typeid(T) == typeid(double))
			return GDT_Float64;
		return GDT_Unknown;
	}
	void setDSInfo(GDAL_DSInfo* info)
	{
		nrows = info->nrows;
		ncols = info->ncols;
		slice = info->slice;
		numbands = info->numbands;
		filename = info->filename;
		pszFormat = info->pszFormat;
		projection = info->projection;
		setGeoTransform(info->adfGeoTransform);
	}
	void crop(OGREnvelope destBB, std::string outfile)
	{
		int startCol = (int)((destBB.MinX - bound.MinX) / adfGeoTransform[1]);
		int endCol = (int)((destBB.MaxX - bound.MinX) / adfGeoTransform[1]);
		int numCols = endCol - startCol + 1;

		int startRow = (int)((destBB.MaxY - bound.MaxY) / adfGeoTransform[5]);
		int endRow = (int)((destBB.MinY - bound.MaxY) / adfGeoTransform[5]);

		int numRows = endRow - startRow + 1;
		destBB.MinX = bound.MinX + startCol * adfGeoTransform[1];
		destBB.MaxY = bound.MaxY + startRow * adfGeoTransform[5];
		GDAL_DS* newdt = new GDAL_DS<T>();

		memcpy(newdt->adfGeoTransform, adfGeoTransform, sizeof(double) * 6);
		newdt->adfGeoTransform[0] = destBB.MinX;
		newdt->adfGeoTransform[3] = destBB.MaxY;
		newdt->bound = destBB;
		newdt->numbands = numbands;
		newdt->ncols = numCols;
		newdt->nrows = numRows;
		newdt->projection = this->projection;
		newdt->pszFormat = this->pszFormat;
		newdt->slice = numCols * numRows;

		std::vector<T*> datasets;
		std::vector<double> nodatas;
		for (size_t nband = 1; nband <= numbands; nband++)
		{
			T* srcdata = this->readData(nband);
			double nodata = this->getNoData(nband);
			nodatas.push_back(nodata);
			T* destdata = new T[numRows*numCols];
			int idx = 0;

			for (int irow = 0; irow < numRows; irow++)
			{
				for (int icol = 0; icol < numCols; icol++)
				{
					destdata[idx++] = nodata;
				}
			}

			for (int irow = 0; irow < numRows; irow++)
			{
				int srcRow = startRow + irow;
				if (srcRow > nrows - 1 || srcRow < 0)
					continue;
				for (int icol = 0; icol < numCols; icol++)
				{
					int srcCol = startCol + icol;
					if (srcCol > ncols - 1 || srcCol < 0)
						continue;
					destdata[icol + irow * numCols] = srcdata[srcCol + srcRow * this->ncols];
				}
			}
			datasets.push_back(destdata);
			delete[] srcdata;
		}

		this->close();
		newdt->create(outfile);

		for (size_t nband = 1; nband <= numbands; nband++)
		{
			double nodata = nodatas[nband - 1];
			T* destdata = datasets[nband-1];
			newdt->writeData(nband, destdata, nodata);
			delete[] destdata;
		}
		delete newdt;
		this->open(filename);
	}
	//void extractByMask(std::string maskfile, std::string outfile)
	//{
	//	GDAL_DS* maskdt = new GDAL_DS<T>();
	//	maskdt->open(maskfile);
	//	T* maskdata = maskdt->readData(1);
	//	GDAL_DS* newdt = new GDAL_DS<T>();
	//	newdt->setDSInfo(this);
	//	newdt->create(outfile);
	//	double masknodata = maskdt->getNoData(1);
	//
	//	for (size_t nband = 1; nband <= numbands; nband++)
	//	{
	//		double nodata = this->getNoData(nband);
	//		T* srcdata = this->readData(nband);
	//		T* destdata = new T[slice];
	//		int idx = 0;
	//		for (int ncell = 0; ncell < slice; ncell++)
	//		{
	//			T valmask = maskdata[ncell];
	//			if (valmask == masknodata)
	//			{
	//				destdata[idx] = nodata;
	//			}
	//			else
	//			{
	//				destdata[idx] = srcdata[ncell];
	//			}
	//			idx++;
	//		}

	//		newdt->writeData(nband, destdata, nodata);
	//		delete[] srcdata;
	//		delete[] destdata;
	//	}
	//	delete maskdt;
	//	delete[] maskdata;
	//	delete newdt;
	//}
	void extractByMask(std::string maskfile, std::string outfile)
	{
		GDAL_DS* maskdt = new GDAL_DS<T>();
		maskdt->open(maskfile);
		T* maskdata = maskdt->readData(1);
		//int startrow, int nrows, int startcol, int ncols;
		double resol = this->adfGeoTransform[1];
		int startrow = (this->bound.MaxY - maskdt->bound.MaxY) / resol;
		int startcol = (maskdt->bound.MinX - this->bound.MinX) / resol;

		double masknodata = maskdt->getNoData(1);
		GDAL_DS* newdt = new GDAL_DS<T>();
		newdt->setDSInfo(maskdt);
		newdt->numbands = this->numbands;
		std::vector<T*> datasets;
		std::vector<double> nodatas;
		for (size_t nband = 1; nband <= numbands; nband++)
		{
			double nodata = this->getNoData(nband);
			T* srcdata = this->readData(nband);
			T* destdata = new T[slice];
			for (int irow = 0; irow < newdt->nrows; irow++)
			{
				for (int icol = 0; icol < newdt->ncols; icol++)
				{
					int srcCell = (startcol + icol) + (startrow + irow) * this->ncols;
					int destCell = icol + irow * newdt->ncols;
					T valmask = maskdata[destCell];
					if (valmask == masknodata)
					{
						destdata[destCell] = nodata;
					}
					else
					{
						destdata[destCell] = srcdata[srcCell];
					}
				}
			}
			datasets.push_back(destdata);
			nodatas.push_back(nodata);
			delete[] srcdata;
		}
		this->close();
		newdt->create(outfile);
		for (size_t nband = 1; nband <= numbands; nband++)
		{
			double nodata = nodatas[nband - 1];
			T* destdata = datasets[nband - 1];
			newdt->writeData(nband, destdata, nodata);
			delete[] destdata;
		}
		delete maskdt;
		delete[] maskdata;
		delete newdt;
		this->open(filename);
	}

public:
	//OGREnvelope bound;
	//int nrows;
	//int ncols;
	//int slice;
	//int numbands;
	//double adfGeoTransform[6];
	//std::string filename;
	GDALDataset* dataset;
private:

};
//typedef enum {
//	/*! Unknown Orange unspecified type */          GDT_Unknown = 0,
//	/*! Eight bit unsigned integer */           GDT_Byte = 1,
//	/*! Sixteen bit unsigned integer */         GDT_UInt16 = 2,
//	/*! Sixteen bit signed integer */           GDT_Int16 = 3,
//	/*! Thirty two bit unsigned integer */      GDT_UInt32 = 4,
//	/*! Thirty two bit signed integer */        GDT_Int32 = 5,
//	/*! Thirty two bit floating point */        GDT_Float32 = 6,
//	/*! Sixty four bit floating point */        GDT_Float64 = 7,
//	/*! Complex Int16 */                        GDT_CInt16 = 8,
//	/*! Complex Int32 */                        GDT_CInt32 = 9,
//	/*! Complex Float32 */                      GDT_CFloat32 = 10,
//	/*! Complex Float64 */                      GDT_CFloat64 = 11,
//	GDT_TypeCount = 12          /* maximum type # + 1 */
//} GDALDataType;
template <class T>
class SVFComputer
{
public:
	GDAL_DS<T>* g_pDS;
	T* g_pData;
	osg::Vec4i skyColor;
	osg::Vec4i greenColor;
	osg::Vec4i outsideColor;
	SVFComputer(GDAL_DS<T>* ds) {
			g_pDS = ds;
			g_pData = ds->readData(1);
	}
	SVFComputer() {
		g_pDS = NULL;
		g_pData = NULL;
	}
	~SVFComputer()
	{
		delete[] g_pData;
	}
	void setData(GDAL_DS<T>* ds)
	{
		if (g_pData)
			delete[] g_pData;
		g_pDS = ds;
		g_pData = ds->readData(1);
	}
	double maxHorizon(osg::Vec3d eye, double azimuthAngle)
	{
		azimuthAngle = azimuthAngle * 3.1415926 / 180;
		osg::Vec2d dir(sin(azimuthAngle), cos(azimuthAngle));
		return maxHorizon(eye, dir);		
	}
	double maxHorizon(osg::Vec3d eye, osg::Vec2d dir)
	{
		osg::Vec2d curpos(eye.x(), eye.y());
		osg::Vec2d eyepos = curpos;
		double step = g_pDS->adfGeoTransform[1] * 0.5;
		OGREnvelope bb = g_pDS->bound;
		double resol = g_pDS->adfGeoTransform[1];
		double maxhorizon = -180;
		double dist = 0;
		double h = 0;
		while (true)
		{
			curpos = curpos + dir * step;
			dist += step;
			if (dist > 1)
			{
				int irow = (int)(bb.MaxY - curpos.y()) / resol;
				int icol = (int)(curpos.x() - bb.MinX) / resol;
				if (irow < 0 || irow > g_pDS->nrows - 1 || icol < 0 || icol > g_pDS->ncols - 1)
					break;
				double run = (eyepos - curpos).length();
				T height = g_pData[icol + irow*g_pDS->ncols];
				double rise = height - eye.z();
				double horizon = atan2(rise, run) * 180 / 3.1415926;
				if (horizon > maxhorizon)
				{
					maxhorizon = horizon;
					h = height;
				}
			}


		}
		if (maxhorizon < 0)
			maxhorizon = 0;
		if (maxhorizon > 90)
			maxhorizon = 90;
		printf("%f,%f,%f,%f\n", dir.x(), dir.y(), h, maxhorizon);
		return maxhorizon;
	}
	std::vector<osg::Vec2d> computeSkymap(osg::Vec3d eye, double numsteps)
	{
		std::vector<osg::Vec2d> coords;
		double step = 360.0 / numsteps;
		for (size_t i = 0; i < numsteps; i++)
		{
			double azimuth = step * i;
			double azimuthInRadians = azimuth * 3.1415926 / 180;
			osg::Vec2d dir(sin(azimuthInRadians), cos(azimuthInRadians));
			double maxhorizon = maxHorizon(eye, dir);
			double len = (90.0 - maxhorizon) / 90.0;
			osg::Vec2d pos(dir.x()*len, dir.y()*len);
			coords.push_back(pos);
		}
		return coords;
	}
	void drawSkymap(std::vector<osg::Vec2d> points, int imagesize, std::string outfile)
	{
		QImage* qImage = new QImage(imagesize, imagesize, QImage::Format::Format_ARGB32);
		//QImage* qImage2 = new QImage(imagesize, imagesize, QImage::Format::Format_ARGB32);
		qImage->fill(qRgba(0, 0, 0, 0));
		//qImage2->fill(qRgba(0, 0, 0, 0));
		QPainter painter(qImage);
		painter.begin(qImage);
		//painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
		//painter.setRenderHint(QPainter::Antialiasing); //make it look nicer
		QPointF center(imagesize*0.5, imagesize*0.5);
	
		//QTransform transform = QTransform().scale(imagesize*0.5, imagesize*0.5).translate(imagesize*0.5, imagesize*0.5);
		//QTransform transform = QTransform().scale(imagesize*0.5, imagesize*0.5)
		QTransform transform = QTransform().scale(imagesize*0.5, imagesize*0.5) * QTransform().translate(imagesize*0.5, imagesize*0.5);
		painter.setTransform(transform);
		painter.setBrush(QBrush(qRgba(255, 0, 0, 255)));
		painter.setPen(Qt::NoPen);
		painter.drawEllipse(QPointF(0, 0), 1, 1);

		QPainterPath qpath;
		QPolygonF qpolygon;
		for (size_t i = 0; i < points.size(); i++)
		{
			osg::Vec2d pt = points[i];
			qpolygon.push_back(QPointF(pt.x(), pt.y()));
		}
		qpolygon.push_back(QPointF(points[0].x(), points[0].y()));
		qpath.addPolygon(qpolygon);
		//painter.setPen(Qt::NoPen);
		QPen pen(qRgba(0, 0, 0, 0));
		pen.setWidth(0);
		painter.setPen(pen);
		painter.fillPath(qpath, QBrush(qRgba(0, 255, 0, 255)));
		painter.end();
		QString str = outfile.data();
		std::string fix = "png";
		if (str.endsWith("jpg"))
			fix = "JPG";
		else if (str.endsWith("bmp"))
			fix = "bmp";
		//osg::ref_ptr<osg::Image> image = new osg::Image;
		//image->allocateImage(qImage->width(), qImage->height(), 1, 4, GL_BGRA, GL_UNSIGNED_BYTE);
		//The pixel format is always RGBA to support transparency
		//memset(image->srcdata(), 0, qImage->byteCount());
		//image->setImage(qImage->width(), qImage->height(), 1,
		//	4,
		//	GL_RGBA, GL_UNSIGNED_BYTE, //Why not GL_RGBA - QGIS bug?
		//	srcdata,
		//	osg::Image::USE_NEW_DELETE, 1);
		//memset(srcdata, 0, qImage->byteCount());
		//image->flipVertical();
		//osgDB::writeImageFile(*image, outfile.srcdata());
		int elementSize = 4;
		uchar* data = new uchar[qImage->width() * qImage->height() * elementSize];
		uchar* pdata = qImage->bits();
		for (size_t i = 0; i < qImage->height(); i++)
		{
			for (size_t j = 0; j < qImage->width(); j++)
			{
				//BGRA
				if (pdata[2] == 255)
					pdata[3] = 255;
				else if (pdata[1] == 255)
					pdata[3] = 128;
				else
					pdata[3] = 0;
				pdata += elementSize;
			}
		}

		for (size_t i = 0; i < qImage->height(); i++)
		{
			int oriLine = qImage->height() - i - 1;
			memcpy(data + i * qImage->width() * elementSize, qImage->bits() + oriLine*qImage->width() * elementSize, qImage->width() * elementSize);
		}
		memcpy(qImage->bits(), data, qImage->width() * qImage->height() * elementSize);
		delete[] data;
		qImage->save(str, fix.data());
		delete qImage;
	}
	T getHeightAt(double x, double y)
	{
		double resol = g_pDS->adfGeoTransform[1];
		int irow = (int)(g_pDS->bound.MaxY - y) / resol;
		int icol = (int)(x - g_pDS->bound.MinX) / resol;
		if (irow < 0 || irow > g_pDS->nrows - 1 || icol < 0 || icol > g_pDS->ncols - 1)
			return g_pDS->getNoData(1);
		return g_pData[icol + irow*g_pDS->ncols];
	}
	//double calSVF(bool applyLambert = false)
	//{
	//	unsigned int skypixels = 0;
	//	unsigned int nonskypixels = 0;
	//	unsigned int ncols = g_pDS->ncols;
	//	unsigned int nrows = g_pDS->nrows;
	//	unsigned int numpixels = ncols * nrows;
	//	if (ncols != nrows)
	//		return 0;
	//	T* srcdata = g_pDS->readData(4);
	//	double resol = 1.0 / nrows;
	//	double totalarea = 0;
	//	double skyarea = 0;
	//	double y = resol * 0.5 - 0.5;
	//	int npixel = 0;
	//	for (unsigned int row = 0; row < nrows; row++)
	//	{
	//		y += resol;
	//		double x = resol * 0.5 - 0.5;
	//		for (unsigned int col = 0; col < ncols; col++)
	//		{
	//			T a = srcdata[npixel];

	//			//printf("%d,", (int)a);
	//			npixel++;
	//			if (outsideColor.a() && r == greenColor.r() && g == greenColor.g() && b == greenColor.b()) {
	//				x += resol;
	//				continue;//outside
	//			}
	//			double zenithD = sqrt(x*x + y*y) * 90.0;//in degrees
	//			if (zenithD <= 0.000000001)
	//				zenithD = 0.000000001;
	//			double zenithR = zenithD * 3.1415926 / 180.0;
	//			double wproj = sin(zenithR) / (zenithD / 90);//weight for equal-areal projection
	//			if (applyLambert)
	//			{
	//				wproj = wproj * cos(zenithR);
	//			}
	//			totalarea += wproj;
	//			if (a == 255)
	//			{
	//				nonskypixels++;
	//			}
	//			else
	//			{
	//				skypixels++;
	//				skyarea += wproj;
	//			}
	//			x += resol;
	//		}

	//	}
	//	double svf = skyarea / totalarea;
	//	delete[] srcdata;
	//	return svf;
	//}
	double calSVF(bool applyLambert = false)
	{
		unsigned int skypixels = 0;
		unsigned int nonskypixels = 0;
		unsigned int ncols = g_pDS->ncols;
		unsigned int nrows = g_pDS->nrows;
		unsigned int numpixels = ncols * nrows;
		if (ncols != nrows)
			return 0;
		T* data = g_pDS->readData(4);
		T* rdata = g_pDS->readData(1);
		T* gdata = g_pDS->readData(2);
		T* bdata = g_pDS->readData(3);

		double resol = 1.0 / nrows;
		double totalarea = 0;
		double skyarea = 0;
		double y = resol * 0.5 - 0.5;
		int npixel = 0;
		for (unsigned int row = 0; row < nrows; row++)
		{
			y += resol;
			double x = resol * 0.5 - 0.5;
			for (unsigned int col = 0; col < ncols; col++)
			{
				T a = data[npixel];
				T r = rdata[npixel];
				T g = gdata[npixel];
				T b = bdata[npixel];
				//printf("%d,", (int)a);
				npixel++;
				if (a == outsideColor.a() && r == outsideColor.r() && g == outsideColor.g() && b == outsideColor.b()) {
					x += resol;
					continue;//outside
				}
				double zenithD = sqrt(x*x + y*y) * 90.0;//in degrees
				if (zenithD <= 0.000000001)
					zenithD = 0.000000001;
				double zenithR = zenithD * 3.1415926 / 180.0;
				double wproj = sin(zenithR) / (zenithD / 90);//weight for equal-areal projection
				if (applyLambert)
				{
					wproj = wproj * cos(zenithR);
				}
				totalarea += wproj;
				if ( g > 120)
				{
					printf("");
				}
				if (a == skyColor.a() && r == skyColor.r() && g == skyColor.g() && b == skyColor.b())
				{
					skypixels++;
					skyarea += wproj;
				}
				else
				{
					nonskypixels++;
				}
				x += resol;
			}

		}
		double svf = skyarea / totalarea;
		delete[] data;
		delete[] rdata;
		delete[] gdata;
		delete[] bdata;
		return svf;
	}
	double calGreenery(bool applyLambert = false)
	{
		unsigned int skypixels = 0;
		unsigned int nonskypixels = 0;
		unsigned int ncols = g_pDS->ncols;
		unsigned int nrows = g_pDS->nrows;
		unsigned int numpixels = ncols * nrows;
		if (ncols != nrows)
			return 0;
		T* data = g_pDS->readData(4);
		T* rdata = g_pDS->readData(1);
		T* gdata = g_pDS->readData(2);
		T* bdata = g_pDS->readData(3);

		double resol = 1.0 / nrows;
		double totalarea = 0;
		double skyarea = 0;
		double y = resol * 0.5 - 0.5;
		int npixel = 0;
		for (unsigned int row = 0; row < nrows; row++)
		{
			y += resol;
			double x = resol * 0.5 - 0.5;
			for (unsigned int col = 0; col < ncols; col++)
			{
				T a = data[npixel];
				T r = rdata[npixel];
				T g = gdata[npixel];
				T b = bdata[npixel];
				//printf("%d,", (int)a);
				npixel++;
				if (a == outsideColor.a() && r == outsideColor.r() && g == outsideColor.g() && b == outsideColor.b()) {
					x += resol;
					continue;//outside
				}
				double zenithD = sqrt(x*x + y*y) * 90.0;//in degrees
				if (zenithD <= 0.000000001)
					zenithD = 0.000000001;
				double zenithR = zenithD * 3.1415926 / 180.0;
				double wproj = sin(zenithR) / (zenithD / 90);//weight for equal-areal projection
				if (applyLambert)
				{
					wproj = wproj * cos(zenithR);
				}
				totalarea += wproj;
				if (a == greenColor.a() && r == greenColor.r() && g == greenColor.g() && b == greenColor.b())
				{
					skypixels++;
					skyarea += wproj;
				}
				else
				{
					nonskypixels++;
				}
				x += resol;
			}

		}
		double svf = skyarea / totalarea;
		delete[] data;
		delete[] rdata;
		delete[] gdata;
		delete[] bdata;
		return svf;
	}

	//double maxHorizon(osg::Vec3d eye, osg::Vec2 dir, double z);
};
//vec3 spherical2Cartisian(double lon, double lat,double )
//{
//	double theta = lon * 0.0174533;
//	double phi =   lat* 0.0174533;
//	return vec3(cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi));
//}
