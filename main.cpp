// PolygonStatistics.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ogrsf_frmts.h"
#include "gdal_priv.h"
#include <QFileinfo>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <qdir.h>
#include "qfileinfo.h"
#include "DbfFile.h"
#include "ShapeFile.h"
#include "BoundManager.h"
#include "Grid.h"
#include "Utils.h"
#include "Preprocessor.h"
#include "ShapeCreator.h"
#include "ScaleByFuel.h"
#include "NonPointProcessor.h"
#include "RoadTimeIDW.h"
#include "EMFAC.h"
#include "TimestructTool.h"
#include "Reallocator.h"
#include "MarionElementa.h"
#include "Geographies.h"
#include "BatchDOE.h";
const std::string SMOKE_INDIR = "B:/NonRoadTemporal/NonRoadTemporal/";
const std::string SMOKE_OUTDIR = "B:/LA_Version2/Time/SCC/";
bool copy_SCC_PROFILE(std::string inname, std::string outname)
{
	std::string infile = SMOKE_INDIR + inname + ".txt";
	if (!QFileInfo(infile.data()).exists())
		return false;
	std::string outfile = SMOKE_OUTDIR + outname + ".txt";
	if (!QFileInfo(outfile.data()).exists())
		return true;
	QFile::copy(infile.data(), outfile.data());
	return true;
}
void verifyRailRoad()
{

	double* sums = new double[1000000];
	memset(sums, 0, sizeof(double) * 1000000);

	ShapeFile rail1("B:/Baltimore/gridPrep_SHP_master/WGS84/OnRoad.shp");
	ShapeFile rail2("B:/Baltimore/gridPrep_SHP_master/WGS84/result/OnRoad.shp");

	OGRFeature *poFeature;
	rail2.poLayer->ResetReading();
	int id = 0;
	while ((poFeature = rail2.poLayer->GetNextFeature()) != NULL)
	{
		//int id = poFeature->GetFieldAsInteger("FID_OnRoad");
		//printf("%d\n", id);
		//sums[id] = sums[id] + poFeature->GetFieldAsDouble("ca");
		OGRFeature::DestroyFeature(poFeature);
	}


	rail1.poLayer->ResetReading();
	id = 0;
	while ((poFeature = rail1.poLayer->GetNextFeature()) != NULL)
	{
		double sum = poFeature->GetFieldAsDouble("ca");
		if (abs(sums[id] - sum) > 0.001)
		{
			printf("%d,%.5f,%.5f\n", id, sum, sums[id]);
		}
		id++;
		OGRFeature::DestroyFeature(poFeature);
	}


	//"FID_Railro"
}
void sum(std::string dir)
{


	std::vector<std::string> files;

	QDir input_dir(dir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	std::string indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		if (fileInfo.fileName().endsWith("fishnet.shp"))
			continue;


		files.push_back(fileInfo.fileName().toLocal8Bit().data());
	}

	double sum = 0;

	for (size_t i = 0; i < files.size(); i++)
	{
		double subtotal = 0;
		ShapeFile shp(indir + files[i]);
		OGRFeature *poFeature;
		shp.poLayer->ResetReading();
		int idx = shp.poLayer->GetLayerDefn()->GetFieldIndex("ca");
		while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
		{
			sum += poFeature->GetFieldAsDouble(idx);
			subtotal += poFeature->GetFieldAsDouble(idx);
			OGRFeature::DestroyFeature(poFeature);
		}
		printf("sector %s: %.5f\n", files[i].data(), subtotal / 1000000 / 1000);
	}



	printf("%.5f\n", sum);
	printf("%.5f\n", sum / 1000000 / 1000);

	getchar();

}
//typedef enum
//{
//	wkbUnknown = 0,         /**< unknown type, non-standard */
//	wkbPoint = 1,           /**< 0-dimensional geometric object, standard WKB */
//	wkbLineString = 2,      /**< 1-dimensional geometric object with linear
//							*   interpolation between Points, standard WKB */
//	wkbPolygon = 3,         /**< planar 2-dimensional geometric object defined
//							*   by 1 exterior boundary and 0 Orange more interior
//							*   boundaries, standard WKB */
//	wkbMultiPoint = 4,      /**< GeometryCollection of Points, standard WKB */
//	wkbMultiLineString = 5, /**< GeometryCollection of LineStrings, standard WKB */
//	wkbMultiPolygon = 6,    /**< GeometryCollection of Polygons, standard WKB */
//	wkbGeometryCollection = 7, /**< geometric object that is a collection of 1
//							   Orange more geometric objects, standard WKB */
//	wkbNone = 100,          /**< non-standard, for pure attribute records */
//	wkbLinearRing = 101,    /**< non-standard, just for createGeometry() */
//	wkbPoint25D = 0x80000001, /**< 2.5D extension as per 99-402 */
//	wkbLineString25D = 0x80000002, /**< 2.5D extension as per 99-402 */
//	wkbPolygon25D = 0x80000003, /**< 2.5D extension as per 99-402 */
//	wkbMultiPoint25D = 0x80000004, /**< 2.5D extension as per 99-402 */
//	wkbMultiLineString25D = 0x80000005, /**< 2.5D extension as per 99-402 */
//	wkbMultiPolygon25D = 0x80000006, /**< 2.5D extension as per 99-402 */
//	wkbGeometryCollection25D = 0x80000007 /**< 2.5D extension as per 99-402 */
//} OGRwkbGeometryType;
void updateFieldAfterIntersection(std::string filename)
{
	ShapeFile input(filename, 1);
	OGRFeature *poFeature;
	input.poLayer->ResetReading();
	std::vector<int> fields;

	for (size_t i = 0; i < input.poLayer->GetLayerDefn()->GetFieldCount(); i++)
	{
		OGRFieldDefn* field = input.poLayer->GetLayerDefn()->GetFieldDefn(i);
		std::string fieldname = field->GetNameRef();
		if (fieldname.size() > 1 && QString(field->GetNameRef()).toLower().indexOf("ca") > -1)
			fields.push_back(input.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data()));
	}

	int idIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("Id");

	OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
	int footprintIndex = -1;
	if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
	{
		footprintIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("length");
	}
	else if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		footprintIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("area");
	}
	else
	{
		return;
	}
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* geo = poFeature->GetGeometryRef();
		OGRwkbGeometryType gtype = geo->getGeometryType();

		double footprintOld = poFeature->GetFieldAsDouble(footprintIndex);
		double footprintNew = 0;

		if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
		{
			footprintNew = Utils::calPolylineLength(poFeature->GetGeometryRef());
		}
		else if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype == wkbPolygon25D)
		{
			footprintNew = Utils::calPolygonArea(poFeature->GetGeometryRef());
		}

		double fraction = 0;
		if (!isinf(footprintOld) && footprintOld > 0)
			fraction = footprintNew / footprintOld;

		for (size_t i = 0; i < fields.size(); i++)
		{
			double ca = poFeature->GetFieldAsDouble(fields[i]);
			if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
			{
				poFeature->SetField(fields[i], ca*fraction);
			}
			else if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype == wkbPolygon25D)
			{
				poFeature->SetField(fields[i], ca*fraction);
			}
		}

		poFeature->SetField(footprintIndex, footprintNew);
		input.poLayer->SetFeature(poFeature);


		OGRFeature::DestroyFeature(poFeature);

	}

}
void updateFracAfterIntersection(std::string filename)
{
	ShapeFile input(filename, 1);
	OGRFeature *poFeature;
	input.poLayer->ResetReading();


	int idIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("Id");

	OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
	int footprintIndex = -1;
	if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
	{
		footprintIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("length");
	}
	else if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
	{
		footprintIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("area");
	}
	else
	{
		return;
	}
	int fracIndex = input.getOrCreateField("frac",OGRFieldType::OFTReal);

	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* geo = poFeature->GetGeometryRef();
		OGRwkbGeometryType gtype = geo->getGeometryType();

		double footprintOld = poFeature->GetFieldAsDouble(footprintIndex);
		double footprintNew = 0;

		if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
		{
			footprintNew = Utils::calPolylineLength(poFeature->GetGeometryRef());
		}
		else if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype == wkbPolygon25D)
		{
			footprintNew = Utils::calPolygonArea(poFeature->GetGeometryRef());
		}

		double fraction = 0;
		if (!isinf(footprintOld) && footprintOld > 0)
			fraction = footprintNew / footprintOld;
		poFeature->SetField(fracIndex, fraction);
		poFeature->SetField(footprintIndex, footprintNew);
		input.poLayer->SetFeature(poFeature);


		OGRFeature::DestroyFeature(poFeature);

	}

}
void updateFieldAfterIntersectionForDir(std::string indir)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (isdigit(name[0]))
			continue;
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		if (fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
		std::string infile = indir + files[files.size() - 1];
		updateFieldAfterIntersection(infile);
	}
}
void intersectWithPolygonForDir(std::string boundary, std::string indir, std::string outdir)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QDir output_dir(outdir.data());
	if (!output_dir.exists())
		output_dir.mkpath(".");
	outdir = (output_dir.absolutePath() + "/").toLocal8Bit().data();


	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (isdigit(name[0]))
			continue;
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		if (fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
		std::string infile = indir + files[files.size() - 1];
		std::string outfile = outdir + files[files.size() - 1];
		//intersectWithArcGIS(infile, boundary, outfile);
		updateFieldAfterIntersection(outfile);
		printf("%s\n", infile.data());
	}



}

void SubsetVulcan(std::string boundary, std::string indir, std::string outdir)
{
	//Utils::updateFootprintForDir(indir,true);
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QDir output_dir(outdir.data());
	if (!output_dir.exists())
		output_dir.mkpath(".");
	outdir = (output_dir.absolutePath() + "/").toLocal8Bit().data();


	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (isdigit(name[0]))
			continue;
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		if (fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
		std::string infile = indir + files[files.size() - 1];
		std::string outfile = outdir + files[files.size() - 1];
		Preprocessor::intersectWithArcGIS(infile, boundary, outfile);
		updateFracAfterIntersection(outfile);
		printf("%s\n", infile.data());
	}



}
OGREnvelope splitIntoTiles(std::string indir, std::string outdir, OGREnvelope bound, double gridsize)
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
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (isdigit(name[0]))
			continue;
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
		std::string infile = indir + files[files.size() - 1];
		//updateFootprint(infile);
	}
	std::string fishnetfile = outdir + "fishnet.shp";
	ShapeFile spatialRefFile(indir + files[0]);
	if (!spatialRefFile.poLayer->GetSpatialRef()->IsGeographic())
	{
		//FOOTPRINT_SCALE_FACTOR = 1;
	}
	else
	{
		//FOOTPRINT_SCALE_FACTOR = getDegreeToMeter(spatialRefFile.poLayer);
		gridsize = gridsize / Utils::getDegreeToMeter(spatialRefFile.poLayer);
	}
	Grid grid(bound, gridsize, 1);
	//if(!QFileInfo(fishnetfile.outlineDS()).exists())
	char wkt[512];
	char* pwkt = wkt;
	if (spatialRefFile.poLayer->GetSpatialRef())
		spatialRefFile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	grid.toShape(wkt, fishnetfile, false);
	for (int n = 0; n < grid.nrows*grid.ncols; n++)
	{
		std::stringstream subboundss;
		subboundss << outdir << n << ".shp";
		//if (QFileInfo(subboundss.str().outlineDS()).exists())
		//	continue;
		if (!QFileInfo(subboundss.str().data()).exists())
		{
			char wkt[512];
			char* pwkt = wkt;
			if (spatialRefFile.poLayer->GetSpatialRef())
				spatialRefFile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
			grid.toShape(wkt, subboundss.str().data(), n);
		}

		std::stringstream subdirss;
		subdirss << outdir << n;
		QDir qsubdir(subdirss.str().data());
		if (!qsubdir.exists())
			qsubdir.mkpath(".");
		else
			continue;
		std::string subdir = (qsubdir.absolutePath() + "/").toLocal8Bit().data();

		for (int i = 0; i < files.size(); i++)
		{

			std::string outputfile = subdir + files[i];
			std::string tagfile = outputfile + ".locked";

			if (QFileInfo(outputfile.data()).exists() || QFileInfo(tagfile.data()).exists())
				continue;

			std::ofstream ofs(tagfile.data());
			ofs.close();
			printf("%s\n", outputfile.data());

			Preprocessor::clipWithArcGIS(indir + files[i], subboundss.str().data(), outputfile);
			updateFieldAfterIntersection(outputfile);
			if (QFileInfo(tagfile.data()).exists())
				QFile::remove(tagfile.data());
		}

		//for (int i = 0; i < files.size(); i++)
		//{

		//	std::string outputfile = subdir + files[i];
		//	std::string tagfile = outputfile + ".locked";

		//	if (QFileInfo(outputfile.outlineDS()).exists() || QFileInfo(tagfile.outlineDS()).exists())
		//		continue;

		//	std::ofstream ofs(tagfile.outlineDS());
		//	ofs.close();
		//	printf("%s\n", outputfile.outlineDS());
		//	clipWithArcGIS(indir + files[i], subboundss.str().outlineDS(), outputfile);
		//	if (QFileInfo(tagfile.outlineDS()).exists())
		//		QFile::remove(tagfile.outlineDS());
		//}


	}
	return grid.bound;
}
OGREnvelope intersectByTiles(std::string indir, std::string outdir, OGREnvelope bound, double gridsize, double tilegridsize)
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
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (isdigit(name[0]))
			continue;
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
	}
	std::string fishnetfile = outdir + "fishnet.shp";
	ShapeFile spatialRefFile(indir + files[0]);
	if (!spatialRefFile.poLayer->GetSpatialRef()->IsGeographic())
	{
		//FOOTPRINT_SCALE_FACTOR = 1;
	}
	else
	{
		//FOOTPRINT_SCALE_FACTOR = getDegreeToMeter(spatialRefFile.poLayer);
		gridsize = gridsize / Utils::getDegreeToMeter(spatialRefFile.poLayer);
	}
	Grid grid(bound, gridsize, 1);
	//if (!QFileInfo(fishnetfile.outlineDS()).exists())
	//	grid.toShape(&spatialRefFile, fishnetfile, false);
	for (int n = 0; n < grid.nrows*grid.ncols; n++)
	{
		std::stringstream subboundss;
		subboundss << outdir << n << ".shp";
		//if (QFileInfo(subboundss.str().outlineDS()).exists())
		//	continue;
		//if (!QFileInfo(subboundss.str().outlineDS()).exists())
		//	grid.toShape(&spatialRefFile, subboundss.str().outlineDS(), n);

		OGREnvelope tilebound = BoundManager::readBoundFromShape(subboundss.str());

		std::stringstream subdirss;
		subdirss << outdir << n;
		QDir qsubdir(subdirss.str().data());
		//if (!qsubdir.exists())
		//	qsubdir.mkpath(".");
		std::string subdir = (qsubdir.absolutePath() + "/").toLocal8Bit().data();


		std::stringstream subdirssout;
		subdirssout << outdir << n << "/" << tilegridsize << "m/";
		QDir qsubdirout(subdirssout.str().data());
		if (!qsubdirout.exists())
			qsubdirout.mkpath(".");
		else
			continue;
		std::string subdirout = (qsubdirout.absolutePath() + "/").toLocal8Bit().data();
		ShapeFile tilespatialRefFile(subboundss.str());
		Grid tilegrid(tilebound, tilegridsize, 0);
		std::string tilefishnetfile = subdirout + "fishnet.shp";
		if (!QFileInfo(tilefishnetfile.data()).exists())
		{
			char wkt[512];
			char* pwkt = wkt;
			if (tilespatialRefFile.poLayer->GetSpatialRef())
				tilespatialRefFile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);

			tilegrid.toShape(pwkt, tilefishnetfile.data());
		}

		for (int i = 0; i < files.size(); i++)
		{

			std::string inputfile = subdir + files[i];
			std::string outputfile = subdirout + files[i];

			//std::string tagfile = outputfile + ".locked";

			//if (QFileInfo(outputfile.outlineDS()).exists() || QFileInfo(tagfile.outlineDS()).exists())
			//continue;

			//std::ofstream ofs(tagfile.outlineDS());
			//ofs.close();
			printf("%s\n", outputfile.data());
			//clipWithArcGIS(indir + files[i], subboundss.str().outlineDS(), outputfile);
			Preprocessor::intersectWithArcGIS(tilefishnetfile, inputfile, outputfile);
			//if (QFileInfo(tagfile.outlineDS()).exists())
			//QFile::remove(tagfile.outlineDS());


		}
		QFile::remove(tilefishnetfile.data());
	}
	return grid.bound;
}
void gridFolderByRaster(std::string indir, std::string outdir, double gridsize, bool skipNonRoad = false)
{
	QDir qoutdir(outdir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	std::vector<std::string> files;
	std::vector<std::string> files2;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
		files2.push_back(input_file);
	}

	std::string nonroadfile = indir + "NonRoad.shp";
	std::string fishnetfile = outdir + "fishnet.shp";
	std::string rasterfile = outdir + "fishnet.tif";

	OGREnvelope bound = BoundManager::readBoundFromShapes(files2);

	ShapeFile boundfile(files2[0].data());
	if (!boundfile.poLayer->GetSpatialRef()->IsGeographic())
	{
		//FOOTPRINT_SCALE_FACTOR = 1;
	}
	else
	{
		//FOOTPRINT_SCALE_FACTOR = getDegreeToMeter(boundfile.poLayer);
		gridsize = gridsize / Utils::getDegreeToMeter(boundfile.poLayer);
	}


	Grid* fishnet = new Grid(bound, gridsize, 1);
	fishnet->reset();

	char wkt[512];
	char* pwkt = wkt;
	if (boundfile.poLayer->GetSpatialRef())
		boundfile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);

	fishnet->toShape(wkt, fishnetfile, false);

	for (size_t i = 0; i < files.size(); i++)
	{
		if (skipNonRoad && files[i] == "NonRoad.shp")
			continue;
		printf("%s\n", (outdir + files[i]).data());
		//intersectWithFishnet(fishnet, fishnetfile, indir + files[i], lidarFisheyeDir + files[i]);
		if (!QFileInfo((outdir + files[i]).data()).exists())
		{
			Preprocessor::intersectWithArcGIS(fishnetfile, indir + files[i], outdir + files[i]);
			updateFieldAfterIntersection(outdir + files[i]);
		}

		ShapeFile input((outdir + files[i]).data());
		fishnet->gatherCells(&input);
	}
	memset(wkt, 0, 512);
	pwkt = wkt;
	if (boundfile.poLayer->GetSpatialRef())
		boundfile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	fishnet->toShape(wkt, fishnetfile, true);
	fishnet->toRaster(rasterfile, wkt);
	boundfile.close();

	//for (size_t i = 0; i < files.size(); i++)
	//{
	//	if (skipNonRoad && files[i] == "NonRoad.Shp")
	//		continue;
	//	ShapeFile input((lidarFisheyeDir + files[i]).outlineDS());
	//	if (input.poLayer->GetGeomType() == wkbPoint || input.poLayer->GetGeomType() == wkbMultiPoint || input.poLayer->GetGeomType() == wkbPoint25D)
	//		continue;
	//	std::string sectorfile = files[i].substr(0, files[i].length() - 4);
	//	std::string sectorfileout = lidarFisheyeDir + sectorfile + ".tif";
	//	printf("%s\n", sectorfile.outlineDS());
	//	fishnet->reset();
	//	gatherCells(fishnet, &input);
	//	fishnet->toRaster(sectorfileout);

	//	sectorfileout = lidarFisheyeDir + sectorfile + "2.Shp";
	//	fishnet->toShape(&input, sectorfileout,true);

	//}

	delete fishnet;
}
void gridFolderByRaster(std::string indir, std::string outdir, std::string intersectFile)
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
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
	}

	//Utils::updateFootprint(intersectFile);
	for (size_t i = 0; i < files.size(); i++)
	{

		//updateFieldAfterIntersection(outdir + files[i]);
		//intersectWithFishnet(fishnet, fishnetfile, indir + files[i], lidarFisheyeDir + files[i]);
		if (!QFileInfo((outdir + files[i]).data()).exists())
		{
			Preprocessor::intersectWithArcGIS(intersectFile, indir + files[i], outdir + files[i]);
			updateFieldAfterIntersection(outdir + files[i]);
		}

	}

}
void findShapes(std::string indir, std::vector<std::string>& files)
{

	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.absoluteFilePath().toLocal8Bit().data());
	}

}
void gridFiles(std::vector<std::string> files, std::string outdir, OGREnvelope bb, double gridsize)
{
	QDir qoutdir(outdir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();

	std::string fishnetfile = outdir + "fishnet.shp";
	std::string rasterfile = outdir + "fishnet.tif";
	ShapeFile boundfile(files[0]);

	Grid* fishnet = new Grid(bb, gridsize, 1);
	fishnet->reset();
	char wkt[512];
	char* pwkt = wkt;
	if (boundfile.poLayer->GetSpatialRef())
		boundfile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	fishnet->toShape(wkt, fishnetfile, false);
	std::vector<std::string> comfiles;
	std::vector<std::string> indfiles;
	std::vector<std::string> resfiles;
	std::vector<std::string> onroadfiles;

	for (size_t i = 0; i < files.size(); i++)
	{
		QFileInfo fileinfo(files[i].data());
		std::string dir = fileinfo.absoluteDir().absolutePath().toLocal8Bit().data();
		std::string dirname = QFileInfo(dir.data()).baseName().toLocal8Bit().data();
		std::string name = fileinfo.fileName().toLocal8Bit().data();
		std::string outfile = outdir + dirname + "_" + name;
		printf("%s\n", outfile.data());
		//intersectWithFishnet(fishnet, fishnetfile, indir + files[i], lidarFisheyeDir + files[i]);
		if (!QFileInfo(outfile.data()).exists())
		{
			Preprocessor::intersectWithArcGIS(fishnetfile, files[i], outfile);
			updateFieldAfterIntersection(outfile);
		}
		QString qoutfile = outfile.data();
		qoutfile = qoutfile.toLower();
		if (qoutfile.endsWith("compoint.shp") || qoutfile.endsWith("comnonpoint.shp"))
			comfiles.push_back(outfile);
		if (qoutfile.endsWith("resnonpoint.shp"))
			resfiles.push_back(outfile);
		if (qoutfile.endsWith("indpoint.shp") || qoutfile.endsWith("indnonpoint.shp"))
			indfiles.push_back(outfile);
		if (qoutfile.endsWith("onroad.shp") && !qoutfile.endsWith("nonroad.shp"))
			onroadfiles.push_back(outfile);
		ShapeFile input(outfile);
		fishnet->gatherCells(&input);
	}

	fishnet->toShape(pwkt, fishnetfile, true);
	fishnet->toRaster(rasterfile, pwkt);
	fishnet->reset();
	for (size_t i = 0; i < onroadfiles.size(); i++)
	{
		ShapeFile input(onroadfiles[i]);
		fishnet->gatherCells(&input);
	}
	fishnet->toShape(pwkt, outdir + "fishnet_onroad.shp", true);
	fishnet->toRaster(outdir + "fishnet_onroad.tif", pwkt);

	fishnet->reset();
	for (size_t i = 0; i < comfiles.size(); i++)
	{
		ShapeFile input(comfiles[i]);
		fishnet->gatherCells(&input);
	}

	fishnet->toShape(wkt, outdir + "fishnet_com.shp", true);
	fishnet->toRaster(outdir + "fishnet_com.tif", wkt);

	fishnet->reset();
	for (size_t i = 0; i < resfiles.size(); i++)
	{
		ShapeFile input(resfiles[i]);
		fishnet->gatherCells(&input);
	}
	fishnet->toShape(pwkt, outdir + "fishnet_res.shp", true);
	fishnet->toRaster(outdir + "fishnet_res.tif", pwkt);

	fishnet->reset();
	for (size_t i = 0; i < indfiles.size(); i++)
	{
		ShapeFile input(indfiles[i]);
		fishnet->gatherCells(&input);
	}
	fishnet->toShape(pwkt, outdir + "fishnet_ind.shp", true);
	fishnet->toRaster(outdir + "fishnet_ind.tif", pwkt);



	boundfile.close();

	//for (size_t i = 0; i < files.size(); i++)
	//{
	//	if (skipNonRoad && files[i] == "NonRoad.Shp")
	//		continue;
	//	ShapeFile input((lidarFisheyeDir + files[i]).outlineDS());
	//	if (input.poLayer->GetGeomType() == wkbPoint || input.poLayer->GetGeomType() == wkbMultiPoint || input.poLayer->GetGeomType() == wkbPoint25D)
	//		continue;
	//	std::string sectorfile = files[i].substr(0, files[i].length() - 4);
	//	std::string sectorfileout = lidarFisheyeDir + sectorfile + ".tif";
	//	printf("%s\n", sectorfile.outlineDS());
	//	fishnet->reset();
	//	gatherCells(fishnet, &input);
	//	fishnet->toRaster(sectorfileout);

	//	sectorfileout = lidarFisheyeDir + sectorfile + "2.Shp";
	//	fishnet->toShape(&input, sectorfileout,true);

	//}

	delete fishnet;
}
#include "GDAL_DS.h"
void aggregrateSingleFile(std::string filenamein, std::string filename)
{
	std::map<int, double> tb;
	std::vector<std::string> files;
	OGRFeature *poFeature;
	int ca11field;
	int id = 0;

	ShapeFile shp(filenamein);
	int blockgroupid = shp.getField("FID_blockg");

	ca11field = shp.getField("ca11");
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		id = poFeature->GetFieldAsInteger(blockgroupid);
		double val = poFeature->GetFieldAsDouble(ca11field);
		std::map<int, double>::iterator iter = tb.find(id);
		if (iter == tb.end())
		{
			tb[id] = val;
		}
		else
		{
			iter->second += val;
		}

		OGRFeature::DestroyFeature(poFeature);
	}


	ShapeFile shpblk(filename, 1);

	ca11field = shpblk.getOrCreateField("ca11", OGRFieldType::OFTReal);

	shpblk.poLayer->ResetReading();
	id = 0;
	std::ofstream ofsCom((filename + ".csv").data());
	ofsCom << "id,ca11" << "\n";
	while ((poFeature = shpblk.poLayer->GetNextFeature()) != NULL)
	{
		std::map<int, double>::iterator iter = tb.find(id);
		if (iter != tb.end())
		{
			poFeature->SetField(ca11field, iter->second);
			ofsCom << id << "," << iter->second << "\n";
		}
		else
		{
			ofsCom << id << "," << 0 << "\n";
		}
		shpblk.poLayer->SetFeature(poFeature);
		id++;
		OGRFeature::DestroyFeature(poFeature);
	}
	shpblk.close();
	ofsCom.close();
}
void gridBaltimoreAggregrateByBlockGroup(std::string blockgroupfile, std::string basefile, std::string outdir)
{

	GDAL_DS<double>* baseds = new GDAL_DS<double>();
	baseds->open(basefile);
	std::vector<int> gridscales;
	for (int i = 1; i < 10; i++)
	{
		gridscales.push_back(i);
	}
	for (int i = 1; i <= 10; i++)
	{
		gridscales.push_back(i * 10);
	}
	char wkt[512];
	memcpy(wkt, baseds->projection.data(), baseds->projection.length());
	char* proj = wkt;
	OGRSpatialReference spatialref;
	spatialref.importFromWkt(&proj);
	double* data = baseds->readData(1);

	//std::ofstream ofs("B:/SpatialGranuality/Baltimore/blockgroup_aggregrates.csv");

	for (int iscale = 0; iscale < gridscales.size(); iscale++)
	{
		int scale = gridscales[iscale];
		std::stringstream ss;
		ss << outdir << scale * 10 << "/";
		if (!QDir(ss.str().data()).exists())
		{
			QDir(ss.str().data()).mkpath(".");
		}
		std::string pointsShpFile = ss.str() + "points.shp";
		std::string pointsBLKShpFile = ss.str() + "pointsBLK.shp";
		std::string BLKShpFile = ss.str() + "BLK.shp";


		if (!QFile(pointsShpFile.data()).exists())
		{
			ShapeFile pointsShp;
			pointsShp.create(pointsShpFile, &spatialref, 0, OGRwkbGeometryType::wkbPoint);
			int fid = pointsShp.getOrCreateField("ca11", OGRFieldType::OFTReal);
			int newcols = (int)(ceil(baseds->ncols / (double)scale));
			int newrows = (int)(ceil(baseds->nrows / (double)scale));


			double* pdata = data;
			double* aggdata = new double[newcols*newrows];
			memset(aggdata, 0, sizeof(double) * newcols*newrows);
			double total = 0;
			double val = 0;
			for (int oldrow = 0; oldrow < baseds->nrows; oldrow++)
			{
				int newrow = oldrow / scale;
				for (int oldcol = 0; oldcol < baseds->ncols; oldcol++)
				{
					int newcol = oldcol / scale;
					val = *pdata++;
					aggdata[newcol + newrow*newcols] += val;
					total += val;
				}
			}
			printf("%f\n", total);
			double resol = baseds->adfGeoTransform[1] * scale;
			double halfresol = resol * 0.5;
			double* paggdata = aggdata;
			total = 0;
			for (int newrow = 0; newrow < newrows; newrow++)
			{
				double y = baseds->adfGeoTransform[3] - resol * newrow - halfresol;
				for (int newcol = 0; newcol < newcols; newcol++)
				{
					double x = baseds->adfGeoTransform[0] + resol * newcol + halfresol;
					val = *paggdata++;
					//if (newrow == 0 && newcol == 249)
					//{
					// printf("");
					//}
					if (val > 0)
					{
						OGRPoint pt;
						pt.setX(x);
						pt.setY(y);
						OGRFeature* poFeature = OGRFeature::CreateFeature(pointsShp.poLayer->GetLayerDefn());
						poFeature->SetField(fid, val);
						poFeature->SetGeometry(&pt);
						pointsShp.poLayer->CreateFeature(poFeature);
						OGRFeature::DestroyFeature(poFeature);
						total += val;
					}

				}
			}

			delete[] aggdata;
		}


		if (!QFile(pointsBLKShpFile.data()).exists())
		{
			Preprocessor::intersectWithArcGIS(pointsShpFile, blockgroupfile, pointsBLKShpFile);
		}

		if (!QFile(BLKShpFile.data()).exists())
		{
			ShapeFile::copy(blockgroupfile, BLKShpFile);
			aggregrateSingleFile(pointsBLKShpFile, BLKShpFile);
		}
	}
	delete baseds;
	delete[] data;

}
void gridAllFolders(std::string indir, double gridsize)
{

	QDir rootdir(indir.data());
	rootdir.setFilter(QDir::Dirs | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	rootdir.setSorting(QDir::Name);
	std::vector<std::string> files;
	QFileInfoList list = rootdir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string dir = fileInfo.absoluteFilePath().toLocal8Bit().data();
		//Utils::updateFootprintForDir(dir,true);
		findShapes(dir, files);
	}

	std::string boundfile = indir + "bound.txt";
	OGREnvelope bound;
	if (QFileInfo(boundfile.data()).exists())
	{
		bound = BoundManager::readBound(boundfile);
	}
	else
	{
		std::vector<std::string> filesForBound;
		for (size_t i = 0; i < files.size(); i++)
		{
			QString qfilename = files[i].data();
			qfilename = qfilename.toLower();
			if (qfilename.endsWith("point.shp") || qfilename.endsWith("elecprod.shp") || qfilename.endsWith("airport.shp"))
				continue;
			filesForBound.push_back(files[i]);
		}
		bound = BoundManager::readBoundFromShapes(filesForBound);
		BoundManager::writeBound(bound, boundfile);
	}

	gridFiles(files, indir, bound, gridsize);

}
//void gridShapeFile(std::string infile, std::string outfile, std::string fishnetfile)
//{
//
//	std::string tmpefile = (QFileInfo(outfile.outlineDS()).dir().absolutePath() + "/" +  QFileInfo(outfile.outlineDS()).baseName() + "_2.Shp").toLocal8Bit().outlineDS();
//	std::string outfiletif = (QFileInfo(outfile.outlineDS()).dir().absolutePath() + "/" + QFileInfo(outfile.outlineDS()).baseName() + ".tif").toLocal8Bit().outlineDS();
//	Utils::updateFootprint(infile);
//	Preprocessor::intersectWithArcGIS(fishnetfile, infile, tmpefile);
//	//updateFieldAfterIntersection(tmpefile);
//	ShapeFile input(tmpefile);
//	Grid fishnet;
//	fishnet.fromFishnetShape(fishnetfile);
//
//	fishnet.gatherCells(&input,"C11");
//	fishnet.toShape(&input, outfile, true);
//	fishnet.toRaster(outfiletif);
//	const char *pszDriverName = "ESRI Shapefile";
//	OGRSFDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
//		pszDriverName);
//	poDriver->DeleteDataSource(tmpefile.outlineDS());
//}
//void gridShapeFile(std::string infile, std::string outfile, double gridsize)
//{
//	OGREnvelope bound = BoundManager::readBoundFromShape(infile.outlineDS());
//	Grid fishnet(bound, gridsize, 0);
//	fishnet.reset();
//	ShapeFile input(infile);
//	std::string tmpefishnetfile = (QFileInfo(outfile.outlineDS()).baseName() + "_1.Shp").toLocal8Bit().outlineDS();
//	std::string tmpefile = (QFileInfo(outfile.outlineDS()).baseName() + "_2.Shp").toLocal8Bit().outlineDS();
//	fishnet.toShape(&input, tmpefishnetfile, false);
//	Preprocessor::intersectWithArcGIS(tmpefishnetfile, infile, tmpefile);
//	input.close();
//	input = ShapeFile(tmpefile);
//	fishnet.gatherCells(&input);
//	fishnet.toShape(&input, outfile, true);
//	const char *pszDriverName = "ESRI Shapefile";
//	OGRSFDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
//		pszDriverName);
//	poDriver->DeleteDataSource(tmpefishnetfile.outlineDS());
//	poDriver->DeleteDataSource(tmpefile.outlineDS());
//}


void gridFolderForSectors(std::string indir, std::string outdir, double gridsize, bool skipNonRoad = false)
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
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
	}

	std::string nonroadfile = indir + "NonRoad.shp";
	std::string fishnetfile = outdir + "fishnet.shp";
	std::string rasterfile = outdir + "fishnet.tif";

	ShapeFile boundfile(nonroadfile.data());
	char wkt[512];
	char* pwkt = wkt;
	if (boundfile.poLayer->GetSpatialRef())
		boundfile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	if (!boundfile.poLayer->GetSpatialRef()->IsGeographic())
	{
		//FOOTPRINT_SCALE_FACTOR = 1;
	}
	else
	{
		//FOOTPRINT_SCALE_FACTOR = getDegreeToMeter(boundfile.poLayer);
		gridsize = gridsize / Utils::getDegreeToMeter(boundfile.poLayer);
	}


	Grid* fishnet = Grid::createFishnet(nonroadfile, gridsize);
	fishnet->reset();
	//fishnet->toShape(&boundfile, fishnetfile, false);

	for (size_t i = 0; i < files.size(); i++)
	{
		if (skipNonRoad && files[i] == "NonRoad.shp")
			continue;
		printf("%s\n", (outdir + files[i]).data());
		Preprocessor::intersectWithFishnet(fishnet, fishnetfile, indir + files[i], outdir + files[i]);
	}

	//fishnet->toShape(&boundfile, fishnetfile, true);
	//fishnet->toRaster(rasterfile);
	boundfile.close();

	for (size_t i = 0; i < files.size(); i++)
	{
		if (skipNonRoad && files[i] == "NonRoad.shp")
			continue;

		ShapeFile input((outdir + files[i]).data());
		if (input.poLayer->GetGeomType() == wkbPoint || input.poLayer->GetGeomType() == wkbMultiPoint || input.poLayer->GetGeomType() == wkbPoint25D)
			continue;
		std::string sectorfile = files[i].substr(0, files[i].length() - 4);
		sectorfile = outdir + sectorfile + ".tif";
		printf("%s\n", sectorfile.data());
		fishnet->reset();
		fishnet->gatherCells(&input);
		fishnet->toRaster(sectorfile, wkt);
	}

	delete fishnet;
}
//OGRErr OGRLayer::SetIgnoredFields(const char ** 	papszFields)
//virtual
//Set which fields can be omitted when retrieving features from the layer.
//
//If the driver supports this functionality(testable using OLCIgnoreFields capability), it will not fetch the specified fields in subsequent calls to GetFeature() / GetNextFeature() and thus save some processing time and / Orange bandwidth.
//
//Besides field names of the layers, the following special fields can be passed : "OGR_GEOMETRY" to ignore geometry and "OGR_STYLE" to ignore layer style.
//
//By default, no fields are ignored.
//
//This method is the same as the C function OGR_L_SetIgnoredFields()
//
//Parameters
//papszFields	an array of field names terminated by NULL item.If NULL is passed, the ignored list is cleared.
//Returns
//OGRERR_NONE if all field names have been resolved(even if the driver does not support this method)
//Reimplemented in GNMGenericLayer, OGRProxiedLayer, OGRUnionLayer, OGRMutexedLayer, and OGRLayerDecorator.

std::vector<std::string> loadAllSectorFile(std::string indir)
{
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	std::vector<std::string> files;

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().toLower().endsWith("fishnet.shp"))
			continue;

		files.push_back(indir + fileInfo.fileName().toLocal8Bit().data());
	}
	return files;
}

//std::vector<std::string> findChildDirctories(std::string lidarFisheyeDir, OGREnvelope bound, double gridsize, std::string resolution)
//{
//	QDir qoutdir(lidarFisheyeDir.outlineDS());
//	if (!qoutdir.exists())
//		qoutdir.mkpath(".");
//	lidarFisheyeDir = (qoutdir.absolutePath() + "/").toLocal8Bit().outlineDS();
//	std::vector<std::string> dirs;
//	Grid grid(bound, gridsize, 1);
//	for (int n = 0; n < grid.nrows*grid.ncols; n++)
//	{
//		std::stringstream subdirss;
//		subdirss << lidarFisheyeDir << n;
//		OGREnvelope bound = BoundManager::readBoundFromShape(subdirss.str() + ".Shp");
//
//		subdirss << "/" << resolution;
//		Utils::updateFishnet(subdirss.str().outlineDS(), bound, 32.8084);
//
//		dirs.push_back(subdirss.str().outlineDS());
//	}
//	return dirs;
//}

void processFolderForGrids(std::string indir, std::string outdir, std::string gridoutdir, double gridsize)
{
	QDir qoutdir(outdir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();


	QDir qgridoutdir(gridoutdir.data());
	if (!qgridoutdir.exists())
		qgridoutdir.mkpath(".");
	gridoutdir = (qgridoutdir.absolutePath() + "/").toLocal8Bit().data();

	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
	}

	std::string nonroadfile = indir + "NonRoad.shp";
	std::string fishnetfile = outdir + "fishnet.shp";
	std::string rasterfile = outdir + "fishnet.tif";

	ShapeFile boundfile(nonroadfile.data());
	if (!boundfile.poLayer->GetSpatialRef()->IsGeographic())
	{
		//FOOTPRINT_SCALE_FACTOR = 1;
	}
	else
	{
		//FOOTPRINT_SCALE_FACTOR = getDegreeToMeter(boundfile.poLayer);
		gridsize = gridsize / Utils::getDegreeToMeter(boundfile.poLayer);
	}

	Grid* fishnet = Grid::createFishnet(nonroadfile, gridsize);

	/*fishnet->toShape(&boundfile, fishnetfile, false);

	for (size_t i = 0; i < files.size(); i++)
	{

	printf("%s\n", (lidarFisheyeDir + files[i]).outlineDS());
	intersectWithArcGIS(fishnetfile, indir + files[i], lidarFisheyeDir + files[i]);
	updateFieldAfterIntersection(lidarFisheyeDir + files[i]);
	ShapeFile input((lidarFisheyeDir + files[i]).outlineDS());
	gatherCells(fishnet, &input);
	}

	fishnet->toShape(&boundfile, fishnetfile, true);
	fishnet->toRaster(rasterfile);*/



	fishnet->reset();
	for (size_t i = 0; i < files.size(); i++)
	{
		ShapeFile input((outdir + files[i]).data());
		if (input.poLayer->GetGeomType() == wkbPoint || input.poLayer->GetGeomType() == wkbMultiPoint || input.poLayer->GetGeomType() == wkbPoint25D)
			continue;
		std::string sectorfile = files[i].substr(0, files[i].length() - 4);
		std::string sectorfileout = gridoutdir + sectorfile + ".tif";
		printf("%s\n", sectorfile.data());
		fishnet->gatherCells(&input);

	}
	fishnet->normalizedByArea();


	char wkt[512];
	char* pwkt = wkt;
	if (boundfile.poLayer->GetSpatialRef())
		boundfile.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	fishnet->toRaster(gridoutdir + "AllSectors.tif", wkt);

	fishnet->toShape(wkt, gridoutdir + "AllSectors.shp", true);

	for (size_t i = 0; i < files.size(); i++)
	{
		ShapeFile input((outdir + files[i]).data());
		if (input.poLayer->GetGeomType() == wkbPoint || input.poLayer->GetGeomType() == wkbMultiPoint || input.poLayer->GetGeomType() == wkbPoint25D)
			continue;
		std::string sectorfile = files[i].substr(0, files[i].length() - 4);
		std::string sectorfileout = gridoutdir + sectorfile + ".tif";
		printf("%s\n", sectorfile.data());
		fishnet->reset();
		fishnet->gatherCells(&input);
		fishnet->toRaster(sectorfileout, wkt);
		fishnet->normalizedByArea();
		sectorfileout = gridoutdir + sectorfile + ".shp";
		fishnet->toShape(wkt, sectorfileout, true);
	}

	fishnet->reset();
	for (size_t i = 0; i < files.size(); i++)
	{
		ShapeFile input((outdir + files[i]).data());
		if (input.poLayer->GetGeomType() == wkbPoint || input.poLayer->GetGeomType() == wkbMultiPoint || input.poLayer->GetGeomType() == wkbPoint25D)
			continue;
		std::string sectorfile = files[i].substr(0, files[i].length() - 4);
		if (!QString(sectorfile.data()).contains("NonPoint"))
			continue;
		fishnet->gatherCells(&input);

	}
	fishnet->normalizedByArea();
	fishnet->toRaster(gridoutdir + "NonPoint.tif", wkt);

	fishnet->toShape(wkt, gridoutdir + "NonPoint.shp", true);
	boundfile.close();

	delete fishnet;
}

std::string double2string(double d)
{
	std::stringstream ss;
	ss << (long)d;
	return ss.str();
}

void createBySector(std::string filename, std::string dir, int fips)
{
	struct PointRecord
	{
		/*2 IND
		3 COM
		6 NonRoad
		7 Railroad
		9 Airport*/
		std::string scc;
		double ca;
		std::string facilityid;
		std::string name;
		int sector;
		double lon;
		double lat;
		int FIPS;
		PointRecord(OGRFeature* fea)
		{
			//name	lat	lon	scc	facilityid	source	sector	ca
			scc = double2string(fea->GetFieldAsDouble("scc"));

			OGRPoint*  p = (OGRPoint*)fea->GetGeometryRef();
			lon = p->getX();
			lat = p->getY();
			ca = fea->GetFieldAsDouble("SR_CO2_tC");

			facilityid = double2string(fea->GetFieldAsDouble("facility_i"));

			sector = (int)fea->GetFieldAsDouble("Sector_ID");
			FIPS = (int)fea->GetFieldAsDouble("FIPS");

			name = fea->GetFieldAsString("facility_n");

		}

		void setFields(OGRFeature* poFeature)
		{
			poFeature->SetField("name", name.data());
			poFeature->SetField("scc", scc.data());
			poFeature->SetField("ca", ca);
			poFeature->SetField("facilityid", facilityid.data());

			OGRPoint pt;
			pt.setX(lon);
			pt.setY(lat);
			poFeature->SetGeometry(&pt);
		}

	};
	QFileInfo fileInfo(filename.data());
	std::vector<ShapeFileBySector> shpBySectors;


	//std::string dir = "Z:/Hestia_Baltimore/gridPrep_SHP_master/";
	ShapeFile IndPoint;
	IndPoint.create((dir + "IndPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&IndPoint, "IND", 2));

	ShapeFile ComPoint;
	ComPoint.create((dir + "ComPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&ComPoint, "COM", 3));

	ShapeFile NonRoadPoint;
	NonRoadPoint.create((dir + "NonRoadPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&NonRoadPoint, "NonRoad", 6));

	ShapeFile createRailroadPoint;
	createRailroadPoint.create((dir + "RailroadPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&createRailroadPoint, "RAIL", 7));

	ShapeFile Airport;
	Airport.create((dir + "Airport.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	shpBySectors.push_back(ShapeFileBySector(&Airport, "AIR", 9));



	for (size_t i = 0; i < shpBySectors.size(); i++)
	{
		OGRFieldDefn fieldname("name", OGRFieldType::OFTString);
		OGRFieldDefn fieldscc("scc", OGRFieldType::OFTString);
		OGRFieldDefn fieldfacilityid("facilityid", OGRFieldType::OFTString);
		OGRFieldDefn fieldca("ca", OGRFieldType::OFTReal);
		shpBySectors[i].SHP->poLayer->CreateField(&fieldname);
		shpBySectors[i].SHP->poLayer->CreateField(&fieldfacilityid);
		shpBySectors[i].SHP->poLayer->CreateField(&fieldscc);
		shpBySectors[i].SHP->poLayer->CreateField(&fieldca);


	}

	ShapeFile shpscc(filename);
	OGRFeature *poFeature;
	shpscc.poLayer->ResetReading();

	while ((poFeature = shpscc.poLayer->GetNextFeature()) != NULL)
	{
		PointRecord record(poFeature);
		if (record.FIPS != fips)
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		for (size_t i = 0; i < shpBySectors.size(); i++)
		{
			ShapeFileBySector& shpBySector = shpBySectors[i];
			if (shpBySector.SecID != record.sector)
				continue;
			OGRFeature* poFeatureNew = OGRFeature::CreateFeature(shpBySector.SHP->poLayer->GetLayerDefn());
			record.setFields(poFeatureNew);
			shpBySector.SHP->poLayer->CreateFeature(poFeatureNew);
			OGRFeature::DestroyFeature(poFeatureNew);
		}

		OGRFeature::DestroyFeature(poFeature);
	}

}

void linkSCC(std::string shapefile, std::string sccFile)
{
	struct SCCRecord
	{
		std::string scc;
		double lon;
		double lat;
		double ca;
		std::string facilityid;
		SCCRecord(OGRFeature* fea)
		{
			scc = double2string(fea->GetFieldAsDouble("scc"));
			lon = fea->GetFieldAsDouble("New_LON");
			lat = fea->GetFieldAsDouble("New_LAT");
			ca = fea->GetFieldAsDouble("SR_CO2_tC");
			facilityid = double2string(fea->GetFieldAsDouble("facility_i"));
		}
	};
	std::vector<SCCRecord> records;
	ShapeFile shpscc(sccFile);
	OGRFeature *poFeature;
	shpscc.poLayer->ResetReading();

	std::map<std::string, std::vector<SCCRecord>> recordsmap;
	while ((poFeature = shpscc.poLayer->GetNextFeature()) != NULL)
	{
		SCCRecord record(poFeature);
		records.push_back(record);
		std::map<std::string, std::vector<SCCRecord>>::iterator iter = recordsmap.find(record.facilityid);
		if (iter == recordsmap.end())
		{
			std::vector<SCCRecord> group;
			group.push_back(record);
			recordsmap[record.facilityid] = group;
		}
		else
		{
			iter->second.push_back(record);
		}
		OGRFeature::DestroyFeature(poFeature);
	}


	ShapeFile shp(shapefile, 1);
	OGRFieldDefn field("scc", OGRFieldType::OFTString);
	shp.poLayer->CreateField(&field);
	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		std::string id = double2string(poFeature->GetFieldAsDouble("facilityid"));
		std::map<std::string, std::vector<SCCRecord>>::iterator iter = recordsmap.find(id);
		if (iter != recordsmap.end())
		{
			std::vector<SCCRecord>& group = iter->second;
			double ca = poFeature->GetFieldAsDouble("ca11");
			bool found = false;
			for (int i = 0; i < group.size(); i++)
			{
				SCCRecord& rc = group[i];
				if ((rc.ca - ca) < 0.00001)
				{
					poFeature->SetField("scc", rc.scc.data());

					found = true;
				}
			}
			if (!found)
			{
				printf("%s,%f", id, ca);
			}
		}
		else
		{
			poFeature->SetField("scc", "-1");
		}
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
}


//	/*2 IND
//	3 COM
//	6 NonRoad
//	7 Railroad
//	9 Airport*/
//	int SectorID;
//	double longitude;
//	double latitude;
//	int FIPS;
//	std::string facility_name;
//	std::string scc;
//	double New_LAT;
//	double New_LON;
//	double COMNG;
//	double COMpetrol;
//	double INDNG;
//	double INDpetrol;
//	double ELE;
//	double AIR;
//	double RAIL;
//	double Nonroad;

//sector BlockGroupID
//longitude	latitude
//facility_name scc New_LAT New_LON
//COM(tC / yr)	
//COM NG	
//COM petrol	
//IND(tC / yr)	
//IND NG	
//IND petrol	
//ELE(tC / yr)	
//ELE NG	ELE petrol	
//AIR(tC / yr)	
//RAIL(tC / yr)	
//Nonroad(tC / yr)

typedef std::pair<std::string, std::string> FieldField;

void createSectorShape(std::string csvfile, std::string outfile, ShapeCreator::SECTORID sector, std::vector<std::pair<std::string, std::string>> fuelfields)
{
	//ValueField(OGRFieldType otp, std::string destfieldname, std::string firstfieldname, std::string secondfieldname = "")
	ValueField  filterfield(OGRFieldType::OFTInteger, "", "Sector ID");
	filterfield.i = (int)sector;
	std::vector<ValueField> fields;

	ValueField  lonfield(OGRFieldType::OFTReal, "", "New_LON", "longitude");
	ValueField  latfield(OGRFieldType::OFTReal, "", "New_LAT", "latitude");

	fields.push_back(ValueField(OGRFieldType::OFTString, "name", "facility_name"));
	fields.push_back(ValueField(OGRFieldType::OFTString, "facilityid", "facility_id"));
	fields.push_back(ValueField(OGRFieldType::OFTString, "scc", "scc"));
	for (size_t i = 0; i < fuelfields.size(); i++)
	{
		fields.push_back(ValueField(OGRFieldType::OFTReal, fuelfields[i].first, fuelfields[i].second));
	}

	ShapeCreator shpCreator;
	shpCreator.create(csvfile, outfile, fields, lonfield, latfield, filterfield);
}
//"COM (tC/yr)	COM NG	COM petrol	IND (tC/yr)	IND NG	IND petrol	ELE (tC/yr)	ELE NG	ELE petrol	AIR (tC/yr)	RAIL (tC/yr)	Nonroad (tC/yr)"
std::vector<ValueField> getFieldsForSector(ShapeCreator::SECTORID sector)
{
	std::vector<ValueField> fuelfields;
	if (sector == ShapeCreator::IND)
	{
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca_ng", "IND NG"));
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca_p", "IND petrol"));
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca", "IND (tC/yr)"));
	}
	else if (sector == ShapeCreator::COM)
	{
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca_ng", "COM NG"));
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca_p", "COM petrol"));
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca", "COM (tC/yr)"));
	}
	else if (sector == ShapeCreator::ELE)
	{
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca_ng", "ELE NG"));
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca_p", "ELE petrol"));
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca", "ELE (tC/yr)"));
	}
	else if (sector == ShapeCreator::NONROAD)
	{
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca", "Nonroad (tC/yr)"));
	}
	else if (sector == ShapeCreator::AIRPORT)
	{
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca", "AIR (tC/yr)"));
	}
	else if (sector == ShapeCreator::RAILROAD)
	{
		fuelfields.push_back(ValueField(OGRFieldType::OFTReal, "ca", "RAIL (tC/yr)"));
	}
	return fuelfields;
}

std::string getFileNameForSector(ShapeCreator::SECTORID sector)
{
	std::vector<std::pair<std::string, std::string>> fuelfields;
	if (sector == ShapeCreator::IND)
	{
		return "IndPoint";
	}
	else if (sector == ShapeCreator::COM)
	{
		return "ComPoint";
	}
	else if (sector == ShapeCreator::ELE)
	{
		return "ElecProd";
	}
	else if (sector == ShapeCreator::NONROAD)
	{
		return "NonRoadPoint";
	}
	else if (sector == ShapeCreator::AIRPORT)
	{
		return "Airport";
	}
	else if (sector == ShapeCreator::RAILROAD)
	{
		return "RailroadPoint";
	}
	return "";
}
void createShapesForSectors(std::string csvfile, std::string outdir, std::vector<ShapeCreator::SECTORID> sectors)
{
	ShapeCreator shpCreator;
	ValueField  lonfield(OGRFieldType::OFTReal, "", "New_LON", "longitude");
	ValueField  latfield(OGRFieldType::OFTReal, "", "New_LAT", "latitude");

	for (size_t i = 0; i < sectors.size(); i++)
	{
		ShapeCreator::SECTORID sector = sectors[i];
		ValueField  filterfield(OGRFieldType::OFTInteger, "", "Sector ID");
		filterfield.i = (int)sector;
		std::string outfilename = outdir + getFileNameForSector(sector) + ".shp";
		std::vector<ValueField> fields = getFieldsForSector(sector);

		fields.push_back(ValueField(OGRFieldType::OFTString, "name", "facility_name"));
		fields.push_back(ValueField(OGRFieldType::OFTString, "facilityid", "facility_id"));
		fields.push_back(ValueField(OGRFieldType::OFTString, "scc", "scc"));

		shpCreator.create(csvfile, outfilename, fields, lonfield, latfield, filterfield);
	}
}
//struct SpatialAllocationParams
//{
//	NonPointProcessor::NonPointFuel fuel;
//	int year;
//	int sector;
//	int division;
//	double total;
//	std::string output_field;
//
//};
void calculateNonPointWeight(std::string dir, std::string rscript, std::string sharedCFG)
{
	//std::string sectorname[] = {"ComNonPoint","ResNonPoint","IndNonPoint"};
	NonPointProcessor::NonPointFuel fueltypes[] = { NonPointProcessor::NG,NonPointProcessor::PETROL,NonPointProcessor::COAL };
	std::string fuelfields[] = { "ca_ng", "ca_p","ca_c" };
	std::string sectors[] = { "ComNonPoint","ResNonPoint","IndNonPoint" };
	NonPointProcessor processor;
	for (size_t i = 0; i < 3; i++)
	{
		std::string shapefilename = dir + sectors[i] + ".dbf";
		NonPointProcessor::SpatialAllocationParams ng(NonPointProcessor::NG, 2020, i + 1, 9, 1000000, "ca_ng");
		NonPointProcessor::SpatialAllocationParams petrol(NonPointProcessor::PETROL, 2020, i + 1, 9, 1000000, "ca_p");
		NonPointProcessor::SpatialAllocationParams coal(NonPointProcessor::COAL, 2020, i + 1, 9, 1000000, "ca_c");
		std::vector<NonPointProcessor::SpatialAllocationParams> nonpointparams;
		nonpointparams.push_back(ng);
		nonpointparams.push_back(petrol);
		nonpointparams.push_back(coal);
		processor.runRScript(rscript, sharedCFG, shapefilename, nonpointparams);
	}


}

void calNonRoadTotals(std::vector<std::string> dirs, std::string outfile)
{
	std::ofstream ofs;
	ofs.open(outfile.data());
	for (size_t i = 0; i < dirs.size(); i++)
	{
		double sum = ShapeFile::getTotal(dirs[i] + "NonRoadPoint.shp", "ca");
		ofs << dirs[i] + "NonRoadPoint.shp" << "," << sum << "\n";
	}
	ofs.close();
}
void caTotals(std::vector<std::string> dirs, std::string field, std::string shapename, std::string outfile)
{
	std::ofstream ofs;
	ofs.open(outfile.data());
	for (size_t i = 0; i < dirs.size(); i++)
	{
		double sum = ShapeFile::getTotal(dirs[i] + shapename, field);
		ofs << dirs[i] + shapename << "," << sum << "\n";
	}
	ofs.close();
}

void copyField(std::vector<std::string> dirs, std::string oldfield, std::string newfield)
{

	for (size_t i = 0; i < dirs.size(); i++)
	{
		ShapeFile::copyField(dirs[i] + "NonRoad.shp", oldfield, newfield);
	}

}
void copyElecProdTimeStruct(std::string destFile, std::string srcFile)
{
	//ShapeFile ElecProd((dir + "ElecProd.Shp").outlineDS(), 1);

}


std::string int2string(int val)
{
	std::stringstream ss;
	ss << val;
	return ss.str();
}

void createColumn(std::string dir)
{

	ShapeFile IndPoint((dir + "IndPoint.shp").data(), 1);
	ShapeFile ComPoint((dir + "ComPoint.shp").data(), 1);
	ShapeFile ElecProd((dir + "ElecProd.shp").data(), 1);
	//ShapeFile Airport((dir + "Airport.Shp").outlineDS(), 1);
	ShapeFile createRailroadPoint((dir + "RailroadPoint.shp").data(), 1);
	ShapeFile NonRoadPoint((dir + "NonRoadPoint.shp").data(), 1);
	ShapeFile ResNonPoint((dir + "ResNonPoint.shp").data(), 1);
	ShapeFile ComNonPoint((dir + "ComNonPoint.shp").data(), 1);
	ShapeFile IndNonPoint((dir + "IndNonPoint.shp").data(), 1);
	ShapeFile Railroad((dir + "Railroad.shp").data(), 1);
	ShapeFile OnRoad((dir + "OnRoad.shp").data(), 1);
	//ShapeFile CMVPort((dir + "Port_Polygons.Shp").outlineDS(), 1);
	//ShapeFile CMVUnderway((dir + "Port_Lanes.Shp").outlineDS(), 1);
	//ShapeFile NonRoad((dir + "NonRoad.Shp").outlineDS(), 1);


	//std::string smokeDir = "B:/NonRoadTemporal/NonRoadTemporal/";
	//std::string sccOutDir = "B:/LA_Version2/Time/SCC/2011/";
	int idx = IndPoint.getOrCreateField("timestruct", OFTString);
	OGRFeatureDefn* featureDef = IndPoint.poLayer->GetLayerDefn();
	OGRFeature* poFeature;
	while ((poFeature = IndPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string scc = poFeature->GetFieldAsString("scc");
		if (copy_SCC_PROFILE(scc, "9" + scc))
		{

			poFeature->SetField(idx, ("9" + scc).data());
		}
		else
		{
			poFeature->SetField(idx, "0");
		}
		IndPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

	idx = ComPoint.getOrCreateField("timestruct", OFTString);
	featureDef = ComPoint.poLayer->GetLayerDefn();
	while ((poFeature = ComPoint.poLayer->GetNextFeature()) != NULL)
	{

		std::string scc = poFeature->GetFieldAsString("scc");
		if (copy_SCC_PROFILE(scc, "9" + scc))
		{
			poFeature->SetField(idx, ("9" + scc).data());
		}
		else
		{
			poFeature->SetField(idx, "0");
		}
		ComPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

	idx = ElecProd.getOrCreateField("timestruct", OFTString);
	featureDef = ElecProd.poLayer->GetLayerDefn();
	while ((poFeature = ElecProd.poLayer->GetNextFeature()) != NULL)
	{

		std::string scc = poFeature->GetFieldAsString("scc");
		std::string timestructElec = "";
		if (scc != "" && atoi(scc.data()) > -1)
		{

			timestructElec = scc;
			if (copy_SCC_PROFILE(scc, "9" + scc))
			{
				//printf("copy smoke");				
				poFeature->SetField(idx, ("9" + scc).data());
			}
			else
			{
				poFeature->SetField(idx, "0");
			}
		}
		else
		{
			std::string facilityid = poFeature->GetFieldAsString("facilityid");
			poFeature->SetField(idx, ("4" + facilityid).data());
			timestructElec = facilityid;
		}

		//smokeDir
		//else
		//{
		//	poFeature->SetField(idx, "0");
		//}


		ElecProd.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

	//idx = Airport.getOrCreateField("timestruct", OFTString);
	//featureDef = Airport.poLayer->GetLayerDefn();
	//while ((poFeature = Airport.poLayer->GetNextFeature()) != NULL)
	//{
	//	Airport.poLayer->SetFeature(poFeature);
	//	poFeature->SetField(idx, "0");
	//	Airport.poLayer->SetFeature(poFeature);
	//	OGRFeature::DestroyFeature(poFeature);
	//}

	idx = createRailroadPoint.getOrCreateField("timestruct", OFTString);
	featureDef = createRailroadPoint.poLayer->GetLayerDefn();
	while ((poFeature = createRailroadPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string scc = poFeature->GetFieldAsString("scc");
		if (copy_SCC_PROFILE(scc, "9" + scc))
		{
			poFeature->SetField(idx, ("9" + scc).data());
		}
		else
		{
			poFeature->SetField(idx, "0");
		}
		createRailroadPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

	idx = NonRoadPoint.getOrCreateField("timestruct", OFTString);
	featureDef = NonRoadPoint.poLayer->GetLayerDefn();
	while ((poFeature = NonRoadPoint.poLayer->GetNextFeature()) != NULL)
	{

		std::string scc = poFeature->GetFieldAsString("scc");
		if (copy_SCC_PROFILE(scc, "9" + scc))
		{
			poFeature->SetField(idx, ("9" + scc).data());
		}
		else
		{
			poFeature->SetField(idx, "0");
		}
		NonRoadPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}




	idx = ResNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = ResNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = ResNonPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string bt = int2string(poFeature->GetFieldAsInteger("bt"));
		if (bt == "" || atof(bt.data()) <= 0)
		{
			poFeature->SetField(idx, "0");
		}
		else
		{
			poFeature->SetField(idx, ("2" + bt).data());
		}
		ResNonPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}


	idx = ComNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = ComNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = ComNonPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string bt = int2string(poFeature->GetFieldAsInteger("bt"));
		if (bt == "" || atof(bt.data()) <= 0)
		{
			poFeature->SetField(idx, "0");
		}
		else
		{
			poFeature->SetField(idx, ("1" + bt).data());
		}
		ComNonPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}


	idx = IndNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = IndNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = IndNonPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string bt = int2string(poFeature->GetFieldAsInteger("bt"));
		if (bt == "" || atof(bt.data()) <= 0)
		{
			poFeature->SetField(idx, "0");
		}
		else
		{
			poFeature->SetField(idx, ("3" + bt).data());
		}
		IndNonPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

	idx = Railroad.getOrCreateField("timestruct", OFTString);
	featureDef = Railroad.poLayer->GetLayerDefn();
	while ((poFeature = Railroad.poLayer->GetNextFeature()) != NULL)
	{
		poFeature->SetField(idx, "0");
		Railroad.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}


	//idx = OnRoad.getOrCreateField("timestruct", OFTString);
	//featureDef = OnRoad.poLayer->GetLayerDefn();
	//int roadID = 0;
	//while ((poFeature = OnRoad.poLayer->GetNextFeature()) != NULL)
	//{
	//	std::stringstream ss;
	//	roadID++;
	//	ss << "600" << roadID;
	//	poFeature->SetField(idx, ss.str().outlineDS());
	//	OnRoad.poLayer->SetFeature(poFeature);
	//	OGRFeature::DestroyFeature(poFeature);
	//}


	/*idx = CMVPort.getOrCreateField("timestruct", OFTString);
	featureDef = CMVPort.poLayer->GetLayerDefn();
	while ((poFeature = CMVPort.poLayer->GetNextFeature()) != NULL)
	{
	poFeature->SetField(idx, "0");
	CMVPort.poLayer->SetFeature(poFeature);
	OGRFeature::DestroyFeature(poFeature);
	}

	idx = CMVUnderway.getOrCreateField("timestruct", OFTString);
	featureDef = CMVUnderway.poLayer->GetLayerDefn();
	while ((poFeature = CMVUnderway.poLayer->GetNextFeature()) != NULL)
	{
	poFeature->SetField(idx, "0");
	CMVUnderway.poLayer->SetFeature(poFeature);
	OGRFeature::DestroyFeature(poFeature);
	}
	*/
	//idx = NonRoad.getOrCreateField("timestruct", OFTString);
	//featureDef = NonRoad.poLayer->GetLayerDefn();
	//while ((poFeature = NonRoad.poLayer->GetNextFeature()) != NULL)
	//{
	//	poFeature->SetField(idx, "0");
	//	NonRoad.poLayer->SetFeature(poFeature);
	//	OGRFeature::DestroyFeature(poFeature);
	//}

}

void updateNonpointTime(std::string dir)
{

	/*ShapeFile ResNonPoint((dir + "ResNonPoint.Shp").outlineDS(), 1);
	ShapeFile ComNonPoint((dir + "ComNonPoint.Shp").outlineDS(), 1);*/
	ShapeFile IndNonPoint((dir + "IndNonPoint.shp").data(), 1);


	int idx;
	OGRFeatureDefn* featureDef;
	OGRFeature* poFeature;

	/*idx = ResNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = ResNonPoint.poLayer->GetLayerDefn();

	while ((poFeature = ResNonPoint.poLayer->GetNextFeature()) != NULL)
	{
	std::string ts = poFeature->GetFieldAsString(idx);
	ts = ts.substr(1, ts.size() - 1);
	poFeature->SetField(idx, ts.outlineDS());
	ResNonPoint.poLayer->SetFeature(poFeature);
	OGRFeature::DestroyFeature(poFeature);
	}


	idx = ComNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = ComNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = ComNonPoint.poLayer->GetNextFeature()) != NULL)
	{
	std::string ts = poFeature->GetFieldAsString(idx);
	ts = ts.substr(1, ts.size() - 1);
	poFeature->SetField(idx, ts.outlineDS());
	ComNonPoint.poLayer->SetFeature(poFeature);
	OGRFeature::DestroyFeature(poFeature);
	}

	*/
	idx = IndNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = IndNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = IndNonPoint.poLayer->GetNextFeature()) != NULL)
	{

		poFeature->SetField(idx, "3000");
		IndNonPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}
void copyFolderDropFields(std::vector<std::string> dirs)
{
	std::vector<std::string> fields2keep;
	fields2keep.push_back("ca11");
	for (size_t i = 0; i < dirs.size(); i++)
	{
		std::string foldername = QDir(dirs[i].data()).dirName().toLocal8Bit().data();
		ShapeFile::copyDir(dirs[i], "B:/LA_Version2/ca11/" + foldername + "/", fields2keep);
		//ShapeFile::scaleCAFields(dirs[i]);
	}

}
void copy2PublicFolder(std::vector<std::string> dirs, std::string outdir)
{
	std::vector<std::string> fields2keep;
	fields2keep.push_back("name"); fields2keep.push_back("timestruct"); fields2keep.push_back("tfs"); fields2keep.push_back("yearbuilt");
	fields2keep.push_back("ca11"); fields2keep.push_back("ca11_ng"); fields2keep.push_back("ca11_p"); fields2keep.push_back("ca11_c");
	fields2keep.push_back("ca12"); fields2keep.push_back("ca12_ng"); fields2keep.push_back("ca12_p"); fields2keep.push_back("ca12_c");
	fields2keep.push_back("ca13"); fields2keep.push_back("ca13_ng"); fields2keep.push_back("ca13_p"); fields2keep.push_back("ca13_c");
	fields2keep.push_back("ca14"); fields2keep.push_back("ca14_ng"); fields2keep.push_back("ca14_p"); fields2keep.push_back("ca14_c");

	for (size_t i = 0; i < dirs.size(); i++)
	{
		std::string foldername = QDir(dirs[i].data()).dirName().toLocal8Bit().data();
		ShapeFile::copyDir(dirs[i], outdir + foldername + "/", fields2keep);
	}

}

void processDirs(std::vector<std::string> dirs)
{
	std::string destwktfile = "B:/LA_Version2/gridPrep_SHP_master/NAD_1983_StatePlane_California_V_FIPS_0405_Feet.txt";
	std::vector<ShapeCreator::SECTORID> sectors;
	sectors.push_back(ShapeCreator::IND); sectors.push_back(ShapeCreator::COM); sectors.push_back(ShapeCreator::RAILROAD); sectors.push_back(ShapeCreator::AIRPORT); sectors.push_back(ShapeCreator::NONROAD);
	for (size_t i = 0; i < dirs.size(); i++)
	{
		std::string dir = dirs[i];

		//Utils::updateFootprintForDir(dir);
		//ShapeFile::copyField(dirs[i] + "OnRoad.Shp", "tC2011", "ca");
		//createShapesForSectors(dir + "NEI.csv", dir, sectors);
		ScaleByFuel scalebyfuel;
		scalebyfuel.scale(dir + "Scaling.csv", dir);
		//std::string foldername = QDir(dirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
		//ShapeFile::reprojectDir(dirs[i], dirs[i] + "../../reprojected/" + foldername + "/", destwktfile);
		//gridFolderByRaster(dir, dir + "grid/", 500 * 3.28084);
		//Utils::updateFootprintForDir(dir);

	}

}
void gridDirs(std::vector<std::string> dirs)
{
	std::string destwktfile = "B:/LA_Version2/gridPrep_SHP_master/NAD_1983_StatePlane_California_V_FIPS_0405_Feet.txt";
	std::vector<ShapeCreator::SECTORID> sectors;
	sectors.push_back(ShapeCreator::IND); sectors.push_back(ShapeCreator::COM); sectors.push_back(ShapeCreator::RAILROAD); sectors.push_back(ShapeCreator::AIRPORT); sectors.push_back(ShapeCreator::NONROAD);
	for (size_t i = 0; i < dirs.size(); i++)
	{
		std::string dir = dirs[i];

		//ShapeFile::copyField(dirs[i] + "OnRoad.Shp", "tC2011", "ca");
		//createShapesForSectors(dir + "NEI.csv", dir, sectors);
		//ScaleByFuel scalebyfuel;
		//scalebyfuel.scale(dir + "Scaling.csv", dir);
		//std::string foldername = QDir(dirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
		//ShapeFile::reprojectDir(dirs[i], dirs[i] + "../reprojected/" + foldername + "/", destwktfile);
		gridFolderByRaster(dir, dir + "grid/", 500 * 3.28084);

	}

}

void createTimeStruct(std::vector<std::string> dirs)
{
	std::vector<std::string> fields2keep;
	fields2keep.push_back("ca11");
	for (size_t i = 0; i < dirs.size(); i++)
	{
		createColumn(dirs[i]);
	}

}


void createTimeStructForNonRoad(std::vector<std::string> dirs, int fips[])
{
	for (size_t idir = 0; idir < dirs.size(); idir++)
	{
		std::vector<std::string> allfilesUnderDir = Utils::findFiles(dirs[idir], ".shp");
		std::vector<std::string> nonroadfiles;
		int fip = fips[idir];
		for (size_t i = 0; i < allfilesUnderDir.size(); i++)
		{
			std::string filename = allfilesUnderDir[i];
			QFileInfo fileinfo(filename.data());
			QString name = fileinfo.baseName();
			if (!name.startsWith("NonRoad", Qt::CaseSensitivity::CaseInsensitive))
				continue;
			if (!name.contains("_", Qt::CaseSensitivity::CaseInsensitive))
				continue;
			nonroadfiles.push_back(filename);
			ShapeFile shp(filename, 1);
			int tsfield = shp.getOrCreateField("timestruct", OGRFieldType::OFTString);
			int surrogateCode = atoi(std::string(name.toLocal8Bit().data()).substr(8, 3).data());
			std::stringstream ss;
			ss << "6" << fip << surrogateCode;
			std::string ts = ss.str();

			OGRFeature *poFeature;
			double sum = 0;
			shp.poLayer->ResetReading();

			while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
			{
				poFeature->SetField(tsfield, ts.data());
				shp.poLayer->SetFeature(poFeature);
				OGRFeature::DestroyFeature(poFeature);
			}
		}
	}

}
void createTimeStructForAirport(std::vector<std::string> dirs)
{
	for (size_t idir = 0; idir < dirs.size(); idir++)
	{
		std::vector<std::string> allfilesUnderDir = Utils::findFiles(dirs[idir], ".shp");
		std::string tsdir = "B:/LA_Version2/Time/Airport/output/";
		std::vector<std::string> tsfiles = Utils::findFiles(tsdir, ".txt");

		struct  TSFileInfo
		{
			std::string name;
			std::string path;
			int idx;
		};
		std::map<std::string, TSFileInfo> tsdic;
		for (size_t i = 0; i < tsfiles.size(); i++)
		{
			QFileInfo fi(tsfiles[i].data());
			TSFileInfo tfi;
			tfi.name = fi.completeBaseName().toLower().toLocal8Bit().data();
			tfi.path = fi.absoluteFilePath().toLocal8Bit().data();
			tfi.idx = i;
			tsdic[tfi.name] = tfi;
		}

		std::ifstream ifs_class;
		ifs_class.open("B:/LA_Version2/Time/Airport/AircraftClass.csv");
		std::string line;
		std::getline(ifs_class, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		std::map<std::string, std::string> classdic;
		while (ifs_class.peek() != -1)
		{
			std::getline(ifs_class, line);
			splits = Utils::splitCSV(',', line);
			classdic[splits[0].data()] = splits[1];
		}
		ifs_class.close();


		std::ifstream ifs_code;
		ifs_code.open("B:/LA_Version2/Time/Airport/AirportCode.csv");
		std::getline(ifs_code, line);
		splits = Utils::splitCSV(',', line);
		std::map<std::string, std::string> codedic;
		while (ifs_code.peek() != -1)
		{
			std::getline(ifs_code, line);
			splits = Utils::splitCSV(',', line);
			codedic[QString(splits[1].data()).toLower().toLocal8Bit().data()] = splits[2];
		}
		ifs_code.close();

		for (size_t i = 0; i < allfilesUnderDir.size(); i++)
		{
			std::string filename = allfilesUnderDir[i];
			QFileInfo fileinfo(filename.data());
			QString name = fileinfo.baseName();
			if (!name.startsWith("Airport", Qt::CaseSensitivity::CaseInsensitive))
				continue;

			ShapeFile shp(filename, 0);
			int tsfield = shp.getOrCreateField("timestruct", OGRFieldType::OFTString);
			int namefield = shp.getOrCreateField("name", OGRFieldType::OFTString);
			int sccfield = shp.getOrCreateField("scc", OGRFieldType::OFTString);

			OGRFeature *poFeature;
			double sum = 0;
			shp.poLayer->ResetReading();

			while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
			{
				std::string name = poFeature->GetFieldAsString(namefield);
				name = QString(name.data()).toLower().toLocal8Bit().data();
				std::string scc = poFeature->GetFieldAsString(sccfield);

				std::string codename = "";
				if (codedic.find(name) != codedic.end())
					codename = codedic[name];

				std::string classname = "";
				if (classdic.find(scc) != classdic.end())
					classname = classdic[scc];


				if (classname == "" || codename == "")
				{
					poFeature->SetField(tsfield, "0");
				}
				else
				{
					std::stringstream sstsfilesrc;
					sstsfilesrc << tsdir << codename << "_" << classname << ".txt";
					std::string idname = codename + "_" + classname;
					idname = QString(idname.data()).toLower().toLocal8Bit().data();
					TSFileInfo tsfi;
					if (tsdic.find(idname) == tsdic.end())
					{
						poFeature->SetField(tsfield, "0");
					}
					else
					{
						tsfi = tsdic[idname];
						std::stringstream ssts;
						ssts << "8" << tsfi.idx;
						std::stringstream sstsfiledest;
						sstsfiledest << "B:/LA_Version2/Time/2011/Airport/" << ssts.str() << ".txt";
						QFile::copy(sstsfilesrc.str().data(), sstsfiledest.str().data());
						poFeature->SetField(tsfield, ssts.str().data());
					}

				}
				shp.poLayer->SetFeature(poFeature);
				OGRFeature::DestroyFeature(poFeature);
			}
		}
	}

}
void createTimeStructForComNonPoint(std::string filename)
{
	ShapeFile ComNonPoint(filename.data(), 1);

	int idx = 0;
	OGRFeatureDefn* featureDef;
	OGRFeature* poFeature;
	idx = ComNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = ComNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = ComNonPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string bt = int2string(poFeature->GetFieldAsInteger("bt"));
		if (bt == "" || atof(bt.data()) <= 0)
		{
			bt = "91";
			poFeature->SetField("bt", 91);
		}
		poFeature->SetField(idx, ("1" + bt).data());
		ComNonPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}
void createTimeStructForResNonPoint(std::string filename)
{
	ShapeFile ResNonPoint(filename.data(), 1);

	int idx = 0;
	OGRFeatureDefn* featureDef;
	OGRFeature* poFeature;
	idx = ResNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = ResNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = ResNonPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string bt = int2string(poFeature->GetFieldAsInteger("bt"));
		if (bt == "" || atof(bt.data()) <= 0)
		{
			poFeature->SetField(idx, "0");
		}
		else
		{
			poFeature->SetField(idx, ("2" + bt).data());
		}
		ResNonPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}


}
void createTimeStructForIndNonPoint(std::string filename)
{
	ShapeFile IndNonPoint(filename.data(), 1);
	int idx;
	OGRFeatureDefn* featureDef;
	OGRFeature* poFeature;
	idx = IndNonPoint.getOrCreateField("timestruct", OFTString);
	featureDef = IndNonPoint.poLayer->GetLayerDefn();
	while ((poFeature = IndNonPoint.poLayer->GetNextFeature()) != NULL)
	{
		std::string bt = int2string(poFeature->GetFieldAsInteger("bt"));
		//if (bt == "" || atof(bt.outlineDS()) <= 0)
		//{
		//	poFeature->SetField(idx, "0");
		//}
		//else
		//{
		//	//poFeature->SetField(idx, ("3" + bt).outlineDS());
		//	poFeature->SetField(idx, "3000");
		//}
		poFeature->SetField(idx, "3000");
		IndNonPoint.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
}
void copyBuildingFiles(std::vector<std::string> dirs, std::string outrootdir)
{
	std::vector<std::string> buildingfiles;
	buildingfiles.push_back("ResNonPoint.shp");
	buildingfiles.push_back("ComNonPoint.shp");
	buildingfiles.push_back("IndNonPoint.shp");
	QDir qrootdir(outrootdir.data());
	qrootdir.mkpath(".");
	std::vector<std::string> fields;
	fields.push_back("bt"); fields.push_back("bc"); fields.push_back("yearbuilt"); fields.push_back("tfs");


	for (size_t i = 0; i < dirs.size(); i++)
	{
		std::string dir = dirs[i];
		QDir qoutdir = qrootdir.absolutePath() + "/" + QDir(dir.data()).dirName();
		qoutdir.mkpath(".");
		std::string outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();

		for (size_t j = 0; j < 3; j++)
		{
			std::string srcfile = dir + buildingfiles[j];
			std::string destcfile = outdir + buildingfiles[j];
			ShapeFile::copy(srcfile, destcfile, fields);
		}

	}
}


void reallocateOnRoad(std::string filename, std::string fieldnames[], double totals[], int numyears, std::string weightfield)
{
	ShapeFile shp(filename, 1);

	int weightIdx = shp.poLayer->GetLayerDefn()->GetFieldIndex(weightfield.data());
	std::vector<int> outfields;

	for (int i = 0; i < numyears; i++)
	{
		int outputIdx = -1;
		if (totals[i] > 0)
			outputIdx = shp.getOrCreateField(fieldnames[i].data(), OGRFieldType::OFTReal);
		outfields.push_back(outputIdx);
	}

	OGRFeature *poFeature;
	double sum = 0;
	shp.poLayer->ResetReading();
	std::vector<double> weights;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		weights.push_back(poFeature->GetFieldAsDouble(weightIdx));
		sum += weights[weights.size() - 1];
		OGRFeature::DestroyFeature(poFeature);
	}

	shp.poLayer->ResetReading();
	int idx = -1;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		idx++;
		for (int i = 0; i < numyears; i++)
		{
			if (outfields[i] > 0)
			{
				double val = weights[idx] / sum * totals[i];
				poFeature->SetField(outfields[i], val);
			}

		}
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}

double calAverageDistance(std::string shapefile)
{

	ANNSearchEngine searchEngine;

	ShapeFile shp(shapefile);
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	std::vector<OGRPoint> points;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGREnvelope env;
		fea->GetGeometryRef()->getEnvelope(&env);
		OGRPoint newp;
		newp.setX((env.MaxX + env.MinX) * 0.5);
		newp.setY((env.MaxY + env.MinY) * 0.5);
		points.push_back(newp);
		OGRFeature::DestroyFeature(fea);
	}
	searchEngine.Create(points, 2);
	double sumdist = 0;

	for (size_t i = 0; i < points.size(); i++)
	{
		int numresults = searchEngine.Select_Nearest_Points(points[i]);
		double dist = 0;
		double maxdist = 0;
		int nearestStationID = -1;
		for (size_t n = 0; n < numresults; n++)
		{
			double dist = 0;
			int selectedIdx = -1;
			searchEngine.Get_Selected_Point(n, selectedIdx, dist);
			if (maxdist < dist)
			{
				maxdist = dist;
			}

		}
		sumdist += sqrt(maxdist);
	}
	double averageDist = sumdist / points.size();
	printf("%s,%f\n", shapefile.data(), averageDist);
	return averageDist;
}

struct HourGrid
{
	float* data;
	double adfTransform[6];
	int nrows;
	int ncols;
	HourGrid(std::string gridfilename)
	{
		GDALDataset* dataset = (GDALDataset*)GDALOpen(gridfilename.data(), GA_ReadOnly);
		nrows = dataset->GetRasterYSize();
		ncols = dataset->GetRasterXSize();
		dataset->GetGeoTransform(adfTransform);
		int num = nrows*ncols;
		data = new float[num];
		GDALRasterBand *band = dataset->GetRasterBand(1);
		band->RasterIO(GF_Read, 0, 0, ncols, nrows, data, ncols, nrows, GDT_Float32, 0, 0);
		GDALClose(dataset);
	}
	~HourGrid()
	{
		delete[] data;
	}

	double getValue(double x, double y)
	{
		int yindex = (int)((y - adfTransform[3]) / adfTransform[5]);
		int xindex = (int)((x - adfTransform[0]) / adfTransform[1]);
		if (yindex < 0) {
			yindex = 0;
		}
		if (yindex > nrows - 1) {
			yindex = nrows - 1;
		}
		if (xindex < 0) {
			xindex = 0;
		}
		if (xindex > ncols - 1) {
			xindex = ncols - 1;
		}

		int idx = yindex*ncols + xindex;
		double val = data[idx];
		if (val == -3.4028234663852886e+38 || val < 0)
			val = 0;
		return val;
	}

};



void normalizeSeries(double dailyprofile[24])
{
	double sum = 0;
	for (size_t i = 0; i < 24; i++)
	{
		sum += dailyprofile[i];
	}

	for (size_t i = 0; i < 24; i++)
	{
		dailyprofile[i] = dailyprofile[i] / sum;
	}
}
void createHourlyProfiles(std::string shapefile, std::string rasterDir, std::string outdir)
{
	std::vector<HourGrid*> surfaces;
	for (int i = 0; i < 24; i++)
	{
		std::stringstream ssfilename;
		ssfilename << rasterDir << i << "/hdr.adf";
		surfaces.push_back(new HourGrid(ssfilename.str()));
	}

	ShapeFile shp(shapefile, 1);
	int timestructField = shp.getOrCreateField("timestruct", OGRFieldType::OFTString);
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	std::vector<OGRPoint> points;
	double dailyprofile[24];
	int id = 0;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGREnvelope env;
		fea->GetGeometryRef()->getEnvelope(&env);
		OGRPoint center;
		double x = (env.MaxX + env.MinX) * 0.5;
		double y = (env.MaxY + env.MinY) * 0.5;

		std::stringstream ssTS;
		ssTS << "600" << id;
		fea->SetField(timestructField, ssTS.str().data());
		std::stringstream ssoutfile;
		ssoutfile << outdir << "600" << id << ".txt";
		id++;

		std::ofstream ofs;
		ofs.open(ssoutfile.str().data());
		for (int i = 0; i < 24; i++)
		{
			double val = surfaces[i]->getValue(x, y);
			dailyprofile[i] = val;

		}
		normalizeSeries(dailyprofile);
		for (int i = 0; i < 24; i++)
		{
			ofs << dailyprofile[i] << std::endl;
		}
		ofs.close();
		shp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}


}
#include "TemporalGridding.h"
std::vector<std::string> constructShapeFilelist(std::string dir, std::string* filearr, int numoffiles)
{
	std::vector<std::string> files;
	for (size_t i = 0; i < numoffiles; i++)
	{
		files.push_back(dir + filearr[i] + ".shp");
	}
	return files;
}

double calSum(const char* filename)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);
	if (!poDataset)
		return 0;
	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	if (poDataset == NULL) {
		return NULL;
	}
	float *data;
	int   xsize = poBand->GetXSize();
	int   ysize = poBand->GetYSize();
	double adftransform[6];
	poDataset->GetGeoTransform(adftransform);
	double cellsize = adftransform[1];
	int   ncells = xsize*ysize;
	data = (float *)malloc(sizeof(float)*xsize*ysize);


	double sum = 0;
	for (size_t i = 0; i < poDataset->GetRasterCount(); i++)
	{
		poBand = poDataset->GetRasterBand(i + 1);
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
			data, xsize, ysize, GDT_Float32,
			0, 0);
		for (float* pdata = data; pdata< data + ncells; pdata++)
		{
			float val = *pdata;
			if (val <= 0 || val != val || val > 1000000000 || isinf(val))
				continue;
			sum += val;
		}


	}
	free(data);
	GDALClose((GDALDatasetH)poDataset);
	return sum;
}



//0.002
//-86.695 -85.5750000000026 
//40.3797990008037 39.3377990008061

//-86.6609878540039 -85.6339029390365
//40.1509819030762 39.3671539416537
//
double calSum(const char* filename, int start, int num)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);
	if (!poDataset)
		return 0;
	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	if (poDataset == NULL) {
		return NULL;
	}
	float *data;
	int   xsize = poBand->GetXSize();
	int   ysize = poBand->GetYSize();
	double adftransform[6];
	poDataset->GetGeoTransform(adftransform);
	double cellsize = adftransform[1];
	int   ncells = xsize*ysize;
	data = (float *)malloc(sizeof(float)*xsize*ysize);

	double ymax = 40.3797990008037;
	double xmin = -86.695;
	double resol = 0.002;
	double sum = 0;
	for (size_t i = start; i < start + num; i++)
	{
		int bandnum = i + 1;
		if (bandnum > poDataset->GetRasterCount())
			break;
		poBand = poDataset->GetRasterBand(bandnum);
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
			data, xsize, ysize, GDT_Float32,
			0, 0);
		float* pdata = data;

		for (int irow = 0; irow < ysize; irow++)
		{
			double y = ymax - irow * resol - resol * 0.5;
			for (int icol = 0; icol < xsize; icol++)
			{
				double x = xmin + icol * resol + resol * 0.5;
				//int idx = irow 
				float val = *pdata;
				if (val <= 0 || val != val || val > 10000000000 || isinf(val) || isnan(val))
				{
					val = 0;
				}
				else
				{
					//-86.6609878540039 -85.6339029390365
					//40.1509819030762 39.3671539416537
					//if (y <= 40.1509819030762 && y >= 39.3671539416537 && x <= -85.6339029390365 && x >= -86.6609878540039)
					sum += val;
				}

				*pdata++;
			}
		}

	}
	free(data);
	GDALClose((GDALDatasetH)poDataset);
	return sum;
}

std::vector<float*> loadFFDASSubset(const char* filename, double*& total)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);

	GDALRasterBand *poBand = poDataset->GetRasterBand(1);

	int   xsize = poBand->GetXSize();
	int   ysize = poBand->GetYSize();
	double adftransform[6];
	poDataset->GetGeoTransform(adftransform);
	double cellsize = adftransform[1];
	int   ncells = xsize*ysize;
	total = new double[ncells];
	std::vector<float*> datasets;
	double sum = 0;
	memset(total, 0, sizeof(double)*ncells);
	for (size_t i = 0; i < poDataset->GetRasterCount(); i++)
	{
		float* data = new float[xsize*ysize];
		poBand = poDataset->GetRasterBand(i + 1);
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
			data, xsize, ysize, GDT_Float32,
			0, 0);
		double* ptotal = total;
		int idx = 0;
		for (float* pdata = data; pdata< data + ncells; pdata++)
		{
			float val = *pdata;
			if (val <= 0 || val != val || val > 1000000000 || isinf(val))
			{
				*pdata = 0;
				val = 0;
			}
			*ptotal = *ptotal + (double)val;
			*ptotal++;
			//if (idx == 1000)
			//{
			//	printf("%f,%f\n", *ptotal, *pdata);
			//}
			idx++;
			//printf("%f,%f\n", *ptotal, *pdata);
		}
		datasets.push_back(data);
	}
	GDALClose((GDALDatasetH)poDataset);
	return datasets;
}

void removeShippingAviationFromFFDAS(const char* filename, float* shipping, float* aviation)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename, GA_Update);
	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	int   xsize = poBand->GetXSize();
	int   ysize = poBand->GetYSize();
	double adftransform[6];
	poDataset->GetGeoTransform(adftransform);
	double cellsize = adftransform[1];
	int   ncells = xsize*ysize;
	float* data = new float[ncells];
	for (size_t iband = 0; iband < poDataset->GetRasterCount(); iband++)
	{
		memset(data, 0, sizeof(float)*ncells);
		poBand = poDataset->GetRasterBand(iband + 1);
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
			data, xsize, ysize, GDT_Float32,
			0, 0);
		for (size_t icell = 0; icell < ncells; icell++)
		{
			float val = data[icell];
			if (val <= 0 || val != val || val > 1000000000 || isinf(val))
			{
				val = 0;
			}
			val = val - shipping[icell] - aviation[icell];
			if (val < 0.00001)
				val = 0;
			data[icell] = val;
		}
		poBand->RasterIO(GF_Write, 0, 0, xsize, ysize,
			data, xsize, ysize, GDT_Float32,
			0, 0);
	}
	GDALClose((GDALDatasetH)poDataset);

}
void removeAviationFromFFDAS(const char* filename, float* aviation)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename, GA_Update);
	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	int   xsize = poBand->GetXSize();
	int   ysize = poBand->GetYSize();
	double adftransform[6];
	poDataset->GetGeoTransform(adftransform);
	double cellsize = adftransform[1];
	int   ncells = xsize*ysize;
	float* data = new float[ncells];
	for (size_t iband = 0; iband < poDataset->GetRasterCount(); iband++)
	{
		memset(data, 0, sizeof(float)*ncells);
		poBand = poDataset->GetRasterBand(iband + 1);
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
			data, xsize, ysize, GDT_Float32,
			0, 0);
		for (size_t icell = 0; icell < ncells; icell++)
		{
			float val = data[icell];
			if (val <= 0 || val != val || val > 1000000000 || isinf(val))
			{
				val = 0;
			}
			val = val - aviation[icell];
			if (val < 0.00001)
				val = 0;
			data[icell] = val;
		}
		poBand->RasterIO(GF_Write, 0, 0, xsize, ysize,
			data, xsize, ysize, GDT_Float32,
			0, 0);
	}
	GDALClose((GDALDatasetH)poDataset);

}

void checkLA()
{

	SubdirManager dirmanager("B:/LA_Version2/gridPrep_SHP_master/");



	std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[1]{ "resnonpoint" }, 1));
	std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[2]{ "RailroadPoint","Railroad" }, 2));
	std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[1]{ "OnRoad" }, 1));
	std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[1]{ "NonRoad" }, 1));
	std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[1]{ "ElecProd" }, 1));
	std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(constructShapeFilelist("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));

	std::vector<std::string> testshapes = comshapefiles;
	std::string cafield = "ca14";
	for (size_t i = 0; i < testshapes.size(); i++)
	{
		double sum = ShapeFile::getTotal(testshapes[i], "ca14");
		printf("%f,%s\n", sum, testshapes[i].data());
	}
	//testshapes = indshapefiles;
	//for (size_t i = 0; i < testshapes.size(); i++)
	//{
	//	double sum = ShapeFile::getTotal(testshapes[i], cafield);
	//	printf("%f,%s\n", sum, testshapes[i].outlineDS());
	//}
	//testshapes = resshapefiles;
	//for (size_t i = 0; i < testshapes.size(); i++)
	//{
	//	double sum = ShapeFile::getTotal(testshapes[i], cafield);
	//	printf("%f,%s\n", sum, testshapes[i].outlineDS());
	//}
	//testshapes = railroadshapefiles;
	//for (size_t i = 0; i < testshapes.size(); i++)
	//{
	//	double sum = ShapeFile::getTotal(testshapes[i], cafield);
	//	printf("%f,%s\n", sum, testshapes[i].outlineDS());
	//}
	testshapes = onroadshapefiles;
	for (size_t i = 0; i < testshapes.size(); i++)
	{
		double sum = ShapeFile::getTotal(testshapes[i], cafield);
		printf("%f,%s\n", sum, testshapes[i].data());
	}

	testshapes = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/").findFilesMatch("NonRoad");
	printf("%f,%s\n", ShapeFile::getTotal(testshapes, cafield), testshapes[0].data());
	testshapes = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/Orange/").findFilesMatch("NonRoad");
	printf("%f,%s\n", ShapeFile::getTotal(testshapes, cafield), testshapes[0].data());
	testshapes = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/").findFilesMatch("NonRoad");
	printf("%f,%s\n", ShapeFile::getTotal(testshapes, cafield), testshapes[0].data());
	testshapes = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/Riverside/").findFilesMatch("NonRoad");
	printf("%f,%s\n", ShapeFile::getTotal(testshapes, cafield), testshapes[0].data());
	testshapes = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/Ventura/").findFilesMatch("NonRoad");
	printf("%f,%s\n", ShapeFile::getTotal(testshapes, cafield), testshapes[0].data());

	//testshapes = elecprodshapefiles;
	//for (size_t i = 0; i < testshapes.size(); i++)
	//{
	//	double sum = ShapeFile::getTotal(testshapes[i], cafield);
	//	printf("%f,%s\n", sum, testshapes[i].outlineDS());
	//}
	//testshapes = marineshapefiles;
	//for (size_t i = 0; i < testshapes.size(); i++)
	//{
	//	double sum = ShapeFile::getTotal(testshapes[i], cafield);
	//	printf("%f,%s\n", sum, testshapes[i].outlineDS());
	//}
	/*int resol = 1000;
	std::stringstream fishnetname;
	fishnetname << "fishnet" << resol << "m";
	std::stringstream resolname;
	resolname << resol << "m";
	std::string lidarFisheyeDir = "C:/HestiaGridding/Baltimore/2011/" + resolname.str() + "/";
	QDir(lidarFisheyeDir.outlineDS()).mkpath(".");
	TemporalGridder gridder(lidarFisheyeDir);
	gridder.fromFishnetRaster("B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");
	gridder.loadtimestruct("B:/Baltimore/Time2011/", "C:/Baltimore/2011timestructs.bin");
	std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";
	std::vector<std::string> comshapefiles = constructShapeFilelist(indir, new std::string[2]{ "comnonpoint","compoint" }, 2);
	std::vector<std::string> indshapefiles = constructShapeFilelist(indir, new std::string[2]{ "indnonpoint","indpoint" }, 2);
	std::vector<std::string> resshapefiles = constructShapeFilelist(indir, new std::string[1]{ "resnonpoint" }, 1);
	std::vector<std::string> railroadshapefiles = constructShapeFilelist(indir, new std::string[2]{ "RailroadPoint","Railroad" }, 2);
	std::vector<std::string> onroadshapefiles = constructShapeFilelist(indir, new std::string[1]{ "OnRoad" }, 1);
	std::vector<std::string> elecprodshapefiles = constructShapeFilelist(indir, new std::string[1]{ "ElecProd" }, 1);
	std::vector<std::string> marineshapefiles = constructShapeFilelist(indir, new std::string[2]{ "CMVPort","CMVUnderway" }, 2);
	gridder.addSectorGrid(comshapefiles, "com_annual.nc", "com_hourly.nc");
	gridder.addSectorGrid(indshapefiles, "ind_annual.nc", "ind_hourly.nc");
	gridder.addSectorGrid(indshapefiles, "res_annual.nc", "res_hourly.nc");
	gridder.addSectorGrid(railroadshapefiles, "railroad_annual.nc", "railroad_hourly.nc");
	gridder.addSectorGrid(onroadshapefiles, "onroad_annual.nc", "onroad_hourly.nc");
	gridder.addSectorGrid(elecprodshapefiles, "elecprod_annual.nc", "elecprod_hourly.nc");
	gridder.addSectorGrid(marineshapefiles, "marine_annual.nc", "marine_hourly.nc");
	gridder.makeHourlyTotal("total_annual.nc", "total_hourly.nc");*/

}
void export2PUBLIC()
{
	//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"timestruct" }, 6);

	std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"timestruct" }, 6);

	std::string fuels[5] = { "g","d","ng","p","" };
	std::vector<std::string> yearnames = Utils::buildVector("", new std::string[5]{ "10","11","12","13","14" }, 5);
	for (size_t iyear = 0; iyear < 5; iyear++)
	{

		for (size_t ifuel = 0; ifuel < 5; ifuel++)
		{
			std::stringstream canamess;
			canamess << "ca" << yearnames[iyear];
			if (fuels[ifuel] != "")
				canamess << "_" << fuels[ifuel];
			std::string fuelname = canamess.str().data();
			printf("%s\n", fuelname.data());
			fields2keep.push_back(fuelname.data());
		}
	}
	fields2keep.push_back("timestruct");
	//ShapeFile::copyDir("B:/Baltimore/gridPrep_SHP_master/StatePlane/", "B:/Baltimore/public/", fields2keep);
	std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	std::vector<std::string> subdirs = Utils::findSubdirectories(indir);

	for (size_t i = 0; i < subdirs.size(); i++)
	{
		std::string foldername = QDir(subdirs[i].data()).dirName().toLocal8Bit().data();
		//Utils::updateFootprintForDir(subdirs[i],true);
		//ShapeFile::copyDir(subdirs[i], "B:/LA_Version2/public/" + foldername + "/", fields2keep);



		ShapeFile::copy(subdirs[i] + "ResNonPoint.shp", "B:/LA_Version2/public/" + foldername + "/" + "ResNonPoint.shp", fields2keep);
		ShapeFile::copy(subdirs[i] + "ComNonPoint.shp", "B:/LA_Version2/public/" + foldername + "/" + "ComNonPoint.shp", fields2keep);
		//ShapeFile::copy(subdirs[i] + "OnRoad.Shp", "Z:/Hestia/LAbasin/LAbasin_v2.0/public_shapefiles/" + foldername + "/" + "OnRoad.Shp", fields2keep);
		//ShapeFile::copy(subdirs[i] + "OnRoad.Shp", "Z:/Hestia/LAbasin/LAbasin_v2.0/gridPrep_SHP_master/" + foldername + "/" + "OnRoad.Shp");


	}
}


struct GridSectorConfig
{
	std::string sectorname;
	std::vector<std::string> shapefiles;
	std::vector<std::string> timestructfiles;
	GridSectorConfig(){}
	GridSectorConfig(std::string _sectorname, std::vector<std::string> _shapefiles, std::vector<std::string> _timestructfiles) {
		sectorname = _sectorname;
		shapefiles = _shapefiles;
		timestructfiles = _timestructfiles;
	}
	GridSectorConfig(std::string _sectorname, std::vector<std::string> _shapefiles) {
		sectorname = _sectorname;
		shapefiles = _shapefiles;
	}
};
void gridLA()
{
	std::string outdir = "B:/LA_Version2/Vulcan_output/GriddedEmissions/";
	std::string shapesindir = "B:/LA_Version2/Vulcan_output/Spatial/";
	std::string intersected = "B:/LA_Version2/Vulcan_output/intersected_tmp/";
	std::string indir = "B:/LA_Version2/Vulcan_output/intersected/";
	std::string timedir = "B:/LA_Version2/Vulcan_output/time/";
	std::string fishnetraster = "B:/LA_Version2/gridPrep_SHP_master/fishnet.tif";
	QDir(outdir.data()).mkdir(".");
	QDir(shapesindir.data()).mkdir(".");
	QDir(intersected.data()).mkdir(".");
	QDir(indir.data()).mkdir(".");

	Utils::updateFootprintForDir(shapesindir, false);
	Preprocessor::gridFolderByRaster(shapesindir, intersected, fishnetraster);
	std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"ca15","timestruct" }, 7);
	ShapeFile::copyDirDropGeometry(intersected, indir, fields2keep);

	std::string years[]{ "2010","2011","2012","2013","2014","2015" };
	std::string cafields[]{ "ca10","ca11","ca12","ca13","ca14","ca15" };
	std::map<std::string, std::vector<std::string>> timestructfile_crosswalk;

	std::vector<GridSectorConfig> sectors;
	GridSectorConfig com("com", Utils::buildVector(new std::string[2]{ "comnonpoint","compoint" }, 2),
		Utils::buildVector(new std::string[2]{ "NonPoint","SMOKE" }, 2));
	GridSectorConfig res("res", Utils::buildVector(new std::string[1]{ "resnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "NonPoint" }, 1));
	GridSectorConfig ind("ind", Utils::buildVector(new std::string[2]{ "indnonpoint","indpoint" }, 2),
		Utils::buildVector(new std::string[2]{ "NonPoint","SMOKE" }, 2));
	GridSectorConfig railroad("railroad", Utils::buildVector(new std::string[1]{ "Railroad" }, 1),
		Utils::buildVector(new std::string[1]{ "SMOKE" }, 1));
	GridSectorConfig nonroad("nonroad", Utils::buildVector(new std::string[1]{ "NonRoad" }, 1),
		Utils::buildVector(new std::string[2]{ "NonRoad","SMOKE" }, 2));
	GridSectorConfig onroad("onroad", Utils::buildVector(new std::string[1]{ "OnRoad.shp" }, 1),
		Utils::buildVector(new std::string[1]{ "OnRoad" }, 1));
	GridSectorConfig cmv("cmv", Utils::buildVector(new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2),
		Utils::buildVector(new std::string[1]{ "cmv" }, 1));
	GridSectorConfig airport("airport", Utils::buildVector(new std::string[1]{ "airport" }, 1),
		Utils::buildVector(new std::string[1]{ "airport" }, 1));
	GridSectorConfig elecprod("elecprod", Utils::buildVector(new std::string[2]{ "ElecProd","ElecNEI" }, 2),
		Utils::buildVector(new std::string[2]{ "ElecProd" ,"SMOKE" }, 2));
	GridSectorConfig cement("cement", Utils::buildVector(new std::string[1]{ "cement" }, 1));
	sectors.push_back(com);
	sectors.push_back(res);
	sectors.push_back(ind);
	sectors.push_back(onroad);
	sectors.push_back(railroad);
	sectors.push_back(nonroad);
	sectors.push_back(elecprod);
	sectors.push_back(airport);
	sectors.push_back(cmv);
	sectors.push_back(cement);

	QDir(outdir.data()).mkpath(".");
	std::string cityname = "LAbasin";
	std::string version = "v2.5";
	Grid grid;
	grid.fromFishnetRaster(fishnetraster);
	for (int i = 0; i < 6; i++)
	{

		std::string year = years[i];
		std::string attName = "ca" + year.substr(2, 2);
		grid.reset(attName);
		int numhous = 8760;
		if (year == "2012")
		{
			numhous = 8784;
		}
		TemporalGridder gridder(outdir + year + "/", numhous);


		gridder.fromFishnetRaster(fishnetraster);
		std::string timestructdir= timedir + year + "/";
		//TimestructTool::normalizeBinary(timestructfile);

		SubdirManager dirmanager(indir);
		for (int nsector = 0; nsector < sectors.size(); nsector++)
		{
			std::vector<std::string> shapefiles = dirmanager.findFilesMatch(sectors[nsector].shapefiles);
			HestiaGrid* sectorGrid = gridder.addSectorGrid(shapefiles, sectors[nsector].sectorname);
			if (QFileInfo(std::string(outdir + cityname + ".total.hourly." + year + "." + version + ".nc").data()).exists())
				continue;
			for (int ntimestruct = 0; ntimestruct < sectors[nsector].timestructfiles.size(); ntimestruct++)
			{
				sectorGrid->loadtimestruct(timestructdir + sectors[nsector].timestructfiles[ntimestruct] + ".bin");
			}
			
		}

		if (!QFileInfo(std::string(outdir + cityname + ".total.hourly." + year + "." + version + ".nc").data()).exists()){
			gridder.loadAttribute(cafields[i]);
		}
		if (!QFileInfo(std::string(outdir + cityname + ".total.annual." + year + "." + version + ".nc").data()).exists()) {
			gridder.makeAnnualTotal(outdir + cityname + ".total.annual." + year + "." + version + ".nc");
		}

		for (size_t isector = 0; isector < gridder.sectorGrids.size(); isector++)
		{
			HestiaGrid* sectorGrid = gridder.sectorGrids[isector];


			printf("%s\n", std::string(outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".shp").data());
			if (!QFileInfo(std::string(outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".shp").data()).exists())
			{	
				sectorGrid->getTotal();
				sectorGrid->toShapefile(shapesindir + "onroad.shp", outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".shp");
			}

			if (QFileInfo(std::string(outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".shp").data()).exists()
				&& !QFileInfo(std::string(outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".tif").data()).exists())
			{
				sectorGrid->getTotal();
				sectorGrid->toShapefile(shapesindir + "onroad.shp", outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".shp");
				grid.reset(attName);
				grid.gatherCells(outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".shp", attName.data());
				grid.toRaster(outdir + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".tif");
			}


		}

	
		if (!QFileInfo(std::string(outdir + cityname + ".total.annual." + year + "." + version + ".shp").data()).exists()){	
			gridder.toShapefile(shapesindir + "onroad.shp", outdir + cityname + ".total.annual." + year + "." + version + ".shp");
		}

		if (QFileInfo(std::string(outdir + cityname + ".total.annual." + year + "." + version + ".shp").data()).exists()
			&& !QFileInfo(std::string(outdir + cityname + ".total.annual." + year + "." + version + ".tif").data()).exists())
		{
			grid.reset(attName);
			grid.gatherCells(outdir + cityname + ".total.annual." + year + "." + version + ".shp", attName.data());
			grid.toRaster(outdir + cityname + ".total.annual." + year + "." + version + ".tif");
		}
		if (!QFileInfo(std::string(outdir + cityname + ".total.hourly." + year + "." + version + ".nc").data()).exists())
		{
			gridder.makeHourlyTotal(outdir + cityname + ".total.hourly." + year + "." + version + ".nc");
			for (size_t isector = 0; isector < gridder.sectorGrids.size(); isector++)
			{
				HestiaGrid* sectorGrid = gridder.sectorGrids[isector];
				sectorGrid->clearTimeStruct();
			}
		}
	}



}
void mergeLA()
{

	//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[3]{ "ca11","length","area"}, 3);

	std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	std::vector<std::string> subdirs;// = Utils::findSubdirectories(indir);
	subdirs.push_back(indir + "Los_Angeles/");
	subdirs.push_back(indir + "Orange/");
	subdirs.push_back(indir + "Riverside/");
	subdirs.push_back(indir + "San_Bernardino/");
	subdirs.push_back(indir + "Ventura/");
	std::string outdir = indir + "Vulcan/";
	//std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"timestruct" }, 6);
	//std::string outdir = "B:/LA_Version2/Merged/";
	//std::vector<std::string> subdirs = Utils::findSubdirectories(indir);
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	std::string foldername = QDir(subdirs[i].data()).dirName().toLocal8Bit().data();
	//	ShapeFile::copyDir(indir + foldername + "/", outdir + foldername + "/", fields2keep);
	//}

	std::map<std::string, std::vector<std::string>> sectors;
	for (size_t i = 0; i < subdirs.size(); i++)
	{
		std::string foldername = QDir(subdirs[i].data()).dirName().toLocal8Bit().data();
		QDir input_dir(subdirs[i].data());
		input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
		input_dir.setSorting(QDir::Name);
		indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
		QFileInfoList list = input_dir.entryInfoList();
		for (int i = 0; i < list.size(); ++i) {
			QFileInfo fileInfo = list.at(i);
			if (!fileInfo.fileName().endsWith(".shp", Qt::CaseInsensitive))
				continue;
			//if (!fileInfo.fileName().endsWith("Railroad.shp", Qt::CaseInsensitive) && !fileInfo.fileName().contains("Nonroad_", Qt::CaseInsensitive))
			//	continue;
			//if (!fileInfo.fileName().endsWith("Railroad.shp", Qt::CaseInsensitive) )
			//	continue;
			if (!fileInfo.fileName().contains("Nonroad_", Qt::CaseInsensitive))
				continue; 
			/*if (!fileInfo.fileName().endsWith("NonPoint.shp", Qt::CaseInsensitive))
				continue;*/
			std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();// outdir + foldername + "/" + fileInfo.fileName().toLocal8Bit().data();
			std::string name = fileInfo.completeBaseName().toLocal8Bit().data();
			std::map<std::string, std::vector<std::string>>::iterator iter = sectors.find(name);
			if (iter == sectors.end())
			{
				std::vector<std::string> newlist;
				newlist.push_back(input_file);
				sectors[name] = newlist;
			}
			else
			{
				iter->second.push_back(input_file);
			}
		}
	}
	std::map<std::string, std::vector<std::string>>::iterator iter = sectors.begin();
	while (iter != sectors.end())
	{

		std::vector<std::string> filelist = iter->second;
		std::string outputfile = outdir + iter->first + ".shp";
		ShapeFile::copy(filelist[0], outputfile);
		ShapeFile master(outputfile, 1);
		for (int n = 0; n < master.poLayer->GetLayerDefn()->GetFieldCount(); n++)
		{
			printf("%s,", master.poLayer->GetLayerDefn()->GetFieldDefn(n)->GetNameRef());
		}
		printf("\n");
		for (size_t i = 1; i < filelist.size(); i++)
		{
			ShapeFile child(filelist[i]);
			if (child.poLayer->GetLayerDefn()->GetFieldCount() != master.poLayer->GetLayerDefn()->GetFieldCount())
			{
				printf("%s,%s\n", outputfile.data(), filelist[i].data());
			}
			for (int n = 0; n < master.poLayer->GetLayerDefn()->GetFieldCount(); n++)
			{
				printf("%s,", child.poLayer->GetLayerDefn()->GetFieldDefn(n)->GetNameRef());
			}
			printf("\n");
			OGRFeature* fea;
			while ((fea = child.poLayer->GetNextFeature()) != NULL)
			{
				OGRFeature* newfea = OGRFeature::CreateFeature(master.poLayer->GetLayerDefn());
				newfea->SetGeometry(fea->GetGeometryRef());
				for (int n = 0; n < master.poLayer->GetLayerDefn()->GetFieldCount(); n++)
				{
					if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTReal)
					{
						newfea->SetField(n, fea->GetFieldAsDouble(n));
					}
					else if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTInteger)
					{
						newfea->SetField(n, fea->GetFieldAsInteger(n));
					}
					else if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTString)
					{
						newfea->SetField(n, fea->GetFieldAsString(n));
					}
					else if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTInteger64)
					{
						newfea->SetField(n, fea->GetFieldAsInteger64(n));
					}
					else
					{
						newfea->SetField(n, fea->GetRawFieldRef(n));
					}
				}
				master.poLayer->CreateFeature(newfea);
				OGRFeature::DestroyFeature(newfea);
				OGRFeature::DestroyFeature(fea);
			}

		}

		iter++;
	}
	
}

void mergeLANonroad()
{

	//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[3]{ "ca11","length","area"}, 3);

	std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	std::vector<std::string> subdirs;// = Utils::findSubdirectories(indir);
	subdirs.push_back(indir + "Los_Angeles/");
	subdirs.push_back(indir + "Orange/");
	subdirs.push_back(indir + "Riverside/");
	subdirs.push_back(indir + "San_Bernardino/");
	subdirs.push_back(indir + "Ventura/");
	std::string outdir = indir + "Vulcan/";
	//std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"timestruct" }, 6);
	//std::string outdir = "B:/LA_Version2/Merged/";
	//std::vector<std::string> subdirs = Utils::findSubdirectories(indir);
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	std::string foldername = QDir(subdirs[i].data()).dirName().toLocal8Bit().data();
	//	ShapeFile::copyDir(indir + foldername + "/", outdir + foldername + "/", fields2keep);
	//}
	std::string SurrogateVulcan[] = { "100","140","260", "300","310","311","350","400","505","510", "520", "525", "525", "850","860","890" };
	std::string SubCountyCode[] = { "STFID","STFID","US_RAIL_", "US_LRES_ID","US_AG2K_ID","US_AG2K_ID","US_LW2K_ID","STFID","BlkGrp_ID","BlkGrp_ID", "BlkGrp_ID", "BlkGrp_ID", "US_GOLF_ID", "US_GOLF_ID","MINE_PT_ID","US_TIMB_ID" };
		



	/*std::ifstream ifs;
	ifs.open(rootdir + "GRID_SCB_0_01.csv");
	std::string headerline;
	std::getline(ifs, headerline);
	int latfield = shp.getOrCreateField("lat", OGRFieldType::OFTReal);
	int lonfield = shp.getOrCreateField("lon", OGRFieldType::OFTReal);
	int idfield = shp.getOrCreateField("Id", OGRFieldType::OFTInteger);
	double halfcellsize = 0.01 * 0.5;
	int numcells = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		double lat = atof(splits[1].data());
		double lon = atof(splits[2].data());*/

	std::map<std::string, std::vector<std::string>> sectors;


	for (size_t sec = 0; sec < 16; sec++)
	{
		std::vector<std::string> filelist;
		std::stringstream ss;
		ss << "NonRoad_" << SurrogateVulcan[sec] << ".shp";
		for (size_t i = 0; i < subdirs.size(); i++)
		{
			std::string foldername = QDir(subdirs[i].data()).dirName().toLocal8Bit().data();
			QDir input_dir(subdirs[i].data());
			input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
			input_dir.setSorting(QDir::Name);
			indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
			filelist.push_back(indir + ss.str());
		}

		std::string outputfile = outdir + ss.str();
		ShapeFile::copy(filelist[0], outputfile);
		ShapeFile master(outputfile, 1);
		for (int n = 0; n < master.poLayer->GetLayerDefn()->GetFieldCount(); n++)
		{
			printf("%s,", master.poLayer->GetLayerDefn()->GetFieldDefn(n)->GetNameRef());
		}
		printf("\n");
		for (size_t i = 1; i < filelist.size(); i++)
		{
			ShapeFile child(filelist[i]);
			if (child.poLayer->GetLayerDefn()->GetFieldCount() != master.poLayer->GetLayerDefn()->GetFieldCount())
			{
				printf("%s,%s\n", outputfile.data(), filelist[i].data());
			}
			for (int n = 0; n < master.poLayer->GetLayerDefn()->GetFieldCount(); n++)
			{
				printf("%s,", child.poLayer->GetLayerDefn()->GetFieldDefn(n)->GetNameRef());
			}
			printf("\n");
			OGRFeature* fea;
			while ((fea = child.poLayer->GetNextFeature()) != NULL)
			{
				OGRFeature* newfea = OGRFeature::CreateFeature(master.poLayer->GetLayerDefn());
				newfea->SetGeometry(fea->GetGeometryRef());
				for (int n = 0; n < master.poLayer->GetLayerDefn()->GetFieldCount(); n++)
				{
					if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTReal)
					{
						newfea->SetField(n, fea->GetFieldAsDouble(n));
					}
					else if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTInteger)
					{
						newfea->SetField(n, fea->GetFieldAsInteger(n));
					}
					else if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTString)
					{
						newfea->SetField(n, fea->GetFieldAsString(n));
					}
					else if (newfea->GetFieldDefnRef(n)->GetType() == OGRFieldType::OFTInteger64)
					{
						newfea->SetField(n, fea->GetFieldAsInteger64(n));
					}
					else
					{
						newfea->SetField(n, fea->GetRawFieldRef(n));
					}
				}
				master.poLayer->CreateFeature(newfea);
				OGRFeature::DestroyFeature(newfea);
				OGRFeature::DestroyFeature(fea);
			}

		}
	}

	


}

#include "GDAL_DS.h"
double checkTotal(std::string filename,double scale,std::string units)
{
	//std::string infilename = "C:/HestiaGridding/Los_Angeles/2010/total_hourly_2010.nc";
	//std::string infilename = "C:/HestiaGridding/Baltimore/2011/total_hourly_2011.nc";

	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename.data(), GA_ReadOnly);
	if (!poDataset)
		return 0;
	int numbands = poDataset->GetRasterCount();
	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	int   xsize = poBand->GetXSize();
	int   ysize = poBand->GetYSize();
	int   ncells = xsize*ysize;
	double *srcData = new double[xsize*ysize];
	double sum = 0;

	//double *destData = new double[xsize*ysize];
	//memset(destData, 0, sizeof(double)*ncells);

	std::ofstream ofs;
	ofs.open(filename + ".totals.csv");
	ofs << "Layer," + units << std::endl;

	memset(srcData, 0, sizeof(double)*ncells);
	for (int i = 1; i < numbands + 1; i++)
	{

		poBand = poDataset->GetRasterBand(i);
		double nodata = poBand->GetNoDataValue();
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
			srcData, xsize, ysize, GDT_Float64,
			0, 0);
		double subtotal = 0;
		for (size_t icell = 0; icell < ncells; icell++)
		{
			if (srcData[icell] == nodata || isnan(srcData[icell]) || srcData[icell] < 0 || srcData[icell] > 10000000000  )
				continue;
			subtotal += srcData[icell];
		}
		subtotal = subtotal *scale;
		ofs << std::fixed << "layer " << i << "," << subtotal << std::endl;
		//printf("%s,%f\n", filename.data(), subtotal);
		sum += subtotal;
		//if(i % 100 == 0)
		//   printf("%d/%d\n", i,numbands);
	}
	//printf("%s,%f,%f\n", filename.data(), sum, sum / 1000000.0);
	ofs.close();
	printf("[%s]\n", filename.data());
	printf("%f\n", sum);
	//delete[] destData;
	delete[] srcData;
	GDALClose((GDALDatasetH)poDataset);
	return sum;

}
void checkTotalsInDir(std::string indir, QString ext, double scale, std::string units,std::string outname)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	//if (QFileInfo(input_dir.absolutePath() + "/totals.csv").exists())
	//	return;
	QFileInfoList list = input_dir.entryInfoList();
	double sum = 0;
	std::ofstream ofs;
	ofs.open(indir + outname);
	ofs << "File," + units << std::endl;

	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(ext, Qt::CaseInsensitive))
			continue;
		//if (fileInfo.fileName().contains("total", Qt::CaseInsensitive))
		//	continue;
		double subtotal = checkTotal(fileInfo.absoluteFilePath().toLocal8Bit().data(), scale, units);
		sum += subtotal;
		//ofs << std::fixed << std::setprecision(15) << fileInfo.completeBaseName().toLocal8Bit().data() << "," << subtotal / 1000000.0 << std::endl;
		ofs << std::fixed << fileInfo.completeBaseName().toLocal8Bit().data() << "," << subtotal << std::endl;
	}
	//ofs << "Total" << "," << sum / 1000000.0 << std::endl;
	ofs.close();
	//printf("sum=%f\n", sum /*/ 1000000.0*/);
}
void checkTotalsInDirs(std::string indir, QString ext, double scale, std::string units, std::string outname)
{
	std::vector<std::string> dirs = Utils::findSubdirectories(indir);
	for (size_t i = 0; i < dirs.size(); i++)
	{
		checkTotalsInDir(dirs[i], ext, scale, units, outname);
	}
}
void scaleRaster(std::string indir, std::string outdir, QString ext,double scale)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(ext, Qt::CaseInsensitive))
			continue;
		QString outfile = QString(outdir.data()) + "/" + fileInfo.fileName();
		QFile::copy(fileInfo.absoluteFilePath(), outfile);
		GDAL_DS<double>* ds = new GDAL_DS<double>();
		ds->open(outfile.toLocal8Bit().data(), GDALAccess::GA_Update);
		ds->scale(1, scale);
		delete ds;
	}
}



void tif2csv(std::string indir)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".tif", Qt::CaseInsensitive))
			continue;
		files.push_back(fileInfo.absoluteFilePath().toLocal8Bit().data());
		std::string infile = fileInfo.absoluteFilePath().toLocal8Bit().data();
		std::string outfile = indir + std::string(fileInfo.completeBaseName().toLocal8Bit().data()) + ".csv";
		GDAL_DS<double>* ds = new GDAL_DS<double>();
		ds->open(infile);
		double* data = ds->readData(1);
		std::ofstream ofs;
		ofs.open(outfile.data());
		for (size_t ncell = 0; ncell < ds->nrows * ds->ncols; ncell++)
		{
			ofs << (int)data[ncell] << std::endl;
		}
		ofs.close();
		delete[] data;
		delete ds;
	}

}
void aggregrateTIFF(int num, std::string indir, std::string outdir)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".tif", Qt::CaseInsensitive))
			continue;
		files.push_back(fileInfo.absoluteFilePath().toLocal8Bit().data());
		std::string infile = fileInfo.absoluteFilePath().toLocal8Bit().data();
		std::string outfile = outdir + std::string(fileInfo.fileName().toLocal8Bit().data());
		GDAL_DS<double>* ds = new GDAL_DS<double>();
		ds->open(infile);
		double* olddata = ds->readData(1);
		int newrows = ds->nrows / num;
		int newcols = ds->ncols / num;
		double* newdata = new double[newrows*newcols];
		memset(newdata, 0, newrows*newcols*sizeof(double));
		int oldidx = 0;
		for (int ioldrow = 0; ioldrow < ds->nrows; ioldrow++)
		{
			int inewrow = ioldrow / num;
			for (int ioldcol = 0; ioldcol < ds->ncols; ioldcol++)
			{
				int inewcol = ioldcol / num;
				int newidx = inewrow*newcols + inewcol;

				if (inewrow < newrows && inewcol < newcols)
					newdata[newidx] = newdata[newidx] + olddata[oldidx];
				oldidx++;
			}
		}
		delete[] olddata;
		ds->ncols = newcols;
		ds->nrows = newrows;

		ds->adfGeoTransform[1] = ds->adfGeoTransform[1] * num;
		ds->adfGeoTransform[5] = ds->adfGeoTransform[5] * num;
		//  adfGeoTransform[0] /* top left x */
		//	adfGeoTransform[1] /* w-e pixel resolution */
		//	adfGeoTransform[2] /* 0 */
		//	adfGeoTransform[3] /* top left y */
		//	adfGeoTransform[4] /* 0 */
		//	adfGeoTransform[5] /* n-s pixel resolution (negative value) */

		ds->bound.MinX = ds->adfGeoTransform[0];
		ds->bound.MaxY = ds->adfGeoTransform[3];
		ds->bound.MaxX = ds->bound.MinX + ds->adfGeoTransform[1] * newcols;
		ds->bound.MinY = ds->bound.MaxY + ds->adfGeoTransform[5] * newrows;
		ds->create(outfile);
		ds->writeData(1, newdata, 0);

		delete[] newdata;
		delete ds;
		//printf("%f,%f\n", checkTotal(infile), checkTotal(outfile));
	}
	tif2csv(indir);
	tif2csv(outdir);
}
void aggregrateTIFF(std::string indir)
{
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".csv", Qt::CaseInsensitive))
			continue;
		std::string infile = fileInfo.absoluteFilePath().toLocal8Bit().data();
		std::string outfile = indir + std::string(fileInfo.completeBaseName().toLocal8Bit().data()) + ".png";
		std::stringstream ss;
		ss << "RScript E:/LASpatialAnalysis/histogram.R" << " " << infile << " " << outfile << " " << "\"" << "Natural Log(X)" << "\"" << " " << fileInfo.completeBaseName().toLocal8Bit().data();
		system(ss.str().data());
	}

}

void mergeFields(std::string indir, std::string mergeShp)
{
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	std::vector<QFileInfo> files;
	std::vector<std::string> names;

	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".tif", Qt::CaseInsensitive) || fileInfo.fileName().endsWith("lspop.tif", Qt::CaseInsensitive))
			continue;
		names.push_back(fileInfo.completeBaseName().toLocal8Bit().data());
		files.push_back(fileInfo);
	}

	names.push_back("lspop");
	files.push_back(QFileInfo(std::string(indir + "lspop.tif").data()));

	Grid grid;
	grid.fromFishnetRaster(files[0].absoluteFilePath().toLocal8Bit().data());
	grid.reset(names);

	for (size_t i = 0; i < names.size(); i++)
	{
		std::string infile = files[i].absoluteFilePath().toLocal8Bit().data();
		GDAL_DS<double>* ds = new GDAL_DS<double>();
		ds->open(infile);
		double* data = ds->readData(1);
		memcpy(grid.cells + i*grid.slice, data, grid.slice*sizeof(double));
		printf("%s\n", infile.data());
		delete ds;
		delete[] data;
	}

	grid.toShape(grid.proj, mergeShp, true);
}

//for (size_t i = 0; i < shapefiles.size(); i++)
//{
//	std::string name = QFileInfo(shapefiles[i].outlineDS()).baseName().toLocal8Bit().outlineDS();
//	std::stringstream ss;
//	std::string infilename = indir + name + ".dbf";
//	ss << "RScript B:/LA_Version2/histogramDBF.R" << " " << infilename << " " << "ca11" << " " << indir + name + ".png" << " " << "\"" << "Natural Log(X)[Carbon = e^X(KgC)]" << "\"" << " " << name;
//	system(ss.str().outlineDS());
//	//"B:/LA_Version2/histogramDBF.R"
//}
void gridLA_TIFF()
{

	std::string fishnetfile = "B:/LA_Version2/gridoutput/fishnet1000m.shp";
	std::string reffileshp = "B:/LA_Version2/gridoutput/fishnet1000m.shp";
	std::string reffiletif = "B:/LA_Version2/gridoutput/fishnet1000m.tif";
	//Grid grid;
	//grid.fromFishnetRaster(reffiletif);
	//grid.toShape(proj, fishnetfile);
	//grid.reset();
	//grid.toRaster("B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");
	//Utils::updateFootprintForDir("B:/Baltimore/gridPrep_SHP_master/ca/");
	std::string indir = "B:/LA_Version2/gridoutput_nogeometry/";
	std::string outroot = "E:/LASpatialAnalysis/1km/";
	//Preprocessor::gridFolderByRaster("B:/Baltimore/gridPrep_SHP_master/ca_1.1/", indir, reffiletif);
	//E:/LanduseRegressionBaltimore/Landsat8/Baltimore_0.tif
	std::string years[]{ "2011" };
	std::string cafields[]{ "ca11" };
	std::string cityname = "LA";
	std::string version = "";
	//Baltimore.total.hourly.2010.v1.1.nc


	std::string year = years[0];
	TemporalGridder gridder(outroot,8760);
	gridder.fromFishnetRaster(reffiletif);



	SubdirManager dirmanager = SubdirManager(indir);
	//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";
	std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "RailroadPoint","Railroad" }, 2));
	std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.shp" }, 1));
	std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "NonRoad" }, 1));
	std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));

	gridder.addSectorGrid(comshapefiles, "com");
	gridder.addSectorGrid(indshapefiles, "ind");
	gridder.addSectorGrid(resshapefiles, "res");
	gridder.addSectorGrid(railroadshapefiles, "railroad");
	gridder.addSectorGrid(onroadshapefiles, "onroad");
	gridder.addSectorGrid(nonroadshapefiles, "nonroad");
	gridder.addSectorGrid(elecprodshapefiles, "elecprod");
	gridder.addSectorGrid(airportshapefiles, "airport");
	gridder.addSectorGrid(marineshapefiles, "marine");
	gridder.loadAttribute(cafields[0]);

	gridder.getTotal();

	GDAL_DS<double>* ds = new GDAL_DS<double>();
	gridder.toShapefile(reffileshp, outroot + "total.shp");
	ds->open(reffiletif);
	ds->create(outroot + "total.tif");
	ds->writeData(1, gridder.cells, 0);
	delete ds;

	for (size_t isector = 0; isector < gridder.sectorGrids.size(); isector++)
	{
		HestiaGrid* sectorGrid = gridder.sectorGrids[isector];
		sectorGrid->getTotal();
		std::string sectorgridfiletif = outroot + sectorGrid->sectorname + ".tif";
		ds = new GDAL_DS<double>();
		ds->open(reffiletif);
		ds->create(sectorgridfiletif);
		ds->writeData(1, sectorGrid->cells, 0);
		delete ds;
		sectorGrid->toShapefile(reffileshp, outroot + sectorGrid->sectorname + ".shp");
	}

}
void gridLA_NoShipping()
{

	std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	std::vector<std::string> subdirs = Utils::findSubdirectories(indir);
	indir = "B:/LA_Version2/gridoutput/";



	std::string fishnetfile = "B:/LA_Version2/gridoutput/fishnet1000m.shp";

	std::string year = "2011";
	TemporalGridder gridder("C:/HestiaGridding/Los_Angeles/",8760);
	gridder.fromFishnetShapeFile(fishnetfile);

	SubdirManager dirmanager = SubdirManager(indir);
	//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";
	std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "RailroadPoint","Railroad" }, 2));
	std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.shp" }, 1));
	std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "NonRoad" }, 1));
	std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));


	gridder.addSectorGrid(comshapefiles, "com");
	gridder.addSectorGrid(indshapefiles, "ind");
	gridder.addSectorGrid(resshapefiles, "res");
	gridder.addSectorGrid(railroadshapefiles, "railroad");
	gridder.addSectorGrid(onroadshapefiles, "onroad");
	gridder.addSectorGrid(nonroadshapefiles, "nonroad");
	gridder.addSectorGrid(elecprodshapefiles, "elecprod");
	gridder.addSectorGrid(marineshapefiles, "marine");
	gridder.addSectorGrid(airportshapefiles, "airport");
	gridder.loadAttribute("ca11");
	gridder.makeAnnualTotal("C:/HestiaGridding/Los_Angeles/Hestia2011.nc");
	gridder.toShapefile(subdirs[0] + "onroad.shp", "C:/HestiaGridding/Los_Angeles/Hestia2011.shp");


}
void regrid_FFDASODIAC_LA()
{
	std::string fishnetfile = "C:/HestiaGridding/Los_Angeles/fishnet2011.shp";
	std::string infileodiac = "C:/HestiaGridding/Los_Angeles/ODIAC_WGS84.shp";
	std::string infileffdas = "C:/HestiaGridding/Los_Angeles/FFDAS_WGS84.shp";

	std::string outfileodiac_tmp = "C:/HestiaGridding/ODIAC_WGS84_intersected.shp";
	std::string outfileffdas_tmp = "C:/HestiaGridding/FFDAS_WGS84_intersected.shp";
	//Utils::updateFootprint(infileodiac, 110000);
	//Utils::updateFootprint(infileffdas, 110000);
	Preprocessor::intersectWithArcGIS(infileodiac, fishnetfile, outfileodiac_tmp);
	Preprocessor::intersectWithArcGIS(infileffdas, fishnetfile, outfileffdas_tmp);
	Preprocessor::updateFieldAfterIntersection(outfileodiac_tmp, 110000);
	Preprocessor::updateFieldAfterIntersection(outfileffdas_tmp, 110000);

	std::string outfileodiac = "C:/HestiaGridding/Los_Angeles/ODIAC2011";
	std::string outfileffdas = "C:/HestiaGridding/Los_Angeles/FFDAS2011";
	Grid grid;
	grid.fromFishnetShape(fishnetfile);
	ShapeFile refshape(fishnetfile);
	ShapeFile outfileodiac_tmp_shp(outfileodiac_tmp);
	grid.reset();
	grid.gatherCells(&outfileodiac_tmp_shp, "ca11", 110000);


	char wkt[512];
	char* pwkt = wkt;
	if (refshape.poLayer->GetSpatialRef())
		refshape.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	grid.toRaster(outfileodiac + ".tif", wkt);
	grid.toShape(wkt, outfileodiac + ".shp", true);
	grid.reset();

	ShapeFile outfileffdas_tmp_shp(outfileffdas_tmp);
	grid.reset();
	grid.gatherCells(&outfileffdas_tmp_shp, "ca11", 110000);
	grid.toRaster(outfileffdas + ".tif", wkt);
	grid.toShape(wkt, outfileffdas + ".shp", true);

	//OGREnvelope regionbound;
	//regionbound.MaxX = -113.7916669; regionbound.MinX = -120.0750002; regionbound.MaxY = 36.12500022; regionbound.MinY = 32.5166669;
	//
	//int startrow = 6465;
	//int nrows = 434;
	//int startcol = 7191;
	//int ncols = 755;
	//double resol = 0.0083333333;
	//double adfTransform[6];// = grid._adfGeoTransform;
	//adfTransform[0] = -180 + startcol * resol;
	//adfTransform[3] = 90 - startrow * resol;
	//adfTransform[1] = resol;
	//adfTransform[5] = -resol;
	//adfTransform[2] = 0;
	//adfTransform[4] = 0;
	////Preprocessor::extractGlobalRaster("B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/ODIAC_2011.tif", "B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/ODIAC_2011_LA.tif", adfTransform, startrow, nrows, startcol, ncols);
	////Preprocessor::extractGlobalRaster("B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/FFDAS_2011.tif", "B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/FFDAS_2011_LA.tif", adfTransform, startrow, nrows, startcol, ncols);


	//Grid grid(adfTransform, ncols,nrows);
	//ShapeFile refshape("B:/FFDAS/urban_boundary/modis.Shp");
	//grid.reset();
	//float* outlineDS = Preprocessor::readData("B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/ODIAC_2011_LA.tif");
	//for (size_t i = 0; i < grid.slice; i++)
	//{
	//	grid.cells[i] = outlineDS[i];
	//}
	//grid.toShape(&refshape, "C:/HestiaGridding/Los_Angeles/ODIAC_WGS84.Shp",true);


	//grid.reset();
	//delete[] outlineDS;
	//outlineDS = Preprocessor::readData("B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/FFDAS_2011_LA.tif");
	//for (size_t i = 0; i < grid.slice; i++)
	//{
	//	grid.cells[i] = outlineDS[i];
	//}
	//grid.toShape(&refshape, "C:/HestiaGridding/Los_Angeles/FFDAS_WGS84.Shp", true);
	//start col	7191		start rows	6465
	//	end col	7945		end rows	6898
	//	ncols	755		nrows	434



	//std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	//std::vector<std::string> subdirs = Utils::findSubdirectories(indir);

	///*subdirs = Utils::findSubdirectories("B:/LA_Version2/gridPrep_SHP_master/");
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//updateNonpointTime(subdirs[i]);
	//}
	//subdirs = Utils::findSubdirectories("B:/LA_Version2/ca/");
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//updateNonpointTime(subdirs[i]);
	//}

	//subdirs = Utils::findSubdirectories("B:/LA_Version2/gridoutput/");
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//updateNonpointTime(subdirs[i]);
	//}*/

	//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[1]{ "ca11"}, 1);
	////for (size_t i = 0; i < subdirs.size(); i++)
	////{
	////	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
	////	Utils::updateFootprintForDir(subdirs[i],true);
	////	ShapeFile::copyDir(subdirs[i], "B:/LA_Version2/ca/" + foldername + "/", fields2keep);
	////}

	////Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Los_Angeles/", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	////Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Orange/", "B:/LA_Version2/gridPrep_SHP_master/Orange/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	////Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Riverside/", "B:/LA_Version2/gridPrep_SHP_master/Riverside/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	////Preprocessor::reprojectDir("B:/LA_Version2/nonroad/San_Bernardino/", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	////Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Ventura/", "B:/LA_Version2/gridPrep_SHP_master/Ventura/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	//std::string years[]{ "2010","2011","2012","2013","2014" };
	//std::string cafields[]{ "ca10","ca11","ca12","ca13","ca14" };
	//indir = "B:/LA_Version2/gridoutput/";

	//SubdirManager dirmanager(indir);
	//subdirs = Utils::findSubdirectories(indir);
	//std::string boundfile = indir + "bound.txt";
	//if (!QFileInfo(boundfile.outlineDS()).exists())
	//{
	//	std::vector<std::string> allfiles = dirmanager.findFilesMatch(".Shp");
	//	OGREnvelope bound = BoundManager::readBoundFromShapes(allfiles);
	//	BoundManager::writeBound(bound, boundfile);
	//}

	//int resol = 1000;
	//std::stringstream fishnetname;
	//fishnetname << "fishnet" << resol << "m";
	//std::stringstream resolname;
	//resolname << resol << "m";

	//OGREnvelope bound = BoundManager::readBound(boundfile);


	//Grid fishnetgrid(bound, resol * 3.28084, 1);
	//ShapeFile refshape(subdirs[0] + "onroad.Shp");
	//if (!QFileInfo((indir + fishnetname.str() + ".tif").outlineDS()).exists())
	//{
	//	fishnetgrid.toShape(&refshape, indir + fishnetname.str() + ".Shp");
	//	fishnetgrid.reset();
	//	fishnetgrid.toRaster(indir + fishnetname.str() + ".tif");
	//}
	////Preprocessor::gridFolderByRaster(lidarFisheyeDir + resolname.str() + "/", lidarFisheyeDir + resolname.str() + "/", "B:/Baltimore/gridPrep_SHP_master/"+ fishnetname.str() + ".tif");
	//subdirs = Utils::findSubdirectories(indir);
	////for (size_t i = 0; i < subdirs.size(); i++)
	////{
	////	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
	////	Preprocessor::gridFolderByRaster(subdirs[i], lidarFisheyeDir + foldername + "/", lidarFisheyeDir + fishnetname.str() + ".tif");

	////}
	////std::vector<std::string> testshapes = comshapefiles;
	//////Preprocessor::gridFolderByRaster("B:/Baltimore/gridPrep_SHP_master/ca11/", "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/", "B:/Baltimore/gridPrep_SHP_master/"+ fishnetname.str() + ".tif");

	//std::string lidarFisheyeDir = "C:/HestiaGridding/Los_Angeles/";
	//QDir(lidarFisheyeDir.outlineDS()).mkpath(".");


	//for (size_t i = 0; i < 5; i++)
	//{
	//	std::string year = years[i];
	//	TemporalGridder gridder(lidarFisheyeDir + year + "/");
	//	gridder.fromFishnetRaster(indir + fishnetname.str() + ".tif");
	//	std::string timestructfile = lidarFisheyeDir + +"timestructs_" + year + ".bin";
	//	//TimestructTool::normalizeBinary(timestructfile);
	//	gridder.loadtimestruct("B:/LA_Version2/Time/" + year + "/", timestructfile);
	//	dirmanager = SubdirManager(indir);
	//	//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";
	//	std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	//	std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	//	std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	//	std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "createRailroadPoint","Railroad" }, 2));
	//	std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.Shp" }, 1));
	//	std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "NonRoad" }, 1));
	//	std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	//	std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	//	std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));


	//	gridder.addSectorGrid(comshapefiles, "com");
	//	gridder.addSectorGrid(indshapefiles, "ind");
	//	gridder.addSectorGrid(resshapefiles, "res");
	//	gridder.addSectorGrid(railroadshapefiles, "railroad");
	//	gridder.addSectorGrid(onroadshapefiles, "onroad");
	//	gridder.addSectorGrid(nonroadshapefiles, "nonroad");
	//	gridder.addSectorGrid(elecprodshapefiles, "elecprod");
	//	gridder.addSectorGrid(marineshapefiles, "marine");
	//	gridder.addSectorGrid(airportshapefiles, "airport");
	//	gridder.loadAttribute(cafields[i]);
	//	//gridder.makeAnnualTotal("total_annual_" + year + ".nc");
	//	gridder.makeHourlyTotal("total_hourly_" + year + ".nc");
	//	//for (size_t i = 0; i < gridder.sectorGrids.size(); i++)
	//	//{
	//	//	HestiaGrid* outputGrid = gridder.sectorGrids[i];
	//	//	outputGrid->getTotal();
	//	//	outputGrid->toShapefile(subdirs[0] + "onroad.Shp", lidarFisheyeDir + year + "/" + outputGrid->sectorname + "_" + year + ".Shp");
	//	//}
	//	//gridder.toShapefile(subdirs[0] + "onroad.Shp", lidarFisheyeDir + year + "/" "total_" + year + ".Shp");
	//}



}

void updateAttributes(std::string indir, std::string outdir)
{
	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();

	std::vector<std::string> fields = Utils::buildVector("", new std::string[5]{ "ca10" , "ca11" , "ca12" , "ca13" , "ca14" }, 5);
	int basefield = 4;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp", Qt::CaseInsensitive))
			continue;
		if (fileInfo.fileName().endsWith("nonpoint.shp", Qt::CaseInsensitive))
			continue;
		if (fileInfo.fileName().endsWith("onroad.shp", Qt::CaseInsensitive))
			continue;
		std::string sectorName = fileInfo.completeBaseName().toLocal8Bit().data();
		std::string out_file = (QDir(outdir.data()).absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();

		ShapeFile shpin(input_file.data());

		OGRwkbGeometryType gtype = shpin.poLayer->GetGeomType();
		bool ispoint = false;
		if (gtype == wkbPoint || gtype == wkbMultiPoint || gtype == wkbPoint25D)
		{
			ispoint = true;
		}
		//if (!ispoint)
		//{
		//	shpin.close();
		//	continue;
		//}


		ShapeFile shpout(out_file.data(), 1);
		std::vector<int> fieldindicesSrc;
		std::vector<int> fieldindicesDest;
		for (size_t ifield = 0; ifield < 5; ifield++)
		{
			fieldindicesSrc.push_back(shpin.poLayer->GetLayerDefn()->GetFieldIndex(fields[ifield].data()));
			fieldindicesDest.push_back(shpout.getOrCreateField(fields[ifield].data(), OGRFieldType::OFTReal));
		}
		std::vector<std::map<int, double>> maps;
		std::map<int, double> footprintmap;
		maps.resize(fieldindicesSrc.size());
		OGRFeature *poFeature;
		int numchars = 6;
		if (sectorName.size() < numchars)
			numchars = sectorName.size();
		std::string destIDFieldName = "FID_" + sectorName.substr(0, numchars);
		int destIDField = shpout.poLayer->GetLayerDefn()->GetFieldIndex(destIDFieldName.data());
		int srcFID = 0;
		while ((poFeature = shpin.poLayer->GetNextFeature()) != NULL)
		{
			for (size_t ifield = 0; ifield < 5; ifield++)
			{
				int srcField = fieldindicesSrc[ifield];
				double ca = poFeature->GetFieldAsDouble(srcField);
				maps[srcField][srcFID] = ca;
			}
			/*	footprintmap[srcFID] = ca;*/
			srcFID++;
			OGRFeature::DestroyFeature(poFeature);
		}
		shpin.close();
		int srcbaseFieldIdx = fieldindicesSrc[basefield];
		int destbaseFieldIdx = fieldindicesDest[basefield];
		while ((poFeature = shpout.poLayer->GetNextFeature()) != NULL)
		{

			srcFID = poFeature->GetFieldAsInteger(destIDField);
			double srcBaseValue = maps[basefield][srcFID];
			double destBaseValue = poFeature->GetFieldAsDouble(destbaseFieldIdx);
			double frac = 1;
			if (!ispoint)
			{
				if (srcBaseValue == 0)
				{
					frac = 1;
				}
				else
				{
					frac = destBaseValue / srcBaseValue;
				}

			}

			for (size_t ifield = 0; ifield < 5; ifield++)
			{
				if (ifield == basefield)
					continue;
				int srcField = fieldindicesSrc[ifield];
				double srcValue = maps[srcField][srcFID];
				int destField = fieldindicesDest[ifield];
				double destValue = srcValue * frac;

				poFeature->SetField(destField, destValue);

			}
			shpout.poLayer->SetFeature(poFeature);
			OGRFeature::DestroyFeature(poFeature);
		}
		shpout.close();
	}


}

void gridLA_With_FFDAS()
{

	//Preprocessor::extractEDGASHourly("C:/FFDAS/aviation/2014.nc", "C:/FFDAS/WRF/aviation_2014.tif", adfTransform, 437, 212, 468, 284);
	//Preprocessor::extractEDGASHourly("C:/FFDAS/shipping/2014.nc", "C:/FFDAS/WRF/shipping_2014.tif", adfTransform, 437, 212, 468, 284);
	//Preprocessor::extractFFDASHourly("C:/FFDAS/2014/", "C:/FFDAS/FFDAS_SUBSET.tif", adfTransform, 437, 212, 468, 284);
	std::string rootdir = "C:/HestiaGridding/Los_Angeles/WRF_FFDAS_Gridding/";
	std::string timedir = "C:/HestiaGridding/Los_Angeles/Time/";
	std::vector<std::string> ffdasversions = Utils::buildVector("", new std::string[3]{ "WithoutShippingAviation" ,"WithShippingOnly","WithShippingAviation" }, 3);
	std::string outroot = "C:/HestiaGridding/Los_Angeles/WRF_FFDAS_Gridding/result/";
	std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[8]{ "ca10" , "ca11" , "ca12" , "ca13" , "ca14" ,"length","area" ,"timestruct" }, 8);
	std::string years[]{ "2010","2014" };
	std::string cafields[]{ "ca10","ca14" };
	std::string inputdirs[]{ rootdir + "grid01/intersected_nogeometry/",rootdir + "grid02/intersected_nogeometry/",rootdir + "grid03/intersected_nogeometry/" };
	std::string gridnames[]{ "geo_em.d01","geo_em.d02","geo_em.d03" };
	//std::string ffdasfile = "C:/FFDAS/WRF/FFDAS_SUBSET_With_Shipping.tif";
	//removeAviationFromFFDAS(ffdasfile.outlineDS(),Preprocessor::readData("C:/FFDAS/WRF/aviation_2014.tif"));

	//std::string indir = "B:/LA_Version2/ca/";
	//std::vector<std::string> subdirs = Utils::findSubdirectories(indir);


	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
	//	//Utils::updateFootprintForDir(subdirs[i],true);
	//	ShapeFile::copyDir(subdirs[i], "B:/LA_Version2/ca14/" + foldername + "/", fields2keep);
	//}


	//int startrow = 437;
	//int nrows = 212;
	//int startcol = 468;
	//int ncols = 284;
	//OGREnvelope ffdasbound;
	//ffdasbound.MinX = -180 + startcol * 0.1;
	//ffdasbound.MaxX = ffdasbound.MinX + ncols * 0.1;
	//ffdasbound.MaxY = 90 - startrow * 0.1;
	//ffdasbound.MinY = ffdasbound.MaxY - nrows * 0.1;
	//Grid grid;
	//grid.fromFishnetShape("B:/LA_Version2/Gridding/FFDASGrid.Shp");
	//grid._adfGeoTransform[0] = -180 + startcol * 0.1;
	//grid._adfGeoTransform[3] = 90 - startrow * 0.1;
	//grid._adfGeoTransform[1] = 0.1;
	//grid._adfGeoTransform[5] = -0.1;
	//grid.resetValue("ca14");
	//grid.toShape(&refshape, "B:/LA_Version2/ca14/FFDAS/ffdas.Shp",true);
	//Utils::updateFootprint("B:/LA_Version2/ca14/FFDAS/ffdas.Shp");
	//Preprocessor::updateFieldAfterIntersection("B:/LA_Version2/ca14/FFDAS/ffdas.Shp");

	//Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Los_Angeles/", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	//Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Orange/", "B:/LA_Version2/gridPrep_SHP_master/Orange/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	//Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Riverside/", "B:/LA_Version2/gridPrep_SHP_master/Riverside/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	//Preprocessor::reprojectDir("B:/LA_Version2/nonroad/San_Bernardino/", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");
	//Preprocessor::reprojectDir("B:/LA_Version2/nonroad/Ventura/", "B:/LA_Version2/gridPrep_SHP_master/Ventura/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectLANonRoad.py");


	//for (size_t i = 0; i < 10000; i++)
	//{
	//	if(i % 100 == 0)
	//	  printf("%f\n", ffdasAnnual[i]);
	//}
	for (size_t iyear = 0; iyear < 2; iyear++)
	{
		std::string year = years[iyear];
		for (int iffdasversion = 0; iffdasversion < 3; iffdasversion++)
		{
			std::string ffdasfile = outroot + year + "/" + "FFDAS_" + ffdasversions[iffdasversion] + ".tif";
			std::string outdir = outroot + year + "/" + ffdasversions[iffdasversion] + "/";
			if (!QDir(outdir.data()).exists())
				QDir(outdir.data()).mkpath(".");
			double* ffdasAnnual;
			std::vector<float*> ffdasHourly;
			for (size_t igriddomain = 0; igriddomain < 3; igriddomain++)
			{
				TemporalGridder gridder("", 8760, "meter");
				//gridder.fromFishnetRaster(indir + fishnetname.str() + ".tif");
				gridder.fromFishnetShapeFile(inputdirs[igriddomain] + "fishnet.shp");
				SubdirManager dirmanager = SubdirManager(inputdirs[igriddomain]);
				//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";
				std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
				std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
				std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
				std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "RailroadPoint","Railroad" }, 2));
				std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.shp" }, 1));
				std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "NonRoad" }, 1));
				std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
				std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
				std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));
				std::vector<std::string> ffdasshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ffdas" }, 1));
				std::string timestructfile = timedir + +"timestructs_" + year + ".bin";
				gridder.loadtimestruct("", timestructfile,8760);
				if (ffdasHourly.size() == 0)
				{
					ffdasHourly = loadFFDASSubset(ffdasfile.data(), ffdasAnnual);
				}
				gridder.addSectorGrid(comshapefiles, "com");
				gridder.addSectorGrid(indshapefiles, "ind");
				gridder.addSectorGrid(resshapefiles, "res");
				gridder.addSectorGrid(railroadshapefiles, "railroad");
				gridder.addSectorGrid(onroadshapefiles, "onroad");
				gridder.addSectorGrid(nonroadshapefiles, "nonroad");
				gridder.addSectorGrid(elecprodshapefiles, "elecprod");
				gridder.addSectorGrid(marineshapefiles, "marine");
				gridder.addSectorGrid(airportshapefiles, "airport");
				gridder.addFFDASGrid(ffdasshapefiles, "ffdas", ffdasAnnual, ffdasHourly);

				gridder.loadAttribute(cafields[iyear]);

				/*			if (!QDir((outdir + "grids/").outlineDS()).exists())
				QDir((outdir + "grids/").outlineDS()).mkpath(".");
				if (!QDir((outdir + "shapes/").outlineDS()).exists())
				QDir((outdir + "shapes/").outlineDS()).mkpath(".");*/
				gridder.makeAnnualTotal(outdir + gridnames[igriddomain] + "." + "total.annual." + year + ".nc");
				gridder.toShapefile(inputdirs[igriddomain] + "fishnet.shp", outdir + gridnames[igriddomain] + "." + "total.annual." + year + ".shp");
				gridder.makeHourlyTotal(outdir + gridnames[igriddomain] + "." + "total.hourly." + year + ".nc");
				//for (size_t i = 0; i < gridder.sectorGrids.size(); i++)
				//{
				//	HestiaGrid* outputGrid = gridder.sectorGrids[i];
				//	outputGrid->getTotal();
				//	outputGrid->toShapefile(subdirs[0] + "onroad.Shp", lidarFisheyeDir + year + "/" + outputGrid->sectorname + "_" + year + ".Shp");
				//}
			}
			for (size_t iffdashours = 0; iffdashours < ffdasHourly.size(); iffdashours++)
			{
				delete[] ffdasHourly[iffdashours];
			}
			delete ffdasAnnual;
		}

	}


}
void testTimeshift()
{
	TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/time/2011/ComNonPoint.bin", "B:/LA_Version2/Vulcan_output/time/2011/ComNonPoint/", ".csv");
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/time/2015/ResNonPoint.bin", "B:/LA_Version2/Vulcan_output/time/2015/ResNonPoint/", ".csv");
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/time/2015/Airport.bin", "B:/LA_Version2/Vulcan_output/time/2015/Airport/", ".csv");
	return;
	//std::string id;
	//double* srcfracs = TimestructTool::readText("B:/LA_Version2/Time/2011/560375237.txt",id,true);
	//TimestructTool::toText("C:/HestiaGridding/560375237_2011.csv", srcfracs);
	//double* destfracs = TimestructTool::timeshift(srcfracs, 2011, 2013);
	//TimestructTool::normalizeFractions(destfracs);
	//TimestructTool::toText("C:/HestiaGridding/560375237_2013.csv", destfracs);
	//5603718112

	//TimestructTool::updateBinary("C:/HestiaGridding/Los_Angeles/2012timestructs.bin", "C:/HestiaGridding/Los_Angeles/2012timestructs.bin", "B:/LA_Version2/Time/2012/");
	//TimestructTool::selectFromBinary2Text("C:/HestiaGridding/Los_Angeles/2011timestructs.bin", "4334", "C:/HestiaGridding/Los_Angeles/4334_2011.csv");
	//TimestructTool::selectFromBinary2Text("C:/HestiaGridding/Los_Angeles/2012timestructs.bin", "4334", "C:/HestiaGridding/Los_Angeles/4334_2012.csv");


	int srcyear = 2011;
	int destyears[]{ 2010,2012,2013,2014,2015 };
	std::string indir = "B:/LA_Version2/Vulcan_output/Time/";
	//std::string indir = "B:/LA_Version2/Time/";
	//std::string testid = "560375237";
	//std::string testid = "4335";
	//TimestructTool::txt2binary("B:/LA_Version2/Time/2011/", "C:/HestiaGridding/Los_Angeles/2011timestructs.bin",true);
	//TimestructTool::selectFromBinary2Text("C:/HestiaGridding/Los_Angeles/2011timestructs.bin", testid, "C:/HestiaGridding/Los_Angeles/" + testid + "_2011.csv");
	//}
	//TimestructTool::txt2binary("B:/Baltimore/Time2011/", "C:/HestiaGridding/Baltimore/2011timestructs.bin");
	//TimestructTool::binary2Text("C:/HestiaGridding/Baltimore/2011timestructs.bin", "C:/HestiaGridding/Baltimore/2011/",".csv");
	std::string baseyeardir = indir + "2011/";
	std::vector<std::string> subdirs;
	QDir rootdir(baseyeardir.data());
	rootdir.setFilter(QDir::Dirs | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	rootdir.setSorting(QDir::Name);
	std::vector<std::string> files;
	QFileInfoList list = rootdir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string dirname = (fileInfo.baseName()).toLocal8Bit().data();
		subdirs.push_back(dirname);
		std::string srcbinaryfile = baseyeardir + dirname + ".bin";
		if (!QFileInfo(srcbinaryfile.data()).exists())
		{
			TimestructTool::txt2binary(8760, baseyeardir + dirname + "/", srcbinaryfile,  true);
		}
		for (int nyear = 0; nyear < 5; nyear++)
		{
			int destyear = destyears[nyear];
			int numhours = TimestructTool::numOfHours(destyear);
			std::stringstream destyeardir;
			destyeardir << indir << destyear << "/";
			std::string destbinaryfile = destyeardir.str() + dirname + ".bin";
			std::string destdir = destyeardir.str() + dirname + "/";

			if (QFileInfo(destbinaryfile.data()).exists())
			{
				continue;
			}
			else if (QDir(destdir.data()).exists())
			{
				TimestructTool::txt2binary(numhours, destdir, destbinaryfile, true);
			}
			if (!QFileInfo(destbinaryfile.data()).exists()) {
				TimestructTool::timeshiftbinary(srcbinaryfile, destbinaryfile, srcyear, destyear, true);
			}
			//TimestructTool::selectFromBinary2Text(destbinaryfiless.str(), testid, testfiless.str());
		}
	}


	//for (size_t i = 0; i < 4; i++)
	//{
	//	int destyear = destyears[i];
	//	std::stringstream srcbinaryfiless;
	//	srcbinaryfiless << outdir << srcyear << "timestructs.bin";
	//	std::stringstream destbinaryfiless; http://www.foxnews.com/science/2016/09/08/nasa-spacecraft-on-way-to-asteroid-to-bring-back-samples.html
	//	destbinaryfiless << outdir << destyear << "timestructs.bin";

	//	std::stringstream destupdatedirss;
	//	destupdatedirss << indir << destyear << "/";
	//	//TimestructTool::timeshiftbinary(srcbinaryfiless.str(), destbinaryfiless.str(), srcyear, destyear);
	//	//TimestructTool::updateBinary(destbinaryfiless.str(), destbinaryfiless.str(), destupdatedirss.str());
	//	std::stringstream testfiless;
	//	testfiless << outdir << testid << "_" << destyear << ".csv";
	//	TimestructTool::selectFromBinary2Text(destbinaryfiless.str(), testid, testfiless.str());
	//}

	//	TimestructTool::timeshiftbinary("C:/HestiaGridding/Los_Angeles/2011timestructs.bin", "C:/HestiaGridding/Los_Angeles/2012timestructs.bin", 2011, 2012);
}
void reallocateLA()
{
	std::string outdir = "E:/LASpatialAnalysis/BlockGroup_Intersected/";
	std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	//std::string rootdir = "B:/LA_Version2/gridPrep_SHP_master/";
	//std::vector<std::string> subdirs = Utils::findSubdirectories(rootdir);
	//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[3]{ "ca10","length","area"}, 3);
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
	//	ShapeFile::copyDir(subdirs[i], indir + foldername + "/", fields2keep);
	//}

	std::vector<std::string> subdirs = Utils::findSubdirectories(indir);
	std::string intersectionShapefile = "E:/LASpatialAnalysis/CensusShapes/pop_income_bg.shp";
	//Utils::updateFootprint(intersectionShapefile);
	for (size_t i = 0; i < subdirs.size(); i++)
	{
		std::string foldername = QDir(subdirs[i].data()).dirName().toLocal8Bit().data();
		std::string srcdir = subdirs[i];
		std::string destdir = outdir + foldername + "/";
		Preprocessor::gridFolderByShape(srcdir, destdir, intersectionShapefile);;
	}

	indir = "E:/LASpatialAnalysis/ca10_BlockGroup_Intersected/";
	outdir = "E:/LASpatialAnalysis/ca10_BlockGroup/";
	subdirs = Utils::findSubdirectories(indir);

	SubdirManager dirmanager(indir);
	subdirs = Utils::findSubdirectories(indir);

	QDir(outdir.data()).mkpath(".");

	Reallocator allocator;
	std::string timestructfile = "";// lidarFisheyeDir + +"timestructs_" + year + ".bin";
									//gridder.loadtimestruct("B:/LA_Version2/Time/" + year + "/", timestructfile);
	dirmanager = SubdirManager(indir);
	//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca10/" + resolname.str() + "/";
	std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "RailroadPoint","Railroad" }, 2));
	std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.shp" }, 1));
	std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "NonRoad" }, 1));
	std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));

	allocator.addSector(comshapefiles, "com");
	allocator.addSector(indshapefiles, "ind");
	allocator.addSector(resshapefiles, "res");
	allocator.addSector(railroadshapefiles, "railroad");
	allocator.addSector(onroadshapefiles, "onroad");
	allocator.addSector(nonroadshapefiles, "nonroad");
	allocator.addSector(elecprodshapefiles, "elecprod");
	allocator.addSector(marineshapefiles, "marine");
	allocator.addSector(airportshapefiles, "airport");
	std::vector<std::string> idfields;
	idfields.push_back("TRACT"); idfields.push_back("BLKGRP");
	allocator.loadAttribute("ca10", idfields);
	allocator.aggregrate(intersectionShapefile, outdir + "total.shp");

}
void reallocateHestiaFFDAS()
{
	//indir = "E:/LASpatialAnalysis/ca11_BlockGroup_Intersected/";
	//lidarFisheyeDir = "E:/LASpatialAnalysis/ca11_BlockGroup/";
	Utils::updateFootprintForDir("E:/LASpatialAnalysis/ca10_BlockGroup_Intersected/");
	std::string intersectionShapefile = "E:/LASpatialAnalysis/CensusShapes/pop_income_bg.shp";
	//Utils::updateFootprintForDir("E:/LASpatialAnalysis/FFDAS_Hestia/");
	//Preprocessor::gridFolderByShape("E:/LASpatialAnalysis/FFDAS_Hestia/", "E:/LASpatialAnalysis/ca11_CensusTract_Intersected/FFDAS_Hestia/", intersectionShapefile);

	std::string indir = "E:/LASpatialAnalysis/ca10_BlockGroup_Intersected/FFDAS_Hestia/";
	std::string outdir = "E:/LASpatialAnalysis/ca10_BlockGroup/";

	SubdirManager dirmanager(indir);

	QDir(outdir.data()).mkpath(".");

	Reallocator allocator;
	std::string timestructfile = "";// lidarFisheyeDir + +"timestructs_" + year + ".bin";
									//gridder.loadtimestruct("B:/LA_Version2/Time/" + year + "/", timestructfile);
									//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";

	std::vector<std::string> ffdasshapes = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "FFDAS2011" }, 1));
	std::vector<std::string> odiacshapes = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ODIAC2011" }, 1));
	std::vector<std::string> hestiashapes = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "Hestia2011" }, 1));

	allocator.addSector(ffdasshapes, "ffdas");
	allocator.addSector(odiacshapes, "odiac");
	allocator.addSector(hestiashapes, "hestia");
	std::vector<std::string> idfields;
	idfields.push_back("TRACT"); idfields.push_back("BLKGRP");
	allocator.loadAttribute("ca10", idfields);
	for (size_t i = 0; i < allocator.sectors.size(); i++)
	{
		allocator.sectors[i]->aggregrate(intersectionShapefile, outdir + allocator.sectors[i]->sectorname + ".shp");
	}

}
void tif2csv(std::string tiffile, std::string outfile)
{
	int xsize, ysize;
	std::vector<float*> data = Preprocessor::readData(tiffile.data(), xsize, ysize);
	int num = xsize * ysize;
	std::ofstream ofs;
	ofs.open(outfile);
	float* pdata = data[0];
	for (size_t i = 0; i < num; i++)
	{
		if (pdata[i] < 0)
			pdata[i] = 0;
		ofs << pdata[i] << std::endl;
	}
	ofs.close();

	for (size_t i = 0; i < data.size(); i++)
	{
		delete[] data[i];
	}
}
void replaceNoData(std::string infile, std::string outfile)
{
	int xsize, ysize;

	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(infile.data(), GA_ReadOnly);


	xsize = poDataset->GetRasterXSize();
	ysize = poDataset->GetRasterYSize();
	double adftransform[6];
	poDataset->GetGeoTransform(adftransform);
	std::string proj = poDataset->GetProjectionRef();
	float* data = new float[xsize * ysize];
	int   ncells = xsize*ysize;

	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
		data, xsize, ysize, GDT_Float32,
		0, 0);

	GDALClose((GDALDatasetH)poDataset);
	int num = xsize * ysize;
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset* pDataset = poDriver->Create(outfile.data(), xsize, ysize, 1, GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(adftransform);
	pDataset->SetProjection(proj.data());
	float _nodata = 0;
	for (size_t i = 0; i < num; i++)
	{
		double val = data[i];
		if (val <= 0)
			data[i] = _nodata;
	}

	GDALRasterBand *pBand = pDataset->GetRasterBand(1);
	pBand->RasterIO(GF_Write, 0, 0, xsize, ysize, data, xsize, ysize, GDT_Float32, 0, 0);
	pBand->SetNoDataValue(_nodata);

	GDALClose((GDALDatasetH)pDataset);

	delete[] data;
}
void calIntensity(std::string indir)
{
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp", Qt::CaseInsensitive))
			continue;
		ShapeFile file(fileInfo.absoluteFilePath().toLocal8Bit().data(), 1);
		OGRFeature *poFeature;
		int caidx = file.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
		int areaidx = file.poLayer->GetLayerDefn()->GetFieldIndex("area");
		int intensityidx = file.getOrCreateField("intensity", OGRFieldType::OFTReal);
		file.poLayer->ResetReading();
		while ((poFeature = file.poLayer->GetNextFeature()) != NULL)
		{
			double area = poFeature->GetFieldAsDouble(areaidx);
			double ca = poFeature->GetFieldAsDouble(caidx);
			double intensity = 0;
			if (ca > 0 && area > 0)
			{
				intensity = ca / area;
			}
			poFeature->SetField(intensityidx, intensity);
			file.poLayer->SetFeature(poFeature);
			OGRFeature::DestroyFeature(poFeature);
		}
	}


}


void updateOnRoadFields()
{
	SubdirManager dirmanager = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/");
	std::vector<std::string> onroadfiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "onroad.shp" }, 1));
	std::vector<std::string> oldnames = Utils::buildVector("", new std::string[5]{ "ca10_p","ca11_p","ca12_p","ca13_p","ca14_p" }, 5);
	std::vector<std::string> newnames = Utils::buildVector("", new std::string[5]{ "ca10_g","ca11_g","ca12_g","ca13_g","ca14_g" }, 5);


	//for (size_t i = 0; i < onroadfiles.size(); i++)
	//{
	//	ShapeFile::changeFieldNames(onroadfiles[i], oldnames, newnames);
	//	//ShapeFile::copyDropFields();
	//}
	std::vector<std::string> subdirs = Utils::findSubdirectories("B:/LA_Version2/gridPrep_SHP_master/");
	std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[16]{
		"ca10","ca11","ca12","ca13","ca14",
		"ca10_g","ca11_g","ca12_g","ca13_g","ca14_g",
		"ca10_d","ca11_d","ca12_d","ca13_d","ca14_d","timestruct" }, 16);
	for (size_t i = 0; i < subdirs.size(); i++)
	{
		std::string foldername = QDir(subdirs[i].data()).dirName().toLocal8Bit().data();
		std::string srcFile = subdirs[i] + "OnRoad.shp";
		std::string destFile = "B:/LA_Version2/Public/" + foldername + "/" + "OnRoad.shp";
		ShapeFile::copy(srcFile, destFile, fields2keep);

		//Utils::updateFootprintForDir(subdirs[i],false);
		//ShapeFile::copyDir(subdirs[i], "B:/LA_Version2/ca/" + foldername + "/", fields2keep);
		//Preprocessor::gridFolderByRaster("B:/LA_Version2/ca/" + foldername + "/", indir + foldername + "/", "B:/LA_Version2/gridoutput/fishnet1000m.tif");
	}
}
#include "EMFAC_FuelSplit.h"


void plotRaster()
{
	//
	std::string indir = "C:/HestiaGridding/Los_Angeles/2011/";
	std::vector<std::string> shapefiles = Utils::findFiles(indir, ".shp");
	//infilename = args[1]
	//fieldname = args[2]
	//figurename = args[3]
	//xtitle = args[4]
	//min = args[5]
	//max = args[6]
	for (size_t i = 0; i < shapefiles.size(); i++)
	{
		std::string name = QFileInfo(shapefiles[i].data()).baseName().toLocal8Bit().data();
		std::stringstream ss;
		std::string filename = indir + name + ".dbf";
		ss << "RScript B:/LA_Version2/histogramDBF.R" << " " << filename << " " << "ca11" << " " << indir + name + ".png" << " " << "\"" << "Natural Log(X)[Carbon = e^X(KgC)]" << "\"" << " " << name;
		system(ss.str().data());
		//"B:/LA_Version2/histogramDBF.R"
	}

}

void plotShapes()
{
	//
	SubdirManager dirmanager = SubdirManager("B:/LA_Version2/Public/");
	std::vector<std::string> shapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ ".shp" }, 1));
	//infilename = args[1]
	//fieldname = args[2]
	//figurename = args[3]
	//xtitle = args[4]
	//min = args[5]
	//max = args[6]
	for (size_t i = 0; i < shapefiles.size(); i++)
	{
		std::string name = QFileInfo(shapefiles[i].data()).baseName().toLocal8Bit().data();
		std::string indir = QFileInfo(shapefiles[i].data()).absoluteDir().absolutePath().toLocal8Bit().data() + std::string("/");
		std::stringstream ss;
		std::string filename = indir + name + ".dbf";
		ss << "RScript B:/LA_Version2/histogramDBF.R" << " " << filename << " " << "ca11" << " " << indir + name + ".png" << " " << "\"" << "Natural Log(X)[Carbon = e^X(KgC)]" << "\"" << " " << name;
		system(ss.str().data());
		//"B:/LA_Version2/histogramDBF.R"
	}

}
void gridBaltimore()
{
	int resol = 200;
	std::stringstream fishnetname;
	fishnetname << "fishnet" << resol << "m";
	std::stringstream resolname;
	resolname << resol << "m";


	//Grid grid(BoundManager::readBound("B:/Baltimore/gridPrep_SHP_master/StatePlane/bound.txt"), resol * 3.28084, 1);
	//ShapeFile refshape("B:/Baltimore/gridPrep_SHP_master/StatePlane/nonroad.Shp");
	//grid.toShape(&refshape, "B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".Shp");
	//grid.reset();
	//grid.toRaster("B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");
	//Utils::updateFootprintForDir("B:/Baltimore/gridPrep_SHP_master/ca/");
	//Preprocessor::gridFolderByRaster("B:/Baltimore/gridPrep_SHP_master/ca_1.1/", "B:/Baltimore/gridPrep_SHP_master/ca_1.1/" + resolname.str() + "/", "B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");


	std::string timestruct2001file = "C:/HestiaGridding/Baltimore/1.1/Baltimore.timestructs.2011.v1.1.bin";
	std::string years[]{ "2010","2011","2012","2013","2014" };
	std::string cafields[]{ "ca10","ca11","ca12","ca13","ca14" };
	std::string cityname = "Baltimore";
	std::string version = "v1.1";
	//Baltimore.total.hourly.2010.v1.1.nc
	std::string outroot = "C:/HestiaGridding/Baltimore/1.1/";
	std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca_1.1/" + resolname.str() + "/";
	for (size_t i = 1; i < 2; i++)
	{
		std::string year = years[i];
		TemporalGridder gridder(outroot,8760);
		gridder.fromFishnetRaster("B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");
		std::string timestructfile = outroot + cityname + ".timestructs." + year + "." + version + ".bin";

		if (!QFileInfo(timestructfile.data()).exists())
		{
			//TimestructTool::timeshiftbinary(timestruct2001file.outlineDS(), timestructfile.outlineDS(),2011, atoi(year.outlineDS()),true);
		}

		//gridder.loadtimestruct("B:/Baltimore/Time2011/", timestructfile);
		std::vector<std::string> comshapefiles = constructShapeFilelist(indir, new std::string[2]{ "comnonpoint","compoint" }, 2);
		std::vector<std::string> indshapefiles = constructShapeFilelist(indir, new std::string[2]{ "indnonpoint","indpoint" }, 2);
		std::vector<std::string> resshapefiles = constructShapeFilelist(indir, new std::string[1]{ "resnonpoint" }, 1);
		std::vector<std::string> railroadshapefiles = constructShapeFilelist(indir, new std::string[2]{ "RailroadPoint","Railroad" }, 2);
		std::vector<std::string> onroadshapefiles = constructShapeFilelist(indir, new std::string[1]{ "OnRoad" }, 1);
		std::vector<std::string> nonroadshapefiles = constructShapeFilelist(indir, new std::string[1]{ "NonRoad" }, 1);
		std::vector<std::string> elecprodshapefiles = constructShapeFilelist(indir, new std::string[1]{ "ElecProd" }, 1);
		std::vector<std::string> marineshapefiles = constructShapeFilelist(indir, new std::string[2]{ "CMVPort","CMVUnderway" }, 2);
		gridder.addSectorGrid(comshapefiles, "com");
		gridder.addSectorGrid(indshapefiles, "ind");
		gridder.addSectorGrid(resshapefiles, "res");
		gridder.addSectorGrid(railroadshapefiles, "railroad");
		gridder.addSectorGrid(onroadshapefiles, "onroad");
		gridder.addSectorGrid(nonroadshapefiles, "nonroad");
		gridder.addSectorGrid(elecprodshapefiles, "elecprod");
		gridder.addSectorGrid(marineshapefiles, "marine");
		gridder.loadAttribute(cafields[i]);

		gridder.makeAnnualTotal(outroot + cityname + ".total.annual." + year + "." + version + ".nc");
		gridder.toShapefile(onroadshapefiles[0], outroot + cityname + ".total.annual." + year + "." + version + ".shp");
		//gridder.makeHourlyTotal(outroot + cityname + ".total.hourly." + year + "." + version + ".nc");
		for (size_t isector = 0; isector < gridder.sectorGrids.size(); isector++)
		{
			HestiaGrid* sectorGrid = gridder.sectorGrids[isector];
			sectorGrid->getTotal();
			//outputGrid->toShapefile(indir + "onroad.Shp", lidarFisheyeDir + outputGrid->sectorname + "_2011.Shp");
			sectorGrid->toShapefile(indir + "onroad.shp", outroot + cityname + "." + sectorGrid->sectorname + "." + year + "." + version + ".shp");
		}
		//printf("%f,%f\n", checkTotal(outroot + cityname + ".total.annual." + year + "." + version + ".nc"), checkTotal(outroot + cityname + ".total.hourly." + year + "." + version + ".nc"));

	}


	//gridder.toShapefile(indir + "onroad.Shp", lidarFisheyeDir + "total_2011.Shp");
}
void gridIndy()
{
	std::string outroot = "C:/HestiaGridding/Indy/";
	std::string inroot = "B:/Indianapolis/";
	std::string year = "2013";

	std::string indirshapes = inroot + "y" + year + "/Shapes/";
	std::string timestructfile = outroot + "timestructs_" + year + ".bin";
	std::string timestructdir = inroot + "y" + year + "/Time/";
	std::string cafield = "CarbonKG" + year.substr(2, 2);



	TemporalGridder gridder(outroot, 8760,"degree");
	gridder.fromFishnetShapeFile("B:/Indianapolis/Grids/g2d_WGS84.shp");
	//gridder.loadtimestruct(timestructdir, timestructfile);
	std::vector<std::string> allsectors = Utils::findFiles(indirshapes, ".shp");
	gridder.addSectorGrid(allsectors, "allsectors");
	gridder.loadAttribute(cafield);
	gridder.makeAnnualTotal(outroot + "total_annual_" + year + ".nc");
	//gridder.makeDailyTotal(outroot + "total_daily_" + year + ".nc");
	gridder.toShapefile("B:/Indianapolis/Grids/g2d_WGS84.shp", outroot + "total_" + year + ".shp");

	//for (size_t i = 0; i < gridder.sectorGrids.size(); i++)
	//{
	//	HestiaGrid* outputGrid = gridder.sectorGrids[i];
	//	outputGrid->getTotal();
	//	outputGrid->toShapefile(indir + "onroad.Shp", lidarFisheyeDir + outputGrid->sectorname + "_2011.Shp");
	//}
	//gridder.toShapefile(indir + "onroad.Shp", lidarFisheyeDir + "total_2011.Shp");
}


void gridLA3()
{
	std::string indir = "H:/HestiaGridding/Los_Angeles_Vulcan/Shapes/";
	std::string outdir = "H:/HestiaGridding/Los_Angeles_Vulcan/Gridded/";
	std::string fishnetrasterfile = "H:/HestiaGridding/Los_Angeles_Vulcan/fishnet.tif";
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
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
	}

	std::string outfishnetfile = outdir + "fishnet.shp";
	std::string outrasterfile = outdir + "fishnet.tif";

	std::string fishnetfile = (QFileInfo(fishnetrasterfile.data()).absoluteDir().absolutePath() + "/fishnet.shp").toLocal8Bit().data();
	std::string rasterfile = fishnetrasterfile;

	Grid fishnet;
	fishnet.fromFishnetRaster(fishnetrasterfile);
	fishnet.reset();

	GDALDataset* pDataset = (GDALDataset*)GDALOpen(fishnetrasterfile.data(), GA_ReadOnly);
	std::string wkt = pDataset->GetProjectionRef();
	GDALClose(pDataset);
	if (!QFileInfo(fishnetfile.data()).exists())
	{
		fishnet.toShape(wkt, fishnetfile, false);
	}
	for (size_t i = 0; i < files.size(); i++)
	{
		/*if (skipNonRoad && files[i] == "NonRoad.shp")
		continue;*/
		printf("%s\n", (outdir + files[i]).data());
		//intersectWithFishnet(fishnet, fishnetfile, indir + files[i], outdir + files[i]);

		if (!QFileInfo((outdir + files[i]).data()).exists())
		{
			Utils::updateFootprint(indir + files[i]);
			Preprocessor::intersectWithArcGIS(indir + files[i], fishnetfile, outdir + files[i]);
			while (!QFileInfo((outdir + files[i]).data()).exists())
			{
				printf("ʧ��!");
				Preprocessor::intersectWithArcGIS(indir + files[i], fishnetfile, outdir + files[i]);
			}
			updateFieldAfterIntersection(outdir + files[i]);
		}

		//ShapeFile input((outdir + files[i]).data());
		//fishnet.gatherCells(&input);
	}

	for (size_t i = 0; i < files.size(); i++)
	{
		fishnet.reset();
		ShapeFile input((outdir + files[i]).data());

		std::string sectorfile = files[i].substr(0, files[i].length() - 4);
		std::string sectorfileout = outdir + sectorfile + ".tif";
		printf("%s\n", sectorfile.data());
		fishnet.gatherCells(&input);
		fishnet.toRaster(sectorfileout);
		sectorfileout = outdir + sectorfile + ".tif";
		fishnet.toRaster(sectorfileout);
	}

}
//void gridMarion()
//{
//
//
//	TemporalGridder gridder("", "feet");
//	gridder.fromFishnetRaster("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridPrep_SHP_intersect/grids/Airport.tif");
//	std::string indir = "Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridPrep_SHP_intersect/grids/";
//	gridder.tif2netcdf(indir + "AllSectors.tif", indir + "AllSectors.nc");
//	gridder.tif2netcdf(indir + "Airport.tif", indir + "Airport.nc");
//	gridder.tif2netcdf(indir + "AllSectors.tif", indir + "AllSectors.nc");
//	gridder.tif2netcdf(indir + "ElecProd.tif", indir + "ElecProd.nc");
//	gridder.tif2netcdf(indir + "NonRoad.tif", indir + "NonRoad.nc");
//	gridder.tif2netcdf(indir + "ResNonPoint.tif", indir + "Res.nc");
//	gridder.tif2netcdf(indir + "OnRoad.tif", indir + "OnRoad.nc");
//	std::vector<std::string> comfiles;
//	comfiles.push_back(indir + "ComNonPoint.tif");
//	comfiles.push_back(indir + "ComPoint.tif");
//	gridder.tif2netcdf(comfiles, indir + "Com.nc");
//	std::vector<std::string> indfiles;
//	comfiles.push_back(indir + "IndNonPoint.tif");
//	comfiles.push_back(indir + "IndPoint.tif");
//	gridder.tif2netcdf(indfiles, indir + "Ind.nc");
//
//	std::vector<std::string> railroadfiles;
//	comfiles.push_back(indir + "Railroad.tif");
//	comfiles.push_back(indir + "RailroadPoint.tif");
//	gridder.tif2netcdf(indfiles, indir + "Railroad.nc");
//}




void gridMarion()
{
	double resol = 200;
	resol = resol * 3.28084;
	Grid grid(BoundManager::readBound("H:/Marion/bound.txt"), resol, 1);
	//
	if (!QFileInfo("B:/Marion/fishnet.tif").exists()) {
		grid.toRaster("B:/Marion/fishnet.tif");
	}
	if (!QFileInfo("B:/Marion/fishnet.shp").exists()) {
		grid.toShape(grid.proj, "B:/Marion/fishnet.shp");
	}
	//######################################################################
	std::string rootdir = "B:/Marion/";
	std::string outdir = rootdir + "GriddedEmissions/";
	std::string shapesindir = "B:/Indianapolis/Marion/Vulcan_output/Spatial/";
	std::string intersected = rootdir + "intersected_tmp/";
	std::string indir = rootdir + "intersected/";
	std::string timedir = "B:/Indianapolis/Time/";
	std::string gridShapeFile = "B:/Marion/fishnet.shp";
	std::string fishnetraster = "B:/Marion/fishnet.tif";


	QDir(outdir.data()).mkdir(".");
	QDir(shapesindir.data()).mkdir(".");
	QDir(intersected.data()).mkdir(".");
	QDir(indir.data()).mkdir(".");

	//Utils::updateFootprintForDir(shapesindir, true);
	Preprocessor::gridFolderByShape(shapesindir, intersected, gridShapeFile);
	std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"ca15","timestruct" }, 7);
	ShapeFile::copyDirDropGeometry(intersected, indir, fields2keep);

	std::string years[]{ "2010","2011","2012","2013","2014","2015" };
	std::string cafields[]{ "ca10","ca11","ca12","ca13","ca14","ca15" };
	std::map<std::string, std::vector<std::string>> timestructfile_crosswalk;

	std::vector<GridSectorConfig> sectors;
	GridSectorConfig com("commercial", Utils::buildVector(new std::string[1]{ "comnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));
	GridSectorConfig res("residential", Utils::buildVector(new std::string[1]{ "resnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));
	GridSectorConfig ind("industrial", Utils::buildVector(new std::string[1]{ "indnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));
	GridSectorConfig onroad("onroad", Utils::buildVector(new std::string[1]{ "onroad" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));

	//sectors.push_back(com);
	//sectors.push_back(res);
	//sectors.push_back(ind);
	sectors.push_back(onroad);
	QDir(outdir.data()).mkpath(".");
	std::string cityname = "";
	std::string version = "";
	//Grid grid;
	grid.fromFishnetRaster(fishnetraster);
	for (int i = 0; i < 6; i++)
	{

		std::string year = years[i];
		std::string attName = "ca" + year.substr(2, 2);
		//std::string timestructdir = timedir + year + "/";
		std::string timestructdir = timedir;
		//TimestructTool::normalizeBinary(timestructfile);
		SubdirManager dirmanager(indir);
		for (int nsector = 0; nsector < sectors.size(); nsector++)
		{
			grid.reset(attName);
			int numhous = 8760;
			if (year == "2012")
			{
				numhous = 8784;
			}
			TemporalGridder gridder(outdir + year + "/", numhous, "deg");
			gridder.fromFishnetRaster(fishnetraster);

			std::vector<std::string> shapefiles = dirmanager.findFilesMatch(sectors[nsector].shapefiles);
			HestiaGrid* sectorGrid = gridder.addSectorGrid(shapefiles, sectors[nsector].sectorname);
			sectorGrid->loadtimestruct(timestructdir + year + ".bin");
			gridder.loadAttribute(cafields[i]);
			printf("%s\n", std::string(outdir + sectorGrid->sectorname + "." + year + ".annual"  + ".shp").data());
			sectorGrid->getTotal();
			sectorGrid->toShapefile(shapesindir + "ResNonPoint.shp", outdir + sectorGrid->sectorname + "." + year + ".annual" + ".shp");
			grid.reset(attName);
			grid.gatherCells(outdir + sectorGrid->sectorname + "." + year + ".annual" + + ".shp", attName.data());
			grid.toRaster(outdir + sectorGrid->sectorname + "." + year + ".annual" + ".tif");
			gridder.makeAnnualTotal(outdir + sectorGrid->sectorname + "." + year + ".annual" + ".nc");
			gridder.makeHourlyTotal(outdir + sectorGrid->sectorname + "." + year + ".hourly" + ".nc");
			
		}
	}
	//TemporalGridder gridder("", "feet");
	//gridder.fromFishnetRaster("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridPrep_SHP_intersect/grids/Airport.tif");
	//int resol = 200;
	//std::stringstream fishnetname;
	//fishnetname << "fishnet" << resol << "m";
	//std::stringstream resolname;
	//resolname << resol << "m";


	//Grid grid(BoundManager::readBound("B:/Baltimore/gridPrep_SHP_master/StatePlane/bound.txt"), resol * 3.28084, 1);
	//ShapeFile refshape("B:/Baltimore/gridPrep_SHP_master/StatePlane/nonroad.Shp");
	//grid.toShape(&refshape, "B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".Shp");
	//grid.reset();
	//grid.toRaster("B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");
	//Utils::updateFootprintForDir("B:/Baltimore/gridPrep_SHP_master/ca/");
	//


	//std::string timestruct2001file = "C:/HestiaGridding/Baltimore/1.1/Baltimore.timestructs.2011.v1.1.bin";
	//std::string years[]{ "2010","2011","2012","2013","2014" };
	//std::string cafields[]{ "ca10","ca11","ca12","ca13","ca14" };
	//std::string cityname = "Marion";
	//std::string version = "v3.0";
	////Baltimore.total.hourly.2010.v1.1.nc
	//std::string outroot = "Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/";
	//std::string indir = "B:/Marion/2011/wgs84/";
	//std::string fishnetshpfile = "B:/Indianapolis/Grids/g2d_WGS84.shp";
	////Grid grid;
	////grid.fromFishnetShape(fishnetshpfile);
	//std::string fishnetrasterfile = "B:/Indianapolis/Grids/g2d_WGS84.tif";
	////grid.toRaster(fishnetrasterfile);

	////Preprocessor::gridFolderByRaster("B:/Marion/2011/", "B:/Marion/2011/wgs84/", fishnetrasterfile);

	//for (size_t i = 0; i < 5; i++)
	//{
	//	std::string year = years[i];
	//	TemporalGridder gridder(outroot, 8760, "degrees");
	//	gridder.fromFishnetRaster(fishnetrasterfile);
	//	//std::string timestructfile = outroot + cityname + ".timestructs." + year + "." + version + ".bin";

	//	//if (!QFileInfo(timestructfile.outlineDS()).exists())
	//	//{
	//	//	//TimestructTool::timeshiftbinary(timestruct2001file.outlineDS(), timestructfile.outlineDS(),2011, atoi(year.outlineDS()),true);
	//	//}

	//	//gridder.loadtimestruct("B:/Baltimore/Time2011/", timestructfile);
	//	std::vector<std::string> comshapefiles = constructShapeFilelist(indir, new std::string[2]{ "comnonpoint","compoint" }, 2);
	//	std::vector<std::string> indshapefiles = constructShapeFilelist(indir, new std::string[2]{ "indnonpoint","indpoint" }, 2);
	//	std::vector<std::string> resshapefiles = constructShapeFilelist(indir, new std::string[1]{ "resnonpoint" }, 1);
	//	std::vector<std::string> railroadshapefiles = constructShapeFilelist(indir, new std::string[2]{ "RailroadPoint","Railroad" }, 2);
	//	std::vector<std::string> onroadshapefiles = constructShapeFilelist(indir, new std::string[1]{ "OnRoad" }, 1);
	//	std::vector<std::string> airportshapefiles = constructShapeFilelist(indir, new std::string[1]{ "Airport" }, 1);
	//	std::vector<std::string> nonroadshapefiles = constructShapeFilelist(indir, new std::string[1]{ "NonRoad" }, 1);
	//	std::vector<std::string> elecprodshapefiles = constructShapeFilelist(indir, new std::string[1]{ "ElecProd" }, 1);
	//	//std::vector<std::string> marineshapefiles = constructShapeFilelist(indir, new std::string[2]{ "CMVPort","CMVUnderway" }, 2);
	//	gridder.addSectorGrid(airportshapefiles, "airport");
	//	gridder.addSectorGrid(comshapefiles, "com");
	//	gridder.addSectorGrid(indshapefiles, "ind");
	//	gridder.addSectorGrid(resshapefiles, "res");
	//	gridder.addSectorGrid(railroadshapefiles, "railroad");
	//	gridder.addSectorGrid(onroadshapefiles, "onroad");
	//	gridder.addSectorGrid(nonroadshapefiles, "nonroad");
	//	gridder.addSectorGrid(elecprodshapefiles, "elecprod");

	//	//gridder.addSectorGrid(marineshapefiles, "marine");
	//	gridder.loadAttribute(cafields[i]);

	//	//gridder.makeAnnualTotal(outroot + cityname + ".total.annual." + year + "." + version + ".nc");
	//	//gridder.toShapefile(onroadshapefiles[0], outroot + cityname + ".total.annual." + year + "." + version + ".Shp");

	//	gridder.makeAnnualTotal(outroot, cityname, year, version);
	//	gridder.toShapefile(onroadshapefiles[0], outroot + cityname + ".total.annual." + year + "." + version + ".shp");
	//	//gridder.makeHourlyTotal(outroot + cityname + ".total.hourly." + year + "." + version + ".nc");
	//	//for (size_t isector = 0; isector < gridder.sectorGrids.size(); isector++)
	//	//{
	//	//	HestiaGrid* outputGrid = gridder.sectorGrids[isector];
	//	//	outputGrid->getTotal();
	//	//	//outputGrid->toShapefile(indir + "onroad.Shp", lidarFisheyeDir + outputGrid->sectorname + "_2011.Shp");
	//	//	outputGrid->toShapefile(indir + "onroad.Shp", outroot + cityname + "." + outputGrid->sectorname + "." + year + "." + version + ".Shp");

	//	//}
	//	//printf("%f,%f\n", checkTotal(outroot + cityname + ".total.annual." + year + "." + version + ".nc"), checkTotal(outroot + cityname + ".total.hourly." + year + "." + version + ".nc"));

	//}
}
#include <qdatetime.h>
#include <osg/Vec3>
#include <osg/Image>
#include <osgDB/readFile>
#include <osgDB/WriteFile>
//template<typename T>
//void renderImage(T* outlineDS, std::vector<T> bins, size_t xsize, size_t ysize, T nodata, const char* output_filename, osg::Image* colorRamp)
//{
//	
//	size_t ncells = xsize*ysize;
//	osg::ref_ptr<osg::Image> image = new osg::Image;
//	double range = bins.size() - 1;
//	image->allocateImage(xsize, ysize, 1, GL_RGBA, GL_BYTE);
//
//	unsigned char* pImageData = image->outlineDS();
//	for (T* pdata = outlineDS; pdata < outlineDS + ncells; pdata++)
//	{
//		T val = *pdata;
//		if (val == nodata)
//		{
//			*pImageData++ = 0; *pImageData++ = 0; *pImageData++ = 0; *pImageData++ = 0;
//			continue;
//		}
//		double pos = 0;
//		for (size_t n = 0; n < bins.size(); n++)
//		{
//			if (val <= bins[n])
//			{
//				pos = n;
//				break;
//			}
//		}
//
//		double normalizedPos = pos / range;
//		if (normalizedPos > 1)
//			normalizedPos = 1;
//		if (normalizedPos < 0)
//			normalizedPos = 0;
//		osg::Vec4 color = colorRamp->getColor(osg::Vec2(normalizedPos, 0.5));
//		*pImageData++ = (unsigned char)(color.r() * 255); *pImageData++ = (unsigned char)(color.g() * 255); *pImageData++ = (unsigned char)(color.b() * 255); *pImageData++ = 255;
//	}
//
//	image->flipVertical();
//	osgDB::writeImageFile(*image, output_filename);
//}
//void createColorRamp(std::string outfile, int xsize,int ysize,osg::Vec3 startColor, osg::Vec3 endColor)
//{
//	osg::Image* img = new osg::Image;
//
//	img->allocateImage(xsize, ysize, 1, GL_RGB, GL_BYTE);
//	unsigned char* outlineDS = img->outlineDS();
//	for (int i = 0; i < xsize; i++)
//	{
//		float rate = (float)i / (xsize - 1);
//		osg::Vec3 color = startColor + (endColor - startColor) * rate;
//		for (int j = 0; j < ysize; j++)
//		{
//			
//		}
//	}
//
//	osgDB::writeImageFile(*img, outfile);
//
//}
void gridMarionForFFDASODIAC()
{


	std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[3]{ "ca11","length","area" }, 3);
	ShapeFile::copyDir("B:/Marion/2011/", "B:/Marion/ca11/", fields2keep);
	Preprocessor::gridFolderByRaster("B:/Marion/ca11/", "B:/Hestia_FFDAS_ODIAC/MarionFFDAS/", "B:/Hestia_FFDAS_ODIAC/Comparison/Marion_FFDAS.tif");
}

void gridForEntropy()
{
	std::string names[] = { "Los_Angeles","Marion","SaltLake","Baltimore" };
	std::string rootdir = "B:/SpatialGranuality/UrbanArea/";
	std::string outdir = "B:/SpatialGranuality/UrbanArea/500m/";
	std::string rasterdir = "B:/SpatialGranuality/UrbanArea/MODIS/";
	std::string resols[] = { "100m","200m","300m","400m","500m" };
	std::ofstream ofs;
	//ofs.open("e:/randomness_comnonpoint.csv");
	ofs.open("e:/randomness_resnonpoint.csv");

	ofs << "Resolution," << names[0] << "," << names[1] << "," << names[2] << "," << names[3] << std::endl;
	for (size_t iresol = 0; iresol < 5; iresol++)
	{
		std::string resol = resols[iresol];
		ofs << resol;
		for (size_t icounty = 0; icounty < 4; icounty++)
		{
			std::string county = names[icounty];
			std::string comnonpointfile = rootdir + county + "/" + "ResNonPoint.shp";
			std::string fishnetFile = outdir + county + "_fishnet" + resol + ".shp";
			std::string outputfile = outdir + county + "_ResNonPoint" + resol + ".shp";
			if (!QFileInfo(outputfile.data()).exists())
			{
				Grid grid;
				grid.fromFishnetRaster(rasterdir + county + resol + ".tif");
				const char *pszFormat = "GTiff";
				char **papszOptions = NULL;
				GDALDataset* pDataset = (GDALDataset*)GDALOpen((rasterdir + county + ".tif").data(), GA_ReadOnly);
				grid.toShape(pDataset->GetProjectionRef(), fishnetFile, false);
				GDALClose((GDALDatasetH)pDataset);
				Preprocessor::gridShapeFile(comnonpointfile, outputfile, fishnetFile, "ca11");
			}
			Grid grid1;
			grid1.fromFishnetRaster(rasterdir + county + resol + ".tif", true);

			Grid grid2;
			grid2.fromFishnetRaster(outdir + county + "_ResNonPoint" + resol + ".tif", true);

			double totalurban = 0;
			double totalcomnonpoint = 0;
			for (size_t i = 0; i < grid1.slice; i++)
			{
				if (grid1.cells[i] <= 0)
					continue;
				totalurban += 1;
				if (grid2.cells[i] > 0)
					totalcomnonpoint += 1;
			}
			ofs << "," << totalcomnonpoint / totalurban * 100;
			printf("%s=%f\n", (county + resol).data(), totalcomnonpoint / totalurban * 100);
		}
		ofs << std::endl;

	}
	ofs.close();
	//std::string names[] = { "Los_Angeles","Marion","SaltLake","Baltimore" };
	//std::string rootdir = "B:/SpatialGranuality/UrbanArea/100m/";
	//for (int i = 0; i < 4; i++)
	//{
	//	std::string comnonpointfile = rootdir + names[i] + "_" + "ComNonPoint" + ".tif";
	//	GDALDataset* pDataset = (GDALDataset*)GDALOpen(comnonpointfile.outlineDS(), GA_ReadOnly);
	//	int ncols = pDataset->GetRasterXSize();
	//	int nrows = pDataset->GetRasterYSize();
	//	int numcells = ncols * nrows;
	//	int* cells = new int[numcells];
	//	const char *pszFormat = "GTiff";
	//	char **papszOptions = NULL;
	//	double _adfGeoTransform[6];
	//	pDataset->GetGeoTransform(_adfGeoTransform);
	//	OGREnvelope bound;
	//	bound.MinX = _adfGeoTransform[0];
	//	bound.MaxY = _adfGeoTransform[3];
	//	bound.MaxX = _adfGeoTransform[0] + _adfGeoTransform[1] * ncols;
	//	bound.MinY = _adfGeoTransform[3] + _adfGeoTransform[5] * nrows;
	//	pDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, ncols, nrows, cells, ncols, nrows, GDT_Int32, 0, 0);
	//	GDALClose((GDALDatasetH)pDataset);

	//	for (int ncell = 0; ncell < numcells; ncell++){
	//		cells[ncell] = 1 - cells[ncell];
	//	}

	//	int numnulls = 0;

	//	for (int ncell = 0; ncell < numcells; ncell++) {
	//		if (cells[ncell] == 0)
	//			numnulls++;
	//	}
	//	printf("%s,%f\n", names[i].outlineDS(), (double)numnulls / (double)numcells * 100);
	//	delete[] cells;

	//}

}



#include "ScaleByFuel.h"
#include "RemoteSensing.h"

void sumSector()
{
	SubdirManager dirmanager = SubdirManager("B:/LA_Version2/gridoutput_nogeometry/Ventura/");
	//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";
	std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "RailroadPoint","Railroad" }, 2));
	std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.shp" }, 1));
	std::vector<std::string> allshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ ".shp" }, 1));



	double sum = ShapeFile::getTotal(allshapefiles, "ca11");
	printf("total=%f\n", sum);
	//std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	//std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	//std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));


}

//
//void updateLAResNonPointTime()
//{
//	//int srcYear = 2011;
//	//std::string id;
//	////double* srcFracs = TimestructTool::readText("C:/HestiaGridding/Los_Angeles/12_2011.csv",id,true);
//	////TimestructTool::toText("C:/HestiaGridding/Los_Angeles/12_2011.csv", srcFracs);
//	////delete[] srcFracs;
//	//for (int iyear = 2010; iyear <= 2014; iyear++)
//	//{
//	//	std::stringstream ssout;
//	//	ssout << "C:/HestiaGridding/Los_Angeles/" << "12_" << iyear << ".csv";
//
//
//	//	std::stringstream ssindir;
//	//	ssindir << "C:/HestiaGridding/Los_Angeles/ori/" << iyear << "/";
//
//	//	//if (iyear != srcYear)
//	//	//{
//	//	//	double* destFracs = TimestructTool::timeshift(srcFracs, srcYear, iyear);
//	//	//	TimestructTool::toText(ssout.str(), destFracs);
//	//	//}
//	//	std::stringstream ss;
//	//	ss << "C:/HestiaGridding/Los_Angeles/timestructs_" << iyear << ".bin";
//	//	TimestructTool::updateBinary(ss.str(), ss.str(), ssindir.str(),true);
//	//	TimestructTool::selectFromBinary2Text(ss.str(), "12", ssout.str());
//	//}
//
//	std::string id = "66111525";
//	//double* srcFracs = TimestructTool::readText("C:/HestiaGridding/Los_Angeles/" + id + "_2011.csv",id,true);
//	//double* destFracs = TimestructTool::timeshift(srcFracs, 2011, 2014);
//	//TimestructTool::toText("C:/HestiaGridding/Los_Angeles/" + id + "_2014.csv", destFracs);
//	//delete[] srcFracs;
//	for (int iyear = 2010; iyear <= 2014; iyear++)
//	{
//		std::stringstream ssout;
//		ssout << "C:/HestiaGridding/Los_Angeles/" << id << "_" << iyear << ".csv";
//		std::stringstream ss;
//		ss << "C:/HestiaGridding/Los_Angeles/timestructs_" << iyear << ".bin";
//		TimestructTool::selectFromBinary2Text(ss.str(), id, ssout.str());
//	}
//
//}

void processOnRoad()
{
	std::string emfacdir = "B:/LA_Version2/OnRoad/";
	std::vector<std::string> yearnames = Utils::buildVector("", new std::string[5]{ "10","11","12","13","14" }, 5);
	std::vector<std::string> emfacfiles;
	emfacfiles.push_back(emfacdir + "fuelsplit_2010.csv");
	emfacfiles.push_back(emfacdir + "fuelsplit_2011.csv");
	emfacfiles.push_back(emfacdir + "fuelsplit_2012.csv");
	emfacfiles.push_back(emfacdir + "fuelsplit_2013.csv");
	emfacfiles.push_back(emfacdir + "fuelsplit_2014.csv");


	EMFAC_FuelSplit emfacFUEL;
	std::string rootdir = "B:/LA_Version2/";
	std::vector<std::string> citynames;
	citynames.push_back("Los_Angeles"); citynames.push_back("Orange"); citynames.push_back("San_Bernardino"); citynames.push_back("Riverside"); citynames.push_back("Ventura");
	std::vector<std::string> citynames2;
	citynames2.push_back("Los Angeles"); citynames2.push_back("Orange"); citynames2.push_back("San Bernardino"); citynames2.push_back("Riverside"); citynames2.push_back("Ventura");
	std::vector<std::string> cafields2keep = Utils::buildVector("", new std::string[8]{ "ca10","ca11","ca12","ca13" ,"ca14","length","area" ,"timestruct" }, 8);
	std::vector<std::string> ca14fields2keep = Utils::buildVector("", new std::string[4]{ "ca14","length","area" ,"timestruct" }, 4);

	std::vector<std::string> publicfields2keep = Utils::buildVector("", new std::string[16]{
		"ca10","ca11","ca12","ca13","ca14",
		"ca10_g","ca11_g","ca12_g","ca13_g","ca14_g",
		"ca10_d","ca11_d","ca12_d","ca13_d","ca14_d","timestruct" }, 16);

	//for (size_t i = 0; i < 5; i++)
	//{
	//	std::string cityname = citynames[i];
	//	std::string cityname2 = citynames2[i];
	//	std::string masterfile = rootdir + "gridPrep_SHP_master/" + cityname + "/OnRoad.Shp";
	//	std::string publicfile = rootdir + "public/" + cityname + "/OnRoad.Shp";
	//	emfacFUEL.allocate(masterfile, emfacfiles, yearnames, cityname2);
	//	std::string cafile = rootdir + "ca/" + cityname + "/OnRoad.Shp";
	//	std::string ca14file = rootdir + "ca14/" + cityname + "/OnRoad.Shp";
	//	ShapeFile::copy(masterfile, publicfile, publicfields2keep);
	//	ShapeFile::copy(masterfile, cafile, cafields2keep);
	//	ShapeFile::copy(masterfile, ca14file, ca14fields2keep);
	//
	//	ShapeFile::copy(masterfile, "Z:/Hestia/LAbasin/LAbasin_v2.0/gridPrep_SHP_master/OnRoad.Shp");
	//	ShapeFile::copy(publicfile, "Z:/Hestia/LAbasin/LAbasin_v2.0/public_shapefiles/OnRoad.Shp");
	//}

	for (size_t i = 0; i < 5; i++)
	{
		std::string cityname = citynames[i];
		std::string cityname2 = citynames2[i];
		std::string masterfile = rootdir + "gridPrep_SHP_master/" + cityname + "/OnRoad.shp";
		std::string publicfile = rootdir + "public/" + cityname + "/OnRoad.shp";

		ShapeFile::copy(masterfile, "Z:/Hestia/LAbasin/LAbasin_v2.0/gridPrep_SHP_master/" + cityname + "/OnRoad.shp");
		ShapeFile::copy(publicfile, "Z:/Hestia/LAbasin/LAbasin_v2.0/public_shapefiles/" + cityname + "/OnRoad.shp");
	}
}
struct BlockGroup
{
	//std::string StateName;
	std::string StateCode;
	//std::string CountyName;
	std::string CountyCode;
	std::string CensusTractCode;
	std::string BlockGroupCode;
	//double Population;
	//double Pcapincome;
	//double Aggincome;
	std::string BlockGroupID;
	std::string CensusTractID;
	std::vector<double> attributes;
	BlockGroup()
	{

	}
	BlockGroup(std::string line)
	{
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		int iarg = 0;
		//StateName = splits[0];
		StateCode = splits[iarg++];
		//CountyName = splits[iarg++];
		CountyCode = splits[iarg++];
		CensusTractCode = splits[iarg++];
		BlockGroupCode = splits[iarg++];
		for (size_t i = iarg; i < splits.size(); i++)
		{
			attributes.push_back(atof(splits[i].data()));
		}
		//Population = atoi(splits[7].data());
		//Pcapincome = atoi(splits[iarg++].data());
		//Aggincome = atoi(splits[iarg++].data());

		while (StateCode.size() < 2)
		{
			StateCode = "0" + StateCode;
		}

		while (CountyCode.size() < 3)
		{
			CountyCode = "0" + CountyCode;
		}

		while (CensusTractCode.size() < 6)
		{
			CensusTractCode = "0" + CensusTractCode;
		}

		BlockGroupID = StateCode + CountyCode + CensusTractCode + BlockGroupCode;
		CensusTractID = StateCode + CountyCode + CensusTractCode;
	}
};
void jointBlockGroup(std::string csvfile, std::string shapefile)
{


	std::ifstream ifs;
	ifs.open(csvfile);
	std::string line;

	std::getline(ifs, line);
	std::vector<std::string> splits = Utils::splitCSV(',', line);

	std::map<std::string, BlockGroup> bkmap;
	ShapeFile shp(shapefile, 1);

	int stateIdx = shp.poLayer->GetLayerDefn()->GetFieldIndex("STATE");
	int blockgroupIdx = shp.poLayer->GetLayerDefn()->GetFieldIndex("BLKGRP");
	int tractIdx = shp.poLayer->GetLayerDefn()->GetFieldIndex("TRACT");
	int countyIdx = shp.poLayer->GetLayerDefn()->GetFieldIndex("COUNTY");
	std::vector<int> attributes;
	for (size_t i = 4; i < splits.size(); i++)
	{
		std::string fname = splits[i];
		if (fname.size() > 10)
			fname = fname.substr(0, 10);
		attributes.push_back(shp.getOrCreateField(fname.data(), OGRFieldType::OFTReal));
	}
	//StateName	StateCode	CountyName	CountyCode	CensusTractCode	BlockGroupCode	Population	Pcapincome	Aggincome
	//	Alabama	1	Autauga County	1	20100	1	499	27411	13678000

	while (ifs.peek() != -1)
	{

		std::getline(ifs, line);
		BlockGroup bk(line);
		//if (bk.StateCode != "6")
		//	continue;
		bkmap[bk.BlockGroupID] = bk;
	}
	ifs.close();



	//int popIdx = Shp.getOrCreateField("Population", OGRFieldType::OFTReal);
	//int pcapincomeIdx = Shp.getOrCreateField("Pcapincome", OGRFieldType::OFTReal);
	//int aggincomeIdx = Shp.getOrCreateField("Aggincome", OGRFieldType::OFTReal);

	OGRFeature *poFeature;
	double sum = 0;
	shp.poLayer->ResetReading();
	int idx = -1;

	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		idx++;

		//int stateIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("STATE");
		//int blockgroupIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("BLKGRP");
		//int tractIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("TRACT");
		//int countyIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("COUNTY");
		std::string statestr = poFeature->GetFieldAsString(stateIdx);
		std::string countystr = poFeature->GetFieldAsString(countyIdx);
		std::string tractstr = poFeature->GetFieldAsString(tractIdx);
		std::string blockgroupstr = poFeature->GetFieldAsString(blockgroupIdx);

		std::string id = statestr + countystr + tractstr + blockgroupstr;
		std::map<std::string, BlockGroup>::iterator iter = bkmap.find(id);
		if (iter == bkmap.end())
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		BlockGroup pbk = iter->second;
		for (size_t i = 0; i < attributes.size(); i++)
		{
			poFeature->SetField(attributes[i], pbk.attributes[i]);
		}
		//poFeature->SetField(popIdx, pbk.Population);
		//poFeature->SetField(pcapincomeIdx, pbk.Pcapincome);
		//poFeature->SetField(aggincomeIdx, pbk.Aggincome);
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();


}
//void jointCensusTract(std::string csvfile, std::string yuyuShapefile)
//{
//
//
//	std::ifstream ifs;
//	ifs.open(csvfile);
//	std::string line;
//
//	std::getline(ifs, line);
//	std::vector<std::string> splits = Utils::splitCSV(',', line);
//	std::vector<std::string> fields;
//	std::map<std::string, BlockGroup> bkmap;
//
//	for (size_t i = 1; i < splits.size(); i++)
//	{
//		fields.push_back(splits[i]);
//	}
//
//	//StateName	StateCode	CountyName	CountyCode	CensusTractCode	BlockGroupCode	Population	Pcapincome	Aggincome
//	//	Alabama	1	Autauga County	1	20100	1	499	27411	13678000
//
//	while (ifs.peek() != -1)
//	{
//
//		std::getline(ifs, line);
//		BlockGroup bk(line);
//		if (bk.StateCode != "6")
//			continue;
//		std::map<std::string, BlockGroup>::iterator iter = bkmap.find(bk.CensusTractID);
//		if (iter == bkmap.end())
//		{
//			bkmap[bk.CensusTractID] = bk;
//		}
//		else
//		{
//			iter->second.Aggincome += bk.Aggincome;
//			iter->second.Population += bk.Population;
//			iter->second.Pcapincome = iter->second.Aggincome/ iter->second.Population;
//		}
//	}
//	ifs.close();
//
//	ShapeFile Shp(yuyuShapefile, 1);
//	int stateIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("STATE");
//	int blockgroupIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("BLKGRP");
//	int tractIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("TRACT");
//	int countyIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("COUNTY");
//
//	int popIdx = Shp.getOrCreateField("Population", OGRFieldType::OFTReal);
//	int pcapincomeIdx = Shp.getOrCreateField("Pcapincome", OGRFieldType::OFTReal);
//	int aggincomeIdx = Shp.getOrCreateField("Aggincome", OGRFieldType::OFTReal);
//
//	OGRFeature *poFeature;
//	double sum = 0;
//	Shp.poLayer->ResetReading();
//	int idx = -1;
//
//	while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
//	{
//		idx++;
//
//		//int stateIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("STATE");
//		//int blockgroupIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("BLKGRP");
//		//int tractIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("TRACT");
//		//int countyIdx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("COUNTY");
//		std::string statestr = poFeature->GetFieldAsString(stateIdx);
//		std::string countystr = poFeature->GetFieldAsString(countyIdx);
//		std::string tractstr = poFeature->GetFieldAsString(tractIdx);
//		std::string blockgroupstr = poFeature->GetFieldAsString(blockgroupIdx);
//
//		std::string id = countystr + tractstr;
//		std::map<std::string, BlockGroup>::iterator iter = bkmap.find(id);
//		if (iter == bkmap.end())
//		{
//			OGRFeature::DestroyFeature(poFeature);
//			continue;
//		}
//		BlockGroup pbk = iter->second;
//		poFeature->SetField(popIdx, pbk.Population);
//		poFeature->SetField(pcapincomeIdx, pbk.Pcapincome);
//		poFeature->SetField(aggincomeIdx, pbk.Aggincome);
//		Shp.poLayer->SetFeature(poFeature);
//		OGRFeature::DestroyFeature(poFeature);
//	}
//	Shp.close();
//
//
//}


void createPolygonGeometry()
{
	OGREnvelope bound;
	bound.MinX = -50;
	bound.MaxX = 50;
	bound.MinY = -50;
	bound.MaxY = 50;
	double xspan = bound.MaxX - bound.MinX;
	double yspan = bound.MaxY - bound.MinY;
	ShapeFile shp;
	shp.create("B:/SpatialGranuality/Synthesized/DenseSmallPolygon.shp", 0, 0, OGRwkbGeometryType::wkbPolygon);
	shp.getOrCreateField("ca11", OGRFieldType::OFTReal);
	int numrows = 10;
	int numcols = 10;
	double cellsizex = xspan / numcols;
	double cellsizey = yspan / numrows;
	srand(time(NULL));
	for (size_t irow = 0; irow < numrows; irow++)
	{
		for (size_t icol = 0; icol < numcols; icol++)
		{
			OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
			OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);

			double cx = bound.MinX + cellsizex * (icol + 0.5);
			double cy = bound.MaxY - cellsizey * (irow + 0.5);
			double trisize = xspan / 4 * 0.5 * 0.3;

			linearRing->addPoint(cx, cy + trisize);
			linearRing->addPoint(cx + trisize, cy - trisize);
			linearRing->addPoint(cx - trisize, cy - trisize);
			linearRing->addPoint(cx, cy + trisize);
			poPolygon->addRing(linearRing);//also crashed
			OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
			poFeaPolygon->SetField("ca11", 100000 + rand() % 100000);
			poFeaPolygon->SetGeometry(poPolygon);
			shp.poLayer->CreateFeature(poFeaPolygon);
			OGRFeature::DestroyFeature(poFeaPolygon);
		}
	}

	shp.close();

}
void mergeCounties()
{
	std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[6]{ "ca11","ca11_p","ca11_ng","ca11_c","tfs","length" }, 6);
	SubdirManager dirmanager = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/");
	//std::string indir = "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/";
	std::vector<std::string> comnonpointshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "comnonpoint" }, 1));
	std::vector<std::string> compointshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "compoint" }, 1));
	std::vector<std::string> indnonpointshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "indnonpoint" }, 1));
	std::vector<std::string> indpointshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "indpoint" }, 1));
	std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "Railroad.shp" }, 1));
	std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.shp" }, 1));
	std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "Airport" }, 1));

	//std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	//std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	//std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));
	std::string outdir = "B:/LA_Version2/Mapping/";
	ShapeFile::mergeFile(comnonpointshapefiles, fields2keep, outdir + "ComNonPoint.shp");
	//ShapeFile::mergeFile(compointshapefiles, fields2keep, lidarFisheyeDir + "ComPoint.Shp");
	//ShapeFile::mergeFile(indnonpointshapefiles, fields2keep, lidarFisheyeDir + "IndNonPoint.Shp");
	//ShapeFile::mergeFile(indpointshapefiles, fields2keep, lidarFisheyeDir + "IndPoint.Shp");
	//ShapeFile::mergeFile(resshapefiles, fields2keep, lidarFisheyeDir + "ResNonPoint.Shp");
	//ShapeFile::mergeFile(onroadshapefiles, fields2keep, lidarFisheyeDir + "OnRoad.Shp");
	//ShapeFile::mergeFile(railroadshapefiles, fields2keep, lidarFisheyeDir + "Railroad.Shp");
	//ShapeFile::mergeFile(elecprodshapefiles, fields2keep, lidarFisheyeDir + "ElecProd.Shp");
	//ShapeFile::mergeFile(airportshapefiles, fields2keep, lidarFisheyeDir + "Airport.Shp");


}

void copy(std::string src, std::string dest)
{
	ShapeFile srcshp(src, 0);
	ShapeFile destshp;
	destshp.create(dest, srcshp.poLayer->GetSpatialRef(), 0, srcshp.poLayer->GetGeomType());

	std::vector<int> oldfields;
	std::vector<int> newfields;
	OGRFeatureDefn* layerdef = srcshp.poLayer->GetLayerDefn();
	int fieldcount = layerdef->GetFieldCount();
	for (size_t i = 0; i < fieldcount; i++)
	{
		int oldidx = i;
		OGRFieldDefn* fdefn = layerdef->GetFieldDefn(oldidx);
		destshp.poLayer->CreateField(fdefn);
		newfields.push_back(destshp.poLayer->GetLayerDefn()->GetFieldIndex(fdefn->GetNameRef()));
		oldfields.push_back(oldidx);
	}
	srcshp.poLayer->ResetReading();
	OGRFeature *poFeatureOld;
	int fipsid = srcshp.poLayer->GetLayerDefn()->GetFieldIndex("FIPS");
	while ((poFeatureOld = srcshp.poLayer->GetNextFeature()) != NULL)
	{

		std::string fips = poFeatureOld->GetFieldAsString(fipsid);
		if (fips.substr(0, 5) != "06037")
		{
			OGRFeature::DestroyFeature(poFeatureOld);
			continue;
		}
		OGRFeature *poFeatureNew = OGRFeature::CreateFeature(destshp.poLayer->GetLayerDefn());
		poFeatureNew->SetGeometry(poFeatureOld->GetGeometryRef());
		for (size_t i = 0; i < fieldcount; i++)
		{
			int oldfieldidx = oldfields[i];
			int newfieldidx = newfields[i];
			if (oldfieldidx == -1 || newfieldidx == -1)
				continue;
			poFeatureNew->SetField(newfieldidx, poFeatureOld->GetRawFieldRef(oldfieldidx));
		}

		destshp.poLayer->CreateFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureOld);
	}

}
void updateFEMA_SQFT(std::string filename)
{
	ShapeFile shp(filename, 1);

	int comsqft = shp.getOrCreateField("COMSQFT", OGRFieldType::OFTReal);
	int ressqft = shp.getOrCreateField("RESSQFT", OGRFieldType::OFTReal);
	int indsqft = shp.getOrCreateField("INDSQFT", OGRFieldType::OFTReal);

	std::vector<int> resfields;
	std::vector<int> comfields;
	std::vector<int> indfields;
	OGRFeatureDefn* layerdef = shp.poLayer->GetLayerDefn();
	int fieldcount = layerdef->GetFieldCount();

	std::vector<std::string> resfieldnames = Utils::buildVector("", new std::string[6]{ "RES1","RES2","RES3","RES4","RES5","RES6" }, 6);
	std::vector<std::string> indfieldnames = Utils::buildVector("", new std::string[6]{ "IND1","IND2","IND3","IND4","IND5","IND6" }, 6);
	std::vector<std::string> comfieldnames = Utils::buildVector("", new std::string[16]{ "COM1","COM2","COM3","COM4","COM5","COM6","COM7","COM8","COM9","COM10","AGR1","REL1","GOV1","GOV2","EDU1","EDU2" }, 16);
	for (size_t i = 0; i < resfieldnames.size(); i++)
	{
		resfields.push_back(layerdef->GetFieldIndex(resfieldnames[i].data()));
	}
	for (size_t i = 0; i < indfieldnames.size(); i++)
	{
		indfields.push_back(layerdef->GetFieldIndex(indfieldnames[i].data()));
	}
	for (size_t i = 0; i < comfieldnames.size(); i++)
	{
		comfields.push_back(layerdef->GetFieldIndex(comfieldnames[i].data()));
	}

	shp.poLayer->ResetReading();
	OGRFeature *poFeature;
	int fipsid = shp.poLayer->GetLayerDefn()->GetFieldIndex("FIPS");

	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		double RES_SQRT = 0;
		double IND_SQRT = 0;
		double COM_SQRT = 0;
		for (size_t i = 0; i < resfields.size(); i++)
		{
			RES_SQRT += poFeature->GetFieldAsDouble(resfields[i]);
		}
		for (size_t i = 0; i < comfields.size(); i++)
		{
			COM_SQRT += poFeature->GetFieldAsDouble(comfields[i]);
		}
		for (size_t i = 0; i < indfields.size(); i++)
		{
			IND_SQRT += poFeature->GetFieldAsDouble(indfields[i]);
		}
		poFeature->SetField(comsqft, COM_SQRT);
		poFeature->SetField(ressqft, RES_SQRT);
		poFeature->SetField(indsqft, IND_SQRT);


		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}
struct BLOCKGROUP_SQFT
{
	double COMSQFT;
	double RESSQFT;
	double INDSQFT;
	double COMSQFT_FEMA;
	double RESSQFT_FEMA;
	double INDSQFT_FEMA;
	int id;
	BLOCKGROUP_SQFT()
	{
		COMSQFT = 0;
		RESSQFT = 0;
		INDSQFT = 0;
		COMSQFT_FEMA = 0;
		RESSQFT_FEMA = 0;
		INDSQFT_FEMA = 0;
	}
};
void aggregrateToFEMA(std::string filename)
{
	ShapeFile shp(filename);

	int comsqft = shp.getField("COMSQFT");
	int ressqft = shp.getField("RESSQFT");
	int indsqft = shp.getField("INDSQFT");
	int bc = shp.getField("bc");
	int tfs = shp.getField("tfs");
	int blockgroupid = shp.getField("toID");

	OGRFeatureDefn* layerdef = shp.poLayer->GetLayerDefn();

	shp.poLayer->ResetReading();
	OGRFeature *poFeature;

	std::map<int, BLOCKGROUP_SQFT> tsf_table;
	std::map<int, BLOCKGROUP_SQFT>::iterator iter;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		int id = poFeature->GetFieldAsInteger(blockgroupid);
		if (id < 0)
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		iter = tsf_table.find(id);
		int sector = poFeature->GetFieldAsInteger(bc);
		double parcelTFS = poFeature->GetFieldAsDouble(tfs);
		if (iter == tsf_table.end())
		{
			BLOCKGROUP_SQFT newBG;
			if (sector == 1) {
				newBG.COMSQFT = parcelTFS;
			}
			else if (sector == 2) {
				newBG.RESSQFT = parcelTFS;
			}
			else {
				newBG.INDSQFT = parcelTFS;
			}
			newBG.COMSQFT_FEMA = poFeature->GetFieldAsDouble(comsqft) * 1000;
			newBG.RESSQFT_FEMA = poFeature->GetFieldAsDouble(ressqft) * 1000;
			newBG.INDSQFT_FEMA = poFeature->GetFieldAsDouble(indsqft) * 1000;
			newBG.id = id;
			tsf_table[id] = newBG;
		}
		else
		{
			if (sector == 1) {
				iter->second.COMSQFT += parcelTFS;
			}
			else if (sector == 2) {
				iter->second.RESSQFT += parcelTFS;
			}
			else {
				iter->second.INDSQFT += parcelTFS;
			}
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();
	std::ofstream ofsCom("B:/LA_Version2/Fema_LA/parcel_vs_blockgroup.csv");
	iter = tsf_table.begin();
	ofsCom << "ID," << "Parcels_COMSQFT," << "FEMA_COMSQFT," << "Parcels_RESSQFT," << "FEMA_RESSQFT," << "Parcels_INDSQFT," << "FEMA_INDSQFT" << std::endl;
	while (iter != tsf_table.end())
	{
		BLOCKGROUP_SQFT& bks = iter->second;
		ofsCom << bks.id << ","
			<< bks.COMSQFT << "," << bks.COMSQFT_FEMA << ","
			<< bks.RESSQFT << "," << bks.RESSQFT_FEMA << ","
			<< bks.INDSQFT << "," << bks.INDSQFT_FEMA << std::endl;
		iter++;
	}

	ofsCom.close();
}
void aggregrate(std::string dir, std::string filename)
{
	std::map<int, double> tb;
	std::vector<std::string> files;
	QDir input_dir(dir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	OGRFeature *poFeature;
	QFileInfoList list = input_dir.entryInfoList();
	int ca11field;
	int id = 0;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		ShapeFile shp(input_file);
		int blockgroupid = shp.getField("FID_blockg");

		ca11field = shp.getField("ca11");
		while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
		{
			id = poFeature->GetFieldAsInteger(blockgroupid);
			double val = poFeature->GetFieldAsDouble(ca11field);
			std::map<int, double>::iterator iter = tb.find(id);
			if (iter == tb.end())
			{
				tb[id] = val;
			}
			else
			{
				iter->second += val;
			}

			OGRFeature::DestroyFeature(poFeature);
		}
	}

	ShapeFile shpblk(filename, 1);

	ca11field = shpblk.getOrCreateField("ca11", OGRFieldType::OFTReal);

	shpblk.poLayer->ResetReading();
	id = 0;
	std::ofstream ofsCom((filename + ".csv").data());
	ofsCom << "id,ca11" << "\n";
	while ((poFeature = shpblk.poLayer->GetNextFeature()) != NULL)
	{
		std::map<int, double>::iterator iter = tb.find(id);
		if (iter != tb.end())
		{
			poFeature->SetField(ca11field, iter->second);
			ofsCom << id << "," << iter->second << "\n";
		}
		else
		{
			ofsCom << id << "," << 0 << "\n";
		}
		shpblk.poLayer->SetFeature(poFeature);
		id++;
		OGRFeature::DestroyFeature(poFeature);
	}
	shpblk.close();
	ofsCom.close();
}

//void sumNLCD(std::string infilename,std::string femaSQFT)
//{
//	ShapeFile Shp(infilename, 1);
//
//	int nlcd_field = Shp.getField("nlcd_score");
//	int nlcdbk_field = Shp.getOrCreateField("nlcd_bk", OGRFieldType::OFTReal);
//	int bk_field = Shp.getField("toID");
//	int sqft_field = Shp.getField("tfs");
//	int fema_sqft_field = Shp.getField(femaSQFT.outlineDS());
//	std::map<int, int> NLCD_BY_BLOCKGROUP;
//	OGRFeature *poFeature;
//	Shp.poLayer->ResetReading();
//	while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
//	{
//	
//		int bkID = poFeature->GetFieldAsInteger(bk_field);
//		int nlcd = poFeature->GetFieldAsInteger(nlcd_field);
//		poFeature->GetFieldAsInteger(bkID);
//		std::map<int, int>::iterator iter = NLCD_BY_BLOCKGROUP.find(bkID);
//		if (iter == NLCD_BY_BLOCKGROUP.end())
//		{
//			NLCD_BY_BLOCKGROUP[bkID] = 0;
//			iter = NLCD_BY_BLOCKGROUP.find(bkID);
//		}
//		iter->second += nlcd;
//		OGRFeature::DestroyFeature(poFeature);
//	}
//
//	Shp.poLayer->ResetReading();
//	while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
//	{
//		int bkID = poFeature->GetFieldAsInteger(bk_field);
//		double sqft_bk = poFeature->GetFieldAsDouble(fema_sqft_field);
//		double nlcd = (double)poFeature->GetFieldAsInteger(nlcd_field);
//		double nlcd_total = 0;
//		std::map<int, int>::iterator iter = NLCD_BY_BLOCKGROUP.find(bkID);
//		if (iter != NLCD_BY_BLOCKGROUP.end())
//		{
//			nlcd_total = (double)iter->second;
//		}
//		else
//		{
//			printf("");
//		}
//		double sqft = sqft_bk * 1000 * (nlcd / nlcd_total);
//		poFeature->SetField(sqft_field, sqft);
//		poFeature->SetField(nlcdbk_field, (int)nlcd_total);
//		Shp.poLayer->SetFeature(poFeature);
//		OGRFeature::DestroyFeature(poFeature);
//	}
//	Shp.close();
//
//}
#include "ExtractByPolygon.h"
struct Parcel
{
	int parcelID;
	std::vector<float> heights;
	int bt;
	int bc;
	std::vector<float> areas;
	float totalBuildingArea;
	float parcelFloorSpace;
	float maxHeight;
	int numbuildings;
	Parcel()
	{
		parcelID = -1;
		totalBuildingArea = 0;
		parcelFloorSpace = 0;
		bt = -1;
		maxHeight = 0;
		numbuildings = 0;
	}
	//This was approached as follows : 
	//if the residential building was less than 30 feet tall, it was considered a one - story building.
	//	If great than 30 feet, the number of floors was the rounded ratio of the height / 15.
	//Commercial buildings used the same approach but the number of floors used height / 20.
	//Industrial buildings were all considered one story tall.
	void addBuilding(int _parcelID, int _bt, float tfs, float height, float area)
	{
		parcelID = _parcelID; bt = _bt; parcelFloorSpace = tfs;
		//if (height > 0)
		//{ 
		int numstory = 1;// (int)(height / 15);
		totalBuildingArea += (area * numstory);
		//}
		heights.push_back(height);
		areas.push_back(area);
		if (height > maxHeight)
			maxHeight = height;
		numbuildings++;
	}
};

void aggregrateToParcels(std::string shapefile, std::string outcsv)
{
	ShapeFile shp(shapefile);
	int parcelIDIdx = shp.getField("FID_Parcel");
	int heightIdx = shp.getField("HEIGHT");
	int areaIdx = shp.getField("Shape_Area");
	int btIdx = shp.getField("CRMcode");
	int tfsIdx = shp.getField("tfs");
	Parcel* parcels = new Parcel[5000000];

	OGRFeatureDefn* layerdef = shp.poLayer->GetLayerDefn();

	shp.poLayer->ResetReading();
	OGRFeature *poFeature;

	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		int id = poFeature->GetFieldAsInteger(parcelIDIdx);
		float height = poFeature->GetFieldAsDouble(heightIdx);
		float area = poFeature->GetFieldAsDouble(areaIdx);
		int bt = poFeature->GetFieldAsInteger(btIdx);
		float tfs = poFeature->GetFieldAsDouble(tfsIdx);
		parcels[id].addBuilding(id, bt, tfs, height, area);
		OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();
	std::ofstream ofs(outcsv.data());

	ofs << "ParcelID," << "BT," << "BuildingArea," << "ParcelFloorSpace" << std::endl;
	for (size_t i = 0; i < 5000000; i++)
	{
		if (parcels[i].parcelID < 0)
			continue;
		Parcel& pc = parcels[i];
		ofs << pc.parcelID << "," << pc.bt << "," << pc.totalBuildingArea << "," << pc.parcelFloorSpace << std::endl;
	}

	ofs.close();
}

void aggregrateToParcelsUpdate(std::string shapefile, std::string outcsv, std::string parcelFile)
{

	ShapeFile shp(shapefile);


	int parcelIDIdx = shp.getField("FID_Parcel");
	int heightIdx = shp.getField("HEIGHT");
	int areaIdx = shp.getField("Shape_Area");
	int btIdx = shp.getField("bt");
	int tfsIdx = shp.getField("tfs");
	Parcel* parcels = new Parcel[5000000];

	OGRFeatureDefn* layerdef = shp.poLayer->GetLayerDefn();

	shp.poLayer->ResetReading();
	OGRFeature *poFeature;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		int id = poFeature->GetFieldAsInteger(parcelIDIdx);
		float height = poFeature->GetFieldAsDouble(heightIdx);
		float area = poFeature->GetFieldAsDouble(areaIdx);
		//int bt = poFeature->GetFieldAsInteger(btIdx);
		float tfs = poFeature->GetFieldAsDouble(tfsIdx);
		parcels[id].addBuilding(id, -1, tfs, height, area);
		if (id == 0)
		{
			printf("%d,%f\n", id, area);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();



	size_t fid = 0;

	ShapeFile shpParcel(parcelFile, 1);
	//int heightIdxParcel = shpParcel.getOrCreateField("HEIGHT", OGRFieldType::OFTReal);
	int areaIdxParcel = shpParcel.getOrCreateField("Shape_Area", OGRFieldType::OFTReal);
	int tfsIdxParcel = shpParcel.getOrCreateField("tfs", OGRFieldType::OFTReal);
	int oritfsIdxParcel = shpParcel.getOrCreateField("ori_tfs", OGRFieldType::OFTReal);
	int numbuildingsIdxParcel = shpParcel.getOrCreateField("numbds", OGRFieldType::OFTInteger);
	int maxheightIdxParcel = shpParcel.getOrCreateField("maxheight", OGRFieldType::OFTReal);
	int bcIdxParcel = shpParcel.getOrCreateField("bc", OGRFieldType::OFTInteger);

	int numfeas = shpParcel.poLayer->GetFeatureCount();
	printf("%d\n", numfeas);
	for (size_t i = 0; i < numfeas; i++)
	{
		poFeature = shpParcel.poLayer->GetFeature(i);
		double area = 0;
		double maxheight = 0;
		int numbuildings = 0;
		fid = i;
		double tfs = poFeature->GetFieldAsInteger(oritfsIdxParcel);
		poFeature->SetField(oritfsIdxParcel, tfs);
		double totaltfs = 0;
		if (parcels[fid].parcelID >= 0)
		{
			Parcel& pc = parcels[fid];
			//area = pc.totalBuildingArea;
			maxheight = pc.maxHeight;
			numbuildings = pc.numbuildings;
			//float tfs = poFeature->GetFieldAsDouble(tfsIdxParcel);
			int bc = poFeature->GetFieldAsInteger(bcIdxParcel);

			//This was approached as follows : 
			//if the residential building was less than 30 feet tall, it was considered a one - story building.
			//	If great than 30 feet, the number of floors was the rounded ratio of the height / 15.
			//Commercial buildings used the same approach but the number of floors used height / 20.
			//Industrial buildings were all considered one story tall.
			//if (pc.totalBuildingArea > 0 && tfs <= 0)
			//{
			for (size_t ibuilding = 0; ibuilding < pc.heights.size(); ibuilding++)
			{
				int numstories = 0;
				double height = pc.heights[ibuilding];
				double footprint = pc.areas[ibuilding];
				if (height > 30)
				{
					if (bc == 1)
					{
						numstories = (int)(height / 20);
					}
					else if (bc == 2)
					{
						numstories = (int)(height / 15);
					}
				}
				else if (height > 0)
				{
					numstories = 1;
				}
				totaltfs += numstories *footprint;
			}
			if (fid == 597103)
				totaltfs = 0;

		}

		if (tfs <= 0)
			tfs = totaltfs;
		tfs = (int)tfs;
		poFeature->SetField(tfsIdxParcel, tfs);
		poFeature->SetField(areaIdxParcel, totaltfs);

		//if (tfs == 0)
		//{
		//	printf("%d,%f\n", fid, area);
		//}
		//if (fid % 1000 == 0)
		//{
		//	printf("%d,%d\n", fid, numfeas);
		//}

		poFeature->SetField(maxheightIdxParcel, maxheight);
		poFeature->SetField(numbuildingsIdxParcel, numbuildings);
		shpParcel.poLayer->SetFeature(poFeature);
		fid++;
		OGRFeature::DestroyFeature(poFeature);
	}
	//while ((poFeature = shpParcel.poLayer->GetNextFeature()) != NULL)
	//{
	//	double area = 0;
	//	if (parcels[fid].parcelID > 0)
	//	{
	//		Parcel& pc = parcels[fid];
	//		area = pc.totalBuildingArea;
	//	}
	//	if (fid == 2065108)
	//	{
	//		printf("%d,%f\n", parcels[fid].parcelID, parcels[fid].totalBuildingArea);
	//	}
	//	if (fid % 1000==0)
	//	{
	//		printf("%d,%d\n", fid, numfeas);
	//	}
	//	//poFeature->SetField(areaIdxParcel, area);
	//	//shpParcel.poLayer->SetFeature(poFeature);
	//	fid++;
	//	OGRFeature::DestroyFeature(poFeature);
	//}
	//printf("%d\n", fid);
	shpParcel.close();


	std::ofstream ofs(outcsv.data());

	ofs << "ParcelID," << "BT," << "BuildingArea," << "ParcelFloorSpace" << std::endl;
	for (size_t i = 0; i < 5000000; i++)
	{
		if (parcels[i].parcelID < 0)
			continue;
		Parcel& pc = parcels[i];
		ofs << pc.parcelID << "," /*<< pc.bt << ","*/ << pc.totalBuildingArea << "," << pc.parcelFloorSpace << std::endl;
	}
	delete[] parcels;
	ofs.close();
}

//void update(std::string yuyuShapefile, std::string outcsv, std::string parcelFile)
//{
//
//	size_t fid = 0;
//	ShapeFile shpParcel(parcelFile, 1);
//	//int heightIdxParcel = shpParcel.getOrCreateField("HEIGHT", OGRFieldType::OFTReal);
//	int areaIdxParcel = shpParcel.getOrCreateField("Shape_Area", OGRFieldType::OFTReal);
//	int numfeas = shpParcel.poLayer->GetFeatureCount();
//	printf("%d\n", numfeas);
//	for (size_t i = 0; i < numfeas; i++)
//	{
//		OGRFeature* poFeature = shpParcel.poLayer->GetFeature(i);
//		
//		poFeature->SetField(areaIdxParcel, area);
//		shpParcel.poLayer->SetFeature(poFeature);
//		fid++;
//		OGRFeature::DestroyFeature(poFeature);
//	}
//
//	shpParcel.close();
//
//}
//void computeSVF_In_Manhattan()
//{
//	std::string dataDir = "E:/GoogleStreetview/Manhattan/Data/";
//	std::string fisheyeDir = "E:/GoogleStreetview/Manhattan/SegnetFisheye2/";
//	std::string outshapefile = "E:/GoogleStreetview/Manhattan/svf.Shp";
//	ShapeFile Shp(outshapefile,1);
//	//Shp.create(outshapefile, 0, 0, OGRwkbGeometryType::wkbPoint);
//	int idIdx = Shp.getOrCreateField("id", OGRFieldType::OFTInteger);
//	int panoidIdx = Shp.getOrCreateField("panoid", OGRFieldType::OFTString);
//	int svfIdx = Shp.getOrCreateField("svf", OGRFieldType::OFTReal);
//
//	SVFComputer<unsigned char>* googlesvfcom = new SVFComputer<unsigned char>();
//	for (int i = 0; i < 12528; i++)
//	{
//		std::stringstream fisheyefiless;
//		fisheyefiless << fisheyeDir << i << ".png";
//		std::string fisheyefile = fisheyefiless.str();
//		if (!QFileInfo(fisheyefile.outlineDS()).exists())
//			continue;
//		//GDAL_DS<unsigned char>* dsImage = new GDAL_DS<unsigned char>();
//		//dsImage->open(fisheyefile);
//		//googlesvfcom->setData(dsImage);
//		//double googleSVF = googlesvfcom->calSVF();
//		//delete dsImage;
//
//		//std::stringstream infofiless;
//		//infofiless << dataDir << i << "/info.txt";
//		//std::ifstream ifs(infofiless.str().outlineDS());
//		//std::string newline;
//		//std::getline(ifs, newline);
//		//ifs.close();
//		////uLTVrBeL8z9SiHLKQ8XbuA, -73.9689696451315, 40.7069884749583
//		//std::vector<std::string> infoAttributes = Utils::splitCSV(',', newline);
//		//double lon = atof(infoAttributes[1].outlineDS());
//		//double lat = atof(infoAttributes[2].outlineDS());
//		//OGRPoint pt;
//		//pt.setX(lon);
//		//pt.setY(lat);
//		OGRFeature *poFeatureNew = Shp.poLayer->GetNextFeature();
//		poFeatureNew->SetField(idIdx, i);
//		Shp.poLayer->SetFeature(poFeatureNew);
//		OGRFeature::DestroyFeature(poFeatureNew);
//		if (i % 100 == 0)
//		{
//			printf("%d\n", i);
//		}
//	}
//	Shp.close();
//}
void computeSVF_In_Manhattan()
{
	std::string dataDir = "E:/GoogleStreetview/Manhattan/Data/";
	std::string fisheyeDir = "E:/GoogleStreetview/Manhattan/SegnetFisheye2/";
	std::string outshapefile = "E:/GoogleStreetview/Manhattan/svf.shp";
	ShapeFile shp;
	shp.create(outshapefile, 0, 0, OGRwkbGeometryType::wkbPoint);
	int idIdx = shp.getOrCreateField("id", OGRFieldType::OFTInteger);
	int panoidIdx = shp.getOrCreateField("panoid", OGRFieldType::OFTString);
	int svfIdx = shp.getOrCreateField("svf", OGRFieldType::OFTReal);

	SVFComputer<unsigned char>* googlesvfcom = new SVFComputer<unsigned char>();
	for (int i = 0; i < 12528; i++)
	{
		std::stringstream fisheyefiless;
		fisheyefiless << fisheyeDir << i << ".png";
		std::string fisheyefile = fisheyefiless.str();
		if (!QFileInfo(fisheyefile.data()).exists())
			continue;
		GDAL_DS<unsigned char>* dsImage = new GDAL_DS<unsigned char>();
		dsImage->open(fisheyefile);
		googlesvfcom->setData(dsImage);
		double googleSVF = googlesvfcom->calSVF();
		delete dsImage;

		std::stringstream infofiless;
		infofiless << dataDir << i << "/info.txt";
		std::ifstream ifs(infofiless.str().data());
		std::string newline;
		std::getline(ifs, newline);
		ifs.close();
		//uLTVrBeL8z9SiHLKQ8XbuA, -73.9689696451315, 40.7069884749583
		std::vector<std::string> infoAttributes = Utils::splitCSV(',', newline);
		double lon = atof(infoAttributes[1].data());
		double lat = atof(infoAttributes[2].data());
		OGRPoint pt;
		pt.setX(lon);
		pt.setY(lat);
		OGRFeature *poFeatureNew = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		poFeatureNew->SetGeometry(&pt);
		poFeatureNew->SetField(panoidIdx, infoAttributes[0].data());
		poFeatureNew->SetField(svfIdx, googleSVF);
		poFeatureNew->SetField(idIdx, i);
		shp.poLayer->CreateFeature(poFeatureNew);
		OGRFeature::DestroyFeature(poFeatureNew);
		if (i % 100 == 0)
		{
			printf("%d\n", i);
		}
	}
	shp.close();
}

//};
void createShpFromCSV(std::string csvfielname = "E:/GoogleStreetview/CardiffGoogleStreet/points.csv", std::string outfilename = "E:/GoogleStreetview/CardiffGoogleStreet/points.shp")
{

	ShapeFile shp;
	shp.create(outfilename.data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);


	std::ifstream ifs;
	ifs.open(csvfielname.data());
	std::string headerline;
	std::getline(ifs, headerline);
	std::vector<std::string> splits = Utils::splitCSV(',', headerline);
	int idfield = shp.getOrCreateField("ID", OGRFieldType::OFTString);
	int numrecords = 0;
	while (ifs.peek() != -1)
	{
		std::string line;
		std::getline(ifs, line);
		splits = Utils::splitCSV(',', line);
		double xcoord = atof(splits[1].data());
		double ycoord = atof(splits[2].data());
		std::string id = splits[0];
		OGRFeature* poFeature = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		poFeature->SetField(idfield, id.data());
		OGRPoint pt;
		pt.setX(xcoord);
		pt.setY(ycoord);
		poFeature->SetGeometry(&pt);

		shp.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);

	}

}
void vintageStats(std::string filename, int& zero, int& pre1980, int& post1979)
{
	ShapeFile shp(filename);
	zero = 0;
	pre1980 = 0;
	post1979 = 0;
	int yearbuiltIdx = shp.getOrCreateField("yearbuilt", OGRFieldType::OFTReal);

	OGRFeatureDefn* layerdef = shp.poLayer->GetLayerDefn();

	shp.poLayer->ResetReading();
	OGRFeature *poFeature;

	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		int yearbuilt = poFeature->GetFieldAsInteger(yearbuiltIdx);
		if (yearbuilt == 0)
		{
			zero++;
		}
		else if (yearbuilt < 1980)
		{
			pre1980++;
		}
		else
		{
			post1979++;
		}
		OGRFeature::DestroyFeature(poFeature);
	}

}


struct VintageComposition
{
	std::string name;
	std::string sect;
	int total;
	int zero;
	int pre1980;
	int post1979;
	int sum()
	{
		total = zero + pre1980 + post1979;
		return total;
	}
};
void checkNCTootals(std::string dir)
{


	std::vector<std::string> files;

	QDir input_dir(dir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	std::string indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		if (!fileInfo.fileName().endsWith(".nc"))
			continue;
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();

		//printf("%f,%s\n", checkTotal(input_file), input_file.data());

	}

}


void grid100m()
{

	std::string rootdir = "B:/LA_Version2/Cumplots/";
	std::string indir = rootdir + "raw/";
	std::string boundfile = rootdir + "bound.txt";
	std::string fishnetshapefile = rootdir + "100m/fishnet.shp";
	std::string fishnetrasterfile = rootdir + "100m/fishnet.tif";
	std::string intersectedDir = rootdir + "100m/intersected/";
	std::string griddedDir = rootdir + "100m/gridded/";
	double meter2feet = 3.28084;
	double resol = meter2feet * 100;
	std::vector<std::string> files;

	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		OGREnvelope bound = BoundManager::readBound(boundfile);
		Grid* fishnetGrid = NULL;
		ShapeFile shp(input_file);
		char wkt[512];
		char* pwkt = wkt;
		if (shp.poLayer->GetSpatialRef())
			shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		shp.close();
		if (!QFileInfo(fishnetrasterfile.data()).exists())
		{
			fishnetGrid = new Grid(bound, resol);
			fishnetGrid->reset();
			fishnetGrid->toShape(pwkt, fishnetshapefile, false);
			fishnetGrid->toRaster(fishnetrasterfile, pwkt);
		}
		else
		{
			fishnetGrid = new Grid();
			fishnetGrid->fromFishnetRaster(fishnetrasterfile);
			fishnetGrid->reset();
		}

		std::string intersectedShapeFile = intersectedDir + fileInfo.fileName().toLocal8Bit().data();
		std::string griddedShapeFile = griddedDir + fileInfo.fileName().toLocal8Bit().data();

		if (!QFileInfo(intersectedShapeFile.data()).exists())
		{
			Utils::updateFootprint(input_file);
			Preprocessor::intersectWithArcGIS(input_file, fishnetshapefile, intersectedShapeFile);
			updateFieldAfterIntersection(intersectedShapeFile);
		}


		if (!QFileInfo(griddedShapeFile.data()).exists())
		{
			ShapeFile inshp(intersectedShapeFile);
			fishnetGrid->gatherCells(&inshp, "ca11");
			fishnetGrid->toShape(pwkt, griddedShapeFile, true);
		}

		delete fishnetGrid;
	}

}
void grid1km()
{

	std::string rootdir = "B:/LA_Version2/Cumplots/";
	std::string indir = rootdir + "raw/";
	std::string boundfile = rootdir + "bound.txt";
	std::string fishnetshapefile = rootdir + "1km/fishnet.shp";
	std::string fishnetrasterfile = rootdir + "1km/fishnet.tif";
	std::string intersectedDir = rootdir + "1km/intersected/";
	std::string griddedDir = rootdir + "1km/gridded/";
	double meter2feet = 3.28084;
	double resol = meter2feet * 1000;
	std::vector<std::string> files;

	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		OGREnvelope bound = BoundManager::readBound(boundfile);
		Grid* fishnetGrid = NULL;
		ShapeFile shp(input_file);
		char wkt[512];
		char* pwkt = wkt;
		if (shp.poLayer->GetSpatialRef())
			shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		shp.close();
		if (!QFileInfo(fishnetrasterfile.data()).exists())
		{
			fishnetGrid = new Grid(bound, resol);
			fishnetGrid->reset();
			fishnetGrid->toShape(pwkt, fishnetshapefile, false);
			fishnetGrid->toRaster(fishnetrasterfile, pwkt);
		}
		else
		{
			fishnetGrid = new Grid();
			fishnetGrid->fromFishnetRaster(fishnetrasterfile);
			fishnetGrid->reset();
		}

		std::string intersectedShapeFile = intersectedDir + fileInfo.fileName().toLocal8Bit().data();
		std::string griddedShapeFile = griddedDir + fileInfo.fileName().toLocal8Bit().data();

		if (!QFileInfo(intersectedShapeFile.data()).exists())
		{
			Utils::updateFootprint(input_file);
			Preprocessor::intersectWithArcGIS(input_file, fishnetshapefile, intersectedShapeFile);
			updateFieldAfterIntersection(intersectedShapeFile);
		}


		if (!QFileInfo(griddedShapeFile.data()).exists())
		{
			ShapeFile inshp(intersectedShapeFile);
			fishnetGrid->gatherCells(&inshp, "ca11");
			fishnetGrid->toShape(pwkt, griddedShapeFile, true);
		}

		delete fishnetGrid;
	}

}
void grid10km()
{

	std::string rootdir = "B:/LA_Version2/Cumplots/";
	std::string indir = rootdir + "raw/";
	std::string boundfile = rootdir + "bound.txt";
	std::string fishnetshapefile = rootdir + "10km/fishnet.shp";
	std::string fishnetrasterfile = rootdir + "10km/fishnet.tif";
	std::string intersectedDir = rootdir + "10km/intersected/";
	std::string griddedDir = rootdir + "10km/gridded/";
	double meter2feet = 3.28084;
	double resol = meter2feet * 10000;
	std::vector<std::string> files;

	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		OGREnvelope bound = BoundManager::readBound(boundfile);
		Grid* fishnetGrid = NULL;
		ShapeFile shp(input_file);
		char wkt[512];
		char* pwkt = wkt;
		if (shp.poLayer->GetSpatialRef())
			shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		shp.close();
		if (!QFileInfo(fishnetrasterfile.data()).exists())
		{
			fishnetGrid = new Grid(bound, resol);
			fishnetGrid->reset();
			fishnetGrid->toShape(pwkt, fishnetshapefile, false);
			fishnetGrid->toRaster(fishnetrasterfile, pwkt);
		}
		else
		{
			fishnetGrid = new Grid();
			fishnetGrid->fromFishnetRaster(fishnetrasterfile);
			fishnetGrid->reset();
		}

		std::string intersectedShapeFile = intersectedDir + fileInfo.fileName().toLocal8Bit().data();
		std::string griddedShapeFile = griddedDir + fileInfo.fileName().toLocal8Bit().data();

		if (!QFileInfo(intersectedShapeFile.data()).exists())
		{
			Utils::updateFootprint(input_file);
			Preprocessor::intersectWithArcGIS(input_file, fishnetshapefile, intersectedShapeFile);
			updateFieldAfterIntersection(intersectedShapeFile);
		}


		if (!QFileInfo(griddedShapeFile.data()).exists())
		{
			ShapeFile inshp(intersectedShapeFile);
			fishnetGrid->gatherCells(&inshp, "ca11");
			fishnetGrid->toShape(pwkt, griddedShapeFile, true);
		}

		delete fishnetGrid;
	}

}
void updateRail(std::string neifile, std::string srcfname, std::string valuefname, std::string railfile, std::string destfname)
{

	std::map<int, double> NEI_Emissions;
	std::ifstream ifs;
	ifs.open(neifile.data());
	std::string headerline;
	std::getline(ifs, headerline);
	std::vector<std::string> splits = Utils::splitCSV(',', headerline);
	splits = Utils::splitCSV(',', headerline);
	int srcid = -1;
	int valueid = -1;

	for (size_t i = 0; i < splits.size(); i++)
	{
		if (srcfname == splits[i])
		{
			srcid = i;
		}
		if (valuefname == splits[i])
		{
			valueid = i;
		}
	}


	while (ifs.peek() != -1)
	{
		std::string line;
		std::getline(ifs, line);
		splits = Utils::splitCSV(',', line);
		if (splits[srcid] == "NA" || splits[srcid] == "")
			continue;
		int id = atoi(splits[srcid].data());
		double emissions = atof(splits[valueid].data());
		NEI_Emissions[id] = emissions;
	}
	ifs.close();

	ShapeFile shp(railfile, 1);
	int idIdx = shp.getOrCreateField(destfname.data(), OGRFieldType::OFTInteger);
	int caIdx = shp.getOrCreateField("ca11", OGRFieldType::OFTReal);
	OGRFeatureDefn* layerdef = shp.poLayer->GetLayerDefn();
	shp.poLayer->ResetReading();
	OGRFeature *poFeature;
	int total = 0;
	int matches = 0;
	int nonzero = 0;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		int id = poFeature->GetFieldAsInteger(idIdx);
		std::map<int, double>::iterator iter = NEI_Emissions.find(id);
		if (iter != NEI_Emissions.end())
		{
			poFeature->SetField(caIdx, iter->second);
			shp.poLayer->SetFeature(poFeature);
			matches++;
			if (iter->second > 0)
				nonzero++;
		}
		else
		{
			//printf("%d\n", id);
		}
		total++;
		OGRFeature::DestroyFeature(poFeature);
	}
	printf("%d,%d,%d,%f\n", matches, nonzero, total, (double)matches / total);
	shp.close();
}
void gridRail(std::string railfile)
{
	std::string name = QFileInfo(railfile.data()).completeBaseName().toLocal8Bit().data();
	std::string dir = QFileInfo(railfile.data()).absoluteDir().absolutePath().toLocal8Bit().data() + std::string("/");
	OGREnvelope bb = BoundManager::readBoundFromShape(railfile);
	Grid grid(bb, 0.1);
	ShapeFile shp(railfile);
	char wkt[512];
	char* pwkt = wkt;
	if (shp.poLayer->GetSpatialRef())
		shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	shp.close();
	grid.reset();
	std::string fishnet = dir + +"fishnet.shp";
	if (!QFileInfo(fishnet.data()).exists())
	{
		grid.toShape(wkt, fishnet);
	}
	std::string outgrid = dir + +"grid.shp";

	std::string intersectedShapeFile = dir + "intersected.shp";

	if (!QFileInfo(intersectedShapeFile.data()).exists())
	{
		Utils::updateFootprint(railfile);
		Preprocessor::intersectWithArcGIS(railfile, fishnet, intersectedShapeFile);
		updateFieldAfterIntersection(intersectedShapeFile);
	}

	ShapeFile inshp(intersectedShapeFile);
	grid.gatherCells(&inshp, "ca11");
	grid.toShape(pwkt, outgrid, true);

}

void calSVF()
{

	GDAL_DS<double>* ds = new GDAL_DS<double>();
	ds->open("F:/Oblique_Photogrammetry/weihai/Weihai_DSM_025_UTM51N.tif");
	OGRSpatialReference oTargetSRS;
	char* buf = new char[ds->projection.size()];
	memcpy(buf, ds->projection.data(), ds->projection.size());
	oTargetSRS.importFromWkt(&buf);
	OGRSpatialReference oSourceSRS;
	oSourceSRS.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
	double x = -3.17970006;
	double y = 51.4789137;

	double destx = x;
	double desty = y;
	poCT->Transform(1, &destx, &desty);
	printf("(%f,%f) -> (%f,%f)\n", x, y, destx, desty);
	SVFComputer<double>* svfcom = new SVFComputer<double>(ds);

	std::string lidarFisheyeDir = "E:/GoogleStreetview/CardiffGoogleStreet2/LiDARFisheye/";
	std::string googleFisheyeDir = "E:/GoogleStreetview/CardiffGoogleStreet2/SegnetFisheye/";
	std::ofstream ssSVF("E:/GoogleStreetview/CardiffGoogleStreet2/svf.csv");
	ssSVF << "ID,PanoID,Lon,Lat,X,Y,GroundElev,Greenery,LiDAR,Google" << std::endl;
	SVFComputer<unsigned char>* googlesvfcom = new SVFComputer<unsigned char>();


	//double x = atof(splits[1].outlineDS());
	//double y = atof(splits[2].outlineDS());
	//double lon = x;
	//double lat = y;
	//poCT->Transform(1, &x, &y);
	////if (QDir((dir + id).outlineDS()).exists())
	////{
	//std::stringstream ssIdx;
	//ssIdx << idx << ".png";
	//std::string lidarfile = lidarFisheyeDir + ssIdx.str();
	//std::string googlefile = googleFisheyeDir + ssIdx.str();
	//double groundelev = svfcomMax->getHeightAt(x, y);
	//double eyeheight = groundelev + 1.4 + 1;

	//osg::Vec3d eye(x, y, eyeheight);
	//std::vector<osg::Vec2d> skymapcoords = svfcomMax->computeSkymap(eye, 3600);
	//svfcomMax->drawSkymap(skymapcoords, 512, lidarfile);
	//
	//GDAL_DS<unsigned char>* dsImage = new GDAL_DS<unsigned char>();
	//dsImage->open(lidarfile);
	//googlesvfcom->setData(dsImage);
	//double lidarSVF = googlesvfcom->calSVF();
}
void calSVF_Kechunag()
{

	GDAL_DS<double>* ds = new GDAL_DS<double>();
	ds->open("E:/ClimateSkyModel64/Kechuang10cm.tif");
	OGRSpatialReference oTargetSRS;
	char* buf = new char[ds->projection.size()];
	memcpy(buf, ds->projection.data(), ds->projection.size());

	std::string dsmFisheyeDir = "E:/SolarComp/bin/DSM/";
	////std::ofstream ssSVF("E:/GoogleStreetview/CardiffGoogleStreet2/svf.csv");
	//ssSVF << "ID,PanoID,Lon,Lat,X,Y,GroundElev,Greenery,LiDAR,Google" << std::endl;
	SVFComputer<double>* svfcom = new SVFComputer<double>(ds);
	ShapeFile shp("E:/SolarComp/bin/points_kechuang.shp");
	OGRFeature *poFeature;
	shp.poLayer->ResetReading();
	int pointid = 0;
	int idx = 0;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{

		OGRPoint* ogrp = (OGRPoint*)poFeature->GetGeometryRef();
		double x = ogrp->getX();
		double y = ogrp->getY();
		double groundelev = svfcom->getHeightAt(x, y);
		double eyeheight = groundelev + 0.1;
		osg::Vec3d pos(x, y, eyeheight);
		std::stringstream ssIdx;
		ssIdx << idx << ".png";
		std::string dsmfile = dsmFisheyeDir + ssIdx.str();
		std::vector<osg::Vec2d> skymapcoords = svfcom->computeSkymap(pos, 3600);
		svfcom->drawSkymap(skymapcoords, 1024, dsmfile);
		idx++;
	}

	//double x = atof(splits[1].outlineDS());
	//double y = atof(splits[2].outlineDS());
	//double lon = x;
	//double lat = y;
	//poCT->Transform(1, &x, &y);
	////if (QDir((dir + id).outlineDS()).exists())
	////{
	//std::stringstream ssIdx;
	//ssIdx << idx << ".png";
	//std::string lidarfile = lidarFisheyeDir + ssIdx.str();
	//std::string googlefile = googleFisheyeDir + ssIdx.str();
	//double groundelev = svfcomMax->getHeightAt(x, y);
	//double eyeheight = groundelev + 1.4 + 1;

	//osg::Vec3d eye(x, y, eyeheight);
	//std::vector<osg::Vec2d> skymapcoords = svfcomMax->computeSkymap(eye, 3600);
	//svfcomMax->drawSkymap(skymapcoords, 512, lidarfile);
	//
	//GDAL_DS<unsigned char>* dsImage = new GDAL_DS<unsigned char>();
	//dsImage->open(lidarfile);
	//googlesvfcom->setData(dsImage);
	//double lidarSVF = googlesvfcom->calSVF();
}
void combineGrids(std::string infile1, std::string infile2, std::string outfile)
{
	GDAL_DS<float>* dsfile1 = new GDAL_DS<float>();
	dsfile1->open(infile1);
	float* data1 = dsfile1->readData(1);
	delete dsfile1;

	GDAL_DS<float>* dsfile2 = new GDAL_DS<float>();
	dsfile2->open(infile2);
	float* data2 = dsfile2->readData(1);
	delete dsfile2;

	GDAL_DS<float>* dsoutfile = new GDAL_DS<float>();
	dsoutfile->open(infile1);
	dsoutfile->create(outfile);
	int numcells = dsoutfile->ncols * dsoutfile->nrows;
	for (int i = 0; i < numcells; i++)
	{
		if (data1[i] < 0)
			data1[i] = 0;
		if (data2[i] < 0)
			data2[i] = 0;
		data1[i] = data1[i] + data2[i];
	}
	dsoutfile->writeData(1, data1, 0);
	delete dsoutfile;
}
void subtractGrids(std::string infile1, std::string infile2, std::string outfile)
{
	GDAL_DS<float>* dsfile1 = new GDAL_DS<float>();
	dsfile1->open(infile1);
	float* data1 = dsfile1->readData(1);
	delete dsfile1;

	GDAL_DS<float>* dsfile2 = new GDAL_DS<float>();
	dsfile2->open(infile2);
	float* data2 = dsfile2->readData(1);
	delete dsfile2;

	GDAL_DS<float>* dsoutfile = new GDAL_DS<float>();
	dsoutfile->open(infile1);
	dsoutfile->create(outfile);
	int numcells = dsoutfile->ncols * dsoutfile->nrows;
	for (int i = 0; i < numcells; i++)
	{
		if (data1[i] < 0)
			data1[i] = 0;
		if (data2[i] < 0)
			data2[i] = 0;
		if (data1[i] > 0)
			data1[i] = (data2[i] - data1[i]) / data1[i];
		if (data1[i] > 1)
			data1[i] = 1;
		else if (data1[i] < -1)
			data1[i] = -1;
	}
	dsoutfile->writeData(1, data1, 0);
	delete dsoutfile;
}
void cleanGrid(std::string infile1, std::string outfile)
{
	GDAL_DS<float>* dsoutfile = new GDAL_DS<float>();
	dsoutfile->open(infile1);
	float* data = dsoutfile->readData(1);
	dsoutfile->create(outfile);
	int numcells = dsoutfile->ncols * dsoutfile->nrows;
	for (int i = 0; i < numcells; i++)
	{
		if (data[i] < 0)
			data[i] = 0;
	}
	dsoutfile->writeData(1, data, 0);
	delete dsoutfile;
}

void toCSV(std::string infile, std::string outcsvfile)
{
	GDAL_DS<float>* dsoutfile = new GDAL_DS<float>();
	dsoutfile->open(infile);
	float* data = dsoutfile->readData(1);
	int numcells = dsoutfile->ncols * dsoutfile->nrows;
	delete dsoutfile;
	std::ofstream ofs(outcsvfile.data());

	for (int i = 0; i < numcells; i++)
	{
		if (data[i] < 0)
			continue;
		ofs << data[i] << std::endl;
	}
	delete[] data;
	ofs.close();

}
void count(std::string filename)
{

	GDAL_DS<float>* dsoutfile = new GDAL_DS<float>();
	dsoutfile->open(filename);
	float* data = dsoutfile->readData(1);
	int numcells = dsoutfile->ncols * dsoutfile->nrows;
	delete dsoutfile;
	int greaterthan0 = 0;
	int greaterthan2000 = 0;
	for (int i = 0; i < numcells; i++)
	{
		if (data[i] < 1)
			continue;
		greaterthan0++;
		if (data[i] > 2000)
			greaterthan2000++;
	}
	delete[] data;
	printf("%d,%d,%d", greaterthan2000, greaterthan0, numcells);

}
void gridFolderByRaster(std::string rootdir, int resolution)
{

	//std::string rootdir = "B:/LA_Version2/Cumplots/";

	std::stringstream ssresol;
	ssresol << resolution << "m";
	std::string strresol = ssresol.str();
	std::string indir = rootdir + "raw/";
	std::string boundfile = rootdir + "bound.txt";
	std::string fishnetshapefile = rootdir + strresol + "/fishnet.shp";
	std::string fishnetrasterfile = rootdir + strresol + "/fishnet.tif";
	std::string intersectedDir = rootdir + strresol + "/intersected/";
	std::string griddedDir = rootdir + strresol + "/gridded/";
	double meter2feet = 3.28084;
	double resol = meter2feet * resolution;
	std::vector<std::string> files;

	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QDir(intersectedDir.data()).mkpath(".");
	//QDir(griddedDir.outlineDS()).mkpath(".");
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		OGREnvelope bound = BoundManager::readBound(boundfile);
		Grid* fishnetGrid = NULL;
		ShapeFile shp(input_file);
		char wkt[512];
		char* pwkt = wkt;
		if (shp.poLayer->GetSpatialRef())
			shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
		shp.close();
		if (!QFileInfo(fishnetrasterfile.data()).exists())
		{
			fishnetGrid = new Grid(bound, resol);
			fishnetGrid->reset();
			fishnetGrid->toShape(pwkt, fishnetshapefile, false);
			fishnetGrid->toRaster(fishnetrasterfile, pwkt);
		}
		else
		{
			fishnetGrid = new Grid();
			fishnetGrid->fromFishnetRaster(fishnetrasterfile);
			fishnetGrid->reset();
		}

		std::string intersectedShapeFile = intersectedDir + fileInfo.fileName().toLocal8Bit().data();
		std::string griddedShapeFile = griddedDir + fileInfo.fileName().toLocal8Bit().data();

		if (!QFileInfo(intersectedShapeFile.data()).exists())
		{
			Utils::updateFootprint(input_file);
			Preprocessor::intersectWithArcGIS(input_file, fishnetshapefile, intersectedShapeFile);
			updateFieldAfterIntersection(intersectedShapeFile);
		}


		//if (!QFileInfo(griddedShapeFile.outlineDS()).exists())
		//{
		//	ShapeFile inshp(intersectedShapeFile);
		//	fishnetGrid->gatherCells(&inshp, "ca11");
		//	fishnetGrid->toShape(pwkt, griddedShapeFile, true);
		//}

		delete fishnetGrid;
	}

}
void IntersectLA1000to100()
{
	//std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	std::vector<std::string> subdirs = Utils::findSubdirectories("B:/LA_Version2/gridPrep_SHP_master/");
	std::string outdir = "E:/HestiaUncertainty/LACounty/Raw/";
	//for (int resol = 100; resol <= 1000; resol+=100)
	//{
	int resol = 1000;
	gridFolderByRaster("E:/HestiaUncertainty/LACounty/", resol);
	//break;
	//}

	//std::string outdir = "C:/HestiaGridding/Los_Angeles/result/";
	//std::string indir = "C:/HestiaGridding/Los_Angeles/Gridding/intersected_nogeometry/";
	//std::string timedir = "C:/HestiaGridding/Los_Angeles/Time/";

	//export2PUBLIC();
	///*subdirs = Utils::findSubdirectories("B:/LA_Version2/gridPrep_SHP_master/");
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	updateNonpointTime(subdirs[i]);
	//}
	//subdirs = Utils::findSubdirectories("B:/LA_Version2/ca/");
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	updateNonpointTime(subdirs[i]);
	//}

	//subdirs = Utils::findSubdirectories("B:/LA_Version2/gridoutput/");
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	updateNonpointTime(subdirs[i]);
	//}*/
	//indir = "B:/LA_Version2/gridoutput/";
	//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[6]{"ca","ca_ng","ca_p","ca_c","area","length"}, 6);
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
	//	//std::string srcFile = subdirs[i] + "OnRoad.Shp";
	//	//std::string destFile = "B:/LA_Version2/ca/" + foldername + "/" + "OnRoad.Shp";
	//	//if (!QDir(("B:/LA_Version2/ca/" + foldername + "/").outlineDS()).exists())
	//	//{
	//	//	QDir(("B:/LA_Version2/ca/" + foldername + "/").outlineDS()).mkpath(".");
	//	//}
	//	ShapeFile::copyDir(subdirs[i], outdir + foldername + "/", fields2keep);
	//	//Utils::updateFootprint(destFile);
	//	// Utils::updateFootprintForDir(subdirs[i], true);
	//	//ShapeFile::copyDir(subdirs[i], "B:/LA_Version2/ca/" + foldername + "/", fields2keep);
	//	//Preprocessor::gridFolderByRaster("B:/LA_Version2/ca/" + foldername + "/", "C:/HestiaGridding/Los_Angeles/gridding/intersected/" + foldername + "/", "C:/HestiaGridding/Los_Angeles/gridding/intersected/fishnet.tif");
	//	//ShapeFile::copyDirDropGeometry("C:/HestiaGridding/Los_Angeles/gridding/intersected/" + foldername + "/", indir + foldername + "/",fields2keep);
	//}
	//return;

}

void cropLA()
{
	char wkt[512];
	char* pwkt = wkt;
	printf("%d, %d,%d\n", 12, 12 / 5, 12 % 5);
	printf("%d, %d,%d\n", 18, 18 / 5, 18 % 5);
	printf("%d, %d,%d\n", 20, 20 / 5, 20 % 5);

	printf("%d, %d,%d\n", 20, 20 / 5, 20 % 5);
	//int maxrow = 2220799 / 1623;
	//int maxcol = 1600;


	Grid grid;
	grid.fromFishnetRaster("E:/HestiaUncertainty/LACounty/100m/fishnet.tif");
	int left = grid.ncols - 1; int right = 0;
	int top = grid.nrows - 1; int bottom = 0;
	int ncols = grid.ncols;
	int rows = grid.nrows;
	QDir input_dir("E:/HestiaUncertainty/LACounty100m/Intersected/");
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		ShapeFile shp(input_file.data());
		OGRFeatureDefn* layerdef = shp.poLayer->GetLayerDefn();
		int idIdx = layerdef->GetFieldIndex("Id");
		shp.poLayer->ResetReading();
		OGRFeature *poFeature;
		while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
		{
			int id = poFeature->GetFieldAsInteger(idIdx);
			int col = id % ncols;
			int row = id / ncols;
			if (row < 2460)
			{
				if (left > col)
					left = col;
				if (right < col)
					right = col;
				if (top > row)
					top = row;
				if (bottom < row)
					bottom = row;
			}

			OGRFeature::DestroyFeature(poFeature);
		}
	}

	grid.ncols = right - left + 1;
	grid.nrows = bottom - top + 1;
	grid.bound.MinX = grid.bound.MinX + left * grid._adfGeoTransform[1];
	grid.bound.MaxY = grid.bound.MaxY - top * grid._adfGeoTransform[1];
	grid._adfGeoTransform[0] = grid.bound.MinX;
	grid._adfGeoTransform[3] = grid.bound.MaxY;
	grid.toShape(grid.proj, "E:/HestiaUncertainty/LACounty100m/fishnet.shp");
	grid.toRaster("E:/HestiaUncertainty/LACounty100m/fishnet.tif", grid.proj);
}

struct OriInfo
{
	double footprint;
	double ca11;
	double ca11_ng;
	double ca11_p;
	double ca11_c;
};

class OriInfoAccessor
{

public:
	int  footprint_Idx;
	int  ca11_Idx;
	int  ca11_ng_Idx;
	int  ca11_p_Idx;
	int  ca11_c_Idx;

	int  footprint_Idx_dest;
	int  id_Idx_dest;
	int  ca11_Idx_dest;
	int  ca11_ng_Idx_dest;
	int  ca11_p_Idx_dest;
	int  ca11_c_Idx_dest;
	OGRFeatureDefn* layerdef;
	std::map<int, OriInfo> attributeTB;
	OGRFeatureDefn* pParentLayerdef;
	OriInfoAccessor(OGRFeatureDefn* layerdef)
	{
		pParentLayerdef = layerdef;
		footprint_Idx = layerdef->GetFieldIndex("area");
		if (footprint_Idx < 0)
		{
			footprint_Idx = layerdef->GetFieldIndex("length");
		}

		if (ca11_Idx < 0)
			ca11_Idx = layerdef->GetFieldIndex("ca11");
		if (ca11_ng_Idx < 0)
			ca11_ng_Idx = layerdef->GetFieldIndex("ca11_ng");
		if (ca11_p_Idx < 0)
			ca11_p_Idx = layerdef->GetFieldIndex("ca11_p");
		if (ca11_c_Idx < 0)
			ca11_c_Idx = layerdef->GetFieldIndex("ca11_c");
		footprint_Idx_dest = -1;
		ca11_Idx_dest = -1;
		ca11_ng_Idx_dest = -1;
		ca11_p_Idx_dest = -1;
		ca11_c_Idx_dest = -1;
	}
	void updateAttributeTB(OGRLayer* layer)
	{
		id_Idx_dest = layer->GetLayerDefn()->GetFieldIndex("Id");
		if (footprint_Idx > -1)
		{
			footprint_Idx_dest = layer->CreateField(pParentLayerdef->GetFieldDefn(footprint_Idx));
		}
		if (ca11_Idx > -1)
		{
			ca11_Idx_dest = layer->CreateField(pParentLayerdef->GetFieldDefn(ca11_Idx));
		}
		if (ca11_ng_Idx > -1)
		{
			ca11_ng_Idx_dest = layer->CreateField(pParentLayerdef->GetFieldDefn(ca11_ng_Idx));
		}
		if (ca11_p_Idx > -1)
		{
			ca11_p_Idx_dest = layer->CreateField(pParentLayerdef->GetFieldDefn(ca11_p_Idx));
		}
		if (ca11_c_Idx > -1)
		{
			ca11_c_Idx_dest = layer->CreateField(pParentLayerdef->GetFieldDefn(ca11_c_Idx));
		}

	}

	void fetchAttributes(int id, OGRFeature *poFeature)
	{
		OriInfo info;
		if (footprint_Idx > 0)
			info.footprint = poFeature->GetFieldAsDouble(footprint_Idx);
		if (ca11_Idx < 0)
			info.ca11 = poFeature->GetFieldAsDouble(ca11_Idx);
		if (ca11_ng_Idx < 0)
			info.ca11_ng = poFeature->GetFieldAsDouble(ca11_ng_Idx);
		if (ca11_p_Idx < 0)
			info.ca11_p = poFeature->GetFieldAsDouble(ca11_p_Idx);
		if (ca11_c_Idx < 0)
			info.ca11_c = poFeature->GetFieldAsDouble(ca11_c_Idx);
		attributeTB[id] = info;
	}
	void setAttributes(OGRFeature *poFeature)
	{
		OriInfo info;
		int id = poFeature->GetFieldAsInteger(id_Idx_dest);
		std::map<int, OriInfo>::iterator iter = attributeTB.find(id);
		if (iter != attributeTB.end())
		{
			if (footprint_Idx > 0)
				poFeature->SetField(footprint_Idx_dest, iter->second.footprint);
			if (ca11_Idx < 0)
				poFeature->SetField(ca11_Idx_dest, iter->second.ca11);
			if (ca11_ng_Idx < 0)
				poFeature->SetField(ca11_ng_Idx_dest, iter->second.ca11_ng);
			if (ca11_p_Idx < 0)
				poFeature->SetField(ca11_p_Idx_dest, iter->second.ca11_p);
			if (ca11_c_Idx < 0)
				poFeature->SetField(ca11_c_Idx_dest, iter->second.ca11_c);
		}
		else
		{
			printf("");
		}
	}
};
//B:\LA_Version2\gridPrep_SHP_master\Los_Angeles
void linkIntersected2Original(std::string intersectedDir, std::string originalDir)
{
	QDir input_dir(intersectedDir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string intersectedfile = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		name = name.substr(3, name.length() - 3);
		std::string orifile = originalDir + name;
		ShapeFile ori_shp(orifile.data());
		OGRFeatureDefn* ori_layerdef = ori_shp.poLayer->GetLayerDefn();
		ShapeFile intersected_shp(intersectedfile.data(), 1);
		OGRFeatureDefn* intersected_layerdef = intersected_shp.poLayer->GetLayerDefn();
		OriInfoAccessor accessor(ori_layerdef);
		accessor.updateAttributeTB(intersected_shp.poLayer);

		ori_shp.poLayer->ResetReading();

		OGRFeature *poFeature;
		int id = 0;
		while ((poFeature = ori_shp.poLayer->GetNextFeature()) != NULL)
		{
			accessor.fetchAttributes(id, poFeature);
			id++;
			OGRFeature::DestroyFeature(poFeature);
		}
		intersected_shp.poLayer->ResetReading();
		while ((poFeature = intersected_shp.poLayer->GetNextFeature()) != NULL)
		{
			accessor.setAttributes(poFeature);
			OGRFeature::DestroyFeature(poFeature);
		}
		intersected_shp.close();
		ori_shp.close();
	}
}

//void segments2points(std::string yuyuShapefile, std::string outshpfile, std::string outcsvfile)
//{
//	ShapeFile Shp(yuyuShapefile, 0);
//	ShapeFile outshp;
//	outshp.create(outshpfile, Shp.poLayer->GetSpatialRef(), 0, wkbPoint);
//	int segidx = outshp.getOrCreateField("oriID", OGRFieldType::OFTInteger);
//	OGRFeature *poFeatureNew;
//	OGRFeature *poFeature;
//	int oriid = 0;
//
//	while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
//	{
//		OGRGeometry* geom = poFeature->GetGeometryRef();
//		OGRPoint pt;
//		OGRwkbGeometryType gtype = poFeature->GetGeometryRef()->getGeometryType();
//
//		if (gtype == wkbLineString || gtype == wkbLineString25D)
//		{
//			OGRLineString* line = (OGRLineString*)geom;
//			int len = line->getNumPoints();
//			for (size_t i = 0; i < len; i++)
//			{
//				OGRPoint p;
//				line->getPoint(i, &p);
//				poFeatureNew = OGRFeature::CreateFeature(outshp.poLayer->GetLayerDefn());
//				poFeatureNew->SetGeometry(&p);
//				poFeatureNew->SetField(segidx, oriid);
//				outshp.poLayer->CreateFeature(poFeatureNew);
//				OGRFeature::DestroyFeature(poFeatureNew);
//			}
//		}
//		else if (gtype == wkbMultiLineString || gtype == wkbMultiLineString25D)
//		{
//
//			OGRMultiLineString *multipoly = (OGRMultiLineString*)geom;
//
//			for (size_t n = 0; n < multipoly->getNumGeometries(); n++)
//			{
//				OGRLineString* line = (OGRLineString*)multipoly->getGeometryRef(n);
//				int len = line->getNumPoints();
//				for (size_t i = 0; i < len; i++)
//				{
//					OGRPoint p;
//					line->getPoint(i, &p);
//					poFeatureNew = OGRFeature::CreateFeature(outshp.poLayer->GetLayerDefn());
//					poFeatureNew->SetGeometry(&p);
//					poFeatureNew->SetField(segidx, oriid);
//					outshp.poLayer->CreateFeature(poFeatureNew);
//					OGRFeature::DestroyFeature(poFeatureNew);
//				}
//			}
//		}
//		oriid++;
//		OGRFeature::DestroyFeature(poFeature);
//	}
//
//}

void segments2points(std::string shapefile, std::string outshpfile, std::string outcsvfile)
{
	ShapeFile shp(shapefile, 0);
	ShapeFile outshp;
	outshp.create(outshpfile, shp.poLayer->GetSpatialRef(), 0, wkbPoint);
	int segidx = outshp.getOrCreateField("oriID", OGRFieldType::OFTInteger);
	OGRFeature *poFeatureNew;
	OGRFeature *poFeature;
	int oriid = 0;

	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* geom = poFeature->GetGeometryRef();
		OGRPoint pt;
		OGRwkbGeometryType gtype = poFeature->GetGeometryRef()->getGeometryType();

		if (gtype == wkbLineString || gtype == wkbLineString25D)
		{
			OGRLineString* line = (OGRLineString*)geom;
			int len = line->getNumPoints();
			for (size_t i = 0; i < len; i++)
			{
				OGRPoint p;
				line->getPoint(i, &p);
				poFeatureNew = OGRFeature::CreateFeature(outshp.poLayer->GetLayerDefn());
				poFeatureNew->SetGeometry(&p);
				poFeatureNew->SetField(segidx, oriid);
				outshp.poLayer->CreateFeature(poFeatureNew);
				OGRFeature::DestroyFeature(poFeatureNew);
			}
		}
		else if (gtype == wkbMultiLineString || gtype == wkbMultiLineString25D)
		{

			OGRMultiLineString *multipoly = (OGRMultiLineString*)geom;

			for (size_t n = 0; n < multipoly->getNumGeometries(); n++)
			{
				OGRLineString* line = (OGRLineString*)multipoly->getGeometryRef(n);
				int len = line->getNumPoints();
				for (size_t i = 0; i < len; i++)
				{
					OGRPoint p;
					line->getPoint(i, &p);
					poFeatureNew = OGRFeature::CreateFeature(outshp.poLayer->GetLayerDefn());
					poFeatureNew->SetGeometry(&p);
					poFeatureNew->SetField(segidx, oriid);
					outshp.poLayer->CreateFeature(poFeatureNew);
					OGRFeature::DestroyFeature(poFeatureNew);
				}
			}
		}
		oriid++;
		OGRFeature::DestroyFeature(poFeature);
	}

}

void exportParcelURL()
{
	QDir input_dir("E:/DC_Corridor/Parcels/Maryland/");
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string outfile = fileInfo.absoluteFilePath().toLocal8Bit().data() + std::string(".csv");
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		ShapeFile shp(fileInfo.absoluteFilePath().toLocal8Bit().data(), 0);
		int id = shp.poLayer->GetLayerDefn()->GetFieldIndex("SDATWEBADR");
		OGRFeature *poFeature;
		std::ofstream ofs(outfile.data());
		while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
		{
			std::string url = poFeature->GetFieldAsString(id);
			ofs << url << std::endl;
		}
		ofs.close();
		OGRFeature::DestroyFeature(poFeature);
	}

}
#include "Vulcan2014.h"
void mergeRoads(std::string indir, std::string outfile, std::string roadtype)
{
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	ShapeFile* master = NULL;
	int srcSEGID = -1;
	int destSEGID = -1;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string fullname = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (!fileInfo.fileName().contains(roadtype.data()))
			continue;

		ShapeFile shp(fullname.data());

		if (master == NULL)
		{
			master = new ShapeFile;
			master->create(outfile, shp.poLayer->GetSpatialRef(), NULL, shp.poLayer->GetGeomType());
			OGRFieldDefn newField("SEGID", OGRFieldType::OFTString);
			destSEGID = master->poLayer->CreateField(&newField);
			master->poLayer->GetLayerDefn()->GetFieldDefn(destSEGID)->SetWidth(12);
		}
		srcSEGID = shp.poLayer->GetLayerDefn()->GetFieldIndex("SEGID");
		shp.poLayer->ResetReading();
		OGRFeature *poFeature;
		int id = 0;
		while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
		{
			OGRFeature *poFeatureNew = OGRFeature::CreateFeature(master->poLayer->GetLayerDefn());
			poFeatureNew->SetGeometry(poFeature->GetGeometryRef());
			poFeatureNew->SetField(destSEGID, poFeature->GetFieldAsString(srcSEGID));
			master->poLayer->CreateFeature(poFeatureNew);
			OGRFeature::DestroyFeature(poFeatureNew);
			OGRFeature::DestroyFeature(poFeature);
		}
		shp.close();
		printf("%s\n", name.data());
	}
	delete master;
}


void reprojectVulcan()
{
	ShapeFile shp("e:/Vulcan/gridPrep_SHP_master/OnroadRuralLocal.shp");
	char wkt[512];
	char* pwkt = wkt;
	if (shp.poLayer->GetSpatialRef())
		shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	std::string destwktfile = "e:/Vulcan/destwktfile.txt";
	std::ofstream ofs;
	ofs.open(destwktfile);
	ofs << std::string(pwkt);
	ofs.close();
	//ShapeFile::reprojectDir("e:/Vulcan/gridPrep_SHP_master/", "e:/Vulcan/gridPrep_SHP_master/reprojected/", destwktfile);

	ShapeFile::reproject("E:/Vulcan/fishnet/ContiguousStates.shp", "E:/Vulcan/fishnet/ContiguousStates_Lambert.shp", destwktfile);

}

void calVulcanExtent()
{
	OGREnvelope initBound = BoundManager::readBoundFromShape("E:/Vulcan/fishnet/fema_blockgroups.shp");

	QDir input_dir("e:/Vulcan/gridPrep_SHP_master/");
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();

	double tolerance = 1000 * 100;
	OGREnvelope bound = initBound;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string fullname = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (fileInfo.fileName().contains("CMV"))
			continue;
		OGREnvelope shpbound = BoundManager::readBoundFromShape(fullname);
		if (shpbound.MaxX > initBound.MaxX && shpbound.MaxX - initBound.MaxX < tolerance)
		{
			if (bound.MaxX < shpbound.MaxX)
				bound.MaxX = shpbound.MaxX;
		}

		if (shpbound.MaxY > initBound.MaxY && shpbound.MaxY - initBound.MaxY < tolerance)
		{
			if (bound.MaxY < shpbound.MaxX)
				bound.MaxY = shpbound.MaxY;
		}

		if (shpbound.MinX < initBound.MinX && abs(shpbound.MinX - initBound.MinX) < tolerance)
		{
			if (bound.MinX > shpbound.MinX)
				bound.MinX = shpbound.MinX;
		}

		if (shpbound.MinY < initBound.MinY && abs(shpbound.MinY - initBound.MinY) < tolerance)
		{
			if (bound.MinY > shpbound.MinY)
				bound.MinY = shpbound.MinY;
		}
	}

	BoundManager::writeBound(bound, "E:/Vulcan/fishnet/VulcanExtent.txt");
	ShapeFile shp("E:/Vulcan/fishnet/fema_blockgroups.shp");
	char wkt[512];
	char* pwkt = wkt;
	if (shp.poLayer->GetSpatialRef())
		shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	Grid* grid = new Grid(bound, 1000, 1);
	grid->toBoundaryShape("E:/Vulcan/fishnet/VulcanExtent.shp", "E:/Vulcan/fishnet/fema_blockgroups.shp");
	grid->toRaster("E:/Vulcan/fishnet/fishnet.tif", pwkt);
	grid->toShape(pwkt, "E:/Vulcan/fishnet/fishnet.shp", false);
}

OGRPolygon* createCircle(int numslices, double radius)
{
	OGRPolygon* rings = new OGRPolygon();
	OGRLinearRing* ring = new OGRLinearRing();
	double step = 1 / (double)numslices * 360.0;
	for (int i = 0; i < numslices; i++)
	{
		double angle = step * i * 0.0174533;
		OGRPoint pt(cos(angle) * radius, sin(angle) * radius);
		printf("%f,%f\n", pt.getX(), pt.getY());
		ring->addPoint(&pt);
	}
	OGRPoint firstp;
	ring->getPoint(0, &firstp);
	ring->addPoint(&firstp);
	rings->addRing(ring);
	return rings;
}

void createCircle(int numslices, double radius, OGRLayer* layer)
{
	OGRPolygon* rings = createCircle(numslices, radius);
	OGRFeature *poFeatureNew = OGRFeature::CreateFeature(layer->GetLayerDefn());
	poFeatureNew->SetGeometry(rings);
	layer->CreateFeature(poFeatureNew);
	OGRFeature::DestroyFeature(poFeatureNew);
}

OGRPolygon* createBox(OGREnvelope bb)
{
	OGRPolygon *poPolygon = new OGRPolygon(); //(OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
	OGRLinearRing  *linearRing = new OGRLinearRing(); //(OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);
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
#include "TemporalGridderByShapes.h"
//OK, let��s go ahead and make the netCDF.We shall call this the ��WRF.VY.*****�� to indicate that this is the WRF grid for Vineet Yadav.
void gridLA_VY()
{

	std::string rootdir = "B:/LA_Version2/Vulcan_output/Vineet/";
	std::string outdir = rootdir + "GriddedEmissions/";
	std::string shapesindir = "B:/LA_Version2/Vulcan_output/Spatial/";
	std::string intersected = rootdir + "intersected_tmp/";
	std::string indir = rootdir + "intersected/";
	std::string timedir = "B:/LA_Version2/Vulcan_output/time/";
	std::string gridShapeFile = rootdir + "fishnet.shp";
	QDir(rootdir.data()).mkdir(".");
	QDir(outdir.data()).mkdir(".");
	QDir(shapesindir.data()).mkdir(".");
	QDir(intersected.data()).mkdir(".");
	QDir(indir.data()).mkdir(".");

	Utils::updateFootprintForDir(shapesindir, false);
	OGRSpatialReference oSRS;
	oSRS.SetWellKnownGeogCS("WGS84");
	ShapeFile shp;
	shp.create(gridShapeFile, &oSRS, NULL, OGRwkbGeometryType::wkbPolygon);
	std::string line;
	std::ifstream ifs;
	ifs.open(rootdir + "GRID_SCB_0_01.csv");
	std::string headerline;
	std::getline(ifs, headerline);
	int latfield = shp.getOrCreateField("lat", OGRFieldType::OFTReal);
	int lonfield = shp.getOrCreateField("lon", OGRFieldType::OFTReal);
	int idfield = shp.getOrCreateField("Id", OGRFieldType::OFTInteger);
	double halfcellsize = 0.01 * 0.5;
	int numcells = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		double lat = atof(splits[1].data());
		double lon = atof(splits[2].data());
		OGRFeature* poFeature = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		OGREnvelope gridcell;
		gridcell.MaxX = lon + halfcellsize;
		gridcell.MaxY = lat + halfcellsize;
		gridcell.MinX = lon - halfcellsize;
		gridcell.MinY = lat - halfcellsize;
		OGRGeometry* cellpoly = createBox(gridcell);
		poFeature->SetGeometry(cellpoly);
		poFeature->SetField(latfield, lat);
		poFeature->SetField(lonfield, lon);
		poFeature->SetField(idfield, numcells);
		shp.poLayer->CreateFeature(poFeature);
		//OGRGeometryFactory::destroyGeometry(cellpoly);
		OGRFeature::DestroyFeature(poFeature);
		numcells++;
	}
	ifs.close();
	shp.close();

	Preprocessor::gridFolderByShape(shapesindir, intersected, gridShapeFile);
	std::vector<std::string> fields2keep;
	ShapeFile::copyDirDropGeometry(intersected, indir, fields2keep);

	std::string years[]{ "2010","2011","2012","2013","2014","2015" };
	std::string cafields[]{ "ca10","ca11","ca12","ca13","ca14","ca15" };
	std::map<std::string, std::vector<std::string>> timestructfile_crosswalk;

	std::vector<GridSectorConfig> sectors;
	GridSectorConfig com("com", Utils::buildVector(new std::string[2]{ "comnonpoint","compoint" }, 2),
		Utils::buildVector(new std::string[2]{ "NonPoint","SMOKE" }, 2));
	GridSectorConfig res("res", Utils::buildVector(new std::string[1]{ "resnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "NonPoint" }, 1));
	GridSectorConfig ind("ind", Utils::buildVector(new std::string[2]{ "indnonpoint","indpoint" }, 2),
		Utils::buildVector(new std::string[2]{ "NonPoint","SMOKE" }, 2));
	GridSectorConfig railroad("railroad", Utils::buildVector(new std::string[1]{ "Railroad" }, 1),
		Utils::buildVector(new std::string[1]{ "SMOKE" }, 1));
	GridSectorConfig nonroad("nonroad", Utils::buildVector(new std::string[1]{ "NonRoad" }, 1),
		Utils::buildVector(new std::string[2]{ "NonRoad","SMOKE" }, 2));
	GridSectorConfig onroad("onroad", Utils::buildVector(new std::string[1]{ "OnRoad.shp" }, 1),
		Utils::buildVector(new std::string[1]{ "OnRoad" }, 1));
	GridSectorConfig cmv("cmv", Utils::buildVector(new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2),
		Utils::buildVector(new std::string[1]{ "cmv" }, 1));
	GridSectorConfig airport("airport", Utils::buildVector(new std::string[1]{ "airport" }, 1),
		Utils::buildVector(new std::string[1]{ "airport" }, 1));
	GridSectorConfig elecprod("elecprod", Utils::buildVector(new std::string[2]{ "ElecProd","ElecNEI" }, 2),
		Utils::buildVector(new std::string[2]{ "ElecProd" ,"SMOKE" }, 2));
	GridSectorConfig cement("cement", Utils::buildVector(new std::string[1]{ "cement" }, 1));


	sectors.push_back(onroad);
	sectors.push_back(com);
	sectors.push_back(res);
	sectors.push_back(ind);
	sectors.push_back(railroad);
	sectors.push_back(nonroad);
	sectors.push_back(elecprod);
	sectors.push_back(airport);
	sectors.push_back(cmv);
	sectors.push_back(cement);

	QDir(outdir.data()).mkpath(".");
	std::string cityname = "LAbasin";
	std::string version = "v2.5";
	std::string projectName = "WRF.VY";

	for (int i = 0; i < 6; i++)
	{
		std::string year = years[i];
		int numhous = 8760;
		if (year == "2012")
		{
			numhous = 8784;
		}
		TemporalGridderByShapes gridder(numhous);
		gridder.numcells = numcells;
		gridder.cells = new double[gridder.numcells];
		memset(gridder.cells, 0, gridder.numcells*sizeof(double));
		std::string timestructdir = timedir + year + "/";
		SubdirManager dirmanager(indir);
		for (int nsector = 0; nsector < sectors.size(); nsector++)
		{
			std::vector<std::string> shapefiles = dirmanager.findFilesMatch(sectors[nsector].shapefiles);
			gridder.addSectorGrid(shapefiles, sectors[nsector].sectorname);
			for (int ntimestruct = 0; ntimestruct < sectors[nsector].timestructfiles.size(); ntimestruct++)
			{
				gridder.loadtimestruct(timestructdir + sectors[nsector].timestructfiles[ntimestruct] + ".bin");
			}

		}

		gridder.loadAttribute(cafields[i]);
		std::string annualFile = outdir + projectName + ".total.annual." + year + ".bin";
		std::string hourlyFile = outdir + projectName + ".total.hourly." + year + ".bin";
		gridder.makeAnnualTotal(annualFile);
		for (size_t isector = 0; isector < gridder.sectors.size(); isector++)
		{
			FFCO2Sector* sectorGrid = gridder.sectors[isector];
			std::string sectorAnnualFile = outdir + projectName + "." + sectorGrid->sectorname + ".annual." + year + ".bin";
			std::string sectorHourlyFile = outdir + projectName + "." + sectorGrid->sectorname + ".hourly." + year + ".bin";
			sectorGrid->makeAnnualTotal(sectorAnnualFile);
			sectorGrid->makeHourlyTotal(sectorHourlyFile);
		}
		gridder.makeHourlyTotal(hourlyFile);
		for (size_t isector = 0; isector < gridder.sectors.size(); isector++)
		{
			FFCO2Sector* sectorGrid = gridder.sectors[isector];
		}
		gridder.clearSectors();
	}


	////��WRF.VY.*****�� to indicate that this is the WRF grid for Vineet Yadav
	////std::string sourceDir = "B:/LA_Version2/gridPrep_SHP_master/";
	////std::string caDir = "C:/HestiaGridding/Los_Angeles/Vineet/ca/";
	////std::string indir = "C:/HestiaGridding/Los_Angeles/Vineet/intersected/";
	//std::string timedir = "C:/HestiaGridding/Los_Angeles/Time/";
	//std::string intersectedDirNoGeometry = "C:/HestiaGridding/Los_Angeles/Vineet/intersected_nogeometry/";
	//std::string outdir = "C:/HestiaGridding/Los_Angeles/Vineet/result/";
	///*QDir(caDir.outlineDS()).mkpath(".");
	//QDir(indir.outlineDS()).mkpath(".");*/
	//QDir(intersectedDirNoGeometry.data()).mkpath(".");
	//QDir(outdir.data()).mkpath(".");
	//std::string gridShapeFile = intersectedDirNoGeometry + "fishnet.shp";
	//std::string projectName = "WRF.VY";
	//std::string year = "2014";
	//std::string cafield = "ca" + year.substr(2, 2);
	//std::string timestructfile = timedir + "timestructs_" + year + ".bin";
	////std::vector<std::string> subdirs = Utils::findSubdirectories(sourceDir);
	////std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[9]{ "ca10","ca11","ca12","ca13","ca14","length","area" ,"timestruct","Id" }, 9);
	////for (size_t i = 0; i < subdirs.size(); i++)
	////{
	////	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
	////	//Utils::updateFootprint(destFile);
	////	// Utils::updateFootprintForDir(subdirs[i], false);
	////	//ShapeFile::copyDir(subdirs[i], caDir + foldername + "/", fields2keep);
	////	//Preprocessor::gridFolderByShape(caDir + foldername + "/", indir + foldername + "/", gridShapeFile);
	////	//ShapeFile::copyDirDropGeometry(indir + foldername + "/", intersectedDirNoGeometry + foldername + "/", fields2keep);
	////}

	//SubdirManager dirmanager = SubdirManager(intersectedDirNoGeometry);
	//std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	//std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	//std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	//std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "RailroadPoint","Railroad" }, 2));
	//std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.shp" }, 1));
	//std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "NonRoad" }, 1));
	//std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	//std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	//std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));




	//TemporalGridderByShapes gridder;
	//gridder.numcells = 1826;
	//gridder.cells = new double[gridder.numcells];
	//memset(gridder.cells, 0, gridder.numcells*sizeof(double));
	//gridder.loadtimestruct("", timestructfile,8760);
	//gridder.addSectorGrid(comshapefiles, "com");
	//gridder.addSectorGrid(indshapefiles, "ind");
	//gridder.addSectorGrid(resshapefiles, "res");
	//gridder.addSectorGrid(railroadshapefiles, "railroad");
	//gridder.addSectorGrid(onroadshapefiles, "onroad");
	//gridder.addSectorGrid(nonroadshapefiles, "nonroad");
	//gridder.addSectorGrid(elecprodshapefiles, "elecprod");
	//gridder.addSectorGrid(marineshapefiles, "marine");
	//gridder.addSectorGrid(airportshapefiles, "airport");
	//gridder.loadAttribute(cafield);


	//std::string annualFile = outdir + projectName + ".total.annual." + year + ".bin";
	//std::string hourlyFile = outdir + projectName + ".total.hourly." + year + ".bin";
	//gridder.makeAnnualTotal(annualFile);

	//for (size_t isector = 0; isector < gridder.sectors.size(); isector++)
	//{
	//	FFCO2Sector* sectorGrid = gridder.sectors[isector];
	//	std::string sectorAnnualFile = outdir + projectName + "." + sectorGrid->sectorname + ".annual." + year + ".bin";
	//	std::string sectorHourlyFile = outdir + projectName + "." + sectorGrid->sectorname + ".hourly." + year + ".bin";
	//	sectorGrid->makeAnnualTotal(sectorAnnualFile);
	//	sectorGrid->makeHourlyTotal(sectorHourlyFile);
	//}

	//gridder.makeHourlyTotal(hourlyFile);
	//Utils::updateFootprint("C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp", 110959.0* 110959.0);
	//ShapeFile* file = new ShapeFile("C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp", 1);
	//OGRFeature *poFeature;
	//std::ofstream ofs;
	//ofs.open("VYgrids.csv");
	//while ((poFeature = file->poLayer->GetNextFeature()) != NULL)
	//{
	//	OGRPolygon* polygon = (OGRPolygon*)poFeature->GetGeometryRef();
	//	OGREnvelope bound;
	//	polygon->getEnvelope(&bound);
	//	printf("%f,%f\n", bound.MaxX - bound.MinX, bound.MaxY - bound.MinY);
	//	ofs << bound.MaxX - bound.MinX << "," << bound.MaxY - bound.MinY << std::endl;
	//	//double area = ;
	//	//poFeature->SetField(idx, area);
	//	//file->poLayer->SetFeature(poFeature);
	//	//OGRFeature::DestroyFeature(poFeature);
	//}
	//ofs.close();
	//delete file;
	//Grid grid;
	//grid.fromFishnetShape("C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp");
	//grid.toShape(grid.proj,"C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp");
}
void gridIndy_Thomas()
{

 //######################################################################
	std::string rootdir = "B:/Indianapolis/Marion/Vulcan_output/";
	std::string outdir = rootdir + "GriddedEmissions/";
	std::string shapesindir = "B:/Indianapolis/Marion/Vulcan_output/Spatial/";
	std::string intersected = rootdir + "intersected_tmp/";
	std::string indir = rootdir + "intersected/";
	std::string timedir = "B:/Indianapolis/Time/";
	std::string gridShapeFile = "B:/Indianapolis/GridDef/IndyGrid.shp";
	std::string fishnetraster = "B:/Indianapolis/GridDef/IndyGrid.tif";


	QDir(outdir.data()).mkdir(".");
	QDir(shapesindir.data()).mkdir(".");
	QDir(intersected.data()).mkdir(".");
	QDir(indir.data()).mkdir(".");

	Utils::updateFootprintForDir(shapesindir, false);
	Preprocessor::gridFolderByShape(shapesindir, intersected, gridShapeFile);
	std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"ca15","timestruct" }, 7);
	ShapeFile::copyDirDropGeometry(intersected, indir, fields2keep);

	std::string years[]{ "2010","2011","2012","2013","2014","2015" };
	std::string cafields[]{ "ca10","ca11","ca12","ca13","ca14","ca15" };
	std::map<std::string, std::vector<std::string>> timestructfile_crosswalk;

	std::vector<GridSectorConfig> sectors;
	GridSectorConfig com("commercial", Utils::buildVector(new std::string[1]{ "comnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));
	GridSectorConfig res("residential", Utils::buildVector(new std::string[1]{ "resnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));
	GridSectorConfig ind("industrial", Utils::buildVector(new std::string[1]{ "indnonpoint" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));
	GridSectorConfig onroad("onroad", Utils::buildVector(new std::string[1]{ "onroad" }, 1),
		Utils::buildVector(new std::string[1]{ "" }, 1));

	//sectors.push_back(com);
	//sectors.push_back(res);
	//sectors.push_back(ind);
	sectors.push_back(onroad);


	QDir(outdir.data()).mkpath(".");
	std::string cityname = "";
	std::string version = "";
	Grid grid;
	grid.fromFishnetRaster(fishnetraster);
	for (int i = 0; i < 6; i++)
	{

		std::string year = years[i];
		std::string attName = "ca" + year.substr(2, 2);
		//std::string timestructdir = timedir + year + "/";
		std::string timestructdir = timedir;
		//TimestructTool::normalizeBinary(timestructfile);
		SubdirManager dirmanager(indir);
		for (int nsector = 0; nsector < sectors.size(); nsector++)
		{
			grid.reset(attName);
			int numhous = 8760;
			if (year == "2012")
			{
				numhous = 8784;
			}
			TemporalGridder gridder(outdir + year + "/", numhous,"deg");
			gridder.fromFishnetRaster(fishnetraster);

			std::vector<std::string> shapefiles = dirmanager.findFilesMatch(sectors[nsector].shapefiles);
			HestiaGrid* sectorGrid = gridder.addSectorGrid(shapefiles, sectors[nsector].sectorname);
		    sectorGrid->loadtimestruct(timestructdir + year + ".bin");
			gridder.loadAttribute(cafields[i]);

			printf("%s\n", std::string(outdir + sectorGrid->sectorname + "." + year + ".annual" + ".shp").data());
			sectorGrid->getTotal();
			sectorGrid->toShapefile(shapesindir + "ResNonPoint.shp", outdir + sectorGrid->sectorname + "." + year + ".annual" + ".shp");
			grid.reset(attName);
			grid.gatherCells(outdir + sectorGrid->sectorname + "." + year + ".annual" + +".shp", attName.data());
			grid.toRaster(outdir + sectorGrid->sectorname + "." + year + ".annual" + ".tif");
			gridder.makeAnnualTotal(outdir + sectorGrid->sectorname + "." + year + ".annual" + ".nc");
			gridder.makeHourlyTotal(outdir + sectorGrid->sectorname + "." + year + ".hourly" + ".nc");

		}
	}

}
double load_LA_VY_To_Shapefile(std::string gridfile, std::string rootdir, std::string name)
{

	OGRSpatialReference oSRS;
	oSRS.SetWellKnownGeogCS("WGS84");
	ShapeFile shp;
	shp.create(rootdir + name + ".shp", &oSRS, NULL, OGRwkbGeometryType::wkbPolygon);
	std::string line;
	std::ifstream ifs;
	ifs.open(rootdir + "GRID_SCB_0_01.csv");
	std::string headerline;
	std::getline(ifs, headerline);
	//int latfield = shp.getOrCreateField("lat", OGRFieldType::OFTReal);
	//int lonfield = shp.getOrCreateField("lon", OGRFieldType::OFTReal);
	//int idfield = shp.getOrCreateField("Id", OGRFieldType::OFTInteger);
	int field = shp.getOrCreateField(name.data(), OGRFieldType::OFTReal);
	double halfcellsize = 0.01 * 0.5;

	std::ifstream ifsAnnual(gridfile.data(), std::ios::binary);
	ifsAnnual.seekg(0, ifsAnnual.end);
	size_t fileSize = ifsAnnual.tellg();
	ifsAnnual.seekg(0, ifsAnnual.beg);
	size_t numcells = fileSize / 8;
	double* data = new double[numcells];
	ifsAnnual.read((char*)data, fileSize);
	ifsAnnual.close();
	int id = 0;
	double totalAnnual = 0;
	double* pdata = data;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		double lat = atof(splits[1].data());
		double lon = atof(splits[2].data());
		OGRFeature* poFeature = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		OGREnvelope gridcell;
		gridcell.MaxX = lon + halfcellsize;
		gridcell.MaxY = lat + halfcellsize;
		gridcell.MinX = lon - halfcellsize;
		gridcell.MinY = lat - halfcellsize;
		OGRGeometry* cellpoly = createBox(gridcell);
		poFeature->SetGeometry(cellpoly);
		//poFeature->SetField(latfield, lat);
		//poFeature->SetField(lonfield, lon);
		poFeature->SetField(field, *pdata);
		shp.poLayer->CreateFeature(poFeature);
		//OGRGeometryFactory::destroyGeometry(cellpoly);
		OGRFeature::DestroyFeature(poFeature);
		totalAnnual += *pdata++;
		id++;
	}
	ifs.close();
	shp.close();
	delete[] data;
	
	//ShapeFile shp(shapefile, 1);

	
	
	/*double totalAnnual = 0;*/
	/*for (size_t i = 0; i < numcells; i++)
	{
		fea = shp.poLayer->GetNextFeature();
		fea->SetField(field, data[i]);
		OGRFeature::DestroyFeature(fea);
		totalAnnual += data[i];
	}*/
	printf("total = %f\n", totalAnnual/1000000/1000);

	//shp.close();
	return totalAnnual;
}

void gridLA_VY_Totals(std::string year)
{

	/*ShapeFile shp("B:/LA_Version2/Vulcan_output/Vineet/Results", 1);
	int caIdx = shp.getOrCreateField("ca", OGRFieldType::OFTReal);
	OGRFeature *poFeature;
	
	shp.poLayer->ResetReading();
	int idx = -1;*/
	int numcells = 17181;
	//int numhours = 8784;
	double totalAnnual = 0;
	double totalHourly = 0;


	std::vector<std::string> files = Utils::findFiles("B:/LA_Version2/Vulcan_output/Vineet/GriddedEmissions/", year + ".bin");
	std::ofstream ofstotals;
	ofstotals.open("B:/LA_Version2/Vulcan_output/Vineet/GriddedEmissions/" + year + "_totals.csv");
	ofstotals << "File,Total(MtC)" << std::endl;
	for (size_t nfile = 0; nfile < files.size(); nfile++)
	{
		std::ifstream ifsAnnual(files[nfile].data(), std::ios::binary);
		ifsAnnual.seekg(0, ifsAnnual.end);
		size_t fileSize = ifsAnnual.tellg();
		ifsAnnual.seekg(0, ifsAnnual.beg);
		size_t numcells = fileSize / 8;
		double* data = new double[numcells];
		ifsAnnual.read((char*)data, fileSize);
		ifsAnnual.close();
		totalAnnual = 0;
		for (size_t i = 0; i < numcells; i++)
		{
			totalAnnual += data[i];
		}
		ofstotals << std::fixed << std::setprecision(2) << QFileInfo(files[nfile].data()).fileName().toLocal8Bit().data() << "," << totalAnnual/1000000.0/1000.0 << std::endl;
		delete[] data;
	}
	ofstotals.close();

	//std::ifstream ifsAnnual("C:/HestiaGridding/Los_Angeles/Vineet/output/WRF.VY.onroad.hourly.2014.bin", std::ios::binary);
	//ifsAnnual.read((char*)outlineDS, sizeof(double) * numcells * 8760);
	//ifsAnnual.close();
	////	memset(outlineDS, 0, sizeof(double) * gridder.numcells * 8760);
	////	std::ifstream ifsHourly(sectorHourlyFile.outlineDS(), std::ios::binary);
	////	ifsHourly.read((char*)outlineDS, sizeof(double) * gridder.numcells * 8760);
	////	ifsHourly.close();
	////	for (size_t j = 0; j < gridder.numcells * 8760; j++)
	////	{
	////		totalHourly += outlineDS[j];
	////	}
	//double* pdata = outlineDS + (numcells * 100);
	//std::ofstream ofs("C:/HestiaGridding/Los_Angeles/Vineet/onroadtime.csv", std::ios::binary);
	//for (size_t i = 0; i < 8760; i++)
	//{
	//	pdata = outlineDS + (numcells * i);
	//	for (size_t icell = 0; icell < 100; icell++)
	//	{
	//		ofs << pdata[icell + 500];
	//		if (icell < 100 - 1)
	//			ofs << ",";
	//		else
	//			ofs << "\n";

	//	}
	//}

	//ofs.close();



	//while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
	//{
	//	idx++;
	//	poFeature->SetField(caIdx, *pdata++);
	//	Shp.poLayer->SetFeature(poFeature);
	//	OGRFeature::DestroyFeature(poFeature);
	//}
	//Shp.close();
	////��WRF.VY.*****�� to indicate that this is the WRF grid for Vineet Yadav
	//std::string sourceDir = "B:/LA_Version2/gridPrep_SHP_master/";
	//std::string caDir = "C:/HestiaGridding/Los_Angeles/Vineet/ca/";
	//std::string indir = "C:/HestiaGridding/Los_Angeles/Vineet/intersected/";
	//std::string intersectedDirNoGeometry = "C:/HestiaGridding/Los_Angeles/Vineet/intersectedNoGeometry/";
	//std::string lidarFisheyeDir = "C:/HestiaGridding/Los_Angeles/Vineet/output/";
	//QDir(caDir.outlineDS()).mkpath(".");
	//QDir(indir.outlineDS()).mkpath(".");
	//QDir(intersectedDirNoGeometry.outlineDS()).mkpath(".");
	//QDir(lidarFisheyeDir.outlineDS()).mkpath(".");
	//std::string gridShapeFile = "C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp";
	//std::string projectName = "WRF.VY";
	//std::string year = "2014";
	//std::string cafield = "ca" + year.substr(2, 2);
	//std::string timestructfile = "C:/HestiaGridding/Los_Angeles/timestructs_" + year + ".bin";
	//std::vector<std::string> subdirs = Utils::findSubdirectories(sourceDir);
	//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[9]{ "ca10","ca11","ca12","ca13","ca14","length","area" ,"timestruct","Id" }, 9);
	//for (size_t i = 0; i < subdirs.size(); i++)
	//{
	//	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
	//	//Utils::updateFootprint(destFile);
	//	// Utils::updateFootprintForDir(subdirs[i], false);
	//	//ShapeFile::copyDir(subdirs[i], caDir + foldername + "/", fields2keep);
	//	//Preprocessor::gridFolderByShape(caDir + foldername + "/", indir + foldername + "/", gridShapeFile);
	//	//ShapeFile::copyDirDropGeometry(indir + foldername + "/", intersectedDirNoGeometry + foldername + "/", fields2keep);
	//}

	//SubdirManager dirmanager = SubdirManager(intersectedDirNoGeometry);
	//std::vector<std::string> comshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "comnonpoint","compoint" }, 2));
	//std::vector<std::string> indshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "indnonpoint","indpoint" }, 2));
	//std::vector<std::string> resshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "resnonpoint" }, 1));
	//std::vector<std::string> railroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "createRailroadPoint","Railroad" }, 2));
	//std::vector<std::string> onroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "OnRoad.Shp" }, 1));
	//std::vector<std::string> nonroadshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "NonRoad" }, 1));
	//std::vector<std::string> elecprodshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "ElecProd" }, 1));
	//std::vector<std::string> marineshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[2]{ "Port_Lanes","Port_Polygons" }, 2));
	//std::vector<std::string> airportshapefiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "airport" }, 1));
	//TemporalGridderByShapes gridder;
	//gridder.numcells = 1826;
	//gridder.cells = new double[gridder.numcells];
	//memset(gridder.cells, 0, gridder.numcells*sizeof(double));
	////gridder.loadtimestruct("", timestructfile);
	//gridder.addSectorGrid(comshapefiles, "com");
	//gridder.addSectorGrid(indshapefiles, "ind");
	//gridder.addSectorGrid(resshapefiles, "res");
	//gridder.addSectorGrid(railroadshapefiles, "railroad");
	//gridder.addSectorGrid(onroadshapefiles, "onroad");
	//gridder.addSectorGrid(nonroadshapefiles, "nonroad");
	//gridder.addSectorGrid(elecprodshapefiles, "elecprod");
	//gridder.addSectorGrid(marineshapefiles, "marine");
	//gridder.addSectorGrid(airportshapefiles, "airport");
	//gridder.loadAttribute(cafield);
	//std::string annualFile = lidarFisheyeDir + projectName + ".total.annual." + year + ".bin";
	//std::string hourlyFile = lidarFisheyeDir + projectName + ".total.hourly." + year + ".bin";
	//for (size_t isector = 0; isector < gridder.sectors.size(); isector++)
	//{
	//	FFCO2Sector* outputGrid = gridder.sectors[isector];
	//	std::string sectorAnnualFile = lidarFisheyeDir + projectName + "." + outputGrid->sectorname + ".annual." + year + ".bin";
	//	std::string sectorHourlyFile = lidarFisheyeDir + projectName + "." + outputGrid->sectorname + ".hourly." + year + ".bin";
	//	double total = 0;
	//	for (size_t j = 0; j < outputGrid->features.size(); j++)
	//	{
	//		total += outputGrid->features[j].total;
	//	}
	//	double totalAnnual = 0;
	//	double totalHourly = 0;
	//	double* outlineDS = new double[gridder.numcells * 8760];
	//	memset(outlineDS, 0, sizeof(double) * gridder.numcells * 8760);
	//	std::ifstream ifsAnnual(sectorAnnualFile.outlineDS(), std::ios::binary);
	//	ifsAnnual.read((char*)outlineDS, sizeof(double) * gridder.numcells);
	//	ifsAnnual.close();
	//	for (size_t j = 0; j < gridder.numcells; j++)
	//	{
	//		totalAnnual += outlineDS[j];
	//	}
	//	memset(outlineDS, 0, sizeof(double) * gridder.numcells * 8760);
	//	std::ifstream ifsHourly(sectorHourlyFile.outlineDS(), std::ios::binary);
	//	ifsHourly.read((char*)outlineDS, sizeof(double) * gridder.numcells * 8760);
	//	ifsHourly.close();
	//	for (size_t j = 0; j < gridder.numcells * 8760; j++)
	//	{
	//		totalHourly += outlineDS[j];
	//	}

	//	printf("%s:total=%f,totalAnnual=%f,totalHourly=%f\n", outputGrid->sectorname.outlineDS(),total, totalAnnual, totalHourly);
	//}


	/*gridder.makeAnnualTotal(annualFile);

	for (size_t isector = 0; isector < gridder.sectors.size(); isector++)
	{
	FFCO2Sector* outputGrid = gridder.sectors[isector];
	std::string sectorAnnualFile = lidarFisheyeDir + projectName + "." + outputGrid->sectorname + ".annual." + year + ".csv";
	std::string sectorHourlyFile = lidarFisheyeDir + projectName + "." + outputGrid->sectorname + ".hourly." + year + ".csv";
	outputGrid->makeAnnualTotal(sectorAnnualFile);
	outputGrid->makeHourlyTotal(sectorHourlyFile);
	}

	gridder.makeHourlyTotal(hourlyFile);*/
	//Utils::updateFootprint("C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp", 110959.0* 110959.0);
	//ShapeFile* file = new ShapeFile("C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp", 1);
	//OGRFeature *poFeature;
	//std::ofstream ofs;
	//ofs.open("VYgrids.csv");
	//while ((poFeature = file->poLayer->GetNextFeature()) != NULL)
	//{
	//	OGRPolygon* polygon = (OGRPolygon*)poFeature->GetGeometryRef();
	//	OGREnvelope bound;
	//	polygon->getEnvelope(&bound);
	//	printf("%f,%f\n", bound.MaxX - bound.MinX, bound.MaxY - bound.MinY);
	//	ofs << bound.MaxX - bound.MinX << "," << bound.MaxY - bound.MinY << std::endl;
	//	//double area = ;
	//	//poFeature->SetField(idx, area);
	//	//file->poLayer->SetFeature(poFeature);
	//	//OGRFeature::DestroyFeature(poFeature);
	//}
	//ofs.close();
	//delete file;
	//Grid grid;
	//grid.fromFishnetShape("C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp");
	//grid.toShape(grid.proj,"C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp");
}
void createBox(OGREnvelope bb, OGRLayer* layer)
{
	OGRPolygon *poPolygon = createBox(bb);
	OGRFeature *poFeatureNew = OGRFeature::CreateFeature(layer->GetLayerDefn());
	poFeatureNew->SetGeometry(poPolygon);
	layer->CreateFeature(poFeatureNew);
	OGRFeature::DestroyFeature(poFeatureNew);
}
#include "MultitheadIntersection.h"
void intersect()
{

	//ShapeFile Shp(fullname.outlineDS());

	//if (master == NULL)
	//{
	//	master = new ShapeFile;
	//	master->create(outfile, Shp.poLayer->GetSpatialRef(), NULL, Shp.poLayer->GetGeomType());
	//	OGRFieldDefn newField("SEGID", OGRFieldType::OFTString);
	//	destSEGID = master->poLayer->CreateField(&newField);
	//	master->poLayer->GetLayerDefn()->GetFieldDefn(destSEGID)->SetWidth(12);
	//}
	//srcSEGID = Shp.poLayer->GetLayerDefn()->GetFieldIndex("SEGID");
	//Shp.poLayer->ResetReading();
	//OGRFeature *poFeature;
	//int id = 0;
	//while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
	//{
	//	OGRFeature *poFeatureNew = OGRFeature::CreateFeature(master->poLayer->GetLayerDefn());
	//	poFeatureNew->SetGeometry(poFeature->GetGeometryRef());
	//	poFeatureNew->SetField(destSEGID, poFeature->GetFieldAsString(srcSEGID));
	//	master->poLayer->CreateFeature(poFeatureNew);
	//	OGRFeature::DestroyFeature(poFeatureNew);
	//	OGRFeature::DestroyFeature(poFeature);
	//}
	//Shp.close();

	//ShapeFile ringshp;
	//ringshp.create("ring.Shp", 0, 0, wkbPolygon);
	//OGRPolygon* rings = createCircle(36,100);

	//OGRFeature *poFeatureNew = OGRFeature::CreateFeature(ringshp.poLayer->GetLayerDefn());
	//poFeatureNew->SetGeometry(rings);
	//ringshp.poLayer->CreateFeature(poFeatureNew);
	//OGRFeature::DestroyFeature(poFeatureNew);
	//ringshp.close();

	//ShapeFile ringshp;
	//ringshp.create("ring.Shp", 0, 0, wkbPolygon);
	//createCircle(36, 100, ringshp.poLayer);
	//ringshp.close();

	//ShapeFile boxshp;
	//boxshp.create("box.Shp", 0, 0, wkbPolygon);
	//OGREnvelope bb;
	//bb.MinX = -120;
	//bb.MaxX = 120;
	//bb.MinY = -80;
	//bb.MaxY = 80;
	//createBox(bb, boxshp.poLayer);

	//boxshp.close();

	//OGRPolygon* p1 = createCircle(36, 100);
	//OGRPolygon* p2 = createBox(bb);
	//OGRPolygon* pIntersection = (OGRPolygon*)p1->Intersection(p2);
	//ShapeFile intersectshp;
	//intersectshp.create("intersection.Shp", 0, 0, wkbPolygon);
	//OGRFeature *poFeatureNew = OGRFeature::CreateFeature(intersectshp.poLayer->GetLayerDefn());
	//poFeatureNew->SetGeometry(pIntersection);
	//intersectshp.poLayer->CreateFeature(poFeatureNew);
	//OGRFeature::DestroyFeature(poFeatureNew);
	//intersectshp.close();
	//OGREnvelope bb;
	//bb.MinX = -100;
	//bb.MaxX = 100;
	//bb.MinY = -100;
	//bb.MaxY = 100;
	//Grid* grid = new Grid(bb, 10);
	//grid->reset();
	//grid->intersect("ring.Shp", "intersection.Shp");

	std::string inshapefile = "C:/VulcanGridding/gridPrep_SHP_master/ComNonPoint.shp";
	//OGREnvelope gridbb = BoundManager::readBoundFromShape(inshapefile);
	std::string outshapefile = "C:/VulcanGridding/Intersected/ComNonPoint.shp";
	MultitheadIntersection gridder;
	//gridder.setupGrid(gridbb, 200 * 3.28084);
	gridder.setupGrid("C:/VulcanGridding/fishnet.tif");
	gridder.intersectParallel(inshapefile, outshapefile, 8);

}

void linkParcels2UtilityArea(std::string parcelfile, std::string utilityAreaFile)
{
	std::pair<std::string, float> emissionRatios[9];
	emissionRatios[0] = std::pair<std::string, float>("Los Angeles Department of Water & Power", 642.911);
	emissionRatios[1] = std::pair<std::string, float>("Southern California Edison", 595.099);
	emissionRatios[2] = std::pair<std::string, float>("Azusa Light & Power", 792.78);
	emissionRatios[3] = std::pair<std::string, float>("City of Industry", 792.78);
	emissionRatios[4] = std::pair<std::string, float>("City of Cerritos", 792.78);
	emissionRatios[5] = std::pair<std::string, float>("Glendale Water & Power", 958.106);
	emissionRatios[6] = std::pair<std::string, float>("City of Vernon Municipal Light Department", 807.89);
	emissionRatios[7] = std::pair<std::string, float>("Burbank Water & Power", 792.78);
	emissionRatios[8] = std::pair<std::string, float>("Pasadena Water & Power", 1063.44);

	GDAL_DS<unsigned char>* utilityArea = new GDAL_DS<unsigned char>();
	utilityArea->open(utilityAreaFile);
	unsigned char* utilityAreaGrid = utilityArea->readData(1);
	ShapeFile shp(parcelfile, 1);
	int utilityAreaNameIdx = shp.getOrCreateField("Utility", OGRFieldType::OFTString);
	int ER11Idx = shp.getOrCreateField("ER11", OGRFieldType::OFTReal);

	OGRFeature* poFeature;
	double  minX = utilityArea->adfGeoTransform[0];
	double  maxY = utilityArea->adfGeoTransform[3];
	double cellsize = utilityArea->adfGeoTransform[1];
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGREnvelope bound;
		poFeature->GetGeometryRef()->getEnvelope(&bound);
		double cx = (bound.MinX + bound.MaxX) * 0.5;
		double cy = (bound.MinY + bound.MaxY) * 0.5;
		int row = (int)((maxY - cy) / cellsize);
		int col = (int)((cx - minX) / cellsize);
		if (row < 0)
			row = 0;
		if (row > utilityArea->nrows - 1)
			row = utilityArea->nrows - 1;
		if (col < 0)
			col = 0;
		if (col > utilityArea->ncols - 1)
			col = utilityArea->ncols - 1;
		unsigned char utilityAreaCode = utilityAreaGrid[col + row * utilityArea->ncols];
		poFeature->SetField(utilityAreaNameIdx, emissionRatios[utilityAreaCode].first.data());
		poFeature->SetField(ER11Idx, emissionRatios[utilityAreaCode].second);
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();

}
void CalScope2Utility(std::string parcelfile)
{

	ShapeFile shp(parcelfile, 1);
	int ER11Idx = shp.getOrCreateField("ER11", OGRFieldType::OFTReal);
	int kwhIdx = shp.getOrCreateField("kwh", OGRFieldType::OFTReal);
	int caIdx = shp.getOrCreateField("ca_ER11", OGRFieldType::OFTReal);
	double pound2kg = 0.453592;
	double co2c = 12.01 / 44.01;
	double conversion = pound2kg * co2c;
	OGRFeature* poFeature;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		double ER11 = poFeature->GetFieldAsDouble(ER11Idx);
		double kwh = poFeature->GetFieldAsDouble(kwhIdx);
		double ca = kwh * ER11 / 1000.0 * conversion;
		poFeature->SetField(caIdx, ca);
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();

}

void parcels2neighborhood()
{

	std::string dir = "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/";
	//CalScope2Utility(dir + "ComNonPoint.Shp");
	//CalScope2Utility(dir + "ResNonPoint.Shp");
	//CalScope2Utility(dir + "IndNonPoint.Shp");

	//reallocateLA();

	std::vector<OGRFieldData> matchFields;

	matchFields.push_back(OGRFieldData("yearbuilt", "yearbuilt", OperatorType::AVERAGE));
	//srcnames.push_back("area");
	//srcnames.push_back("ca_ng");
	//srcnames.push_back("ca_p");
	////srcnames.push_back("elec");
	//srcnames.push_back("ca10_ng");
	//srcnames.push_back("ca10_p");
	//srcnames.push_back("ca11_ng");
	//srcnames.push_back("ca11_p");

	//srcnames.push_back("kwh");
	//srcnames.push_back("ca_ER11");

	//Geographies geographies;
	//std::string destFile = dir + "LA_Neighborhoods/LA_Neighborhoods.Shp";
	//////Preprocessor::gridFolderByShape(dir , dir + "LA_Neighborhoods/", destFile);
	//geographies.setup(dir + "LA_Neighborhoods/ComNonPoint.Shp", matchFields, dir + "LA_Neighborhoods/ComNeighborhoods.Shp");
	//geographies.gather(dir + "LA_Neighborhoods/ComNonPoint.Shp", "FID_LA_Nei");
	//geographies.update();
	//
	//geographies.setup(dir + "LA_Neighborhoods/ResNonPoint.Shp", matchFields, dir + "LA_Neighborhoods/ResNeighborhoods.Shp");
	//geographies.gather(dir + "LA_Neighborhoods/ResNonPoint.Shp", "FID_LA_Nei");
	//geographies.update();

	//geographies.setup(dir + "LA_Neighborhoods/IndNonPoint.Shp", srcnames, dir + "LA_Neighborhoods/IndNeighborhoods.Shp");
	//geographies.gather(dir + "LA_Neighborhoods/IndNonPoint.Shp", "FID_LA_Nei");
	//geographies.update();
}

void blockgroups2neighborhood()
{

	std::string dir = "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/";
	//CalScope2Utility(dir + "ComNonPoint.Shp");
	//CalScope2Utility(dir + "ResNonPoint.Shp");
	//CalScope2Utility(dir + "IndNonPoint.Shp");

	//reallocateLA();

	std::vector<OGRFieldData> matchFields;

	matchFields.push_back(OGRFieldData("yearbuilt", "yearbuilt", OperatorType::AVERAGE));
	//srcnames.push_back("area");
	//srcnames.push_back("ca_ng");
	//srcnames.push_back("ca_p");
	////srcnames.push_back("elec");
	//srcnames.push_back("ca10_ng");
	//srcnames.push_back("ca10_p");
	//srcnames.push_back("ca11_ng");
	//srcnames.push_back("ca11_p");

	//srcnames.push_back("kwh");
	//srcnames.push_back("ca_ER11");

	//Geographies geographies;
	//std::string destFile = dir + "LA_Neighborhoods/LA_Neighborhoods.Shp";
	//////Preprocessor::gridFolderByShape(dir , dir + "LA_Neighborhoods/", destFile);
	//geographies.setup(dir + "LA_Neighborhoods/ComNonPoint.Shp", matchFields, dir + "LA_Neighborhoods/ComNeighborhoods.Shp");
	//geographies.gather(dir + "LA_Neighborhoods/ComNonPoint.Shp", "FID_LA_Nei");
	//geographies.update();
	//
	//geographies.setup(dir + "LA_Neighborhoods/ResNonPoint.Shp", matchFields, dir + "LA_Neighborhoods/ResNeighborhoods.Shp");
	//geographies.gather(dir + "LA_Neighborhoods/ResNonPoint.Shp", "FID_LA_Nei");
	//geographies.update();

	//geographies.setup(dir + "LA_Neighborhoods/IndNonPoint.Shp", srcnames, dir + "LA_Neighborhoods/IndNeighborhoods.Shp");
	//geographies.gather(dir + "LA_Neighborhoods/IndNonPoint.Shp", "FID_LA_Nei");
	//geographies.update();
}
void linkBlockgroups2Neighborhoods(std::string blockgroupfile, std::string neighborhoodFile)
{

	GDAL_DS<unsigned char>* neighborhoods = new GDAL_DS<unsigned char>();
	neighborhoods->open(neighborhoodFile);
	unsigned char* utilityAreaGrid = neighborhoods->readData(1);
	ShapeFile shp(blockgroupfile, 1);
	int NeihoodIDIdx = shp.getOrCreateField("NeihoodID", OGRFieldType::OFTInteger);

	OGRFeature* poFeature;
	double  minX = neighborhoods->adfGeoTransform[0];
	double  maxY = neighborhoods->adfGeoTransform[3];
	double cellsize = neighborhoods->adfGeoTransform[1];
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGREnvelope bound;
		poFeature->GetGeometryRef()->getEnvelope(&bound);
		double cx = (bound.MinX + bound.MaxX) * 0.5;
		double cy = (bound.MinY + bound.MaxY) * 0.5;
		int row = (int)((maxY - cy) / cellsize);
		int col = (int)((cx - minX) / cellsize);
		if (row < 0)
			row = 0;
		if (row > neighborhoods->nrows - 1)
			row = neighborhoods->nrows - 1;
		if (col < 0)
			col = 0;
		if (col > neighborhoods->ncols - 1)
			col = neighborhoods->ncols - 1;
		unsigned char neihoodCode = utilityAreaGrid[col + row * neighborhoods->ncols];
		poFeature->SetField(NeihoodIDIdx, neihoodCode);
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();

}

void updateID()
{
	ShapeFile inshp("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/LA_Neighborhoods/ResNonPoint.shp");
	ShapeFile outshp("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint.shp", 1);
	OGRFeature* fea;
	int srcfield = inshp.poLayer->GetLayerDefn()->GetFieldIndex("FID_ResNon");
	int destfield = inshp.poLayer->GetLayerDefn()->GetFieldIndex("FID_LA_Nei");
	std::map<int, double> areaMap;
	int numdestfeas = outshp.poLayer->GetFeatureCount();

	int* idmap = new int[numdestfeas];
	printf("%d\n", numdestfeas);
	printf("%d\n", srcfield);
	printf("%d\n", destfield);
	for (size_t i = 0; i < numdestfeas; i++)
	{
		idmap[i] = -1;
	}
	inshp.poLayer->ResetReading();
	while ((fea = inshp.poLayer->GetNextFeature()) != NULL)
	{
		double area = Utils::calPolygonArea(fea->GetGeometryRef());
		int srcid = fea->GetFieldAsInteger(srcfield);
		int destid = fea->GetFieldAsInteger(destfield);
		std::map<int, double>::iterator iter = areaMap.find(srcid);
		if (iter == areaMap.end())
		{
			areaMap[srcid] = area;
			idmap[srcid] = destid;
		}
		else
		{
			if (iter->second > area)
			{
				areaMap[srcid] = area;
				idmap[srcid] = destid;
			}
		}
		OGRFeature::DestroyFeature(fea);
	}
	inshp.close();
	destfield = outshp.poLayer->GetLayerDefn()->GetFieldIndex("nei_id");
	outshp.poLayer->ResetReading();
	int curid = 0;
	while ((fea = outshp.poLayer->GetNextFeature()) != NULL)
	{
		fea->SetField(destfield, idmap[curid]);
		curid++;
		outshp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}
	outshp.close();

}

void MODIS_LST_2_NEIGHBORHOODS(std::string resnonpointfile, std::string neighborhoodfile, std::string lstfile)
{
	GDAL_DS<unsigned short>* resnonpoint_ds = new GDAL_DS<unsigned short>();
	resnonpoint_ds->open(resnonpointfile);
	unsigned short* resnonpoint_buf = resnonpoint_ds->readData(1);
	unsigned short resnonpoint_nodata = (unsigned short)resnonpoint_ds->getNoData(1);

	GDAL_DS<unsigned short>* neighborhood_ds = new GDAL_DS<unsigned short>();
	neighborhood_ds->open(neighborhoodfile);
	unsigned short* neighborhood_buf = neighborhood_ds->readData(1);
	unsigned short neighborhood_nodata = (unsigned short)neighborhood_ds->getNoData(1);

	//resnonpoint_ds->crop(neighborhood_ds->bound, "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/Neighborhood2.tif");

	GDAL_DS<unsigned short>* lst_ds = new GDAL_DS<unsigned short>();
	lst_ds->open(lstfile);
	unsigned short* lst_buf = lst_ds->readData(2);
	unsigned short lst_nodata = (unsigned short)lst_ds->getNoData(1);

	OGRSpatialReference oTargetSRS;
	char* destprojbuf = new char[lst_ds->projection.size()];
	memcpy(destprojbuf, lst_ds->projection.data(), lst_ds->projection.size());
	oTargetSRS.importFromWkt(&destprojbuf);

	OGRSpatialReference oSourceSRS;
	char* srcprojbuf = new char[neighborhood_ds->projection.size()];
	memcpy(srcprojbuf, neighborhood_ds->projection.data(), neighborhood_ds->projection.size());
	oSourceSRS.importFromWkt(&srcprojbuf);

	OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
	double x = -3.17970006;
	double y = 51.4789137;

	double destx = x;
	double desty = y;
	poCT->Transform(1, &destx, &desty);
	double* lstSum = new double[272];
	double* lstCount = new double[272];
	for (size_t i = 0; i < 272; i++)
	{
		lstSum[i] = 0;
		lstCount[i] = 0;
	}
	int srcidx = 0;
	for (int irow = 0; irow < neighborhood_ds->nrows; irow++)
	{
		double y = neighborhood_ds->bound.MaxY + irow * neighborhood_ds->adfGeoTransform[5];
		for (int icol = 0; icol < neighborhood_ds->ncols; icol++)
		{
			double x = neighborhood_ds->bound.MinX + icol * neighborhood_ds->adfGeoTransform[1];
			unsigned int neiidx = neighborhood_buf[srcidx];
			unsigned int residx = resnonpoint_buf[srcidx];
			srcidx++;
			if (neiidx == neighborhood_nodata)
			{
				//printf("%d\n", neiidx);
				continue;
			}
			if (residx == resnonpoint_nodata)
			{
				//printf("%d\n", residx);
				continue;
			}
			double destx = x;
			double desty = y;
			poCT->Transform(1, &destx, &desty);
			int destCol = (int)((destx - lst_ds->bound.MinX) / lst_ds->adfGeoTransform[1]);
			int destRow = (int)((desty - lst_ds->bound.MaxY) / lst_ds->adfGeoTransform[5]);

			if (destx < lst_ds->bound.MinX || destx > lst_ds->bound.MaxX || desty < lst_ds->bound.MinY || desty > lst_ds->bound.MaxY)
				continue;

			if (destCol < 0 || destCol > lst_ds->ncols - 1 || destRow < 0 || destRow > lst_ds->nrows - 1)
				continue;
			unsigned int lst = lst_buf[destCol + destRow*lst_ds->ncols];
			if (lst < 10000 || lst > 20000)
				continue;
			lstCount[neiidx] = lstCount[neiidx] + 1;
			lstSum[neiidx] += lst;
		}
	}


	std::ofstream ofs;
	ofs.open((lstfile + "_NIGHT.csv").data());

	ofs << "NeighborhoodID,LST_DAY" << std::endl;

	for (size_t i = 0; i < 272; i++)
	{
		ofs << i << "," << (lstSum[i] / lstCount[i]) * 0.02 - 273.15 << std::endl;
	}

	ofs.close();

}

#include "NeiborhoodOptimizationMatrixPrep.h"
void TwinCitiesComparison()
{

	GDAL_DS<double>* dsMax = new GDAL_DS<double>();
	dsMax->open("E:/TwinCitiesComparison/max.tif");
	GDAL_DS<double>* dsMin = new GDAL_DS<double>();
	dsMin->open("E:/TwinCitiesComparison/min.tif");

	OGRSpatialReference oTargetSRS;
	char* buf = new char[dsMax->projection.size()];
	memcpy(buf, dsMax->projection.data(), dsMax->projection.size());
	oTargetSRS.importFromWkt(&buf);
	OGRSpatialReference oSourceSRS;
	oSourceSRS.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
	double x = -3.17970006;
	double y = 51.4789137;

	double destx = x;
	double desty = y;
	poCT->Transform(1, &destx, &desty);
	printf("(%f,%f) -> (%f,%f)\n", x, y, destx, desty);
	SVFComputer<double>* svfcomMax = new SVFComputer<double>(dsMax);
	SVFComputer<double>* svfcomMin = new SVFComputer<double>(dsMin);
	std::ifstream ifs;
	ifs.open("E:/TwinCitiesComparison/GSV_Data/points.csv");
	std::string line;
	std::getline(ifs, line);
	std::string lidarFisheyeDir = "E:/TwinCitiesComparison/GSV_Data/LiDARFisheye/";
	std::string googleFisheyeDir = "E:/TwinCitiesComparison/GSV_Data/SegnetFisheye/";
	std::ofstream ssSVF("E:/TwinCitiesComparison/GSV_Data/svf3.csv");
	ssSVF << "ID,Lon,Lat,X,Y,GroundElev,Greenery,LiDAR,Google" << std::endl;
	SVFComputer<unsigned char>* googlesvfcom = new SVFComputer<unsigned char>();

	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		std::string panoid = splits[0];
		if (atoi(panoid.data()) != 2815)
		{
			continue;
		}
		double x = atof(splits[1].data());
		double y = atof(splits[2].data());
		double lon = x;
		double lat = y;
		poCT->Transform(1, &x, &y);

		//x = 480913.715;
		//y = 4984213.878;
		//if (QDir((dir + id).outlineDS()).exists())
		//{
		std::stringstream ssIdx;
		ssIdx << panoid << ".png";
		std::string lidarfile = lidarFisheyeDir + ssIdx.str();
		std::string googlefile = googleFisheyeDir + ssIdx.str();
		double groundelev = svfcomMin->getHeightAt(x, y);
		double eyeheight = groundelev + 1.4 + 1;
		printf("%f\n", eyeheight);
		GDAL_DS<unsigned char>* dsImage = new GDAL_DS<unsigned char>();
		googlesvfcom->skyColor = osg::Vec4i(128, 128, 128, 255);
		googlesvfcom->greenColor = osg::Vec4i(128, 128, 192, 255);
		googlesvfcom->outsideColor = osg::Vec4i(255, 255, 255, 0);
		//-93.2537	44.9421

		double googleSVF = -1;
		double greenery = -1;
		if (dsImage->open(googlefile))
		{
			googlesvfcom->setData(dsImage);
			googleSVF = googlesvfcom->calSVF();
			greenery = googlesvfcom->calGreenery();
		}

		googlesvfcom->skyColor = osg::Vec4i(0, 255, 0, 128);
		googlesvfcom->outsideColor = osg::Vec4i(0, 0, 0, 0);
		double lidarSVF = -1;
		bool isInside = dsMax->isInside(x, y);

		//if (!QFileInfo(lidarfile.outlineDS()).exists() && isInside)
		//{
		osg::Vec3d eye(x, y, eyeheight);
		std::vector<osg::Vec2d> skymapcoords = svfcomMax->computeSkymap(eye, 360);
		svfcomMax->drawSkymap(skymapcoords, 512, lidarfile);
		//}

		if (isInside && dsImage->open(lidarfile))
		{
			googlesvfcom->setData(dsImage);
			lidarSVF = googlesvfcom->calSVF();
		}

		printf("%s,%f,%f\n", panoid.data(), googleSVF, lidarSVF);

		delete dsImage;
		ssSVF << panoid << "," << lon << "," << lat << "," << x << "," << y << "," << groundelev << "," << greenery << "," << lidarSVF << "," << googleSVF << std::endl;
		//idx++;
	}

}
void dissolvePolyline(const OGRLineString  *polyline, double& sampleDist, std::vector<OGRPoint>& points)
{

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
			OGRPoint newp(px, py);
			points.push_back(newp);
		}

	}

}
void polyline2points(double LINE_SAMPLE_DIST, std::string infilename, std::string outfilename)
{

	ShapeFile input(infilename.data());
	ShapeFile output;
	output.create(outfilename, input.poLayer->GetSpatialRef(), 0, wkbPoint);
	int srcFieldIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
	if (srcFieldIndex < 0)
		return;
	int destIDFieldIndex = output.getOrCreateField("rfid", OGRFieldType::OFTInteger);
	int destFieldIndex = output.getOrCreateField("ca11", OGRFieldType::OFTReal);

	OGRFeature* poFeature;
	int fid = -1;
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{
		fid++;

		OGRGeometry* poGeometry = poFeature->GetGeometryRef();
		if (!poGeometry)
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		double value = poFeature->GetFieldAsDouble(srcFieldIndex);
		std::vector<OGRPoint> splitPoints;

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
		double val = poFeature->GetFieldAsDouble(srcFieldIndex);
		for (int igeom = 0; igeom < geomlist.size(); igeom++)
		{
			OGRGeometry* poGeometry = geomlist[igeom];
			OGREnvelope geobb;

			if (dynamic_cast<OGRLineString*>(poGeometry))
			{
				OGRLineString* linestr = (OGRLineString *)poGeometry;
				dissolvePolyline(linestr, LINE_SAMPLE_DIST, splitPoints);
			}
			else if (dynamic_cast<OGRCircularString*>(poGeometry))
			{
				OGRCircularString* circularstr = (OGRCircularString *)poGeometry;
				OGRLineString* linestr = dynamic_cast<OGRLineString*>(circularstr->getLinearGeometry());
				dissolvePolyline(linestr, LINE_SAMPLE_DIST, splitPoints);
			}

		}
		double allocatedCA = value / splitPoints.size();
		for (size_t i = 0; i < splitPoints.size(); i++)
		{
			OGRFeature* newfea = OGRFeature::CreateFeature(output.poLayer->GetLayerDefn());
			newfea->SetGeometry(&splitPoints[i]);
			newfea->SetField(destIDFieldIndex, fid);
			newfea->SetField(destFieldIndex, allocatedCA);

			output.poLayer->CreateFeature(newfea);
			OGRFeature::DestroyFeature(newfea);
		}
		OGRFeature::DestroyFeature(poFeature);
		if (fid % 100 == 0)
			printf("%d\n", fid);
	}

}



class UrbanOutline
{
public:
	int destIDFieldIndex;
	int destFieldIndex;
	std::vector<std::string> dict;
	ShapeFile output;
	GDAL_DS<unsigned int>* maskds;
	unsigned int* mask;
	unsigned int nodata;
	UrbanOutline(std::string rasterfilename, std::string outputfile, ShapeFile* refshape)
	{
		/*std::ifstream ifs;
		ifs.open(rasterfilename + ".txt");
		std::string line;
		std::getline(ifs, line);
		while (ifs.peek() != -1)
		{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		dict.push_back(splits[splits.size() - 1]);

		}*/
		output.create(outputfile, refshape->poLayer->GetSpatialRef(), 0, wkbPoint);
		QFileInfo info(rasterfilename.data());
		std::string name = info.fileName().toLocal8Bit().data();
		int numchars = 6;
		if (name.size() < numchars)
			numchars = name.size();
		std::string destIDFieldName = "FID_" + name.substr(0, numchars);
		destIDFieldIndex = output.getOrCreateField(destIDFieldName.data(), OGRFieldType::OFTInteger);
		destFieldIndex = output.getOrCreateField("ca11", OGRFieldType::OFTReal);
		maskds = new GDAL_DS<unsigned int>();
		maskds->open(rasterfilename);
		mask = maskds->readData(1);
		nodata = (unsigned int)maskds->getNoData(1);
	}
	void addPoint(OGRPoint* pt, double& ca)
	{

		if (pt->getX() < maskds->bound.MinX || pt->getX() > maskds->bound.MaxX ||
			pt->getY() < maskds->bound.MinY || pt->getY() > maskds->bound.MaxY)
			return;

		unsigned int xindex = (unsigned int)((pt->getX() - maskds->bound.MinX) / maskds->adfGeoTransform[1]);
		unsigned int yindex = (unsigned int)((maskds->bound.MaxY - pt->getY()) / maskds->adfGeoTransform[1]);

		if (xindex < 0 || xindex > maskds->ncols - 1 ||
			yindex < 0 || yindex > maskds->nrows - 1)
			return;
		int fid = mask[xindex + yindex * maskds->ncols];
		if (fid == nodata)
			return;
		OGRFeature* newfea = OGRFeature::CreateFeature(output.poLayer->GetLayerDefn());
		newfea->SetGeometry(pt);
		newfea->SetField(destIDFieldIndex, fid - 1);
		newfea->SetField(destFieldIndex, ca);
		output.poLayer->CreateFeature(newfea);
		OGRFeature::DestroyFeature(newfea);
	}
	~UrbanOutline()
	{
		delete[] mask;
		delete maskds;
	}
private:

};


void polyline2points(double LINE_SAMPLE_DIST, std::string infilename, std::vector<std::string> masks, std::vector<std::string> outputfiles)
{

	ShapeFile input(infilename.data());

	int srcFieldIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
	if (srcFieldIndex < 0)
		return;
	std::vector<UrbanOutline*> urbanOutlines;
	for (size_t i = 0; i < masks.size(); i++)
	{
		urbanOutlines.push_back(new UrbanOutline(masks[i], outputfiles[i], &input));
	}


	OGRFeature* poFeature;
	int fid = -1;
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{
		fid++;

		OGRGeometry* poGeometry = poFeature->GetGeometryRef();
		if (!poGeometry)
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		double value = poFeature->GetFieldAsDouble(srcFieldIndex);
		std::vector<OGRPoint> splitPoints;

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
		double val = poFeature->GetFieldAsDouble(srcFieldIndex);
		for (int igeom = 0; igeom < geomlist.size(); igeom++)
		{
			OGRGeometry* poGeometry = geomlist[igeom];
			OGREnvelope geobb;

			if (dynamic_cast<OGRLineString*>(poGeometry))
			{
				OGRLineString* linestr = (OGRLineString *)poGeometry;
				dissolvePolyline(linestr, LINE_SAMPLE_DIST, splitPoints);
			}
			else if (dynamic_cast<OGRCircularString*>(poGeometry))
			{
				OGRCircularString* circularstr = (OGRCircularString *)poGeometry;
				OGRLineString* linestr = dynamic_cast<OGRLineString*>(circularstr->getLinearGeometry());
				dissolvePolyline(linestr, LINE_SAMPLE_DIST, splitPoints);
			}

		}
		double allocatedCA = value / splitPoints.size();
		for (size_t i = 0; i < splitPoints.size(); i++)
		{
			for (int n_urban = 0; n_urban < urbanOutlines.size(); n_urban++)
			{
				urbanOutlines[n_urban]->addPoint(&splitPoints[i], allocatedCA);
			}

		}
		OGRFeature::DestroyFeature(poFeature);
		if (fid % 100 == 0)
			printf("%d\n", fid);
	}
	for (size_t i = 0; i < masks.size(); i++)
	{
		delete urbanOutlines[i];
	}
}
void dropFields(std::string srcdir, std::vector<std::string> fields2keep)
{

	QDir input_dir(srcdir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	srcdir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".dbf"))
		{
			continue;
		}
		std::string shpfile = (input_dir.absolutePath() + "/" + fileInfo.completeBaseName()).toLocal8Bit().data() + std::string(".dbf");
		std::vector<std::string> newfields = fields2keep;
		ShapeFile shp(shpfile);
		for (int n = 0; n < shp.poLayer->GetLayerDefn()->GetFieldCount(); n++)
		{
			std::string name = shp.poLayer->GetLayerDefn()->GetFieldDefn(n)->GetNameRef();
			if (name.substr(0, 3) == "FID")
			{
				newfields.push_back(name);
			}
		}
		std::string infile = fileInfo.absoluteFilePath().toLocal8Bit().data();
		std::stringstream ss;
		ss << "RScript B:/Shapefiles2Grid/copydbf.R " << infile << " " << infile;
		for (size_t j = 0; j < newfields.size(); j++)
		{
			ss << " " << newfields[j];
		}
		system(ss.str().data());
	}
}
void copyFields(std::string srcdir, std::string destdir, std::vector<std::string> fields2keep)
{

	QDir input_dir(srcdir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	srcdir = (input_dir.absolutePath() + "/").toLocal8Bit().data();

	QDir output_dir(destdir.data());
	destdir = (output_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string name = fileInfo.fileName().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".dbf"))
		{
			continue;
		}
		//FID_ResNon
		std::string destID = "FID_";
		if (name.size() < 6)
			destID = destID + name;
		else
			destID = destID + name.substr(0, 6);
		std::string srcID = "FID";
		std::string infile = fileInfo.absoluteFilePath().toLocal8Bit().data();
		std::string outfile = (output_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		std::stringstream ss;
		ss << "RScript B:/Shapefiles2Grid/copyfieldmatch.R " << infile << " " << outfile << " " << srcID << " " << destID;
		for (size_t j = 0; j < fields2keep.size(); j++)
		{
			ss << " " << fields2keep[j];
		}
		printf("%s\n", ss.str().data());
		system(ss.str().data());
	}
}
void checkVulcanFiles()
{
	std::ifstream ifs;
	ifs.open("E:/Vulcan/gridPrep_SHP_master/Vulcan2SpatialCrosswalk.csv");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		ShapeFile dbf("E:/Vulcan/gridPrep_SHP_master/emissions/" + splits[1] + ".dbf");
		ShapeFile shp("E:/Vulcan/gridPrep_SHP_master/spatial/" + splits[3] + ".dbf");
		if (dbf.poLayer->GetFeatureCount() != shp.poLayer->GetFeatureCount())
		{
			printf("%s\n", line.data());
		}
	}
}
void linkUrbanOutline(std::string outline_file, std::string griddedDir)
{

	GDAL_DS<unsigned short>* ds = new GDAL_DS<unsigned short>();
	GDAL_DS<double>* dsArea = new GDAL_DS<double>();
	ds->open(outline_file);
	dsArea->setDSInfo(ds);
	dsArea->create(griddedDir + "Area.tif");

	double* areadata = new double[ds->slice];
	unsigned short* outlineDS = ds->readData(1);
	unsigned short min = 65535;
	unsigned short max = 0;
	double cellarea = dsArea->adfGeoTransform[1] * dsArea->adfGeoTransform[1] / 1000000.0;
	for (int i = 0; i < ds->slice; i++)
	{
		unsigned short val = outlineDS[i];
		areadata[i] = 0;
		if (val == 65535)
			continue;
		if (min > val)
			min = val;
		if (max < val)
			max = val;
		areadata[i] = cellarea;
	}
	dsArea->writeData(1,areadata, 0);
	delete dsArea;
	delete[] areadata;
	delete ds;
	//std::string outname = QFileInfo(outline_file.data()).completeBaseName().toLocal8Bit().data();
	std::string outname = QFileInfo(griddedDir.substr(0, griddedDir.size()-1).data()).completeBaseName().toLocal8Bit().data();

	std::ofstream ofs;
	ofs.open(griddedDir + outname + ".csv");

	std::vector<std::string> sectornames;
	std::vector<double*> totals;

	QDir input_dir(griddedDir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		if (!fileInfo.fileName().endsWith(".tif", Qt::CaseInsensitive))
			continue;
		sectornames.push_back(fileInfo.completeBaseName().toLocal8Bit().data());
		GDAL_DS<double>* dsFFCO2 = new GDAL_DS<double>();
		dsFFCO2->open(fileInfo.absoluteFilePath().toLocal8Bit().data());
		double* ffco2DS = dsFFCO2->readData(1);
		double nodata = dsFFCO2->getNoData(1);
		double* sumByCity = new double[max + 1];
		for (int n = 0; n < max + 1; n++)
		{
			sumByCity[n] = 0;
		}
		for (int n = 0; n < dsFFCO2->slice; n++)
		{
			double val = ffco2DS[n];
			if (val == nodata)
				continue;
			unsigned short urbanID = outlineDS[n];
			if (urbanID == 65535)
				continue;
			sumByCity[urbanID] += val;
		}
		totals.push_back(sumByCity);
		delete dsFFCO2;
		delete ffco2DS;
	}

	ofs << "Id";
	for (int i = 0; i < sectornames.size(); i++)
	{
		ofs << "," + sectornames[i];
	}
	ofs << std::endl;


	for (int ncity = 0; ncity < max + 1; ncity++)
	{
		ofs << ncity;
		for (int sector = 0; sector < totals.size(); sector++)
		{
			ofs << std::fixed << std::setprecision(5) << "," << totals[sector][ncity];
		}
		ofs << std::endl;
	}

	for (int i = 0; i < sectornames.size(); i++)
	{
		delete[] totals[i];
	}
	ofs.close();
	delete[] outlineDS;
}

void raster2csv(std::string indir)
{

	std::vector<std::string> files;
	QDir input_dir(indir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	double sum = 0;
	std::vector<std::string> csvfiles;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".tif", Qt::CaseInsensitive))
			continue;

		std::string outfile = input_file + ".csv";
		std::string name = fileInfo.completeBaseName().toLocal8Bit().data();
		std::ofstream ofs;
		ofs.open(outfile);

		GDAL_DS<double>* ds = new GDAL_DS<double>();
		ds->open(input_file);
		double* data = ds->readData(1);
		double nodata = ds->getNoData(1);
		double* pdata = data;
		ofs << name << std::endl;
		for (int n = 0; n < ds->slice; n++)
		{
			double val = *pdata;
			if (val == nodata || val < 0)
				val = 0;
			ofs << std::fixed << std::setprecision(5) << val << std::endl;
			pdata++;
		}
		ofs.close();
		csvfiles.push_back(outfile);
		delete ds;
		delete[] data;
		printf("%s\n", input_file);
		//ofs << std::fixed << std::setprecision(5) << fileInfo.completeBaseName().toLocal8Bit().data() << "," << subtotal /*/ 1000000.0*/ << std::endl;
	}


}
struct CumsumObject
{
	double a;
	double b;
};
bool sortCumsumObject(CumsumObject i, CumsumObject j) { return i.a < j.a; }

void cumsum(std::string infile1, std::string infile2, std::string outfile)
{

	GDAL_DS<double>* ds1 = new GDAL_DS<double>();
	ds1->open(infile1);
	double* data1 = ds1->readData(1);
	double nodata1 = ds1->getNoData(1);

	GDAL_DS<double>* ds2 = new GDAL_DS<double>();
	ds2->open(infile2);
	double* data2 = ds2->readData(1);
	double nodata2 = ds2->getNoData(1);
	std::vector<CumsumObject> cumsumobjs;
	double suma = 0;
	double sumb = 0;
	size_t num = ds1->slice;
	for (size_t n = 0; n < num; n++)
	{
		double a = data1[n];
		if (a == nodata1 || a < 0)
			a = 0;
		double b = data2[n];
		if (b == nodata2 || b < 0)
			b = 0;
		suma += a;
		sumb += b;
		CumsumObject obj;
		obj.a = a;
		obj.b = b;
		cumsumobjs.push_back(obj);
	}
	sort(cumsumobjs.begin(), cumsumobjs.end(), sortCumsumObject);

	CumsumObject cumsumpercent[1000];
	for (size_t n = 0; n < 1000; n++)
	{
		cumsumpercent[n].a = -1;
	}
	double cuma = 0;
	double cumb = 0;
	printf("%f,%f\n", suma, sumb);
	for (size_t k = 0; k < ds1->slice; k++)
	{
		double a = cumsumobjs[k].a;
		double b = cumsumobjs[k].b;
		cuma += a;
		cumb += b;
		CumsumObject per;
		per.a = cuma / suma;
		per.b = cumb / sumb;
		int thousper = (int)(per.a * 1000);
		if (cumsumpercent[thousper].a < 0)
		{
			cumsumpercent[thousper] = per;
		}
	}
	delete ds1;
	delete[] data1;
	delete ds2;
	delete[] data2;
	std::ofstream ofs;
	ofs.open(outfile);

	ofs << QFileInfo(infile1.data()).completeBaseName().toLocal8Bit().data() << "," << QFileInfo(infile1.data()).completeBaseName().toLocal8Bit().data() << std::endl;
	for (size_t n = 0; n < 1000; n++)
	{
		ofs << std::fixed << std::setprecision(5) << cumsumpercent[n].a << "," << cumsumpercent[n].b << std::endl;
	}
	ofs.close();

}

double calUrbanArea(std::string rasterfile)
{
	GDAL_DS<unsigned short>* ds1 = new GDAL_DS<unsigned short>();
	ds1->open(rasterfile);
	unsigned short* data1 = ds1->readData(1);
	unsigned short nodata1 = (unsigned short)ds1->getNoData(1);
	double cellarea = ds1->adfGeoTransform[1] * ds1->adfGeoTransform[1] / 1000000;
	double area1 = 0;
	for (size_t i = 0; i < ds1->slice; i++)
	{
		unsigned short val = data1[i];
		if (val == nodata1)
			continue;
		area1 += cellarea;
	}
	printf("%f,%s\n", area1, rasterfile.data());
	return area1;
}
void linkYuyu2ACS(std::string yuyuFile, std::string acsFile, std::string yuyuShapefile, std::string acsShapefile)
{
	GDAL_DS<unsigned short>* ds1 = new GDAL_DS<unsigned short>();
	ds1->open(yuyuFile);
	unsigned short* data1 = ds1->readData(1);
	unsigned short nodata1 = (unsigned short)ds1->getNoData(1);
	unsigned short max1 = 0;
	double cellarea = ds1->adfGeoTransform[1] * ds1->adfGeoTransform[1] / 1000000;
	double area1 = 0;
	for (size_t i = 0; i < ds1->slice; i++)
	{
		unsigned short val = data1[i];
		if (val == nodata1)
			continue;
		area1 += cellarea;
		if (max1 < val)
			max1 = val;
	}

	GDAL_DS<unsigned short>* ds2 = new GDAL_DS<unsigned short>();
	ds2->open(acsFile);
	unsigned short* data2 = ds2->readData(1);
	unsigned short nodata2 = (unsigned short)ds2->getNoData(1);
	double area2 = 0;
	for (size_t i = 0; i < ds2->slice; i++)
	{
		unsigned short val = data2[i];
		if (val == nodata2)
			continue;
		area2 += cellarea;

	}
	//printf("%f,%f\n", area1, area2);
	std::vector<std::map<int, int>> matchlist;
	max1 = max1 + 1;
	matchlist.resize(max1);
	for (size_t i = 0; i < ds1->slice; i++)
	{
		unsigned short val1 = data1[i];
		if (val1 == nodata1)
			continue;
		unsigned short val2 = data2[i];
		if (val2 == nodata2)
			continue;
		//if (val1 <= matchlist.size() - 1)
		//	continue;
		std::map<int, int>& map1 = matchlist[val1];

		std::map<int, int>::iterator iter = map1.find(val2);
		if (iter == map1.end())
		{
			map1[val2] = 1;
		}
		else
		{
			iter->second++;
		}
	}
	std::map<int, int> yuyu2acs;
	for (size_t i = 0; i < ds1->slice; i++)
	{
		unsigned short val1 = data1[i];
		if (val1 == nodata1)
			continue;
		std::map<int, int>& map1 = matchlist[val1];
		std::map<int, int>::iterator iter = map1.begin();
		int max_occurance = 0;
		int max_occurance_id = nodata2;
		while (iter != map1.end())
		{
			if (iter->second > max_occurance)
			{
				max_occurance = iter->second;
				max_occurance_id = iter->first;
			}

			iter++;
		}
		yuyu2acs[val1] = max_occurance_id;
	}

	ShapeFile acsShp(acsShapefile);

	std::vector<int> destFieldIds;
	ShapeFile Shp(yuyuShapefile, 1);
	std::vector<OGRFieldDefn*> fdn;
	for (size_t i = 0; i < acsShp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
	{
		fdn.push_back(acsShp.poLayer->GetLayerDefn()->GetFieldDefn(i));
		//Shp.poLayer->CreateField(fdn[i]);
		int idx = Shp.poLayer->GetLayerDefn()->GetFieldIndex(fdn[i]->GetNameRef());
		if (idx < 0)
		{
			Shp.poLayer->CreateField(fdn[i]);
			idx = Shp.poLayer->GetLayerDefn()->GetFieldIndex(fdn[i]->GetNameRef());
		}
		destFieldIds.push_back(idx);
	}
	OGRFeature* poFeature;
	int acsId = Shp.getOrCreateField("acsId", OGRFieldType::OFTInteger);
	unsigned fid = -1;
	while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
	{
		fid++;
		int id = nodata1;
		std::map<int, int>::iterator iter = yuyu2acs.find(fid);
		if (iter != yuyu2acs.end())
		{
			id = iter->second;
		}
		poFeature->SetField(acsId, id);
		if (id != nodata1)
		{
			OGRFeature* poFeature2 = acsShp.poLayer->GetFeature(id);
			for (size_t i = 0; i < destFieldIds.size(); i++)
			{
				poFeature->SetField(destFieldIds[i], poFeature2->GetRawFieldRef(i));
			}
			OGRFeature::DestroyFeature(poFeature2);
		}
		Shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
	Shp.close();
}
void remapYuyu(std::string inshpfile, std::string outshpfile, std::string intiffile, std::string outtiffile)
{

	ShapeFile inshp(inshpfile);
	ShapeFile outshp;
	outshp.create(outshpfile, inshp.poLayer->GetSpatialRef(), inshp.poLayer->GetLayerDefn(), OGRwkbGeometryType::wkbMultiPolygon);
	OGRFeature* poFeature;
	int acsIdField = inshp.getOrCreateField("acsId", OGRFieldType::OFTInteger);
	unsigned fid = -1;
	std::vector<int> acsIds[5000];
	std::vector<int> singleIds;
	while ((poFeature = inshp.poLayer->GetNextFeature()) != NULL)
	{
		fid++;
		int acsId = poFeature->GetFieldAsInteger(acsIdField);
		if (acsId != 65535)
		{
			acsIds[acsId].push_back(fid);
		}
		else
		{
			singleIds.push_back(fid);
		}
		OGRFeature::DestroyFeature(poFeature);
	}
	fid = -1;
	std::map<int, int> fidMap;
	for (size_t i = 0; i < 5000; i++)
	{
		std::vector<int>& idgroup = acsIds[i];
		if (idgroup.size() == 0)
			continue;
		fid++;
		OGRFeature* mergedfea = OGRFeature::CreateFeature(inshp.poLayer->GetLayerDefn());

		OGRMultiPolygon* multipoly = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(OGRwkbGeometryType::wkbMultiPolygon);
		for (size_t j = 0; j< idgroup.size(); j++)
		{
			OGRFeature* fea = inshp.poLayer->GetFeature(idgroup[j]);
			fidMap[idgroup[j]] = fid;
			OGRGeometry* srcGeom = fea->GetGeometryRef();
			if (srcGeom->getGeometryType() == OGRwkbGeometryType::wkbPolygon)
			{
				multipoly->addGeometry(srcGeom);
				//printf("wkbPolygon\n");
			}
			else if (srcGeom->getGeometryType() == OGRwkbGeometryType::wkbMultiPolygon)
			{
				OGRMultiPolygon* srcMultipoly = (OGRMultiPolygon*)fea->GetGeometryRef();
				for (size_t n = 0; n < srcMultipoly->getNumGeometries(); n++)
				{
					multipoly->addGeometry(srcMultipoly->getGeometryRef(n));
				}
				//printf("wkbMultiPolygon\n");
			}

			if (j == 0)
			{
				for (size_t n = 0; n < fea->GetFieldCount(); n++)
				{
					mergedfea->SetField(n, fea->GetRawFieldRef(n));
				}
			}
			OGRFeature::DestroyFeature(fea);
		}

		mergedfea->SetGeometry(multipoly);
		outshp.poLayer->CreateFeature(mergedfea);
		OGRFeature::DestroyFeature(mergedfea);
	}

	for (size_t i = 0; i < singleIds.size(); i++)
	{
		fid++;
		fidMap[singleIds[i]] = fid;
		OGRFeature* mergedfea = OGRFeature::CreateFeature(inshp.poLayer->GetLayerDefn());
		OGRMultiPolygon* multipoly = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(OGRwkbGeometryType::wkbMultiPolygon);
		OGRFeature* fea = inshp.poLayer->GetFeature(singleIds[i]);
		OGRGeometry* srcGeom = fea->GetGeometryRef();
		if (srcGeom->getGeometryType() == OGRwkbGeometryType::wkbPolygon)
		{
			multipoly->addGeometry(srcGeom);
			//printf("wkbPolygon\n");
		}
		else if (srcGeom->getGeometryType() == OGRwkbGeometryType::wkbMultiPolygon)
		{
			OGRMultiPolygon* srcMultipoly = (OGRMultiPolygon*)fea->GetGeometryRef();
			for (size_t n = 0; n < srcMultipoly->getNumGeometries(); n++)
			{
				multipoly->addGeometry(srcMultipoly->getGeometryRef(n));
			}
			//printf("wkbMultiPolygon\n");
		}

		for (size_t n = 0; n < fea->GetFieldCount(); n++)
		{
			mergedfea->SetField(n, fea->GetRawFieldRef(n));
		}
		mergedfea->SetGeometry(multipoly);
		OGRFeature::DestroyFeature(fea);
		outshp.poLayer->CreateFeature(mergedfea);
		OGRFeature::DestroyFeature(mergedfea);
	}


	GDAL_DS<unsigned short>* dsin = new GDAL_DS<unsigned short>();
	dsin->open(intiffile);
	unsigned short nodata = (unsigned short)dsin->getNoData(1);
	unsigned short* data1 = dsin->readData(1);
	dsin->create(outtiffile);
	for (size_t i = 0; i < dsin->slice; i++)
	{
		unsigned short val = data1[i];
		if (val == nodata)
			continue;
		std::map<int, int>::iterator iter = fidMap.find(val);
		data1[i] = iter->second;
	}
	dsin->writeData(1, data1, nodata);
	delete[] data1;
	delete dsin;
	/*QDir input_dir("C:/VulcanGridding/urban_outlines/Yuyu/");
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	std::string srcdir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
	QFileInfo fileInfo = list.at(i);
	std::string name = fileInfo.fileName().toLocal8Bit().data();
	if (!fileInfo.fileName().endsWith(".dbf"))
	{
	continue;
	}
	ShapeFile shp(fileInfo.absoluteFilePath().toLocal8Bit().data(), 1);
	printf("%s\n", fileInfo.absoluteFilePath().toLocal8Bit().data());
	int gridField = shp.poLayer->GetLayerDefn()->GetFieldIndex("gridId");
	int feaField = shp.poLayer->GetLayerDefn()->GetFieldIndex("feaId");
	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
	int gridId = poFeature->GetFieldAsInteger(gridField);
	unsigned short feaId = data1[gridId];
	poFeature->SetField(feaField, feaId);
	shp.poLayer->SetFeature(poFeature);
	OGRFeature::DestroyFeature(poFeature);
	}
	shp.close();
	}*/

	//int gridField = output.getOrCreateField("gridId", OGRFieldType::OFTInteger);

}
#include "SparseFractionGrid.h"
std::vector<double> loadAttribute(std::string dbffile, std::string attribute)
{
	std::vector<double> result;
	ShapeFile shp(dbffile);
	OGRFeature* poFeature;
	int fid = shp.poLayer->GetLayerDefn()->GetFieldIndex(attribute.data());
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		double val = poFeature->GetFieldAsDouble(fid);
		if (val != val || val < 0 || isinf((float)val) || !isfinite(val))
		{
			val = 0;
		}
		result.push_back(val);
		OGRFeature::DestroyFeature(poFeature);
	}
	return result;
}
void reallocateEmissionsByPop(std::string urban_outline_file, std::string co2_sparse_grid_file, std::string pop_sparse_grid_file, std::string co2_sparse_dbf_file, std::string pop_sparse_dbf_file, std::string outtif)
{
	Grid outputGrid;
	outputGrid.fromFishnetRaster(urban_outline_file);
	outputGrid.reset();
	double _NODATA = 3.40282306074e+038;
	for (int icell = 0; icell < outputGrid.slice; icell++)
	{
		outputGrid.cells[icell] = _NODATA;
	}


	std::string outDir = QFileInfo(pop_sparse_grid_file.data()).absoluteDir().absolutePath().toLocal8Bit().data() + std::string("/");
	//std::string outtif = outDir + ".tif";


	std::vector<double> popDBF = loadAttribute(pop_sparse_dbf_file, "P0010001");
	std::vector<double> co2DBF = loadAttribute(co2_sparse_dbf_file, "Total");
	SparseFractionGrid co2grid;
	co2grid.open(co2_sparse_grid_file);
	SparseFractionGrid popgrid;
	popgrid.open(pop_sparse_grid_file);

	Grid popGrid;
	popGrid.fromFishnetRaster(urban_outline_file);
	popGrid.reset();

	for (int icell = 0; icell < popGrid.slice; icell++)
	{
		popGrid.cells[icell] = 0;
	}
	for (size_t n = 0; n < popgrid.cells.size(); n++)
	{
		SparseGridCell& cell = popgrid.cells[n];
		int gridId = cell.gridId;
		int feaId = cell.feaId;
		double frac = cell.fraction;
		double pop = popDBF[feaId] * frac;
		if (pop > 0 && pop < 1000000000)
		{
			popGrid.cells[gridId] += pop;
		}
	}

	//P0010001
	//B19127_001
	//Total
	std::vector<std::vector<int>> co2_fid2gridid;
	co2_fid2gridid.resize(co2DBF.size());
	for (size_t n = 0; n < co2grid.cells.size(); n++)
	{
		SparseGridCell& cell = co2grid.cells[n];
		co2_fid2gridid[cell.feaId].push_back(cell.gridId);
	}


	for (size_t fid = 0; fid < co2DBF.size(); fid++)
	{

		double totalco2 = co2DBF[fid];
		double totalpop = 0;
		std::vector<int> cellvec = co2_fid2gridid[fid];
		for (size_t ipopcell = 0; ipopcell < cellvec.size(); ipopcell++)
		{
			int cellid = cellvec[ipopcell];
			totalpop += popGrid.cells[cellid];
		}
		for (size_t ipopcell = 0; ipopcell < cellvec.size(); ipopcell++)
		{
			int cellid = cellvec[ipopcell];
			double frac = 1.0 / cellvec.size();
			if (totalpop > 0)
			{
				frac = popGrid.cells[cellid] / totalpop;
			}
			double co2 = totalco2 * frac;
			if (outputGrid.cells[cellid] == _NODATA)
				outputGrid.cells[cellid] = 0;
			outputGrid.cells[cellid] += co2;
		}

	}

	GDAL_DS<double>* ds = new GDAL_DS<double>();
	ds->open(urban_outline_file);
	ds->create(outtif);
	ds->writeData(1, outputGrid.cells, _NODATA);
	delete ds;

}
void reallocateEmissionsByPop(std::string outline_name)
{
	std::string rootdir = "C:/VulcanGridding/urban_outlines/";
	reallocateEmissionsByPop(rootdir + outline_name + ".tif",
		rootdir + outline_name + "/fe_2007_us_zcta500.bin",
		rootdir + outline_name + "/us_census_bk.bin",
		"E:/Vulcan/gridPrep_SHP_master/emissions/fe_2007_us_zcta500.dbf",
		"E:/Vulcan/gridPrep_SHP_master/emissions/us_census_bk.dbf",
		rootdir + outline_name + "/Jones_Kammen_Total.tif");
	//linkUrbanOutline(rootdir + outline_name + ".tif",
	//	rootdir + outline_name);
}

void makeCommonGrid(std::vector<std::string> rasterFiles)
{
	GDAL_DS<unsigned short>* ds = new  GDAL_DS<unsigned short>();
	ds->open(rasterFiles[1]);
	OGREnvelope bound = ds->bound;
	for (size_t i = 0; i < rasterFiles.size(); i++)
	{
		GDAL_DS<unsigned short>* subds = new  GDAL_DS<unsigned short>();
		subds->open(rasterFiles[i]);
		bound.Merge(subds->bound);
		delete subds;
	}
	double* adftransform = ds->adfGeoTransform;
	double resol = 300;
	int nrows = (int)(ceil((bound.MaxY - bound.MinY) / resol));
	int ncols = (int)(ceil((bound.MaxX - bound.MinX) / resol));
	bound.MinY = bound.MaxY - resol * nrows;
	bound.MaxX = bound.MinX + resol * ncols;
	adftransform[1] = resol;
	adftransform[5] = -resol;
	adftransform[0] = bound.MinX;
	adftransform[3] = bound.MaxY;
	ds->nrows = nrows;
	ds->ncols = ncols;
	ds->setGeoTransform(adftransform);
	ds->create("C:/VulcanGridding/urban_outlines/commongrd.tif");
	delete ds;
}
//GDAL_DS<double>* ffdas = new GDAL_DS<double>();
//ffdas->open("B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/FFDAS_2011.tif");
//OGREnvelope bound;
//bound.MinX = -124.8;
//bound.MaxX = -66.9;
//bound.MinY = 24.5;
//bound.MaxY = 49.5;
//ffdas->crop(bound, "E:/Vulcan/gridPrep_SHP_master/spatial/FFDAS2011.tif");
//delete ffdas;
//
//GDAL_DS<double>* odiac = new GDAL_DS<double>();
//odiac->open("B:/Hestia_FFDAS_ODIAC/FFDAS_ODIAC/ODIAC_2011.tif");
//
//bound.MinX = -124.8;
//bound.MaxX = -66.9;
//bound.MinY = 24.5;
//bound.MaxY = 49.5;
//odiac->crop(bound, "E:/Vulcan/gridPrep_SHP_master/spatial/ODIAC2011.tif");
//delete odiac;

void makeTIFF(std::string indir, std::string outdir, std::string fishnetrasterfile)
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
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp") || fileInfo.fileName().endsWith("fishnet.shp"))
			continue;
		files.push_back(fileInfo.fileName().toLocal8Bit().data());
	}

	std::string outfishnetfile = outdir + "total.shp";
	std::string outrasterfile = outdir + "total.tif";

	std::string fishnetfile = (QFileInfo(fishnetrasterfile.data()).absoluteDir().absolutePath() + "/fishnet.shp").toLocal8Bit().data();
	std::string rasterfile = fishnetrasterfile;


	Grid fishnet;
	fishnet.fromFishnetRaster(fishnetrasterfile);
	fishnet.reset();

	GDALDataset* pDataset = (GDALDataset*)GDALOpen(fishnetrasterfile.data(), GA_ReadOnly);
	std::string wkt = pDataset->GetProjectionRef();
	GDALClose(pDataset);

	//if (!QFileInfo(fishnetfile.data()).exists())
	//{
	//	fishnet.toShape(wkt, fishnetfile, false);
	//}
	for (size_t i = 0; i < files.size(); i++)
	{
		/*if (skipNonRoad && files[i] == "NonRoad.shp")
		continue;*/
		printf("%s\n", (outdir + files[i]).data());
		//intersectWithFishnet(fishnet, fishnetfile, indir + files[i], outdir + files[i]);

		//if (!QFileInfo((outdir + files[i]).data()).exists())
		//{
		//	Utils::updateFootprint(indir + files[i]);
		//	Preprocessor::intersectWithArcGIS(indir + files[i], fishnetfile, outdir + files[i]);
		//	while (!QFileInfo((outdir + files[i]).data()).exists())
		//	{
		//		printf("ʧ��!");
		//		Preprocessor::intersectWithArcGIS(indir + files[i], fishnetfile, outdir + files[i]);
		//	}
		//	updateFieldAfterIntersection(outdir + files[i]);
		//}

		ShapeFile input((outdir + files[i]).data());
		fishnet.gatherCells(&input);
	}
	fishnet.toShape(wkt, outfishnetfile, true);
	fishnet.toRaster(outrasterfile, wkt);
}

void loadCities(std::string templateFile,std::string dir,std::string outfile)
{
	GDAL_DS<unsigned int>* ds = new GDAL_DS<unsigned int>();
	ds->open(templateFile);
	ds->create(outfile);
	unsigned int* cities = new unsigned int[ds->slice];
	for (size_t i = 0; i < ds->slice; i++)
	{
		cities[i] = 65535;
	}
	//memset(cities, 0, sizeof(unsigned int) * ds->slice);
	QDir input_dir(dir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	dir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	double resol = ds->adfGeoTransform[1];
	QFileInfoList list = input_dir.entryInfoList();
	ShapeFile attributes;
	attributes.create(outfile + ".shp",NULL,NULL, wkbNone);
	int popfield = attributes.getOrCreateField("pop", OGRFieldType::OFTReal);

	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		if (!fileInfo.fileName().endsWith(".dat"))
			continue;
		std::ifstream ifs;
		ifs.open(fileInfo.absoluteFilePath().toLocal8Bit().data());
		int rows, cols, elements;
		double minx, maxy;
		int row, col;
		double pop, area;
		std::string line;
		std::getline(ifs, line);
		std::stringstream(line) >> cols >> rows >> elements >> minx >> maxy;
		int cityid = atoi(fileInfo.completeBaseName().toLocal8Bit().data()) - 1;
		double popsum = 0;
		for (size_t j = 0; j < elements; j++)
		{
			std::getline(ifs, line);
			std::stringstream(line) >> col >> row >> pop >> area;
			popsum += pop;
			col = col - 1;
			row = row - 1;
			double x = minx + col * resol + resol * 0.5;
			double y = maxy - row * resol - resol * 0.5;
			int nx = (int)((x - ds->bound.MinX) / resol);
			int ny = (int)((ds->bound.MaxY - y) / resol);
			cities[nx + ny * ds->ncols] = cityid;
		}
		OGRFeature* fea = OGRFeature::CreateFeature(attributes.poLayer->GetLayerDefn());
		fea->SetField(popfield, popsum);
		attributes.poLayer->CreateFeature(fea);

		ifs.close();
		
	}
	attributes.close();
	ds->writeData(1, cities, 65535);
	delete ds;
	delete[] cities;
}

void loadGrump(std::string infile, std::string outfile)
{

	std::ifstream ifs;
	ifs.open(infile.data());
	int rows, cols, elements;
	double minx, maxy;
	int row, col;
	float pop;
	std::string line;
	std::getline(ifs, line);
	std::stringstream(line) >> cols >> rows >> elements >> minx >> maxy;
	GDAL_DS<unsigned int>* ds = new GDAL_DS<unsigned int>();
	double resol = 0.00833333;
	ds->nrows = rows;
	ds->ncols = cols;
	ds->slice = rows * cols;
	ds->numbands = 1;
	ds->adfGeoTransform[0] = minx;
	ds->adfGeoTransform[3] = maxy;
	ds->adfGeoTransform[1] = resol;
	ds->adfGeoTransform[5] = -resol;
	ds->setGeoTransform(ds->adfGeoTransform);
	//bound.MinX = adfGeoTransform[0];
	//bound.MaxY = adfGeoTransform[3];
	//bound.MaxX = adfGeoTransform[0] + adfGeoTransform[1] * ncols;
	//bound.MinY = adfGeoTransform[3] + adfGeoTransform[5] * nrows;

	unsigned int* pops = new unsigned int[ds->slice];
	memset(pops, 0, sizeof(unsigned int) * ds->slice);
	ds->create(outfile);
	double sum = 0;
	for (size_t j = 0; j < elements; j++)
	{
		std::getline(ifs, line);
		std::stringstream(line) >> col >> row >> pop;
		sum += pop;
		col = col - 1;
		row = row - 1;
		//double x = minx + col * resol + resol * 0.5;
		//double y = maxy - row * resol - resol * 0.5;
		//int nx = (int)((x - ds->bound.MinX) / resol);
		//int ny = (int)((ds->bound.MaxY - y) / resol);
		pops[col + row * ds->ncols] = pops[col + row * ds->ncols] + (unsigned int)pop;
	}
	printf("%f\n", sum);
	ds->writeData(1, pops, 0);
	delete ds;
	delete[] pops;

}
void scaleRaster(std::string filename, double scale)
{
	GDAL_DS<double>* ds = new GDAL_DS<double>();
	ds->open(filename,GDALAccess::GA_Update);
	ds->scale(1, scale);
	delete ds;
}

//void groupCCA()
//{
//	ShapeFile ccaShp("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.shp");
//
//	std::vector<int> destFieldIds;
//	ShapeFile Shp(yuyuShapefile, 1);
//	std::vector<OGRFieldDefn*> fdn;
//	for (size_t i = 0; i < acsShp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
//	{
//		fdn.push_back(acsShp.poLayer->GetLayerDefn()->GetFieldDefn(i));
//		//Shp.poLayer->CreateField(fdn[i]);
//		int idx = Shp.poLayer->GetLayerDefn()->GetFieldIndex(fdn[i]->GetNameRef());
//		if (idx < 0)
//		{
//			Shp.poLayer->CreateField(fdn[i]);
//			idx = Shp.poLayer->GetLayerDefn()->GetFieldIndex(fdn[i]->GetNameRef());
//		}
//		destFieldIds.push_back(idx);
//	}
//	OGRFeature* poFeature;
//	int acsId = Shp.getOrCreateField("acsId", OGRFieldType::OFTInteger);
//	unsigned fid = -1;
//	while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
//	{
//		fid++;
//		int id = nodata1;
//		std::map<int, int>::iterator iter = yuyu2acs.find(fid);
//		if (iter != yuyu2acs.end())
//		{
//			id = iter->second;
//		}
//		poFeature->SetField(acsId, id);
//		if (id != nodata1)
//		{
//			OGRFeature* poFeature2 = acsShp.poLayer->GetFeature(id);
//			for (size_t i = 0; i < destFieldIds.size(); i++)
//			{
//				poFeature->SetField(destFieldIds[i], poFeature2->GetRawFieldRef(i));
//			}
//			OGRFeature::DestroyFeature(poFeature2);
//		}
//		Shp.poLayer->SetFeature(poFeature);
//		OGRFeature::DestroyFeature(poFeature);
//	}
//	Shp.close();
//
//}
//

void filterWeatherStations()
{
	struct StationInfo
	{
		int idx;
		int cls;
		StationInfo()
		{

		}
		StationInfo(int _idx, std::string clsstr)
		{
			idx = _idx;
			if (clsstr == "I")
				cls = 1;
			else if (clsstr == "II")
				cls = 2;
			else
				cls = 3;
		}
	};
	ShapeFile inshp("E:/Vulcan/time/TMY3/TMY3_Vulcan_Proj.shp",1);
	OGRFeature* fea;
	int clsfield = inshp.poLayer->GetLayerDefn()->GetFieldIndex("Class");
	int filterfield = inshp.getOrCreateField("Filter",OGRFieldType::OFTInteger);
	
	OGREnvelope bound;
	inshp.poLayer->GetExtent(&bound);
	double resol = 150 * 1000;
	int numcols = (int)(ceil((bound.MaxX - bound.MinX) / resol));
	std::map<long, std::vector<StationInfo>> grid;
	int idx = 0;
	inshp.poLayer->ResetReading();
	while ((fea = inshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
		long row = (long)((bound.MaxY - p->getY()) / resol);
		long col = (long)((p->getX() - bound.MinX) / resol);
		long gridid = col + row * numcols;
		std::map<long, std::vector<StationInfo>>::iterator iter = grid.find(gridid);
		StationInfo info(idx,fea->GetFieldAsString(clsfield));
		if (iter == grid.end())
		{
			std::vector<StationInfo> newgroup;
			newgroup.push_back(info);
			grid[gridid] = newgroup;
		}
		else
		{
			iter->second.push_back(info);
		}
		idx++;
		OGRFeature::DestroyFeature(fea);
	}


	std::vector<int> selflags;
	selflags.resize(idx);
	memset(&selflags[0], 0, sizeof(int)*idx);
	std::map<long, std::vector<StationInfo>>::iterator iter = grid.begin();
	while (iter != grid.end())
	{
		std::vector<StationInfo> group = iter->second;
		int selCls = 4;
		int selID = 0;
		for (size_t i = 0; i < group.size(); i++)
		{
			if (group[i].cls < selCls)
			{
				selID = group[i].idx;
				selCls = group[i].cls;
			}
		}	
		selflags[selID] = 1;
		iter++;
	}
	idx = 0;
	inshp.poLayer->ResetReading();
	while ((fea = inshp.poLayer->GetNextFeature()) != NULL)
	{
		fea->SetField(filterfield, selflags[idx]);
		inshp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
		idx++;
	}
	inshp.close();

}

void matchAirports(std::string neifile, std::string nationalfile)
{
	ShapeFile neishp(neifile, 0);
	ShapeFile nationalshp(nationalfile, 1);
	int nei_airport_field = neishp.poLayer->GetLayerDefn()->GetFieldIndex("Facility_n");
	int national_airport_field = nationalshp.poLayer->GetLayerDefn()->GetFieldIndex("AIRPT_NAME");
	int facility_field = nationalshp.getOrCreateField("NEI_NAME",OGRFieldType::OFTString);
	int match_dist_field = nationalshp.getOrCreateField("NEI_DIST", OGRFieldType::OFTReal);

	OGRFeature* fea;
	nationalshp.poLayer->ResetReading();
	neishp.poLayer->ResetReading();
	std::vector<OGRPoint> points;
	std::vector<std::string> nei_airport_names;
	int fid2 = -1;
	while ((fea = neishp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
		OGRPoint newp;
		newp.setX(p->getX() * 100000);
		newp.setY(p->getY() * 100000);
		points.push_back(newp);
		nei_airport_names.push_back(fea->GetFieldAsString(nei_airport_field));
		fid2++;
		if(fid2 == 7535)
			printf(fea->GetFieldAsString(nei_airport_field));
		if (fea->GetFieldAsString(nei_airport_field) == std::string("Nantucket Memorial"))
			printf("%f,%f\n", p->getX(), p->getY());
		OGRFeature::DestroyFeature(fea);
	}
	neishp.close();
	ANNSearchEngine ann;
	int numresults = 30;
	ann.Create(points, numresults);
	std::vector<int> matchString;
	matchString.resize(numresults);
	std::vector<std::string> names;
	names.resize(numresults);


	OGRPoint queryP;
	int fid = 640;
	fea = nationalshp.poLayer->GetFeature(640);
	OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
	printf("%s\n", fea->GetFieldAsString(national_airport_field));
	OGRFeature::DestroyFeature(fea);
	double lon = p->getX();
	double lat = p->getY();
	queryP.setX(lon * 100000);
	queryP.setY(lat * 100000);
	ann.Select_Nearest_Points(queryP);
	for (int i = 0; i < numresults; i++)
	{
		double dist = 0;
		int selectedIdx = -1;
		ann.Get_Selected_Point(i, selectedIdx, dist);
		dist = sqrt(dist);
		std::string nei_name = nei_airport_names[selectedIdx];
		printf("%s,%f\n", nei_name.data(), dist);
	}

	int numunmatched = 0;
	while ((fea = nationalshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
		OGRPoint newp;
		newp.setX(p->getX() * 100000);
		newp.setY(p->getY() * 100000);
		 ann.Select_Nearest_Points(newp);
		std::string national_name = fea->GetFieldAsString(national_airport_field);
		national_name = QString(national_name.data()).toUpper().toLocal8Bit().data();
		if (national_name.size() > 7 && national_name.substr(national_name.size() - 7, 7) == "AIRPORT")
			national_name = national_name.substr(0, national_name.size() - 7);

		if (national_name == "NANTUCKET MEMORIAL ")
		{
			printf("");
		}
		int nearestStationID = -1;
		int maxmatches = -1;
		int matchid = -1;
		double seldist = 0;
		for (int i = 0; i < 3; i++)
		{
			double dist = 0;
			int selectedIdx = -1;
			ann.Get_Selected_Point(i, selectedIdx, dist);
			dist = sqrt(dist);
			matchString[i] = 0;
			std::string nei_name = nei_airport_names[selectedIdx];
			names[i] = nei_name;
			nei_name = QString(nei_name.data()).toUpper().toLocal8Bit().data();
			int numchars = nei_name.size();
			if (numchars > national_name.size())
				numchars = national_name.size();

			for (int nchar = 0; nchar < numchars; nchar++)
			{
				if (nei_name[nchar] != national_name[nchar])
					break;
				matchString[i]++;
			}

			if (matchString[i] > maxmatches)
			{
				maxmatches = matchString[i];
				matchid = i;
				seldist = dist;
			}
		}
		if (maxmatches < 2 /*|| seldist > 20*/)
		{
			maxmatches = -1;
			matchid = -1;
			seldist = 0;
			for (int i = 0; i < 3; i++)
			{
				double dist = 0;
				int selectedIdx = -1;
				ann.Get_Selected_Point(i, selectedIdx, dist);
				dist = sqrt(dist);
				matchString[i] = 0;
				std::string nei_name = nei_airport_names[selectedIdx];
				names[i] = nei_name;
				nei_name = QString(nei_name.data()).toUpper().toLocal8Bit().data();
				int numchars = nei_name.size();
				if (numchars > national_name.size())
					numchars = national_name.size();
		
				for (int m = 0; m < numchars; m++)
				{
					int matchcount = -1;
					for (int nchar = m; nchar < numchars; nchar++)
					{
						if (nei_name[nchar] != national_name[nchar])
							break;
						matchcount++;
					}
					if (matchString[i] < matchcount)
						matchString[i] = matchcount;
				}

				if (matchString[i] > maxmatches)
				{
					maxmatches = matchString[i];
					matchid = i;
					seldist = dist;
				}
			}
			if (maxmatches < 2 /*|| seldist > 0.3*/)
			{
				numunmatched++;
				fea->SetField(facility_field, "");
				fea->SetField(match_dist_field, 0);
			}
			else
			{
				fea->SetField(facility_field, names[matchid].data());
				fea->SetField(match_dist_field, -1);
			}

		}
		else
		{
			fea->SetField(facility_field, names[matchid].data());
			fea->SetField(match_dist_field, seldist);
		}

		nationalshp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}
	nationalshp.close();
	printf("%d\n", numunmatched);
}
void matchAirNavAirports(std::string neifile, std::string airnavfile)
{
	ShapeFile neishp(neifile, 0);
	ShapeFile nationalshp(airnavfile, 1);
	int nei_airport_field = neishp.poLayer->GetLayerDefn()->GetFieldIndex("Facility_n");
	int national_airport_field = nationalshp.poLayer->GetLayerDefn()->GetFieldIndex("AirportNam");
	int facility_field = nationalshp.getOrCreateField("NEI_NAME", OGRFieldType::OFTString);
	int match_dist_field = nationalshp.getOrCreateField("NEI_DIST", OGRFieldType::OFTReal);

	OGRFeature* fea;
	nationalshp.poLayer->ResetReading();
	neishp.poLayer->ResetReading();
	std::vector<OGRPoint> points;
	std::vector<std::string> nei_airport_names;
	int fid2 = -1;
	while ((fea = neishp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
		OGRPoint newp;
		newp.setX(p->getX() * 100000);
		newp.setY(p->getY() * 100000);
		points.push_back(newp);
		nei_airport_names.push_back(fea->GetFieldAsString(nei_airport_field));
		fid2++;
		//if (fid2 == 7535)
		//	printf(fea->GetFieldAsString(nei_airport_field));
		//if (fea->GetFieldAsString(nei_airport_field) == std::string("Nantucket Memorial"))
		//	printf("%f,%f\n", p->getX(), p->getY());
		OGRFeature::DestroyFeature(fea);
	}
	neishp.close();
	ANNSearchEngine ann;
	int numresults = 30;
	ann.Create(points, numresults);
	std::vector<int> matchString;
	matchString.resize(numresults);
	std::vector<std::string> names;
	names.resize(numresults);


	OGRPoint queryP;
	int fid = 640;
	fea = nationalshp.poLayer->GetFeature(640);
	OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
	printf("%s\n", fea->GetFieldAsString(national_airport_field));
	OGRFeature::DestroyFeature(fea);
	double lon = p->getX();
	double lat = p->getY();
	queryP.setX(lon * 100000);
	queryP.setY(lat * 100000);
	ann.Select_Nearest_Points(queryP);
	for (int i = 0; i < numresults; i++)
	{
		double dist = 0;
		int selectedIdx = -1;
		ann.Get_Selected_Point(i, selectedIdx, dist);
		dist = sqrt(dist);
		std::string nei_name = nei_airport_names[selectedIdx];
		printf("%s,%f\n", nei_name.data(), dist);
	}

	int numunmatched = 0;
	while ((fea = nationalshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
		OGRPoint newp;
		newp.setX(p->getX() * 100000);
		newp.setY(p->getY() * 100000);
		ann.Select_Nearest_Points(newp);
		std::string national_name = fea->GetFieldAsString(national_airport_field);
		bool ishe = false;
		//if (national_name == "Charles M Schulz - Sonoma County Airport")
		//{
		//	ishe = true;
		//	printf("");
		//}
		national_name = QString(national_name.data()).toUpper().toLocal8Bit().data();

		if (national_name.size() > 7 && national_name.substr(national_name.size() - 7, 7) == "AIRPORT")
			national_name = national_name.substr(0, national_name.size() - 8);

		//if (national_name == "NANTUCKET MEMORIAL ")
		//{
		//	printf("");
		//}
		int nearestStationID = -1;
		int maxmatches = -1;
		int matchid = -1;
		double seldist = 0;

		for (int i = 0; i < 3; i++)
		{
			double dist = 0;
			int selectedIdx = -1;
			ann.Get_Selected_Point(i, selectedIdx, dist);
			dist = sqrt(dist);
			//if (dist > 10000)
			//	continue;
			matchString[i] = 0;
			std::string nei_name = nei_airport_names[selectedIdx];
			names[i] = nei_name;
			nei_name = QString(nei_name.data()).toUpper().toLocal8Bit().data();
			//if (ishe)
			//{
			//	printf("%s,%s\n", nei_name.data(), national_name.data());
			//}
			int numchars = nei_name.size();
			if (numchars > national_name.size())
				numchars = national_name.size();

			for (int nchar = 0; nchar < numchars; nchar++)
			{
				if (nei_name[nchar] != national_name[nchar])
					break;
				matchString[i]++;
			}
			if (nei_name.find(national_name) != std::string::npos || national_name.find(nei_name) != std::string::npos)
			{
				matchString[i] = numchars;
			}
			if (matchString[i] > maxmatches)
			{
				maxmatches = matchString[i];
				matchid = i;
				seldist = dist;
			}
		}
		//if (national_name == QString("Hillman").toUpper().toLocal8Bit().data())
		//{
		//	printf("");
		//}

		//if (fea->GetFieldAsInteger("ID") == 2060)
		//{
		//	printf("");
		//}
		int minlen = min(names[matchid].size(), national_name.size());
		if (maxmatches < minlen / 2 /*|| seldist > 20*/)
		{
			//maxmatches = -1;
			//matchid = -1;
			//seldist = 0;
			//for (int i = 0; i < 3; i++)
			//{
			//	double dist = 0;
			//	int selectedIdx = -1;
			//	ann.Get_Selected_Point(i, selectedIdx, dist);
			//	dist = sqrt(dist);
			//	matchString[i] = 0;
			//	std::string nei_name = nei_airport_names[selectedIdx];
			//	names[i] = nei_name;
			//	nei_name = QString(nei_name.data()).toUpper().toLocal8Bit().data();
			//	int numchars = nei_name.size();
			//	if (numchars > national_name.size())
			//		numchars = national_name.size();

			//	for (int m = 0; m < numchars; m++)
			//	{
			//		int matchcount = -1;
			//		for (int nchar = m; nchar < numchars; nchar++)
			//		{
			//			if (nei_name[nchar] != national_name[nchar])
			//				break;
			//			matchcount++;
			//		}
			//		if (matchString[i] < matchcount)
			//			matchString[i] = matchcount;
			//	}

			//	if (matchString[i] > maxmatches)
			//	{
			//		maxmatches = matchString[i];
			//		matchid = i;
			//		seldist = dist;
			//	}
			//}
			//if (maxmatches < 2 /*|| seldist > 0.3*/)
			//{
			//	numunmatched++;
			//	fea->SetField(facility_field, "");
			//	fea->SetField(match_dist_field, 0);
			//}
			//else
			//{
			//	fea->SetField(facility_field, names[matchid].data());
			//	fea->SetField(match_dist_field, -1);
			//}
			fea->SetField(facility_field, "");
			fea->SetField(match_dist_field, -1);
		}
		else
		{
			fea->SetField(facility_field, names[matchid].data());
			fea->SetField(match_dist_field, seldist);
		}

		nationalshp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}
	nationalshp.close();
	printf("%d\n", numunmatched);
}
void matchInternationalAirports(std::string neifile, std::string airnavfile)
{
	ShapeFile neishp(neifile, 0);
	ShapeFile nationalshp(airnavfile, 1);
	int nei_airport_field = neishp.poLayer->GetLayerDefn()->GetFieldIndex("Facility_n");
	int nei_airport_id = neishp.poLayer->GetLayerDefn()->GetFieldIndex("Facility_I");
	int national_airport_field = nationalshp.poLayer->GetLayerDefn()->GetFieldIndex("FacilityNa");
	int facility_field = nationalshp.getOrCreateField("NEI_NAME", OGRFieldType::OFTString);
	int match_dist_field = nationalshp.getOrCreateField("NEI_DIST", OGRFieldType::OFTReal);

	OGRFeature* fea;
	nationalshp.poLayer->ResetReading();
	neishp.poLayer->ResetReading();
	std::vector<OGRPoint> points;
	std::vector<std::string> nei_airport_names;
	int fid2 = -1;
	while ((fea = neishp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
		OGRPoint newp;
		newp.setX(p->getX() * 100000);
		newp.setY(p->getY() * 100000);
		points.push_back(newp);
		nei_airport_names.push_back(fea->GetFieldAsString(nei_airport_field));
		fid2++;
		//if (fid2 == 7535)
		//	printf(fea->GetFieldAsString(nei_airport_field));
		//if (fea->GetFieldAsString(nei_airport_field) == std::string("Nantucket Memorial"))
		//	printf("%f,%f\n", p->getX(), p->getY());
		OGRFeature::DestroyFeature(fea);
	}
	neishp.close();
	ANNSearchEngine ann;
	int numresults = 30;
	ann.Create(points, numresults);
	std::vector<int> matchString;
	matchString.resize(numresults);
	std::vector<std::string> names;
	names.resize(numresults);


	OGRPoint queryP;
	//int fid = 640;
	//fea = nationalshp.poLayer->GetFeature(640);
	//OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
	//printf("%s\n", fea->GetFieldAsString(national_airport_field));
	//OGRFeature::DestroyFeature(fea);
	//double lon = p->getX();
	//double lat = p->getY();
	//queryP.setX(lon * 100000);
	//queryP.setY(lat * 100000);
	//ann.Select_Nearest_Points(queryP);
	//for (int i = 0; i < numresults; i++)
	//{
	//	double dist = 0;
	//	int selectedIdx = -1;
	//	ann.Get_Selected_Point(i, selectedIdx, dist);
	//	dist = sqrt(dist);
	//	std::string nei_name = nei_airport_names[selectedIdx];
	//	printf("%s,%f\n", nei_name.data(), dist);
	//}

	int numunmatched = 0;
	while ((fea = nationalshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* p = (OGRPoint*)fea->GetGeometryRef();
		OGRPoint newp;
		newp.setX(p->getX() * 100000);
		newp.setY(p->getY() * 100000);
		ann.Select_Nearest_Points(newp);
		std::string national_name = fea->GetFieldAsString(national_airport_field);
		bool ishe = false;
		//if (national_name == "Charles M Schulz - Sonoma County Airport")
		//{
		//	ishe = true;
		//	printf("");
		//}
		national_name = QString(national_name.data()).toUpper().toLocal8Bit().data();

		if (national_name.size() > 7 && national_name.substr(national_name.size() - 7, 7) == "AIRPORT")
			national_name = national_name.substr(0, national_name.size() - 8);

		//if (national_name == "NANTUCKET MEMORIAL ")
		//{
		//	printf("");
		//}
		int nearestStationID = -1;
		int maxmatches = -1;
		int matchid = -1;
		double seldist = 0;

		for (int i = 0; i < 3; i++)
		{
			double dist = 0;
			int selectedIdx = -1;
			ann.Get_Selected_Point(i, selectedIdx, dist);
			dist = sqrt(dist);
			//if (dist > 10000)
			//	continue;
			matchString[i] = 0;
			std::string nei_name = nei_airport_names[selectedIdx];
			names[i] = nei_name;
			nei_name = QString(nei_name.data()).toUpper().toLocal8Bit().data();
			//if (ishe)
			//{
			//	printf("%s,%s\n", nei_name.data(), national_name.data());
			//}
			int numchars = nei_name.size();
			if (numchars > national_name.size())
				numchars = national_name.size();

			for (int nchar = 0; nchar < numchars; nchar++)
			{
				if (nei_name[nchar] != national_name[nchar])
					break;
				matchString[i]++;
			}
			if (nei_name.find(national_name) != std::string::npos || national_name.find(nei_name) != std::string::npos)
			{
				matchString[i] = numchars;
			}
			if (matchString[i] > maxmatches)
			{
				maxmatches = matchString[i];
				matchid = i;
				seldist = dist;
			}
		}
		//if (national_name == QString("Hillman").toUpper().toLocal8Bit().data())
		//{
		//	printf("");
		//}

		//if (fea->GetFieldAsInteger("ID") == 2060)
		//{
		//	printf("");
		//}
		int minlen = min(names[matchid].size(), national_name.size());
		if (maxmatches < minlen / 2 /*|| seldist > 20*/)
		{
			//maxmatches = -1;
			//matchid = -1;
			//seldist = 0;
			//for (int i = 0; i < 3; i++)
			//{
			//	double dist = 0;
			//	int selectedIdx = -1;
			//	ann.Get_Selected_Point(i, selectedIdx, dist);
			//	dist = sqrt(dist);
			//	matchString[i] = 0;
			//	std::string nei_name = nei_airport_names[selectedIdx];
			//	names[i] = nei_name;
			//	nei_name = QString(nei_name.data()).toUpper().toLocal8Bit().data();
			//	int numchars = nei_name.size();
			//	if (numchars > national_name.size())
			//		numchars = national_name.size();

			//	for (int m = 0; m < numchars; m++)
			//	{
			//		int matchcount = -1;
			//		for (int nchar = m; nchar < numchars; nchar++)
			//		{
			//			if (nei_name[nchar] != national_name[nchar])
			//				break;
			//			matchcount++;
			//		}
			//		if (matchString[i] < matchcount)
			//			matchString[i] = matchcount;
			//	}

			//	if (matchString[i] > maxmatches)
			//	{
			//		maxmatches = matchString[i];
			//		matchid = i;
			//		seldist = dist;
			//	}
			//}
			//if (maxmatches < 2 /*|| seldist > 0.3*/)
			//{
			//	numunmatched++;
			//	fea->SetField(facility_field, "");
			//	fea->SetField(match_dist_field, 0);
			//}
			//else
			//{
			//	fea->SetField(facility_field, names[matchid].data());
			//	fea->SetField(match_dist_field, -1);
			//}
			fea->SetField(facility_field, "");
			fea->SetField(match_dist_field, -1);
		}
		else
		{
			fea->SetField(facility_field, names[matchid].data());
			fea->SetField(match_dist_field, seldist);
		}

		nationalshp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}
	nationalshp.close();
	printf("%d\n", numunmatched);
}
void linkUrbanOutline2States(std::string outlineFile, std::string outlineShpFile, std::string stateFile, std::string mapfile)
{
	GDAL_DS<unsigned short>* ds1 = new GDAL_DS<unsigned short>();
	ds1->open(outlineFile);
	unsigned short* data1 = ds1->readData(1);
	unsigned short nodata1 = (unsigned short)ds1->getNoData(1);
	unsigned short max1 = 0;
	double cellarea = ds1->adfGeoTransform[1] * ds1->adfGeoTransform[1] / 1000000;
	double area1 = 0;
	for (size_t i = 0; i < ds1->slice; i++)
	{
		unsigned short val = data1[i];
		if (val == nodata1)
			continue;
		area1 += cellarea;
		if (max1 < val)
			max1 = val;
	}

	GDAL_DS<unsigned short>* ds2 = new GDAL_DS<unsigned short>();
	ds2->open(stateFile);
	unsigned short* data2 = ds2->readData(1);
	unsigned short nodata2 = (unsigned short)ds2->getNoData(1);
	double area2 = 0;
	for (size_t i = 0; i < ds2->slice; i++)
	{
		unsigned short val = data2[i];
		if (val == nodata2)
			continue;
		area2 += cellarea;

	}
	//printf("%f,%f\n", area1, area2);
	std::vector<std::map<int, int>> matchlist;
	max1 = max1 + 1;
	matchlist.resize(max1);
	OGRSpatialReference oSourceSRS;
	char wkt[512];
	memcpy(wkt, ds1->projection.data(), ds1->projection.size());
	char* pwkt = wkt;
	oSourceSRS.importFromWkt(&pwkt);
	OGRSpatialReference oTargetSRS;
	oTargetSRS.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *poCT = OGRCreateCoordinateTransformation(&oSourceSRS, &oTargetSRS);
	double x = -3.17970006;
	double y = 51.4789137;

	double destx = x;
	double desty = y;
	poCT->Transform(1, &destx, &desty);
	for (size_t i = 0; i < ds1->slice; i++)
	{
		unsigned short val1 = data1[i];
		if (val1 == nodata1)
			continue;
		unsigned short val2 = data2[i];
		if (val2 == nodata2)
			continue;
		//if (val1 <= matchlist.size() - 1)
		//	continue;
		std::map<int, int>& map1 = matchlist[val1];

		std::map<int, int>::iterator iter = map1.find(val2);
		if (iter == map1.end())
		{
			map1[val2] = 1;
		}
		else
		{
			iter->second++;
		}
	}
	ShapeFile shp(outlineShpFile.data(), 1);
	int fipsField = shp.getOrCreateField("StateFIPS", OGRFieldType::OFTInteger);
	//printf("(%f,%f) -> (%f,%f)\n", x, y, destx, desty);
	std::map<int, int> yuyu2acs;
	//std::ofstream ofs(mapfile.data());
	//ofs << "ClusterID,StateFIPS,CentroidLongi,CentroidLati" << std::endl;
	for (int val1 = 0; val1 < matchlist.size(); val1++)
	{
		std::map<int, int>& map1 = matchlist[val1];
		std::map<int, int>::iterator iter = map1.begin();
		int max_occurance = 0;
		int max_occurance_id = nodata2;
		while (iter != map1.end())
		{
			if (iter->second > max_occurance && iter->first != 255)
			{
				max_occurance = iter->second;
				max_occurance_id = iter->first;
			}
			iter++;
		}
		//ofs << val1 << "," << max_occurance_id << std::endl;
		OGRFeature* fea = shp.poLayer->GetFeature(val1);
		fea->SetField(fipsField, max_occurance_id);
		shp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}
	shp.close();
	//ofs.close();
	delete ds1;
	delete ds2;
}

void remapCCA(std::string tiffile, std::string shpfile)
{
	GDAL_DS<unsigned short>* ds = new GDAL_DS<unsigned short>();
	ds->open(tiffile,GDALAccess::GA_Update);
	unsigned short* data = ds->readData(1);
	unsigned short nodata = (unsigned short)ds->getNoData(1);
	ShapeFile shp(shpfile.data());
	int grididfield = shp.getOrCreateField("gridcode",OGRFieldType::OFTInteger);
	OGRFeature* fea ;
	std::map<int, int> gridcodemap;
	int id = 0;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		gridcodemap[fea->GetFieldAsInteger(grididfield)] = id;
		id++;
		OGRFeature::DestroyFeature(fea);
	}

	for (size_t i = 0; i < ds->slice; i++)
	{
		unsigned short gridid = data[i];
		if (gridid == nodata)
			continue;
		if (gridcodemap.find(gridid) == gridcodemap.end())
		{
			printf("%d\n", gridid);
		}
		else
		{
			data[i] = gridcodemap[gridid];
		}
	}
	ds->writeData(1, data, nodata);
	delete[] data;
	shp.close();
	delete ds;
}

void gridVulcanLAOnroad()
{
	Grid fishnetGrid;
	std::string fishnetrasterfile = "B:/LA_Version2/Vulcan/fishnet100m.tif";
	std::string fishnetshapefile = "B:/LA_Version2/Vulcan/fishnet100m.shp";
	std::string intersectedShapeFile = "B:/LA_Version2/Vulcan/gridded/Vulcan_Onroad100m.shp";
	std::string vulcanonroadfile = "B:/LA_Version2/Vulcan/Vulcan_Onroad100m.shp";
	fishnetGrid.fromFishnetRaster(fishnetrasterfile);
	
	//Utils::updateFootprint(vulcanonroadfile);
	Preprocessor::intersectWithArcGIS(vulcanonroadfile, fishnetshapefile, intersectedShapeFile);
	updateFieldAfterIntersection(intersectedShapeFile);
	ShapeFile inshp(intersectedShapeFile);
	//ShapeFile inshp("B:/LA_Version2/LA_Basin_Bound.shp");

	fishnetGrid.gatherCells(&inshp, "ca11");
	fishnetGrid.toRaster("B:/LA_Version2/Vulcan/gridded/Vulcan_Onroad100m.tif", fishnetGrid.proj);
}
void gridHestiaLAOnroad()
{
	Grid fishnetGrid;
	std::string fishnetrasterfile = "B:/LA_Version2/Vulcan/fishnet100m.tif";
	std::string fishnetshapefile = "B:/LA_Version2/Vulcan/fishnet100m.shp";
	std::string intersectedShapeFile = "B:/LA_Version2/Vulcan/gridded/Hestia_Onroad100m.shp";
	std::string vulcanonroadfile = "B:/LA_Version2/Vulcan/Hestia_Onroad100m.shp";
	fishnetGrid.fromFishnetRaster(fishnetrasterfile);

	//Utils::updateFootprint(vulcanonroadfile);
	Preprocessor::intersectWithArcGIS(vulcanonroadfile, fishnetshapefile, intersectedShapeFile);
	updateFieldAfterIntersection(intersectedShapeFile);
	ShapeFile inshp(intersectedShapeFile);
	//ShapeFile inshp("B:/LA_Version2/LA_Basin_Bound.shp");

	fishnetGrid.gatherCells(&inshp, "ca11");
	fishnetGrid.toRaster("B:/LA_Version2/Vulcan/gridded/Hestia_Onroad100m.tif", fishnetGrid.proj);
}

void calculatePointCoordinates(std::string shpfile)
{
	
	ShapeFile shp(shpfile.data(),1);
	int fLatitude = shp.getOrCreateField("Latitude", OGRFieldType::OFTReal);
	int fLongitude = shp.getOrCreateField("Longitude", OGRFieldType::OFTReal);
	OGRFeature* fea;
	std::map<int, int> gridcodemap;
	int id = 0;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		fea->SetField(fLatitude, pt->getY());
		fea->SetField(fLongitude, pt->getX());
		shp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}

	shp.close();
	
}
#include "LosAngeles.h"
void copySCCProfiles(std::string filename)
{
	ShapeFile ElecProd(filename, 0);

	int idx = ElecProd.getOrCreateField("timestruct", OFTString);
	OGRFeatureDefn* featureDef = ElecProd.poLayer->GetLayerDefn();
	OGRFeature* poFeature;
	while ((poFeature = ElecProd.poLayer->GetNextFeature()) != NULL)
	{
		std::string scc = poFeature->GetFieldAsString("SCC");

		if (scc == "")
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		if (copy_SCC_PROFILE(scc, "9" + scc))
		{

			//poFeature->SetField(idx, ("9" + scc).data());
		}
		else
		{
			//poFeature->SetField(idx, "0");
		}
		//ElecProd.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
}



void findNearestRoadSegment(std::string srcfile,std::string destfile)
{

	std::vector<OGRPoint> points;
	std::vector<std::string> timestructs;
	ShapeFile srcroadshapes(srcfile, 0);
	OGRFeature *poFeature;
	srcroadshapes.poLayer->ResetReading();
	double maxDist = 10000 * 3.28084;
	int minnumpoints = 0;
	int maxnumpoints = 1;

	int timestructField = srcroadshapes.getOrCreateField("timestruct", OGRFieldType::OFTString);
	while ((poFeature = srcroadshapes.poLayer->GetNextFeature()) != NULL)
	{
		std::vector<OGRPoint> pointOnSegment;
		RoadTimeIDW::getPoints(poFeature->GetGeometryRef(), pointOnSegment);
		int centerIdx = pointOnSegment.size() / 2;
		OGRPoint center = pointOnSegment[centerIdx];
		points.push_back(center);
		timestructs.push_back(poFeature->GetFieldAsString(timestructField));
		OGRFeature::DestroyFeature(poFeature);
	}
	srcroadshapes.close();
	ANNSearchEngine searchEngine;
	searchEngine.Create(points, maxnumpoints);
	ShapeFile destroadshapes(destfile, 1);
	destroadshapes.poLayer->ResetReading();
	timestructField = destroadshapes.getOrCreateField("timestruct", OGRFieldType::OFTString);
	ShapeFile shp;
	shp.create("e:/idw_points.shp", destroadshapes.poLayer->GetSpatialRef(), 0, OGRwkbGeometryType::wkbPoint);
	int timestructField2 = shp.getOrCreateField("timestruct", OGRFieldType::OFTString);
	int idxField2 = shp.getOrCreateField("idx", OGRFieldType::OFTInteger);

	while ((poFeature = destroadshapes.poLayer->GetNextFeature()) != NULL)
	{
		std::vector<OGRPoint> pointOnSegment;
		RoadTimeIDW::getPoints(poFeature->GetGeometryRef(), pointOnSegment);
		int centerIdx = pointOnSegment.size() / 2;
		OGRPoint center = pointOnSegment[centerIdx];
		int numresults = searchEngine.Select_Nearest_Points(center);
		int selectedIdx = -1;
		double dist;
		searchEngine.Get_Selected_Point(0, selectedIdx, dist);
		poFeature->SetField(timestructField, timestructs[selectedIdx].data());
		OGRFeature* fea = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		OGRPoint pt;
		pt.setX(center.getX());
		pt.setY(center.getY());
		fea->SetField(idxField2, selectedIdx);
		fea->SetField(timestructField2, timestructs[selectedIdx].data());
		fea->SetGeometry(&pt);
		destroadshapes.poLayer->SetFeature(poFeature);
		shp.poLayer->CreateFeature(fea);
		OGRFeature::DestroyFeature(fea);
		OGRFeature::DestroyFeature(poFeature);
	}
	destroadshapes.close();

}

void calVulcanExtents() {
	//OGREnvelope bb = BoundManager::readBoundFromRaster("C:/VulcanGridding/VulcanGridContiguousUS/grid.tif");
	//bb.Merge(BoundManager::readBoundFromShape("E:/Vulcan/support_data/US_Boundaries/Vulcan/Contiguous_US_NOAA_Maritime_Limits.shp"));
	//Grid grid(bb, 1000, 0);
	//grid.resetValue("ca11");
	//grid.toBoundaryShape("E:/Vulcan/support_data/US_Boundaries/Vulcan/Contiguous/bound.shp", "E:/Vulcan/support_data/US_Boundaries/Vulcan/Contiguous_US_NOAA_Maritime_Limits.shp");
	//GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	//ds->ncols = grid.ncols;
	//ds->nrows = grid.nrows;
	//ds->numbands = 1;
	//ds->setGeoTransform(grid._adfGeoTransform);
	//ds->bound = grid.bound;
	//ds->create("E:/Vulcan/support_data/US_Boundaries/Vulcan/Contiguous/grid.tif");

	//OGREnvelope bb = BoundManager::readBoundFromRaster("C:/VulcanGridding/VulcanGridAlaska/grid.tif");
	//bb.Merge(BoundManager::readBoundFromShape("E:/Vulcan/support_data/US_Boundaries/Vulcan/Alaska_NOAA_Maritime_Limits.shp"));
	//Grid grid(bb, 1000, 0);
	//grid.resetValue("ca11");
	//grid.toBoundaryShape("E:/Vulcan/support_data/US_Boundaries/Vulcan/Alaska/bound.shp", "E:/Vulcan/support_data/US_Boundaries/Vulcan/Alaska_NOAA_Maritime_Limits.shp");
	//GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	//ds->ncols = grid.ncols;
	//ds->nrows = grid.nrows;
	//ds->numbands = 1;
	//ds->setGeoTransform(grid._adfGeoTransform);
	//ds->bound = grid.bound;
	//ds->create("E:/Vulcan/support_data/US_Boundaries/Vulcan/Alaska/grid.tif");

	//OGREnvelope bb = BoundManager::readBoundFromRaster("H:/VulcanGridding/SphericalGridContiguousUS/grid.tif");
	//bb.Merge(BoundManager::readBoundFromShape("E:/Vulcan/support_data/US_Boundaries/Spherical/Contiguous_US_NOAA_Maritime_Limits.shp"));
	//
	//bb.MinX = -125.1;
	//bb.MaxX = -66.9;
	//bb.MinY = 24.2;
	//bb.MaxY = 49.5;

	//Grid grid(bb, 0.01, 0);
	//grid.resetValue("ca11");
	//grid.toBoundaryShape("E:/Vulcan/support_data/US_Boundaries/Spherical/Contiguous/bound.shp", "E:/Vulcan/support_data/US_Boundaries/Spherical/Contiguous_US_NOAA_Maritime_Limits.shp");
	//GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	//ds->ncols = grid.ncols;
	//ds->nrows = grid.nrows;
	//ds->numbands = 1;
	//ds->setGeoTransform(grid._adfGeoTransform);
	//ds->bound = grid.bound;
	//ds->create("E:/Vulcan/support_data/US_Boundaries/Spherical/Contiguous/grid.tif");
	//delete ds;
	OGREnvelope bb = BoundManager::readBoundFromRaster("H:/VulcanGridding/SphericalGridAlaska/grid.tif");
	bb.Merge(BoundManager::readBoundFromShape("E:/Vulcan/support_data/US_Boundaries/Spherical/Alaska_NOAA_Maritime_Limits.shp"));
	bb.MinX = -179.5;
	bb.MaxX = -129.9;
	bb.MinY = 51.0000;
	bb.MaxY = 71.6000;
	Grid grid(bb, 0.1, 0);
	grid.bound = bb;
	grid.resetValue("ca11");
	grid.toBoundaryShape("E:/Vulcan/support_data/US_Boundaries/Spherical/Alaska/bound.shp", "E:/Vulcan/support_data/US_Boundaries/Spherical/Alaska_NOAA_Maritime_Limits.shp");
	GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	ds->ncols = grid.ncols;
	ds->nrows = grid.nrows;
	ds->numbands = 1;
	ds->setGeoTransform(grid._adfGeoTransform);
	ds->bound = grid.bound;
	ds->create("E:/Vulcan/support_data/US_Boundaries/Spherical/Alaska/grid10km.tif");
	delete ds;
}

void OnroadCenterPoints()
{
	std::string indir = "H:/Vulcan_2014/output_data/onroad/spatial/old/";
	std::vector<std::string> fields;
	ShapeFile::copyDirDropGeometry("H:/Vulcan_2014/output_data/onroad/spatial/old/", "H:/Vulcan_2014/support_data/onroad/jianming/", fields);
}

void bin2txt(std::string name)
{
	std::vector<double> fracs;
	fracs.resize(8760);
	std::string outdir = "H:/Vulcan_2014/support_data/onroad/jianming/time/";
	std::ifstream ifs;
	ifs.open(outdir + name + ".bin", std::ios::binary);
	ifs.read((char*)(&fracs[0]), sizeof(double)* 8760);
	ifs.close();

	std::ofstream ofs;
	ofs.open(outdir + name + ".csv");

	for (size_t i = 0; i < 8760; i++)
	{
		ofs << fracs[i] << std::endl;
	}
	ofs.close();
}

void calTimeProfileCoords()
{
	std::string shapedir = "H:/Vulcan_2014/support_data/onroad/jianming/";
	std::map<int, TimeProfileCoords> timeprofiles;
	ShapeFile shp(shapedir + "CCS_Locations_prj.shp");
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	int idfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("TemporalID");

	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
	
		int pid = fea->GetFieldAsInteger(idfield);
		if (pid  > -1) {
			TimeProfileCoords coords;
			coords.X = fea->GetFieldAsDouble("X");
			coords.Y = fea->GetFieldAsDouble("Y");
			std::stringstream ss;
			ss << "5" << pid;
			timeprofiles[atoi(ss.str().data())] = coords;
		}

		OGRFeature::DestroyFeature(fea);
	}
	shp.close();



	BoundManager manager;
	
	std::vector<std::string> shapenames;
	shapenames.push_back(shapedir + "OnroadRuralLocal.shp");
	shapenames.push_back(shapedir + "OnroadRuralNonLocal.shp");
	shapenames.push_back(shapedir + "OnroadUrbanLocal.shp");
	shapenames.push_back(shapedir + "OnroadUrbanNonLocal.shp");
	OGREnvelope bound = manager.readBoundFromShapes(shapenames);
	double resol = 10 * 1000;
	bound.MinX = bound.MinX - resol;
	bound.MinY = bound.MinY - resol;
	bound.MaxX = bound.MaxX + resol;
	bound.MaxY = bound.MaxY + resol;
	

	int ncol = (int)((bound.MaxX - bound.MinX) / resol) + 1;
	int nrow = (int)((bound.MaxY - bound.MinY) / resol) + 1;

	for (size_t row = 0; row < nrow; row++)
	{
		for (size_t col = 0; col < ncol; col++)
		{
			int gridID = col + row * ncol + 1;
			double cellx = bound.MinX + resol * (col + 0.5);
			double celly = bound.MaxY - resol * (row + 0.5);
			TimeProfileCoords coords;
			coords.X = cellx;
			coords.Y = celly;
			for (int cls = 1; cls <= 4; cls++)
			{
				std::stringstream ss;
				ss << cls << gridID;
				timeprofiles[atoi(ss.str().data())] = coords;
			}
		}
	}

	std::ofstream ofs;
	ofs.open("H:/Vulcan_2014/support_data/onroad/jianming/TimeProfileCoords.csv");
	ofs << "Index,X,Y" << std::endl;
	QDir input_dir((shapedir + "/time").data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	std::vector<QFileInfo> binfiles;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		if (!fileInfo.fileName().endsWith(".bin"))
			continue;
		std::string filename = fileInfo.absoluteFilePath().toLocal8Bit().data();
		binfiles.push_back(list.at(i));
		int id = atoi(QFileInfo(filename.data()).baseName().toLocal8Bit().data());
		TimeProfileCoords coords = timeprofiles[id];
		ofs << std::fixed << id << "," << coords.X << "," << coords.Y << std::endl;
	}
	
	ofs.close();
}


void scaleGridResol(std::string infile, std::string outfile,double scalingFac)
{
	GDAL_DS<unsigned char>* ingrid = new GDAL_DS<unsigned char>();
	ingrid->open(infile,GDALAccess::GA_ReadOnly);
	ingrid->adfGeoTransform[1] = ingrid->adfGeoTransform[1] *scalingFac;
	ingrid->adfGeoTransform[5] = ingrid->adfGeoTransform[5] *scalingFac;;
	double resol = ingrid->adfGeoTransform[1];
	ingrid->ncols = 1;
	while ((ingrid->bound.MinX + ingrid->ncols * resol) < ingrid->bound.MaxX)
	{
		ingrid->ncols++;
	}

	ingrid->nrows = 1;
	while ((ingrid->bound.MaxY - ingrid->nrows * resol) > ingrid->bound.MinY)
	{
		ingrid->nrows++;
	}
	ingrid->close();
	ingrid->create(outfile);
	delete ingrid;
}


void fixOnroad()
{
	ShapeFile shp("H:/Vulcan_2014/output_data/onroad/spatial/OnroadRuralLocal.shp");
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	int idxfield = shp.getOrCreateField("idx", OGRFieldType::OFTString);
	int id = 0;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		printf("%d\n", id);
		std::string idx = fea->GetFieldAsString(idxfield);
		OGRFeature::DestroyFeature(fea);
		id++;
	}
	shp.close();
}

void grid_LA_For_ODIAC()
{

	std::string outdir = "C:/Hestia_ODIAC/LA/";
	std::string shapesindir = "B:/LA_Version2/Vulcan_output/spatial/";
	std::string fishnetraster = "B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/LA.tif";
	QDir(outdir.data()).mkdir(".");
	Utils::updateFootprintForDir(shapesindir, true);
	Preprocessor::gridFolderByRaster(shapesindir, outdir, fishnetraster);
	//std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"timestruct" }, 6);
	//ShapeFile::copyDirDropGeometry(intersected, indir, fields2keep);

}
void grid_LA_Onroad_For_ODIAC()
{

	std::string outdir = "C:/Hestia_ODIAC/LA_Onroad/";
	std::string shapesindir = "B:/LA_Version2/gridPrep_SHP_master/Onroad/";
	std::string fishnetraster = "B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/LA.tif";
	QDir(outdir.data()).mkdir(".");
	//Utils::updateFootprintForDir(shapesindir, true);
	Preprocessor::gridFolderByRaster(shapesindir, outdir, fishnetraster);
	//std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"timestruct" }, 6);
	//ShapeFile::copyDirDropGeometry(intersected, indir, fields2keep);
}
void grid_For_ODIAC(std::string shapesindir, std::string outdir, std::string fishnetraster)
{

	//std::string outdir = "C:/Hestia_ODIAC/LA_Onroad/";
	//std::string shapesindir = "B:/LA_Version2/gridPrep_SHP_master/Onroad/";
	//std::string fishnetraster = "B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/LA.tif";
	QDir(outdir.data()).mkdir(".");
	//Utils::updateFootprintForDir(shapesindir, true);
	//Preprocessor::gridFolderByRaster(shapesindir, outdir, fishnetraster);

	Grid grid;
	grid.fromFishnetRaster(fishnetraster);
	grid.reset();
	grid.gatherCells(outdir + "Onroad.shp");
	grid.toRaster(shapesindir + "Onroad.tif", grid.proj);

	//std::vector<std::string> fields2keep;// = Utils::buildVector("", new std::string[6]{ "ca10","ca11","ca12","ca13" ,"ca14" ,"timestruct" }, 6);
	//ShapeFile::copyDirDropGeometry(intersected, indir, fields2keep);
}



void mergeTiles()
{
	GDAL_DS<unsigned char>* destDS = new GDAL_DS<unsigned char>();
	destDS->open("H:/BaiduYunDownload/Guangzhou3D/Production_OSGB/Guangzhou.tif");
	printf("%d\n", destDS->nrows);
	std::string indir = "H:/BaiduYunDownload/Guangzhou3D/Production_OSGB/Data/DOM/";
	std::vector<unsigned char*> pixelData;
	for (size_t band = 1; band <= 4; band++)
	{
		unsigned char* data = new unsigned char[destDS->slice];
		pixelData.push_back(data);
	}

	for (size_t tilerow = 0; tilerow <= 17; tilerow++)
	{
		for (size_t tilecol = 0; tilecol <= 25; tilecol++)
		{
			std::stringstream ss;
			ss << indir << tilerow << "_" << tilecol << ".tif";
			if (!QFileInfo(ss.str().data()).exists())
				continue;
			printf("%s\n", ss.str().data());
			GDAL_DS<unsigned char>* srcDS = new GDAL_DS<unsigned char>();
			srcDS->open(ss.str().data());
	
			for (size_t band = 1; band <= 4; band++)
			{
				unsigned char* descdata = pixelData[band-1];
				unsigned char* srcdata = srcDS->readData(band);
				unsigned char* psrcdata = srcdata;
				for (size_t row = 0; row < srcDS->nrows; row++)
				{
					size_t startpos = (tilerow*srcDS->nrows + row) * destDS->ncols + tilecol * srcDS->ncols;
					memcpy(descdata + startpos, psrcdata, srcDS->ncols);
					psrcdata += srcDS->ncols;
				}
				delete srcdata;
			}

			delete srcDS;
		}
	}
	destDS->close();
	destDS->create("H:/BaiduYunDownload/Guangzhou3D/Production_OSGB/Guangzhou2.tif");
	for (size_t band = 1; band <= 4; band++)
	{
		destDS->writeData(band, pixelData[band-1]);
		delete[] pixelData[band-1];
	}
	delete destDS;
	
}

void makeAlphaChannel(std::string infile, std::string outfile)
{

	GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	ds->open(infile);
	unsigned char* r = ds->readData(1);
	unsigned char* g = ds->readData(2);
	unsigned char* b = ds->readData(3);
	unsigned char* a = new unsigned char[ds->slice];
	for (size_t i = 0; i < ds->slice; i++)
	{
		unsigned short val = r[i]+g[i]+b[i];
		if (val == 0)
			a[i] = 0;
		else
			a[i] = 255;
	}
	ds->close();
	ds->numbands = 4;
	ds->create(outfile);
	ds->writeData(1,r);ds->writeData(2,g);ds->writeData(3,b);ds->writeData(4,a);
	ds->close();
	delete[] r; delete[] g; delete[] b; delete[] a;
	delete ds;
}

int timeshift(int srcyear, int destyear)
{
	QDate srcDate(srcyear, 1, 1);
	QDate destDate(destyear, 1, 1);
	int lagdays = 0;
	int hoursinweek = 24 * 7;
	double* lastweek = new double[hoursinweek];

	while (destDate.dayOfWeek() != srcDate.dayOfWeek())
	{
		srcDate = srcDate.addDays(1);
		lagdays++;
	}
	int laghours = lagdays * 24;
	return laghours;
}

#include "RasterizationGridder.h"
void integrateGridToShapes(std::string shapefile, std::string rasterfile, std::string fieldname)
{
	Grid grid;
	grid.fromFishnetRaster(rasterfile, true);
	printf("sumofgridcells = %f\n", grid.sum());
	Grid* pGrid = &grid;
	OGREnvelope bound = pGrid->bound;
	ShapeFile input(shapefile.data(), 1);
	int findex = input.getOrCreateField(fieldname.data(), OGRFieldType::OFTReal);
	OGRFeature* poFeature;
	int fid = -1;
	TileRenderer renderer;
	double cellsize = grid._adfGeoTransform[1];
	double cellsize_scale = 1.0 / 20;
	double RASTERIZATION_RESOLUTION = cellsize * cellsize_scale;
	double cellval_scale = cellsize_scale * cellsize_scale;
	RasterizationGridder rasterizer;
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{
		fid++;
		//if (fid != 62141) {
		//	OGRFeature::DestroyFeature(poFeature);
		//	continue;
		//}
		OGRGeometry* poGeometry = poFeature->GetGeometryRef();
		if (!poGeometry)
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		OGRwkbGeometryType gtype = poGeometry->getGeometryType();
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
		double integral = 0;
		for (int igeom = 0; igeom < geomlist.size(); igeom++)
		{
			OGRGeometry* poGeometry = geomlist[igeom];
			OGREnvelope geobb;
			poGeometry->getEnvelope(&geobb);
			if (!bound.Intersects(geobb))
				continue;
			std::vector<unsigned long> feamap;
			renderer.init(RASTERIZATION_RESOLUTION, geobb);
			for (int itile = 0; itile < renderer.numtiles; itile++)
			{
				RasterMask* mask = renderer.drawTile(itile, poGeometry, true);
				rasterizer.gatherCellsFromPolygonMap(mask, pGrid, feamap);
				delete mask;
			}
			double scalefac = renderer.cellsize / cellsize;
			scalefac = scalefac * scalefac;
			for (int n = 0; n < feamap.size(); n++)
			{
				double gridval = pGrid->cells[feamap[n]];
				integral += gridval*scalefac;
			}
		}

		poFeature->SetField(findex, integral);
		input.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
		if (fid % 100 == 0)
			printf("%d\n", fid);
	}
	printf("number of features = %d\n", fid);

}


void subsetShapeByFIPS(QString infile, QString outfile,int fips) {
	QFileInfo fileInfo(infile);
	if (!fileInfo.fileName().endsWith(".shp"))
		return;
	ShapeFile shpin(fileInfo.absoluteFilePath().toLocal8Bit().data());
	int fipsID = shpin.getOrCreateField("FIPS", OGRFieldType::OFTInteger);
	ShapeFile shpout;
	shpout.create(outfile.toLocal8Bit().data(), shpin.poLayer->GetSpatialRef(), shpin.poLayer->GetLayerDefn(), shpin.poLayer->GetGeomType());
	OGRFeature* fea;
	while ((fea = shpin.poLayer->GetNextFeature()) != NULL)
	{
		if (fea->GetFieldAsInteger(fipsID) == fips) {
			shpout.poLayer->CreateFeature(fea);
		}
		OGRFeature::DestroyFeature(fea);
	}
	
}
void subsetShapesByFIPS(QString indir, QString outdir, int fips) {
	QDir input_dir(indir);
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	indir = input_dir.absolutePath() + "/";
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		subsetShapeByFIPS(fileInfo.absoluteFilePath(), outdir + fileInfo.fileName(), fips);
	}
}


void makeMaskChannel(std::string infile, std::string outfile)
{

	GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	ds->open(infile);
	unsigned char* r = ds->readData(1);
	unsigned char* g = ds->readData(2);
	unsigned char* b = ds->readData(3);
	unsigned char* a = new unsigned char[ds->slice];
	for (size_t i = 0; i < ds->slice; i++)
	{
		unsigned short val = r[i] + g[i] + b[i];
		if (r[i] == 255 && g[i] == 0 && b[i]==0)
			a[i] = 255;
		else
			a[i] = 0;
	}
	ds->close();
	ds->numbands = 1;
	ds->create(outfile);
	ds->writeData(1, a, 0);
	ds->close();
	delete[] r; delete[] g; delete[] b; delete[] a;
	delete ds;
}
#include "ogrsf_frmts.h"
struct pos2d
{
	double x;
	double y;
	pos2d(double _x, double _y) :x(_x), y(_y) {

	}
	pos2d() {
	}
};

OGRPolygon* toOGRPolygon(OGRLayer* layer, const OGREnvelope& bb)
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

OGRPolygon* toOGRPolygon(OGRLayer* layer, std::vector<pos2d> posarr)
{

	OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
	OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);
	for (size_t i = 0; i < posarr.size(); i++)
	{
		linearRing->addPoint(posarr[i].x, posarr[i].y);
	}

	poPolygon->addRing(linearRing);//also crashed
	return poPolygon;
}
void createGridFromCenters(std::string outputfile)
{
	char* srcwktbuf = new char[10000];
	OGRSpatialReference oSRS;
	oSRS.SetWellKnownGeogCS("WGS84");
	oSRS.exportToWkt(&srcwktbuf);
	std::string wkt = srcwktbuf;
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
	int idIdx = fishnet.getOrCreateField("Id", OGRFieldType::OFTInteger);
	std::vector<int> fields;
	GDAL_DS<float>* lon_tif = new GDAL_DS<float>();
	lon_tif->open("B:/Indianapolis/GridDef/LON.tif");
	float* lon = lon_tif->readData(1);
	GDAL_DS<float>* lat_tif = new GDAL_DS<float>();
	lat_tif->open("B:/Indianapolis/GridDef/LAT.tif");
	float* lat = lat_tif->readData(1);
	int nrows = lon_tif->nrows;
	int ncols = lon_tif->ncols;
	int id = 0;
	int slice = nrows * ncols;
	int flon = fishnet.getOrCreateField("lon", OGRFieldType::OFTReal);
	int flat = fishnet.getOrCreateField("lat", OGRFieldType::OFTReal);

	std::vector<pos2d> posarray;
	std::vector<std::vector<double>> xarray;
	std::vector<std::vector<double>> yarray;
	xarray.resize((nrows + 1) * (ncols + 1));
	yarray.resize((nrows + 1) * (ncols + 1));
	posarray.resize((nrows + 1) * (ncols + 1));

	int nrows_plus1 = lon_tif->nrows + 1;
	int ncols_plus1 = lon_tif->ncols + 1;
	for (size_t i = 0; i < nrows; i++)
	{
		for (size_t j = 0; j < ncols; j++)
		{
			double centerx = lon[id];
			double centery = lat[id];
			double xresol = 0;
			double yresol = 0;
			if (i < nrows - 1) {
				yresol = centery - lat[(i + 1) * ncols + j];
			}
			else {
				yresol = lat[(i - 1) * ncols + j] - centery;
			}

			if (j < ncols - 1) {
				xresol = lon[i * ncols + (j + 1)] - centerx;
			}
			else {
				xresol = centerx - lon[i * ncols + (j-1)];
			}
			OGREnvelope bb;
			bb.MaxX = centerx + xresol * 0.5;
			bb.MaxY = centery + yresol* 0.5;
			bb.MinX = centerx - xresol * 0.5;
			bb.MinY = centery - yresol* 0.5;
			xarray[i * ncols_plus1 + (j + 1)].push_back(bb.MaxX);
			xarray[(i + 1) * ncols_plus1 + (j + 1)].push_back(bb.MaxX);
			xarray[i * ncols_plus1 + j].push_back(bb.MinX);
			xarray[(i + 1) * ncols_plus1 + j].push_back(bb.MinX);

			yarray[i * ncols_plus1 + j].push_back(bb.MaxY);
			yarray[i * ncols_plus1 + (j + 1)].push_back(bb.MaxY);

			yarray[(i + 1) * ncols_plus1 + j].push_back(bb.MinY);
			yarray[(i + 1) * ncols_plus1 + (j + 1)].push_back(bb.MinY);
			id++;
		}
	}
	id = 0;
	for (size_t i = 0; i < nrows_plus1; i++)
	{
		for (size_t j = 0; j < ncols_plus1; j++)
		{
			double sumx = 0;
			double sumy = 0;
			for (size_t n = 0; n < xarray[id].size(); n++)
			{
				sumx += xarray[id][n];
			}
			for (size_t n = 0; n < yarray[id].size(); n++)
			{
				sumy += yarray[id][n];
			}
			posarray[id].x = sumx / xarray[id].size();
			posarray[id].y = sumy / yarray[id].size();
			
			id++;
		}
	}
	id = 0;
	for (size_t i = 0; i < nrows; i++)
	{
		for (size_t j = 0; j < ncols; j++)
		{
			double centerx = lon[id];
			double centery = lat[id];
			double xresol = 0;
			double yresol = 0;
			if (i < nrows - 1) {
				yresol = centery - lat[(i + 1) * ncols + j];
			}
			else{
				yresol = lat[(i - 1) * ncols + j] - centery;
			}

			if (j < ncols - 1) {
				xresol = lon[i * ncols + (j + 1)] - centerx;
			}
			else {
				xresol = centerx - lon[i * ncols + j];
			}

			OGREnvelope bb;
			bb.MaxX = centerx + xresol * 0.5;
			bb.MaxY = centery + yresol* 0.5;
			bb.MinX = centerx - xresol * 0.5;
			bb.MinY = centery - yresol* 0.5;

			std::vector<pos2d> posarr;
			posarr.push_back(posarray[(i + 1) * ncols_plus1 + j]);
			posarr.push_back(posarray[i * ncols_plus1 + j]);
			posarr.push_back(posarray[i * ncols_plus1 + (j + 1)]);
			posarr.push_back(posarray[(i + 1) * ncols_plus1 + (j + 1)]);
			posarr.push_back(posarray[(i + 1) * ncols_plus1 + j]);

			//linearRing->addPoint(bb.MinX, bb.MinY);
			//linearRing->addPoint(bb.MinX, bb.MaxY);
			//linearRing->addPoint(bb.MaxX, bb.MaxY);
			//linearRing->addPoint(bb.MaxX, bb.MinY);
			//linearRing->addPoint(bb.MinX, bb.MinY);

			OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(fishnet.poLayer->GetLayerDefn());
			poFeaPolygon->SetField(idIdx, id);
			poFeaPolygon->SetField(flon, centerx);
			poFeaPolygon->SetField(flat, centery);
			OGRPolygon *poPolygon = toOGRPolygon(fishnet.poLayer, posarr);
			poFeaPolygon->SetGeometry(poPolygon);

			fishnet.poLayer->CreateFeature(poFeaPolygon);

			OGRFeature::DestroyFeature(poFeaPolygon);

			id++;
		}

	}
	fishnet.close();
}

void createIndyGrid() {

	OGREnvelope bound;
	bound.MinX = -86.660984;	bound.MaxX = -85.636086;	bound.MinY = 39.367186;	bound.MaxY = 40.150980;
	int nrows = 87; int ncols = 87;
	double resolX = (bound.MaxX - bound.MinX) / ncols;
	double resolY = (bound.MinY - bound.MaxY) / nrows;
	double adfGeoTransform[6];
	adfGeoTransform[0] = bound.MinX;
	adfGeoTransform[1] = resolX;
	adfGeoTransform[2] = 0;
	adfGeoTransform[3] = bound.MaxY;
	adfGeoTransform[4] = 0;
	adfGeoTransform[5] = resolY;

	GDAL_DS<char>* ds = new GDAL_DS<char>();
	ds->nrows = nrows; ds->ncols = ncols; ds->numbands = 1;
	ds->slice = ds->nrows * ds->ncols;
	ds->setGeoTransform(adfGeoTransform);
	ds->create("B:/Indianapolis/GridDef/IndyGrid.tif");
	srand(time(NULL));
	char* data = new char[ds->slice];
	for (size_t i = 0; i < ds->slice; i++)
	{
		data[i] = (char)(rand() % 255 + 1);
	}
	ds->writeData(1, data);
	ds->close();
	delete[] data;
}


void local2UTC()
{
	int years[] = { 2010,2011,2012,2013,2014,2015 };
	int hoursinyrs[]{ 8760,8760,8784,8760,8760,8760 };
	int utcoffet = 8;
	int hoursinweek = 24 * 7;
	struct HourRange {
		int year;
		int starthour;
		int endhour;
		HourRange() {

		}
		HourRange(int _year, int _starthour, int _endhour) 
		: year(_year),starthour(_starthour),endhour(_endhour){

		}
	};
	std::vector<std::vector<HourRange>> matches;
	for (int n = 0; n < 6; n++)
	{
		int year = years[n];
		std::vector<HourRange> parts;
		int hoursinyr = hoursinyrs[n];
		parts.push_back(HourRange(year,utcoffet, hoursinyr - 1));
		if (n < 5) {
			parts.push_back(HourRange(year + 1, 0, utcoffet-1));
		}
		else
		{
			parts.push_back(HourRange(year, hoursinyr - hoursinweek, hoursinyr - hoursinweek + (utcoffet - 1)));
		}
		matches.push_back(parts);
	}
	std::string cityname = "LAbasin";
	std::string version = "v2.5";
	std::string outdir = "B:/LA_Version2/Vulcan_output/LABasin_v2.5_utc/upload/";
	std::string indir = "B:/LA_Version2/Vulcan_output/GriddedEmissions/upload/";
	std::string fishnetraster = "B:/LA_Version2/gridPrep_SHP_master/fishnet.tif";

	for (int n = 0; n < 6; n++)
	{
		if (n != 2)
			continue;
		TemporalGridder gridder(outdir, hoursinyrs[n]);
		std::string hourlyFilename = outdir + cityname + ".total.hourly." + Utils::int2string(years[n]) + "." + version + ".nc";
		gridder.fromFishnetRaster(fishnetraster);
		gridder.initializeHourlyNetCDF(hourlyFilename);
		double* xcoords = new double[gridder.ncols];
		double* ycoords = new double[gridder.nrows];
		double* tcoords = new double[gridder.numhours];
		for (int i = 0; i < gridder.ncols; i++)
		{
			xcoords[i] = gridder._adfGeoTransform[0] + gridder._adfGeoTransform[1] * i + gridder._adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < gridder.nrows; i++)
		{
			ycoords[i] = gridder._adfGeoTransform[3] + gridder._adfGeoTransform[5] * i + gridder._adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < gridder.numhours; i++)
		{
			tcoords[i] = i;
		}

		GDAL_DS<double>* ds = NULL;
		int activeYear = -1;

		std::vector<HourRange> parts = matches[n];
		printf("%d: {\n", years[n]);
		int curhour = 0;
		double* slice = new double[gridder.ncols*gridder.nrows];
		double* totalgrid = new double[gridder.ncols*gridder.nrows];
		memset(totalgrid, 0, gridder.slice*sizeof(double));
		for (int m = 0; m < parts.size(); m++)
		{
			if (parts[m].year != activeYear) {
				if (ds) {
					delete ds;
				}
				std::string srchourlyFilename = indir + cityname + ".total.hourly." + Utils::int2string(parts[m].year) + "." + version + ".nc";
				ds = new GDAL_DS<double>();
				ds->open(srchourlyFilename);
			}
			for (int ihour = parts[m].starthour; ihour <= parts[m].endhour; ihour++)
			{
				printf("%d,%d,%d\n", parts[m].year, ihour, curhour);
				ds->readData(ihour+1, slice);
				gridder.addGridCells(slice, totalgrid, gridder.slice);
				gridder.hourlyNCFile.writeSlice(curhour, slice);
				curhour++;
			}
			printf("%d,%d,%d\n", parts[m].year, parts[m].starthour, parts[m].endhour);
		}
		printf("}\n");
		if (ds) {
			delete ds;
		}

		gridder.hourlyNCFile.write(0, xcoords);
		gridder.hourlyNCFile.write(1, ycoords);
		gridder.hourlyNCFile.write(2, tcoords);
		gridder.hourlyNCFile.close();


		TemporalGridder annualGridder(outdir, hoursinyrs[n]);
		std::string annualFilename = outdir + cityname + ".total.annual." + Utils::int2string(years[n]) + "." + version + ".nc";
		annualGridder.fromFishnetRaster(fishnetraster);
		annualGridder.initializeAnnualNetCDF(annualFilename);
		double total = 0;
		for (int i = 0; i < gridder.slice; i++)
		{
			total = total + totalgrid[i];
		}
		printf("total=%f\n", total);
		annualGridder.totalNCFile.writeSlice(0, totalgrid);
		annualGridder.totalNCFile.write(0, xcoords);
		annualGridder.totalNCFile.write(1, ycoords);
		annualGridder.totalNCFile.write(2, tcoords);
		annualGridder.totalNCFile.close();

		delete[] totalgrid;
		delete[] slice;
		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;

	}

}

void local2UTCTotals()
{
	int years[] = { 2010,2011,2012,2013,2014,2015 };
	int hoursinyrs[]{ 8760,8760,8784,8760,8760,8760 };
	int utcoffet = 8;
	int hoursinweek = 24 * 6;
	struct HourRange {
		int year;
		int starthour;
		int endhour;
		HourRange() {

		}
		HourRange(int _year, int _starthour, int _endhour)
			: year(_year), starthour(_starthour), endhour(_endhour) {

		}
	};
	std::vector<std::vector<HourRange>> matches;
	for (int n = 0; n < 6; n++)
	{
		int year = years[n];
		std::vector<HourRange> parts;
		int hoursinyr = hoursinyrs[n];
		parts.push_back(HourRange(year, utcoffet, hoursinyr - 1));
		if (n < 5) {
			parts.push_back(HourRange(year + 1, 0, 7));
		}
		else
		{
			parts.push_back(HourRange(year, hoursinyr - hoursinweek, hoursinyr + 7 - hoursinweek));
		}
		matches.push_back(parts);
	}
	std::string cityname = "LAbasin";
	std::string version = "v2.5";
	std::string outdir = "H:/LABasin_v2.5/";
	std::string indir = "B:/LA_Version2/Vulcan_output/GriddedEmissions/upload/";
	std::string fishnetraster = "B:/LA_Version2/gridPrep_SHP_master/fishnet.tif";

	for (int n = 0; n < 6; n++)
	{
		if (n != 2)
			continue;
		TemporalGridder gridder(outdir, hoursinyrs[n]);
		std::string hourlyFilename = indir + cityname + ".total.hourly." + Utils::int2string(years[n]) + "." + version + ".nc";
		std::string annualFilename = outdir + cityname + ".total.annual." + Utils::int2string(years[n]) + "." + version + ".nc";
		gridder.fromFishnetRaster(fishnetraster);
		gridder.initializeAnnualNetCDF(annualFilename);
		double* xcoords = new double[gridder.ncols];
		double* ycoords = new double[gridder.nrows];
		double* tcoords = new double[gridder.numhours];
		for (int i = 0; i < gridder.ncols; i++)
		{
			xcoords[i] = gridder._adfGeoTransform[0] + gridder._adfGeoTransform[1] * i + gridder._adfGeoTransform[1] * 0.5;
		}
		for (int i = 0; i < gridder.nrows; i++)
		{
			ycoords[i] = gridder._adfGeoTransform[3] + gridder._adfGeoTransform[5] * i + gridder._adfGeoTransform[5] * 0.5;
		}
		for (int i = 0; i < gridder.numhours; i++)
		{
			tcoords[i] = i;
		}

		GDAL_DS<double>* ds = new GDAL_DS<double>();
		ds->open(indir + cityname + ".total.annual." + Utils::int2string(years[n]) + "." + version + ".nc");
		double* totalgrid = ds->readData(1);
		double total = 0;
		total = 0;
		for (int i = 0; i < gridder.slice; i++)
		{
			total = total + totalgrid[i];
		}
		printf("total=%f\n", total);
		delete ds;
		std::vector<HourRange> parts = matches[n];
		printf("%d: {\n", years[n]);
		HourRange lastpart = parts[1];
		ds = new GDAL_DS<double>();
		ds->open(indir + cityname + ".total.hourly." + Utils::int2string(lastpart.year) + "." + version + ".nc");
		for (int ihour = lastpart.starthour; ihour <= lastpart.endhour; ihour++)
		{
			double* slice = ds->readData(ihour + 1);
			gridder.addGridCells(slice, totalgrid, gridder.slice);
			delete[] slice;
		}
		delete ds;

		ds = new GDAL_DS<double>();
		ds->open(indir + cityname + ".total.hourly." + Utils::int2string(years[n]) + "." + version + ".nc");
		for (int ihour = 0; ihour < 8; ihour++)
		{
			double* slice = ds->readData(ihour +1);
			gridder.subtractGridCells(slice, totalgrid, gridder.slice);
			delete[] slice;
		}
		delete ds;
		total = 0;
		for (int i = 0; i < gridder.slice; i++)
		{
			total = total + totalgrid[i];
		}
		printf("total=%f\n", total);
		gridder.totalNCFile.writeSlice(0, totalgrid);
		gridder.totalNCFile.write(0, xcoords);
		gridder.totalNCFile.write(1, ycoords);
		gridder.totalNCFile.write(2, tcoords);
		gridder.totalNCFile.close();
		delete[] xcoords;
		delete[] ycoords;
		delete[] tcoords;
		delete[] totalgrid;
	}

}



void writeDoubleBinary(double* data, size_t num, std::string filename) {
	std::ofstream fileout(filename.data(), std::ios::binary);
	fileout.write((char*)data, num * 8);
	fileout.close();
}
void non_urban_gridcells(std::string maskfile,std::string dir = "E:/Vulcan/gridPrep_SHP_master/urban_gridding/") {
	std::string sectors[11] = { "Airport","CMV","Cement","Commercial" ,"ElecProd" ,"Industrial" ,"Nonroad" ,"Onroad" ,"Railroad" ,"Residential","Total"};
	GDAL_DS<unsigned int>* mask_ds = new GDAL_DS<unsigned int>();
	mask_ds->open(maskfile);
	unsigned int* mask = mask_ds->readData(1);
	unsigned int mask_nodata = mask_ds->getNoData(1);

	GDAL_DS<double>* total_ds = new GDAL_DS<double>();
	total_ds->open(dir + "Total.tif");
	double* total = total_ds->readData(1);
	double total_nodata = total_ds->getNoData(1);

	GDAL_DS<double>* pop_ds = new GDAL_DS<double>();
	pop_ds->open(dir + "Population.tif");
	double* pop = pop_ds->readData(1);
	double pop_nodata = pop_ds->getNoData(1);

	std::vector<unsigned int> indices;
	std::vector<double> pop_data;
	std::vector<double> urban_mask;
	QString qdirname = QString(dir.data()) + QFileInfo(maskfile.data()).baseName() + "/";
	QDir(qdirname).mkpath(".");
	std::string outdir = qdirname.toLocal8Bit().data();
	double sum_total = 0;
	for (size_t ncell = 0; ncell < pop_ds->slice; ncell++)
	{
		double pop_val = pop[ncell];
		double total_val = total[ncell];
		if (total_val != total_nodata) {
			sum_total += total_val;
		}
		unsigned int mask_val = mask[ncell];
		if (pop_val < 1 || pop_val == pop_nodata || total_val == total_nodata ) {
			continue;
		}

		if (mask_val == mask_nodata) {
				mask_val = 0;
		}
		else{
			mask_val = 1;
		}
		pop_data.push_back(pop_val);
		urban_mask.push_back(mask_val);
		indices.push_back(ncell);
	}
	printf("%f\n", sum_total);
	writeDoubleBinary(&urban_mask[0], urban_mask.size(), outdir + "UrbanMask.data");
	writeDoubleBinary(&pop_data[0], pop_data.size(), outdir + "Population.data");
	delete mask_ds;
	delete[] mask;

	delete total_ds;
	delete[] total;

	delete pop_ds;
	delete[] pop;

	for each (std::string sector in sectors)
	{
		GDAL_DS<double>* ffco2_ds = new GDAL_DS<double>();
		ffco2_ds->open(dir + sector + ".tif");
		double* ffco2 = ffco2_ds->readData(1);
		double ffco2_nodata = ffco2_ds->getNoData(1);
		std::vector<double> ffco2_data;
		for each (unsigned int idx in indices) {
			double ffco2_val = ffco2[idx];
			if (ffco2_val == ffco2_nodata) {
				ffco2_val = 0;
			}
			ffco2_data.push_back(ffco2_val);
		}
		writeDoubleBinary(&ffco2_data[0], ffco2_data.size(), outdir + sector + ".data");
		delete[] ffco2;
		delete ffco2_ds;
	}


}


void non_urban_gridcells(std::vector<std::string> maskfiles, std::string dir = "E:/Vulcan/gridPrep_SHP_master/urban_gridding/") {
	std::string sectors[16] = { "Jones_Kammen_Total","Scope2ComtC","Scope2IndtC","Scope2RestC","Scope2TotaltC", "Airport","CMV","Cement","Commercial" ,"ElecProd" ,"Industrial" ,"Nonroad" ,"Onroad" ,"Railroad" ,"Residential","Total" };


	GDAL_DS<double>* total_ds = new GDAL_DS<double>();
	total_ds->open(dir + "Total.tif");
	double* total = total_ds->readData(1);
	double total_nodata = total_ds->getNoData(1);

	GDAL_DS<double>* pop_ds = new GDAL_DS<double>();
	pop_ds->open(dir + "Population.tif");
	double* pop = pop_ds->readData(1);
	double pop_nodata = pop_ds->getNoData(1);

	std::vector<unsigned int> indices;
	std::vector<double> pop_data;

	std::string outdir = dir;
	double sum_total = 0;
	for (size_t ncell = 0; ncell < pop_ds->slice; ncell++)
	{
		double pop_val = pop[ncell];
		double total_val = total[ncell];
		if (total_val != total_nodata) {
			sum_total += total_val;
		}

		if (pop_val < 1 || pop_val == pop_nodata || total_val == total_nodata) {
			continue;
		}

		pop_data.push_back(pop_val);
		indices.push_back(ncell);
	}
	writeDoubleBinary(&pop_data[0], pop_data.size(), outdir + "Population.data");
	for each (std::string maskfile in maskfiles) {
		GDAL_DS<unsigned int>* mask_ds = new GDAL_DS<unsigned int>();
		mask_ds->open(maskfile);
		unsigned int* mask = mask_ds->readData(1);
		unsigned int mask_nodata = mask_ds->getNoData(1);
		std::vector<double> urban_mask;
		for each (unsigned int idx in indices) {
			unsigned int mask_val = mask[idx];
			if (mask_val == mask_nodata) {
				mask_val = 0;
			}
			else {
				mask_val = 1;
			}
			urban_mask.push_back(mask_val);
		}
		delete mask_ds;
		delete[] mask;
		printf("%f\n", sum_total);
		std::string maskname = QFileInfo(maskfile.data()).baseName().toLocal8Bit().data();
		writeDoubleBinary(&urban_mask[0], urban_mask.size(), outdir + maskname + ".data");
	}

	delete total_ds;
	delete[] total;

	delete pop_ds;
	delete[] pop;

	for each (std::string sector in sectors)
	{
		GDAL_DS<double>* ffco2_ds = new GDAL_DS<double>();
		ffco2_ds->open(dir + sector + ".tif");
		double* ffco2 = ffco2_ds->readData(1);
		double ffco2_nodata = ffco2_ds->getNoData(1);
		std::vector<double> ffco2_data;
		for each (unsigned int idx in indices) {
			double ffco2_val = ffco2[idx];
			if (ffco2_val == ffco2_nodata) {
				ffco2_val = 0;
			}
			ffco2_data.push_back(ffco2_val);
		}
		writeDoubleBinary(&ffco2_data[0], ffco2_data.size(), outdir + sector + ".data");
		delete[] ffco2;
		delete ffco2_ds;
	}


}
int main(int argc, char** argv)
{
	OGRRegisterAll();
	GDALAllRegister();
	//Utils::updateFootprint("B:/LA_Version2/Vulcan/shapefiles/gz_2010_us_050_00_5m.shp");
	Preprocessor::gridFolderByRaster("C:/WRF/shapefiles/Hestia/", "C:/WRF/shapefile_grid_intersection_fractions/Hestia/geo_em_d01/", "C:/WRF/geo_em_d01.tif");
	Preprocessor::gridFolderByRaster("C:/WRF/shapefiles/Hestia/", "C:/WRF/shapefile_grid_intersection_fractions/Hestia/geo_em_d02/", "C:/WRF/geo_em_d02.tif");
	Preprocessor::gridFolderByRaster("C:/WRF/shapefiles/Hestia/", "C:/WRF/shapefile_grid_intersection_fractions/Hestia/geo_em_d03/", "C:/WRF/geo_em_d03.tif");

	Preprocessor::gridFolderByRaster("C:/WRF/shapefiles/Vulcan/", "C:/WRF/shapefile_grid_intersection_fractions/Vulcan/geo_em_d01/", "C:/WRF/geo_em_d01.tif");
	Preprocessor::gridFolderByRaster("C:/WRF/shapefiles/Vulcan/", "C:/WRF/shapefile_grid_intersection_fractions/Vulcan/geo_em_d02/", "C:/WRF/geo_em_d02.tif");
	Preprocessor::gridFolderByRaster("C:/WRF/shapefiles/Vulcan/", "C:/WRF/shapefile_grid_intersection_fractions/Vulcan/geo_em_d03/", "C:/WRF/geo_em_d03.tif");
	return 0;
	//Preprocessor::gridFolderByRaster("D:/LA_Version2/WRF/shapefiles/Hestia/", "B:/LA_Version2/Vulcan/shapefile_grid_intersection_fractions/", "B:/LA_Version2/Vulcan/fishnet.tif");
	//Preprocessor::gridFolderByRaster("B:/LA_Version2/Vulcan/shapefiles/", "B:/LA_Version2/Vulcan/shapefile_grid_intersection_fractions/", "B:/LA_Version2/Vulcan/fishnet.tif");
	//Grid grid;
	//grid.fromFishnetRaster("B:/basegrid/grid.tif");
	//grid.toShape(grid.proj, "B:/basegrid/grid.shp",false);
	//grid_For_ODIAC("B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/LA/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/LA/Output/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/LA/fishnet.tif");
	//grid_For_ODIAC("B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/Indianapolis/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/Indianapolis/Output/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/Indianapolis/fishnet.tif");
	//grid_For_ODIAC("B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/SaltLake/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/SaltLake/Output/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/SaltLake/fishnet.tif");
	//grid_For_ODIAC("B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/Baltimore/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/Baltimore/Output/", "B:/Hestia_FFDAS_ODIAC/Comparison/hestia/shapes/Baltimore/fishnet.tif");
	//std::vector<std::string> maskfiles;
	//maskfiles.push_back("H:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif");
	//maskfiles.push_back("H:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif");
	//maskfiles.push_back("H:/VulcanGridding/urban_outlines/Yuyu.tif");
	//maskfiles.push_back("H:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.tif");
	//maskfiles.push_back("H:/VulcanGridding/urban_outlines/GRUMP.tif");
	//maskfiles.push_back("H:/VulcanGridding/urban_outlines/CBSA.tif");
	//non_urban_gridcells(maskfiles);

	/*non_urban_gridcells("H:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif");
	non_urban_gridcells("H:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif");
	non_urban_gridcells("H:/VulcanGridding/urban_outlines/Yuyu.tif");
	non_urban_gridcells("H:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.tif");
	non_urban_gridcells("H:/VulcanGridding/urban_outlines/GRUMP.tif");
	non_urban_gridcells("H:/VulcanGridding/urban_outlines/CBSA.tif");*/

	//Grid grid;
	//grid.fromFishnetRaster("H:/Vulcan_2014/support_data/gridding/LambertContiguous1km.tif", false);
	//grid.toShape(grid.proj, "H:/Vulcan_2014/support_data/gridding/LambertContiguous1km.shp", false);
	//gridIndy_Thomas();
	//gridMarion();
	//TimestructTool::txt2binary(8760, "B:/Indianapolis/Indy.2011.Timestructurefiles/", "B:/Indianapolis/Time/2011.bin", ".txt");
	//TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2010.bin", 2011, 2010);
	//TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2012.bin", 2011, 2012);
	//TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2013.bin", 2011, 2013);
	//TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2014.bin", 2011, 2014);
	//TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2015.bin", 2011, 2015);
	//TimestructTool::binary2Text("B:/Indianapolis/Time/2012.bin", "B:/Indianapolis/Time/2012/", ".csv");

	//checkTotalsInDir("H:/VulcanGridding/urban_outlines/Yuyu/", ".tif", 1, "KgC/year", "totals.tif.csv");
	return 0;
	//Grid grid;
	//grid.fromFishnetRaster("B:/Marion/fishnet.tif");
	//ShapeFile shp("B:/Marion/fishnet.shp");
	//char buf[255];
	//char* pbuf = buf;
	//shp.poLayer->GetSpatialRef()->exportToProj4(&pbuf);
	//std::string str = pbuf;
	//printf("%s\n", pbuf);
	//printf("%s\n", pbuf);
	//checkTotalsInDir("B:/Marion/", ".nc", 1, "KgC/year", "totals.nc.csv");

	//Preprocessor::gridFolderByShape("E:/Vulcan/gridPrep_SHP_master/spatial/", "B:/Marion/VulcanSpatial/", "B:/Marion/fishnet.shp");
	//gridMarion();
	//checkTotalsInDir("B:/Indianapolis/Marion/Vulcan_output/GriddedEmissions/", ".tif", 1, "KgC/year", "totals.tif.csv");
	//checkTotalsInDir("B:/Indianapolis/Marion/Vulcan_output/GriddedEmissions/", ".nc", 1, "KgC/year", "totals.nc.csv");
	//gridIndy_Thomas();
	//TimestructTool::txt2binary(8760, "B:/Indianapolis/Time/2011/", "B:/Indianapolis/Time/2011.bin", ".txt");
	/*TimestructTool::binary2Text("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2011/", ".csv");

	TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2010.bin", 2011, 2010);
	TimestructTool::binary2Text("B:/Indianapolis/Time/2010.bin", "B:/Indianapolis/Time/2010/", ".csv");

	TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2012.bin", 2011, 2012);
	TimestructTool::binary2Text("B:/Indianapolis/Time/2012.bin", "B:/Indianapolis/Time/2012/", ".csv");

	TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2013.bin", 2011, 2013);
	TimestructTool::binary2Text("B:/Indianapolis/Time/2013.bin", "B:/Indianapolis/Time/2013/", ".csv");

	TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2014.bin", 2011, 2014);
	TimestructTool::binary2Text("B:/Indianapolis/Time/2014.bin", "B:/Indianapolis/Time/2014/", ".csv");

	TimestructTool::timeshiftbinary("B:/Indianapolis/Time/2011.bin", "B:/Indianapolis/Time/2015.bin", 2011, 2015);
	TimestructTool::binary2Text("B:/Indianapolis/Time/2015.bin", "B:/Indianapolis/Time/2015/", ".csv");*/
	//Grid grid;
	//grid.fromFishnetRaster("B:/ACES_Vulcan/ACES/ACES_tC_2011_TOTAL.tif", true);
	//grid.toShape(grid.proj, "B:/ACES_Vulcan/ACES_tC_2011_TOTAL.shp",true);
	//checkTotalsInDir("B:/ACES_Vulcan/Vulcan/", ".tif", 1, "tC/year", "totals.csv");
	//createIndyGrid();
	//local2UTC();
	//local2UTCTotals();
	//checkTotalsInDir("B:/LA_Version2/Vulcan_output/LABasin_v2.5_utc/upload/", ".nc", 1.0 / 1000000.0 / 1000.0, "MtC/year", "totals.csv");
	//
	//Utils::updateFootprintForDir("E:/Vulcan/gridPrep_SHP_master/spatial/", true);
	//Preprocessor::gridFolderByShape("E:/Vulcan/gridPrep_SHP_master/spatial/", "B:/Indianapolis/GridDef/IndyGrid/", "B:/Indianapolis/GridDef/IndyGrid.shp");
	return 0;
	//OGRSpatialReference sr;
	//sr.importFromProj4("+proj=lcc +lat_1=30 +lat_2=60 +lat_0=41.8389129639 +lon_0=-77 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs");
	//char srbuf[255];
	//char* psrbuf = srbuf;
	//sr.exportToWkt(&psrbuf);
	//printf("%s\n", psrbuf);
	//std::ofstream ofs("IndyLLC.txt");
	//ofs << psrbuf;
	//ofs.close();
	//linkUrbanOutline("H:/VulcanGridding/2010/Place_2010Census_DP1.tif", "H:/VulcanGridding/2010/Place_2010Census_DP1/");
	//checkTotalsInDir("H:/VulcanGridding/2010/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");

	//linkUrbanOutline("H:/VulcanGridding/2011/Place_2010Census_DP1.tif", "H:/VulcanGridding/2011/Place_2010Census_DP1/");
	//checkTotalsInDir("H:/VulcanGridding/2011/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");

	//linkUrbanOutline("H:/VulcanGridding/2012/Place_2010Census_DP1.tif", "H:/VulcanGridding/2012/Place_2010Census_DP1/");
	//checkTotalsInDir("H:/VulcanGridding/2012/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");

	//linkUrbanOutline("H:/VulcanGridding/2013/Place_2010Census_DP1.tif", "H:/VulcanGridding/2013/Place_2010Census_DP1/");
	//checkTotalsInDir("H:/VulcanGridding/2013/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");

	//linkUrbanOutline("H:/VulcanGridding/2014/Place_2010Census_DP1.tif", "H:/VulcanGridding/2014/Place_2010Census_DP1/");
	//checkTotalsInDir("H:/VulcanGridding/2014/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");

	//linkUrbanOutline("H:/VulcanGridding/2015/Place_2010Census_DP1.tif", "H:/VulcanGridding/2015/Place_2010Census_DP1/");
	//checkTotalsInDir("H:/VulcanGridding/2015/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");



	//for (int year = 2010; year <= 2015; year++) {
	//	std::stringstream ssyear;
	//	ssyear << year;
	//	std::stringstream ss;
	//	ss << "H:/VulcanGridding/" << year << "/Place_2010Census_DP1/";
	//	printf("%s\n", ss.str().data());
	//	std::vector<std::string> files = Utils::findFiles(ss.str().data(),".bin");
	//	for each (std::string filename  in files)
	//	{
	//		if (filename.substr(filename.size() - 8, 4) == ssyear.str()) {
	//			printf("%s\n", filename.data());
	//			QFile::remove(filename.data());
	//		}
	//	}
	//}
	//gridLA_VY_Totals("2010");
	//gridLA_VY_Totals("2011");
	//gridLA_VY_Totals("2012");
	//gridLA_VY_Totals("2013");
	//gridLA_VY_Totals("2014");
	//gridLA_VY_Totals("2015");
	//

	//gridLA();
	//checkTotalsInDir("B:/LA_Version2/Vulcan_output/GriddedEmissions/", ".nc", 1.0 / 1000000.0 / 1000.0, "MtC/year", "totals.csv");
	//gridLA_VY();
	//linkUrbanOutline("H:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif", "H:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/");
	//linkUrbanOutline("H:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif", "H:/VulcanGridding/urban_outlines/Place_2010Census_DP1/");
	//linkUrbanOutline("H:/VulcanGridding/urban_outlines/Yuyu.tif", "H:/VulcanGridding/urban_outlines/Yuyu/");
	//linkUrbanOutline("H:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.tif", "H:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5/");
	//linkUrbanOutline("H:/VulcanGridding/urban_outlines/GRUMP.tif", "H:/VulcanGridding/urban_outlines/GRUMP/");
	//linkUrbanOutline("H:/VulcanGridding/urban_outlines/CBSA.tif", "H:/VulcanGridding/urban_outlines/CBSA/");
	//linkUrbanOutline("H:/VulcanGridding/2015/Place_2010Census_DP1.tif", "H:/VulcanGridding/2015/Place_2010Census_DP1/");

	//checkTotalsInDir("H:/VulcanGridding/2015/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/", ".tif", 1, "tC/year", "totals.csv");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1/", ".tif", 1, "tC/year", "totals.csv");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/Yuyu/", ".tif", 1, "tC/year", "totals.csv");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5/", ".tif", 1, "tC/year", "totals.csv");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/GRUMP/", ".tif", 1, "tC/year", "totals.csv");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/CBSA/", ".tif", 1, "tC/year", "totals.csv");
	//testTimeshift();
	//return 0;
	//TimestructTool::txt2binary(8760, "B:/LA_Version2/Vulcan_output/Time/2010/ElecProd", "B:/LA_Version2/Vulcan_output/Time/2010/ElecProd.bin");
	//TimestructTool::txt2binary(8760, "B:/LA_Version2/Vulcan_output/Time/2011/ElecProd", "B:/LA_Version2/Vulcan_output/Time/2011/ElecProd.bin");
	//TimestructTool::txt2binary(8784, "B:/LA_Version2/Vulcan_output/Time/2012/ElecProd", "B:/LA_Version2/Vulcan_output/Time/2012/ElecProd.bin");
	//TimestructTool::txt2binary(8760, "B:/LA_Version2/Vulcan_output/Time/2013/ElecProd", "B:/LA_Version2/Vulcan_output/Time/2013/ElecProd.bin");
	//TimestructTool::txt2binary(8760, "B:/LA_Version2/Vulcan_output/Time/2014/ElecProd", "B:/LA_Version2/Vulcan_output/Time/2014/ElecProd.bin");
	//TimestructTool::txt2binary(8760, "B:/LA_Version2/Vulcan_output/Time/2015/ElecProd", "B:/LA_Version2/Vulcan_output/Time/2015/ElecProd.bin");
	return 0;
	//linkUrbanOutline("B:/ACES_Vulcan/CountyBound.tif", "B:/ACES_Vulcan/ACES/");
	//checkTotalsInDir("B:/ACES_Vulcan/ACES/", ".tif", 1, "tC/year", "totals.csv");

	//linkUrbanOutline("B:/ACES_Vulcan/CountyBound.tif", "B:/ACES_Vulcan/Vulcan/");
	//checkTotalsInDir("B:/ACES_Vulcan/Vulcan/", ".tif", 1, "tC/year", "totals.csv");

	//mergeLA();

	//GDAL_DS<float>* lon_tif = new GDAL_DS<float>();
	//lon_tif->open("B:/Indianapolis/GridDef/LON.tif");
	//float* lon = lon_tif->readData(1);
	//GDAL_DS<float>* lat_tif = new GDAL_DS<float>();
	//lat_tif->open("B:/Indianapolis/GridDef/LAT.tif");
	//float* lat = lat_tif->readData(1);
	
	//createGridFromCenters("B:/Indianapolis/GridDef/IndyPoleGrid.shp");
	return 0;
	//checkTotalsInDir("B:/ACES_Vulcan/ACES/", ".tif",1.0,"tC","total.csv");
	//Grid grid;
	//grid.fromFishnetRaster("H:/BaiduYunDownload/Guangzhou3D/CO2_25cm.tif",true);
	//Grid* newgrid = grid.upscale(4);
	//newgrid->toRaster("H:/BaiduYunDownload/Guangzhou3D/CO2_1m.tif");
	//delete newgrid;

	//makeMaskChannel("H:/BaiduYunDownload/Guangzhou3D/Production_OSGB/results/Classification/classification2.v1.tif", "H:/BaiduYunDownload/Guangzhou3D/Production_OSGB/results/Classification/buildings.tif");
	//checkTotalsInDir("B:/Indianapolis/Marion/output_data/gridded/residential.2011.hourly/", ".nc",1.0,"KgC","totals.txt");


	//Utils::updateFootprintForDir("B:/Indianapolis/Marion/support_data/spatial/test/", true, 1000000000.0);
	//Preprocessor::gridFolderByRaster("B:/Indianapolis/Marion/support_data/spatial/test/",
	//	"B:/Indianapolis/Marion/support_data/gridding/wgs84_ArcGIS/",
	//	"B:/Indianapolis/Marion/support_data/gridding/grid_wgs84.tif", 1000000000.0);

	//GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	//ds->adfGeoTransform[0] = -86.695000;
	//ds->adfGeoTransform[1] = 0.002;
	//ds->adfGeoTransform[2] = 0;
	//ds->adfGeoTransform[3] = 40.379799;
	//ds->adfGeoTransform[4] = 0;
	//ds->adfGeoTransform[5] = -0.002;
	//ds->ncols= 560;
	//ds->nrows = 521;
	//ds->numbands = 1;
	//char* srcwktbuf = new char[10000];
	//OGRSpatialReference oSRS;
	//oSRS.SetWellKnownGeogCS("WGS84");
	//oSRS.exportToWkt(&srcwktbuf);
	//ds->projection = srcwktbuf;
	//ds->create("B:/Indianapolis/IndyGrid.tif");
	//delete ds;
	//adfGeoTransform[0] /* top left x */
	//adfGeoTransform[1] /* w-e pixel resolution */
	//adfGeoTransform[2] /* 0 */
	//adfGeoTransform[3] /* top left y */
	//adfGeoTransform[4] /* 0 */
	//adfGeoTransform[5] /* n-s pixel resolution (negative value) */
	//subsetShapesByFIPS("B:/Indianapolis/nonroad.shapes/","B:/Indianapolis/Marion/Vulcan_output/vulcan_spatial/",18097);
	//subsetShapeByFIPS("E:/Vulcan/gridPrep_SHP_master/spatial/OnroadRuralLocal.shp", "B:/Indianapolis/Marion/OnroadRuralLocal.shp", 18097);
	//subsetShapeByFIPS("E:/Vulcan/gridPrep_SHP_master/spatial/OnroadUrbanLocal.shp", "B:/Indianapolis/Marion/OnroadUrbanLocal.shp", 18097);
	//subsetShapeByFIPS("E:/Vulcan/gridPrep_SHP_master/spatial/OnroadRuralNonLocal.shp", "B:/Indianapolis/Marion/OnroadRuralNonLocal.shp", 18097);
	//subsetShapeByFIPS("E:/Vulcan/gridPrep_SHP_master/spatial/OnroadUrbanNonLocal.shp", "B:/Indianapolis/Marion/OnroadUrbanNonLocal.shp", 18097);



	return 0;
	//integrateGridToShapes(
	//	"H:/Vulcan_2014/urban_scaling/LandCast/us_census_bk.shp",
	//	"H:/Vulcan_2014/urban_scaling/LandCast/landscan2010.tif",
	//	"lsp2010");

	//integrateGridToShapes(
	//	"H:/Vulcan_2014/urban_scaling/LandCast/us_census_bk.shp",
	//	"H:/Vulcan_2014/urban_scaling/LandCast/LandCast2030.tif",
	//	"lsp2030");

	//integrateGridToShapes(
	//	"H:/Vulcan_2014/urban_scaling/LandCast/us_census_bk.shp",
	//	"H:/Vulcan_2014/urban_scaling/LandCast/LandCast2050.tif",
	//	"lsp2050");

	//printf("%d\n", timeshift(2011, 2010));
	//printf("%d\n", timeshift(2011, 2012));

	/*std::ifstream filein("B:/LA_Version2/Vulcan_output/time/2014/ComNonpoint.bin", std::ios::binary);
	filein.seekg(0, filein.end);
	size_t fileSize = filein.tellg();
	filein.seekg(0, filein.beg);
	int numstructs = 0;
	filein.read((char*)&numstructs, sizeof(int) * 1);
	int numhours = (fileSize - (numstructs * 20) - sizeof(int)) / (sizeof(double)*numstructs);
	char* idbuf = new char[numstructs * 20];
	memset(idbuf, 0, (size_t)numstructs * 20);
	double* fractionbuf = new double[numstructs * numhours];
	filein.read(idbuf, numstructs * 20);
	char* pstr = idbuf;
	for (size_t i = 0; i < numstructs; i++)
	{
		printf("%s\n", pstr);
		pstr += 20;
	}
	filein.read((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	filein.close();



	gridLA_VY();
	checkTotalsInDir("B:/LA_Version2/Vulcan_output/GriddedEmissions/", ".nc",0.001,"Total(tC)","nc.totals.csv");
	checkTotalsInDir("B:/LA_Version2/Vulcan_output/GriddedEmissions/", ".tif", 0.001, "Total(tC)", "tif.totals.csv");
	gridLA_VY_Totals("2010");
	gridLA_VY_Totals("2011");
	gridLA_VY_Totals("2012");
	gridLA_VY_Totals("2013");
	gridLA_VY_Totals("2014");*/
	//checkTotalsInDir("B:/LA_Version2/Vulcan_output/Vineet/GriddedEmissions/", ".nc");
	//copySCCProfiles("B:/LA_Version2/Vulcan_output/ElecProd.shp");
	//copySCCProfiles("B:/LA_Version2/Vulcan_output/RailroadPoint.shp");
	//copySCCProfiles("B:/LA_Version2/Vulcan_output/NonroadPoint.shp");
	//copySCCProfiles("B:/LA_Version2/Vulcan_output/IndPoint.shp");
	//copySCCProfiles("B:/LA_Version2/Vulcan_output/ComPoint.shp");

	//TimestructTool ts;
	//ts.binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/ElecProd.bin", "B:/LA_Version2/Vulcan_output/Time_2012/", ".csv");
	//checkTotal("H:/Vulcan_2014/output_data/gridded/d001.nc");
	//checkTotalsInDir("H:/Vulcan_2014/output_data/gridded/Vulcan1km/ContiguousUS/",".nc");
	//std::vector<std::string> fields;
	//Grid grid(BoundManager::readBoundFromShape("C:/VulcanGrids/Coconino/bound.shp"), 100);
	//grid.proj = Utils::getProjFromShapefile("C:/VulcanGrids/Coconino/bound.shp");
	//grid.toRaster("C:/VulcanGrids/Coconino100m.tif", grid.proj);
	//GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	//ds->open("C:/VulcanGrids/Coconino100m.tif");
	//ds->create("C:/VulcanGrids/Coconino100m.tif");
	//delete ds;

	//Grid grid(BoundManager::readBoundFromShape("C:/VulcanGrids/Flagstaff/bound.shp"), 100,0);
	//grid.proj = Utils::getProjFromShapefile("C:/VulcanGrids/Flagstaff/bound.shp");
	//grid.toRaster("C:/VulcanGrids/Flagstaff100m.tif", grid.proj);
	//GDAL_DS<unsigned char>* ds = new GDAL_DS<unsigned char>();
	//ds->open("C:/VulcanGrids/Flagstaff100m.tif");
	//ds->create("C:/VulcanGrids/Flagstaff100m.tif");
	//delete ds;


	//ShapeFile::copyDropGeometry("C:/GSV2SVF/Data/FinancialDistrictRoads.shp", "C:/GSV2SVF/Data/FinancialDistrictPoints.shp", fields);
	//makeAlphaChannel("C:/Tokyo/20170710_B10.tif", "C:/Tokyo/20170710_B10_alpha.tif");
	return 0;

	//checkTotalsInDir("B:/ACES_Vulcan/ACES/", ".tif");
	//checkTotalsInDir("B:/ACES_Vulcan/Vulcan/GriddedFFCO2/", ".tif");
	//scaleRaster("B:/ACES_Vulcan/ACES/", "B:/ACES_Vulcan/ACES/Scaled/", ".tif", 0.001);

	//Preprocessor::reprojectWithArcGIS("E:/Vulcan/gridPrep_SHP_master/spatial/rail/railway.shp", "E:/Vulcan/gridPrep_SHP_master/spatial/railway.shp", "E:/Vulcan/gridPrep_SHP_master/WGS84ToVulcan.py");

	//mergeTiles();
	//gridLA_For_ODIAC();
	//grid_LA_Onroad_For_ODIAC();
	//Preprocessor::reprojectDir("E:/Vulcan/gridPrep_SHP_master/spatial/", "B:/ACES_Vulcan/Vulcan/Spatial/", "E:/Vulcan/gridPrep_SHP_master/VulcanToACES.py");
	return 0;
	//GDAL_DS<int>* ds = new GDAL_DS<int>();
	//ds->open("B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/LA.tif", GA_Update);
	//ds->scaleTotal(1, 33393.25952 * 1000000);
	//delete ds;


	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/Alaska/onroad.2011.hourly.old/d365.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/Alaska/onroad.2011.hourly/d365.nc");
	//printf("");
	//fixOnroad();
	//scaleGridResol("H:/Vulcan_2014/support_data/gridding/LambertContiguous.tif", "H:/Vulcan_2014/support_data/gridding/LambertContiguous10km.tif",10.0);
	//scaleGridResol("H:/Vulcan_2014/support_data/gridding/LambertAlaska.tif", "H:/Vulcan_2014/support_data/gridding/LambertAlaska10km.tif", 10.0);
	//calTimeProfileCoords();


	//std::vector<std::string> fields;
	//ShapeFile::copyDirDropGeometry("H:/Vulcan_2014/output_data/onroad/spatial/old/", "H:/Vulcan_2014/support_data/onroad/jianming/", fields);

	//RoadTimeIDW roadidw;
	//roadidw.m_maxDist = 1000 * 1000;
	//roadidw.m_minnumpoints = 1;
	//roadidw.m_maxnumpoints = 30;
	//roadidw.loadCCSStations("H:/Vulcan_2014/support_data/onroad/jianming/CCS_Locations_prj.shp", "H:/Vulcan_2014/support_data/onroad/jianming/CCS_Temporal_US.bin");
	//roadidw.createIDW();
	//roadidw.m_outdir = "H:/Vulcan_2014/support_data/onroad/jianming/time/";
	/*BoundManager manager;
	std::string shapedir = "H:/Vulcan_2014/support_data/onroad/jianming/";
	std::vector<std::string> shapenames;
	shapenames.push_back(shapedir + "OnroadRuralLocal.shp");
	shapenames.push_back(shapedir + "OnroadRuralNonLocal.shp");
	shapenames.push_back(shapedir + "OnroadUrbanLocal.shp");
	shapenames.push_back(shapedir + "OnroadUrbanNonLocal.shp");
	OGREnvelope bound = manager.readBoundFromShapes(shapenames);
	double resol = 10 * 1000;

	bound.MinX = bound.MinX - resol;
	bound.MinY = bound.MinY - resol;
	bound.MaxX = bound.MaxX + resol;
	bound.MaxY = bound.MaxY + resol;
	int ncol = (int)((bound.MaxX - bound.MinX) / resol) + 1;
	int nrow = (int)((bound.MaxY - bound.MinY) / resol) + 1;
	bound.MaxX = bound.MinX + ncol * resol;
	bound.MinY = bound.MaxY - nrow * resol;
	std::ofstream ofs;
	ofs.open("extent.txt"); 
	ofs<< std::fixed  << "nrol=" << ncol << std::endl;
	ofs << "nrow=" << nrow << std::endl;
	ofs << std::fixed << "MaxX=" << bound.MaxX << std::endl;
	ofs << std::fixed << "MinX=" << bound.MinX << std::endl;
	ofs << std::fixed << "MaxY=" << bound.MaxY << std::endl;
	ofs << std::fixed << "MinY=" << bound.MinY << std::endl;
	ofs.close();*/
	/*int ncol = (int)((bound.MaxX - bound.MinX) / resol) + 1;
	int col = (int)((-51608.047 - bound.MinX) / resol);
	int row = (int)((bound.MaxY - (-1561498.613)) / resol);
	int gridID = col + row * ncol + 1;
	//printf("%d\n", gridID);*/
	//roadidw.idw_ccs(shapenames[0], bound, resol);
	//roadidw.idw_ccs(shapenames[1], bound, resol);
	//roadidw.idw_ccs(shapenames[2], bound, resol);
	//roadidw.idw_ccs(shapenames[3], bound, resol);
	//std::ofstream ofs;
	//ofs.open("H:/Vulcan_2014/support_data/onroad/jianming/TimeProfileCoords.csv");
	//ofs << "index,X,Y" << std::endl;
	//std::map<int, TimeProfileCoords>::iterator iter = roadidw.m_TimeProfileCoords.begin();
	//while(iter != roadidw.m_TimeProfileCoords.end()){
	//	ofs << std::fixed << iter->first << "," << iter->second.X << "," << iter->second.Y << std::endl;
	//	iter++;
	//}
	//ofs.close();
	//roadidw.m_maxDist = 100000 * 100000;
	//roadidw.fillgap_ccsp(shapenames[0], "");
	//roadidw.fillgap_ccsp(shapenames[1], "");
	//roadidw.fillgap_ccsp(shapenames[2], "");
	//roadidw.fillgap_ccsp(shapenames[3], ""); 
	//TimestructTool::mergeBinary(8760, "H:/Vulcan_2014/support_data/onroad/jianming/time/", "H:/Vulcan_2014/support_data/onroad/jianming/OnroadTemporal.bin", true);
	//bin2txt("1107334");
	//bin2txt("5975");
	//bin2txt("1107335");

	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/onroad.2011.annual.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/onroad.2011.hourly/d005.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/onroad.2011.hourly/d010.nc");
	//checkTotalsInDirs("C:/VulcanOutput/WGS84_Spherical_10km/",".nc");
	/////checkTotalsInDir("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/total.2011.hourly/",".nc");
	//checkTotalsInDirs("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/", ".nc");
	//checkTotalsInDirs("C:/VulcanOutput/WGS84_Spherical_10km/Alaska/", ".nc");
	////checkTotalsInDir("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/onroad.2011.hourly/", ".nc");


	//checkTotalsInDir("C:/VulcanOutput/WGS84_Spherical_10km/Alaska/total.2011.hourly/", ".nc");
	//checkTotalsInDir("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/total.2011.hourly/", ".nc");
	////checkTotalsInDir("C:/VulcanOutput/nonroad.point.2011.hourly/",".nc");

	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/commercial.2011.annual.nc");

	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/onroad.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/elec_prod.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/cement.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/commercial.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/residential.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/nonroad.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/airport.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/cmv.2011.hourly/d001.nc");
	//checkTotal("C:/VulcanOutput/WGS84_Spherical_10km/ContiguousUS/industrial.2011.hourly/d001.nc");


	return 0;

	/*calVulcanExtents();*/
	//checkTotal("B:/com.co2.d001.nc");
	//Preprocessor::reprojectDir("E:/Vulcan/gridPrep_SHP_master/spatial/", "E:/Vulcan/gridPrep_SHP_master/spatial_wgs84/", "E:/Vulcan/gridPrep_SHP_master/VulcanToWGS84.py");
	//gridLA_VY();
	//gridLA_VY_Totals("2012");
	//gridLA_VY_Totals("2010");
	//gridLA_VY_Totals("2011");
	//gridLA_VY_Totals("2013");
	//gridLA_VY_Totals("2014");

	//load_LA_VY_To_Shapefile("B:/LA_Version2/Vulcan_output/Vineet/GriddedEmissions/WRF.VY.total.annual.2012.bin","B:/LA_Version2/Vulcan_output/Vineet/", "total2012");
	//load_LA_VY_To_Shapefile("B:/LA_Version2/Vulcan_output/Vineet/GriddedEmissions/WRF.VY.total.annual.2011.bin", "B:/LA_Version2/Vulcan_output/Vineet/", "total2011");
	//load_LA_VY_To_Shapefile("B:/LA_Version2/Vulcan_output/Vineet/GriddedEmissions/WRF.VY.onroad.annual.2012.bin", "B:/LA_Version2/Vulcan_output/Vineet/", "onroad2012");
	//load_LA_VY_To_Shapefile("B:/LA_Version2/Vulcan_output/Vineet/GriddedEmissions/WRF.VY.onroad.annual.2011.bin", "B:/LA_Version2/Vulcan_output/Vineet/", "onroad2011");
	//testTimeshift();
	//std::map<std::string, double*> timestructMap;
	//TimeStructContainer container;
	//container.load("B:/LA_Version2/Vulcan_output/Time/2012/ElecProd.bin", timestructMap);
	//checkTotalsInDir("B:/LA_Version2/Vulcan_output/GriddedEmissions/",".nc");
	//calVulcanExtents();
	///TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/Onroad.bin", "B:/LA_Version2/Vulcan_output/Time/2012/Onroad/", ".csv");
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/Nonpoint.bin", "B:/LA_Version2/Vulcan_output/Time/2012/Nonpoint/", ".csv");
	//
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/ElecProd.bin", "B:/LA_Version2/Vulcan_output/Time/2012/ElecProd2/", ".csv");
	//
	//checkTotal("B:/LA_Version2/Vulcan_output/GriddedEmissions/LAbasin.total.hourly.2012.v2.4.nc");
	//gridLA_VY();
	
	//Utils::updateFootprint("B:/LA_Version2/Vulcan_output/OnRoad.Shp", true);
	//mergeLA();
	//findNearestRoadSegment("B:/LA_Version2/gridPrep_SHP_master/Vulcan/OnroadVersions2_2.shp", "B:/LA_Version2/Vulcan_output/OnRoad.Shp");
	//testTimeshift();
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/Airport.bin", "B:/LA_Version2/Vulcan_output/Time/2012/Airport/",".csv");
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/ElecProd.bin", "B:/LA_Version2/Vulcan_output/Time/2012/ElecProd/", ".csv");
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/NonPoint.bin", "B:/LA_Version2/Vulcan_output/Time/2012/NonPoint/", ".csv");
	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2012/ComPoint.bin", "B:/LA_Version2/Vulcan_output/Time/2012/ComPoint/", ".csv");

	//TimestructTool::binary2Text("B:/LA_Version2/Vulcan_output/Time/2014/NonPoint.bin", "B:/LA_Version2/Vulcan_output/Time/2014/NonPoint/", ".csv");

	//return 0;
	//RoadTimeIDW roadidw;
	//roadidw.m_maxDist = 10000 * 3.28084;
	//roadidw.m_minnumpoints = 0;
	//roadidw.m_maxnumpoints = 1;
	//roadidw.loadStations("B:/LA_Version2/OnRoadTime/BinaryStationCounts/","E:/LA_OnRoad_Time/");
	//roadidw.coordinatesFromShapefile("B:/LA_Version2/OnRoadTime/BinaryStationCounts/stations.Shp","");
	//roadidw.createIDW();
	////roadidw.stations2shapes("B:/LA_Version2/OnRoadTime/BinaryStationCounts/stations.Shp");
	//roadidw.nearest("B:/LA_Version2/Vulcan_output/OnRoad.Shp", "timestruct", "5", "", 0);
	//roadidw.nearest("B:/LA_Version2/gridPrep_SHP_master/Vulcan/OnRoad.Shp", "timestruct", "5", "", 0);
	//Utils::updateFootprint("B:/OnroadSubset/Onroad.shp");
	//Preprocessor::intersectWithArcGIS("B:/OnroadSubset/Onroad.shp", "B:/LA_Version2/gridPrep_SHP_master/bound.shp", "B:/LA_Version2/gridPrep_SHP_master/Vulcan/Onroad.shp");
	//updateFieldAfterIntersection("B:/LA_Version2/gridPrep_SHP_master/Vulcan/Onroad.shp");

	//RoadTimeIDW roadidw;
	//roadidw.m_maxDist = 10000 * 3.28084;
	//roadidw.m_minnumpoints = 0;
	//roadidw.loadStations("B:/LA_Version2/OnRoadTime/BinaryStationCounts/","E:/LA_OnRoad_Time/");
	//roadidw.coordinatesFromShapefile("B:/LA_Version2/OnRoadTime/BinaryStationCounts/stations.Shp","");
	//roadidw.createIDW();
	////roadidw.stations2shapes("B:/LA_Version2/OnRoadTime/BinaryStationCounts/stations.Shp");
	//roadidw.idw("B:/LA_Version2/Vulcan_output/OnRoad.Shp", "timestruct", "5", 0);

	//mergeLA();

	return 0;
	//gridLA3();
	
	/*std::string indir = "H:/HestiaGridding/Los_Angeles_Vulcan/Shapes/";
	std::string outdir = "H:/HestiaGridding/Los_Angeles_Vulcan/Gridded/";
	std::string fishnetrasterfile = "H:/HestiaGridding/Los_Angeles_Vulcan/fishnet.tif";
	
	std::string outfishnetfile = outdir + "fishnet.shp";
	std::string outrasterfile = outdir + "fishnet.tif";

	std::string fishnetfile = (QFileInfo(fishnetrasterfile.data()).absoluteDir().absolutePath() + "/fishnet.shp").toLocal8Bit().data();
	std::string rasterfile = fishnetrasterfile;

	Grid fishnet;
	fishnet.fromFishnetRaster(fishnetrasterfile);
	fishnet.reset();

	
	fishnet.gatherCells("H:/HestiaGridding/Los_Angeles/output/LAbasin.onroad.2011.v2.0.shp");
	fishnet.toRaster("H:/HestiaGridding/Comparison/2.1/onroad.tif");*/


	//
	return 0;
	

	//OGREnvelope bound;
	//bound.MinX = 6116566.259531;
	//bound.MaxX = 6976304.701139;
	//bound.MinY = 1598198.757620;
	//bound.MaxY = 2020106.129355;
	//GDAL_DS<int>* ds = new GDAL_DS<int>();
	//ds->open("B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/LA.tif", GA_Update);
	//OGREnvelope bound;
	//bound.MinX = -119.478265;
	//bound.MaxX = -116.631772;
	//bound.MinY = 33.382467;
	//bound.MaxY = 34.543147;
	//GDAL_DS<int>* ds = new GDAL_DS<int>();
	//ds->open("B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/LA.tif", GA_Update);
	//ds->scaleTotal(1, 33393.25952 * 1000000);
	//delete ds;
	//ds = new GDAL_DS<int>();
	//ds->open("B:/Hestia_FFDAS_ODIAC/Comparison/ffdas_scaled/LA.tif", GA_Update);
	//ds->scaleTotal(1, 33393.25952 * 1000000);
	//delete ds;

	//checkTotalsInDir("B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/", ".tif");
	//checkTotalsInDir("B:/Hestia_FFDAS_ODIAC/Comparison/ffdas_scaled/", ".tif");
	//std::vector<std::string> files;
	//files.push_back("B:/Hestia_FFDAS_ODIAC/Comparison/hestia/LA.tif");
	//files.push_back("B:/Hestia_FFDAS_ODIAC/Comparison/odiac/LA.tif");
	//files.push_back("B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/LA.tif");
	//files.push_back("B:/Hestia_FFDAS_ODIAC/Comparison/ffdas/LA.tif");
	//files.push_back("B:/Hestia_FFDAS_ODIAC/Comparison/ffdas_scaled/LA.tif");
	//for (size_t i = 0; i < files.size(); i++)
	//{
	//	GDAL_DS<int>* ds = new GDAL_DS<int>();
	//	ds->open(files[i]);
	//	ds->extractByMask("B:/Hestia_FFDAS_ODIAC/Comparison/Grump_LA.tif",files[i]);
	//	delete ds;
	//}
	//checkTotalsInDir("B:/Hestia_FFDAS_ODIAC/Comparison/ffdas/",".tif");
	//LosAngeles LAS;
	//LAS.ChangeDataType();
	//BatchDOE batchDOE;
	//batchDOE.init("C:/doe23/eQUEST_Prototypes/","C:/doe23/",  "E:/Vulcan/time/TMY3/TMY3.dbf");
	//batchDOE.gapfill();
	//batchDOE.findNearestStations("E:/Vulcan/time/BuildingTimeCombined/fema_bsf_2002bnd_centroid.shp");
	//calculatePointCoordinates("E:/Vulcan/time/BuildingTimeCombined/fema_bsf_2002bnd_centroid.shp");
	//BatchDOE batchDOE;
	//batchDOE.init("C:/doe23/eQUEST_Prototypes/","C:/doe23/",  "E:/Vulcan/time/TMY3/TMY3.dbf");
	//batchDOE.runAll("2011");
	//batchDOE.gapfill();
	//batchDOE.output2shapes();
	//gridVulcanLAOnroad();
	//gridHestiaLAOnroad();
	/*Grid fishnetGrid;
	std::string fishnetrasterfile = "B:/LA_Version2/gridPrep_SHP_master/fishnet.tif";
	std::string fishnetshapefile = "B:/LA_Version2/gridPrep_SHP_master/fishnet.shp";
	std::string intersectedShapeFile = "B:/LA_Version2/Vulcan/gridded/Vulcan_Onroad.shp";
	std::string vulcanonroadfile = "B:/LA_Version2/Vulcan/Vulcan_Onroad.shp";
	fishnetGrid.fromFishnetRaster(fishnetrasterfile);
	fishnetGrid.ncols = fishnetGrid.ncols * 10;
	fishnetGrid.nrows = fishnetGrid.nrows * 10;
	fishnetGrid._adfGeoTransform[1] = fishnetGrid._adfGeoTransform[1] / 10;
	fishnetGrid._adfGeoTransform[5] = fishnetGrid._adfGeoTransform[5] / 10;
	fishnetGrid.reset();
	for (size_t i = 0; i < fishnetGrid.ncols*fishnetGrid.nrows; i++)
	{
		fishnetGrid.cells[i] = 1;
	}
	fishnetGrid.toRaster("B:/LA_Version2/Vulcan/fishnet100m.tif", fishnetGrid.proj);*/
	//fishnetGrid.toShape(fishnetGrid.proj,"B:/LA_Version2/Vulcan/fishnet100m.shp",false);
	////Utils::updateFootprint(vulcanonroadfile);
	////Preprocessor::intersectWithArcGIS(vulcanonroadfile, fishnetshapefile, intersectedShapeFile);
	////updateFieldAfterIntersection(intersectedShapeFile);
	////ShapeFile inshp(intersectedShapeFile);
	//ShapeFile inshp("B:/LA_Version2/LA_Basin_Bound.shp");
	//
	//fishnetGrid.gatherCells(&inshp, "ca11");
	//fishnetGrid.toRaster("B:/LA_Version2/Vulcan/gridded/LAMask.tif", fishnetGrid.proj);
	//Grid fishnetGrid;
	//std::string fishnetrasterfile = "B:/LA_Version2/gridPrep_SHP_master/fishnet.tif";
	//std::string fishnetshapefile = "B:/LA_Version2/gridPrep_SHP_master/fishnet.shp";
	//std::string intersectedShapeFile = "B:/LA_Version2/Vulcan/gridded/Hestia_Onroad.shp";
	//std::string vulcanonroadfile = "B:/LA_Version2/Vulcan/Hestia_Onroad.shp";
	//fishnetGrid.fromFishnetRaster(fishnetrasterfile);
	//fishnetGrid.reset();
	//Utils::updateFootprint(vulcanonroadfile);
	//Preprocessor::intersectWithArcGIS(vulcanonroadfile, fishnetshapefile, intersectedShapeFile);
	//updateFieldAfterIntersection(intersectedShapeFile);
	//ShapeFile inshp(intersectedShapeFile);
	//fishnetGrid.gatherCells(&inshp, "ca11");
	//fishnetGrid.toRaster("B:/LA_Version2/Vulcan/gridded/Hestia_Onroad.tif", fishnetGrid.proj);
	//BatchDOE batchDOE;
	//batchDOE.init("C:/doe23/eQUEST_Prototypes/","C:/doe23/",  "E:/Vulcan/time/TMY3/TMY3.dbf");
	//batchDOE.runAll("2011");


	//std::vector<std::string> fields;
	//fields.push_back("SEGID");
	//ShapeFile::copyDir("E:/Vulcan/output_data/onroad/", "E:/Vulcan/gridPrep_SHP_master/spatial/", fields);

	//matchInternationalAirports("E:/Vulcan/airports/nei_airports.shp", "E:/Vulcan/airports/InternationalAirports.shp");
	//remapCCA("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.tif", "C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.shp");
	//linkUrbanOutline2States("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif","C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.shp", "C:/VulcanGridding/cb_2016_us_state_500k/States.tif", "C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/StateFIPS.csv");
	//linkUrbanOutline2States("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.tif","C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5.shp", "C:/VulcanGridding/cb_2016_us_state_500k/States.tif", "C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5/StateFIPS.csv");
	//linkUrbanOutline2States("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif", "C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.shp", "C:/VulcanGridding/cb_2016_us_state_500k/States.tif", "C:/VulcanGridding/urban_outlines/Place_2010Census_DP1/StateFIPS.csv");
	//linkUrbanOutline2States("C:/VulcanGridding/urban_outlines/Yuyu.tif", "C:/VulcanGridding/urban_outlines/Yuyu.shp", "C:/VulcanGridding/cb_2016_us_state_500k/States.tif", "C:/VulcanGridding/urban_outlines/Yuyu/StateFIPS.csv");
	//linkUrbanOutline2States("C:/VulcanGridding/urban_outlines/Grump.tif", "C:/VulcanGridding/urban_outlines/Grump.shp", "C:/VulcanGridding/cb_2016_us_state_500k/States.tif", "C:/VulcanGridding/urban_outlines/Grump/StateFIPS.csv");
	//linkUrbanOutline2States("C:/VulcanGridding/urban_outlines/CBSA.tif", "C:/VulcanGridding/urban_outlines/CBSA.shp", "C:/VulcanGridding/cb_2016_us_state_500k/States.tif", "C:/VulcanGridding/urban_outlines/CBSA/StateFIPS.csv");
	//reallocateEmissionsByPop("ACS_2015_5YR_UA");
	//reallocateEmissionsByPop("GRUMP_CCA_d1500_l5");
	//reallocateEmissionsByPop("Place_2010Census_DP1");
	//reallocateEmissionsByPop("Yuyu");
	//reallocateEmissionsByPop("Grump");
	//reallocateEmissionsByPop("CBSA");

	//checkTotalsInDir("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/", ".tif");
	//return 0;

	//SubsetVulcan("B:/Baltimore/VulcanSubset/BaltimoreCityBound_Vulcan.shp", "E:/Vulcan/gridPrep_SHP_master/spatial/", "B:/Baltimore/VulcanSubset/gridPrep_SHP_master/");
	//matchAirports("E:/Vulcan/airports/nei_airports.shp", "E:/Vulcan/airports/national_airports.shp");
	//matchAirNavAirports("E:/Vulcan/airports/nei_airports.shp", "E:/Vulcan/airports/airnav_airports.shp");
	//BatchDOE batchDOE;
	//batchDOE.init("C:/doe23/eQUEST_Prototypes/","C:/doe23/",  "E:/Vulcan/time/TMY3/TMY3.dbf");
	//batchDOE.runAll("2011");
	//return 0;

	//filterWeatherStations();

	//Vulcan2014 VULCAN2;
	//VULCAN2.createAirport();
	//VULCAN2.createElecProd();
	//VULCAN2.createCementPoint();
	//VULCAN2.createNonroadPoint();
	//VULCAN2.createRailroadPoint();
	
	//Preprocessor::reprojectDir("E:/Vulcan/gridPrep_SHP_master/points_wgs84/", "E:/Vulcan/gridPrep_SHP_master/spatial/", "E:/Vulcan/gridPrep_SHP_master/WGS84ToVulcan.py");

	//std::vector<std::string> LA_DIRS;
	//LA_DIRS.push_back("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/");eateComPoint();
	//VULCAN2.createIndPoint(); LA_DIRS.push_back("B:/LA_Version2/gridPrep_SHP_master/Orange/");
	//LA_DIRS.push_back("B:/LA_Version2/gridPrep_SHP_master/Riverside/"); LA_DIRS.push_back("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/");
	//LA_DIRS.push_back("B:/LA_Version2/gridPrep_SHP_master/Ventura/");

	//createTimeStructForAirport(LA_DIRS);

	/*return 0;
	std::ifstream ifs;
	ifs.open("B:/LA_Version2/gridded_WRF_VY/WRF.VY.res.annual.2014.bin", std::ios::binary);
	double* data = new double[1826];
	ifs.read((char*)data, sizeof(double) * 1826);
	double total = 0;
	for (size_t i = 0; i < 1826; i++)
	{
		total += data[i];
	}
	ifs.close();
	total = total / 1000000;
	printf("%f,%f\n", total, total * 12 );
	printf("");*/
	//loadGrump("C:/VulcanGridding/urban_outlines/CCA/GRUMP_218/218.dat", "C:/VulcanGridding/urban_outlines/CCA/GRUMP_218/218.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/CCA/GRUMP_218/218.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/Grump/gl_grumpv1_pdens_00_grid_30/gluds00ag.tif");
	//loadCities("C:/VulcanGridding/urban_outlines/CCA/gluds90g1.tif",
	//	"C:/VulcanGridding/urban_outlines/CCA/d1500_l5/",
	//	"C:/VulcanGridding/urban_outlines/CCA/d1500_l5/d1500_l5.tif");

	//loadCities("C:/VulcanGridding/urban_outlines/CCA/gluds90g1.tif",
	//	"C:/VulcanGridding/urban_outlines/CCA/d1500_l5/",
	//	"C:/VulcanGridding/urban_outlines/CCA/GRUM_CCA_d1500_l5_2.tif");
	//loadCities("C:/VulcanGridding/urban_outlines/CCA/gluds90g1.tif",
	//	"C:/VulcanGridding/urban_outlines/CCA/d1500_l5/",
	//	"C:/VulcanGridding/urban_outlines/CCA/d1500_l5/d1500_l5.tif");
	/*GDAL_DS<unsigned int>* ffdas = new GDAL_DS<unsigned int>();
	ffdas->open("C:/VulcanGridding/urban_outlines/CCA/GRUM_CCA_d1500_l5.tif");
	OGREnvelope bound;
	bound.MinX = -124.8;
	bound.MaxX = -66.9;
	bound.MinY = 24.5;
	bound.MaxY = 49.5;
	ffdas->crop(bound, "C:/VulcanGridding/urban_outlines/GRUM_CCA_d1500_l5.tif");
	delete ffdas;*/
	//return 0;
	//mergeLA();

	//GDAL_DS<double>* ffdas = new GDAL_DS<double>();
	//ffdas->open("B:/FFDAS/urban_boundary/grump.tif");
	//OGREnvelope bound;
	//bound.MinX = -124.8;
	//bound.MaxX = -66.9;
	//bound.MinY = 24.5;
	//bound.MaxY = 49.5;
	//ffdas->crop(bound, "E:/grump.tif");
	//delete ffdas;
	//Preprocessor::gridFolderByRaster("B:/LA_Version2/Merged/", "B:/LA_Version2/Merged/gridded/","B:/Hestia_FFDAS_ODIAC/Comparison/LA_Hestia.tif");
	//makeTIFF("B:/LA_Version2/Merged/", "B:/LA_Version2/Merged/gridded/", "B:/Hestia_FFDAS_ODIAC/Comparison/LA_Hestia.tif");
	//checkTotalsInDir("B:/Hestia_FFDAS_ODIAC/Comparison/", ".tif");
	//checkTotal("B:/Hestia_FFDAS_ODIAC/Comparison/LA_Hestia.tif");
	//double total = checkTotal("B:/LA_Version2/Merged/gridded/total.tif");
	//checkTotal("B:/Hestia_FFDAS_ODIAC/Comparison/LA_FFDAS.tif");
	//checkTotal("B:/Hestia_FFDAS_ODIAC/Comparison/LA_ODIAC.tif");
	//GDAL_DS<int>* ffdas = new GDAL_DS<int>();
	//ffdas->open("B:/Hestia_FFDAS_ODIAC/Comparison/LA_FFDAS.tif",GDALAccess::GA_Update);
	//ffdas->scaleTotal(1, total);
	//delete ffdas;

	//GDAL_DS<int>* odiac = new GDAL_DS<int>();
	//odiac->open("B:/Hestia_FFDAS_ODIAC/Comparison/LA_ODIAC.tif", GDALAccess::GA_Update);
	//odiac->scaleTotal(1, total);
	//delete odiac;

	//scaleRaster("C:/VulcanGridding/urban_outlines/Yuyu/ODIAC2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/Yuyu/FFDAS2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/ODIAC2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/FFDAS2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1/ODIAC2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1/FFDAS2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5/ODIAC2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5/FFDAS2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/GRUMP/ODIAC2011.tif", 0.001);
	//scaleRaster("C:/VulcanGridding/urban_outlines/GRUMP/FFDAS2011.tif", 0.001);
	//
	//checkTotal("C:/VulcanGridding/urban_outlines/Yuyu/FFDAS2011.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/Yuyu/Total.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/FFDAS2011.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/Total.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1/FFDAS2011.tif");

	//scaleRaster("E:/Vulcan/gridPrep_SHP_master/spatial/ODIAC2011.tif", 0.001);
	//scaleRaster("E:/Vulcan/gridPrep_SHP_master/spatial/FFDAS2011.tif", 0.001);
	//scaleRaster("E:/Vulcan/gridPrep_SHP_master/emissions/ODIAC2011.tif", 0.001);
	//scaleRaster("E:/Vulcan/gridPrep_SHP_master/emissions/FFDAS2011.tif", 0.001);

	//checkTotal("E:/Vulcan/gridPrep_SHP_master/emissions/FFDAS2011.tif");
	//checkTotal("E:/Vulcan/gridPrep_SHP_master/emissions/ODIAC2011.tif");
	//GDAL_DS<double>* odiac = new GDAL_DS<double>();
	//odiac->open("C:/VulcanGridding/urban_outlines/Yuyu/ODIAC2011.tif",GDALAccess::GA_Update);
	//odiac->multiply(0.001, 1);
	//delete odiac;
	//odiac = new GDAL_DS<double>();
	//odiac->open("C:/VulcanGridding/urban_outlines/Yuyu/FFDAS2011.tif", GDALAccess::GA_Update);
	//odiac->multiply(0.001, 1);
	//delete odiac;
	//checkTotal("E:/Vulcan/gridPrep_SHP_master/emissions/FFDAS2011.tif");

	

	//
	//checkTotal("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/population.nc");
	//checkTotal("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/onroad.tif");
	//checkTotalsInDir("C:/VulcanGridding/VulcanGridAlaska/GriddedFFCO2/", ".tif");
	//checkTotalsInDir("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/", ".tif");
	//checkTotalsInDir("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/", ".tif");
	//checkTotalsInDir("C:/VulcanGridding/SphericalGridContiguousUS/GriddedFFCO2/", ".tif");
    //std::vector<std::string> strs;
    //ShapeFile::copyDir("E:/Vulcan/output_data/onroad/Shapes/", "E:/Vulcan/output_data/onroad/", strs);
	return 0;
	
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5",".tif");

	
	//reallocateEmissionsByPop("GRUMP_CCA_d1500_l5");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/GRUMP_CCA_d1500_l5", ".tif");
	

	//reallocateEmissionsByPop("GRUMP");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/GRUMP",".tif");
	//
	//std::string indir = "B:/Hestia_FFDAS_ODIAC/Comparison/";
	//std::string outdir = "B:/Hestia_FFDAS_ODIAC/Comparison/before_scaling/";
	//std::string names[4] = {"LA","SaltLake" ,"Baltimore" ,"Indianapolis"};
	//std::string srcname = "ODIAC";
	//for (size_t i = 0; i < 4; i++)
	//{
	//	std::string srcDT = indir + srcname + ".tif";
	//	std::string maskDT = indir + names[i] + "_" + srcname + ".tif";
	//	std::string outDT = outdir + names[i] + "_" + srcname + ".tif";
	//	GDAL_DS<int>* ds = new GDAL_DS<int>();
	//	ds->open(srcDT);
	//	ds->extractByMask(maskDT, outDT);
	//	delete ds;
	//}
	//checkTotal("B:/Hestia_FFDAS_ODIAC/Comparison/before_scaling/LA_ODIAC.tif");
	//checkTotal("B:/Hestia_FFDAS_ODIAC/Comparison/after_scaling/LA_ODIAC.tif");
	//checkTotal("B:/Hestia_FFDAS_ODIAC/Comparison/raw/LA_ODIAC.tif");
//checkTotalsInDir("B:/Hestia_FFDAS_ODIAC/Comparison/odiac/", ".tif");
//checkTotalsInDir("B:/Hestia_FFDAS_ODIAC/Comparison/odiac_scaled/", ".tif");
//checkTotalsInDir("B:/Hestia_FFDAS_ODIAC/Comparison/hestia/", ".tif");
	return 0;


	/*GDAL_DS<float>* ffdas = new GDAL_DS<float>();
	ffdas->open("E:/Vulcan/gridPrep_SHP_master/emissions/FFDAS2011.tif");
	ffdas->extractByMask("C:/VulcanGridding/SphericalGridContiguousUS/ContiguousUS.tif",
	"E:/Vulcan/gridPrep_SHP_master/spatial/FFDAS2011.tif");
	delete ffdas;

	GDAL_DS<float>* odiac = new GDAL_DS<float>();
	odiac->open("E:/Vulcan/gridPrep_SHP_master/emissions/ODIAC2011.tif");

	odiac->extractByMask("C:/VulcanGridding/SphericalGridContiguousUS/ContiguousUS.tif",
	"E:/Vulcan/gridPrep_SHP_master/spatial/ODIAC2011.tif");
	delete odiac;*/
	//double odiac = checkTotal("E:/Vulcan/gridPrep_SHP_master/spatial/ODIAC2011.tif");
	//odiac = checkTotal("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/ODIAC2011.tif");
	//double ffdas = checkTotal("E:/Vulcan/gridPrep_SHP_master/spatial/FFDAS2011.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/FFDAS2011.tif");


	//std::vector<std::string> rasterFiles;
	//rasterFiles.push_back("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif");
	//rasterFiles.push_back("C:/VulcanGridding/urban_outlines/Yuyu.tif");
	//rasterFiles.push_back("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif");
	//rasterFiles.push_back("C:/VulcanGridding/urban_outlines/Grump.tif");
	//makeCommonGrid(rasterFiles);

	/*	reallocateEmissionsByPop("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif",
	"C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA_FULL/fe_2007_us_zcta500.dbf.bin",
	"C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA_FULL/us_census_bk.dbf.bin",
	"E:/Vulcan/gridPrep_SHP_master/emissions/fe_2007_us_zcta500.dbf",
	"E:/Vulcan/gridPrep_SHP_master/emissions/us_census_bk.dbf",
	"C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA_FULL/Jones_Kammen_Total.tif");*/
	//checkTotal("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA_FULL/Jones_Kammen_Total.tif");
	//checkTotal("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/Jones_Kammen_Total.tif");
	//linkUrbanOutline("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif", 
	//	"C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA_FULL/");



	//
	//
	//return 0;
	//linkYuyu2ACS("C:/VulcanGridding/urban_outlines/Yuyu.tif", 
	//	"C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif",
	//	"C:/VulcanGridding/urban_outlines/Yuyu.dbf",
	//	"C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.dbf");
	//remapYuyu
	//	("C:/VulcanGridding/urban_outlines/Yuyu.shp",
	//	"C:/VulcanGridding/urban_outlines/linkedYuyu/Yuyu.shp",
	//	"C:/VulcanGridding/urban_outlines/Yuyu.tif",
	//	"C:/VulcanGridding/urban_outlines/linkedYuyu/Yuyu.tif");
	//return 0;


	//calUrbanArea("C:/VulcanGridding/urban_outlines/Yuyu.tif");
	//calUrbanArea("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif");
	//calUrbanArea("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif");
	//calUrbanArea("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif");

	//std::vector<std::string> dropfields;
	//raster2csv("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/");
	/*cumsum("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/population.tif",
	"C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/total.tif",
	"C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/pop_total.csv");*/
	//jointBlockGroup("E:/LASpatialAnalysis/FetchCensusData/census_data.csv", "E:/LASpatialAnalysis/FetchCensusData/national_bk_shapes.Shp");
	//ShapeFile::copyDirDropGeometry("C:/VulcanGridding/urban_outlines/Urban results/", "C:/VulcanGridding/urban_outlines/Urban results/points/", dropfields);
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA/", ".tif");

	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1/", ".tif");
	//linkUrbanOutline("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif", "C:/VulcanGridding/urban_outlines/Place_2010Census_DP1/");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/Grump/", ".tif");
	//linkUrbanOutline("C:/VulcanGridding/urban_outlines/Grump.tif", "C:/VulcanGridding/urban_outlines/Grump/");
	//Preprocessor::reprojectDir("E:/Vulcan/gridPrep_SHP_master/spatial/reprojected/", "E:/Vulcan/gridPrep_SHP_master/spatial_wgs84/", "E:/Vulcan/gridPrep_SHP_master/VulcanToWGS84.py");
	//linkUrbanOutline("C:/VulcanGridding/urban_outlines/Yuyu.tif", "C:/VulcanGridding/urban_outlines/Yuyu/");
	//checkTotalsInDir("C:/VulcanGridding/urban_outlines/Yuyu/", ".tif");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/ElecProd.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/OnRoad.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Commercial.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Residential.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Industrial.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Airport.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Cement.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Railroad.nc");
	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Nonroad.nc");

	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Total.nc");

	//checkVulcanFiles();
	//checkTotalsInDir("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/", ".tif");
	//checkTotalsInDir("C:/VulcanGridding/SphericalGridContiguousUS/GriddedFFCO2/", ".tif");
	//checkTotalsInDir("C:/VulcanGridding/VulcanGridAlaska/GriddedFFCO2/", ".tif");
	//checkTotalsInDir("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/", ".tif");

	//checkTotal("C:/VulcanGridding/SphericalGridAlaska/GriddedFFCO2/Total.tif");
	//checkTotal("C:/VulcanGridding/VulcanGridAlaska/GriddedFFCO2/Total.tif");
	//checkTotal("C:/VulcanGridding/SphericalGridContiguousUS/GriddedFFCO2/Total.tif");
	//checkTotal("C:/VulcanGridding/VulcanGridContiguousUS/GriddedFFCO2/Total.tif");
	//
	return 0;
	//ShapeFile Shp("E:/Cement.dbf");
	//OGRFeature* poFeature;
	//int fid = -1;
	//while ((poFeature = Shp.poLayer->GetNextFeature()) != NULL)
	//{
	//	OGRGeometry* poGeometry = poFeature->GetGeometryRef();
	//	if (!poGeometry)
	//	{
	//		double value = poFeature->GetFieldAsDouble("ca11");
	//		printf("%f\n", value);
	//	}
	//}
	/*TwinCitiesComparison();
	return 0;
	std::string srcdir = "E:/Vulcan/gridPrep_SHP_master - Copy/";
	QDir input_dir(srcdir.outlineDS());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	srcdir = (input_dir.absolutePath() + "/").toLocal8Bit().outlineDS();

	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
	QFileInfo fileInfo = list.at(i);
	std::string name = fileInfo.fileName().toLocal8Bit().outlineDS();
	std::string name2 = name;
	bool found = false;
	if (fileInfo.fileName().endsWith(".xml"))
	{
	continue;
	}
	if (!fileInfo.fileName().startsWith("Nonroad."))
	{
	continue;
	}
	for (size_t n = 0; n < name2.size() - 4; n++)
	{
	if (name2[n] == '.')
	{
	name2[n] = '_';
	found = true;
	}
	}
	printf("%s,%s\n", name.outlineDS(), name2.outlineDS());
	QFile::rename((srcdir + name).outlineDS(), (srcdir + name2).outlineDS());

	}*/

	//std::vector<std::string> fields;
	//fields.push_back("ca11");
	//fields.push_back("area");	
	//fields.push_back("AREA");
	//fields.push_back("length");
	//fields.push_back("LENGTH");
	//fields.push_back("ca11");
	//Utils::renameFieldsForDir("C:/VulcanGridding/gridPrep_SHP_master/", srcfields, destfields);
	//Utils::dropFieldsForDir("C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/", fields);
	//Utils::dropFieldsForDir("C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/", fields);
	//Utils::dropFieldsForDir("C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/", fields);
	//Utils::copyFolder("E:/Vulcan/gridPrep_SHP_master/", "E:/Vulcan/gridPrep_SHP_master/clean/", fields);
	//dropFields("C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/", fields);
	//dropFields("C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/", fields);
	//dropFields("C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/", fields);
	//Preprocessor::gridFolderByShape("C:/VulcanGridding/gridPrep_SHP_master/", "C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/", "C:/VulcanGridding/urban_outlines/cb_2012_us_uac10_500k.shp");
	//Preprocessor::gridFolderByShape("C:/VulcanGridding/gridPrep_SHP_master/", "C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/", "C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.shp");
	//Preprocessor::gridFolderByShape("C:/VulcanGridding/gridPrep_SHP_master/", "C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/", "C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.shp");

	//copyFields("C:/VulcanGridding/gridPrep_SHP_master/","C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/", fields);
	//copyFields("C:/VulcanGridding/gridPrep_SHP_master/", "C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/", fields);
	//copyFields("C:/VulcanGridding/gridPrep_SHP_master/", "C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/", fields);
	//Utils::dropFieldsForDir("C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/", fields);
	//Utils::dropFieldsForDir("C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/", fields);
	//Utils::dropFieldsForDir("C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/", fields);
	//Preprocessor::updateFieldAfterIntersectionForDir("C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/");
	//Preprocessor::updateFieldAfterIntersectionForDir("C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/");
	//Preprocessor::updateFieldAfterIntersectionForDir("C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/");
	Vulcan2014 VULCAN;
	//VULCAN.createAirport();
	//VULCAN.createElecProd();
	//VULCAN.createComPoint();
	//VULCAN.createIndPoint();
	//VULCAN.createCementPoint();
	//VULCAN.createNonroadPoint();
	//VULCAN.createRailroadPoint();
	//vulcan.sumBySectorRaster("C:/VulcanGridding/Gridded/", "C:/VulcanGridding/Gridded/totals.csv");
	//system("RScript B:/Shapefiles2Grid/cal_dbf_total.R C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/ ca11");
	//system("RScript B:/Shapefiles2Grid/cal_dbf_total.R C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/ ca11");
	//system("RScript B:/Shapefiles2Grid/cal_dbf_total.R C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/ ca11");
	//system("RScript B:/Shapefiles2Grid/cal_dbf_total.R E:/Vulcan/gridPrep_SHP_master/ ca11");
	//VULCAN.sumBySectorShape("C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/", "C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA_totals.csv");
	//VULCAN.sumBySectorShape("C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/", "C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1_totals.csv");
	//VULCAN.sumBySectorShape("C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/", "C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k_totals.csv");
	VULCAN.sumBySectorShape("E:/Vulcan/gridPrep_SHP_master/emissions/", "E:/Vulcan/gridPrep_SHP_master/totals.csv");
	return 0;
	//GDAL_DS<unsigned int>* ds = new GDAL_DS<unsigned int>();
	//ds->open("C:/VulcanGridding/urban_outlines/cb_2012_us_uac10_500k.tif");
	//GDALRasterAttributeTable* tb = ds->dataset->GetRasterBand(1)->GetDefaultRAT();
	//for (int i = 0; i < 100; i++)
	//{
	//	std::string val = tb->GetValueAsString(i, 0);
	//	printf("%s\n", val.outlineDS());
	//}

	std::vector<std::string> masks;
	masks.push_back("C:/VulcanGridding/urban_outlines/cb_2012_us_uac10_500k.tif");
	masks.push_back("C:/VulcanGridding/urban_outlines/Place_2010Census_DP1.tif");
	masks.push_back("C:/VulcanGridding/urban_outlines/ACS_2015_5YR_UA.tif");

	std::vector<std::string> outdirs;
	outdirs.push_back("C:/VulcanGridding/Intersected_urban/cb_2012_us_uac10_500k/");
	outdirs.push_back("C:/VulcanGridding/Intersected_urban/Place_2010Census_DP1/");
	outdirs.push_back("C:/VulcanGridding/Intersected_urban/ACS_2015_5YR_UA/");

	std::vector<std::string> onroadfiles;
	onroadfiles.push_back("C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadRuralLocal.shp");
	onroadfiles.push_back("C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadRuralNonLocal.shp");
	onroadfiles.push_back("C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadUrbanLocal.shp");
	onroadfiles.push_back("C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadUrbanNonLocal.shp");
	//for (size_t i = 0; i < onroadfiles.size(); i++)
	//{
	//	QFileInfo info(onroadfiles[i].outlineDS());
	//	std::vector<std::string> outputfiles;
	//	for (size_t j = 0; j < masks.size(); j++)
	//	{
	//		outputfiles.push_back(outdirs[j] + info.fileName().toLocal8Bit().outlineDS());
	//	}
	//	polyline2points(328.084, onroadfiles[i], masks, outputfiles);
	//}
	//polyline2points(328.084, "C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadRuralLocal.Shp", masks);
	//polyline2points(328.084,"C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadRuralLocal.Shp", "C:/VulcanGridding/gridPrep_SHP_master/OnroadRuralLocal.Shp");
	//polyline2points(328.084, "C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadRuralNonLocal.Shp", "C:/VulcanGridding/gridPrep_SHP_master/OnroadRuralNonLocal.Shp");

	//polyline2points(328.084, "C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadUrbanLocal.Shp", "C:/VulcanGridding/gridPrep_SHP_master/OnroadUrbanLocal.Shp");
	//polyline2points(328.084, "C:/VulcanGridding/gridPrep_SHP_master/OnRoad/OnroadUrbanNonLocal.Shp", "C:/VulcanGridding/gridPrep_SHP_master/OnroadUrbanNonLocal.Shp");

	//MODIS_LST_2_NEIGHBORHOODS("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint100ft.tif",
	//	"E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/Neighborhood100ft.tif",
	//	"E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/MODIS_LST_2010_12-2.tif");
	//Utils::updateFootprintForDir("C:/VulcanGridding/urban_outlines/");

	//TwinCitiesComparison();
	return 0;
	//updateID();
	//return 0;
	//std::vector<std::string> fields;
	//ExtractByPolygon<char>::link("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint.Shp", "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/bg06037.Shp", "bk_id", fields);
	//ExtractByPolygon<char>::link("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint.Shp", "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNeighborhoods.Shp", "nei_id", fields);

	//ExtractByPolygon<char>::link2("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/bg06037.Shp", "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint.Shp", "bk_id", fields);
	//ExtractByPolygon<char>::link2("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/CAResNeighborhoods.Shp", "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint.Shp", "nei_id", fields);
	//return 0;
	//parcels2neighborhood();
	//linkBlockgroups2Neighborhoods("E:/LASpatialAnalysis/FetchCensusData/bg06037.Shp", "E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/LA_Neighborhoods/Neighborhoods.tif");

	//linkParcels2UtilityArea(dir + "ComNonPoint.Shp", dir + "LA_Neighborhoods/LA_Utility_Areas.tif");
	//linkParcels2UtilityArea(dir + "ResNonPoint.Shp", dir + "LA_Neighborhoods/LA_Utility_Areas.tif");
	//linkParcels2UtilityArea(dir + "IndNonPoint.Shp", dir + "LA_Neighborhoods/LA_Utility_Areas.tif");

	//return 0;
	//std::vector<std::string> CBECS12NEEUI_files;
	//CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI (Meged with California Survey)/CBECS12NEEUI.csv");
	////CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Pre-1980.csv");
	////CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Post-1979.csv");
	/*std::vector<std::string> fieldnames;
	fieldnames.push_back("ca10_ng");
	fieldnames.push_back("ca10");
	fieldnames.push_back("kwh");
	fieldnames.push_back("tfs");
	Preprocessor::updateFieldAfterIntersection("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint_Neiborhoods.Shp", fieldnames);*/
	//std::vector<std::string> RECS09NEEUI_files;
	//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Pre-1980.csv");
	//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Post-1979.csv");
	//std::vector<int> breaks2;
	//breaks2.push_back(1980);	
	//int censusDivision = 9;
	//NeiborOptimMatrixPrep neiborprep;
	//neiborprep.updateNEEUIClassCode("E:/LASpatialAnalysis/LA_Comparison/Los_Angeles/ResNonPoint.Shp", RECS09NEEUI_files, breaks2, censusDivision);
	return 0;
	//std::vector<std::string> MECSNEEUI_files;
	//MECSNEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/Mecs10NEEUI.csv");
	//std::vector<std::string> CBECS12NEEUI_files;
	//CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI.csv");

	//std::vector<std::string> RECS09NEEUI_files;
	//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI.csv");

	/*std::vector<int> breaks;
	breaks.push_back(10000);
	int censusDivision = 9;
	std::vector<int> breaks1;
	breaks1.push_back(10000);

	std::vector<int> breaks2;
	breaks2.push_back(1980);

	NonPointProcessor::updateNEEUI(dir + "ComNonPoint.Shp", CBECS12NEEUI_files, breaks1, censusDivision);
	NonPointProcessor::updateNEEUI(dir + "ResNonPoint.Shp", RECS09NEEUI_files, breaks2, censusDivision);
	NonPointProcessor::updateNEEUI(dir + "IndNonPoint.Shp", MECSNEEUI_files, breaks2, censusDivision);
	*/

	//FID_LA_Nei ComElec.Shp geographies.setup("E:/LASpatialAnalysis/LA_Comparison/ComElecNeighborhoods2.Shp", );
	//
	//return 0;
	//gridMarion();
	//double total = 0;

	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.airport.annual.2011.v3.0.nc");
	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.com.annual.2011.v3.0.nc");
	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.res.annual.2011.v3.0.nc");
	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.ind.annual.2011.v3.0.nc");
	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.nonroad.annual.2011.v3.0.nc");
	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.onroad.annual.2011.v3.0.nc");
	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.elecprod.annual.2011.v3.0.nc");
	//total += checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.railroad.annual.2011.v3.0.nc");


	//printf("%f\n", total);
	//total = checkTotal("Z:/Hestia/Indy/INDIANAPOLIS_NEI2011/gridded/Marion.total.annual.2011.v3.0.nc");
	//printf("%f\n", total);

	//jointBlockGroup("E:/LASpatialAnalysis/CensusShapes/pop_income.csv","E:/LASpatialAnalysis/CensusShapes/gz_2010_06_150_00_500k.Shp");
	//jointCensusTract("E:/LASpatialAnalysis/CensusShapes/pop_income.csv", "E:/LASpatialAnalysis/ca11_CensusTract/res.Shp");
	//jointCensusTract("E:/LASpatialAnalysis/CensusShapes/pop_income.csv", "E:/LASpatialAnalysis/ca11_CensusTract/com.Shp");
	//jointCensusTract("E:/LASpatialAnalysis/CensusShapes/pop_income.csv", "E:/LASpatialAnalysis/ca11_CensusTract/hestia.Shp");
	//jointCensusTract("E:/LASpatialAnalysis/CensusShapes/pop_income.csv", "E:/LASpatialAnalysis/ca11_CensusTract/total.Shp");
	//jointCensusTract("E:/LASpatialAnalysis/CensusShapes/pop_income.csv", "E:/LASpatialAnalysis/ca11_CensusTract/ffdas.Shp");
	//jointCensusTract("E:/LASpatialAnalysis/CensusShapes/pop_income.csv", "E:/LASpatialAnalysis/ca11_CensusTract/odiac.Shp");
	//aggregrateTIFF("E:/LASpatialAnalysis/1km/");
	//aggregrateTIFF("E:/LASpatialAnalysis/10km/");
	//mergeFields("E:/LASpatialAnalysis/10km/", "E:/LASpatialAnalysis/LA_10km_grid.Shp");
	//aggregrateTIFF(10,"E:/LASpatialAnalysis/1km/", "E:/LASpatialAnalysis/10km/");
	//gridLA_TIFF();

	//-119.68896007 -114.058285362
	//35.832632771 32.689068972
	//char wkt[512];
	//char* pwkt = wkt;
	//if (input.poLayer->GetSpatialRef())
	//	input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	/*Grid grid;
	grid.fromFishnetRaster("E:/LASpatialAnalysis/lspop2011.tif", true);
	OGREnvelope bound;
	bound.MinX = -119.68896007; bound.MaxX = -114.058285362; bound.MinY = 32.689068972; bound.MaxY = 35.832632771;*/
	//grid.toShape(grid.proj, "E:/LASpatialAnalysis/lspop2011_LA_clip.Shp", bound, true);
	//Utils::updateFootprint("E:/LASpatialAnalysis/lspop2011_LA_clip.Shp", 10000);
	//Preprocessor::intersectWithArcGIS("E:/LASpatialAnalysis/lspop2011_LA_clip.Shp", "B:/LA_Version2/gridoutput/fishnet1000m.Shp", "E:/LASpatialAnalysis/lspop2011_LA.Shp");
	//Preprocessor::updateFieldAfterIntersection("E:/LASpatialAnalysis/lspop2011_LA.Shp",10000);
	//Grid grid;
	//grid.fromFishnetShape("B:/LA_Version2/gridoutput/fishnet1000m.Shp");
	//grid.reset();
	//ShapeFile fishnet("E:/LASpatialAnalysis/lspop2011_LA.Shp");
	//grid.gatherCells(&fishnet,"ca11", 10000);
	//grid.toRaster("E:/LASpatialAnalysis/lspopLA.tif");
	//grid.toShape(grid.proj,"E:/LASpatialAnalysis/lspopLA.Shp",true);
	//calIntensity("E:/LASpatialAnalysis/ca11_CensusTract/");
	//gridLA_With_FFDAS();
	//Utils::updateFootprintForDir("E:/LASpatialAnalysis/ca11_CensusTract/", true,9.2903e-8);

	//replaceNoData("E:/LASpatialAnalysis/FFDAS_Hestia/FFDAS2011.tif", "E:/LASpatialAnalysis/FFDAS_Hestia/FFDAS.tif");
	//replaceNoData("E:/LASpatialAnalysis/FFDAS_Hestia/ODIAC2011.tif", "E:/LASpatialAnalysis/FFDAS_Hestia/ODIAC.tif");
	//replaceNoData("E:/LASpatialAnalysis/FFDAS_Hestia/Hestia2011.tif", "E:/LASpatialAnalysis/FFDAS_Hestia/Hestia.tif");

	//tif2csv("E:/LASpatialAnalysis/FFDAS_Hestia/FFDAS.tif", "E:/LASpatialAnalysis/FFDAS_Hestia/FFDAS.csv");
	//tif2csv("E:/LASpatialAnalysis/FFDAS_Hestia/ODIAC.tif", "E:/LASpatialAnalysis/FFDAS_Hestia/ODIAC.csv");
	//tif2csv("E:/LASpatialAnalysis/FFDAS_Hestia/Hestia.tif", "E:/LASpatialAnalysis/FFDAS_Hestia/Hestia.csv");
	//reallocateLA();
	//reallocateHestiaFFDAS();
	//gridLA_NoShipping();
	//regrid_FFDASODIAC_LA();
	//Utils::updateFootprintForDir("E:/LASpatialAnalysis/FFDAS_Hestia/");
	//

	//return 0;
	//intersect();
	//OGREnvelope bb = BoundManager::readBoundFromShape("E:/Rasterization_gridding/ComNonPoint.Shp");
	//Grid grid(bb, 656.168, 0);
	//ShapeFile Shp("E:/Rasterization_gridding/ComNonPoint.Shp");
	//char wkt[512];
	//char* pwkt = wkt;
	//if (Shp.poLayer->GetSpatialRef())
	//	Shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	//Shp.close();
	//grid.reset();
	//std::string fishnet = "E:/Rasterization_gridding/fishnet";
	//grid.toShape(wkt, fishnet + ".Shp");
	//grid.toRaster(fishnet + ".tif", wkt);
	//Preprocessor::updateFieldAfterIntersection("E:/Rasterization_gridding/Intersection/ComNonPoint_Intersect.Shp");

	//ShapeFile Shp("B:/Baltimore/gridPrep_SHP_master/ca_1.2/200m/OnRoad.Shp");
	//char wkt[512];
	//char* pwkt = wkt;
	//if (Shp.poLayer->GetSpatialRef())
	//	Shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
	//Grid grid;
	//grid.fromFishnetRaster("B:/Baltimore/gridPrep_SHP_master/ca_1.2/200m/fishnet.tif", true);
	//grid.reset();
	//grid.gatherCells(&Shp);
	//grid.toShape(pwkt, "E:/Rasterization_gridding/Rasterization/OnRoad2.Shp",true);

	//Preprocessor::intersectWithArcGIS("E:/Rasterization_gridding/ComNonPoint.Shp", "E:/Rasterization_gridding/fishnet.Shp", "E:/Rasterization_gridding/Intersection/ComNonPoint.Shp");

	//Preprocessor::gridShapeFile("E:/Vulcan/fishnet.Shp", "E:/Vulcan/gridPrep_SHP_master/OnroadUrbanNonLocal.Shp", "E:/Vulcan/gridPrep_SHP_master/OnroadUrbanNonLocal.Shp");

	//Preprocessor::gridFolderByRaster("C:/VulcanGridding/gridPrep_SHP_master/", "C:/VulcanGridding/Intersected/", "C:/VulcanGridding/fishnet.tif");

	//reprojectVulcan();
	//calVulcanExtent();
	//mergeRoads("E:/Vulcan/onroad/LocalRoads/", "E:/Vulcan/onroad/OnroadUrbanLocal.Shp", "Urban");
	//mergeRoads("E:/Vulcan/onroad/LocalRoads/", "E:/Vulcan/onroad/OnroadRuralLocal.Shp", "Rural");

	//mergeRoads("E:/Vulcan/onroad/GapfilledRoads_Raw/", "E:/Vulcan/onroad/OnroadUrbanNonLocal.Shp", "u");
	//mergeRoads("E:/Vulcan/onroad/GapfilledRoads_Raw/", "E:/Vulcan/onroad/OnroadRuralNonLocal.Shp", "r");

	//E:\Vulcan\onroad\GapfilledRoads_Raw

	//return 0;
	//calSVF_Kechunag();
	//exportParcelURL();
	//reallocateLA();
	//	aggregrate("B:/SpatialGranuality/Baltimore/blockgroup/intersected/", "B:/SpatialGranuality/Baltimore/blockgroup/blockgroup.Shp");
	//gridFolderByRaster("B:/SpatialGranuality/Baltimore/", "B:/SpatialGranuality/Baltimore/blockgroup/intersected/", "B:/SpatialGranuality/Baltimore/blockgroup/blockgroup.Shp");
	//gridFolderByRaster("B:/SpatialGranuality/Baltimore/", "B:/SpatialGranuality/Baltimore/100m/intersected/", 328.084, false);
	//aggregrateSingleFile("B:/SpatialGranuality/Baltimore/100m/points_blockgroup.Shp", "B:/SpatialGranuality/Baltimore/100m/blockgroup.Shp");
	//gridFolderByRaster("B:/SpatialGranuality/Baltimore/", "B:/SpatialGranuality/Baltimore/200m/intersected/", 328.084*2, false);
	//gridFolderByRaster("B:/SpatialGranuality/Baltimore/", "B:/SpatialGranuality/Baltimore/500m/intersected/", 328.084*5, false);

	//gridFolderByRaster("B:/SpatialGranuality/Baltimore/", "B:/SpatialGranuality/Baltimore/300m/intersected/", 328.084 * 3, false);
	//gridFolderByRaster("B:/SpatialGranuality/Baltimore/", "B:/SpatialGranuality/Baltimore/400m/intersected/", 328.084 * 4, false);

	//aggregrateSingleFile("B:/SpatialGranuality/Baltimore/200m/points_blockgroup.Shp", "B:/SpatialGranuality/Baltimore/200m/blockgroup.Shp");
	//aggregrateSingleFile("B:/SpatialGranuality/Baltimore/500m/points_blockgroup.Shp", "B:/SpatialGranuality/Baltimore/500m/blockgroup.Shp");
	//gridBaltimoreAggregrateByBlockGroup("B:/SpatialGranuality/Baltimore/blockgroups/blockgroup.Shp", "B:/SpatialGranuality/Baltimore/blockgroups/fishnet10m.tif", "B:/SpatialGranuality/Baltimore/blockgroups/");
	//return 0;
	//std::vector<std::string> fields;
	//fields.push_back("Length");

	//ShapeFile::copyDropGeometry("E:/GoogleStreetview/Manhattan/StreetSegment_Clip.Shp", "E:/GoogleStreetview/Manhattan/StreetSegmentPoints.Shp", fields);
	//ShapeFile::copyDropGeometry("E:/TwinCitiesComparison/streets_downtown_(no_highway).Shp", "E:/TwinCitiesComparison/streets_points.Shp", fields);

	//segments2points("E:/TwinCitiesComparison/streets_downtown_(no_highway).Shp", "E:/TwinCitiesComparison/streets_points.Shp", "E:/TwinCitiesComparison/streets_points.csv");
	//linkIntersected2Original("E:/HestiaUncertainty/LACounty100m/Intersected/","B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/");
	//IntersectLA1000to100();
	//combineGrids("B:/FFDAS/withroad_vesus_withoutroad/withroad_other_2010.nc", "B:/FFDAS/withroad_vesus_withoutroad/ffdas_flux_road_2010.nc","B:/FFDAS/withroad_vesus_withoutroad/withroad_other_2010.tif");
	//cleanGrid("B:/FFDAS/withroad_vesus_withoutroad/withoutroad_other_2010.nc", "B:/FFDAS/withroad_vesus_withoutroad/withoutroad_other_2010.tif");
	//toCSV("B:/FFDAS/withroad_vesus_withoutroad/non_zero/ffdas_flux_road_2010.tif", "B:/FFDAS/withroad_vesus_withoutroad/non_zero/ffdas_flux_road_2010.csv");
	//toCSV("B:/FFDAS/withroad_vesus_withoutroad/non_zero/withoutroad_other_2010.tif", "B:/FFDAS/withroad_vesus_withoutroad/non_zero/withoutroad_other_2010.csv");
	//toCSV("B:/FFDAS/withroad_vesus_withoutroad/non_zero/withroad_other_2010.tif", "B:/FFDAS/withroad_vesus_withoutroad/non_zero/withroad_other_2010.csv");
	//count("C:/FFDAS/1km/lspop2011.tif");
	//toCSV("C:/FFDAS/1km/lspop2011.tif", "C:/FFDAS/1km/lspop2011.csv");
	//grid1km();
	//grid10km();
	//subtractGrids("B:/FFDAS/withroad_vesus_withoutroad/non_zero/withoutroad_other_2010.tif", "B:/FFDAS/withroad_vesus_withoutroad/non_zero/withroad_other_2010.tif", "B:/FFDAS/withroad_vesus_withoutroad/non_zero/percent_diff.tif");

	//gridLA();
	//updateRail("E:/rail/rail2011.csv","shape_id","final.CO2.tC","E:/rail/rail_lines.Shp","FRAARCID");
	///updateRail("E:/rail/S260FRAC.csv", "US_RAIL_", "FIPS", "E:/rail/us_rail2k.Shp","US_RAIL_");
	//Utils::dbf2csv("E:/rail/rail_lines.Shp", "E:/rail/railway.csv");
	//Utils::dbf2csv("E:/rail/us_rail2k.Shp", "E:/rail/us_rail2k.csv");
	//Utils::dbf2csv("E:/CMV_shapes/Ports_030216.Shp", "E:/CMV_shapes/Ports_030216.csv");
	//Utils::dbf2csv("E:/CMV_shapes/Ports_091913.Shp", "E:/CMV_shapes/Ports_091913.csv");
	//Utils::dbf2csv("E:/CMV_shapes/ShippingLanes_080515.Shp", "E:/CMV_shapes/ShippingLanes_080515.csv");
	//Utils::dbf2csv("E:/CMV_shapes/ShippingLanes_112812_FINAL.Shp", "E:/CMV_shapes/ShippingLanes_112812_FINAL.csv");
	//updateRail("E:/CMV_shapes/Nonpoint.rail_CMV.em.CO2.final.2011 (marine only).csv", "shape_id", "final.CO2.tC", "E:/CMV_shapes/Ports_030216.Shp", "ShapeID");
	//updateRail("E:/CMV_shapes/Nonpoint.rail_CMV.em.CO2.final.2011 (marine only).csv", "shape_id", "final.CO2.tC", "E:/CMV_shapes/Ports_091913.Shp", "ShapeID");
	//updateRail("E:/CMV_shapes/Nonpoint.rail_CMV.em.CO2.final.2011 (marine only).csv", "shape_id", "final.CO2.tC", "E:/CMV_shapes/ShippingLanes_080515.Shp", "ShapeID");
	//updateRail("E:/CMV_shapes/Nonpoint.rail_CMV.em.CO2.final.2011 (marine only).csv", "shape_id", "final.CO2.tC", "E:/CMV_shapes/ShippingLanes_112812_FINAL.Shp", "ShapeID");

	//Ports_030216.Shp
	//	Ports_091913.Shp
	//	ShippingLanes_080515.Shp
	//	ShippingLanes_112812_FINAL.Shp
	//return 0;
	//gridRail("E:/rail/rail_lines.Shp");

	//ShapeFile::copyDir("B:/LA_Version2/Public/Los_Angeles/", "Z:/Hestia/LAbasin/LAbasin_v2.0/gridPrep_SHP_master/Los_Angeles/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/Public/Los_Angeles/", "Z:/Hestia/LAbasin/LAbasin_v2.0/public_shapefiles/Los_Angeles/", fds);
	//grid100m();

	//std::vector<std::string> fds;
	//ShapeFile::copyDir("B:/LA_Version2/Public/Los_Angeles/", "B:/LA_Version2/Public2/Los_Angeles/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/Public/Orange/", "B:/LA_Version2/Public2/Orange/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/Public/San_Bernardino/", "B:/LA_Version2/Public2/San_Bernardino/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/Public/Ventura/", "B:/LA_Version2/Public2/Ventura/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/Public/Riverside/", "B:/LA_Version2/Public2/Riverside/", fds);


	//ShapeFile::copyDir("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/", "B:/LA_Version2/gridPrep_SHP_master2/Los_Angeles/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/gridPrep_SHP_master/Orange/", "B:/LA_Version2/gridPrep_SHP_master2/Orange/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", "B:/LA_Version2/gridPrep_SHP_master2/San_Bernardino/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/gridPrep_SHP_master/Ventura/", "B:/LA_Version2/gridPrep_SHP_master2/Ventura/", fds);
	//ShapeFile::copyDir("B:/LA_Version2/gridPrep_SHP_master/Riverside/", "B:/LA_Version2/gridPrep_SHP_master2/Riverside/", fds);

	//checkNCTootals("C:/HestiaGridding/Los_Angeles/result/");
	//printf("")
	//double calSum(const char* infilename)
	//gridLA();
	//std::vector<std::string> CBECS12NEEUI_files;
	//CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI (Meged with California Survey)/CBECS12NEEUI.csv");
	////CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Pre-1980.csv");
	////CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Post-1979.csv");

	//std::vector<std::string> RECS09NEEUI_files;
	//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Pre-1980.csv");
	//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Post-1979.csv");

	//std::vector<std::string> MECSNEEUI_files;
	//MECSNEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/Mecs10NEEUI.csv");
	//std::vector<std::string> CBECS12NEEUI_files;
	//CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI.csv");

	//std::vector<std::string> RECS09NEEUI_files;
	//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI.csv");

	//std::vector<int> breaks;
	//breaks.push_back(10000);
	//int censusDivision = 9;
	//std::vector<int> breaks1;
	//breaks1.push_back(10000);

	//std::vector<int> breaks2;
	//breaks2.push_back(1980);
	//std::vector<std::string> LA_DIRS;
	//LA_DIRS.push_back("Los_Angeles"); 
	//LA_DIRS.push_back("Orange");
	//LA_DIRS.push_back("San_Bernardino");
	//LA_DIRS.push_back("Ventura");
	//LA_DIRS.push_back("Riverside");
	//std::string rootdir = "B:/LA_Version2/gridPrep_SHP_master/";
	//ScaleByFuel scaleByFuel;
	//for (size_t i = 0; i < LA_DIRS.size(); i++)
	//{
	//	std::string dir = rootdir + LA_DIRS[i] + "/";
	//	//NonPointProcessor::updateNEEUI(dir + "ComNonPoint.Shp", CBECS12NEEUI_files, breaks1, censusDivision);
	//	//NonPointProcessor::updateNEEUI(dir + "ResNonPoint.Shp", RECS09NEEUI_files, breaks2, censusDivision);
	//	//createTimeStructForComNonPoint(dir + "ComNonPoint.Shp");
	//	//createTimeStructForResNonPoint(dir + "ResNonPoint.Shp");
	//	//createTimeStructForIndNonPoint(dir + "IndNonPoint.Shp");
	//	//printf("%s\n", dir);
	//	//scaleByFuel.scale(dir + "Scaling.csv", dir);
	//}

	//return 0;
	//for (size_t i = 0; i < LA_DIRS.size(); i++)
	//{
	//	std::string dir = LA_DIRS[i];
	//	NonPointProcessor::updateNEEUI(dir + "IndNonPoint.Shp", MECSNEEUI_files, breaks, -1);
	//}



	//for (size_t i = 0; i < LA_DIRS.size(); i++)
	//{
	//	std::string dir = rootdir + LA_DIRS[i];
	//	NonPointProcessor::updateNEEUI(rootdir + "IndNonPoint.Shp", MECSNEEUI_files, breaks, -1);
	//}
	/*std::vector<VintageComposition> vintages;
	for (size_t i = 0; i < LA_DIRS.size(); i++)
	{
	std::string dir = rootdir + LA_DIRS[i] + "/";
	VintageComposition comVintage;
	VintageComposition resVintage;
	VintageComposition indVintage;
	comVintage.name = LA_DIRS[i]; comVintage.sect = "COM";
	resVintage.name = LA_DIRS[i]; resVintage.sect = "RES";
	indVintage.name = LA_DIRS[i]; indVintage.sect = "IND";
	vintageStats(dir + "ComNonPoint.Shp", comVintage.zero, comVintage.pre1980, comVintage.post1979);
	vintageStats(dir + "ResNonPoint.Shp", resVintage.zero, resVintage.pre1980, resVintage.post1979);
	vintageStats(dir + "IndNonPoint.Shp", indVintage.zero, indVintage.pre1980, indVintage.post1979);
	comVintage.sum();
	resVintage.sum();
	indVintage.sum();
	vintages.push_back(comVintage);
	vintages.push_back(resVintage);
	vintages.push_back(indVintage);
	}
	std::ofstream ofs("vintage.csv");

	for (size_t i = 0; i < vintages.size(); i++)
	{
	ofs << "," << vintages[i].name << "," << vintages[i].sect << "," << vintages[i].total << "," << vintages[i].zero << "," << vintages[i].pre1980 << "," << vintages[i].post1979 << std::endl;
	}

	ofs.close();*/

	//std::vector<VintageComposition> compos;
	//
	//int zero, int& pre1980, int& post1979
	//createShpFromCSV();

	//GDAL_DS<float>* ncfds = new GDAL_DS<float>();
	////ncfds->open("E:/INDY.Total.2013.annual.nc");
	//ncfds->open("C:/Users/jliang41/Downloads/INDY.Total.2013.hourly.nc");

	//int numcells = ncfds->ncols*ncfds->nrows;
	//std::ofstream ofs;
	//ofs.open("indy2.csv");
	//for (size_t i = 1; i <= ncfds->numbands; i++)
	//{
	//	float* outlineDS = ncfds->readData(i);
	//	float minval = FLT_MAX;
	//	float maxval = -FLT_MAX;
	//	for (size_t n = 0; n < ncfds->ncols; n++)
	//	{
	//		float val = outlineDS[n];
	//		if (minval > val)
	//			minval = val;
	//		if (maxval < val)
	//			maxval = val;
	//	}
	//	ofs << i << "," << minval << "," << maxval << std::endl;
	//	printf("%d,%f,%f\n", i, minval, maxval);
	//	delete[] outlineDS;
	//}
	//ofs.close();
	//return 0;
	//computeSVF_In_Manhattan();
	//std::vector<std::string> fields;
	///ShapeFile::copyDropGeometry("B:/LA_Version2/LARIAC4Buildings/Footprints.Shp", "B:/LA_Version2/LARIAC4Buildings/Footprints_NoGeometry.Shp", fields);
	//aggregrateToParcelsUpdate("B:/LA_Version2/LARIAC4Buildings/Footprints_NoGeometry.Shp", "B:/LA_Version2/LARIAC4Buildings/Footprints.csv", "B:/LA_Version2/LARIAC4Buildings/Parcels_06037_LA.Shp");
	//aggregrateToParcels("B:/LA_Version2/LARIAC4Buildings/ResFootprints.Shp", "B:/LA_Version2/LARIAC4Buildings/ResFootprints.csv");
	//aggregrateToParcels("B:/LA_Version2/LARIAC4Buildings/ComFootprints.Shp", "B:/LA_Version2/LARIAC4Buildings/ComFootprints.csv");
	//aggregrateToParcels("B:/LA_Version2/LARIAC4Buildings/IndFootprints.Shp", "B:/LA_Version2/LARIAC4Buildings/IndFootprints.csv");

	//fields.push_back("FID_Parcel");
	//fields.push_back("FID_Export");
	//fields.push_back("HEIGHT");
	//fields.push_back("Sector");
	//fields.push_back("CRMcode");
	//fields.push_back("tfs");
	//fields.push_back("Shape_Area");


	//
	//aggregrateToFEMA("B:/LA_Version2/Fema_LA/parcels_trimmed.Shp");
	//
	//Preprocessor::intersectWithArcGIS("B:/LA_Version2/LARIAC4Buildings/Export_Output.Shp", "B:/LA_Version2/LARIAC4Buildings/Parcels_06037_LA.Shp", "B:/LA_Version2/LARIAC4Buildings/Footprints_06037_LA.Shp");


	//return 0;
	/*	gridLA_With_FFDAS();
	gridLA();
	gridLA_VY()*/;



	//Preprocessor::gridFolderByShape("B:/LA_Version2/ca/Riverside/", "B:/LA_Version2/gridoutput/Riverside/", "B:/LA_Version2/gridoutput/fishnet1000m.Shp");
	//ShapeFile::copyDirDropGeometry("B:/LA_Version2/gridoutput/Riverside/", "B:/LA_Version2/gridoutput_nogeometry/Riverside/", fields2keep);

	//Preprocessor::gridFolderByShape("B:/LA_Version2/ca/Riverside/", "C:/HestiaGridding/Los_Angeles/Vineet/intersected/Riverside/", "C:/HestiaGridding/Los_Angeles/Vineet/fishnet.Shp");
	//ShapeFile::copyDropGeometry("C:/HestiaGridding/Los_Angeles/Vineet/intersected/Riverside/IndNonPoint.Shp", "C:/HestiaGridding/Los_Angeles/Vineet/intersectedNoGeometry/Riverside/IndNonPoint.Shp", fields2keep);
	//ShapeFile::copyDropGeometry("C:/HestiaGridding/Los_Angeles/Vineet/intersected/Riverside/ComNonPoint.Shp", "C:/HestiaGridding/Los_Angeles/Vineet/intersectedNoGeometry/Riverside/ComNonPoint.Shp", fields2keep);
	//return 0;
	//Preprocessor::gridFolderByShape("B:/LA_Version2/ca/Riverside/", "B:/LA_Version2/WRF_FFDAS_GRIDDING/grid01/Riverside/", "B:/LA_Version2/WRF_FFDAS_GRIDDING/grid01/fishnet.Shp");
	//Preprocessor::gridFolderByShape("B:/LA_Version2/ca/Riverside/", "B:/LA_Version2/WRF_FFDAS_GRIDDING/grid02/Riverside/", "B:/LA_Version2/WRF_FFDAS_GRIDDING/grid02/fishnet.Shp");
	//Preprocessor::gridFolderByShape("B:/LA_Version2/ca/Riverside/", "B:/LA_Version2/WRF_FFDAS_GRIDDING/grid03/Riverside/", "B:/LA_Version2/WRF_FFDAS_GRIDDING/grid03/fishnet.Shp");

	//for (size_t i = 0; i < foldernames.size(); i++)
	//{
	//	ShapeFile::copyDirDropGeometry("C:/HestiaGridding/Los_Angeles/WRF_FFDAS_GRIDDING/grid01/"+ foldernames[i] + "/", "C:/HestiaGridding/Los_Angeles/WRF_FFDAS_GRIDDING/grid01_nogeometry/" + foldernames[i] + "/", fields2keep);
	//	ShapeFile::copyDirDropGeometry("C:/HestiaGridding/Los_Angeles/WRF_FFDAS_GRIDDING/grid02/" + foldernames[i] + "/", "C:/HestiaGridding/Los_Angeles/WRF_FFDAS_GRIDDING/grid02_nogeometry/" + foldernames[i] + "/", fields2keep);
	//	ShapeFile::copyDirDropGeometry("C:/HestiaGridding/Los_Angeles/WRF_FFDAS_GRIDDING/grid03/" + foldernames[i] + "/", "C:/HestiaGridding/Los_Angeles/WRF_FFDAS_GRIDDING/grid03_nogeometry/" + foldernames[i] + "/", fields2keep);
	//}

	//B:\LA_Version2\WRF_FFDAS_GRIDDING\grid01


	//sumNLCD("B:/LA_Version2/ComNonPoint_RI/ComNonPoint.Shp","COMSQFT");
	//sumNLCD("B:/LA_Version2/ComNonPoint_RI/IndNonPoint.Shp", "INDSQFT");
	//createTimeStructForComNonPoint("B:/LA_Version2/ComNonPoint_RI/ComNonPoint.Shp");
	//createTimeStructForIndNonPoint("B:/LA_Version2/ComNonPoint_RI/IndNonPoint.Shp");
	//GDAL_DS<unsigned char>* dsNLCD = new GDAL_DS<unsigned char>();
	//dsNLCD->open("B:/LA_Version2/Mapping/NLCD_Urban_LABasin_30m_percentage.tif");
	//OGREnvelope cropbound;
	//cropbound.MinX = 6659839.260581;
	//cropbound.MaxX = 7640995.917358;
	//cropbound.MinY = 1615202.045136;
	//cropbound.MaxY = 1870488.236334;
	//dsNLCD->crop(cropbound, "B:/LA_Version2/ComNonPoint_RI/NLCD.tif");

	//GDAL_DS<char>* dsMax = new GDAL_DS<char>();
	//dsMax->open("B:/LA_Version2/Mapping/NLCD_Urban_LABasin_30m_percentage.tif");
	//ExtractByPolygon<char>* extractor = new ExtractByPolygon<char>(dsMax);
	//extractor->extract("B:/LA_Version2/ComNonPoint_RI/Parcels.Shp", "nlcd_score");

	//Preprocessor::intersectWithArcGIS( "B:/LA_Version2/ComNonPoint_RI/Parcels.Shp","B:/LA_Version2/ComNonPoint_RI/fema_SQFT_RI.Shp", "B:/LA_Version2/ComNonPoint_RI/ParcelsByBlockgroup.Shp");
	//21 0 less than 20 percent of total cover
	//	22 35 20 - 49 percent of total cover
	//	23 65 50 - 79 percent of the total cover
	//	24 90 80 to 100 percent of the total cover

	//GDAL_DS<char>* dsMax = new GDAL_DS<char>();
	//dsMax->open("B:/LA_Version2/Mapping/NLCD_Urban_LABasin_30m.tif");
	//char nodata = dsMax->getNoData(1);
	//char* outlineDS = dsMax->readData(1);
	//char* pdata = outlineDS;

	//for (size_t i = 0; i < dsMax->slice; i++)
	//{
	//	char val = *pdata;
	//	if (val == nodata || val == 21)
	//		val = 0;
	//	else if (val == 22)
	//		val = 35;
	//	else if (val == 23)
	//		val = 65;
	//	else if (val == 24)
	//		val = 90;
	//	*pdata++ = val;
	//}
	//dsMax->create("B:/LA_Version2/Mapping/NLCD_Urban_LABasin_30m_percentage.tif");
	//dsMax->writeData(1, outlineDS, 0);
	//delete dsMax;

	/*GDAL_DS<double>* dsMax = new GDAL_DS<double>();
	dsMax->open("E:/GoogleStreetview/CardiffLiDAR/Mosaic.tif");*/

	//OGREnvelope newbb;
	//newbb.MinX = 316732.601214033;
	//newbb.MaxX = 320056.314149747;
	//newbb.MinY = 175258.342692511;
	//newbb.MaxY = 177680.081485819;
	//dsMax->crop(newbb, "E:/GoogleStreetview/CardiffLiDAR/subset.tif");

	//copy("B:/LA_Version2/ComNonPoint_RI/spatial/FEMA_shapfile/fema_bsf_2002bnd.Shp", "B:/LA_Version2/Fema_LA/fema_RI.Shp");

	//updateFEMA_SQFT("B:/LA_Version2/ComNonPoint_RI/fema_SQFT_RI.Shp");
	//updateFEMA_SQFT("B:/LA_Version2/ComNonPoint_RI/spatial/FEMA_shapfile/fema_bsf_2002bnd.Shp");
	//std::vector<std::string> fields;
	////fields.push_back("COMSQFT"); fields.push_back("INDSQFT"); fields.push_back("RESSQFT");
	//fields.push_back("yearbuilt"); fields.push_back("bc"); fields.push_back("toID");
	//fields.push_back("COMSQFT"); fields.push_back("INDSQFT"); fields.push_back("RESSQFT");
	//ShapeFile::copyDropGeometry("B:/LA_Version2/Fema_LA/parcels.Shp", "B:/LA_Version2/Fema_LA/parcels_trimmed.Shp", fields);

	//ExtractByPolygon<char>::link("B:/LA_Version2/Fema_LA/parcels.Shp", "B:/LA_Version2/Fema_LA/fema_RI.Shp", fields);

	//std::vector<std::string> fields;
	//fields.push_back("COMSQFT"); fields.push_back("INDSQFT"); fields.push_back("RESSQFT");



	//mergeCounties();

	//createPolygonGeometry();
	return 0;
	//reallocateLA();
	//reallocateHestiaFFDAS();
	//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithoutShippingAviation/grids/geo_em.d01.total.annual.2010.nc"));
	//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithoutShippingAviation/grids/geo_em.d02.total.annual.2010.nc"));
	//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithoutShippingAviation/grids/geo_em.d03.total.annual.2010.nc"));
	//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithoutShippingAviation/grids/geo_em.d01.total.hourly.2010.nc"));

	//geo_em.d03.total.annual.2010:   43.43 MtC


	//Grid grid;468284
	//QFile::copy("C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_SUBSET.tif", "C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_WithoutShippingAviation.tif");
	//removeShippingAviationFromFFDAS("C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_WithoutShippingAviation.tif", Preprocessor::readData("C:/HestiaGridding/Los_Angeles/WRF/2010/shipping_2010.tif"), Preprocessor::readData("C:/HestiaGridding/Los_Angeles/WRF/2010/aviation_2010.tif"));
	//gridLA_With_FFDAS();

	//std::ofstream ofs("geo_em.2010.csv");
	//ofs << "WithoutShippingAviation" << std::endl;
	//ofs << "geo_em.d03.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithoutShippingAviation/grids/geo_em.d03.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;
	//ofs << "geo_em.d02.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithoutShippingAviation/grids/geo_em.d02.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;
	//ofs << "geo_em.d01.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithoutShippingAviation/grids/geo_em.d01.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;

	//ofs << "WithShippingAviation" << std::endl;
	//ofs << "geo_em.d03.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithShippingAviation/grids/geo_em.d03.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;
	//ofs << "geo_em.d02.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithShippingAviation/grids/geo_em.d02.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;
	//ofs << "geo_em.d01.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithShippingAviation/grids/geo_em.d01.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;

	//ofs << "WithShippingOnly" << std::endl;
	//ofs << "geo_em.d03.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithShippingOnly/grids/geo_em.d03.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;
	//ofs << "geo_em.d02.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithShippingOnly/grids/geo_em.d02.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;
	//ofs << "geo_em.d01.total.annual.2010: " << checkTotal("C:/HestiaGridding/Los_Angeles/WRF/2010/WithShippingOnly/grids/geo_em.d01.total.annual.2010.nc") / 1000000000.0 << " MtC" << std::endl;

	//ofs.close();
	//QFile::copy("C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_SUBSET.tif", "C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_WithShippingOnly.tif");
	//removeAviationFromFFDAS("C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_WithShippingOnly.tif", Preprocessor::readData("C:/HestiaGridding/Los_Angeles/WRF/2010/aviation_2010.tif"));

	//QFile::copy("C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_SUBSET.tif", "C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_WithShippingAviation.tif");
	//grid.fromFishnetRaster()


}

//updateAttributes("B:/LA_Version2/Public/Los_Angeles/",    "B:/LA_Version2/Gridding/grid01/Los_Angeles/");
//updateAttributes("B:/LA_Version2/Public/Orange/",         "B:/LA_Version2/Gridding/grid01/Orange/");
//updateAttributes("B:/LA_Version2/Public/Riverside/",      "B:/LA_Version2/Gridding/grid01/Riverside/");
//updateAttributes("B:/LA_Version2/Public/San_Bernardino/", "B:/LA_Version2/Gridding/grid01/San_Bernardino/");
//updateAttributes("B:/LA_Version2/Public/Ventura/",        "B:/LA_Version2/Gridding/grid01/Ventura/");

//updateAttributes("B:/LA_Version2/Public/Los_Angeles/",    "B:/LA_Version2/Gridding/grid02/Los_Angeles/");
//updateAttributes("B:/LA_Version2/Public/Orange/",         "B:/LA_Version2/Gridding/grid02/Orange/");
//updateAttributes("B:/LA_Version2/Public/Riverside/",      "B:/LA_Version2/Gridding/grid02/Riverside/");
//updateAttributes("B:/LA_Version2/Public/San_Bernardino/", "B:/LA_Version2/Gridding/grid02/San_Bernardino/");
//updateAttributes("B:/LA_Version2/Public/Ventura/",        "B:/LA_Version2/Gridding/grid02/Ventura/");

//updateAttributes("B:/LA_Version2/Public/Los_Angeles/",    "B:/LA_Version2/Gridding/grid03/Los_Angeles/");
//updateAttributes("B:/LA_Version2/Public/Orange/",         "B:/LA_Version2/Gridding/grid03/Orange/");
//updateAttributes("B:/LA_Version2/Public/Riverside/",      "B:/LA_Version2/Gridding/grid03/Riverside/");
//updateAttributes("B:/LA_Version2/Public/San_Bernardino/", "B:/LA_Version2/Gridding/grid03/San_Bernardino/");
//updateAttributes("B:/LA_Version2/Public/Ventura/",        "B:/LA_Version2/Gridding/grid03/Ventura/");
//gridLA_With_FFDAS();
//gridLA_VY();
//gridLA_VY_Totals();


//processOnRoad();
//
//sumSector();
//updateLAResNonPointTime();


//scalebyfuel.scale("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/Scaling.csv", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/");
//scalebyfuel.scale("B:/LA_Version2/gridPrep_SHP_master/Orange/Scaling.csv", "B:/LA_Version2/gridPrep_SHP_master/Orange/");
//scalebyfuel.scale("B:/LA_Version2/gridPrep_SHP_master/Riverside/Scaling.csv", "B:/LA_Version2/gridPrep_SHP_master/Riverside/");
//scalebyfuel.scale("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/Scaling.csv", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/");
//scalebyfuel.scale("B:/LA_Version2/gridPrep_SHP_master/Ventura/Scaling.csv", "B:/LA_Version2/gridPrep_SHP_master/Ventura/");
//RemoteSensing rs;
//rs.processBaltimore();
//gridBaltimore();
//export2PUBLIC();
//gridLA();
//gridLA_With_FFDAS();


//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/WRF/WithoutShippingAviation/grids/geo_em.d01.total.annual.2014.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/WRF/WithoutShippingAviation/grids/geo_em.d02.total.annual.2014.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/WRF/WithoutShippingAviation/grids/geo_em.d03.total.annual.2014.nc"));


//double dif01 = checkTotal("C:/HestiaGridding/Los_Angeles/WRF//WithoutShippingAviation/grids/geo_em.d01.total.annual.2014.nc") - checkTotal("e:/geo_em.d01.total.annual.2014.nc");
//double dif02 = checkTotal("C:/HestiaGridding/Los_Angeles/WRF//WithoutShippingAviation/grids/geo_em.d02.total.annual.2014.nc") - checkTotal("e:/geo_em.d02.total.annual.2014.nc");
//double dif03 = checkTotal("C:/HestiaGridding/Los_Angeles/WRF//WithoutShippingAviation/grids/geo_em.d03.total.annual.2014.nc") - checkTotal("e:/geo_em.d03.total.annual.2014.nc");
//printf("%f,%f,%f\n", dif01/1000000/1000, dif02 / 1000000 / 1000, dif03 / 1000000 / 1000);


//printf("%f\n", checkTotal("C:/FFDAS/WRF/WithShippingAviation/grids/geo_em.d01.total.annual.2014.nc"));
//printf("%f\n", checkTotal("C:/FFDAS/WRF/WithShippingAviation/grids/geo_em.d02.total.annual.2014.nc"));
//printf("%f\n", checkTotal("C:/FFDAS/WRF/WithShippingAviation/grids/geo_em.d03.total.annual.2014.nc"));

//printf("%f\n", checkTotal("C:/FFDAS/WRF/WithShippingOnly/grids/geo_em.d01.total.annual.2014.nc"));
//printf("%f\n", checkTotal("C:/FFDAS/WRF/WithShippingOnly/grids/geo_em.d02.total.annual.2014.nc"));
//printf("%f\n", checkTotal("C:/FFDAS/WRF/WithShippingOnly/grids/geo_em.d03.total.annual.2014.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.annual.2010.v2.0.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.annual.2011.v2.0.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.annual.2012.v2.0.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.annual.2013.v2.0.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.annual.2014.v2.0.nc"));

//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.hourly.2011.v2.0.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.hourly.2014.v2.0.nc"));
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/LAbasin.total.annual.2011.v2.0.nc"));
//getchar();
//gridForEntropy();
/*SubdirManager dirmanager = SubdirManager("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/");

std::vector<std::string> onroadfiles = dirmanager.findFilesMatch(Utils::buildVector("", new std::string[1]{ "onroad" }, 1));

EMFAC emfac;
emfac.HPMS2SCAG("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/OnRoad.Shp", "B:/LA_Version2/OnRoad/Rural.csv", "B:/LA_Version2/OnRoad/Urban.csv");
emfac.HPMS2SCAG("B:/LA_Version2/gridPrep_SHP_master/Orange/OnRoad.Shp", "B:/LA_Version2/OnRoad/Rural.csv", "B:/LA_Version2/OnRoad/Urban.csv");
emfac.HPMS2SCAG("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/OnRoad.Shp", "B:/LA_Version2/OnRoad/Rural.csv", "B:/LA_Version2/OnRoad/Urban.csv");
emfac.HPMS2SCAG("B:/LA_Version2/gridPrep_SHP_master/Riverside/OnRoad.Shp", "B:/LA_Version2/OnRoad/Rural.csv", "B:/LA_Version2/OnRoad/Urban.csv");
emfac.HPMS2SCAG("B:/LA_Version2/gridPrep_SHP_master/Ventura/OnRoad.Shp", "B:/LA_Version2/OnRoad/Rural.csv", "B:/LA_Version2/OnRoad/Urban.csv");*/



//gridLA();
//gridBaltimore();
//C:\HestiaGridding\Baltimore\timestructs / timestructs_2010

//std::string indir = "C:/HestiaGridding/Baltimore/compress/";
//std::string tsid = "4001553";
//for (size_t iyear = 2010; iyear <= 2014; iyear++)
//{
//	std::stringstream ssin;
//	ssin << indir << "timestructs_" << iyear << ".bin";
//	std::stringstream ssout;
//	ssout << indir << iyear << "_" + tsid << ".txt";
//	std::stringstream txtindir;
//	txtindir << "B:/Baltimore/Points/Powerplants_Time/" << iyear << "/";
//	//TimestructTool::updateBinary(ssin.str(), ssin.str(), txtindir.str(), true);
//	TimestructTool::selectFromBinary2Text(ssin.str(), tsid, ssout.str());
//}

//gridMarionForFFDASODIAC();

//float val = 330.963184;
//float* outlineDS = Preprocessor::readData("B:/INDY.Total.2012.annual.nc");
//for (size_t i = 0; i < 521; i++)
//{
//	for (size_t j = 0; j < 560; j++)
//	{
//		/*if(!isnan(*outlineDS))
//		{
//			printf("%f\n", *outlineDS);
//		}*/
//		if (abs(*outlineDS - 330.963184) < 0.01)
//		{
//			printf("%d,%d,%f\n", i, j, *outlineDS);
//		}
//		outlineDS++;
//	}
//}

//for the period September 1, 2012 to April 30, 2013
//std::ofstream ofs;
//ofs.open("hestia_version2_totals.csv");
//int startday = QDate(2012, 1, 1).daysTo(QDate(2012, 9, 1));
//int numdays  = QDate(2012, 9, 1).daysTo(QDate(2012, 12, 31)) + 1;
////printf("2012=%f\n", calSum("C:/HestiaGridding/Indy/total_daily_2012.nc", startday, numdays));
////printf("2012=%f\n", calSum("H:/INDY.Total.2012.hourly.nc", startday * 24, numdays * 24));
//double sum1 = calSum("C:/HestiaGridding/Indy/total_daily_2012.nc", startday, numdays);
//startday = 0;
//numdays = QDate(2013, 1, 1).daysTo(QDate(2013, 4, 30)) + 1;
//double sum2 = calSum("C:/HestiaGridding/Indy/total_daily_2013.nc", startday, numdays);
//ofs << sum1 << "," << sum2 << "," << sum1 + sum2;
//ofs.close();
//printf("2013=%f\n", calSum("C:/HestiaGridding/Indy/total_daily_2013.nc", startday, numdays));
//printf("2013=%f\n", calSum("H:/INDY.Total.2013.hourly.nc", startday * 24, numdays * 24));
//printf("prior=%f\n", calSum("H:/INDY.Total.2013.hourly.nc", startday * 24, numdays * 24));
//printf("2012=%f\n", calSum("H:/INDY.Total.2012.hourly.nc", 0, 24*10));
//printf("2012=%f\n", calSum("C:/HestiaGridding/Indy/total_daily_2012.nc", 0, 10));



//std::ofstream ofs;
//ofs.open("prior_post_totals.csv");
//std::string indir = "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/";
//ofs << "Case" << "," << "prior" << "," << "prior" << std::endl;
//for (size_t i = 0; i < 4; i++)
//{
//	std::string casename;
//	std::stringstream sscase;
//	sscase << "CASE" << i;
//	casename = sscase.str();
//	std::string infilename = indir + casename + "_Inversion_Thomas.tif";
//	ofs << casename << "," << calSum(infilename.outlineDS(), 0, 1) << "," <<calSum(infilename.outlineDS(), 1, 1) << std::endl;
//}


//ofs.close();

//printf("CASE0:prior=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE0_Inversion_Thomas.tif", 0, 1));
//printf("CASE0:post=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE0_Inversion_Thomas.tif", 1, 1));

//printf("CASE1:prior=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE1_Inversion_Thomas.tif", 0, 1));
//printf("CASE1:post=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE1_Inversion_Thomas.tif", 1, 1));

//printf("CASE2:prior=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE2_Inversion_Thomas.tif", 0, 1));
//printf("CASE2:post=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE2_Inversion_Thomas.tif", 1, 1));

//printf("CASE3:prior=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE3_Inversion_Thomas.tif", 0, 1));
//printf("CASE3:post=%f\n", calSum("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE3_Inversion_Thomas.tif", 1, 1));



/*getchar();*/
//gridIndy();
//MarionElementa elementa;
//elementa.nc2tif("B:/Hestia_FFDAS_ODIAC/Influx/CASE0_Inversion_Thomas.nc", "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE0_Inversion_Thomas.tif");
//elementa.nc2tif("B:/Hestia_FFDAS_ODIAC/Influx/CASE1_Inversion_Thomas.nc", "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE1_Inversion_Thomas.tif");
//elementa.nc2tif("B:/Hestia_FFDAS_ODIAC/Influx/CASE2_Inversion_Thomas.nc", "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE2_Inversion_Thomas.tif");
//elementa.nc2tif("B:/Hestia_FFDAS_ODIAC/Influx/CASE3_Inversion_Thomas.nc", "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE3_Inversion_Thomas.tif");
//std::string maskfile = "B:/Hestia_FFDAS_ODIAC/Influx/MarionMask.tif";
//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE0_Inversion_Thomas.tif", maskfile, 0, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE0_outside.tif");
//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE0_Inversion_Thomas.tif", maskfile, 1, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE0_inside.tif");

//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE1_Inversion_Thomas.tif", maskfile, 0, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE1_outside.tif");
//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE1_Inversion_Thomas.tif", maskfile, 1, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE1_inside.tif");


//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE2_Inversion_Thomas.tif", maskfile, 0, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE2_outside.tif");
//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE2_Inversion_Thomas.tif", maskfile, 1, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE2_inside.tif");

//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE3_Inversion_Thomas.tif", maskfile, 0, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE3_outside.tif");
//elementa.crop("B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE3_Inversion_Thomas.tif", maskfile, 1, "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/CASE3_inside.tif");


//std::string indir = "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/";
//std::string datasets[] = {
//	"CASE0_outside","CASE0_inside",
// "CASE1_outside","CASE1_inside",
// "CASE2_outside","CASE2_inside",
// "CASE3_outside","CASE3_inside"
//};
//int numdatasets = 8;
//for (size_t i = 0; i < numdatasets; i++)
//{
//	std::string infile = indir + datasets[i] + ".tif";
//	std::string outfile = indir + datasets[i] + "_dif.tif";
//	std::string outcsffile = indir + datasets[i] + "_dif.csv";
//	elementa.dif(infile, outfile, outcsffile);
//	//std::string lidarFisheyeDir = indir + datasets[i];
//	//if (!QDir(lidarFisheyeDir.outlineDS()).exists())
//	//{
//	//	QDir(lidarFisheyeDir.outlineDS()).mkpath(".");
//	//}
//}
//std::string indir = "B:/Hestia_FFDAS_ODIAC/Influx/prior_post/";
//std::string datasets[] = {
//	"CASE0_Inversion_Thomas","CASE1_Inversion_Thomas",
// "CASE2_Inversion_Thomas","CASE3_Inversion_Thomas"};
//int numdatasets = 4;
//for (size_t i = 0; i < numdatasets; i++)
//{
//	std::string infile = indir + datasets[i] + ".tif";
//	std::string outfile = indir + datasets[i] + "_dif.tif";
//	std::string outcsffile = indir + datasets[i] + "_dif.csv";
//	elementa.dif(infile, outfile, outcsffile);
//	//std::string lidarFisheyeDir = indir + datasets[i];
//	//if (!QDir(lidarFisheyeDir.outlineDS()).exists())
//	//{
//	//	QDir(lidarFisheyeDir.outlineDS()).mkpath(".");
//	//}
//}
//plotRaster();
//plotShapes();
//updateOnRoadFields();
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/compress/LABasin.total.hourly.2010.v2.0.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/compress/LABasin.total.annual.2010.v2.0.nc"));
//getchar();
//gridLA();





//std::vector<std::string> fields2drop= Utils::buildVector("", new std::string[2]{ "AREA" ,"LENGTH" }, 2);
//ShapeFile Shp("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/NonRoad_510.Shp");
//int fieldidx = Shp.poLayer->GetLayerDefn()->GetFieldIndex("area");
//printf("%d\n", fieldidx);
//std::vector<std::string> subdirs = Utils::findSubdirectories(indir);
//for (size_t i = 0; i < subdirs.size(); i++)
//{
//	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
//	std::string srcdir = subdirs[i];
//	std::string destdir = "B:/LA_Version2/all/" + foldername + "/";
//	ShapeFile::copyDirDropFields(subdirs[i], "B:/LA_Version2/all/" + foldername + "/", fields2drop);
//	Utils::updateFootprintForDir(destdir, true);
//}
//
//std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
//std::vector<std::string> subdirs = Utils::findSubdirectories(indir);
//std::vector<std::string> fields2keep = Utils::buildVector("", new std::string[4]{ "ca11","area","length","timestruct" }, 4);
//for (size_t i = 0; i < subdirs.size(); i++)
//{
//	std::string foldername = QDir(subdirs[i].outlineDS()).dirName().toLocal8Bit().outlineDS();
//	ShapeFile::copyDir(subdirs[i], "B:/LA_Version2/ca11/" + foldername + "/", fields2keep);
//	//Utils::updateFootprintForDir(subdirs[i], true);
//}
//printf("%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/2014/total_annual_2014.nc"));

//gridLA_With_FFDAS();

//getchar();

//Utils::updateFootprint("B:/LA_Version2/Gridding/FFDASGrid.Shp");
//Utils::updateFootprint("B:/LA_Version2/ca/FFDAS/ffdas.Shp");
//OGREnvelope worldbound;
//worldbound.MaxX = 180; worldbound.MinX = -180; worldbound.MaxY = 90; worldbound.MinY = -90;
//Grid grid(worldbound,0.1,0);

////grid.toShape(&refshape, "B:/LA_Version2/Gridding/FFDASGrid.Shp");
//int startrow = 437;
//int nrows = 212;
//int startcol = 468;
//int ncols = 284;
//double* adfTransform = grid._adfGeoTransform;
//adfTransform[0] = -180 + startcol * 0.1;
//adfTransform[3] = 90 - startrow * 0.1;
//adfTransform[1] = 0.1;
//adfTransform[5] = -0.1;
/*Preprocessor::extractEDGASHourly("C:/FFDAS/aviation/2010.nc", "C:/HestiaGridding/Los_Angeles/WRF/2010/aviation_2010.tif", adfTransform, 437, 212, 468, 284);
Preprocessor::extractEDGASHourly("C:/FFDAS/shipping/2010.nc", "C:/HestiaGridding/Los_Angeles/WRF/2010/shipping_2010.tif", adfTransform, 437, 212, 468, 284);
Preprocessor::extractFFDASHourly("C:/FFDAS/2010/", "C:/HestiaGridding/Los_Angeles/WRF/2010/FFDAS_SUBSET.tif", adfTransform, 437, 212, 468, 284);*/


//gridBaltimore();
//export2PUBLIC();
//testTimeshift();
//gridBaltimore();
//gridLA();
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/2011/total_annual_2011.nc"), checkTotal("C:/HestiaGridding/Los_Angeles/2011/total_hourly_2011.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/2012/total_annual_2012.nc"), checkTotal("C:/HestiaGridding/Los_Angeles/2012/total_hourly_2012.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/2013/total_annual_2013.nc"), checkTotal("C:/HestiaGridding/Los_Angeles/2013/total_hourly_2013.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Los_Angeles/2014/total_annual_2014.nc"), checkTotal("C:/HestiaGridding/Los_Angeles/2014/total_hourly_2014.nc"));

//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Baltimore/total_annual_2010.nc"), checkTotal("C:/HestiaGridding/Baltimore/total_hourly_2010.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Baltimore/total_annual_2011.nc"), checkTotal("C:/HestiaGridding/Baltimore/total_hourly_2011.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Baltimore/total_annual_2012.nc"), checkTotal("C:/HestiaGridding/Baltimore/total_hourly_2012.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Baltimore/total_annual_2013.nc"), checkTotal("C:/HestiaGridding/Baltimore/total_hourly_2013.nc"));
//printf("%f,%f\n", checkTotal("C:/HestiaGridding/Baltimore/total_annual_2014.nc"), checkTotal("C:/HestiaGridding/Baltimore/total_hourly_2014.nc"));





//getchar();
//gridLA();
//double onroad_hourly = calSum("C:/Baltimore/200m/onroad_hourly.nc");
//double onroad_annual = calSum("C:/Baltimore/200m/onroad_annual.nc");

//printf("%f\n%f\n", onroad_annual, onroad_hourly);
//createHourlyProfiles("B:/Baltimore/OnRoadTimeStructureCreation/OnRoadWGS84.Shp","B:/Baltimore/OnRoadTimeStructureCreation/Raster_Extents/","B:/Baltimore/OnRoadTimeStructureCreation/Timefiles/WeekDays/");
//Preprocessor::reprojectDir("B:/Baltimore/gridPrep_SHP_master/StatePlane", "B:/Baltimore/gridPrep_SHP_master/reprojected/", "B:/Hestia_FFDAS_ODIAC/Reproject/ReprojectBaltimoreToFeet.py");
/*calAverageDistance("B:/SpatialGranuality/UrbanArea/Baltimore/ComNonPoint.Shp");
calAverageDistance("B:/SpatialGranuality/UrbanArea/Los_Angeles/ComNonPoint.Shp");
calAverageDistance("B:/SpatialGranuality/UrbanArea/Marion/ComNonPoint.Shp");
calAverageDistance("B:/SpatialGranuality/UrbanArea/SaltLake/ComNonPoint.Shp");*/
//Utils::updateFootprintForDir("B:/Baltimore/gridPrep_SHP_master/reprojected/",true);
//Utils::updateFootprint("B:/Baltimore/gridPrep_SHP_master/StatePlane/onroad.Shp");
//TemporalGridder gridder2;
//int nrows = 50;
//int ncols = 100;
//int times = 24;
//int ncount = nrows * ncols * times;
//double* outlineDS = new double[ncount];
//double* pdata = outlineDS;
//double* xcoords = new double[ncols];
//double* ycoords = new double[nrows];
//for (int i = 0; i < ncols; i++)
//{
//	xcoords[i] = i;
//}
//for (int i = 0; i < nrows; i++)
//{
//	ycoords[i] = nrows-i-1;
//}

//for (int ntime = 0; ntime< times; ntime++)
//{
//	for (int nrow = 0; nrow < nrows; nrow++)
//	{
//		for (int ncol = 0; ncol < ncols; ncol++)
//		{
//			*pdata++ = ntime;
//		}
//	}
//}
//NCFile* nc = gridder2.createNC("e:/timegrid.nc", ncols, nrows, times);
//printf("%f\n", outlineDS[1000]);
//double* slicedata = new double[nrows*ncols];
//for (int i = 0; i < nrows*ncols; i++)
//{
//	slicedata[i] = i;
//}
//nc->writeSlice(2, slicedata);
////nc->write(3, outlineDS);
//nc->write(0, xcoords);
//nc->write(1, ycoords);
//nc->write(2, xcoords);
//return 0;
//Grid grid;
//std::vector<std::string> shapefiles;
//shapefiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/nonroad.Shp");
//shapefiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/onroad.Shp");
//shapefiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/resnonpoint.Shp");
//shapefiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/comnonpoint.Shp");
//shapefiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/indnonpoint.Shp");
//shapefiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/CMVPort.Shp");
//shapefiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/CMVUnderway.Shp");
//BoundManager::writeBound(BoundManager::readBoundFromShapes(shapefiles), "B:/Baltimore/gridPrep_SHP_master/StatePlane/bound.txt");


//Grid grid(BoundManager::readBound("B:/Baltimore/gridPrep_SHP_master/StatePlane/bound.txt"), resol * 3.28084,1);
//ShapeFile refshape("B:/Baltimore/gridPrep_SHP_master/StatePlane/nonroad.Shp");
//grid.toShape(&refshape, "B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".Shp");
//grid.reset();
//grid.toRaster("B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");
//Preprocessor::gridFolderByRaster("B:/Baltimore/gridPrep_SHP_master/ca/", "B:/Baltimore/gridPrep_SHP_master/ca11/" + resolname.str() + "/", "B:/Baltimore/gridPrep_SHP_master/"+ fishnetname.str() + ".tif");

/*std::vector<std::string> fields;
fields.push_back("ca11");
fields.push_back("length");
fields.push_back("area");
fields.push_back("timestruct");
ShapeFile::copyDir("B:/Baltimore/gridPrep_SHP_master/StatePlane/", "B:/Baltimore/gridPrep_SHP_master/ca11/", fields);*/
//Utils::updateFootprintForDir("B:/Baltimore/gridPrep_SHP_master/ca11/");
//grid.fromFishnetRaster("B:/LA_Version2/Gridding/geo_em.d01_VAR.tif",true);
//grid.resetValue("ca", 1.0);
//grid.toShape(NULL, "B:/LA_Version2/Gridding/geo_em.d01.fishnet.Shp", true);
//grid.fromFishnetRaster("B:/LA_Version2/Gridding/geo_em.d02_VAR.tif", true);
//grid.resetValue("ca", 1.0);
//grid.toShape(NULL, "B:/LA_Version2/Gridding/geo_em.d02.fishnet.Shp", true);
//grid.fromFishnetRaster("B:/LA_Version2/Gridding/geo_em.d03_VAR.tif", true);
//grid.resetValue("ca", 1.0);
//grid.toShape(NULL, "B:/LA_Version2/Gridding/geo_em.d03.fishnet.Shp", true);

/*grid.fromFishnetRaster("B:/LA_Version2/Gridding/geo_em.d01_VAR.tif",true);
grid.toBoundaryShape("B:/LA_Version2/Gridding/geo_em.d01.bound.Shp", "B:/LA_Version2/Gridding/geo_em.d01.fishnet.Shp");
grid.fromFishnetRaster("B:/LA_Version2/Gridding/geo_em.d02_VAR.tif", true);
grid.toBoundaryShape("B:/LA_Version2/Gridding/geo_em.d02.bound.Shp", "B:/LA_Version2/Gridding/geo_em.d02.fishnet.Shp");
grid.fromFishnetRaster("B:/LA_Version2/Gridding/geo_em.d03_VAR.tif", true);
grid.toBoundaryShape("B:/LA_Version2/Gridding/geo_em.d03.bound.Shp", "B:/LA_Version2/Gridding/geo_em.d03.fishnet.Shp");
return 0;*/
//"Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Post-1979.csv"
//"Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Pre-1980.csv"
//"Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Post-1979.csv"
//"Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Pre-1980.csv"
//std::vector<std::string> CBECS12NEEUI_files;
//CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Pre-1980.csv");
//CBECS12NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/CBECS12NEEUI_Post-1979.csv");

//std::vector<std::string> RECS09NEEUI_files;
//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Pre-1980.csv");
//RECS09NEEUI_files.push_back("Z:/Hestia/Precalculated NEEUI/RECS09NEEUI_Post-1979.csv");


//std::vector<int> breaks;
//breaks.push_back(1980);
//int censusDivision = 5;
//NonPointProcessor::updateNEEUI("B:/Baltimore/gridPrep_SHP_master/StatePlane/ComNonPoint.Shp", CBECS12NEEUI_files, breaks, censusDivision);
//NonPointProcessor::updateNEEUI("B:/Baltimore/gridPrep_SHP_master/StatePlane/ResNonPoint.Shp", RECS09NEEUI_files, breaks, censusDivision);


//ScaleByFuel scaleByFuel;
//scaleByFuel.scaleNonRoad("B:/LA_Version2/gridPrep_SHP_master","B:/LA_Version2/Processing/NonRoad.csv");
//RoadTimeIDW roadidw;
//roadidw.m_maxDist = 10000 * 3.28084;
//roadidw.m_minnumpoints = 0;
//roadidw.loadStations("B:/LA_Version2/OnRoadTime/BinaryStationCounts/","E:/LA_OnRoad_Time/");
//roadidw.coordinatesFromShapefile("B:/LA_Version2/OnRoadTime/BinaryStationCounts/stations.Shp","");
//roadidw.createIDW();
////roadidw.stations2shapes("B:/LA_Version2/OnRoadTime/BinaryStationCounts/stations.Shp");
//roadidw.idw("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/OnRoad.Shp", "BT", "timestruct", "5", "6037", 0);
//roadidw.idw("B:/LA_Version2/gridPrep_SHP_master/Orange/OnRoad.Shp", "BT", "timestruct", "5", "6059", 0);
//roadidw.idw("B:/LA_Version2/gridPrep_SHP_master/Riverside/OnRoad.Shp", "BT", "timestruct", "5", "6065", 0);
//roadidw.idw("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/OnRoad.Shp", "BT", "timestruct", "5", "6071", 0);
//roadidw.idw("B:/LA_Version2/gridPrep_SHP_master/Ventura/OnRoad.Shp", "BT", "timestruct", "5", "6111", 0);

//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/", 6037);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Orange/", 6059);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Riverside/", 6065);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", 6071);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Ventura/", 6111);



//copyBuildingFiles(LA_DIRS,"B:/LA_Version2/VintageExperiments/NoVintage/");
//int fips[] = { 6037,6059,6065,6071,6111 };
//Los_Angeles 6037
//Orange 6059
//Riverside 6065
//San_Bernardino 6071
//Ventura 6111
//createTimeStructForNonRoad(LA_DIRS, fips);
//createTimeStruct(LA_DIRS);
//processDirs(LA_DIRS);
//gridAllFolders("B:/LA_Version2/ca11/", 1000 * 3.28084);
//gridDirs(LA_DIRS);

//copyFolderDropFields(LA_DIRS);

//calNonRoadTotals(LA_DIRS, "nonroadpoint.csv");
//caTotals(LA_DIRS, "ca", "createRailroadPoint.Shp","railroadpoint.csv");
//caTotals(LA_DIRS, "ca", "createRailroadPoint.Shp", "railroadpoint.csv");
//copyField(LA_DIRS, "ca11", "ca");
//return 0;
//ShapeCreator shpCreator;
//shpCreator.create("B:/LA_Version2/Processing/Los_Angeles/Powerplants.txt", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/ElecProd.Shp");
//shpCreator.create("B:/LA_Version2/Processing/Orange/Powerplants.txt", "B:/LA_Version2/gridPrep_SHP_master/Orange/ElecProd.Shp");
//shpCreator.create("B:/LA_Version2/Processing/Riverside/Powerplants.txt", "B:/LA_Version2/gridPrep_SHP_master/Riverside/ElecProd.Shp");
//shpCreator.create("B:/LA_Version2/Processing/San_Bernardino/Powerplants.txt", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/ElecProd.Shp");
//shpCreator.create("B:/LA_Version2/Processing/Ventura/Powerplants.txt", "B:/LA_Version2/gridPrep_SHP_master/Ventura/ElecProd.Shp");
//linkSCC("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/ElecProd.Shp", "B:/LA_Version2/Processing/NEI_PointSources/LA_basin_ElecProd.Shp");
//linkSCC("B:/LA_Version2/gridPrep_SHP_master/Orange/ElecProd.Shp", "B:/LA_Version2/Processing/NEI_PointSources/LA_basin_ElecProd.Shp");
//linkSCC("B:/LA_Version2/gridPrep_SHP_master/Riverside/ElecProd.Shp", "B:/LA_Version2/Processing/NEI_PointSources/LA_basin_ElecProd.Shp");
//linkSCC("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/ElecProd.Shp", "B:/LA_Version2/Processing/NEI_PointSources/LA_basin_ElecProd.Shp");
//linkSCC("B:/LA_Version2/gridPrep_SHP_master/Ventura/ElecProd.Shp", "B:/LA_Version2/Processing/NEI_PointSources/LA_basin_ElecProd.Shp");

//////////NonPointProcessor npProcessor;
//////////npProcessor.exportBySector("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/Parcels_06037_LA.Shp", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/");
//////////npProcessor.exportBySector("B:/LA_Version2/gridPrep_SHP_master/Orange/Parcels_06059_OR.Shp", "B:/LA_Version2/gridPrep_SHP_master/Orange/");
//////////npProcessor.exportBySector("B:/LA_Version2/gridPrep_SHP_master/Riverside/Parcels_06065_RI.Shp", "B:/LA_Version2/gridPrep_SHP_master/Riverside/");
//////////npProcessor.exportBySector("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/Parcels_06071_SB.Shp", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/");
//////////npProcessor.exportBySector("B:/LA_Version2/gridPrep_SHP_master/Ventura/Parcels_06111_VE.Shp", "B:/LA_Version2/gridPrep_SHP_master/Ventura/");
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/", 6037);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Orange/", 6059);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Riverside/", 6065);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", 6071);
//createBySector("B:/LA_Version2/Processing/NEI_PointSources/LA_basin_NonElecProd.Shp", "B:/LA_Version2/gridPrep_SHP_master/Ventura/", 6111);

//std::string nonpoint_RSCRIPT = "B:/LA_Version2/Parcel_Carbon_Allocation/SpatialAllocation.R";
//std::string sharedCFG = "B:/LA_Version2/Parcel_Carbon_Allocation/shared.csv";
//calculateNonPointWeight("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/", nonpoint_RSCRIPT, sharedCFG);
//calculateNonPointWeight("B:/LA_Version2/gridPrep_SHP_master/Orange/", nonpoint_RSCRIPT, sharedCFG);
//calculateNonPointWeight("B:/LA_Version2/gridPrep_SHP_master/Riverside/", nonpoint_RSCRIPT, sharedCFG);
//calculateNonPointWeight("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", nonpoint_RSCRIPT, sharedCFG);
//calculateNonPointWeight("B:/LA_Version2/gridPrep_SHP_master/Ventura/", nonpoint_RSCRIPT, sharedCFG);
//tC2011
//std::vector<ShapeCreator::SECTORID> sectors;
//sectors.push_back(ShapeCreator::IND); sectors.push_back(ShapeCreator::COM); sectors.push_back(ShapeCreator::RAILROAD); sectors.push_back(ShapeCreator::AIRPORT); sectors.push_back(ShapeCreator::NONROAD);
//createShapesForSectors("B:/LA_Version2/Processing/Los_Angeles/NEI.csv", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/",sectors);
//createShapesForSectors("B:/LA_Version2/Processing/Orange/NEI.csv", "B:/LA_Version2/gridPrep_SHP_master/Orange/", sectors);
//createShapesForSectors("B:/LA_Version2/Processing/Riverside/NEI.csv", "B:/LA_Version2/gridPrep_SHP_master/Riverside/", sectors);
//createShapesForSectors("B:/LA_Version2/Processing/San_Bernardino/NEI.csv", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", sectors);
//createShapesForSectors("B:/LA_Version2/Processing/Ventura/NEI.csv", "B:/LA_Version2/gridPrep_SHP_master/Ventura/", sectors);



//gridFolderByRaster("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/", "B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/grid/", 500 * 3.28084);
//gridFolderByRaster("B:/LA_Version2/gridPrep_SHP_master/Orange/", "B:/LA_Version2/gridPrep_SHP_master/Orange/grid/", 500 * 3.28084);
//gridFolderByRaster("B:/LA_Version2/gridPrep_SHP_master/Riverside/", "B:/LA_Version2/gridPrep_SHP_master/Riverside/grid/", 500 * 3.28084);
//gridFolderByRaster("B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/", "B:/LA_Version2/gridPrep_SHP_master/San_Bernardino/grid/", 500 * 3.28084);
//gridFolderByRaster("B:/LA_Version2/gridPrep_SHP_master/Ventura/", "B:/LA_Version2/gridPrep_SHP_master/Ventura/grid/", 500 * 3.28084);

// Preprocessor::gridShapeFile("B:/Shapefiles2Grid/OSM_ONROAD/Baltimore_OnRoad_OSM.Shp", "B:/Shapefiles2Grid/OSM_ONROAD/Baltimore_OnRoad_OSM_fishnet.Shp", "Z:/Hestia/BALTIMORE/gridPrep_SHP_intersect/fishnet.Shp","C11");
//Preprocessor::gridShapeFile("B:/Shapefiles2Grid/OSM_ONROAD/Marion_Onroad_OSM.Shp", "B:/Shapefiles2Grid/OSM_ONROAD/Marion_Onroad_OSM_fishnet.Shp", "Z:/Hestia/INDIANAPOLIS/INDIANAPOLIS_NEI2011/gridPrep_SHP_intersect/fishnet.Shp","C11");

//grid.fromFishnetRaster("Z:/Hestia/BALTIMORE/gridPrep_SHP_intersect/fishnet.tif");
//Preprocessor::intersectWithFishnet
// Preprocessor::convertPolyline("E:/Hestia_Onroad_Subsets/src/All_Roads_Baltimore.Shp", "E:/Hestia_Onroad_Subsets/src/Baltimore/OnRoad.Shp");
//Preprocessor::convertPolyline("E:/Hestia_Onroad_Subsets/src/All_Roads_Marion.Shp", "E:/Hestia_Onroad_Subsets/src/Marion/OnRoad.Shp");
//Preprocessor::convertPolyline("E:/Hestia_Onroad_Subsets/src/All_Roads_SLC.Shp", "E:/Hestia_Onroad_Subsets/src/SaltLake/OnRoad.Shp");


//std::string onroadfields[] = {"ca02","ca10","ca11","ca12","ca13","ca14" };

//double onroadtotalsBaltimore[] = { 0,	512420524.7,	501115185.2,	489019464.1,	478998999.6,	468816182 };//Baltimore
//double onroadtotalsMarion[] = { 0,	1743373570,	1700826830,	1702462348,	1718010167,	1745519149 };//Marion
//double onroadtotalsSLC[] = { 1080345357,	1363071468,	1361110800,	1471271909,	0,	0 };//SaltLake
//reallocateOnRoad("E:/Hestia_Onroad_Subsets/src/Baltimore/OnRoad.Shp", onroadfields, onroadtotalsBaltimore, 6, "CO2_tC");
//reallocateOnRoad("E:/Hestia_Onroad_Subsets/src/Marion/OnRoad.Shp", onroadfields, onroadtotalsMarion, 6, "CO2_tC");
//reallocateOnRoad("E:/Hestia_Onroad_Subsets/src/SaltLake/OnRoad.Shp", onroadfields, onroadtotalsSLC, 6, "CO2_tC");

//Utils::updateFootprint("B:/SpatialGranuality/Baltimore/reprejected/OnRoad.Shp");

//   std::vector<std::string> fields; fields.push_back("ca11");
//ShapeFile::copy("B:/LA_Version2/gridPrep_SHP_master/Los_Angeles/OnRoad.Shp", "B:/SpatialGranuality/Los_Angeles/OnRoad.Shp", fields);
//Utils::updateFootprint("B:/SpatialGranuality/Los_Angeles/OnRoad.Shp");
//Utils::updateFootprint("E:/Hestia_Onroad_Subsets/Baltimore/OnRoad.Shp");
//Utils::updateFootprint("E:/Hestia_Onroad_Subsets/Marion/OnRoad.Shp");
//Utils::updateFootprint("E:/Hestia_Onroad_Subsets/SaltLake/OnRoad.Shp");

//   ShapeFile::copy("E:/Hestia_Onroad_Subsets/Baltimore/OnRoad.Shp", "B:/SpatialGranuality/Baltimore/OnRoad.Shp", fields);
//ShapeFile::copy("E:/Hestia_Onroad_Subsets/Marion/OnRoad.Shp", "B:/SpatialGranuality/Marion/OnRoad.Shp", fields);
//ShapeFile::copy("E:/Hestia_Onroad_Subsets/SaltLake/OnRoad.Shp", "B:/SpatialGranuality/SaltLake/OnRoad.Shp", fields);

//   City	2002	2010	2011	2012	2013	2014
//Baltimore	0	512420524.7	501115185.2	489019464.1	478998999.6	468816182
//Marion	0	1743373570	1700826830	1702462348	1718010167	1745519149
//SaltLake	1080345357	1363071468	1361110800	1471271909	0	0

//DbfFile_c dbf("B:/Gridded_Cities/INDIANAPOLIS/0/10m/OnRoad.dbf");
////dbf.DumpAll("B:/Gridded_Cities/INDIANAPOLIS/0/OnRoad.txt");
////dbf.DumpAll("B:/Gridded_Cities/INDIANAPOLIS/0/OnRoad.txt");
//std::vector<const char*> fields2;
//std::vector<std::string> fields;
//fields.push_back("Id");
//fields.push_back("ca11");

//fields2.push_back("Id");
//fields2.push_back("ca11");
////dbf.DumpFields("test.txt",&fields2[0],2);
//std::vector<std::vector<std::string>> table = dbf.ReadFields(fields);
//for (size_t i = 0; i < table.size(); i++)
//{
//	std::vector<std::string>& row = table[i];
//	int id = atoi(row[0].outlineDS());
//	double ca = atof(row[1].outlineDS());
//	printf("%d,%d,%f\n", i, id, ca);
//}
//return 0;
//unit is in foot
//BoundManager::writeBound(BoundManager::readBoundFromShape("B:/Gridded_Cities/SaltLake/NonRoad.Shp"),"B:/Gridded_Cities/SaltLake/bound.txt");
//splitIntoTiles("B:/Gridded_Cities/SaltLake", "B:/Gridded_Cities/SaltLake", BoundManager::readBound("B:/Gridded_Cities/SaltLake/bound.txt"),16000*3.28084);
//intersectByTiles("B:/Gridded_Cities/SaltLake", "B:/Gridded_Cities/SaltLake", BoundManager::readBound("B:/Gridded_Cities/SaltLake/bound.txt"), 16000 * 3.28084, 32.8084);
//gridShapeFile("B:/Baltimore/Gridding/Gridding/circle.Shp", "B:/Baltimore/Gridding/Gridding/fishnet_hehe.Shp", 0.1);
//BoundManager::writeBound(BoundManager::readBoundFromShape("B:/Gridded_Cities/INDIANAPOLIS/NonRoad.Shp"), "B:/Gridded_Cities/INDIANAPOLIS/bound.txt");
//updateFootprintForDir("B:/Gridded_Cities/INDIANAPOLIS");
//splitIntoTiles("B:/Gridded_Cities/INDIANAPOLIS", "B:/Gridded_Cities/INDIANAPOLIS", BoundManager::readBound("B:/Gridded_Cities/INDIANAPOLIS/bound.txt"),16000*3.28084);
//intersectByTiles("B:/Gridded_Cities/INDIANAPOLIS", "B:/Gridded_Cities/INDIANAPOLIS", BoundManager::readBound("B:/Gridded_Cities/INDIANAPOLIS/bound.txt"), 16000 * 3.28084, 32.8084);
//gridFolderByRaster("B:/Gridded_Cities/INDIANAPOLIS", "B:/Gridded_Cities/INDIANAPOLIS/10", 32.8084);
//updateFootprintForDir("B:/Gridded_Cities/Los_Angeles/", true);
//intersectWithPolygonForDir("B:/Gridded_Cities/LA_Boundary_Continental.Shp", "B:/Gridded_Cities/LA_Basin/", "B:/Gridded_Cities/LA_County/");
//OGREnvelope bound = BoundManager::readBound("B:/Baltimore/gridPrep_SHP_master/StatePlane/bound.txt");
//Grid grid(bound, 200, 1);
//ShapeFile ref("B:/Baltimore/gridPrep_SHP_master/StatePlane/NonRoad.Shp");
//grid.toShape(&ref,"B:/Baltimore/Hestia_gridding/input/Baltimore/Grids/fishnet.Shp");
//gridFolderByRaster("B:/Baltimore/gridPrep_SHP_master/StatePlane/", "B:/Baltimore/gridPrep_SHP_master/StatePlane/200_time/", 200);
//processFolderForGrids("B:/Baltimore/gridPrep_SHP_master/StatePlane/", "B:/Baltimore/gridPrep_SHP_master/StatePlane/200_time/", "B:/Baltimore/gridPrep_SHP_master/StatePlane/200_time/grids/", 200);
//upscaleFishnetBySector("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/", 20, "B:/Baltimore/gridPrep_SHP_master/StatePlane/200/");
//updateFieldAfterIntersectionForDir("B:/Gridded_Cities/LA_County/");
//updateFootprintForDir("B:/Gridded_Cities/LA_County/trimmed/");

//BoundManager::writeBound(BoundManager::readBoundFromShape("B:/Gridded_Cities/LA_County/NonRoad.Shp"),"B:/Gridded_Cities/LA_County/bound.txt");
//splitIntoTiles("B:/Gridded_Cities/LA_County", "B:/Gridded_Cities/LA_County", BoundManager::readBound("B:/Gridded_Cities/LA_County/bound.txt"),16000*3.28084);
//intersectByTiles("B:/Gridded_Cities/LA_County", "B:/Gridded_Cities/LA_County", BoundManager::readBound("B:/Gridded_Cities/LA_County/bound.txt"), 16000 * 3.28084, 32.8084);

//std::vector<std::string> dirs;
//dirs.push_back("B:/Gridded_Cities/LA_County/");
//dirs.push_back("B:/Gridded_Cities/INDIANAPOLIS/");
//dirs.push_back("B:/Gridded_Cities/SaltLake/");
//for (size_t i = 0; i < dirs.size(); i++)
//{
//
//	std::string boundfile = dirs[i] + "bound.txt";
//	mergeTiles(dirs[i], dirs[i] + "merged/", BoundManager::readBound(boundfile), 16000 * 3.28084, "10m");
//	//std::vector<std::string> subdirs = findChildDirctories(dirs[i], BoundManager::readBound(boundfile), 16000 * 3.28084,"10m");
//	//for (size_t j = 0; j < subdirs.size(); j++)
//	//{
//	//	printf("%s\n", subdirs[j].outlineDS());
//	//	std::stringstream ss;
//	//	OGREnvelope bound = BoundManager::readBoundFromShape(subdirs[j] + ".Shp");
//	//	updateFieldAfterIntersectionForDir(subdirs[j]);
//	//	updateFishnet(subdirs[j].outlineDS(), bound, 32.8084);
//	//	
//	//}
//}
//return 0;
//std::string indir = argv[1];
//std::string lidarFisheyeDir = argv[2];
//double resol = atof(argv[3]);
//bool skipNonRoad = false;
//if (argc > 4)
//{
//	if (argv[4] == "true" || atoi(argv[4]) == 1)
//		skipNonRoad = true;
//}

//gridFolderByRaster("B:/Baltimore/gridPrep_SHP_master/StatePlane", "B:/Baltimore/gridPrep_SHP_master/StatePlane/500", 500);
//gridFolderByRaster(indir, lidarFisheyeDir, resol, skipNonRoad);

//OGREnvelope bound = BoundManager::readBoundFromShape("B:/Baltimore/gridPrep_SHP_master/StatePlane/NonRoad.Shp");
////splitIntoTiles("B:/Baltimore/gridPrep_SHP_master/StatePlane", "B:/Baltimore/gridPrep_SHP_master/StatePlane/tiles", bound,1000);
//intersectByTiles("B:/Baltimore/gridPrep_SHP_master/StatePlane", "B:/Baltimore/gridPrep_SHP_master/StatePlane/tiles", bound, 1000,10);

//sumField("B:/Baltimore/gridPrep_SHP_master/WGS84/result");
//verifyRailRoad();

//updateFishnet("B:/Baltimore/gridPrep_SHP_master/StatePlane/10");
//upscaleFishnet("B:/Baltimore/gridPrep_SHP_master/StatePlane", "B:/Baltimore/gridPrep_SHP_master/StatePlane/10",2,"ca11");
//updateFishnet("B:/Baltimore/gridPrep_SHP_master/StatePlane/10", "B:/Baltimore/gridPrep_SHP_master/StatePlane");
//upscaleFishnet("B:/Baltimore/gridPrep_SHP_master/StatePlane/10",10,"ca11", 
//"B:/Baltimore/gridPrep_SHP_master/StatePlane//10/100/fishnet.Shp", "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/100/fishnet.tif");
//BoundManager::writeBound("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/fishnet.Shp", "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/bound.txt");
//BoundManager::writeBound("B:/Baltimore/gridPrep_SHP_master/StatePlane/NonRoad.Shp", "B:/Baltimore/gridPrep_SHP_master/StatePlane/bound.txt");
//BoundManager::writeBound("B:/Baltimore/gridPrep_SHP_master/WGS84/NonRoad.Shp", "B:/Baltimore/gridPrep_SHP_master/WGS84/bound.txt");

//int scales[] = {1,2,4,5,8,10,20,40,50,100};
//int num = sizeof(scales) / sizeof(int);
//for (size_t i = 0; i < num; i++)
//{
//	int scale = scales[i];
//	std::stringstream ssShapeFile, ssGridFile;

//	std::string dir = "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/";
//	ssShapeFile << dir  << "fishnet/" << scale * 10 << ".Shp";
//	ssGridFile <<  dir  << "fishnet/" << scale * 10 << ".tif";
//	if (QFileInfo(ssShapeFile.str().outlineDS()).exists())
//		continue;
//	printf("%s\n", ssShapeFile.str().outlineDS());
//	upscaleFishnet("B:/Baltimore/gridPrep_SHP_master/StatePlane/10", scale, "ca11",
//		ssShapeFile.str(), ssGridFile.str());
//}
//int scales[] = {1};
//int num = sizeof(scales) / sizeof(int);
//for (size_t i = 0; i < num; i++)
//{
//	int scale = scales[i];
//	std::stringstream lidarFisheyeDir;

//	std::string dir = "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/";
//	lidarFisheyeDir << dir  << "sectors" << scale * 10 << "/";
//	upscaleFishnetBySector("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/", scale, lidarFisheyeDir.str());
//	
//}



//computeRMSE("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/","B:/Baltimore/gridPrep_SHP_master/StatePlane/10/rmse.csv");

/*std::vector<std::string> sectorFiles;
sectorFiles.clear();
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/OnRoad.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ComNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ResNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/IndNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/Railroad.Shp");
computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/rmse_all.csv");

sectorFiles.clear();
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/OnRoad.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ComNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ResNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/IndNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/Railroad.Shp");
computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/rmse_all_n.csv", true);

sectorFiles.clear();
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/OnRoad.Shp");
computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/rmse_onroad.csv");

sectorFiles.clear();
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ComNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ResNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/IndNonPoint.Shp");
computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/rmse_building.csv");




sectorFiles.clear();
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/OnRoad.Shp");
computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/rmse_onroad_n.csv",true);

sectorFiles.clear();
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ComNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ResNonPoint.Shp");
sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/IndNonPoint.Shp");
computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/rmse_building_n.csv", true);*/


//std::vector<std::string> sectorFiles = loadAllSectorFile("B:/Baltimore/gridPrep_SHP_master/StatePlane/10");
//computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/RMSE_ALL.csv");

//sectorFiles.clear();
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ComNonPoint.Shp");
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ResNonPoint.Shp");
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/IndNonPoint.Shp");
//computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/RMSE_NonPoint.csv");

//sectorFiles.clear();
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/OnRoad.Shp");
//computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/RMSE_OnRoad.csv");

//sectorFiles.clear();
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/Railroad.Shp");
//computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/RMSE_Railroad.csv");

//sectorFiles.clear();
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ComPoint.Shp");
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/IndPoint.Shp");
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/ElecProd.Shp");
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/createRailroadPoint.Shp");
//sectorFiles.push_back("B:/Baltimore/gridPrep_SHP_master/StatePlane/10/Airport.Shp");

//computeRMSE_BySector(sectorFiles, "B:/Baltimore/gridPrep_SHP_master/StatePlane/10/RMSE_Point.csv");


// 	return 0;
//
//
//
//}
