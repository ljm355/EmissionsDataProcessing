// PolygonStatistics.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "ogrsf_frmts.h"
//#include "gdal_priv.h"
#include <QFileinfo>
#include <map>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <qdir.h>
#include "qfileinfo.h"
#include "DbfFile.h"
#include "osgDB\ReadFile"
#include "osgDB\WriteFile"
#include "osg\Image"


void mapFolder(std::string indir1, std::string indir2, std::string outdir, int numClasses = 32)
{
	if (!QDir(indir1.data()).exists() || !QDir(indir2.data()).exists())
		return;
	osg::ref_ptr<osg::Image> colorRamp = osgDB::readImageFile("B:/FFDAS/color_ramps/colorramp.bmp");
	QDir qoutdir(outdir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	std::vector<std::string> files;

	QDir input_dir1(indir1.data());
	input_dir1.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir1.setSorting(QDir::Name);
	indir1 = (input_dir1.absolutePath() + "/").toLocal8Bit().data();

	indir2 = (QDir(indir2.data()).absolutePath() + "/").toLocal8Bit().data();

	QFileInfoList list = input_dir1.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		std::string input_file = (input_dir1.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".tif"))
			continue;
		files.push_back(fileInfo.baseName().toLocal8Bit().data());
	}



}

void reprojectWithArcGIS(std::string inputfile, std::string outputfile, std::string pythonScriptFile)
{
	QFileInfo info(outputfile.data());

	std::ifstream pyifs(pythonScriptFile.data());
	pyifs.seekg(0, std::ios::end);
	size_t size = pyifs.tellg();
	std::string buffer(size, ' ');
	pyifs.seekg(0);
	pyifs.read(&buffer[0], size);

	int inputStart = buffer.find("INPUT", 0);
	std::string pyscript = buffer.substr(0, inputStart) + inputfile + buffer.substr(inputStart + 5, buffer.size() - inputStart + 5);
	int outputStart = pyscript.find("OUTPUT", 0);
	pyscript = pyscript.substr(0, outputStart) + outputfile + pyscript.substr(outputStart + 6, pyscript.size() - outputStart + 6);

	std::ofstream ofs;
	std::string scriptFile = (info.absoluteDir().absolutePath() + "/" + info.completeBaseName() + ".py").toLocal8Bit().data();
	ofs.open(scriptFile.data());
	ofs << pyscript;
	ofs.close();
	system(scriptFile.data());

	QFile::remove(scriptFile.data());
}

void move(std::string src, std::string dest)
{
	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
		pszDriverName);
	poDriver->Delete(dest.data());
	GDALDataset* poDSSrc = (GDALDataset*)GDALOpenEx(src.data(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	GDALDataset* poDSDest = (GDALDataset*)poDriver->CreateCopy(dest.data(), poDSSrc, 0, 0, 0, 0);
	GDALClose(poDSSrc);
	GDALClose(poDSDest);
	poDriver->Delete(src.data());
}
void readText2Buf(std::string filename, char*& buf, int& length)
{
	std::ifstream is(filename.data(), std::ifstream::binary);
	is.seekg(0, is.end);
	length = is.tellg();
	is.seekg(0, is.beg);

	buf = new char[length];
	// read data as a block:
	is.read(buf, length);

	is.close();
}
void reprojectWithOGR(std::string inputfile, std::string outputfile, std::string destwktfile)
{


	GDALDataset* poDS = (GDALDataset*)GDALOpenEx(inputfile.data(), GDAL_OF_VECTOR, NULL, NULL, NULL);
	OGRLayer* poLayer = poDS->GetLayer(0);
	int srcEPSG = 4326;
	char* buf = NULL;
	int length = 0;
	readText2Buf(destwktfile, buf, length);
	OGRSpatialReference oSRS;
	oSRS.importFromWkt(&buf);

	bool isWGS84 = false;
	std::string srcwktfile = QFileInfo(destwktfile.data()).absoluteDir().absolutePath().toLocal8Bit().data() + std::string("/srcwkt.txt");
	
	if (!poLayer->GetSpatialRef() || poLayer->GetSpatialRef()->GetEPSGGeogCS() == srcEPSG)
	{
		//OGRSpatialReference oSRS;
		//oSRS.SetWellKnownGeogCS("WGS84");
		printf("%s: spatial reference not found.\n", inputfile.data());
		isWGS84 = true;
	}
	else if (poLayer->GetSpatialRef()->IsSame(&oSRS))
	{
		//move(inputfile, outputfile);
		return;
	}
	else
	{
		static char* srcwktbuf = new char[10000];
		poLayer->GetSpatialRef()->exportToWkt(&srcwktbuf);
		//QFileInfo fileinfo(inputfile.)
		std::string srcwktstr = srcwktbuf;
		std::ofstream ofs;
		ofs.open(srcwktfile);
		ofs << srcwktstr;
		ofs.close();
		//delete[] srcwktbuf;
	}
	GDALClose(poDS);
	//static char* wktsrc = new char[100000];
	//static char* wktdest = new char[100000];

	//dSRS->exportToWkt(&wktdest);
	//oSRS.exportToWkt(&wktsrc);



	std::stringstream ss;
	if (isWGS84)
	{
		ss << "ogr2ogr -f" << " " << "\"" << "ESRI Shapefile" << "\"" << " " << outputfile << " " << inputfile
			<< " " << "-s_srs" << " " << "EPSG:" << srcEPSG
			<< " " << "-t_srs" << " " << destwktfile << "\n";
	}
	else
	{
		ss << "ogr2ogr -f" << " " << "\"" << "ESRI Shapefile" << "\"" << " " << outputfile << " " << inputfile
			<< " " << "-s_srs" << " " << srcwktfile
			<< " " << "-t_srs" << " " << destwktfile << "\n";
	
	}
	//std::stringstream ss;
	//ss << "ogr2ogr -f" << " " << "\"" << "ESRI Shapefile" << "\"" << " " << outputfile << " " << inputfile
	//	<< " " << "-s_srs" << " " << "EPSG:" << srcEPSG
	//	<< " " << "-t_srs" << " " << "EPSG:" << destEPSG << "\n";


	//"ogr2ogr -f "ESRI Shapefile" original.shp wgs84.shp -s_srs EPSG:27700 -t_srs EPSG:4326"
	printf("reprojecting: %s.\n", inputfile.data());
	std::string commandline = ss.str().data();
	system(commandline.data());
	QFile::remove(srcwktfile.data());
}
void reprojectDir(std::string indir,std::string outdir, std::string destwktfile)
{
	if (!QDir(indir.data()).exists())
		return;
	
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
		std::string input_file = fileInfo.absoluteFilePath().toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string output_file = outdir + fileInfo.fileName().toLocal8Bit().data();
		//reprojectWithOGR(input_file, output_file, pythonFile);
		reprojectWithOGR(input_file, output_file, destwktfile);
		printf("%s\n", output_file.data());
	}

}


int main(int argc, char** argv)
{
	OGRRegisterAll();
	GDALAllRegister();
	std::string destwktfile = "B:/LA_Version2/gridPrep_SHP_master/NAD_1983_StatePlane_California_V_FIPS_0405_Feet.txt";
	reprojectDir("B:/LA_Version2/gridPrep_SHP_master/LA/", "B:/LA_Version2/gridPrep_SHP_master/LA/reprojected", destwktfile);
	reprojectDir("B:/LA_Version2/gridPrep_SHP_master/OR/", "B:/LA_Version2/gridPrep_SHP_master/OR/reprojected", destwktfile);
	reprojectDir("B:/LA_Version2/gridPrep_SHP_master/RI/", "B:/LA_Version2/gridPrep_SHP_master/RI/reprojected", destwktfile);
	reprojectDir("B:/LA_Version2/gridPrep_SHP_master/SB/", "B:/LA_Version2/gridPrep_SHP_master/SB/reprojected", destwktfile);
	reprojectDir("B:/LA_Version2/gridPrep_SHP_master/VE/", "B:/LA_Version2/gridPrep_SHP_master/VE/reprojected", destwktfile);
	return 0;
}
