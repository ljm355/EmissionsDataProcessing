#include "RemoteSensing.h"
#include "gdal_priv.h"
#include "Utils.h"
#include "QDirIterator"
#include <map>
#include <sstream>
#include "GDAL_DS.h"
#include "TemporalGridding.h"
#include "Grid.h"
#include "Preprocessor.h"
RemoteSensing::RemoteSensing()
{
}


RemoteSensing::~RemoteSensing()
{

}

std::string RemoteSensing::getProj(std::string filename)
{
	if (filename.substr(filename.length() - 3, 3) == "shp")
	{
		ShapeFile shp(filename);
		char wkt[512];
		char* pwkt = wkt;
		if (shp.poLayer->GetSpatialRef())
		{
			shp.poLayer->GetSpatialRef()->exportToWkt(&pwkt);
			return std::string(pwkt);
		}
		return "";
	}

	GDALDataset* pDataset = (GDALDataset*)GDALOpen(filename.data(), GA_ReadOnly);
	std::string proj = pDataset->GetProjectionRef();
	GDALClose(pDataset);

	return proj;
}

OGRSpatialReference RemoteSensing::create(std::string wkt)
{
	OGRSpatialReference spatialref;
	char* proj = (char*)wkt.data();
	OGRErr err = spatialref.importFromWkt(&proj);
	return spatialref;
}

OGRSpatialReference RemoteSensing::getWGS84()
{
	OGRSpatialReference oSRS;
	oSRS.SetWellKnownGeogCS("WGS84");
	//oSRS.SetWellKnownGeogCS("EPSG:4326")
	return oSRS;
}
void RemoteSensing::writeBound(OGREnvelope bound, OGRSpatialReference* srs, std::string outshapefile)
{
	ShapeFile shp;
	shp.create(outshapefile.data(), srs);
	OGRFeature* poFeaPolygon = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
	OGRPolygon *poPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(wkbPolygon);
	OGRLinearRing  *linearRing = (OGRLinearRing  *)OGRGeometryFactory::createGeometry(wkbLinearRing);
	linearRing->addPoint(bound.MinX, bound.MinY);
	linearRing->addPoint(bound.MinX, bound.MaxY);
	linearRing->addPoint(bound.MaxX, bound.MaxY);
	linearRing->addPoint(bound.MaxX, bound.MinY);
	linearRing->addPoint(bound.MinX, bound.MinY);
	poPolygon->addRing(linearRing);//also crashed

	poFeaPolygon->SetGeometry(poPolygon);
	shp.poLayer->CreateFeature(poFeaPolygon);
	OGRFeature::DestroyFeature(poFeaPolygon);
	shp.close();
}
OGREnvelope RemoteSensing::readBound(std::string filename, OGRSpatialReference * targetSRS)
{
	ShapeFile shp(filename);
	OGREnvelope env;
    shp.poLayer->GetExtent(&env);
	env = transformBound(env, shp.poLayer->GetSpatialRef(), targetSRS);
	return env;
}

OGREnvelope RemoteSensing::transformBound(OGREnvelope srcBound, OGRSpatialReference* sourceSRS, OGRSpatialReference* targetSRS)
{
	OGRCoordinateTransformation* poCT = OGRCreateCoordinateTransformation(sourceSRS,targetSRS);
	poCT->Transform(1, &(srcBound.MinX), &(srcBound.MaxY));
	poCT->Transform(1, &(srcBound.MaxX), &(srcBound.MinY));
	OCTDestroyCoordinateTransformation(poCT);
	return srcBound;
}

//void RemoteSensing::extractUrbanBound(const char* boundShapefile, const char * rasterDir, const char * outroot)
//{
//	double total = 0;
//	ShapeFile shp(boundShapefile);
//	OGRFeature *poFeature;
//	shp.poLayer->ResetReading();
//	std::vector<OGREnvelope> polygonbounds;
//	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
//	{
//		OGREnvelope env;
//		poFeature->GetGeometryRef()->getEnvelope(&env);
//		polygonbounds.push_back(env);
//		OGRFeature::DestroyFeature(poFeature);
//	}
//	shp.close();
//	std::vector<std::string> rasterfiles;
//	QDirIterator it(rasterDir, QDirIterator::Subdirectories);
//	while (it.hasNext()) {
//		std::vector<std::string> subset = Utils::findFiles(it.next().toLocal8Bit().data(),".tif");
//		for (size_t i = 0; i < subset.size(); i++)
//		{
//			rasterfiles.push_back(subset[i]);
//		}
//	}
//
//	std::map<int, std::vector<std::string>> bounddict;
//	double adftransform[6];
//	for (int i = 0; i < rasterfiles.size(); i++)
//	{
//		GDALDataset* pDataset = (GDALDataset*)GDALOpen(rasterfiles[i].data(), GA_ReadOnly);
//		pDataset->GetGeoTransform(adftransform);
//		int ncols = pDataset->GetRasterXSize();
//		int nrows = pDataset->GetRasterYSize();
//		OGREnvelope rasterbound;
//		rasterbound.MinX = adftransform[0];
//		rasterbound.MaxY = adftransform[3];
//		rasterbound.MaxX = adftransform[0] + adftransform[1] * ncols;
//		rasterbound.MinY = adftransform[3] + adftransform[5] * nrows;
//		GDALClose(pDataset);
//		for (int npolygon = 0; npolygon < polygonbounds.size(); npolygon++)
//		{
//			OGREnvelope polygonbound = polygonbounds[npolygon];
//			if (!polygonbound.Intersects(rasterbound))
//				continue;
//
//			std::map<int, std::vector<std::string>>::iterator iter = bounddict.find(npolygon);
//			if (iter == bounddict.end())
//			{
//				std::vector<std::string> interectedRasterBounds;
//				interectedRasterBounds.push_back(rasterfiles[i]);
//				bounddict[npolygon] = interectedRasterBounds;
//			}
//			else
//			{
//				iter->second.push_back(rasterfiles[i]);
//			}
//		}
//	}
//
//
//
//}

void RemoteSensing::extractUrbanBound(OGREnvelope targetBound, const char * rasterDir, const char* outdir)
{

	std::vector<std::string> rasterfiles;
	QDirIterator it(rasterDir, QDirIterator::Subdirectories);
	while (it.hasNext()) {
		std::vector<std::string> subset = Utils::findFiles(it.next().toLocal8Bit().data(), ".tif");
		for (size_t i = 0; i < subset.size(); i++)
		{
			rasterfiles.push_back(subset[i]);
		}
	}
	QDir qoutdir(outdir);
	if (!qoutdir.exists())
		qoutdir.mkpath(".");

	std::vector<std::string> interectedRasterBounds;
	double adftransform[6];
	QString qoutdirstring = qoutdir.absolutePath() + "/";
	QString qoutdirname = QDir(outdir).dirName();
	QString buildvrt_command = "gdalbuildvrt " + qoutdirstring + qoutdirname + ".vrt " + qoutdirstring + "*.tif";
	printf("%s\n", buildvrt_command.toLocal8Bit().data());
	//OGRSpatialReference destSRS;// = create(getProj("E:/UrbanExtent/GHSL/BETA/FULL/MT/12/0/1664_mt.tif"));
	//destSRS.SetWellKnownGeogCS("EPSG:3857");
	//OGRSpatialReference destSRS = create(getProj("E:/UrbanExtent/GHSL/BETA/FULL/MT/1664_mt.tif"));
	//char wkt[512];
	//char* pwkt = wkt;
	//destSRS.exportToWkt(&pwkt);
	std::string wkt = getProj("E:/UrbanExtent/GHSL/BETA/FULL/MT/1664_mt.tif");
	for (int i = 0; i < rasterfiles.size(); i++)
	{
		GDALDataset* pDataset = (GDALDataset*)GDALOpen(rasterfiles[i].data(), GA_ReadOnly);
		pDataset->GetGeoTransform(adftransform);
		int ncols = pDataset->GetRasterXSize();
		int nrows = pDataset->GetRasterYSize();
		OGREnvelope rasterbound;
		rasterbound.MinX = adftransform[0];
		rasterbound.MaxY = adftransform[3];
		rasterbound.MaxX = adftransform[0] + adftransform[1] * ncols;
		rasterbound.MinY = adftransform[3] + adftransform[5] * nrows;
		GDALClose(pDataset);	
		if (i % 100 == 0)
			printf("%d/%d\n", i, rasterfiles.size()-1);
		if (!targetBound.Intersects(rasterbound) && !targetBound.Contains(rasterbound) && !rasterbound.Contains(targetBound))
			continue;
		std::vector<std::string> interectedRasterBounds;
		interectedRasterBounds.push_back(rasterfiles[i]);
		std::stringstream ss;
		ss << qoutdir.absolutePath().toLocal8Bit().data() << "/" << i << "." << QFileInfo(rasterfiles[i].data()).completeSuffix().toLocal8Bit().data();
		QString destfile = ss.str().data();// qoutdir.absolutePath() + "/" + QFileInfo(rasterfiles[i].data()).fileName();
		QFile::copy(rasterfiles[i].data(), destfile);

		pDataset = (GDALDataset*)GDALOpen(destfile.toLocal8Bit().data(), GA_Update);
		pDataset->SetProjection(wkt.data());
		GDALClose(pDataset);

	}
	system(buildvrt_command.toLocal8Bit().data());
}

//18011,4720124,"Boone"
//18057,4487263,"Hamilton"
//18059,3423507,"Hancock"
//18063,4558700,"Hendricks"
//18081,3588222,"Johnson"
//18095,5050192,"Madison"
//18097,4418912,"Marion"
//18109,4565755,"Morgan"
//18145,4602585,"Shelby"

struct CountyBuiltup
{
	unsigned short fips;
	unsigned int numlandpixels;
	unsigned int numlandpixels_check;
	unsigned int numbuiltuppixels;
	std::string name;
	CountyBuiltup(unsigned int _fips, unsigned int _numlandpixels, std::string _name,unsigned int _numbuiltuppixels = 0)
	{
		fips = _fips;
		numlandpixels = _numlandpixels;
		numbuiltuppixels = _numbuiltuppixels;
		name = _name;
		numlandpixels_check = 0;
	}

};
#include <fstream>
void countyStats()
{
	std::vector<CountyBuiltup> allcounties;
	allcounties.push_back(CountyBuiltup(18011, 4720124, "Boone"));
	allcounties.push_back(CountyBuiltup(18057, 4487263, "Hamilton"));
	allcounties.push_back(CountyBuiltup(18059, 3423507, "Hancock"));
	allcounties.push_back(CountyBuiltup(18063, 4558700, "Hendricks"));
	allcounties.push_back(CountyBuiltup(18081, 3588222, "Johnson"));
	allcounties.push_back(CountyBuiltup(18095, 5050192, "Madison"));
	allcounties.push_back(CountyBuiltup(18097, 4418912, "Marion"));
	allcounties.push_back(CountyBuiltup(18109, 4565755, "Morgan"));
	allcounties.push_back(CountyBuiltup(18145, 4602585, "Shelby"));


	CountyBuiltup* allcounties_arr[20000];
	for (int ncounty = 0; ncounty < allcounties.size(); ncounty++)
	{
		allcounties_arr[allcounties[ncounty].fips] = &(allcounties[ncounty]);
	}

	GDALDataset* pGHSLDataset = (GDALDataset*)GDALOpen("C:/MOD11A2.005_seasonal/Indianapolis/GHSL_ALL.tif", GA_ReadOnly);
	int ncols = pGHSLDataset->GetRasterXSize();
	int nrows = pGHSLDataset->GetRasterYSize();
	char* builtupdata = new char[ncols*nrows];
	pGHSLDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, ncols, nrows, builtupdata, ncols, nrows, GDT_Byte, 0, 0);
	GDALClose(pGHSLDataset);

	GDALDataset* pCountyDataset = (GDALDataset*)GDALOpen("C:/MOD11A2.005_seasonal/Indianapolis/MarionAndAll_Bound.tif", GA_ReadOnly);
	unsigned short* countymaskdata = new unsigned short[ncols*nrows];
	pCountyDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, ncols, nrows, countymaskdata, ncols, nrows, GDT_UInt16, 0, 0);
	GDALClose(pCountyDataset);

	int numcells = nrows*ncols;

	for (int ncell = 0; ncell < numcells; ncell++)
	{
		unsigned short countymask = countymaskdata[ncell];
		if (countymask == 0)
			continue;
		CountyBuiltup* county = allcounties_arr[countymask];
		county->numlandpixels_check++;
		char builtuptype = builtupdata[ncell];
		if (builtuptype < 3)
			continue;
		county->numbuiltuppixels++;
	}

	std::ofstream ofs("e:/indy_builtup.csv");
	ofs << "county,fips,numlandpixels,numlandpixels_check,numbuiltuppixels,numbuiltuppixels/numlandpixels" << std::endl;

	for (int ncounty = 0; ncounty < allcounties.size(); ncounty++)
	{
		CountyBuiltup& county = allcounties[ncounty];
		ofs << county.name << "," << county.fips << "," << county.numlandpixels << "," << county.numlandpixels_check << "," << county.numbuiltuppixels << "," << (double)county.numbuiltuppixels / (double)county.numlandpixels_check << std::endl;
	}

	ofs.close();
	delete[] countymaskdata;
	delete[] builtupdata;

}

//Bands
//Wavelength
//Resolution
//Band 1 每 Coastal aerosol
//0.43 每 0.45
//30
//Band 2 每 Blue
//0.45 每 0.51
//30
//Band 3 每 Green
//0.53 每 0.59
//30
//Band 4 每 Red
//0.64 每 0.67
//30
//Band 5 每 Near Infrared(NIR)
//0.85 每 0.88
//30
//Band 6 每 SWIR 1
//1.57 每 1.65
//30
//Band 7 每 SWIR 2
//2.11 每 2.29
//30
//Band 8 每 Panchromatic
//0.50 每 0.68
//15
//Band 9 每 Cirrus
//1.36 每 1.38
//30
//Band 10 每 Thermal Infrared
//(TIRS) 1
//10.60 每 11.19
//100 * (30)
//Band 11 每 Thermal Infrared
//(TIRS) 2
//11.50 每 12.51
//100 * (30)

std::vector<std::string> RemoteSensing::constructShapeFilelist(std::string dir, std::string* filearr, int numoffiles)
{
	std::vector<std::string> files;
	for (size_t i = 0; i < numoffiles; i++)
	{
		files.push_back(dir + filearr[i] + ".shp");
	}
	return files;
}

void RemoteSensing::gridBaltimore()
{
	std::string rootdir = "E:/LanduseRegressionBaltimore/";

	int resol = 30;
	std::stringstream ssoutroot;
	ssoutroot << rootdir << "Hestia" << resol << "m" <<"/";
	std::string outroot = ssoutroot.str();

	QDir(outroot.data()).mkpath(".");
	
	std::string fishnetfile = outroot + "fishnet.shp";
	std::string reffileshp = "E:/LanduseRegressionBaltimore/bound.shp";
	std::string reffiletif = "E:/LanduseRegressionBaltimore/bound.tif";
	std::string proj = getProj(reffiletif);
	Grid grid;
	grid.fromFishnetRaster(reffiletif);
	//grid.toShape(proj, fishnetfile);
	//grid.reset();
	//grid.toRaster("B:/Baltimore/gridPrep_SHP_master/" + fishnetname.str() + ".tif");
	//Utils::updateFootprintForDir("B:/Baltimore/gridPrep_SHP_master/ca/");
	std::string indir = outroot + "intersected/";
	//Preprocessor::gridFolderByRaster("B:/Baltimore/gridPrep_SHP_master/ca_1.1/", indir, reffiletif);
	//E:/LanduseRegressionBaltimore/Landsat8/Baltimore_0.tif
	std::string years[]{"2014" };
	std::string cafields[]{"ca14" };
	std::string cityname = "Baltimore";
	std::string version = "";
	//Baltimore.total.hourly.2010.v1.1.nc


	std::string year = years[0];
	TemporalGridder gridder(outroot,8760);
	gridder.fromFishnetRaster(reffiletif);
	std::vector<std::string> comshapefiles = constructShapeFilelist(indir, new std::string[2]{ "comnonpoint","compoint" }, 2);
	std::vector<std::string> comnonpointshapefiles = constructShapeFilelist(indir, new std::string[2]{ "comnonpoint"}, 1);
	std::vector<std::string> indshapefiles = constructShapeFilelist(indir, new std::string[2]{ "indnonpoint","indpoint" }, 2);
	std::vector<std::string> indnonpointshapefiles = constructShapeFilelist(indir, new std::string[2]{ "indnonpoint"}, 1);
	std::vector<std::string> resshapefiles = constructShapeFilelist(indir, new std::string[1]{ "resnonpoint" }, 1);
	std::vector<std::string> railroadshapefiles = constructShapeFilelist(indir, new std::string[2]{ "RailroadPoint","Railroad" }, 2);
	std::vector<std::string> onroadshapefiles = constructShapeFilelist(indir, new std::string[1]{ "OnRoad" }, 1);
	std::vector<std::string> nonroadshapefiles = constructShapeFilelist(indir, new std::string[1]{ "NonRoad" }, 1);
	std::vector<std::string> elecprodshapefiles = constructShapeFilelist(indir, new std::string[1]{ "ElecProd" }, 1);
	//gridder.addSectorGrid(comshapefiles, "com");
	gridder.addSectorGrid(comnonpointshapefiles, "comnonpoint");
	//gridder.addSectorGrid(indshapefiles, "ind");
	gridder.addSectorGrid(indnonpointshapefiles, "indnonpoint");
	gridder.addSectorGrid(resshapefiles, "res");
	gridder.addSectorGrid(railroadshapefiles, "railroad");
	gridder.addSectorGrid(onroadshapefiles, "onroad");
	//gridder.addSectorGrid(nonroadshapefiles, "nonroad");
	//gridder.addSectorGrid(elecprodshapefiles, "elecprod");
	gridder.loadAttribute(cafields[0]);
	std::string gridfiletif = outroot + "total_nonpoint.tif";

	gridder.getTotal();

	GDAL_DS<double>* ds = new GDAL_DS<double>();
	gridder.toShapefile(reffileshp, outroot + "total_nonpoint.shp");
	ds->open(reffiletif);
	ds->create(gridfiletif);
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

void RemoteSensing::processBaltimore()
{

	//gridBaltimore();
	//return;
	//processLandsat("E:/LanduseRegressionBaltimore/Landsat8/", "E:/LanduseRegressionBaltimore/Landsat8/de-clouded/");
	processLandsat("E:/LanduseRegressionBaltimore/NewYork/", "E:/LanduseRegressionBaltimore/NewYork/de-clouded/");

}
void RemoteSensing::test()
{
	countyStats();
	//OGRSpatialReference sourceSRS = RemoteSensing::getWGS84();
	//OGRSpatialReference destSRS = create(getProj("E:/UrbanExtent/GHSL/BETA/FULL/MT/1664_mt.tif"));
	////destSRS.SetWellKnownGeogCS("EPSG:3857");
	////EPSG 3857
	//OGREnvelope LA_BOUND;
	////LA_BOUND.MinX = -119.68896007;
	////LA_BOUND.MaxX = -114.058285362;
	////LA_BOUND.MinY = 32.689068972;
	////LA_BOUND.MaxY = 35.832632771;
	//LA_BOUND.MinX = -13323714.088562895;
	//LA_BOUND.MaxX = -12696910.247228336;
	//LA_BOUND.MinY = 3854105.42089237;
	//LA_BOUND.MaxY = 4277616.270804314;
	////LA_BOUND.MinX = -13505704.550869996;
	////LA_BOUND.MaxX = -12286190.021100197;
	////LA_BOUND.MinY = 3708396.653806514;
	////LA_BOUND.MaxY = 4441082.285844444;
	//
	////LA_BOUND = transformBound(LA_BOUND, &sourceSRS, &destSRS);

	//extractUrbanBound(LA_BOUND, "E:/UrbanExtent/GHSL/BETA/FULL/MT/12/", "C:/MOD11A2.005_seasonal/LosAngeles/GHSL/");
	//OGRSpatialReference sourceSRS = RemoteSensing::getWGS84();


	OGRSpatialReference destSRS = create(getProj("E:/UrbanExtent/GHSL/BETA/FULL/MT/1664_mt.tif"));
	//OGREnvelope outerbound = readBound("B:/Indianapolis/Grids/g2d_WGS84.shp", &destSRS);
	//writeBound(outerbound, &destSRS, "C:/MOD11A2.005_seasonal/Indianapolis/outerbound.shp");

	//OGREnvelope innerbound = readBound("B:/Hestia_FFDAS_ODIAC/Marion_Influx/fishnet.shp", &destSRS);
	//writeBound(innerbound, &destSRS, "C:/MOD11A2.005_seasonal/Indianapolis/innerbound.shp");

	//destSRS.SetWellKnownGeogCS("EPSG:3857");
	//EPSG 3857
	OGREnvelope LA_BOUND;
	//LA_BOUND.MinX = -119.68896007;
	//LA_BOUND.MaxX = -114.058285362;
	//LA_BOUND.MinY = 32.689068972;
	//LA_BOUND.MaxY = 35.832632771;
	LA_BOUND.MinX = -9650843.254322851;
	LA_BOUND.MaxX = -9526165.424634675;
	LA_BOUND.MinY = 4770174.537267407;
	LA_BOUND.MaxY = 4921288.105826897;
	//LA_BOUND.MinX = -13505704.550869996;
	//LA_BOUND.MaxX = -12286190.021100197;
	//LA_BOUND.MinY = 3708396.653806514;
	//LA_BOUND.MaxY = 4441082.285844444;

	//LA_BOUND = transformBound(LA_BOUND, &sourceSRS, &destSRS);

	LA_BOUND.MinX = 13115705.616129244;
	LA_BOUND.MaxX = 13852307.089332189;
	LA_BOUND.MinY = 3401660.329833442;
	LA_BOUND.MaxY = 3877911.28233534;

	LA_BOUND = readBound("B:/Baltimore/gridPrep_SHP_master/fishnet200m.shp", &destSRS);
	extractUrbanBound(LA_BOUND, "E:/UrbanExtent/GHSL/BETA/FULL/MT/12/", "C:/MOD11A2.005_seasonal/Baltimore/GHSL/");

	LA_BOUND = readBound("B:/SpatialGranuality/SaltLake/fishnet.shp", &destSRS);
	extractUrbanBound(LA_BOUND, "E:/UrbanExtent/GHSL/BETA/FULL/MT/12/", "C:/MOD11A2.005_seasonal/SaltLake/GHSL/");


	
}

void RemoteSensing::processLandsat(std::string indir, std::string outdir)
{
	std::vector<std::string> files = Utils::findFiles(indir, ".tif");
	QDir(outdir.data()).mkpath(".");
	for (size_t i = 0; i < files.size(); i++)
	{
		std::string name = QFileInfo(files[i].data()).completeBaseName().toLocal8Bit().data();
		GDAL_DS<float>* ds = new GDAL_DS<float>();
		ds->open(files[i].data(), GA_ReadOnly);
		std::vector<float*> refdata;
		for (size_t i = 0; i < 13; i++)
		{
			refdata.push_back(ds->readData(i + 1));
		}
		float* band1 = refdata[0];
		float* band2 = refdata[1];
		float* band3 = refdata[2];
		float* band4 = refdata[3];
		float* band5 = refdata[4];
		float* band6 = refdata[5];
		float* band7 = refdata[6];
		float* band8 = refdata[7];
		float* band9 = refdata[8];
		float* band10 = refdata[9];
		float* band11 = refdata[10];
		float* band12 = refdata[11];
		float* band13 = refdata[12];
		float* nirband = band5;
		float* redband = band4;
		float* greenband = band3;
		float* blueband = band2;
		float* fmaskband = band13;
		float* lstband = band10;

		GDAL_DSInfo dsinfo = *((GDAL_DSInfo*)ds);

		dsinfo.numbands = 1;
		char* cmaskband = new char[dsinfo.slice];
		int numcloudpixs = 0;
		for (size_t npixel = 0; npixel < dsinfo.slice; npixel++)
		{
			cmaskband[npixel] = (char)(fmaskband[npixel]);
			if (cmaskband[npixel] == 4)
				numcloudpixs++;
		}
		if ((float)numcloudpixs / (float)dsinfo.slice > 0.1)
		{
			for (size_t i = 0; i < 13; i++)
			{
				delete[] refdata[i];
			}
			delete[] cmaskband;
			delete ds;
			continue;
		}
		GDAL_DS<char>* dscloud = new GDAL_DS<char>();
		dscloud->setDSInfo(&dsinfo);
		dscloud->create(outdir + name + "_cloud.tif");
		dscloud->writeData(1, cmaskband, 255);
		delete[] cmaskband;

		float* fndviband = new float[dsinfo.slice];
		GDAL_DS<float>* dsNDVI = new GDAL_DS<float>();
		dsNDVI->setDSInfo(&dsinfo);
		dsNDVI->create(outdir + name + "_NDVI.tif");

		float* feviband = new float[dsinfo.slice];
		GDAL_DS<float>* dsEVI = new GDAL_DS<float>();
		dsEVI->setDSInfo(&dsinfo);
		dsEVI->create(outdir + name + "_EVI.tif");

		float* fbrightnessband = new float[dsinfo.slice];
		GDAL_DS<float>* dsBrightness = new GDAL_DS<float>();
		dsBrightness->setDSInfo(&dsinfo);
		dsBrightness->create(outdir + name + "_Brightness.tif");

		float* fgreennessband = new float[dsinfo.slice];
		GDAL_DS<float>* dsGreennessband = new GDAL_DS<float>();
		dsGreennessband->setDSInfo(&dsinfo);
		dsGreennessband->create(outdir + name + "_Greenness.tif");

		float* fwetnessband = new float[dsinfo.slice];
		GDAL_DS<float>* dsWetnessband = new GDAL_DS<float>();
		dsWetnessband->setDSInfo(&dsinfo);
		dsWetnessband->create(outdir + name + "_Wetness.tif");


		GDAL_DS<float>* dsThermalbandband = new GDAL_DS<float>();
		dsThermalbandband->setDSInfo(&dsinfo);
		dsThermalbandband->create(outdir + name + "_lst.tif");

		// tasseled cap index Muhammad Hasan Ali Baig et al(2014)
		for (size_t npixel = 0; npixel < dsinfo.slice; npixel++)
		{
			float ndvi = (nirband[npixel] - redband[npixel]) / (nirband[npixel] + redband[npixel]);
			if (ndvi < -1 || ndvi > 1 || ndvi != ndvi || isinf(ndvi) || isnan(ndvi))
				ndvi = 0;
			float evi = 2.5 * (nirband[npixel] - redband[npixel]) / (nirband[npixel] + 6 * redband[npixel] - 7.5 * blueband[npixel] + 1);

			//'2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
			//	'NIR': image.select('B5'),
			//	'RED' : image.select('B4'),
			//	'BLUE' : image.select('B2')
			float brightness = band2[npixel] * 0.3029 + band3[npixel] * 0.2786 + band4[npixel] * 0.4733 + band5[npixel] * 0.5599 + band6[npixel] * 0.508 + band7[npixel] * 0.1872;
			float greenness = band2[npixel] * -0.2941 + band3[npixel] * -0.243 + band4[npixel] * -0.5424 + band5[npixel] * 0.7276 + band6[npixel] * 0.0713 + band7[npixel] * -0.1608;
			float wetness = band2[npixel] * 0.1511 + band3[npixel] * 0.1973 + band4[npixel] * 0.3283 + band5[npixel] * 0.3407 + band6[npixel] * -0.7117 + band7[npixel] * -0.4559;
			fndviband[npixel] = ndvi;
			feviband[npixel] = evi;

			fbrightnessband[npixel] = brightness;
			fgreennessband[npixel] = greenness;
			fwetnessband[npixel] = wetness;
			lstband[npixel] = lstband[npixel] - 273.15;
		}


		//-273.15
		dsNDVI->writeData(1, fndviband, 0);
		dsEVI->writeData(1, feviband, 0);
		dsBrightness->writeData(1, fbrightnessband, 0);
		dsGreennessband->writeData(1, fgreennessband, 0);
		dsWetnessband->writeData(1, fwetnessband, 0);
		dsThermalbandband->writeData(1, lstband, 0);
		delete[] fndviband;
		delete[] feviband;
		delete[] fbrightnessband;
		delete[] fgreennessband;
		delete[] fwetnessband;
		delete dsNDVI;
		delete dsEVI;
		delete dsBrightness;
		delete dsGreennessband;
		delete dsWetnessband;
		delete dsThermalbandband;

		//// Compute the EVI using an expression.
		//var evi = image.expression(
		//	'2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
		//		'NIR': image.select('B5'),
		//		'RED' : image.select('B4'),
		//	'BLUE' : image.select('B2')
		//	});


		//dsinfo.numbands = 3;

		//std::vector<float*> refdata;
		//refdata.push_back(ds->readData(4));
		//refdata.push_back(ds->readData(3));
		//refdata.push_back(ds->readData(2));

		//GDAL_DS<char>* dsrgb = new GDAL_DS<char>();
		//dsrgb->setDSInfo(&dsinfo);
		//dsrgb->create(outdir + name + "_rgb.tif");
		//for (size_t nband = 0; nband < 3; nband++)
		//{
		//	char* rgbband = new char[dsinfo.slice];
		//	float* refband = refdata[nband];
		//	float fmin = 1;
		//	float fmax = 0;
		//	for (size_t npixel = 0; npixel < dsinfo.slice; npixel++)
		//	{
		//		if (fmin > refband[npixel])
		//			fmin = refband[npixel];
		//		if (fmax < refband[npixel])
		//			fmax = refband[npixel];
		//	}

		//	for (size_t npixel = 0; npixel < dsinfo.slice; npixel++)
		//	{
		//		rgbband[npixel] = (char)((refband[npixel] - fmin) / (fmax-fmin) * 255);
		//	}

		//	dsrgb->writeData(nband + 1, rgbband, -1);
		//	delete[] rgbband;
		//	delete[] refband;
		//}
		//delete dsrgb;
		for (size_t i = 0; i < 13; i++)
		{
			delete[] refdata[i];
		}
		delete ds;

	}

	//std::vector<std::string> files = Utils::findFiles("E:/LanduseRegressionBaltimore/Landsat8/");
}
