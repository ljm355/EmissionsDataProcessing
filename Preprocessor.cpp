#include "Preprocessor.h"
#include "qfileinfo.h"
#include <fstream>
#include <sstream>
#include <qdir.h>
#include <qstring>
//#include "archive_entry.h"
//#include "archive.h"

Preprocessor::Preprocessor()
{

}
Preprocessor::~Preprocessor()
{

}
void Preprocessor::updateFieldAfterIntersection(std::string filename,double scale)
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
	int fracField = input.getOrCreateField("fraction", OGRFieldType::OFTReal);
	while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
	{
		OGRGeometry* geo = poFeature->GetGeometryRef();
		OGRwkbGeometryType gtype = geo->getGeometryType();

		double footprintOld = poFeature->GetFieldAsDouble(footprintIndex);
		double footprintNew = 0;

		if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
		{
			footprintNew = Utils::calPolylineLength(poFeature->GetGeometryRef()) * scale;
		}
		else if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype == wkbPolygon25D)
		{
			footprintNew = Utils::calPolygonArea(poFeature->GetGeometryRef()) * scale;
		}

		double fraction = 0;
		if (!isinf(footprintOld) && footprintOld > 0)
			fraction = footprintNew / footprintOld;
		poFeature->SetField(fracField, fraction);
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
void Preprocessor::updateFieldAfterIntersection(std::string filename, std::vector<std::string> fieldnames, double scale)
{
	ShapeFile input(filename, 1);
	OGRFeature *poFeature;
	input.poLayer->ResetReading();
	std::vector<int> fields;
	for (size_t i = 0; i < fieldnames.size(); i++)
	{
		int fieldidx = input.poLayer->GetLayerDefn()->GetFieldIndex(fieldnames[i].data());
		if(fieldidx > -1)
		   fields.push_back(fieldidx);
	}

	int units_idx = input.getOrCreateField("units",OGRFieldType::OFTReal);
	int bt_idx = input.poLayer->GetLayerDefn()->GetFieldIndex("bt");

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
			footprintNew = Utils::calPolylineLength(poFeature->GetGeometryRef()) * scale;
		}
		else if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype == wkbPolygon25D)
		{
			footprintNew = Utils::calPolygonArea(poFeature->GetGeometryRef()) * scale;
		}

		int bt = poFeature->GetFieldAsInteger(bt_idx);
		float units = 1;
		if (bt == 4)
			units = 3;
		else if (bt == 5)
			units = 5;

		double fraction = 0;
		if (!isinf(footprintOld) && footprintOld > 0)
			fraction = footprintNew / footprintOld;
		if (fraction >  0.999999)
			fraction = 1;
		else if (fraction <  0.000001)
		{
			fraction = 0;
		}
		units = fraction * units;
		poFeature->SetField(units_idx, units);
		poFeature->SetField(footprintIndex, footprintNew);
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
void Preprocessor::updateFieldAfterIntersectionForDir(std::string indir,double scale)
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
		updateFieldAfterIntersection(infile,scale);
	}
}

void Preprocessor::gridFolderByShape(std::string indir, std::string outdir, std::string fishneshapefile, double scale)
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

	for (size_t i = 0; i < files.size(); i++)
	{
		printf("%s\n", (outdir + files[i]).data());

		if (QFileInfo((outdir + files[i]).data()).exists())
		{
			continue;
			//ShapeFile::remove(outdir + files[i]);
		}
		//if (!QFileInfo((outdir + files[i]).data()).exists())
		//{
			Preprocessor::intersectWithArcGIS(indir + files[i], fishneshapefile, outdir + files[i]);
			updateFieldAfterIntersection(outdir + files[i], scale);
		//}
	}

}

void Preprocessor::gridFolderByRaster(std::string indir, std::string outdir, std::string fishnetrasterfile, double scale)
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

	std::string outfishnetfile = outdir + "fishnet.shp";
	std::string outrasterfile = outdir + "fishnet.tif";

	std::string fishnetfile = (QFileInfo(fishnetrasterfile.data()).absoluteDir().absolutePath() + "/fishnet.shp").toLocal8Bit().data();
	std::string rasterfile = fishnetrasterfile;

	//if (!boundfile.poLayer->GetSpatialRef()->IsGeographic())
	//{
	//	//FOOTPRINT_SCALE_FACTOR = 1;
	//}
	//else
	//{
	//	//FOOTPRINT_SCALE_FACTOR = getDegreeToMeter(boundfile.poLayer);
	//	gridsize = gridsize / Utils::getDegreeToMeter(boundfile.poLayer);
	//}


	Grid fishnet;
	fishnet.fromFishnetRaster(fishnetrasterfile);
	fishnet.reset();

	GDALDataset* pDataset = (GDALDataset*)GDALOpen(fishnetrasterfile.data(), GA_ReadOnly);
	std::string wkt = pDataset->GetProjectionRef();
	GDALClose(pDataset);
	bool needCreateFishnet = false;
	//if (!QFileInfo(fishnetfile.data()).exists())
	//{
		needCreateFishnet = true;
		fishnet.toShape(wkt, fishnetfile, false);
	//}
	for (size_t i = 0; i < files.size(); i++)
	{
		/*if (skipNonRoad && files[i] == "NonRoad.shp")
			continue;*/
		printf("%s\n", (outdir + files[i]).data());
		//intersectWithFishnet(fishnet, fishnetfile, indir + files[i], outdir + files[i]);
		
		if (!QFileInfo((outdir + files[i]).data()).exists())
		{
			Utils::updateFootprint(indir + files[i], scale);
			Preprocessor::intersectWithArcGIS(indir + files[i], fishnetfile, outdir + files[i]);
			while (!QFileInfo((outdir + files[i]).data()).exists())
			{
				printf("Ê§°Ü!");
				Preprocessor::intersectWithArcGIS(indir + files[i], fishnetfile, outdir + files[i]);
			}
			updateFieldAfterIntersection(outdir + files[i], scale);
		}
		if (needCreateFishnet) {
			ShapeFile input((outdir + files[i]).data());
			fishnet.gatherCells(&input);
		}

	}

	if (needCreateFishnet) {
		fishnet.toShape(wkt, outfishnetfile, true);
		fishnet.toRaster(outrasterfile, wkt);
	}
	//for (size_t i = 0; i < files.size(); i++)
	//{
	//	if (skipNonRoad && files[i] == "NonRoad.shp")
	//		continue;
	//	ShapeFile input((outdir + files[i]).data());
	//	if (input.poLayer->GetGeomType() == wkbPoint || input.poLayer->GetGeomType() == wkbMultiPoint || input.poLayer->GetGeomType() == wkbPoint25D)
	//		continue;
	//	std::string sectorfile = files[i].substr(0, files[i].length() - 4);
	//	std::string sectorfileout = outdir + sectorfile + ".tif";
	//	printf("%s\n", sectorfile.data());
	//	fishnet->reset();
	//	gatherCells(fishnet, &input);
	//	fishnet->toRaster(sectorfileout);

	//	sectorfileout = outdir + sectorfile + "2.shp";
	//	fishnet->toShape(&input, sectorfileout,true);

	//}

	//delete fishnet;
}
void Preprocessor::reprojectWithArcGIS(std::string inputfile, std::string outputfile, std::string pythonScriptFile)
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
void Preprocessor::reprojectDir(std::string indir, std::string outdir, std::string pythonFile)
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
		std::string input_file = (input_dir.absolutePath() + "/" + fileInfo.fileName()).toLocal8Bit().data();
		if (!fileInfo.fileName().endsWith(".shp"))
			continue;
		std::string output_file = outdir + fileInfo.fileName().toLocal8Bit().data();
		if (QFileInfo(output_file.data()).exists())
			continue;
			//ShapeFile::remove(output_file);
		reprojectWithArcGIS(input_file, output_file, pythonFile);
		printf("%s\n", output_file.data());
	}

}
void Preprocessor::intersectWithArcGIS(std::string inputfile, std::string fishnet, std::string outputfile)
{
	QFileInfo info(outputfile.data());

	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName);
	poDriver->Delete(outputfile.data());
	std::ofstream ofs;
	std::string scriptFile = (info.absoluteDir().absolutePath() + "/" + info.completeBaseName() + ".py").toLocal8Bit().data();
	ofs.open(scriptFile.data());
	std::stringstream script;
	script << "import arcpy" << "\n";
	////script << "fishnetfile = " << "\"" << fishnetfile << "\"" << "\n";
	////script << "inputfile = " << "\"" << inputfile << "\"" << "\n";
	//script << "outputfile = " << "\"" << outputfile << "\"" << "\n";
	script << "arcpy.Intersect_analysis(";
	//script << "fishnetfile" << ";" << "inputfile" << "," << "outputfile" << "," << "\"ALL\"" << "," << "\"\"" << ","  "\"INPUT\")";
	//script << "\"" << fishnetfile << ";" << inputfile << "\"" << "," << "outputfile" << "," << "\"NO_FID\"" << "," << "\"\"" << ","  "\"INPUT\")";
	script << "\"" << inputfile << ";" << fishnet << "\"" << "," << "\"" << outputfile << "\"" << "," << "\"ALL\"" << "," << "\"\"" << ","  "\"INPUT\")";
	ofs << script.str();
	ofs.close();
	//printf("%s\n", script.str().data());
	system(scriptFile.data());

	//QFile::remove(scriptFile.data());
}
void Preprocessor::clipWithArcGIS(std::string fishnetfile, std::string inputfile, std::string outputfile)
{
	QFileInfo info(outputfile.data());

	const char *pszDriverName = "ESRI Shapefile";
	GDALDriver *poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName);
	poDriver->Delete(outputfile.data());
	std::ofstream ofs;
	std::string scriptFile = (info.absoluteDir().absolutePath() + "/" + info.completeBaseName() + ".py").toLocal8Bit().data();
	ofs.open(scriptFile.data());
	std::stringstream script;
	script << "import arcpy" << "\n";
	script << "arcpy.Intersect_analysis(";
	script << "\"" << fishnetfile << ";" << inputfile << "\"" << "," << "\"" << outputfile << "\"" << "," << "\"\")";
	ofs << script.str();
	ofs.close();
	system(scriptFile.data());
	QFile::remove(scriptFile.data());
}
void Preprocessor::intersectWithFishnet(Grid* grid, const std::string& fishnetfile, const std::string& inputfile, const std::string& outputfile)
{
	bool hasFile = false;
	if (QFileInfo(outputfile.data()).exists())
	{
		hasFile = true;
	}
	else
	{
		intersectWithArcGIS(fishnetfile, inputfile, outputfile);
	}

}
std::string Preprocessor::getFootprintFieldName(OGRwkbGeometryType gtype)
{
	if (gtype == wkbLineString || gtype == wkbLineString25D || gtype == wkbMultiLineString)
	{
		return "length";
	}
	else if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype == wkbPolygon25D)
	{
		return "area";
	}
	return "";
}
bool Preprocessor::has25D(OGRGeometry* geom)
{
	OGRwkbGeometryType gtype = geom->getGeometryType();
	if (gtype == wkbLineString25D || gtype == wkbMultiLineString25D)
	{
		return true;
	}
	else if (gtype == wkbMultiLineString)
	{
		OGRMultiLineString* multipoly = (OGRMultiLineString*)geom;
		for (size_t i = 0; i < multipoly->getNumGeometries(); i++)
		{

			if (multipoly->getGeometryRef(i)->getGeometryType() == geom->getGeometryType() == wkbLineString25D)
				return true;
		}
	}
}
OGRLineString* Preprocessor::convertOGRLineStringFrom25D(OGRGeometry* geom)
{
	OGRLineString* oldline = (OGRLineString*)geom;
	int len = oldline->getNumPoints();
	OGRLineString *newline = (OGRLineString*)OGRGeometryFactory::createGeometry(wkbLineString);
	for (size_t i = 0; i < len; i++)
	{
		OGRPoint p;
		oldline->getPoint(i, &p);
		newline->addPoint(p.getX(), p.getY());
	}
	return oldline;
}
OGRMultiLineString* Preprocessor::convertOGRMultiLineStringFrom25D(OGRGeometry* geom)
{
	OGRMultiLineString *oldmultipoly = (OGRMultiLineString*)geom;
	OGRMultiLineString *newmultipoly = (OGRMultiLineString*)OGRGeometryFactory::createGeometry(wkbMultiLineString);

	for (size_t i = 0; i < oldmultipoly->getNumGeometries(); i++)
	{
		newmultipoly->addGeometry(convertOGRLineStringFrom25D(oldmultipoly->getGeometryRef(i)));
	}
	return newmultipoly;
}
bool Preprocessor::convertPolylineGeometry(OGRGeometry* geom, OGRGeometry*& newgeom)
{
	//if (!has25D(geom))
	//	return false;
	double len = 0;
	OGRwkbGeometryType gtype = geom->getGeometryType();
	if (geom->getGeometryType() == wkbLineString25D || geom->getGeometryType() == wkbLineString)
	{
		newgeom = convertOGRLineStringFrom25D(geom);
	}
	else if (geom->getGeometryType() == wkbMultiLineString25D || geom->getGeometryType() == wkbMultiLineString)
	{
		newgeom = convertOGRMultiLineStringFrom25D(geom);
	}
	return true;
}
bool Preprocessor::convertPolyline(std::string infile, std::string outfile)
{
	ShapeFile oldShp(infile);
	OGRwkbGeometryType gtype = oldShp.poLayer->GetGeomType();
	if (gtype != wkbLineString && gtype != wkbLineString25D && gtype != wkbMultiLineString25D)
		return false;
	OGRwkbGeometryType newgtype = wkbLineString;
	if (gtype == wkbMultiLineString25D)
		newgtype = wkbMultiLineString;
	ShapeFile newShp;
	newShp.create(outfile.data(), oldShp.poLayer->GetSpatialRef(), oldShp.poLayer->GetLayerDefn(), newgtype);
	OGRFeature *oldFeature;
	oldShp.poLayer->ResetReading();
	int id = 0;
	while ((oldFeature = oldShp.poLayer->GetNextFeature()) != NULL)
	{

		OGRFeature* newFeature = OGRFeature::CreateFeature(newShp.poLayer->GetLayerDefn());
		for (size_t i = 0; i < oldShp.poLayer->GetLayerDefn()->GetFieldCount(); i++)
		{
			newFeature->SetField(i, oldFeature->GetRawFieldRef(i));
		}
		OGRGeometry* newGeom = NULL;
		if (convertPolylineGeometry(oldFeature->GetGeometryRef(), newGeom)) {
			newFeature->SetGeometry(newGeom);
		}
		else {
			newFeature->SetGeometry(oldFeature->GetGeometryRef());
		}
		newShp.poLayer->CreateFeature(newFeature);
		OGRFeature::DestroyFeature(newFeature);
		OGRFeature::DestroyFeature(oldFeature);
	}
	newShp.close();
	oldShp.close();
	return true;
}
std::vector<std::string> parseLine(char*& data, char& delimiter)
{
	char* pdata = data;
	char buf[1024];
	size_t num = 0;
	size_t start = 0;
	if (*pdata == delimiter)
	{
		num++;
		start++;
		pdata++;
	}
	std::vector<std::string> strs;
	while (true)
	{
		char c = *pdata;

		if (c == '\n' || c == delimiter)
		{
			//*pdata = '\0';
			//pdata = data + start;
			memcpy(buf, data + start, num);
			buf[num] = '\0';
			std::string str(buf);
			//*splits++ = str;
			strs.push_back(str);
			start += num + 1;
			num = 0;
			pdata++;
			if (c == '\n')
			{
				data = pdata;
				break;
			}
			continue;
		}

		pdata++;
		num++;
	}
	return strs;
}
bool parseFloats(char*& data, char& delimiter, float* values, int count)
{
	char* pdata = data;
	size_t num = 0;
	size_t start = 0;
	if (*pdata == delimiter)
	{
		num++;
		start++;
		pdata++;
	}
	std::vector<std::string> strs;
	while (true)
	{
		char c = *pdata;
		if (c == '\n' || c == delimiter)
		{

			*pdata = '\0';
			*values++ = atof(data + start);
			*pdata = delimiter;
			start += num + 1;
			pdata++;
			num = 0;
			count--;
			if (count == 0 || c == '\n')
			{
				data = pdata;
				break;
			}
			continue;
		}

		pdata++;
		num++;
	}
	return true;
}
bool parseInt(char*& data, char& delimiter, int* values, int count)
{
	char* pdata = data;
	size_t num = 0;
	size_t start = 0;
	if (*pdata == delimiter)
	{
		num++;
		start++;
		pdata++;
	}
	std::vector<std::string> strs;
	while (true)
	{
		char c = *pdata;
		if (c == '\n' || c == delimiter)
		{

			*pdata = '\0';
			*values++ = atoi(data + start);
			*pdata = delimiter;
			start += num + 1;
			pdata++;
			num = 0;
			count--;
			if (count == 0 || c == '\n')
			{
				data = pdata;
				break;
			}
			continue;
		}

		pdata++;
		num++;
	}
	return true;
}
void Preprocessor::gridShapeFile(std::string infile, std::string outfile, std::string fishnetfile,const char* fieldname)
{

	std::string tmpefile = (QFileInfo(outfile.data()).dir().absolutePath() + "/" + QFileInfo(outfile.data()).baseName() + "_2.shp").toLocal8Bit().data();
	std::string outfiletif = (QFileInfo(outfile.data()).dir().absolutePath() + "/" + QFileInfo(outfile.data()).baseName() + ".tif").toLocal8Bit().data();
	Utils::updateFootprint(infile);
	Preprocessor::intersectWithArcGIS(infile,fishnetfile, tmpefile);
	//updateFieldAfterIntersection(tmpefile);
	ShapeFile input(tmpefile);
	Grid fishnet;
	fishnet.fromFishnetShape(fishnetfile);
	fishnet.reset(fieldname);
	if (fieldname != "")
	{
		fishnet.gatherCells(&input, fieldname);
	}

	char wkt[512];
	char* pwkt = wkt;
	if (input.poLayer->GetSpatialRef())
		input.poLayer->GetSpatialRef()->exportToWkt(&pwkt);

	fishnet.toShape(wkt, outfile, true);
	fishnet.toRaster(outfiletif, wkt);
	const char *pszDriverName = "ESRI Shapefile";
	OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(pszDriverName)->Delete(tmpefile.data());
}
//void Preprocessor::readFromZipFile(const char* zipFile, std::vector<float*> dataArray, int xsize)
//{
//	int r;
//
//	struct archive_entry *entry;
//	struct archive *a = archive_read_new();
//	//archive_read_support_compression_all(a);
//	//archive_read_support_format_raw(a);
//	archive_read_support_format_all(a);
//	archive_read_support_filter_all(a);
//	r = archive_read_open_filename(a, zipFile, 8192);
//	r = archive_read_next_header(a, &entry);
//	//printf("%s\\n",archive_entry_pathname(entry));
//
//	size_t total = archive_entry_size(entry);
//	size_t size;
//	size_t bufsize = total;
//	char* buf = new char[bufsize];
//	size = archive_read_data(a, buf, bufsize);
//	std::string line;
//	size_t offset = 0;
//	char delimiter = ',';
//	char* pbuf = buf;
//	std::vector<std::string> splits = parseLine(pbuf, delimiter);
//	std::vector<float> values;
//	values.resize(24);
//	std::vector<int> index;
//	index.resize(2);
//
//	char* bufEnd = buf + bufsize;
//	int i = 0;
//	int numlines = 0;
//	int numfields = splits.size() - 2;
//	while (pbuf < bufEnd)
//	{
//		parseInt(pbuf, delimiter, &index[0], 2);
//		parseFloats(pbuf, delimiter, &values[0], 24);
//		int ncol = atoi(splits[0].data()) - 1;
//		int nrow = atoi(splits[1].data()) - 1;
//
//		for (int j = 0; j < numfields; j++)
//		{
//			//float val = values[j];
//			//if (val < 0)
//			//	val = 0;
//			//float val2 = dataArray[j][(index[1]-1) * xsize + (index[0]-1)];
//			//if (val != val2)
//			//printf("%f,%f\n",val,val2);
//			dataArray[j][(index[1] - 1) * xsize + (index[0] - 1)] = values[j];
//		}
//		numlines++;
//	}
//	delete[] buf;
//	archive_read_free(a);
//}

void Preprocessor::readHourlyFFDASNC(const char* filename, std::vector<float*> dataArrays)
{
	int status;
	int ncidSrc;
	status = nc_open(filename, NC_WRITE | NC_SHARE, &ncidSrc);

	int numdims = 0;
	nc_inq_ndims(ncidSrc, &numdims);
	std::vector<NetCDFDim> dims;
	for (int i = 0; i < numdims; i++)
	{
		NetCDFDim dim;
		int err = nc_inq_dim(ncidSrc, i, dim.name, &dim.len);
		err = nc_inq_dimid(ncidSrc, dim.name, &dim.id);
		dims.push_back(dim);
	}

	size_t start[] = { 0, 0 }; // start at first value 
    std::vector<size_t> latlondim; latlondim.push_back(dims[0].len); latlondim.push_back(dims[2].len);
	size_t* count = &latlondim[0];

	int numvars = 0;
	nc_inq_nvars(ncidSrc, &numvars);

	std::vector<NetCDFVar> vars;
	for (int i = 4; i < numvars; i++)
	{
		NetCDFVar var;
		nc_inq_var(ncidSrc, i, var.name, &var.vtype,
			&var.ndims, var.dimids, &var.varnatt);
		nc_inq_varid(ncidSrc, var.name, &var.id);

		for (size_t j = 0; j < var.varnatt; j++)
		{
			NetCDFAttribute attr;
			nc_inq_attname(ncidSrc, var.id, j, attr.name);
			nc_inq_att(ncidSrc, var.id, attr.name, &attr.vtype, &attr.len);
			nc_get_att_text(ncidSrc, var.id, attr.name, attr.text);
			var.attributes.push_back(attr);
		}
		vars.push_back(var);

	}

	//float* tmpData = new float[dims[0].len + dims[2].len];
	//for (size_t i = 0; i < 4; i++)
	//{
	//	NetCDFVar var = vars[i];
	//	std::vector<size_t> dimlens;
	//	for (size_t j = 0; j < var.ndims; j++)
	//	{
	//		dimlens.push_back(dims[var.dimids[j]].len);
	//	}
	//	nc_get_vara_float(ncidSrc, var.id, start, &dimlens[0], tmpData);
	//}
	//delete[] tmpData;
	numvars = vars.size();
	for (size_t i = 0; i < numvars; i++)
	{
		nc_get_vara_float(ncidSrc, i + 4, start, count, dataArrays[i]);
	}


	nc_close(ncidSrc);
	
}
void Preprocessor::extractFFDASHourly(std::string dir, std::string outfile, double* adftransform, int startrow, int nrows, int startcol, int ncols)
{
	std::vector<float*> dataArrays;
	
	int xsize = 3600;
	int ysize = 1800;
	int numcells = xsize*ysize;
	dataArrays.resize(24);
	for (int j = 0; j < dataArrays.size(); j++)
	{
		dataArrays[j] = new float[xsize*ysize];
	}
	float* buf = new float[nrows*ncols];
	QDir input_dir(dir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	//poDriver->Delete(outfile.data());
	GDALDataset* pDataset = poDriver->Create(outfile.data(), ncols, nrows, 8760, GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(adftransform);

	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		QString input_file = input_dir.absolutePath() + "/" + fileInfo.fileName();
		if (!input_file.endsWith(".nc"))
			continue;
		std::string name = fileInfo.baseName().toLocal8Bit().data();
		name = name.substr(1, 3);
		int day = atoi(name.data());
		int starthour = (day - 1) * 24;
		printf("day=%d,hour=%d,%s\n", day, starthour, name.data());
		readHourlyFFDASNC(input_file.toLocal8Bit().data(), dataArrays);
		for (int ihour = 0; ihour  < 24; ihour++)
		{
			float* pdestdata = buf;
			float* psrcdata = dataArrays[ihour];
			for (size_t irow = startrow; irow < startrow+nrows; irow++)
			{
				for (size_t icol = startcol; icol < startcol+ncols; icol++)
				{
					*pdestdata++ = psrcdata[irow*xsize+icol];
				}
			}
			pDataset->GetRasterBand(starthour + ihour + 1)->SetNoDataValue(0);
			pDataset->GetRasterBand(starthour + ihour + 1)->RasterIO(GF_Write, 0, 0, ncols, nrows, buf, ncols, nrows, GDT_Float32, 0, 0);
		}
	}
	GDALClose((GDALDatasetH)pDataset);
	for (int j = 0; j < dataArrays.size(); j++)
	{
		delete[] dataArrays[j];
	}
	delete[] buf;


}

void Preprocessor::extractEDGASHourly(std::string infile, std::string outfile, double* adftransform, int startrow, int nrows, int startcol, int ncols)
{
	int xsize = 3600;
	int ysize = 1800;
	int numcells = xsize*ysize;
	float* data = readData(infile.data());
	float* buf = new float[nrows*ncols];
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	//poDriver->Delete(outfile.data());
	GDALDataset* pDataset = poDriver->Create(outfile.data(), ncols, nrows, 1, GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(adftransform);
	float* pdestdata = buf;
	float* psrcdata = data;
	for (size_t irow = startrow; irow < startrow + nrows; irow++)
	{
		for (size_t icol = startcol; icol < startcol + ncols; icol++)
		{
			*pdestdata++ = psrcdata[irow*xsize + icol];
		}
	}
	pDataset->GetRasterBand(1)->SetNoDataValue(0);
	pDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, ncols, nrows, buf, ncols, nrows, GDT_Float32, 0, 0);

	GDALClose((GDALDatasetH)pDataset);
	delete[] data;
	delete[] buf;
}

void Preprocessor::extractGlobalRaster(std::string infile, std::string outfile, double* adftransform, int startrow, int nrows, int startcol, int ncols)
{
	int xsize;
	int ysize;
	std::vector<float*>  data = readData(infile.data(),xsize,ysize);
	int numcells = xsize*ysize;
	float* buf = new float[nrows*ncols];
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	//poDriver->Delete(outfile.data());
	GDALDataset* pDataset = poDriver->Create(outfile.data(), ncols, nrows, data.size(), GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(adftransform);

	for (size_t iband = 0; iband < data.size(); iband++)
	{
		float* pdestdata = buf;
		float* psrcdata = data[iband];
		for (size_t irow = startrow; irow < startrow + nrows; irow++)
		{
			for (size_t icol = startcol; icol < startcol + ncols; icol++)
			{
				*pdestdata++ = psrcdata[irow*xsize + icol];
			}
		}
		pDataset->GetRasterBand(iband+1)->SetNoDataValue(0);
		pDataset->GetRasterBand(iband + 1)->RasterIO(GF_Write, 0, 0, ncols, nrows, buf, ncols, nrows, GDT_Float32, 0, 0);
		delete[] data[iband];
	}
	
	GDALClose((GDALDatasetH)pDataset);

	delete[] buf;
}

std::vector<float*> Preprocessor::readData(const char * filename,int& xsize,int& ysize)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);

	std::vector<float*> results;
	xsize = poDataset->GetRasterXSize();
	ysize = poDataset->GetRasterYSize();

	int   ncells = xsize*ysize;
	for (size_t i = 0; i < poDataset->GetRasterCount(); i++)
	{
		GDALRasterBand *poBand = poDataset->GetRasterBand(i+1);
		float* srcData = new float[xsize*ysize];
		poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
			srcData, xsize, ysize, GDT_Float32,
			0, 0);
		results.push_back(srcData);
	}
	
	GDALClose((GDALDatasetH)poDataset);
	return results;
}

float* Preprocessor::readData(const char * filename)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);
	GDALRasterBand *poBand = poDataset->GetRasterBand(1);
	if (poDataset == NULL) {
		return NULL;
	}
	float *srcData;
	int   xsize = poBand->GetXSize();
	int   ysize = poBand->GetYSize();
	int   ncells = xsize*ysize;
	srcData = new float[xsize*ysize];
	poBand->RasterIO(GF_Read, 0, 0, xsize, ysize,
		srcData, xsize, ysize, GDT_Float32,
		0, 0);
	GDALClose((GDALDatasetH)poDataset);
	return srcData;
}
