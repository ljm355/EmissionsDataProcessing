#include "ScaleByFuel.h"
#include <fstream>
#include <sstream>
#include "Utils.h"
#include <map>
#include "qfileinfo.h"
struct CAInstance
{
	std::string sector;
	std::string fuel;
	std::string weightField;
	std::string shapeFileName;
	std::vector<double> totals;
};

ScaleByFuel::ScaleByFuel()
{
}


ScaleByFuel::~ScaleByFuel()
{
}


void ScaleByFuel::spatialAllocationTotal(ShapeFile* shp, std::vector<std::string> fuelfieldnames, std::string outputfieldname)
{

	std::vector<int> fuelfields;
	for (size_t i = 0; i < fuelfieldnames.size(); i++)
	{
		int fuelfieldidx = shp->poLayer->GetLayerDefn()->GetFieldIndex(fuelfieldnames[i].data());
		fuelfields.push_back(fuelfieldidx);
	}

	int outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	if (outputfIdx < 0)
	{
		OGRFieldDefn def(outputfieldname.data(), OGRFieldType::OFTReal);
		OGRErr er = shp->poLayer->CreateField(&def);
		outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	}
	//else
	//{
	//	return;
	//}
	OGRFeature *poFeature;

	shp->poLayer->ResetReading();
	while ((poFeature = shp->poLayer->GetNextFeature()) != NULL)
	{
		double sum = 0;
		for (size_t i = 0; i < fuelfieldnames.size(); i++)
		{
			if (fuelfields[i] > -1)
				sum += poFeature->GetFieldAsDouble(fuelfields[i]);
		}
		poFeature->SetField(outputfIdx, sum);
		shp->poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}



}
void ScaleByFuel::spatialAllocationForNonPoint(ShapeFile* shp, int baseyear, std::string weightField, std::string fuel, double total, std::string outputfieldname)
{

	std::vector<double> weights;
	int weightIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(weightField.data());
	if (weightIdx == -1)
	{
		printf("%s,%s\n", "weight field not found", weightField.data());
		return;
	}

	int yearIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex("yearbuilt");
	//int neeguiIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(("neeui_" + fuel).data());
	//int tfsIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex("tfs");

	int outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	if (outputfIdx < 0)
	{
		OGRFieldDefn def(outputfieldname.data(), OGRFieldType::OFTReal);
		OGRErr er = shp->poLayer->CreateField(&def);
		outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	}
	//else
	//{
	//	return;
	//}


	OGRFeature *poFeature;
	double sum = 0;
	shp->poLayer->ResetReading();
	int idx = -1;
	while ((poFeature = shp->poLayer->GetNextFeature()) != NULL)
	{
		idx++;

		if (poFeature->GetFieldAsInteger(yearIdx) > baseyear)
		{
			weights.push_back(0);
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		double weight = poFeature->GetFieldAsDouble(weightIdx);
		//if (neeguiIdx > -1)
		//weight *= poFeature->GetFieldAsDouble(neeguiIdx);
		weights.push_back(weight);
		sum += weight;
		OGRFeature::DestroyFeature(poFeature);
	}

	shp->poLayer->ResetReading();
	idx = -1;
	while ((poFeature = shp->poLayer->GetNextFeature()) != NULL)
	{
		idx++;

		if (poFeature->GetFieldAsInteger(yearIdx) > baseyear)
		{
			poFeature->SetField(outputfIdx, 0.0);
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}
		double val = weights[idx] / sum * total;
		if (sum == 0 || total == 0)
			val = 0;
		poFeature->SetField(outputfIdx, val);
		shp->poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

}
void ScaleByFuel::spatialAllocation(ShapeFile* shp, std::string weightField, double total, std::string outputfieldname)
{

	std::vector<double> weights;

	int weightIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(weightField.data());
	if (weightIdx == -1)
	{
		printf("%s,%s\n", "weight field not found", weightField.data());
		return;
	}
	int outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	if (outputfIdx < 0)
	{
		OGRFieldDefn def(outputfieldname.data(), OGRFieldType::OFTReal);
		OGRErr er = shp->poLayer->CreateField(&def);
		outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	}

	//const char* ignorefields[] = { "OGR_GEOMETRY" };
	////shp->poLayer->SetIgnoredFields(&ignorefields[0]);
	//shp->poLayer->SetIgnoredFields(ignorefields);
	OGRFeature *poFeature;
	double sum = 0;
	shp->poLayer->ResetReading();
	while ((poFeature = shp->poLayer->GetNextFeature()) != NULL)
	{
		weights.push_back(poFeature->GetFieldAsDouble(weightIdx));
		sum += weights[weights.size() - 1];
		OGRFeature::DestroyFeature(poFeature);
	}


	shp->poLayer->ResetReading();
	int idx = -1;
	while ((poFeature = shp->poLayer->GetNextFeature()) != NULL)
	{
		idx++;
		double val = weights[idx] / sum * total;
		if (sum == 0 || total == 0)
			val = 0;
		poFeature->SetField(outputfIdx, val);
		shp->poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
}
void ScaleByFuel::spatialAllocationElecProd(ShapeFile* shp, std::string weightField, double total, std::string outputfieldname)
{
	if (weightField == outputfieldname)
		return;
	std::vector<double> weights;

	int weightIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(weightField.data());
	if (weightIdx == -1)
	{
		printf("%s,%s\n", "weight field not found", weightField.data());
		return;
	}
	int outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	int sourceIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex("source");
	if (outputfIdx < 0)
	{
		OGRFieldDefn def(outputfieldname.data(), OGRFieldType::OFTReal);
		OGRErr er = shp->poLayer->CreateField(&def);
		outputfIdx = shp->poLayer->GetLayerDefn()->GetFieldIndex(outputfieldname.data());
	}
	OGRFeature *poFeature;
	double sum = 0;
	double sumNonNEI = 0;
	shp->poLayer->ResetReading();
	while ((poFeature = shp->poLayer->GetNextFeature()) != NULL)
	{
		std::string source = poFeature->GetFieldAsString(sourceIdx);
		if (source[0] != 'N' && source[0] != 'n')
		{
			sumNonNEI += poFeature->GetFieldAsDouble(outputfIdx);
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}

		weights.push_back(poFeature->GetFieldAsDouble(weightIdx));
		sum += weights[weights.size() - 1];
		OGRFeature::DestroyFeature(poFeature);
	}


	shp->poLayer->ResetReading();
	int idx = -1;
	total = total - sumNonNEI;
	while ((poFeature = shp->poLayer->GetNextFeature()) != NULL)
	{
		std::string source = poFeature->GetFieldAsString(sourceIdx);
		if (source[0] != 'N' && source[0] != 'n')
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}

		idx++;
		double val = weights[idx] / sum * total;
		if (sum == 0 || total == 0)
			val = 0;
		poFeature->SetField(outputfIdx, val);
		shp->poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
}

void ScaleByFuel::scale(std::string configfile, std::string dir)
{
	

	std::ifstream ifs;
	ifs.open(configfile);
	std::string line;
	std::string weightField = "ca";
	int baseyear = 2011;
	std::getline(ifs, line);
	std::vector<std::string> splits = Utils::splitCSV(',', line);
	std::vector<std::string> years;
	std::vector<int> yearsidx;
	int sectoridx = -1;
	int fuelidx = -1;
	int weightidx = -1;
	int shapefileidx = -1;
	for (size_t i = 0; i < splits.size(); i++)
	{
		std::string str = splits[i];
		if (str == "")
			continue;
		if (str == "weight"){
			weightidx = i;continue;
		}
		if (str == "fuel") {
			fuelidx = i; continue;
		}
		if (str == "sector") {
			sectoridx = i; continue;
		}
		if (str == "shapefile") {
			shapefileidx = i; continue;
		}
		bool isyear = true;
		for (size_t j = 0; j < str.size(); j++)
		{
			if(isalpha(str[j]))
			{
				isyear = false;
				break;
			}
		}
		if (isyear)
		{
		   years.push_back(splits[i]);
		   yearsidx.push_back(i);
		}

	}
	std::vector<CAInstance> initialGroup;
	int linenum = 0;
	while (ifs.peek() != -1)
	{
		if (linenum == 27)
		{
			printf("%d\n", linenum);
		}

		linenum++;
		std::getline(ifs, line);
		splits = Utils::splitCSV(',', line);
		if (splits.size() < 9)
		{
			splits = Utils::splitCSV(',', line);
			printf("%s\n", line.data());
		}
		CAInstance instance;
		if(sectoridx > -1)
		   instance.sector = splits[sectoridx];


		instance.fuel = splits[fuelidx];
		instance.shapeFileName = splits[shapefileidx];
		if (instance.shapeFileName == "")
			continue;
		std::string shapefilename = dir + instance.shapeFileName + ".shp";
		if (!QFileInfo(shapefilename.data()).exists())
		{
			continue;
		}
		instance.weightField = splits[weightidx];
		for (size_t i = 0; i < yearsidx.size(); i++)
		{
			instance.totals.push_back(atof(splits[yearsidx[i]].data()) * 1000000 * 1000);
		}
		initialGroup.push_back(instance);
	}


	std::map<std::string, std::vector<CAInstance>> groups;

	for (size_t i = 0; i < initialGroup.size(); i++)
	{
		CAInstance& instance = initialGroup[i];
		if (groups.find(instance.shapeFileName) == groups.end())
		{
			std::vector<CAInstance> newGroup;
			groups[instance.shapeFileName] = newGroup;
		}
		groups[instance.shapeFileName].push_back(instance);
	}


	std::map<std::string, std::vector<CAInstance>>::iterator iter = groups.begin();
	while (iter != groups.end())
	{
		std::vector<CAInstance>& shapefileGroup = iter->second;
		std::string shapefilename = dir + shapefileGroup[0].shapeFileName + ".shp";
		if (!QFileInfo(shapefilename.data()).exists())
		{
			iter++; continue;
		}
		if (!QString(shapefileGroup[0].shapeFileName.data()).endsWith("NonPoint"))
		{
			iter++; continue;
		}
		ShapeFile shp(shapefilename, 1);
		//bool canignore = shp.poLayer->TestCapability("OLCIgnoreFields ");
		//const char* ignorefields[] = { "ca","OGR_GEOMETRY" };
		//shp->poLayer->SetIgnoredFields(&ignorefields[0]);
		//shp.poLayer->SetIgnoredFields(ignorefields);
		for (size_t iyear = 0; iyear < years.size(); iyear++)
		{
			printf("%s,%s\n", shapefileGroup[0].shapeFileName.data(), years[iyear].data());
			std::vector<std::string> fuelfields;
			for (size_t i = 0; i < shapefileGroup.size(); i++)
			{
				CAInstance& instance = shapefileGroup[i];
				std::string fuel = instance.fuel;
				std::string outfieldname = "ca" + years[iyear].substr(2, 2);
				if (fuel != "")
				{
					outfieldname = outfieldname  + "_" + fuel;
				}

				double scaledTotal = instance.totals[iyear];
				//printf("%d,%d\n", iyear, i);
				weightField = instance.weightField;
				if (QString(shapefileGroup[0].shapeFileName.data()).endsWith("NonPoint"))
				{
					spatialAllocationForNonPoint(&shp, atoi(years[iyear].data()), weightField, fuel, scaledTotal, outfieldname);
				}
				else if (QString(shapefileGroup[0].shapeFileName.data()).endsWith("ElecProd"))
				{
					spatialAllocationElecProd(&shp, weightField, scaledTotal, outfieldname);
				}
				else
				{
					spatialAllocation(&shp, weightField, scaledTotal, outfieldname);

				}
				fuelfields.push_back(outfieldname);
			}
			if(fuelfields.size() > 1)
			   spatialAllocationTotal(&shp, fuelfields, "ca" + years[iyear].substr(2, 2));
		}
		iter++;
	}
}
#include <qdir.h>
void ScaleByFuel::scaleNonRoad(std::string rootdir, std::string cfgfile,std::string weightfield)
{
	std::ifstream ifs;
	ifs.open(cfgfile);
	std::string line;
	
	std::getline(ifs, line);
	std::vector<std::string> splits = Utils::splitCSV(',', line);
	std::vector<std::string> cafields;
	for (size_t i = 1; i < splits.size(); i++)
	{
		cafields.push_back(splits[i]);
	}

	rootdir = (QDir(rootdir.data()).absolutePath() + "/").toLocal8Bit().data();
	while (ifs.peek() != -1)
	{
		
		std::getline(ifs, line);
		splits = Utils::splitCSV(',', line);

		std::string subdirname = splits[0];
		std::vector<double> totals;
		for (size_t i = 1; i < splits.size(); i++)
		{
			totals.push_back(atof(splits[i].data()) * 1000000 * 1000);
		}

		std::vector<std::string> allfilesUnderDir = Utils::findFiles(rootdir + subdirname,".shp");
		std::vector<std::string> nonroadfiles;
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
			std::string dbffile = rootdir + subdirname + "/" + name.toLocal8Bit().data() +".dbf";
			std::string commandline = "RScript B:/LA_Version2/rename_field.R " + dbffile + " ca11 ca";
			//system(commandline.data());
		}
		printf("%s\n", subdirname.data());
		scaleNonRoad(nonroadfiles, cafields, totals, weightfield);

	}
	ifs.close();
}

void ScaleByFuel::scaleNonRoad(std::vector<std::string> shapefiles, std::vector<std::string>  fields, std::vector<double> totals,std::string weightField)
{

	std::vector<std::vector<double>> valuesArr;
	std::vector<std::vector<double>> subfractionsArr;
	std::vector<double> fractions;
	std::vector<double> sums;
	double total = 0;
	for (size_t i = 0; i < shapefiles.size(); i++)
	{
		std::vector<double> values;
		double sum = sumField(shapefiles[i], weightField, values);
		fractions.push_back(sum);
		sums.push_back(sum);
		valuesArr.push_back(values);
		total += sum;
		std::vector<double> subfractions;
		for (size_t j = 0; j < values.size(); j++)
		{
			subfractions.push_back(values[j] / sum);
		}
		subfractionsArr.push_back(subfractions);
	}
	for (size_t i = 0; i < shapefiles.size(); i++)
	{
		fractions[i] = fractions[i] / total;
	}
	for (size_t i = 0; i < shapefiles.size(); i++)
	{
		ShapeFile shp(shapefiles[i], 1);
		printf("%s\n", shapefiles[i].data());

		//int tsfield = shp.getOrCreateField("timestruct", OGRFieldType::OFTString);

		std::vector<int> fieldindices;
		for (size_t ifield = 0; ifield < fields.size(); ifield++)
		{
			std::string cafield = fields[ifield];
			int calFieldIdx = shp.getOrCreateField(cafield.data(), OGRFieldType::OFTReal);
			fieldindices.push_back(calFieldIdx);
		}
		
		OGRFeature *poFeature;
		double sum = 0;
		shp.poLayer->ResetReading();
		std::vector<double>& subfractions = subfractionsArr[i];

		double fraction = fractions[i];
		int idx = 0;
		while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
		{
			for (size_t ifield = 0; ifield < fields.size(); ifield++)
			{
				double val = totals[ifield] * fraction * subfractions[idx];
				poFeature->SetField(fieldindices[ifield], val);
			}
			shp.poLayer->SetFeature(poFeature);
			idx++;
			OGRFeature::DestroyFeature(poFeature);
		}
	}
}

double ScaleByFuel::sumField(std::string shapefile, std::string fieldname, std::vector<double>& values)
{
	ShapeFile shp(shapefile);
	int weightIdx = shp.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data());
	if (weightIdx == -1)
	{
		printf("weight field not found:%s\n", fieldname.data());
		return 0;
	}

	OGRFeature *poFeature;
	double sum = 0;
	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		double val = poFeature->GetFieldAsDouble(weightIdx);
		sum += val;
		values.push_back(val);
		OGRFeature::DestroyFeature(poFeature);
	}
	return sum;

}
