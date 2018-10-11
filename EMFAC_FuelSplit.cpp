#include "EMFAC_FuelSplit.h"
#include <fstream>
#include <map>
#include <sstream>
#include "Utils.h"
EMFAC_FuelSplit::EMFAC_FuelSplit()
{

}


EMFAC_FuelSplit::~EMFAC_FuelSplit()
{

}

void EMFAC_FuelSplit::loadSCAG2HPMS(std::string crosswalkfile,std::string postfix,std::map<std::string,std::string>& crosswalk)
{
	std::string line;
	std::ifstream infile(crosswalkfile.data());
	int HPMSCode = 1;
	while (std::getline(infile, line))
	{
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		std::stringstream ss;
		ss << HPMSCode << postfix;
		for (size_t i = 1; i < splits.size(); i++)
		{
			if (splits[i] == "")
				continue;
			crosswalk[splits[i]] = ss.str();
		}
		HPMSCode++;
	}
	infile.close();
}


void EMFAC_FuelSplit::HPMS2SCAG(std::string onroadfile, std::string ruralCrosswalk, std::string urbanCrosswalk)
{
	std::map<std::string, std::string> crosswalk;
	loadSCAG2HPMS(ruralCrosswalk, "R", crosswalk);
	loadSCAG2HPMS(urbanCrosswalk, "U", crosswalk);

	ShapeFile shp(onroadfile,1);
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	int id = 0;
	int hpmsfield = shp.getOrCreateField("hpms", OGRFieldType::OFTString);
	int scagfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("rdTypeUR");
	int vmtfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("VMT");
	int cafield = shp.poLayer->GetLayerDefn()->GetFieldIndex("ca");
	std::map<int, OGRPoint> pointmaps;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		std::string scagType = fea->GetFieldAsString(scagfield);
		std::string hpmsType = crosswalk[scagType];
		fea->SetField(hpmsfield, hpmsType.data());
		if (hpmsType == "7U" || hpmsType == "7R")
		{
			fea->SetField(vmtfield, fea->GetFieldAsDouble(cafield));
		}
		shp.poLayer->SetFeature(fea);
		id++;
		OGRFeature::DestroyFeature(fea);
	}
	shp.close();
	printf("num of points:%d\n", id);
}

//Passenger cars - gasoline
//Passenger cars - diesel
//Passenger cars - electric
//Light Trucks - gasoline
//Light Trucks - diesel
//Light Trucks - electric
//Single - unit trucks - gasoline
//Single - unit trucks - diesel
//motorcycles - gasoline
//buses - gasoline
//buses - diesel
//combination trucks - gasoline
//combination trucks - diesel
//Total


#include <algorithm>
#include <qstring.h>
EMFAC_FuelSpecificEmissions EMFAC_FuelSplit::loadEMFAC(std::string EMFACTableFILE,std::string cityname)
{
	EMFAC_FuelSpecificEmissions table;
	//std::vector<EMFAC_FuelSpecificEmissions> table;
	std::string line;
	std::ifstream infile(EMFACTableFILE.data());
	std::vector<std::string> splits;
	//if (QString(EMFACTableFILE.data()).contains("2011"))
	//{
	//	printf("");
	//}
	while (true)
	{
		std::getline(infile, line);
		if (QString(line.data()).contains(QString(cityname.data()), Qt::CaseInsensitive))
			break;
	}

	//skip three lines.
	for (size_t i = 0; i < 3; i++)
	{
		std::getline(infile, line);
	}
	for (size_t ihpmsType = 0; ihpmsType < 7; ihpmsType++)
	{
		for (size_t ifuel = 0; ifuel < 3; ifuel++)
		{
			table.emissionsByHPMSRural[ihpmsType][ifuel] = 0;
			table.emissionsByHPMSUrban[ihpmsType][ifuel] = 0;
		}

	}
	printf("\n");
	for (size_t iline = 0; iline < 13; iline++)
	{

		std::getline(infile, line);
		//printf("%s\n", line.data());
		splits = Utils::splitCSV(',', line);
		//printf("%d,%s\n", iline, line.data());
		std::vector<std::string> splits2 = Utils::splitCSV('-', splits[0]);
		std::string fuelname = splits2[splits2.size() - 1];
		fuelname.erase(std::remove(fuelname.begin(), fuelname.end(), ' '), fuelname.end());
		int fuelcode = (int)getFuelCode(fuelname);
		int ruralstart = 1;
		int urbanstart = 9;
		if (fuelcode == FuelTypeCode::unknown)
			continue;
		for (size_t ihpmsType = 0; ihpmsType < 7; ihpmsType++)
		{
			table.emissionsByHPMSRural[ihpmsType][fuelcode] = table.emissionsByHPMSRural[ihpmsType][fuelcode] + atof(splits[ihpmsType + ruralstart].data()) * 1000;
			table.emissionsByHPMSUrban[ihpmsType][fuelcode] = table.emissionsByHPMSUrban[ihpmsType][fuelcode] + atof(splits[ihpmsType + urbanstart].data()) * 1000;
		}

	}
	for (size_t ihpmsType = 0; ihpmsType < 7; ihpmsType++)
	{
		int totalIdx = (int)FuelTypeCode::total;
		for (size_t ifuel = 0; ifuel < 2; ifuel++)
		{
			table.emissionsByHPMSRural[ihpmsType][totalIdx] = table.emissionsByHPMSRural[ihpmsType][totalIdx] + table.emissionsByHPMSRural[ihpmsType][ifuel];
			table.emissionsByHPMSUrban[ihpmsType][totalIdx] = table.emissionsByHPMSUrban[ihpmsType][totalIdx] + table.emissionsByHPMSUrban[ihpmsType][ifuel];
		}

	}
	infile.close();
	return table;
}
void EMFAC_FuelSplit::allocate(std::string onroadfile, std::vector<std::string> EMFACTableFILEs, std::vector<std::string> yearnames, std::string cityname)
{
	
	std::vector<EMFAC_FuelSpecificEmissions> EMFACTables;
	for (size_t i = 0; i < EMFACTableFILEs.size(); i++)
	{
		EMFACTables.push_back(loadEMFAC(EMFACTableFILEs[i],cityname));
	}
	std::string fuels[3] = { "g","d","" };
	ShapeFile shp(onroadfile, 1);
	OGRFeature *fea;
	//std::vector<
	int id = 0;
	int hpmsfield = shp.getOrCreateField("hpms", OGRFieldType::OFTString);
	int scagfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("rdTypeUR");
	int vmtfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("VMT");
	int cafield = shp.poLayer->GetLayerDefn()->GetFieldIndex("ca");
	std::vector<int> fieldindices;
	for (size_t iyear = 0; iyear < yearnames.size(); iyear++)
	{
		EMFAC_FuelSpecificEmissions& yeartable = EMFACTables[iyear];
		

		for (size_t ifuel = 0; ifuel < 3; ifuel++)
		{
			std::stringstream canamess;
			canamess << "ca" << yearnames[iyear];
			if (fuels[ifuel] != "")
				canamess << "_" << fuels[ifuel];
			std::string fuelname = canamess.str().data();
			printf("%s\n", fuelname.data());
			fieldindices.push_back(shp.getOrCreateField(fuelname.data(), OGRFieldType::OFTReal));
		}
	}




	double ruralVMTSum[7];
	double urbanVMTSum[7];
	for (size_t i = 0; i < 7; i++)
	{
		ruralVMTSum[i] = 0;
		urbanVMTSum[i] = 0;
	}
	std::map<int, OGRPoint> pointmaps;
	shp.poLayer->ResetReading();
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		std::string hpmsType = fea->GetFieldAsString(hpmsfield);
		double	vmt = fea->GetFieldAsDouble(vmtfield);
		int hpmsId = atoi(hpmsType.substr(0, 1).data()) - 1;
		if (hpmsType.substr(1, 1) == "U"){
			urbanVMTSum[hpmsId] += vmt;
		}
		else{
			ruralVMTSum[hpmsId] += vmt;
		}
		OGRFeature::DestroyFeature(fea);
	}

	shp.poLayer->ResetReading();
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		std::string hpmsType = fea->GetFieldAsString(hpmsfield);
		double vmt = fea->GetFieldAsDouble(vmtfield);
		int hpmsId = atoi(hpmsType.substr(0, 1).data()) - 1;
		double frac = 1;
		int fieldidx = 0;
		if (hpmsType.substr(1, 1) == "U") {
			frac = vmt / urbanVMTSum[hpmsId];

			fieldidx = 0;
			for (size_t iyear = 0; iyear < yearnames.size(); iyear++)
			{
				EMFAC_FuelSpecificEmissions& yeartable = EMFACTables[iyear];

				for (size_t ifuel = 0; ifuel < 3; ifuel++)
				{
					fea->SetField(fieldindices[fieldidx], frac*yeartable.emissionsByHPMSUrban[hpmsId][ifuel]);
					fieldidx++;
				}

			}
		}
		else
		{
			fieldidx = 0;
			frac = vmt / ruralVMTSum[hpmsId];
			for (size_t iyear = 0; iyear < yearnames.size(); iyear++)
			{
				EMFAC_FuelSpecificEmissions& yeartable = EMFACTables[iyear];

				for (size_t ifuel = 0; ifuel < 3; ifuel++)
				{
					fea->SetField(fieldindices[fieldidx], frac*yeartable.emissionsByHPMSRural[hpmsId][ifuel]);
					fieldidx++;
				}

			}
		}
		shp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}

	shp.close();
	printf("num of points:%d\n", id);
}
void EMFAC_FuelSplit::compare(std::string onroadfile1, std::string onroadfile2, std::string resultfile)
{


}
