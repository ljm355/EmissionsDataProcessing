#include "EMFAC.h"
#include <fstream>
#include <map>
#include <sstream>
#include "Utils.h"
EMFAC::EMFAC()
{

}


EMFAC::~EMFAC()
{

}

void loadSCAG2HPMS(std::string crosswalkfile,std::string postfix,std::map<std::string,std::string>& crosswalk)
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


void EMFAC::HPMS2SCAG(std::string onroadfile, std::string ruralCrosswalk, std::string urbanCrosswalk)
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

struct EMFAC_Emissions
{

	std::string fieldname;
	int fieldidx;
	int year;
	double emissionsByHPMSRural[7];
	double emissionsByHPMSUrban[7];
	EMFAC_Emissions(std::string name)
	{
		fieldname = "ca" + name.substr(name.size() - 2, 2);
		year = atoi(name.substr(5, name.size() - 5).data());
	}


};
std::vector<EMFAC_Emissions> loadEMFAC(std::string EMFACTableFILE)
{
	std::vector<EMFAC_Emissions> table;
	std::string line;
	std::ifstream infile(EMFACTableFILE.data());
	std::getline(infile, line);
	std::vector<std::string> splits = Utils::splitCSV(',', line);
	for (size_t i = 0; i < (splits.size() - 1) / 2; i++)
	{
		std::string header = splits[i * 2+1];
		EMFAC_Emissions emfacYear(header);
		table.push_back(emfacYear);
	}
	int hpmsType = 0;
	while (std::getline(infile, line))
	{
		splits = Utils::splitCSV(',', line);
		for (size_t i = 0; i < (splits.size() - 1) / 2; i++)
		{
			table[i].emissionsByHPMSRural[hpmsType] = atof(splits[i * 2 + 1].data()) * 1000;
			table[i].emissionsByHPMSUrban[hpmsType] = atof(splits[i * 2 + 1 + 1].data()) * 1000;
		}
		hpmsType++;
	}
	infile.close();
	return table;
}
void EMFAC::allocate(std::string onroadfile, std::string EMFACTableFILE)
{
	
	std::vector<EMFAC_Emissions> EMFACTable = loadEMFAC(EMFACTableFILE);

	ShapeFile shp(onroadfile, 1);
	OGRFeature *fea;

	int id = 0;
	int hpmsfield = shp.getOrCreateField("hpms", OGRFieldType::OFTString);
	int scagfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("rdTypeUR");
	int vmtfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("VMT");
	int cafield = shp.poLayer->GetLayerDefn()->GetFieldIndex("ca");
	for (size_t i = 0; i < EMFACTable.size(); i++)
	{
		EMFACTable[i].fieldidx = shp.getOrCreateField(EMFACTable[i].fieldname.data(), OGRFieldType::OFTReal);
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
		if (hpmsType.substr(1, 1) == "U") {
			frac = vmt / urbanVMTSum[hpmsId];
			for (size_t iyear = 0; iyear < EMFACTable.size(); iyear++)
			{
				fea->SetField(EMFACTable[iyear].fieldidx, frac*EMFACTable[iyear].emissionsByHPMSUrban[hpmsId]);
			}
		}
		else
		{
			frac = vmt / ruralVMTSum[hpmsId];
			for (size_t iyear = 0; iyear < EMFACTable.size(); iyear++)
			{
				fea->SetField(EMFACTable[iyear].fieldidx, frac*EMFACTable[iyear].emissionsByHPMSRural[hpmsId]);
			}
		}
		shp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}

	shp.close();
	printf("num of points:%d\n", id);
}
void EMFAC::compare(std::string onroadfile1, std::string onroadfile2, std::string resultfile)
{


}
