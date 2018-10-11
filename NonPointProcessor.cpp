#include "NonPointProcessor.h"
#include "qdir.h"
#include <fstream>
#include <sstream>
#include "Utils.h"
NonPointProcessor::NonPointProcessor()
{
}


NonPointProcessor::~NonPointProcessor()
{

}

std::map<int, NEEUIRecord> NonPointProcessor::loadNEEUITable(std::string filename, int censusdivision)
{

	std::map<int, NEEUIRecord> table;
	std::string line;
	std::ifstream infile(filename.data());
	std::getline(infile, line);
	std::vector<std::string> splits = Utils::splitCSV(',', line);
	int DivisionFidx = -1;
	int BTFidx = -1;
	int NEEUINGF = -1;
	int NEEUIFuelOilF = -1;
	int NEEUIElec = -1;
	int NEEUIElecKWH = -1;
	for (size_t i = 0; i < splits.size(); i++)
	{
		std::string split = splits[i];
		if (split == "Division")
			DivisionFidx = i;
		else if (split == "BT")
			BTFidx = i;
		else if (split == "NEEUING")
			NEEUINGF = i;
		else if (split == "NEEUIFuelOil")
			NEEUIFuelOilF = i;
		else if (split == "NEEUIElec")
			NEEUIElec = i;
		else if (split == "NEEUIElecKWH")
			NEEUIElecKWH = i;

	}
	int hpmsType = 0;
	int index = 1;
	while (std::getline(infile, line))
	{
		splits = Utils::splitCSV(',', line);
		if (DivisionFidx > -1 && atoi(splits[DivisionFidx].data()) != censusdivision)
			continue;
		int bt = atoi(splits[BTFidx].data());
		NEEUIRecord record;
		record.NEEUING = atof(splits[NEEUINGF].data());
		record.NEEUIFuelOil = atof(splits[NEEUIFuelOilF].data());
		if(NEEUIElec > -1)
		   record.NEEUIElec = atof(splits[NEEUIElec].data());
		if(NEEUIElecKWH > -1)
		   record.NEEUIElecKWH = atof(splits[NEEUIElecKWH].data());
		record.index = index;
		table[bt] = record;
		index++;
	}
	infile.close();
	return table;
}

void NonPointProcessor::updateNEEUI(std::string filename,std::vector<std::string> neeui_filenames,std::vector<int> yearbreaks,int censusdivision)
{
	std::vector<std::map<int, NEEUIRecord>> tables;
	bool calElec = true;
	for (size_t i = 0; i < neeui_filenames.size(); i++)
	{
		tables.push_back(loadNEEUITable(neeui_filenames[i], censusdivision));
	}
	ShapeFile shp(filename, 1);
	int neeui_ng_idx = shp.getOrCreateField("neeui_ng",OGRFieldType::OFTReal);
	int neeui_p_idx = shp.getOrCreateField("neeui_p", OGRFieldType::OFTReal);

	//int neeui_c_idx = shp.getOrCreateField("neeui_c", OGRFieldType::OFTReal);
	int ca_ng_idx = shp.getOrCreateField("ca_ng", OGRFieldType::OFTReal);
	int ca_p_idx = shp.getOrCreateField("ca_p", OGRFieldType::OFTReal);
	int ca_c_idx = shp.getOrCreateField("ca_c", OGRFieldType::OFTReal);

	int neeui_elec_idx = -1;
	int neeui_kwh_idx = -1;
	int elec_idx = -1;
	int kwh_idx = -1;
	if (calElec)
	{
		neeui_elec_idx = shp.getOrCreateField("neeui_elec", OGRFieldType::OFTReal);
		neeui_kwh_idx = shp.getOrCreateField("neeui_kwh", OGRFieldType::OFTReal);
		elec_idx = shp.getOrCreateField("elec", OGRFieldType::OFTReal);
		kwh_idx = shp.getOrCreateField("kwh", OGRFieldType::OFTReal);

	}

	int tfs_idx = shp.getOrCreateField("tfs", OGRFieldType::OFTReal);
	int year_idx = shp.getOrCreateField("yearbuilt", OGRFieldType::OFTReal);
	int bt_idx = shp.getOrCreateField("bt", OGRFieldType::OFTReal);
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		double tfs = fea->GetFieldAsDouble(tfs_idx);

		int year = fea->GetFieldAsDouble(year_idx);
		int bt = fea->GetFieldAsDouble(bt_idx);
		int breakidx = 0;
		for (size_t i = 0; i < yearbreaks.size(); i++)
		{
			if (year >= yearbreaks[i])
			{
				breakidx = i+1;
				break;
			}
		}
		if (breakidx > tables.size() - 1) {
			breakidx = tables.size() - 1;
		}
		std::map<int, NEEUIRecord>& table = tables[breakidx];
		NEEUIRecord& record = table[bt];

		double neeui_ng = record.NEEUING; 
		double neeui_p = record.NEEUIFuelOil;

		double neeui_elec = record.NEEUIElec;
		double neeui_kwh = record.NEEUIElecKWH;

		fea->SetField(neeui_ng_idx, neeui_ng);
		fea->SetField(neeui_p_idx, neeui_p);


		//fea->SetField(neeui_c_idx, tfs);

		fea->SetField(ca_ng_idx, neeui_ng * tfs);
		fea->SetField(ca_p_idx, neeui_p * tfs);
		fea->SetField(ca_c_idx, tfs);
		if (calElec)
		{
			fea->SetField(neeui_elec_idx, neeui_elec);
			fea->SetField(neeui_kwh_idx, neeui_kwh);
			fea->SetField(elec_idx, neeui_elec * tfs);
			fea->SetField(kwh_idx, neeui_kwh * tfs);
		}
		shp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}

	shp.close();

}
//
//fuel	year	sector	total	carbon_conversion	building	division	output_field
//ng	2020	2	204878287	0.014481818	B: / LA_Version2 / gridPrep_SHP_master / Los_Angeles / ResNonPoint.dbf	9	ca_ng
//petrol	2020	2	40092210	0.019773792	B : / LA_Version2 / gridPrep_SHP_master / Los_Angeles / ResNonPoint.dbf	9	ca_p
//coal	2020	2	20330	1	B : / LA_Version2 / gridPrep_SHP_master / Los_Angeles / ResNonPoint.dbf	9	ca_c
void NonPointProcessor::toCSV(std::vector<SpatialAllocationParams> params, std::string shapefilename, std::string outfile)
{
	std::ofstream ofs;
	ofs.open(outfile);
	ofs << "fuel," << "year," << "sector," << "total," << "building," << "division," << "output_field" << "\n";
	std::string fuels[] = {"ng","petrol","coal","elec"};
	for (size_t i = 0; i < params.size(); i++)
	{
		SpatialAllocationParams p = params[i];
		ofs << fuels[p.fuel-1] << "," << p.year << "," << p.sector << "," << p.total << "," << shapefilename << "," << p.division << "," << p.output_field << "\n";
	}
	ofs.close();
}
void NonPointProcessor::runRScript(std::string scriptfile, std::string sharedCFG, std::string shapefilename, std::vector<SpatialAllocationParams> params)
{
	std::string tmpfile = sharedCFG + ".csv";
	toCSV(params, shapefilename,tmpfile);
	std::stringstream ss;
	ss << "RScript " << scriptfile << " " << sharedCFG << " " << tmpfile <<"\n";
	system(ss.str().data());
	QFile::remove(tmpfile.data());
}
//void NonPointProcessor::spatialAllocation(std::string sharedCFG, std::string shapefile, std::string sectorField)
//{
//
//}
void NonPointProcessor::exportBySector(std::string infile, std::string outdir, std::string sectorField)
{

	ShapeFile NonPoint(infile);
	outdir = QDir(outdir.data()).absolutePath().toLocal8Bit().data() + std::string("/");
	ShapeFile ComNonPoint;
	ComNonPoint.create(outdir + "ComNonPoint.shp", NonPoint.poLayer->GetSpatialRef(), NonPoint.poLayer->GetLayerDefn(), NonPoint.poLayer->GetGeomType());

	ShapeFile ResNonPoint;
	ResNonPoint.create(outdir + "ResNonPoint.shp", NonPoint.poLayer->GetSpatialRef(), NonPoint.poLayer->GetLayerDefn(), NonPoint.poLayer->GetGeomType());

	ShapeFile IndNonPoint;
	IndNonPoint.create(outdir + "IndNonPoint.shp", NonPoint.poLayer->GetSpatialRef(), NonPoint.poLayer->GetLayerDefn(), NonPoint.poLayer->GetGeomType());

	OGRLayer* layers[] = { ComNonPoint.poLayer,ResNonPoint.poLayer,IndNonPoint.poLayer };
	//IndPoint.create((dir + "IndPoint.shp").data(), NULL, NULL, OGRwkbGeometryType::wkbPoint);
	//shpBySectors.push_back(ShapeFileBySector(&IndPoint, "IND", 2));
	OGRFeature *poFeature;
	double sum = 0;
	int bcidx = NonPoint.poLayer->GetLayerDefn()->GetFieldIndex(sectorField.data());
	NonPoint.poLayer->ResetReading();
	while ((poFeature = NonPoint.poLayer->GetNextFeature()) != NULL)
	{
		int bc = poFeature->GetFieldAsInteger(bcidx);
		if (bc < 1 || bc > 3)
		{
			continue;
			OGRFeature::DestroyFeature(poFeature);
		}
		OGRLayer* layer = layers[bc - 1];
		layer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}


}

