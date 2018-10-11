#include "NeiborhoodOptimizationMatrixPrep.h"



NeiborOptimMatrixPrep::NeiborOptimMatrixPrep()
{


}


NeiborOptimMatrixPrep::~NeiborOptimMatrixPrep()
{


}
void NeiborOptimMatrixPrep::updateNEEUIClassCode(std::string filename, std::vector<std::string> neeui_filenames, std::vector<int> yearbreaks, int censusdivision)
{
	std::vector<std::map<int, NEEUIRecord>> tables;
	bool calElec = true;
	for (size_t i = 0; i < neeui_filenames.size(); i++)
	{
		tables.push_back(NonPointProcessor::loadNEEUITable(neeui_filenames[i], censusdivision));
	}
	ShapeFile shp(filename, 1);
	int neeui_ng_idx = shp.getOrCreateField("neeui_ng", OGRFieldType::OFTReal);
	int neeui_kwh_idx = shp.getOrCreateField("neeui_kwh", OGRFieldType::OFTReal);
	int kwh_idx = shp.getOrCreateField("kwh", OGRFieldType::OFTReal);;
	int tfs_idx = shp.getOrCreateField("tfs", OGRFieldType::OFTReal);
	int year_idx = shp.getOrCreateField("yearbuilt", OGRFieldType::OFTReal);
	int bt_idx = shp.getOrCreateField("bt", OGRFieldType::OFTReal);
	int neeui_type = shp.getOrCreateField("neeui_tp", OGRFieldType::OFTInteger);
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
				breakidx = i + 1;
				break;
			}
		}
		if (breakidx > tables.size() - 1) {
			breakidx = tables.size() - 1;
		}
		std::map<int, NEEUIRecord>& table = tables[breakidx];
		NEEUIRecord& record = table[bt];
		double neeui_ng = record.NEEUING;
		double neeui_kwh = record.NEEUIElecKWH;
		double neeui_ng2 = fea->GetFieldAsDouble(neeui_ng_idx);
		double neeui_kwh2 = fea->GetFieldAsDouble(neeui_kwh_idx);
		fea->SetField(neeui_kwh_idx, neeui_kwh);
		fea->SetField(kwh_idx, neeui_kwh * tfs);
		fea->SetField(neeui_type, record.index);

		shp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}

	shp.close();

}