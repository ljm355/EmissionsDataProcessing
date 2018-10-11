#pragma once
#include "ShapeFile.h"
#include <vector>
#include "gdal_priv.h"
#include "Utils.h"
enum OperatorType
{
	SUM = 0,
	AVERAGE = 1
};



class OGRFieldData
{
public:
	std::vector<double> data;
	OGRFieldType tp;
	std::string srcName;
	std::string destName;
	int srcIdx;
	int destIdx;
	OperatorType operatorType;
	std::vector<int> count;
	OGRFieldData()
	{
		operatorType = OperatorType::SUM;
	}
	OGRFieldData(std::string _srcname,std::string _destname, OperatorType _tp)
	{
		srcName = _srcname; destName = _destname; operatorType = _tp;
	}
	~OGRFieldData()
	{

	}
	void begin()
	{
		count.resize(data.size());
		memset(&count[0], 0, count.size() * sizeof(int));
	}
	void set(double val,int idx)
	{
		if (val <= 0)
			return;
		data[idx] += val;
		count[idx] ++;
	}
	void end()
	{
		if (operatorType == OperatorType::AVERAGE)
		{
			for (int i = 0; i < data.size(); i++)
			{
				data[i] = data[i] / count[i];
			}
		}
	}
	void setField(OGRFeature* fea,int idx)
	{
		if (tp == OGRFieldType::OFTReal)
			fea->SetField(destIdx,data[idx]);
		else if	(tp == OGRFieldType::OFTInteger64)
			fea->SetField(destIdx, (long)data[idx]);
		else
			fea->SetField(destIdx, (int)data[idx]);
	}
};


class Geographies
{
public:
	void setup(std::string srcShpFile, std::vector<OGRFieldData> fields, std::string destShpFile);
	void gather(std::string srcShpFile, std::string fidName, OperatorType tp = OperatorType::SUM);
	void update();
private:
	double calFraction(OGRFeature* fea, const int& footprintIdx);
	std::vector<OGRFieldData> m_AttributeData;
	std::string m_destFileName;
public:

};

