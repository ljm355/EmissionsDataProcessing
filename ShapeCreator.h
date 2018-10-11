#pragma once
#include "ShapeFile.h"
#include <string>
#include <qstring.h>
#include <vector>

struct ValueField
{
	int index;
	int index2;
	std::string s;
	double d;
	int i;
	OGRFieldType tp;
	std::string srcname;
	std::string srcname2;
	std::string destname;
	ValueField(OGRFieldType otp, std::string destfieldname, std::string firstfieldname, std::string secondfieldname = "")
	{
		srcname = firstfieldname;
		srcname2 = secondfieldname;
		tp = otp;
		index = -1;
		index2 = -1;
		destname = destfieldname;
	}


	bool set(std::string& val)
	{
		if (val == "")
		{
			s = "";
			d = 0;
			i = 0;
			return false;
		}

		if (tp == OGRFieldType::OFTString)
			s = val;
		else if (tp == OGRFieldType::OFTReal)
			d = atof(val.data());
		else if (tp == OGRFieldType::OFTInteger)
			i = atoi(val.data());
		return true;
	}
	void set(std::vector<std::string>& vals)
	{
		if (!set(vals[index]) && index2 > -1)
		{
			set(vals[index2]);
		}
	}
	void get(OGRFeature* fea)
	{
		if (tp == OGRFieldType::OFTString)
			fea->SetField(destname.data(), s.data());
		else if (tp == OGRFieldType::OFTReal)
			fea->SetField(destname.data(), d);
		else if (tp == OGRFieldType::OFTInteger)
			fea->SetField(destname.data(), i);
	}
	void findIndex(std::vector<std::string>& vals)
	{
		for (size_t i = 0; i < vals.size(); i++)
		{
			if (srcname == vals[i])
			{
				index = i;
				break;
			}
		}
		if (srcname2 != "")
		{
			for (size_t i = 0; i < vals.size(); i++)
			{
				if (vals[i] != "" && srcname2 == vals[i])
				{
					index2 = i;
					break;
				}
			}
		}


	}
	void createField(OGRLayer* layer)
	{
		OGRFieldDefn field(destname.data(), tp);
		layer->CreateField(&field);
	}
	bool equal(ValueField* field)
	{
		if (tp == OGRFieldType::OFTString)
			return s == field->s;
		else if (tp == OGRFieldType::OFTReal)
			return d == field->d;
		else if (tp == OGRFieldType::OFTInteger)
			return i == field->i;
		return false;
	}

};
struct ShapeFileBySector
{
	ShapeFile* SHP;
	std::string BC;
	int SecID;
public:
	ShapeFileBySector(ShapeFile* shp, std::string bc,int id = -1)
	{
		SHP = shp;
		BC = bc;
		SecID = id;
	}
};
class ShapeCreator
{
public:
	ShapeCreator();
	~ShapeCreator();
	enum SECTORID
	{
		ELE = 1,
		IND = 2,
		COM = 3,
		NONROAD = 6,
		RAILROAD = 7,
		AIRPORT = 9,
	};

	std::vector<std::string> split(const char & delimiter, const std::string & line);
	void split(const char & delimiter, const std::string & line, std::string * psplits);
	void createBySector(std::string configFile, std::string dir);
	void create(std::string configFile, std::string outfile);
	void create(std::string csvfielname, std::string outfilename, std::vector<ValueField> fields, ValueField lonfield, ValueField latfield, ValueField filterfield);
};

