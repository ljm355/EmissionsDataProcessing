#pragma once
#include <string>
class MarionElementa
{
public:
	MarionElementa();
	~MarionElementa();
	void getSpatialRefInfo(std::string reffile, std::string & proj, double* adftransform);
	float * readAttribute(std::string filename, std::string atrname);
	void crop(std::string filename, std::string maskfile, int mask,std::string outfile);
	void dif(std::string filename, std::string outfile, std::string outcsv);
	void nc2tif(std::string ncfile,std::string tiffile,std::string reffile = "B:/Hestia_FFDAS_ODIAC/Influx/MarionMask.tif");
	float* readdata(std::string filename, int band);
};

