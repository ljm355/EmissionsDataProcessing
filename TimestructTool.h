#pragma once
#include <string>
#include <vector>
class TimestructTool
{
public:
	TimestructTool();
	~TimestructTool();
	static bool isLeapYear(int year);
	static int numOfHours(int year);
	static std::vector<double> readText(std::string filename, std::string& id, bool normalize = true);
	static std::vector<double> readBin(std::string filename, std::string & id, bool normalize);
	static void readText(std::string filename, std::vector<double>& fracs, std::string& id, bool normalize = true);
	static void normalizeFractions(double *& fracs,int numhours);
	static void toText(std::string filename, double* fracs, int numhours);
	static void binary2Text(std::string binaryfilename, std::string outdir, std::string ext);
	static void updateBinary(std::string infile, std::string outfile, std::string txtdir, bool normalize = true);
	static void updateFileInBinary(std::string infile, std::string outfile, std::string txtfile, bool normalize = true);
	static void normalizeBinary(std::string infile);
	static void timeshift(double* srcfracs,int srcyear, double*& destfracs,int destyear, int destNumHours);
	static double* timeshift(double* srcfracs, int srcyear, int destyear);
	static void timeshiftOverwrite(double * srcfracs, int srcyear, int destyear);
	static void txt2binary(int numhours, std::string txtdir, std::string binaryfilename,  bool normalize = true);
	static void mergeBinary(int numhours, std::string bindir, std::string binaryfilename, bool normalize);
	static void timeshiftbinary(std::string infile, std::string outfile, int srcyear, int destyear, bool normalize = true);
	static double* readBinary(std::string binaryfilename, std::string id, int& numhours);
	static void selectFromBinary2Text(std::string binaryfilename, std::string id, std::string txtfilename);
	static void padLastDay(double*& fracs);
private:

};

