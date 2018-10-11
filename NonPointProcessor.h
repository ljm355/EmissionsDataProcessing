#pragma once
#include "ScaleByFuel.h"

struct NEEUIRecord
{
	double NEEUING;
	double NEEUIFuelOil;
	double NEEUIElec;
	double NEEUIElecKWH;
	int index;
};


class NonPointProcessor
{
public:

	NonPointProcessor();
	~NonPointProcessor();
	enum NonPointSector
	{
		COM = 1,
		RES = 2,
		IND = 3
	};

	enum NonPointFuel
	{
		NG = 1,
		PETROL = 2,
		COAL = 3,
		ELEC = 4
	};
	struct SpatialAllocationParams
	{
		NonPointFuel fuel;
		int year;
		int sector;
		int division;
		double total;
		std::string output_field;
		SpatialAllocationParams()
		{

		}
		SpatialAllocationParams(NonPointFuel _fuel, int _year, int _sector, int _division, double _total, std::string _fieldname)
			:fuel(_fuel),year(_year),sector(_sector),division(_division),total(_total), output_field(_fieldname)
		{

		}
	};
	static std::map<int, NEEUIRecord> NonPointProcessor::loadNEEUITable(std::string filename, int censusdivision);
	void exportBySector(std::string infile, std::string outdir, std::string sectorField = "bc");
	static void updateNEEUI(std::string filename, std::vector<std::string> neeui_filenames, std::vector<int> yearbreaks, int censusdivision);
	void toCSV(std::vector<SpatialAllocationParams> params, std::string shapefilename, std::string outfile);
	void runRScript(std::string scriptfile, std::string sharedCFG, std::string shapefilename, std::vector<SpatialAllocationParams> params);
};

