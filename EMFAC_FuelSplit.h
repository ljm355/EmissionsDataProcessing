#pragma once
#include "ShapeFile.h"
struct EMFAC_FuelSpecificEmissions
{

	std::string fieldname;
	int fieldidx;
	int year;
	//double emissionsByHPMSRural[7];
	//double emissionsByHPMSUrban[7];

	std::vector<double> emissionsByHPMSRural[7];
	std::vector<double> emissionsByHPMSUrban[7];
	EMFAC_FuelSpecificEmissions()
	{
		for (size_t i = 0; i < 7; i++)
		{
			for (int ifuel = 0; ifuel < 4; ifuel++)
			{
				emissionsByHPMSRural[i].push_back(0);
				emissionsByHPMSUrban[i].push_back(0);
			}
		}
	}
	EMFAC_FuelSpecificEmissions(std::string name)
	{
		fieldname = "ca" + name.substr(name.size() - 2, 2);
		year = atoi(name.substr(5, name.size() - 5).data());
	}


};
class EMFAC_FuelSplit
{
public:
	EMFAC_FuelSplit();
	~EMFAC_FuelSplit();
	void HPMS2SCAG(std::string onroadfile, std::string ruralCrosswalk, std::string urbanCrosswalk);
	void allocate(std::string onroadfile, std::vector<std::string> EMFACTableFILEs, std::vector<std::string> yearnames, std::string cityname);
	void compare(std::string onroadfile1, std::string onroadfile2, std::string resultfile);
	void loadSCAG2HPMS(std::string crosswalkfile, std::string postfix, std::map<std::string, std::string>& crosswalk);
	EMFAC_FuelSpecificEmissions loadEMFAC(std::string EMFACTableFILE, std::string cityname = "");
	enum FuelTypeCode
	{

		gasoline = 0,
		diesel = 1,
		//electric = 2,
		total = 2,
		unknown = 1000
	};
	FuelTypeCode getFuelCode(std::string name)
	{
		if (name == "gasoline")
			return FuelTypeCode::gasoline;
		if (name == "diesel")
			return FuelTypeCode::diesel;
		//if (name == "electric")
		//	return -1;// FuelTypeCode::electric;
		if (name == "total")
			return FuelTypeCode::total;
		return FuelTypeCode::unknown;
	}
};

