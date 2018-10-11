#pragma once
#include <Windows.h>
#include <string>
#include <fstream>
class PeMSStation
{
public:
	PeMSStation();
	~PeMSStation();
	double* m_Fractions;
	int   m_NumHours;
	int m_ID;
	std::string m_Header;
	double m_X;
	double m_Y;
	std::string m_HestiaType;
	int roadclass;
	void fromStream(std::ifstream& ifs);
	void toBin(std::string filename);
	void toTxt(std::string filename);
};


