#pragma once
#include "ShapeFile.h"
class EMFAC
{
public:
	EMFAC();
	~EMFAC();
	void HPMS2SCAG(std::string onroadfile, std::string ruralCrosswalk, std::string urbanCrosswalk);
	void allocate(std::string onroadfile, std::string EMFACTableFILE);
	void compare(std::string onroadfile1, std::string onroadfile2, std::string resultfile);
};

