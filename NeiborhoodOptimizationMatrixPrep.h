#pragma once
#include "NonPointProcessor.h"

class NeiborOptimMatrixPrep
{
public:
	NeiborOptimMatrixPrep();
	~NeiborOptimMatrixPrep();
	void updateNEEUIClassCode(std::string filename, std::vector<std::string> neeui_filenames, std::vector<int> yearbreaks, int censusdivision);
};

