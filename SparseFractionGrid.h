#pragma once
#include <fstream>
#include <map>
#include <vector>
struct SparseGridCell
{
	int feaId;
	int gridId;
	double fraction;
	SparseGridCell()
	{

	}
	SparseGridCell(int fid,int gid,double frac)
	{
		feaId = fid; gridId = gid; fraction = frac;
	}
};

class SparseFractionGrid
{
public:
	std::vector<SparseGridCell> cells;
	SparseFractionGrid();
	~SparseFractionGrid();
	void save(std::string filename);
	void open(std::string filename);
};

