#include "SparseFractionGrid.h"

SparseFractionGrid::SparseFractionGrid()
{

}


SparseFractionGrid::~SparseFractionGrid()
{

}

void SparseFractionGrid::save(std::string filename)
{
	std::ofstream ofs;
	ofs.open(filename, std::ofstream::binary);
	int num = cells.size();
	ofs.write((char*)(&num), sizeof(int));
	size_t bufsize = sizeof(SparseGridCell) * cells.size();
	ofs.write((char*)(&cells[0]), bufsize);
	ofs.close();
}

void SparseFractionGrid::open(std::string filename)
{
	std::ifstream ifs;
	ifs.open(filename, std::ofstream::binary);
	int num = 0;
	ifs.read((char*)(&num), sizeof(int));
	cells.resize(num);
	size_t bufsize = sizeof(SparseGridCell) * cells.size();
	ifs.read((char*)(&cells[0]), bufsize);
	ifs.close();
}
