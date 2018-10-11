
#include "ANNSearchEngine.h"
#include "ann\ANN.h"					// ANN declarations

ANNSearchEngine::ANNSearchEngine(void)
{
	valueArray = NULL;
	nnIdx = NULL;
	dists = NULL;
	kdTree = NULL;
	eps = 0.01;
}
ANNSearchEngine::~ANNSearchEngine(void)
{
	Destroy();
}
bool ANNSearchEngine::Create(std::vector<OGRPoint>& points, int maxPoints)
{
	Destroy();
	m_nPoints_Max = maxPoints;
	nPts = points.size();
	queryPt = annAllocPt(2);					// allocate query point
	dataPts = annAllocPts(nPts, 2);			// allocate data points
	nnIdx = new ANNidx[maxPoints];						// allocate near neigh indices
	dists = new ANNdist[maxPoints];						// allocate near neighbor dists
	/*if (valueArray)
	{
		delete[] valueArray;
		valueArray = NULL;
	}*/
	//valueArray = new double[points.size()];
	for (int i = 0; i<nPts; i++)
	{
		OGRPoint& shp = points[i];
		dataPts[i][0] = shp.getX();
		dataPts[i][1] = shp.getY();
		//valueArray[i] = shp.z;

	}

	kdTree = new ANNkd_tree(					// build search structure
		dataPts,					// the data points
		nPts,						// number of points
		2);

	return true;
}
int ANNSearchEngine::Select_Nearest_Points(const OGRPoint &p)
{
	queryPt[0] = p.getX();
	queryPt[1] = p.getY();
	kdTree->annkSearch(						// search
		queryPt,						// query point
		m_nPoints_Max,								// number of near neighbors
		nnIdx,							// nearest neighbors (returned)
		dists,							// distance (returned)
		eps);							// error bound
	return m_nPoints_Max;
}
bool ANNSearchEngine::Get_Selected_Point(int i, int& idx, double &dist)
{
	int id = nnIdx[i];
	idx = id;
	double x = dataPts[id][0];
	double y = dataPts[id][1];
	//z = valueArray[id];
	dist = dists[i];
	return true;
}
void	ANNSearchEngine::Destroy(void)
{
	if (nnIdx)
	{
		delete[] nnIdx;							// clean things up
		nnIdx = NULL;
	}
	if (dists)
	{
		delete[] dists;
		dists = NULL;
	}
	/*if (valueArray)
	{
		delete[]valueArray;
		valueArray = NULL;
	}*/
	if (kdTree)
	{
		delete kdTree;
		kdTree = NULL;
	}
	annClose();
}
