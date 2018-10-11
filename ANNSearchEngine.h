#pragma once
#include <Windows.h>
#include <vector>
#include "ann/ANN.h"
#include "io.h"
#include "ogrsf_frmts.h"

//邻近点快速搜索方法
class ANNSearchEngine
{
public:
     ANNSearchEngine(void);
     ~ANNSearchEngine(void);

	 //************************************
	 // Method:    Create初始化
	 // FullName:  ANNSearchEngine::Create
	 // Access:    public 
	 // Returns:   bool
	 // Qualifier:
	 // Parameter: std::vector<GeoPoint> & points输入点数组
	 // Parameter: int maxReturns每次搜索返回邻近点个数
	 //************************************
	 bool Create(std::vector<OGRPoint>& points,int maxReturns);
	 int Select_Nearest_Points(const OGRPoint &p);
     //************************************
     // Method:    Get_Selected_Point获取第i个邻近点
     // FullName:  ANNSearchEngine::Get_Selected_Point
     // Access:    public 
     // Returns:   bool
     // Qualifier:
     // Parameter: int i索引
     // Parameter: double & x坐标X
     // Parameter: double & y坐标Y
     // Parameter: double & z属性值
     // Parameter: double & dist距离的平方
     //************************************
     bool Get_Selected_Point(int i, double &x, double &y,double &dist);
	 bool Get_Selected_Point(int i, int& idx,double &dist);
private:
	 //************************************
	 // Method:    Destroy
	 // FullName:  ANNSearchEngine::Destroy
	 // Access:    public 
	 // Returns:   void
	 // Qualifier:
	 // Parameter: void
	 //************************************
	 void	Destroy(void);
private:
	int m_nPoints_Max;
	int					nPts;					// actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances
	ANNkd_tree*			kdTree;					// search structure
	//std::vector<double>             valueArray;
	double*          valueArray;
	int              arraySize;
	double			eps;			// error bound
	bool m_isCreated;

};

