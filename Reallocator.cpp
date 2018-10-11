#include "Reallocator.h"
#include "qfileinfo.h"
#include <fstream>
#include <sstream>
#include "TimestructTool.h"
#include "qdir.h"
void CASector::loadAttribute(std::string _cafieldname, std::vector<std::string> _idfieldnames)
{
	cafieldname = _cafieldname;
	idfieldnames = _idfieldnames;
	features.clear();
	for (int ishapefile = 0; ishapefile < shapefiles.size(); ishapefile++)
	{
		std::string shapefile = shapefiles[ishapefile];
		ShapeFile input(shapefile.data());
		OGRFeature *poFeature;
		input.poLayer->ResetReading();

		int caIndex = input.poLayer->GetLayerDefn()->GetFieldIndex(cafieldname.data());
		std::vector<int> idIndices;
		for (size_t i = 0; i < _idfieldnames.size(); i++)
		{
			int idIndex = input.poLayer->GetLayerDefn()->GetFieldIndex(_idfieldnames[i].data());
			idIndices.push_back(idIndex);
		}

		int timestructIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("timestruct");
		OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
		while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
		{
			double ca = poFeature->GetFieldAsDouble(caIndex);
			if (ca > 0)
			{
				std::stringstream ss;
				for (size_t i = 0; i < idIndices.size(); i++)
				{
					ss << poFeature->GetFieldAsInteger(idIndices[i]);
				}
			
		
				IntersectedFeature fc;
				fc.total = ca;
				fc.id = atoi(ss.str().data());
				if (timestructIndex > -1 && timestructMap.size() > 0)
				{
					std::string tsid = poFeature->GetFieldAsString(timestructIndex);
					if (timestructMap.find(tsid) != timestructMap.end())
					{
						fc.timestructure = timestructMap[tsid];

					}
					else
					{
						fc.timestructure = timestructMap["0"];
						if (timestructMap.size() > 2)
						{
							printf("%s\n", tsid.data());
						}
					}
				}

				features.push_back(fc);
			}
			OGRFeature::DestroyFeature(poFeature);
		}
	}
}

std::vector<double> CASector::aggregrate(std::string basemapfile, std::string outfile)
{

	std::map<int, double> featuremap;
	for (size_t i = 0; i < features.size(); i++)
	{
		IntersectedFeature& fea = features[i];
		std::map<int, double>::iterator iter = featuremap.find(fea.id);
		if (iter != featuremap.end())
		{
			iter->second = iter->second + fea.total;
		}
		else
		{
			featuremap[fea.id] = 0;
		}
	}
	ShapeFile::copy(basemapfile, outfile);
	ShapeFile shp(outfile, 1);
	OGRFeature *poFeature;
	std::vector<double> cacolumn;
	int caidx = shp.getOrCreateField(cafieldname.data(),OGRFieldType::OFTReal);
	std::vector<int> idIndices;
	for (size_t i = 0; i < idfieldnames.size(); i++)
	{
		int idIndex = shp.poLayer->GetLayerDefn()->GetFieldIndex(idfieldnames[i].data());
		idIndices.push_back(idIndex);
	}
	shp.poLayer->ResetReading();
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		std::stringstream ss;
		for (size_t i = 0; i < idIndices.size(); i++)
		{
			ss << poFeature->GetFieldAsInteger(idIndices[i]);
		}
		int id = atoi(ss.str().data());
		std::map<int, double>::iterator iter = featuremap.find(id);
		double ca = 0;
		if (iter != featuremap.end())
		{
			ca = iter->second;
		}
		poFeature->SetField(caidx, ca);
		cacolumn.push_back(ca);
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}

	shp.close();
	return cacolumn;
}

Reallocator::Reallocator()
{
	fractionbuf = NULL;
	idbuf = NULL;
}
void Reallocator::loadAttribute(std::string _cafieldname, std::vector<std::string> _idfieldnames)
{
	cafieldname = _cafieldname;
	idfieldnames = _idfieldnames;
	for (int isector = 0; isector < sectors.size(); isector++)
	{
		CASector* sector = sectors[isector];
		sector->loadAttribute(_cafieldname, _idfieldnames);
	}
}


void Reallocator::addSector(std::vector<std::string> shapefiles, std::string sectorname)
{
	CASector* sector= new CASector;
	sector->timestructMap = timestructMap;
	sector->shapefiles = shapefiles;
	sector->sectorname = sectorname;
	sectors.push_back(sector);
}

std::vector<double> Reallocator::aggregrate(std::string basemapfile, std::string outfile)
{

	QDir qdir = QFileInfo(outfile.data()).absoluteDir();
	if (!qdir.exists())
		qdir.mkpath(".");
	std::string dir = qdir.absolutePath().toLocal8Bit().data();
	dir = dir + "/";
	std::vector<double> sum;

	for (size_t i = 0; i < sectors.size(); i++)
	{
		CASector*  sector = sectors[i];
		std::string sectoroutfile = dir + sector->sectorname + ".shp";
		std::vector<double> sectorsum = sector->aggregrate(basemapfile, sectoroutfile);
		if (sum.size() < 1)
		{
			sum = sectorsum;
		}
		else
		{
			for (size_t j = 0; j < sum.size(); j++)
			{
				sum[j] = sum[j] + sectorsum[j];
			}
		}
	}

	ShapeFile::copy(basemapfile, outfile);
	ShapeFile shp(outfile, 1);
	OGRFeature *poFeature;
	std::vector<double> cacolumn;
	int caidx = shp.getOrCreateField(cafieldname.data(), OGRFieldType::OFTReal);
	//std::vector<int> idIndices;
	//for (size_t i = 0; i < idfieldnames.size(); i++)
	//{
	//	int idIndex = shp.poLayer->GetLayerDefn()->GetFieldIndex(idfieldnames[i].data());
	//	idIndices.push_back(idIndex);
	//}
	shp.poLayer->ResetReading();
	int ifeature = 0;
	while ((poFeature = shp.poLayer->GetNextFeature()) != NULL)
	{
		poFeature->SetField(caidx, sum[ifeature]);
		shp.poLayer->SetFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
		ifeature++;
	}

	shp.close();
	return sum;

}

void Reallocator::loadtimestruct(std::string txtdir, std::string binaryfile,int numofhours)
{
	if (!QFileInfo(binaryfile.data()).exists())
	{
		TimestructTool::txt2binary(numofhours,txtdir, binaryfile);
	}

	if (idbuf)
	{
		delete[] idbuf;
		delete[] fractionbuf;
		timestructMap.clear();
	}
	std::ifstream filein(binaryfile.data(), std::ios::binary);
	int numstructs = 0;
	filein.read((char*)&numstructs, sizeof(int) * 1);
	idbuf = new char[numstructs * 20];
	memset(idbuf, 0, (size_t)numstructs * 20);
	fractionbuf = new double[numstructs * numofhours];
	filein.read(idbuf, numstructs * 20);
	filein.read((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numofhours);
	filein.close();
	//std::ofstream ofs("B:/LA_Version2/Time/2011.csv");
	char* pidbuf = idbuf;
	double* pfractionbuf = fractionbuf;
	//std::string outputdir = "C:/HestiaGridding/Baltimore/";
	for (size_t i = 0; i < numstructs; i++)
	{
		//ofs << pidbuf << std::endl;
		std::string tsid = pidbuf;
		/*double sum = 0;
		std::ofstream ofsts((outputdir + tsid + ".csv").data());
		for (size_t ihour = 0; ihour < 8760; ihour++)
		{
		sum += pfractionbuf[ihour];
		ofsts << pfractionbuf[ihour] << std::endl;
		}
		ofsts.close();
		printf("%f\n",sum);*/

		double sum = 0;

		for (size_t ihour = 0; ihour < numofhours; ihour++)
		{
			sum += pfractionbuf[ihour];
		}
		if ((sum - 1) > 0.0001)
			printf("%f\n", sum);
		timestructMap[tsid] = pfractionbuf;
		pidbuf += 20;
		pfractionbuf += numofhours;

	}
	//ofs.close();

}

void Reallocator::clearSectors()
{
	for (size_t i = 0; i < sectors.size(); i++)
	{
		delete sectors[i];
	}
	sectors.clear();
}
Reallocator::~Reallocator()
{
	timestructMap.clear();
	if (idbuf)
	{
		delete[] idbuf;
	}

	if (fractionbuf)
	{
		delete[] fractionbuf;
	}
}

