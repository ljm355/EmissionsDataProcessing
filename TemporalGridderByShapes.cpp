#include "TemporalGridderByShapes.h"


void FFCO2Sector::loadAttribute(std::string _fieldname)
{
	fieldname = _fieldname;
	features.clear();
	for (int ishapefile = 0; ishapefile < shapefiles.size(); ishapefile++)
	{
		std::string shapefile = shapefiles[ishapefile];
		ShapeFile input(shapefile.data());
		OGRFeature *poFeature;
		input.poLayer->ResetReading();

		int caIndex = input.poLayer->GetLayerDefn()->GetFieldIndex(fieldname.data());
		int idIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("Id");

		int timestructIndex = input.poLayer->GetLayerDefn()->GetFieldIndex("timestruct");
		bool isfloat = false;
		if (input.poLayer->GetLayerDefn()->GetFieldDefn(timestructIndex)->GetType() == OGRFieldType::OFTReal)
			isfloat = true;
		OGRwkbGeometryType gtype = input.poLayer->GetGeomType();
		while ((poFeature = input.poLayer->GetNextFeature()) != NULL)
		{
			double ca = poFeature->GetFieldAsDouble(caIndex);
			if (ca > 0)
			{
				int id = poFeature->GetFieldAsInteger(idIndex);
				std::string tsid;
				if (!isfloat)
				{
					tsid = poFeature->GetFieldAsString(timestructIndex);
				}
				else
				{
					std::stringstream ss;
					ss << (long)poFeature->GetFieldAsDouble(timestructIndex);
					tsid = ss.str();

				}


				FeatureCell fc;
				fc.total = ca;
				fc.gridid = id;
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
				features.push_back(fc);
			}
			OGRFeature::DestroyFeature(poFeature);
		}
	}
}
void FFCO2Sector::getTimeSlice(double* gridcells, int houridx)
{
	//memset(gridcells, 0, sizeof(double)*slice);
	for (size_t i = 0; i < features.size(); i++)
	{
		FeatureCell& fc = features[i];
		if (fc.total <= 0)
			continue;
		//printf("%d\n", fc.gridid);
		gridcells[fc.gridid] += (fc.total * fc.timestructure[houridx]);
	}
}
void FFCO2Sector::getTotal(double* gridcells)
{
	memset(gridcells, 0, sizeof(double)*numcells);
	for (size_t i = 0; i < features.size(); i++)
	{
		FeatureCell& fc = features[i];
		gridcells[fc.gridid] += fc.total;
	}
}
void FFCO2Sector::getTimeSlice(int houridx)
{
	memset(cells, 0, sizeof(double)*numcells);
	for (size_t i = 0; i < features.size(); i++)
	{
		FeatureCell& fc = features[i];
		if (fc.total <= 0)
			continue;
		//printf("%d\n", fc.gridid);
		cells[fc.gridid] += (fc.total * fc.timestructure[houridx]);
	}
}
void FFCO2Sector::getTotal()
{
	memset(cells, 0, sizeof(double)*numcells);
	for (size_t i = 0; i < features.size(); i++)
	{
		FeatureCell& fc = features[i];
		cells[fc.gridid] += fc.total;
	}
}
void FFCO2Sector::addCells(double* grid1, double* grid2, int numcells)
{
	for (size_t i = 0; i < numcells; i++)
	{
		grid2[i] = grid1[i] + grid2[i];
	}
}
void FFCO2Sector::makeHourlyTotal(std::string hourlyFilename)
{
	if (timestructMap.size() == 0)
	{
		printf("Time structures not loaded!");
		return;
	}

	//std::ofstream ofs;
	//ofs.open(hourlyFilename.data());
	std::ofstream ofs(hourlyFilename.data(), std::ios::binary);
	double* sectorCells = new double[numcells];
	double sum = 0;
	for (int ihour = 0; ihour < numhours; ihour++)
	{
		memset(cells, 0, numcells*sizeof(double));
		getTimeSlice(cells, ihour);
		//for (int ncell = 0; ncell < numcells; ncell++)
		//{
		//	ofs << cells[ncell] << std::endl;
		//}
		ofs.write((char*)&cells[0], (size_t)numcells * sizeof(double));
		for (size_t icell = 0; icell < numcells; icell++)
		{
			sum += cells[icell];
		}
		if (ihour % 100 == 0)
		{
			printf("%d/%d\n", ihour, numhours);
		}
	}
	printf("total=%f\n", sum);
	memset(cells, 0, numcells*sizeof(double));
	ofs.close();

}
void FFCO2Sector::makeDailyTotal(std::string dailyFilename)
{
	if (timestructMap.size() == 0)
	{
		printf("Time structures not loaded!");
		return;
	}
	//std::ofstream ofs;
	//ofs.open(dailyFilename.data());
	std::ofstream ofs(dailyFilename.data(), std::ios::binary);
	double* hourlyCells = new double[numcells];
	double sum = 0;
	int ihour = 0;
	for (int iday = 0; iday < numdays; iday++)
	{
		memset(cells, 0, numcells*sizeof(double));
		for (size_t hourofday = 0; hourofday < 24; hourofday++)
		{
			memset(hourlyCells, 0, numcells*sizeof(double));
			getTimeSlice(hourlyCells, ihour);
			addCells(hourlyCells, cells, numcells);
			ihour++;
		}
		//for (int ncell = 0; ncell < numcells; ncell++)
		//{
		//	ofs << cells[ncell] << std::endl;
		//}	
		ofs.write((char*)&cells[0], (size_t)numcells * sizeof(double));
		for (size_t icell = 0; icell < numcells; icell++)
		{
			sum += cells[icell];
		}
		printf("%d/%d\n", ihour, numhours);
	}

	delete[] hourlyCells;
	memset(cells, 0, numcells*sizeof(double));
	ofs.close();

}
void FFCO2Sector::makeAnnualTotal(std::string totalFilename)
{
	//std::ofstream ofs;
	//ofs.open(totalFilename.data());
	std::ofstream ofs(totalFilename.data(), std::ios::binary);
	memset(cells, 0, numcells*sizeof(double));
	getTotal(cells);
	//for (int ncell = 0; ncell < numcells; ncell++)
	//{
	//	ofs << cells[ncell] << std::endl;
	//}
	ofs.write((char*)&cells[0], (size_t)numcells * sizeof(double));
	memset(cells, 0, numcells*sizeof(double));
	ofs.close();
}
TemporalGridderByShapes::~TemporalGridderByShapes()
{
	timestructMap.clear();
	if (idbuf)
	{
		delete[] idbuf;
		delete[] fractionbuf;
	}
	clearSectors();
	delete cells;
}
void TemporalGridderByShapes::normalizeFractions(double*& fractions)
{
	double sum = 0;
	for (size_t i = 0; i < numhours; i++)
	{
		sum += fractions[i];
	}
	for (size_t i = 0; i < numhours; i++)
	{
		fractions[i] = fractions[i] / sum;
	}
}
void TemporalGridderByShapes::clearSectors()
{
	for (size_t i = 0; i < sectors.size(); i++)
	{
		delete sectors[i];
	}
	sectors.clear();
}
void TemporalGridderByShapes::loadtimestruct(std::string txtdir, std::string binaryfile)
{

	if (!QFileInfo(binaryfile.data()).exists())
	{
		TimestructTool::txt2binary(numhours, txtdir, binaryfile, true);
	}
	loadtimestruct(binaryfile);

}
void TemporalGridderByShapes::loadtimestruct(std::string binaryfile)
{
	TimeStructContainer* ts = new TimeStructContainer;
	ts->load(binaryfile, this->timestructMap);
	timestructs.push_back(ts);
}

void TemporalGridderByShapes::clearTimeStruct()
{
	timestructMap.clear();
	for (int i = 0; i < timestructs.size(); i++)
	{
		delete timestructs[i];
	}
	timestructs.clear();
}
void TemporalGridderByShapes::addSectorGrid(std::vector<std::string> shapefiles, std::string sectorname)
{
	FFCO2Sector* sectorGrid = new FFCO2Sector(numhours);
	sectorGrid->timestructMap = timestructMap;
	sectorGrid->shapefiles = shapefiles;
	sectorGrid->sectorname = sectorname;
	sectorGrid->numcells = numcells;
	sectorGrid->numhours = numhours;
	sectorGrid->numdays = numdays;
	sectorGrid->cells = new double[numcells];
	memset(sectorGrid->cells, 0, numcells*sizeof(double));
	sectors.push_back(sectorGrid);
}
void TemporalGridderByShapes::loadAttribute(std::string fieldname)
{
	this->fieldname = fieldname;
	for (int isector = 0; isector < sectors.size(); isector++)
	{
		FFCO2Sector* sectorGrid = sectors[isector];
		sectorGrid->timestructMap = this->timestructMap;
		sectorGrid->loadAttribute(fieldname);
	}
}

void TemporalGridderByShapes::getTotal()
{
	
	double* sectorCells = new double[numcells];
	memset(cells, 0, numcells*sizeof(double));
	for (int isector = 0; isector < sectors.size(); isector++)
	{
		memset(sectorCells, 0, numcells*sizeof(double));
		FFCO2Sector* sectorGrid = sectors[isector];
		sectorGrid->getTotal();
		addCells(sectorGrid->cells, cells, numcells);
	}
	delete[] sectorCells;
}
void TemporalGridderByShapes::makeHourlyTotal(std::string hourlyFilename)
{
	if (timestructMap.size() == 0)
	{
		printf("Time structures not loaded!");
		return;
	}

	std::ofstream ofs(hourlyFilename.data(), std::ios::binary);
	double* sectorCells = new double[numcells];
	double sum = 0;
	for (int ihour = 0; ihour < numhours; ihour++)
	{
		memset(cells, 0, numcells*sizeof(double));
		for (int isector = 0; isector < sectors.size(); isector++)
		{
			memset(sectorCells, 0, numcells*sizeof(double));
			FFCO2Sector* sectorGrid = sectors[isector];
			sectorGrid->getTimeSlice(sectorCells, ihour);
			addCells(sectorCells, cells, numcells);
			//sectorGrid->hourlyNCFile.writenumcells(ihour, sectorCells);
		}
		//for (int ncell = 0; ncell < numcells; ncell++)
		//{
		//	ofs << cells[ncell] << std::endl;
		//}
		ofs.write((char*)&cells[0], (size_t)numcells * sizeof(double));
		for (size_t icell = 0; icell < numcells; icell++)
		{
			sum += cells[icell];
		}
		if (ihour % 100 == 0)
		{
			printf("%d/%d\n", ihour, numhours);
		}
	}
	printf("total=%f\n", sum);

	for (int isector = 0; isector < sectors.size(); isector++)
	{
		FFCO2Sector* sectorGrid = sectors[isector];
		memset(sectorGrid->cells, 0, numcells*sizeof(double));
	}
	memset(cells, 0, numcells*sizeof(double));
	ofs.close();

}
void TemporalGridderByShapes::makeDailyTotal(std::string dailyFilename)
{
	if (timestructMap.size() == 0)
	{
		printf("Time structures not loaded!");
		return;
	}
	std::ofstream ofs(dailyFilename.data(), std::ios::binary);
	//ofs.open(dailyFilename.data());
	double* sectorCells = new double[numcells];
	double sum = 0;
	int ihour = 0;
	for (int iday = 0; iday < numdays; iday++)
	{
		memset(cells, 0, numcells*sizeof(double));
		for (size_t hourofday = 0; hourofday < 24; hourofday++)
		{
			for (int isector = 0; isector < sectors.size(); isector++)
			{
				memset(sectorCells, 0, numcells*sizeof(double));
				FFCO2Sector* sectorGrid = sectors[isector];
				sectorGrid->getTimeSlice(sectorCells, ihour);
				addCells(sectorCells, cells, numcells);
			}
			ihour++;
		}
		//for (int ncell = 0; ncell < numcells; ncell++)
		//{
		//	ofs << cells[ncell] << std::endl;
		//}
		ofs.write((char*)&cells[0], (size_t)numcells * sizeof(double));
		for (size_t icell = 0; icell < numcells; icell++)
		{
			sum += cells[icell];
		}
		printf("%d/%d\n", ihour, numhours);
	}

	for (int isector = 0; isector < sectors.size(); isector++)
	{
		FFCO2Sector* sectorGrid = sectors[isector];
		memset(sectorGrid->cells, 0, numcells*sizeof(double));
	}
	memset(cells, 0, numcells*sizeof(double));
	ofs.close();

}
void TemporalGridderByShapes::makeAnnualTotal(std::string totalFilename)
{
	//std::ofstream ofs;
	//ofs.open(totalFilename.data());
	std::ofstream ofs(totalFilename.data(), std::ios::binary);
	double* sectorCells = new double[numcells];

	memset(cells, 0, numcells*sizeof(double));
	for (int isector = 0; isector < sectors.size(); isector++)
	{
		memset(sectorCells, 0, numcells*sizeof(double));
		FFCO2Sector* sectorGrid = sectors[isector];
		sectorGrid->getTotal();
		addCells(sectorGrid->cells, cells, numcells);
	}
	//for (int ncell = 0; ncell < numcells; ncell++)
	//{
	//	ofs << cells[ncell] << std::endl;
	//}
	ofs.write((char*)&cells[0], (size_t)numcells * sizeof(double));
	ofs.close();
}