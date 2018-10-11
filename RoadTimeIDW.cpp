#include "RoadTimeIDW.h"
#include "Utils.h"
#include "Utils.h"
#include <qdir.h>
#include <sstream>
RoadTimeIDW::RoadTimeIDW()
{
	m_maxnumpoints = 10;
	m_minnumpoints = 0;
	m_maxDist = 0.0083333333 * 10;
}


RoadTimeIDW::~RoadTimeIDW()
{
	for (size_t i = 0; i < m_stations.size(); i++)
	{
		delete m_stations[i];
	}
}


void RoadTimeIDW::getPoints(OGRGeometry* geom, std::vector<OGRPoint>& points)
{
	OGRwkbGeometryType gtype = geom->getGeometryType();
	if (gtype == wkbLineString || gtype == wkbLineString25D)
	{
		OGRLineString* linestr = (OGRLineString*)geom;
		int len = linestr->getNumPoints();
		for (size_t i = 0; i < len; i++)
		{
			OGRPoint p;
			linestr->getPoint(i, &p);
			points.push_back(p);

		}
	}
	else
	{
		OGRMultiLineString * multilinestr = (OGRMultiLineString*)geom;
		for (size_t i = 0; i < multilinestr->getNumGeometries(); i++)
		{
			getPoints(multilinestr->getGeometryRef(i), points);
		}
	}

}

void RoadTimeIDW::loadStations(std::string stationdir, std::string outputdir)
{

	QDir qdir(outputdir.data());
	if (!qdir.exists())
		qdir.mkpath(".");
	m_outdir = (qdir.absolutePath() + "/").toLocal8Bit().data();
	std::vector<std::string> stationfiles = Utils::findFiles(stationdir, ".data");
	for (size_t i = 0; i < m_stations.size(); i++)
	{
		delete m_stations[i];
	}
	m_stations.clear();
	std::vector<OGRPoint> points;
	for (size_t i = 0; i < stationfiles.size(); i++)
	{
		std::ifstream ifs(stationfiles[i].data(), std::ios::binary);
		int numstations = 0;
		ifs.read(((char*)&numstations), sizeof(int));
		for (size_t j = 0; j < numstations; j++)
		{
			PeMSStation* pems = new PeMSStation;
			pems->fromStream(ifs);
			m_stations.push_back(pems);
		}

	}

}



int loadTimeStruct(std::string filename,double*& fractions)
{
	std::string line;
	std::ifstream infile(filename.data());
	std::vector<std::string> lines;
	while (std::getline(infile, line))
	{
		lines.push_back(line);
	}
	infile.close();
	fractions = new double[lines.size()];
	for (size_t i = 0; i < lines.size(); i++)
	{
		fractions[i] = atof(lines[i].data());
	}
	return lines.size();
}

void RoadTimeIDW::loadStations(std::string shapefile, std::string timestructField,std::string timestructDir)
{

	ShapeFile shp(shapefile);
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	int id = 0;
	int tsfield = shp.poLayer->GetLayerDefn()->GetFieldIndex(timestructField.data());
	std::map<int, OGRPoint> pointmaps;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		int tsid = fea->GetFieldAsInteger(tsfield);
		std::stringstream ssfile;
		ssfile << timestructDir << tsid << ".txt";
		PeMSStation* pems = new PeMSStation;
		pems->m_ID = id;
		pems->m_NumHours = loadTimeStruct(ssfile.str(), pems->m_Fractions);
		pems->m_X = pt->getX();
		pems->m_Y = pt->getY();
		pems->m_Header = "";
		pems->m_HestiaType = "";
		m_stations.push_back(pems);
		id++;
		OGRFeature::DestroyFeature(fea);
	}
	printf("num of points:%d\n", id);

}

void RoadTimeIDW::loadCCSStations(std::string shapefile, std::string temporalFile)
{


	std::ifstream ifs;
	ifs.open(temporalFile, std::ios::binary);
	ifs.seekg(0, ifs.end);
	size_t fileSize = ifs.tellg();
	ifs.seekg(0, ifs.beg);
	std::vector<std::vector<double>> fractions;
	size_t numstations = fileSize / (sizeof(double) * 8760);
	printf("%d,%d", fileSize,numstations);
	for (size_t i = 0; i < numstations; i++)
	{
		std::vector<double> fracs;
		fracs.resize(8760);
		ifs.read((char*)&fracs[0], sizeof(double) * 8760);
		fractions.push_back(fracs);
	}
	ifs.close();

	ShapeFile shp(shapefile);
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	int id = 0;
	int idfield = shp.poLayer->GetLayerDefn()->GetFieldIndex("TemporalID");
	std::map<int, OGRPoint> pointmaps;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		int pid = fea->GetFieldAsInteger(idfield);
		if (pid  > -1) {
			PeMSStation* pems = new PeMSStation;
			pems->m_ID = pid;
			pems->m_Fractions = new double[8760];
			std::vector<double>& srcFractions = fractions[pems->m_ID - 1];
			memcpy(pems->m_Fractions, &srcFractions[0], 8760 * sizeof(double));
			pems->m_NumHours = 8760;
			//if (pems->m_ID == 5734) {

			//	pems->toTxt("ljm.txt");
			//	pems->toBin("ljm.bin");
			//	srcFractions = fractions[pems->m_ID - 1];
			//}
			pems->m_X = fea->GetFieldAsDouble("X");
			pems->m_Y = fea->GetFieldAsDouble("Y");
			pems->m_Header = "";
			pems->m_HestiaType = "";
			pems->roadclass = fea->GetFieldAsInteger("class");
			m_stations.push_back(pems);

		}

		id++;
		OGRFeature::DestroyFeature(fea);
	}
	printf("num of points:%d\n", id);

}

void RoadTimeIDW::createIDW()
{
	std::vector<OGRPoint> points;
	for (size_t i = 0; i < m_stations.size(); i++)
	{
		PeMSStation* pems = m_stations[i];
		OGRPoint newp;
		newp.setX(pems->m_X);
		newp.setY(pems->m_Y);
		points.push_back(newp);
	}
	m_searchEngine.Create(points, m_maxnumpoints);
}

void RoadTimeIDW::stations2shapes(std::string shapefile)
{
	ShapeFile shp;
	shp.create(shapefile, 0, 0, OGRwkbGeometryType::wkbPoint);
	int stationidfield = shp.getOrCreateField("StationID", OGRFieldType::OFTInteger);
	for (size_t i = 0; i < m_stations.size(); i++)
	{
		OGRFeature* fea = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		OGRPoint pt;
		pt.setX(m_stations[i]->m_X);
		pt.setY(m_stations[i]->m_Y);
		fea->SetField(stationidfield, m_stations[i]->m_ID);
		fea->SetGeometry(&pt);
		shp.poLayer->CreateFeature(fea);
		OGRFeature::DestroyFeature(fea);
	}
}

void RoadTimeIDW::coordinatesFromShapefile(std::string shapefile,std::string idField)
{
	ShapeFile shp(shapefile);
	OGRFeature *fea;
	shp.poLayer->ResetReading();
	int id = 0;
	int stationidfield = shp.poLayer->GetLayerDefn()->GetFieldIndex(idField.data());
	std::map<int, OGRPoint> pointmaps;
	while ((fea = shp.poLayer->GetNextFeature()) != NULL)
	{

		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();
		if (stationidfield > -1)
		{
			int stationid = fea->GetFieldAsInteger(stationidfield);
			pointmaps[stationidfield] = *pt;
		}
		else
		{
			m_stations[id]->m_X = pt->getX();
			m_stations[id]->m_Y = pt->getY();
		}
		id++;
		OGRFeature::DestroyFeature(fea);
	}
	printf("num of points:%d\n", id);
}
void RoadTimeIDW::idw(std::string roadshapefile, std::string fieldname, std::string profilePrefixSector, int idoffset)
{


	ShapeFile roadshapes(roadshapefile, 1);
	OGRFeature *poFeature;
	roadshapes.poLayer->ResetReading();
	double maxdisx = m_maxDist * m_maxDist;
	ShapeFile shp;
	shp.create("e:/idw_points.shp", roadshapes.poLayer->GetSpatialRef(), 0, OGRwkbGeometryType::wkbPoint);
	int npointField = shp.getOrCreateField("NPoints", OGRFieldType::OFTInteger);
	int id = 0;

	int timestructFieldNew = roadshapes.getOrCreateField(fieldname.data(), OGRFieldType::OFTString);

	while ((poFeature = roadshapes.poLayer->GetNextFeature()) != NULL)
	{
		std::vector<OGRPoint> pointOnSegment;
		getPoints(poFeature->GetGeometryRef(), pointOnSegment);
		int centerIdx = pointOnSegment.size() / 2;
		OGRPoint center = pointOnSegment[centerIdx];
		int numresults = m_searchEngine.Select_Nearest_Points(center);
		std::vector<int> selectedPoints;
		std::vector<double> selectedPointsDist;


		OGRFeature* fea = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		OGRPoint pt;
		pt.setX(center.getX());
		pt.setY(center.getY());

		fea->SetGeometry(&pt);

		int nearestStationID = -1;
		for (size_t i = 0; i < numresults; i++)
		{
			double dist = 0;
			int selectedIdx = -1;
			m_searchEngine.Get_Selected_Point(i, selectedIdx, dist);
			if (dist > maxdisx)
				continue;
			selectedPoints.push_back(selectedIdx);
			selectedPointsDist.push_back(dist);
			//if (dist < m_maxDist / 10000)
			//{
			//	//nearestStationID = selectedIdx;
			//	break;
			//}

		}

		numresults = selectedPoints.size();
		fea->SetField(npointField, numresults);
		int profileID = id + idoffset;
		std::stringstream ss;
		std::stringstream ssProfileID;

		if (numresults == 1 /*|| nearestStationID > 0*/)
		{
			std::map<int, int>::iterator iter = m_stationProfileIDs.find(selectedPoints[0]);
			if (iter == m_stationProfileIDs.end())
			{
				ss << m_outdir << profilePrefixSector << profileID << ".txt";
				m_stations[selectedPoints[0]]->toTxt(ss.str().data());
				m_stationProfileIDs[selectedPoints[0]] = profileID;
			}
			else
			{
				profileID = m_stationProfileIDs[selectedPoints[0]];
			}
			ssProfileID << profilePrefixSector << profileID;
		}
		else if (numresults > 1)
		{
			//if (!QFileInfo(ss.str().data()).exists())
			//{
			ssProfileID << profilePrefixSector << profileID;
			ss << m_outdir << ssProfileID.str().data() << ".txt";
			PeMSStation idwStation;
			idwStation.m_NumHours = m_stations[0]->m_NumHours;
			idwStation.m_Fractions = new double[idwStation.m_NumHours];
			double* fractions = idwStation.m_Fractions;
			for (size_t ihour = 0; ihour < idwStation.m_NumHours; ihour++)
			{
				//反距离权重插值
				double sum1 = 0;
				double sum2 = 0;
				for (int ipoint = 0; ipoint < selectedPoints.size(); ipoint++)
				{

					int sampleID = selectedPoints[ipoint];
					double dist = selectedPointsDist[ipoint];
					double z = m_stations[sampleID]->m_Fractions[ihour];
					dist = 1 / dist;
					sum1 += dist * z;//权重*属性值
					sum2 += dist;//权重
				}
				fractions[ihour] = sum1 / sum2;
			}

			idwStation.toTxt(ss.str().data());
			//}

		}
		else
		{
			ssProfileID << profilePrefixSector << profileID;

		}


		poFeature->SetField(timestructFieldNew, ssProfileID.str().data());
		roadshapes.poLayer->SetFeature(poFeature);
		shp.poLayer->CreateFeature(fea);
		OGRFeature::DestroyFeature(fea);
		id++;
		if (id % 100 == 0)
			printf("segment:%d\n", id);
		OGRFeature::DestroyFeature(poFeature);
	}
	roadshapes.close();


}
void RoadTimeIDW::idw(std::string roadshapefile, std::string fieldnameold, std::string fieldnamenew, std::string profilePrefixSector, std::string profilePrefixFip,int idoffset)
{
	

	ShapeFile roadshapes(roadshapefile, 1);
	OGRFeature *poFeature;
	roadshapes.poLayer->ResetReading();
	double maxdisx = m_maxDist * m_maxDist;
	ShapeFile shp;
	shp.create("e:/idw_points.shp", roadshapes.poLayer->GetSpatialRef(), 0, OGRwkbGeometryType::wkbPoint);
	int npointField = shp.getOrCreateField("NPoints", OGRFieldType::OFTInteger);
	int id = 0;
	int timestructField = roadshapes.getOrCreateField(fieldnameold.data(),OGRFieldType::OFTInteger);
	int timestructFieldNew = roadshapes.getOrCreateField(fieldnamenew.data(), OGRFieldType::OFTString);

	while ((poFeature = roadshapes.poLayer->GetNextFeature()) != NULL)
	{
		std::vector<OGRPoint> pointOnSegment;
		getPoints(poFeature->GetGeometryRef(), pointOnSegment);
		int centerIdx = pointOnSegment.size() / 2;
		OGRPoint center = pointOnSegment[centerIdx];
		int numresults = m_searchEngine.Select_Nearest_Points(center);
		std::vector<int> selectedPoints;
		std::vector<double> selectedPointsDist;


		OGRFeature* fea = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		OGRPoint pt;
		pt.setX(center.getX());
		pt.setY(center.getY());

		fea->SetGeometry(&pt);

		int nearestStationID = -1;
		for (size_t i = 0; i < numresults; i++)
		{
			double dist = 0;
			int selectedIdx = -1;
			m_searchEngine.Get_Selected_Point(i, selectedIdx, dist);
			if (dist > maxdisx)
				continue;
			selectedPoints.push_back(selectedIdx);
			selectedPointsDist.push_back(dist);
			//if (dist < m_maxDist / 10000)
			//{
			//	//nearestStationID = selectedIdx;
			//	break;
			//}

		}

		numresults = selectedPoints.size();
		fea->SetField(npointField, numresults);
		int profileID = id + idoffset;
		std::stringstream ss;
		std::stringstream ssProfileID;

		if (numresults == 1 /*|| nearestStationID > 0*/)
		{
			std::map<int, int>::iterator iter = m_stationProfileIDs.find(selectedPoints[0]);
			if (iter == m_stationProfileIDs.end())
			{
				ss << m_outdir << profilePrefixSector << profilePrefixFip << profileID << ".txt";
				m_stations[selectedPoints[0]]->toTxt(ss.str().data());
				m_stationProfileIDs[selectedPoints[0]] = profileID;
			}
			else
			{
				profileID = m_stationProfileIDs[selectedPoints[0]];
			}
			ssProfileID << profilePrefixSector << profilePrefixFip<< profileID;
		}
		else if (numresults > 1)
		{
			//if (!QFileInfo(ss.str().data()).exists())
			//{
			ssProfileID << profilePrefixSector << profilePrefixFip << profileID;
			ss << m_outdir << ssProfileID .str().data() << ".txt";
			PeMSStation idwStation;
			idwStation.m_NumHours = m_stations[0]->m_NumHours;
			idwStation.m_Fractions = new double[idwStation.m_NumHours];
			double* fractions = idwStation.m_Fractions;
			for (size_t ihour = 0; ihour < idwStation.m_NumHours; ihour++)
			{
				//反距离权重插值
				double sum1 = 0;
				double sum2 = 0;
				for (int ipoint = 0; ipoint < selectedPoints.size(); ipoint++)
				{

					int sampleID = selectedPoints[ipoint];
					double dist = selectedPointsDist[ipoint];
					double z = m_stations[sampleID]->m_Fractions[ihour];
					dist = 1 / dist;
					sum1 += dist * z;//权重*属性值
					sum2 += dist;//权重
				}
				fractions[ihour] = sum1 / sum2;
			}

			idwStation.toTxt(ss.str().data());
			//}

		}
		else
		{
			std::stringstream ssOri;
			int oriID = poFeature->GetFieldAsInteger(timestructField);
			profileID = oriID;// -6000;
			//if (profileID > 9)
			//{
				//printf("error\n");
			//}
			ssProfileID << profilePrefixSector << profileID;
		
		}


		poFeature->SetField(timestructFieldNew, ssProfileID.str().data());
		roadshapes.poLayer->SetFeature(poFeature);
		shp.poLayer->CreateFeature(fea);
		OGRFeature::DestroyFeature(fea);
		id++;
		if (id % 100 == 0)
			printf("segment:%d\n", id);
		OGRFeature::DestroyFeature(poFeature);
	}
	roadshapes.close();


}
void RoadTimeIDW::nearest(std::string roadshapefile, std::string fieldname, std::string profilePrefixSector, std::string profilePrefixFip, int idoffset)
{
	std::map<int, PeMSStation*> stdmap;
	for (size_t i = 0; i < m_stations.size(); i++)
	{
		PeMSStation* idwStation = m_stations[i];
		std::stringstream ss;
		ss << m_outdir << profilePrefixSector << profilePrefixFip << idwStation->m_ID << ".txt";
		if (!QFile::exists(ss.str().data()))
		{
			idwStation->toTxt(ss.str().data());
		}
		std::map<int, PeMSStation*>::iterator iter = stdmap.find(idwStation->m_ID);
		if (iter != stdmap.end())
		{
			printf("%d,%d\n", idwStation->m_ID, iter->first);
			printf("(%f,%f),(%f,%f)\n", idwStation->m_X, idwStation->m_Y, iter->second->m_X, iter->second->m_Y);
		}
		stdmap[idwStation->m_ID] = idwStation;
	}

	ShapeFile roadshapes(roadshapefile, 1);
	OGRFeature *poFeature;
	roadshapes.poLayer->ResetReading();
	double maxdisx = m_maxDist * m_maxDist;
	ShapeFile shp;
	shp.create("e:/idw_points.shp", roadshapes.poLayer->GetSpatialRef(), 0, OGRwkbGeometryType::wkbPoint);
	int npointField = shp.getOrCreateField("StationID", OGRFieldType::OFTInteger);

	int id = 0;
	int timestructField = roadshapes.getOrCreateField(fieldname.data(), OGRFieldType::OFTInteger);

	while ((poFeature = roadshapes.poLayer->GetNextFeature()) != NULL)
	{
		std::vector<OGRPoint> pointOnSegment;
		getPoints(poFeature->GetGeometryRef(), pointOnSegment);
		int centerIdx = pointOnSegment.size() / 2;
		OGRPoint center = pointOnSegment[centerIdx];
		int numresults = m_searchEngine.Select_Nearest_Points(center);
		int selectedIdx = -1;
		double dist;
		m_searchEngine.Get_Selected_Point(0, selectedIdx, dist);

		OGRFeature* fea = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		OGRPoint pt;
		pt.setX(center.getX());
		pt.setY(center.getY());
		fea->SetGeometry(&pt);
		fea->SetField(npointField, m_stations[selectedIdx]->m_ID);

		std::stringstream ss;
		ss << profilePrefixSector << profilePrefixFip << m_stations[selectedIdx]->m_ID;
		poFeature->SetField(timestructField, ss.str().data());
		roadshapes.poLayer->SetFeature(poFeature);
		shp.poLayer->CreateFeature(fea);
		OGRFeature::DestroyFeature(fea);
		id++;
		if (id % 100 == 0)
			printf("segment:%d\n", id);
		OGRFeature::DestroyFeature(poFeature);
	}
	roadshapes.close();


}


void RoadTimeIDW::idw_ccs(std::string roadshapefile, OGREnvelope bound, double resol)
{
	int ncol = (int)((bound.MaxX - bound.MinX) / resol) + 1;
	ShapeFile roadshapes(roadshapefile, 1);
	OGRFeature *poFeature;
	roadshapes.poLayer->ResetReading();
	double maxdisx = m_maxDist * m_maxDist;
	//ShapeFile shp;
	//shp.create("e:/idw_points.shp", roadshapes.poLayer->GetSpatialRef(), 0, OGRwkbGeometryType::wkbPoint);
	//int npointField = shp.getOrCreateField("NPoints", OGRFieldType::OFTInteger);
	int id = 0;
	int f_time = roadshapes.getOrCreateField("index", OGRFieldType::OFTInteger);
	int f_class = roadshapes.getOrCreateField("class", OGRFieldType::OFTInteger);
	while ((poFeature = roadshapes.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint center = *((OGRPoint*)poFeature->GetGeometryRef());
		int classid = poFeature->GetFieldAsInteger(f_class);
		if (classid <= 0)
			classid = 1;
		int col = (int)((center.getX() - bound.MinX) / resol);
		int row = (int)((bound.MaxY - center.getY()) / resol);
		int gridID = col + row * ncol + 1;
		double cellx = bound.MinX + resol * (col + 0.5);
		double celly = bound.MaxY - resol * (row + 0.5);
		center.setX(cellx);
		center.setY(celly);
		TimeProfileCoords coords;
		int numresults = m_searchEngine.Select_Nearest_Points(center);
		std::vector<int> selectedPoints;
		std::vector<double> selectedPointsDist;
		//OGRFeature* fea = OGRFeature::CreateFeature(shp.poLayer->GetLayerDefn());
		//OGRPoint pt;
		//pt.setX(center.getX());
		//pt.setY(center.getY());
		//fea->SetGeometry(&pt);
		int nearestStationID = -1;
		for (size_t i = 0; i < numresults; i++)
		{
			double dist = 0;
			int selectedIdx = -1;
			m_searchEngine.Get_Selected_Point(i, selectedIdx, dist);
			PeMSStation* sel = m_stations[selectedIdx];
			if (classid != sel->roadclass)
				continue;
			dist = sqrt((sel->m_X - cellx) * (sel->m_X - cellx) + (sel->m_Y - celly) * (sel->m_Y - celly));
			dist = dist*dist;
			if (dist > maxdisx)
				continue;
			selectedPoints.push_back(selectedIdx);
			selectedPointsDist.push_back(dist);
		}

		numresults = selectedPoints.size();
		//fea->SetField(npointField, numresults);

		std::stringstream ss;
		std::stringstream ssProfileID;

		bool fileexist = true;
		if (numresults == 1)
		{
			int sampleID = selectedPoints[0];
			gridID = m_stations[sampleID]->m_ID;
			ssProfileID << 5 << gridID;
			ss << m_outdir << ssProfileID.str().data() << ".bin";
			//printf("%d,%d\n", sampleID, sampleID);
			if (!QFileInfo(ss.str().data()).exists()) {
				m_stations[sampleID]->toBin(ss.str().data());
				fileexist = false;
			}
			coords.X = m_stations[sampleID]->m_X;
			coords.Y = m_stations[sampleID]->m_Y;
			m_TimeProfileCoords[atoi(ssProfileID.str().data())] = coords;
		}
		else if (numresults > 1)
		{
			ssProfileID << classid << gridID;
			ss << m_outdir << ssProfileID.str().data() << ".bin";
			if (!QFileInfo(ss.str().data()).exists()) {
				fileexist = false;
				PeMSStation idwStation;
				idwStation.m_NumHours = m_stations[0]->m_NumHours;
				idwStation.m_Fractions = new double[idwStation.m_NumHours];
				double* fractions = idwStation.m_Fractions;
				for (size_t ihour = 0; ihour < idwStation.m_NumHours; ihour++)
				{
					//反距离权重插值
					double sum1 = 0;
					double sum2 = 0;
					int maxnumpoints = 10;
					if (selectedPoints.size() < 10)
						maxnumpoints = selectedPoints.size();
					for (int ipoint = 0; ipoint < maxnumpoints; ipoint++)
					{

						int sampleID = selectedPoints[ipoint];
						double dist = selectedPointsDist[ipoint];
						double z = m_stations[sampleID]->m_Fractions[ihour];
						dist = 1 / dist;
						sum1 += dist * z;//权重*属性值
						sum2 += dist;//权重
					}
					fractions[ihour] = sum1 / sum2;
				}
				idwStation.toBin(ss.str().data());
			}
			coords.X = cellx;
			coords.Y = celly;
			m_TimeProfileCoords[atoi(ssProfileID.str().data())] = coords;
		}
		else
		{
			ssProfileID << "-1";
			fileexist = false;

		}
		if (!fileexist)
		{
			poFeature->SetField(f_time, atoi(ssProfileID.str().data()));
			roadshapes.poLayer->SetFeature(poFeature);
		}


		//shp.poLayer->CreateFeature(fea);
		//OGRFeature::DestroyFeature(fea);
		id++;
		if (id % 100 == 0)
			printf("segment:%d\n", id);
		OGRFeature::DestroyFeature(poFeature);
	}
	roadshapes.close();


}
void RoadTimeIDW::fillgap_ccsp(std::string roadshapefile, std::string samplefile)
{
	
	std::vector<OGRPoint> points;
	std::vector<int> ids;

	std::ifstream ifs;
	ifs.open("H:/Vulcan_2014/support_data/onroad/jianming/TimeProfileCoords.csv");
	std::string header;
	std::getline(ifs, header);

	std::string line;
	while (std::getline(ifs, line))
	{
		std::vector<std::string > splits = Utils::splitCSV(',', line);
		ids.push_back(atoi(splits[0].data()));
		OGRPoint pt;
		pt.setX(atof(splits[1].data()));
		pt.setY(atof(splits[2].data()));
		points.push_back(pt);
	}
	ifs.close();
	ANNSearchEngine searchEngine;
	searchEngine.Create(points, m_maxnumpoints);

	ShapeFile roadshapes(roadshapefile, 1);
	OGRFeature *poFeature;
	roadshapes.poLayer->ResetReading();
	double maxdisx = m_maxDist * m_maxDist;

	int id = 0;
	int f_time = roadshapes.getOrCreateField("index", OGRFieldType::OFTInteger);
	int f_class = roadshapes.getOrCreateField("class", OGRFieldType::OFTInteger);
	while ((poFeature = roadshapes.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint center = *((OGRPoint*)poFeature->GetGeometryRef());
		int classid = poFeature->GetFieldAsInteger(f_class);
		int tid = poFeature->GetFieldAsInteger(f_time);

		if (tid < 0) {
			TimeProfileCoords coords;
			int numresults = searchEngine.Select_Nearest_Points(center);
			std::vector<int> selectedPoints;
			std::vector<double> selectedPointsDist;
			double dist = 0;
			int selectedIdx = -1;
			searchEngine.Get_Selected_Point(0, selectedIdx, dist);
			tid = ids[selectedIdx];
			poFeature->SetField(f_time, tid);
			roadshapes.poLayer->SetFeature(poFeature);
		}
		
		id++;
		if (id % 100 == 0)
			printf("segment:%d\n", id);
		OGRFeature::DestroyFeature(poFeature);
	}
	roadshapes.close();


}