#include "Vulcan2014.h"



Vulcan2014::Vulcan2014()
{
	m_inputdir = "Z:/Vulcan/Vulcan_2014/output_data/";
	m_inputdir = "E:/Vulcan/output_data/";
	m_outputdir = "E:/Vulcan/gridPrep_SHP_master/";

}


Vulcan2014::~Vulcan2014()
{
}
void Vulcan2014::createAirport()
{
	ShapeFile points;
	points.create(m_outputdir + "Airport.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
	std::ifstream ifs;
	ifs.open(m_inputdir + "Airport/Airport.em.CO2.2011.csv");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int latidx = -1;
	int lonidx = -1;
	int co2idx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("latitude").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("longitude").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("final.CO2.tC").toLower())
			co2idx = i;
	}
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		poFeature->SetGeometry(&pt);
		std::string co2str = splits[co2idx];
		if (co2str != "NA")
		{
			poFeature->SetField(co2fidx, atof(co2str.data()));
		}
		points.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);

	}
}
void Vulcan2014::createRailroadPoint()
{
	ShapeFile points;
	points.create(m_outputdir + "RailroadPoint.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
	std::ifstream ifs;
	ifs.open(m_inputdir + "Rail/Point.RAIL.em.CO2.2011.smpl.csv");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int latidx = -1;
	int lonidx = -1;
	int co2idx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("latitude").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("longitude").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("CO2 emit (tC)").toLower())
			co2idx = i;
	}
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		poFeature->SetGeometry(&pt);
		std::string co2str = splits[co2idx];
		if (co2str != "NA")
		{
			poFeature->SetField(co2fidx, atof(co2str.data()));
		}
		points.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);

	}
}
void Vulcan2014::createNonroadPoint()
{
	ShapeFile points;
	points.create(m_outputdir + "NonroadPoint.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
	std::ifstream ifs;
	ifs.open(m_inputdir + "Nonroad/Point.NNRD.em.CO2.2011.smpl.csv");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int latidx = -1;
	int lonidx = -1;
	int co2idx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("latitude").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("longitude").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("CO2 emit (tC)").toLower())
			co2idx = i;
	}
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		poFeature->SetGeometry(&pt);
		std::string co2str = splits[co2idx];
		if (co2str != "NA")
		{
			poFeature->SetField(co2fidx, atof(co2str.data()));
		}
		points.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);

	}
}
//void Vulcan2014::createElecProd()
//{
//
//	ShapeFile points;
//	points.create(m_outputdir + "ElecProd.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
//	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
//	int sourcefidx = points.getOrCreateField("source", OGRFieldType::OFTString);
//	std::ifstream ifs;
//	ifs.open(m_inputdir + "Electricity_production/Point.WF.ELEC.em.CO2.2011.csv");
//	std::string line;
//	std::getline(ifs, line);
//	std::vector<std::string> headerline = Utils::splitCSV(',', line);
//	int latidx = -1;
//	int lonidx = -1;
//	int co2idx = -1;
//	for (size_t i = 0; i < headerline.size(); i++)
//	{
//		if (headerline[i] == "latitude")
//			latidx = i;
//		if (headerline[i] == "longitude")
//			lonidx = i;
//		if (headerline[i] == "SR.CO2.tC")
//			co2idx = i;
//	}
//	int idx = 0;
//	while (ifs.peek() != -1)
//	{
//		std::getline(ifs, line);
//		std::vector<std::string> splits = Utils::splitCSV(',', line);
//		idx++;
//		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
//		std::string latstr = splits[latidx];
//		std::string lonstr = splits[lonidx];
//		OGRPoint pt;
//		pt.setX(atof(lonstr.data()));
//		pt.setY(atof(latstr.data()));
//		if (pt.getY() < -6)
//		{
//			printf("");
//		}
//		poFeature->SetGeometry(&pt);
//		std::string co2str = splits[co2idx];
//		if (co2str != "NA")
//		{
//			poFeature->SetField(co2fidx, atof(co2str.data()));
//		}
//		poFeature->SetField(sourcefidx, "NEI");
//		points.poLayer->CreateFeature(poFeature);
//		OGRFeature::DestroyFeature(poFeature);
//	}
//	ifs.close();
//
//	ifs.open(m_inputdir + "Electricity_production/Electricity.production.em.2011.csv");
//	line = "";
//	//std::getline(ifs, line, '\r');
//	std::getline(ifs, line);
//	headerline = Utils::splitCSV(',', line);
//	latidx = -1;
//	lonidx = -1;
//	co2idx = -1;
//	int sourceidx = -1;
//	for (size_t i = 0; i < headerline.size(); i++)
//	{
//		if (headerline[i] == "longitude")
//			lonidx = i;
//		if (headerline[i] == "latitude")
//			latidx = i;
//		if (headerline[i] == "CO2_MASS.tC")
//			co2idx = i;
//		if (headerline[i] == "emit source")
//			sourceidx = i;
//	}
//	idx = 0;
//	while (ifs.peek() != -1)
//	{
//		//std::getline(ifs, line,'\r');
//		std::getline(ifs, line);
//		std::vector<std::string> splits = Utils::splitCSV(',', line);
//		idx++;
//		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
//		std::string latstr = splits[latidx];
//		std::string lonstr = splits[lonidx];
//		std::string sourcestr = splits[sourceidx];
//		if (isdigit(latstr[latstr.size()-1]) && isdigit(lonstr[lonstr.size()-1]))
//		{
//			OGRPoint pt;
//			pt.setX(atof(lonstr.data()));
//			pt.setY(atof(latstr.data()));
//			if (pt.getY() < -6)
//			{
//				printf("");
//			}
//			poFeature->SetGeometry(&pt);
//			std::string co2str = splits[co2idx];
//			if (co2str != "#N/A")
//			{
//				poFeature->SetField(co2fidx, atof(co2str.data()));
//			}
//			else
//			{
//				//printf("%f,%f,%f\n", pt.getX(), pt.getY(), atof(co2str.data()));
//				printf("%f,%f,%s\n", pt.getX(), pt.getY(), "NA");
//			}
//			poFeature->SetField(sourcefidx, sourcestr.data());
//
//			points.poLayer->CreateFeature(poFeature);
//		}
//		else
//		{
//			printf("%s\n", line.data());
//		}
//
//		OGRFeature::DestroyFeature(poFeature);
//	}
//	ifs.close();
//
//
//}
void Vulcan2014::createElecProd()
{

	ShapeFile points;
	points.create(m_outputdir + "ElecProd.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
	int sourcefidx = points.getOrCreateField("source", OGRFieldType::OFTString);
	std::ifstream ifs;
	ifs.open(m_inputdir + "Elec_prod/Point.ELEC.em.CO2.2011.smpl.csv");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int latidx = -1;
	int lonidx = -1;
	int co2idx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("latitude").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("longitude").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("CO2 emit (tC)").toLower())
			co2idx = i;
	}
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		if (pt.getY() < -6)
		{
			printf("");
		}
		poFeature->SetGeometry(&pt);
		std::string co2str = splits[co2idx];
		if (co2str != "NA")
		{
			poFeature->SetField(co2fidx, atof(co2str.data()));
		}
		poFeature->SetField(sourcefidx, "NEI");
		points.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);
	}
	ifs.close();

	ifs.open(m_inputdir + "Elec_prod/Electricity.production.em.nice.2011.csv");
	line = "";
	//std::getline(ifs, line, '\r');
	std::getline(ifs, line);
	headerline = Utils::splitCSV(',', line);
	latidx = -1;
	lonidx = -1;
	co2idx = -1;
	int sourceidx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("latitude").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("longitude").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("Total fossil (tC)").toLower())
			co2idx = i;
		if (QString(headerline[i].data()).toLower() == QString("emit.source").toLower())
			sourceidx = i;
	}
	idx = 0;
	while (ifs.peek() != -1)
	{
		//std::getline(ifs, line,'\r');
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		std::string sourcestr = splits[sourceidx];
		if (isdigit(latstr[latstr.size() - 1]) && isdigit(lonstr[lonstr.size() - 1]))
		{
			OGRPoint pt;
			pt.setX(atof(lonstr.data()));
			pt.setY(atof(latstr.data()));
			if (pt.getY() < -6)
			{
				printf("");
			}
			poFeature->SetGeometry(&pt);
			std::string co2str = splits[co2idx];
			if (co2str != "#N/A")
			{
				poFeature->SetField(co2fidx, atof(co2str.data()));
			}
			else
			{
				//printf("%f,%f,%f\n", pt.getX(), pt.getY(), atof(co2str.data()));
				printf("%f,%f,%s\n", pt.getX(), pt.getY(), "NA");
			}
			poFeature->SetField(sourcefidx, sourcestr.data());

			points.poLayer->CreateFeature(poFeature);
		}
		else
		{
			printf("%s\n", line.data());
		}

		OGRFeature::DestroyFeature(poFeature);
	}
	ifs.close();


}
void Vulcan2014::createIndPoint()
{
	ShapeFile points;
	points.create(m_outputdir + "IndPoint.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
	std::ifstream ifs;
	ifs.open(m_inputdir + "Industrial/Point.INDS.em.CO2.2011.smpl.csv ");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int latidx = -1;
	int lonidx = -1;
	int co2idx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("latitude").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("longitude").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("CO2 emit (tC)").toLower())
			co2idx = i;
	}
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		poFeature->SetGeometry(&pt);
		std::string co2str = splits[co2idx];
		if (co2str != "NA")
		{
			poFeature->SetField(co2fidx, atof(co2str.data()));
		}
		points.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);

	}
}
void Vulcan2014::createComPoint()
{
	ShapeFile points;
	points.create(m_outputdir + "ComPoint.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
	std::ifstream ifs;
	ifs.open(m_inputdir + "Commercial/Point.COMM.em.CO2.2011.smpl.csv");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int latidx = -1;
	int lonidx = -1;
	int co2idx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("latitude").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("longitude").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("CO2 emit (tC)").toLower())
			co2idx = i;
	}
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		poFeature->SetGeometry(&pt);
		std::string co2str = splits[co2idx];
		if (co2str != "NA")
		{
			poFeature->SetField(co2fidx, atof(co2str.data()));
		}
		points.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);

	}
}
void Vulcan2014::createCementPoint()
{
	ShapeFile points;
	points.create(m_outputdir + "Cement.shp", NULL, NULL, OGRwkbGeometryType::wkbPoint);
	int co2fidx = points.getOrCreateField("ca11", OGRFieldType::OFTReal);
	std::ifstream ifs;
	ifs.open(m_inputdir + "Cement/Cement.em.CO2.2011.csv");
	std::string line;
	std::getline(ifs, line);
	std::vector<std::string> headerline = Utils::splitCSV(',', line);
	int latidx = -1;
	int lonidx = -1;
	int co2idx = -1;
	for (size_t i = 0; i < headerline.size(); i++)
	{
		if (QString(headerline[i].data()).toLower() == QString("lat").toLower())
			latidx = i;
		if (QString(headerline[i].data()).toLower() == QString("lon").toLower())
			lonidx = i;
		if (QString(headerline[i].data()).toLower() == QString("CO2.final.tC").toLower())
			co2idx = i;
	}
	int idx = 0;
	while (ifs.peek() != -1)
	{
		std::getline(ifs, line);
		std::vector<std::string> splits = Utils::splitCSV(',', line);
		idx++;
		OGRFeature* poFeature = OGRFeature::CreateFeature(points.poLayer->GetLayerDefn());
		std::string latstr = splits[latidx];
		std::string lonstr = splits[lonidx];
		OGRPoint pt;
		pt.setX(atof(lonstr.data()));
		pt.setY(atof(latstr.data()));
		poFeature->SetGeometry(&pt);
		std::string co2str = splits[co2idx];
		if (co2str != "NA")
		{
			poFeature->SetField(co2fidx, atof(co2str.data()));
		}
		points.poLayer->CreateFeature(poFeature);
		OGRFeature::DestroyFeature(poFeature);

	}
}

void Vulcan2014::sumBySectorRaster(std::string dir, std::vector<std::vector<std::string>> sectorfiles, std::vector<std::string> sectornames, std::string outcsvfile)
{
	std::ofstream ofs;
	ofs.open(outcsvfile);
	double total = 0;
	for (size_t i = 0; i < sectorfiles.size(); i++)
	{
		double subtotal = 0;
		for (size_t j = 0; j < sectorfiles[i].size(); j++)
		{
			GDAL_DS<double>* ds = new GDAL_DS<double>();
			ds->open(dir + sectorfiles[i][j]);
			double filesum = ds->sum(1);
			subtotal += filesum;
			printf("%s,%f\n", sectorfiles[i][j].data(), filesum);
		}
		total += subtotal;
		ofs << sectornames[i] << "," << subtotal << std::endl;
	}
	ofs << "Total" << "," << total << std::endl;
	ofs.close();
}
void Vulcan2014::sumBySectorRaster(std::string dir, std::string outcsvfile)
{
	std::vector<std::vector<std::string>> sectorfiles;
	std::vector<std::string> sectornames;
	std::vector<std::string> files;

	sectornames.push_back("Electricity");
	files.clear();
	files.push_back("ElecProd.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("OnRoad");
	files.clear();
	files.push_back("OnRoadRuralLocal.tif");
	files.push_back("OnRoadRuralNonLocal.tif");
	files.push_back("OnRoadUrbanLocal.tif");
	files.push_back("OnRoadUrbanNonLocal.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("Commercial");
	files.clear();
	files.push_back("ComNonPoint.tif");
	files.push_back("ComPoint.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("Residential");
	files.clear();
	files.push_back("ResNonPoint.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("Industrial");
	files.clear();
	files.push_back("IndNonPoint.tif");
	files.push_back("IndPoint.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("Airport");
	files.clear();
	files.push_back("Airport.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("Cement");
	files.clear();
	files.push_back("Cement.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("CMV");
	files.clear();
	files.push_back("CMVPort.tif");
	files.push_back("CMVShippingLane.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("Railroad");
	files.clear();
	files.push_back("Railroad.tif");
	files.push_back("RailroadPoint.tif");
	sectorfiles.push_back(files);

	sectornames.push_back("Nonroad");
	files.clear();
	files.push_back("NonroadPoint.tif");
	files.push_back("Nonroad_County.tif");
	files.push_back("Nonroad_S100FRAC.tif");
	files.push_back("Nonroad_S140FRAC.tif");
	files.push_back("Nonroad_S260FRAC.tif");
	files.push_back("Nonroad_S300FRAC.tif");
	files.push_back("Nonroad_S310FRAC.tif");
	files.push_back("Nonroad_S311FRAC.tif");
	files.push_back("Nonroad_S350FRAC.tif");
	files.push_back("Nonroad_S400FRAC.tif");
	files.push_back("Nonroad_S505FRAC.tif");
	files.push_back("Nonroad_S510FRAC.tif");
	files.push_back("Nonroad_S520FRAC.tif");
	files.push_back("Nonroad_S525FRAC_Com_Ind_Inst.tif");
	files.push_back("Nonroad_S525FRAC_Golf.tif");
	files.push_back("Nonroad_S850FRAC.tif");
	files.push_back("Nonroad_S860FRAC.tif");
	files.push_back("Nonroad_S890FRAC.tif");

	sectorfiles.push_back(files);


	sumBySectorRaster(dir, sectorfiles, sectornames, outcsvfile);
}
#include <qfileinfo.h>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision
void Vulcan2014::sumBySectorShape(std::string dir, std::vector<std::vector<std::string>> sectorfiles, std::vector<std::string> sectornames, std::string outcsvfile)
{
	std::ofstream ofs;
	ofs.open(outcsvfile);
	double total = 0;
	for (size_t i = 0; i < sectorfiles.size(); i++)
	{
		double subtotal = 0;
		std::ofstream ofs_sector;
		if (sectorfiles[i].size() > 1)
		{
			ofs_sector.open(dir + sectornames[i] + "_total.csv");
			ofs_sector << "sector_file,total" << std::endl;
		}
		for (size_t j = 0; j < sectorfiles[i].size(); j++)
		{
			double filesum = 0;
			std::string shapefilename = sectorfiles[i][j];
			std::string dbffilename = dir + shapefilename.substr(0, shapefilename.size() - 4) + ".dbf";
			std::string csvfilename = dbffilename + ".csv";
			if (QFileInfo(csvfilename.data()).exists())
			{
				std::ifstream ifs;
				ifs.open(csvfilename.data());
				std::string line;
				getline(ifs, line);
				getline(ifs, line);
				filesum = atof(line.data());
				ifs.close();
			//else
			//{
			//	std::stringstream ss;
			//	ss << "RScript B:/Shapefiles2Grid/copyfieldmatch.R " << dbffilename << " " << "ca11";
			//	for (size_t j = 0; j < fields2keep.size(); j++)
			//	{
			//		ss << " " << fields2keep[j];
			//	}
			//}
			}
			else
			{
				ShapeFile shps(dir + shapefilename);
				int caidx = shps.poLayer->GetLayerDefn()->GetFieldIndex("ca11");
				OGRwkbGeometryType gtype = shps.poLayer->GetGeomType();
				int footprintIndex = -1;
				if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
				{
					footprintIndex = shps.poLayer->GetLayerDefn()->GetFieldIndex("length");
				}
				else if (gtype == wkbPolygon || gtype == wkbPolygon25D || gtype == wkbMultiPolygon)
				{
					footprintIndex = shps.poLayer->GetLayerDefn()->GetFieldIndex("area");
				}
				OGRFeature *poFeature;
				shps.poLayer->ResetReading();

				while ((poFeature = shps.poLayer->GetNextFeature()) != NULL)
				{
					OGRGeometry* geo = poFeature->GetGeometryRef();
					OGRwkbGeometryType gtype = geo->getGeometryType();

					double footprintOld = 0;
					if (footprintIndex > -1)
						footprintOld = poFeature->GetFieldAsDouble(footprintIndex);
					double footprintNew = 0;
					if (gtype == wkbLineString || gtype == wkbMultiLineString || gtype == wkbLineString25D)
					{
						footprintNew = Utils::calPolylineLength(poFeature->GetGeometryRef());
					}
					else if (gtype == wkbPolygon || gtype == wkbMultiPolygon || gtype == wkbPolygon25D)
					{
						footprintNew = Utils::calPolygonArea(poFeature->GetGeometryRef());
					}
					double fraction = 1;
					if (!isinf(footprintOld) && footprintOld > 0 && footprintIndex > -1)
					{
						fraction = footprintNew / footprintOld;
					}
					double ca = poFeature->GetFieldAsDouble(caidx) * fraction;
					filesum += ca;
				}
			}

	
			subtotal += filesum;
			printf("%s,%f\n", sectorfiles[i][j].data(), filesum);
			if (sectorfiles[i].size() > 1)
			{
				ofs_sector << std::fixed << std::setprecision(5) << shapefilename.substr(0, shapefilename.size() - 4) << "," << filesum << std::endl;
			}
		}
		if (sectorfiles[i].size() > 1)
		{
			ofs_sector.close();
		}
		total += subtotal;
		ofs << sectornames[i] << "," << subtotal / 1000000.0 << std::endl;
	}
	ofs << "Total" << "," << total / 1000000.0 << std::endl;
	ofs.close();
}
void Vulcan2014::sumBySectorShape(std::string dir, std::string outcsvfile)
{
	std::vector<std::vector<std::string>> sectorfiles;
	std::vector<std::string> sectornames;
	std::vector<std::string> files;

	sectornames.push_back("Electricity");
	files.clear();
	files.push_back("ElecProd.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("OnRoad");
	files.clear();
	files.push_back("OnRoadRuralLocal.dbf");
	files.push_back("OnRoadRuralNonLocal.dbf");
	files.push_back("OnRoadUrbanLocal.dbf");
	files.push_back("OnRoadUrbanNonLocal.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("Commercial");
	files.clear();
	files.push_back("ComNonPoint.dbf");
	files.push_back("ComPoint.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("Residential");
	files.clear();
	files.push_back("ResNonPoint.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("Industrial");
	files.clear();
	files.push_back("IndNonPoint.dbf");
	files.push_back("IndPoint.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("Airport");
	files.clear();
	files.push_back("Airport.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("Cement");
	files.clear();
	files.push_back("Cement.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("CMV");
	files.clear();
	files.push_back("CMVPort.dbf");
	files.push_back("CMVShippingLane.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("Railroad");
	files.clear();
	files.push_back("Railroad.dbf");
	files.push_back("RailroadPoint.dbf");
	sectorfiles.push_back(files);

	sectornames.push_back("Nonroad");
	files.clear();
	files.push_back("NonroadPoint.dbf");
	files.push_back("Nonroad_County.dbf");
	files.push_back("Nonroad_S100FRAC.dbf");
	files.push_back("Nonroad_S140FRAC.dbf");
	files.push_back("Nonroad_S260FRAC.dbf");
	files.push_back("Nonroad_S300FRAC.dbf");
	files.push_back("Nonroad_S310FRAC.dbf");
	files.push_back("Nonroad_S311FRAC.dbf");
	files.push_back("Nonroad_S350FRAC.dbf");
	files.push_back("Nonroad_S400FRAC.dbf");
	files.push_back("Nonroad_S505FRAC.dbf");
	files.push_back("Nonroad_S510FRAC.dbf");
	files.push_back("Nonroad_S520FRAC.dbf");
	files.push_back("Nonroad_S525FRAC_Com_Ind_Inst.dbf");
	files.push_back("Nonroad_S525FRAC_Golf.dbf");
	files.push_back("Nonroad_S850FRAC.dbf");
	files.push_back("Nonroad_S860FRAC.dbf");
	files.push_back("Nonroad_S890FRAC.dbf");

	sectorfiles.push_back(files);


	sumBySectorShape(dir, sectorfiles, sectornames, outcsvfile);
}
