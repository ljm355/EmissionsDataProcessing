#include "BatchDOE.h"
#include "Utils.h"
#include <qdir.h>
#include <fstream>
#include <sstream>
#include "ANNSearchEngine.h"
BatchDOE::BatchDOE()
{

}


BatchDOE::~BatchDOE()
{

}

void BatchDOE::init(std::string bdindir, std::string doedir, std::string stationShapeFile)
{
	_bdindir = bdindir;
	_doedir = doedir;
	_stationShapeFile = stationShapeFile;

	_bdPrototypes.clear();
	QDir input_dir(_bdindir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	_bdindir = (input_dir.absolutePath() + "/").toLocal8Bit().data();
	QFileInfoList list = input_dir.entryInfoList();
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		if (!fileInfo.fileName().endsWith(".inp"))
			continue;
		std::string oriprototype = fileInfo.completeBaseName().toLocal8Bit().data();
		std::string prototype = oriprototype;
		formatName(prototype);
		if (prototype != oriprototype)
			QFile::rename(fileInfo.absoluteFilePath(), (bdindir + prototype + ".inp").data());
		_bdPrototypes.push_back(prototype);
	}

	ShapeFile inshp(_stationShapeFile, 1);
	int f_filename = inshp.getOrCreateField("Filename", OGRFieldType::OFTString);
	OGRFeature* fea;
	_weatherStations.clear();
	int idx = 0;
	inshp.poLayer->ResetReading();
	while ((fea = inshp.poLayer->GetNextFeature()) != NULL)
	{
		WeatherStation sta;
		sta.State = fea->GetFieldAsString("State");
		sta.Site_Name = fea->GetFieldAsString("Site_Name");
		int elev = fea->GetFieldAsInteger("Elev");
		std::stringstream ss;
		ss << elev;
		sta.Elev = ss.str();
		//sta.Filter = fea->GetFieldAsInteger("Filter");
		sta.DataFileName = fea->GetFieldAsString("Station"); //sta.State + "_" + sta.Site_Name;
		//formatName(sta.DataFileName);
		sta.Lon = fea->GetFieldAsDouble("Longitude");
		sta.Lat = fea->GetFieldAsDouble("Latitude");
		sta.ID = idx;
		bool found = false;
		//for (size_t i = 0; i < sta.DataFileName.size() - 1; i++)
		//{
		//	int len = sta.DataFileName.size() - i;
		//	std::string weatherFileName = sta.DataFileName.substr(0, len);
		//	std::string weatherFile = _doedir + "weather/" + weatherFileName + ".bin";
		//	if (QFile::exists(weatherFile.data()))
		//	{
		//		found = true;
		//		sta.DataFileName = weatherFileName;
		//		break;
		//	}
		//}
		//if (!found)
		//{
		//	sta.DataFileName = "";
		//	fea->SetField(f_filename, sta.DataFileName.data());
		//}
	
		//inshp.poLayer->SetFeature(fea);
		_weatherStations.push_back(sta);
		idx++;
		OGRFeature::DestroyFeature(fea);
	}
	inshp.close();
}

std::vector<std::string> BatchDOE::readTextFile(std::string filename)
{

	std::ifstream ifs;
	ifs.open(filename);
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(ifs, line))
	{
		lines.push_back(line);
	}
	ifs.close();
	return lines;
}

void BatchDOE::writeTextFile(std::vector<std::string>& lines, std::string filename)
{
	std::ofstream ofs;
	ofs.open(filename);
	for (size_t i = 0; i < lines.size(); i++)
	{
		ofs << lines[i] << std::endl;
	}
	ofs.close();
}

void BatchDOE::formatName(std::string & name)
{
	//int idxof = -1;AK_FAIRBANKS/EIELSON_A
	std::vector<char> chars;
	for (size_t i = 0; i < name.size(); i++)
	{
		if (name[i] == '\'')
			continue;
		if (name[i] == ' ')
			name[i] = '_';
		if (name[i] == '/')
			name[i] = '_';
		chars.push_back(name[i]);
		//if (name[i] == '[')
		//{
		//	idxof = i - 1;
		//	break;
		//}
	}
	chars.push_back('\0');
	name = &chars[0];

	//if (idxof > -1)
	//{
	//	name = name.substr(0, idxof);
	//}
}

void BatchDOE::copyPrototype(std::string bdPrototype, std::string year, WeatherStation& weatherStation)
{
	std::string outdir = _doedir + "output/" + weatherStation.DataFileName + "/";
	std::vector<std::string> lines = readTextFile(_bdindir + bdPrototype + ".inp");
	
	for (size_t i = 0; i < 100; i++)
	{
		if (lines[i].find("END-YEAR") != std::string::npos)
		{
			lines[i] = "   END-YEAR         = " + year;
		}
		else if (lines[i].find("BEGIN-YEAR") != std::string::npos)
		{
			lines[i] = "   BEGIN-YEAR       = " + year;
		}
		//else if (lines[i].find("ALTITUDE") != std::string::npos)
		//{
		//	lines[i] = "ALTITUDE = " + weatherStation.Elev;
		//}
	}

	//std::stringstream  hourlyreporting;
	//hourlyreporting << "\"" << "FM1 Total Hourly Report Block" << "\"" << " = REPORT - BLOCK" << std::endl;
	//hourlyreporting << "	LIBRARY - ENTRY " << "\"" << "FM1 Hourly Report Block" << "\"" << std::endl;
	//hourlyreporting << "	.." << std::endl;
	//hourlyreporting << std::endl;
	//hourlyreporting << "\"" << "Hourly Report 1" << "\"" << " = HOURLY - REPORT" << std::endl;
	//hourlyreporting << "	REPORT - SCHEDULE = " << "\"" << "Schedule ON/OFF" << "\"" << std::endl;
	//hourlyreporting << "	REPORT - BLOCK = " << "(" << "\"" << "FM1 Total Hourly Report Block" << "\"" << ")" << std::endl;
	//hourlyreporting << "	.." << std::endl;

	std::vector<std::string> hourlyReportlines;
	hourlyReportlines.push_back("$              Hourly Reporting");
	hourlyReportlines.push_back("$ ---------------------------------------------------------");
	hourlyReportlines.push_back("");
	hourlyReportlines.push_back(std::string("\"") + "FM1 Total Hourly Report Block" + "\"" + " = REPORT-BLOCK" );
	hourlyReportlines.push_back("   LIBRARY-ENTRY " + std::string("\"") + "FM1 Total Hourly Report Block" + std::string("\"") );
	hourlyReportlines.push_back("   ..");
	hourlyReportlines.push_back("");
	hourlyReportlines.push_back(std::string("\"") + "Hourly Report 1" + std::string("\"") + " = HOURLY-REPORT");
	hourlyReportlines.push_back("   REPORT-SCHEDULE  = " + std::string("\"") + "Schedule ON/OFF" + std::string("\""));
	hourlyReportlines.push_back("   REPORT-BLOCK     = " + std::string("(") + std::string("\"") + "FM1 Total Hourly Report Block" + std::string("\"") + ")");
	hourlyReportlines.push_back("   ..");
	hourlyReportlines.push_back("");
	hourlyReportlines.push_back("");
	hourlyReportlines.push_back("$ ---------------------------------------------------------");
	hourlyReportlines.push_back("$              THE END");
	hourlyReportlines.push_back("$ ---------------------------------------------------------");
	hourlyReportlines.push_back("");
	hourlyReportlines.push_back("END ..");
	hourlyReportlines.push_back("COMPUTE ..");
	hourlyReportlines.push_back("STOP ..");
	
	int hourlyReportStart = -1;
	int hourlyReportEnd = -1;
	for (size_t i = lines.size() - 1 ; i > lines.size() - 100; i--)
	{
		if (lines[i].find(" Hourly Reporting") != std::string::npos)
		{
			hourlyReportStart = i;
			break;
		}
	}

	std::vector<std::string> newlines;
	for (size_t i = 0; i < hourlyReportStart; i++)
	{
		newlines.push_back(lines[i]);
	}

	for (size_t i = 0; i < hourlyReportlines.size(); i++)
	{
		newlines.push_back(hourlyReportlines[i]);
	}


	writeTextFile(newlines, outdir + bdPrototype + ".inp");
	//$ ---------------------------------------------------------
	//$              Hourly Reporting
	//$-------------------------------------------------------- -
	//"FM1 Hourly Report Block" = REPORT - BLOCK
	//	LIBRARY - ENTRY "FM1 Hourly Report Block"
	//	..

	//	"Hourly Report 1" = HOURLY - REPORT
	//	REPORT - SCHEDULE = "Schedule ON/OFF"
	//	REPORT - BLOCK = ("FM1 Hourly Report Block")
	//	..
	//$-------------------------------------------------------- -
	//$              THE END
	//$-------------------------------------------------------- -

}

bool BatchDOE::parseSimResults(std::string filename, std::vector<std::string>& results)
{
	size_t bufsize = 461129;
	size_t filelen = QFileInfo(filename.data()).size();
	//std::ifstream f(filename, std::ios::ate);
	//size_t size = f.tellg();
	if (filelen /1000000.0 < 3.0)
		return false;
    //printf("%f,%f\n", bufsize/1024.0,filelen /1024.0);
	std::ifstream f(filename);
	char* buf = new char[bufsize];
	f.seekg(filelen - bufsize, std::ios::beg);
	f.read(buf, bufsize);
	f.close();
	//buf[bufsize] = '\0';
	std::stringstream ss;
	ss << buf;
	std::string txtbuf = ss.str();
	std::vector<int> linebreaks;
	std::vector<std::string> lines;
	char* pbuf = buf;
	int lastpos = 0;
	for (int pos = 0; pos <= bufsize; pos++) {
		char c = *pbuf;
		if (c == '\n' || pos + 1 == bufsize)
		{
			std::string newline = txtbuf.substr(lastpos, pos - lastpos);
			lines.push_back(newline);
			lastpos = pos + 1;
			if (newline.size() > 12 && newline.substr(0, 13) == "HOURLY REPORT")
			{
				//printf("%s\n", newline.data());
				linebreaks.push_back(lines.size() - 1);
			}
		}
		if (pos + 1 == bufsize /*|| linebreaks.size() >= 365*/)
		{	
			printf("%s\n", lines[lines.size()-1].data());
			break;
		}
		pbuf++;
	}

	for (size_t i = 0; i < linebreaks.size(); i++)
	{
		int linenum = linebreaks[i] + 9;
		//if (i == 364)
		//{
		//	printf("%d\n");
		//}
		for (size_t j = 0; j < 24; j++)
		{
			if (linenum > lines.size() - 1)
			{
				printf("%d%d\n", linenum, lines.size());
			}
			std::string& line = lines[linenum];
			std::string num = line.substr(line.size() - 11, 10);
			int totalfuel = atoi(num.data());
			//if (totalfuel > 100000)
			//{
			//	printf("");
			//}
			results.push_back(num);
			linenum++;
		}
	}
	delete[] buf;
	if (results.size() >= 8760)
		return true;
	return false;
}
void BatchDOE::run(std::string bdPrototype, std::string year, WeatherStation& sta)
{
	std::string outdir = _doedir + "output/" + sta.DataFileName + "/";
	if (QFile::exists((outdir + bdPrototype + ".sim").data()) && !QFile::exists((outdir + bdPrototype + ".csv").data()))
	{
		std::vector<std::string> results;
		if (parseSimResults(outdir + bdPrototype + ".sim", results))
		{
			writeTextFile(results, outdir + bdPrototype + ".csv");
		}
		return;
	}

	if (QFile::exists((outdir + bdPrototype + ".csv").data()) || QFile::exists((outdir + bdPrototype + ".sim").data()))
		return;
	copyPrototype(bdPrototype, year, sta);
	std::string command = _doedir + "DOE23" + " " + "EXE49m" + " ";
	command += outdir + bdPrototype + " " + sta.DataFileName;
	system(command.data());
	std::vector<std::string> results;
	if (parseSimResults(outdir + bdPrototype + ".sim", results))
	{
		writeTextFile(results, outdir + bdPrototype + ".csv");
		//QFile::remove((outdir + bdPrototype + ".sim").data());
	}
	//static std::string exts[] = {"BDL","CRT","dsn","inp","PLO","STD","ctr" };
	//for (size_t i = 0; i < 7; i++)
	//{
	//	QString outfile = (outdir + bdPrototype + "." + exts[i]).data();
	//	if(QFile::exists(outfile))
	//		QFile::remove(outfile);
	//}
}
void BatchDOE::cleanup(std::string dir)
{
	static std::string exts[] = {"CRT","dsn","inp","PLO","STD","ctr" };
	for (int itype = 0; itype < _bdPrototypes.size(); itype++)
	{
		std::string bt = _bdPrototypes[itype];
		if (QFile::exists((dir + bt + ".csv").data()))
		{
			QFile::remove((dir + bt + ".sim").data());
			QFile::remove((dir + bt + ".BDL").data());
		}
		for (size_t i = 0; i < 6; i++)
		{
			QString outfile = (dir + bt + "." + exts[i]).data();
			if (QFile::exists(outfile))
				QFile::remove(outfile);
		}
	}
}

#include <time.h> 
void BatchDOE::runAll(std::string year)
{
	time_t timer;
	time(&timer);  /* get current time; same as: timer = time(NULL)  */
	std::stringstream logfile;
	logfile << timer << ".log.txt";
	std::ofstream ofs;
	ofs.open(logfile.str().data());

	for (int istation = 0; istation < _weatherStations.size(); istation++)
	{

		WeatherStation sta = _weatherStations[istation];
		bool hasWeatherData = false;
		
		//for (size_t i = 0; i < sta.DataFileName.size() - 1; i++)
		//{
			//int len = sta.DataFileName.size() - i;
			//std::string weatherFileName = sta.DataFileName.substr(0, len);
			//int cutoff = 0;
			//while (cutoff < sta.DataFileName.size() / 2)
			//{
				std::string weatherFile = _doedir + "weather/" + sta.DataFileName + ".bin";
				if (QFile::exists(weatherFile.data()))
				{
					//sta.DataFileName = weatherFileName;
					hasWeatherData = true;
					//break;
				}
			//	cutoff++;
			//}
			
		//}
		if (!hasWeatherData)
			continue;
		std::string outdir = _doedir + "output/" + sta.DataFileName + "/";
		QDir dir(outdir.data());
		//if (dir.exists())
		//	continue;
		dir.mkpath(".");
		for (int itype = 0; itype < _bdPrototypes.size(); itype++)
		{
			std::string bt = _bdPrototypes[itype];
			//if (QFile::exists((outdir + bt + ".sim").data()))
			//	continue;
			//if (QFile::exists((outdir + bt + ".csv").data()))
			//	continue;
			//std::string outdir = _doedir + "output/" + weatherStation + "/";
			//if (QFile::exists((outdir + bdPrototype + ".csv").data()) || QFile::exists((outdir + bdPrototype + ".sim").data()))
			//	return;
			ofs << outdir << "," << bt << std::endl;
			ofs.flush();
			run(bt, year, sta);
		}
		cleanup(outdir);

	}
	ofs.close();
}
///// <summary>
///// 计算球面距离，以米为单位
///// </summary>
///// <param name="lat1"></param>
///// <param name="lon1"></param>
///// <param name="lat2"></param>
///// <param name="lon2"></param>
///// <returns></returns>
double getDistanceFromLatLonInMeter(double lat1, double lon1, double lat2, double lon2)
{
	double R = 6371.0; // Radius of the earth in km
	double dLat = (lat2 - lat1) * 3.14159 / 180;  // deg2rad below
	double dLon = (lon2 - lon1) * 3.14159 / 180;
	double a =
		sin(dLat / 2) * sin(dLat / 2) +
		cos(lat1* 3.14159 / 180) * cos(lat2* 3.14159 / 180) *
		sin(dLon / 2) * sin(dLon / 2);
	double c = 2 * atan2(sqrt(a), sqrt(1 - a));
	double d = R * c; // Distance in km
	return d;
}
BatchDOE::WeatherStation search_point;
struct less_than_key
{
	//static BatchDOE::WeatherStation point;
	inline bool operator() (const BatchDOE::WeatherStation& struct1, const  BatchDOE::WeatherStation& struct2)
	{
		double d1 = getDistanceFromLatLonInMeter(struct1.Lat, struct1.Lon, search_point.Lat, search_point.Lon);
		double d2 = getDistanceFromLatLonInMeter(struct2.Lat, struct2.Lon, search_point.Lat, search_point.Lon);
		return (d1 < d2);
	}
};
void BatchDOE::gapfill()
{

	ANNSearchEngine searchEngine;
	std::vector<OGRPoint> points;
	int numreturns = 400;
	for (size_t i = 0; i < _weatherStations.size(); i++)
	{
		WeatherStation sta = _weatherStations[i];
		OGRPoint newp;
		newp.setX(sta.Lon * 1000);
		newp.setY(sta.Lat * 1000);
		points.push_back(newp);
	}
	searchEngine.Create(points, numreturns);
	std::vector<WeatherStation> sortedStations = _weatherStations;
	for (int i = 0; i < _weatherStations.size(); i++)
	{

		WeatherStation sta = _weatherStations[i];
		OGRPoint center;
		center.setX(sta.Lon * 1000);
		center.setY(sta.Lat * 1000);
		int numresults = searchEngine.Select_Nearest_Points(center);
		search_point = sta;
		std::sort(sortedStations.begin(), sortedStations.end(), less_than_key());
		for (int itype = 0; itype < _bdPrototypes.size(); itype++)
		{
			std::string bt = _bdPrototypes[itype];
			std::string destfile = _doedir + "output/" + sta.DataFileName + "/" + bt + ".csv";
			if (!QFile::exists(destfile.data()))
			{
				bool found = false;
				for (int istation = 0; istation < numreturns; istation++)
				{
					double dist = 0;
					//int selectedIdx = -1;
					//searchEngine.Get_Selected_Point(istation, selectedIdx, dist);
					//std::string srcfile = _doedir + "output/" + _weatherStations[selectedIdx].DataFileName + "/" + bt + ".csv";
					std::string srcfile = _doedir + "output/" + sortedStations[istation].DataFileName + "/" + bt + ".csv";
	
					if (QFile::exists(srcfile.data()))
					{
				
						//if (sqrt(dist) > 1000)
						//{
						dist = getDistanceFromLatLonInMeter(sta.Lat, sta.Lon, sortedStations[istation].Lat, sortedStations[istation].Lon);
						if (dist > 1000)
						{
							printf("%d, %s\n,%s\n", (int)(dist), srcfile.data(), destfile.data());
							break;
						}
						//}
						found = true;
						QFile::copy(srcfile.data(), destfile.data());
						printf("%s\n,%s\n", srcfile.data(), destfile.data());
						break;
					}
				}
				if (!found)
				{
					printf("%s\n", sta.DataFileName + "/" + bt + ".csv");
				}
			}
		}
	}

}
void BatchDOE::findNearestStations(std::string bkFile)
{
	std::vector<OGRPoint> points;

	for (size_t i = 0; i < _weatherStations.size(); i++)
	{
		if (_weatherStations[i].Lon > 0)
		{
			printf("%s,%f,%f\n", _weatherStations[i].DataFileName.data(), _weatherStations[i].Lon, _weatherStations[i].Lat);
			_weatherStations[i].Lon = -180 - (180 - _weatherStations[i].Lon);
			printf("%s,%f,%f\n", _weatherStations[i].DataFileName.data(), _weatherStations[i].Lon, _weatherStations[i].Lat);
		}
		//WeatherStation sta = _weatherStations[i];
		//OGRPoint newp;
		//newp.setX(sta.Lon * 1000);
		//newp.setY(sta.Lat * 1000);
		//points.push_back(newp);
	}
	std::vector<WeatherStation> sortedStations = _weatherStations;

	ShapeFile inshp(bkFile, 1);
	int stationNameField = inshp.getOrCreateField("Station", OGRFieldType::OFTString);
	int stationIDField = inshp.getOrCreateField("StationID", OGRFieldType::OFTInteger);
	OGRFeature* fea;
	inshp.poLayer->ResetReading();
	int idx = 0;
	while ((fea = inshp.poLayer->GetNextFeature()) != NULL)
	{
		OGRPoint* pt = (OGRPoint*)fea->GetGeometryRef();

		search_point.Lat = pt->getY();
		search_point.Lon = pt->getX();
		if (pt->getX() > 0)
		{
			search_point.Lon = (-180 - (180 - pt->getX()));
		}
		std::sort(sortedStations.begin(), sortedStations.end(), less_than_key());
		std::string srcfile = sortedStations[0].DataFileName;
		fea->SetField(stationNameField, sortedStations[0].DataFileName.data());
		fea->SetField(stationIDField, sortedStations[0].ID);
		inshp.poLayer->SetFeature(fea);
		OGRFeature::DestroyFeature(fea);
		printf("%d\n", idx); idx++;
	}
	inshp.close();
}
void BatchDOE::output2shapes()
{
	for (int itype = 0; itype < _bdPrototypes.size(); itype++)
	{
		std::string bt = _bdPrototypes[itype];
		std::string outshapefile = _doedir + "output/" + bt + ".shp";
		ShapeFile inshp(_stationShapeFile, 0);
		ShapeFile outshp;
		OGRSpatialReference oSRS;
		oSRS.SetWellKnownGeogCS("WGS84");
		outshp.create(outshapefile, &oSRS, inshp.poLayer->GetLayerDefn(), OGRwkbGeometryType::wkbPoint);
		OGRFeature* fea;
		inshp.poLayer->ResetReading();
		while ((fea = inshp.poLayer->GetNextFeature()) != NULL)
		{
			WeatherStation sta;
			sta.State = fea->GetFieldAsString("State");
			sta.Site_Name = fea->GetFieldAsString("Site_Name");
			int elev = fea->GetFieldAsInteger("Elev");
			std::stringstream ss;
			ss << elev;
			sta.Elev = ss.str();
			sta.DataFileName = fea->GetFieldAsString("Filename"); //sta.State + "_" + sta.Site_Name;
																  //formatName(sta.DataFileName);
			sta.Lon = fea->GetFieldAsDouble("Longitude");
			sta.Lat = fea->GetFieldAsDouble("Latitude");

			std::string outdir = _doedir + "output/" + sta.DataFileName + "/";

			if (QFile::exists((outdir + bt + ".csv").data()))
			{
				outshp.poLayer->CreateFeature(fea);
			}
			OGRFeature::DestroyFeature(fea);

		}
		outshp.close();
		inshp.close();
	}

}
