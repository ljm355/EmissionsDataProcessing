#include "TimestructTool.h"
#include <fstream>
#include <sstream>
#include <qfileinfo.h>
#include <qdatetime.h>
#include <qdir.h>
#include <vector>
#include <map>
#include <fstream>

TimestructTool::TimestructTool()
{
}
TimestructTool::~TimestructTool()
{
}
bool TimestructTool::isLeapYear(int year)
{
	if ((year % 400 == 0 || year % 100 != 0) && (year % 4 == 0))
		return true;
	return false;
}
int TimestructTool::numOfHours(int year)
{
	if (isLeapYear(year))
		return 8784;
	return 8760;
}
std::vector<double> TimestructTool::readText(std::string filename, std::string& id, bool normalize)
{
	std::vector<double> fracs;
	readText(filename, fracs, id, normalize);
	return fracs;
}
std::vector<double> TimestructTool::readBin(std::string filename, std::string& id, bool normalize)
{

	std::ifstream ifs;
	ifs.open(filename, std::ios::binary);
	ifs.seekg(0, ifs.end);
	size_t fileSize = ifs.tellg();
	ifs.seekg(0, ifs.beg);
	int numhours = fileSize/sizeof(double);
	std::vector<double> fracs;
	fracs.resize(numhours);
	ifs.read((char*)(&fracs[0]), sizeof(double) * numhours);
	ifs.close();

	id = QFileInfo(filename.data()).baseName().toLocal8Bit().data();
	//normalizeFractions(fracbegin);
	if (normalize)
	{
		double sum = 0;
		for (size_t ihour = 0; ihour < numhours; ihour++)
		{
			sum += fracs[ihour];
		}
		if (sum == 0)
		{
			double frac = 1.0 / numhours;
			for (size_t ihour = 0; ihour < numhours; ihour++)
			{
				fracs[ihour] = frac;
			}
		}
		else
		{
			for (size_t ihour = 0; ihour < numhours; ihour++)
			{
				fracs[ihour] = fracs[ihour] / sum;
			}
		}
	}
	return fracs;
}
void TimestructTool::readText(std::string filename, std::vector<double>& fracs, std::string& id, bool normalize)
{
	std::ifstream ifs(filename.data());
	std::string line;
	while(std::getline(ifs, line))
	{
		fracs.push_back(atof(line.data()));
	}
	int numhours = fracs.size();
	id = QFileInfo(filename.data()).baseName().toLocal8Bit().data();
	//normalizeFractions(fracbegin);
	if (normalize)
	{
		double sum = 0;
		for (size_t ihour = 0; ihour < numhours; ihour++)
		{
			sum += fracs[ihour];
		}
		if (sum == 0)
		{
			double frac = 1.0 / numhours;
			for (size_t ihour = 0; ihour < numhours; ihour++)
			{
				fracs[ihour] = frac;
			}
		}
		else
		{
			for (size_t ihour = 0; ihour < numhours; ihour++)
			{
				fracs[ihour] = fracs[ihour] / sum;
			}
		}
	}


}
void TimestructTool::normalizeFractions(double *& fracs, int numhours)
{
	double sum = 0;
	for (size_t ihour = 0; ihour < numhours; ihour++)
	{
		sum += fracs[ihour];
	}
	if (sum == 0)
	{
		double frac = 1.0 / numhours;
		for (size_t ihour = 0; ihour < numhours; ihour++)
		{
			fracs[ihour] = frac;
		}
	}
	else
	{
		for (size_t ihour = 0; ihour < numhours; ihour++)
		{
			fracs[ihour] = fracs[ihour] / sum;
		}
	}
}
void TimestructTool::toText(std::string txtfilename, double* fracs, int numhours)
{
	std::ofstream ofsts(txtfilename.data());
	for (size_t ihour = 0; ihour < numhours; ihour++)
	{
		ofsts << fracs[ihour] << std::endl;
	}
	ofsts.close();

}
void TimestructTool::binary2Text(std::string binaryfilename,std::string outdir,std::string ext)
{
	QDir qoutdir(outdir.data());
	if (!qoutdir.exists())
		qoutdir.mkpath(".");
	outdir = (qoutdir.absolutePath() + "/").toLocal8Bit().data();
	std::ifstream filein(binaryfilename.data(), std::ios::binary);
	filein.seekg(0, filein.end);
	size_t fileSize = filein.tellg();
	filein.seekg(0, filein.beg);
	int numstructs = 0;
	filein.read((char*)&numstructs, sizeof(int) * 1);
	int numhours = (fileSize - (numstructs * 20) - sizeof(int)) / (sizeof(double)*numstructs);
	char* idbuf = new char[numstructs * 20];
	double* fractionbuf = new double[numstructs * numhours];
	filein.read(idbuf, numstructs * 20);
	filein.read((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	filein.close();
	char* pidbuf = idbuf;
	double* pfractionbuf = fractionbuf;
	for (size_t i = 0; i < numstructs; i++)
	{
		std::string tsid = pidbuf;
		std::ofstream ofsts((outdir + tsid + ext).data());
		double sumfracs = 0;
		for (size_t ihour = 0; ihour < numhours; ihour++)
		{
			sumfracs += pfractionbuf[ihour];
		   ofsts << pfractionbuf[ihour] << std::endl;
		}
		ofsts.close();
		if (i % 100 == 0)
		{
			printf("%d/%d,%f\n", i, numstructs, sumfracs);
		}
		pidbuf += 20;
		pfractionbuf += numhours;
	}

	delete[] idbuf;
	delete[] fractionbuf;
}
void TimestructTool::updateBinary(std::string infile, std::string outfile, std::string txtdir, bool normalize)
{
	QDir input_dir(txtdir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	std::vector<QFileInfo> txtfiles;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		//if (!fileInfo.fileName().endsWith(".txt"))
		//	continue;
		txtfiles.push_back(list.at(i));
	}

	FILE *filein = fopen(infile.data(), "rb");
	fseek(filein, 0, SEEK_END);
	size_t fileSize = ftell(filein);
	rewind(filein);
	size_t numstructs = 0;
	fread(&numstructs, sizeof(int), 1, filein);

	int numhours = (fileSize - (numstructs * 20) - sizeof(int)) / (sizeof(double)*numstructs);
	char* idbuf = new char[numstructs * 20];
	double* fractionbuf = new double[numstructs * numhours];
	fread(idbuf, numstructs * 20, 1, filein);
	fread(fractionbuf, numstructs * sizeof(double) * numhours, 1, filein);
	fclose(filein);
	char* pidbuf = idbuf;
	double* pfractionbuf = fractionbuf;
	std::map<std::string, double*> timestructMap;
	for (size_t i = 0; i < numstructs; i++)
	{
		std::string tsid = pidbuf;
		timestructMap[tsid] = pfractionbuf;
		pidbuf += 20;
		pfractionbuf += numhours;
	}

	//double* tmpbuf = new double[numhours];

	for (int i = 0; i < txtfiles.size(); ++i) {
		QFileInfo& fileInfo = txtfiles[i];
		std::string tsid;
		std::ifstream ifs(fileInfo.absoluteFilePath().toLocal8Bit().data());
		std::vector<double> tmpbuf;
		readText(fileInfo.absoluteFilePath().toLocal8Bit().data(), tmpbuf, tsid);
		std::map<std::string, double*>::iterator iter = timestructMap.find(tsid);
		if (iter != timestructMap.end())
		{
			printf("%s\n", tsid.data());
			//memcpy(iter->second, tmpbuf, numhours * sizeof(double));
			double* pdestfracs = iter->second;
			//double* psrcfracs = tmpbuf;
			for (size_t ihour = 0; ihour < numhours; ihour++)
			{
				*pdestfracs++ = tmpbuf[ihour];
			}
		}
	}
	std::ofstream fileout(outfile.data(), std::ios::binary);
	fileout.write((char*)&numstructs, sizeof(int));
	fileout.write(idbuf, numstructs * 20);
	fileout.write((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	fileout.close();
	//delete[] tmpbuf;
	delete[] idbuf;
	delete[] fractionbuf;
}
void TimestructTool::updateFileInBinary(std::string infile, std::string outfile, std::string txtfile, bool normalize)
{

	FILE *filein = fopen(infile.data(), "rb");
	fseek(filein, 0, SEEK_END);
	size_t fileSize = ftell(filein);
	rewind(filein);
	size_t numstructs = 0;
	fread(&numstructs, sizeof(int), 1, filein);

	int numhours = (fileSize - (numstructs * 20) - sizeof(int)) / (sizeof(double)*numstructs);
	char* idbuf = new char[numstructs * 20];
	double* fractionbuf = new double[numstructs * numhours];
	fread(idbuf, numstructs * 20, 1, filein);
	fread(fractionbuf, numstructs * sizeof(double) * numhours, 1, filein);
	fclose(filein);
	char* pidbuf = idbuf;
	double* pfractionbuf = fractionbuf;
	std::map<std::string, double*> timestructMap;
	for (size_t i = 0; i < numstructs; i++)
	{
		std::string tsid = pidbuf;
		timestructMap[tsid] = pfractionbuf;
		pidbuf += 20;
		pfractionbuf += numhours;
	}

	//double* tmpbuf = new double[numhours];
	//for (int i = 0; i < binfiles.size(); ++i) {
		QFileInfo fileInfo(txtfile.data());
		std::string tsid;
		std::ifstream ifs(fileInfo.absoluteFilePath().toLocal8Bit().data());
		std::vector<double> tmpbuf;
		readText(fileInfo.absoluteFilePath().toLocal8Bit().data(), tmpbuf, tsid);
		std::map<std::string, double*>::iterator iter = timestructMap.find(tsid);
		if (iter != timestructMap.end())
		{
			//memcpy(iter->second, tmpbuf, numhours * sizeof(double));
			double* pdestfracs = iter->second;
			//double* psrcfracs = tmpbuf;
			for (size_t ihour = 0; ihour < numhours; ihour++)
			{
				*pdestfracs++ = tmpbuf[ihour];
			}
		}
	//}
	std::ofstream fileout(outfile.data(), std::ios::binary);
	fileout.write((char*)&numstructs, sizeof(int));
	fileout.write(idbuf, numstructs * 20);
	fileout.write((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	fileout.close();
	//delete[] tmpbuf;
	delete[] idbuf;
	delete[] fractionbuf;
}
void TimestructTool::normalizeBinary(std::string infile)
{
	

	FILE *filein = fopen(infile.data(), "rb");
	fseek(filein, 0, SEEK_END);
	size_t fileSize = ftell(filein);
	rewind(filein);
	size_t numstructs = 0;
	fread(&numstructs, sizeof(int), 1, filein);
	int numhours = (fileSize - (numstructs * 20) - sizeof(int)) / (sizeof(double)*numstructs);
	char* idbuf = new char[numstructs * 20];
	double* fractionbuf = new double[numstructs * numhours];
	fread(idbuf, numstructs * 20, 1, filein);
	fread(fractionbuf, numstructs * sizeof(double) * numhours, 1, filein);
	fclose(filein);
	char* pidbuf = idbuf;
	double* pfractionbuf = fractionbuf;
	for (size_t i = 0; i < numstructs; i++)
	{
		std::string tsid = pidbuf;
		double sum = 0;
		for (size_t ihour = 0; ihour < numhours; ihour++)
		{
			sum += pfractionbuf[ihour];
		}
		if (sum != 1)
		{
			if (sum == 0)
			{
				double frac = 1.0 / numhours;
				for (size_t ihour = 0; ihour < numhours; ihour++)
				{
					pfractionbuf[ihour] = frac;
				}
			}
			else
			{
				for (size_t ihour = 0; ihour < numhours; ihour++)
				{
					pfractionbuf[ihour] = pfractionbuf[ihour] / sum;
				}
			}
		}

		pidbuf += 20;
		pfractionbuf += numhours;
	}
	std::ofstream fileout(infile.data(), std::ios::binary);
	fileout.write((char*)&numstructs, sizeof(int));
	fileout.write(idbuf, numstructs * 20);
	fileout.write((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	fileout.close();
	delete[] idbuf;
	delete[] fractionbuf;
}
void TimestructTool::timeshift(double* srcfracs, int srcyear, double *& destfracs, int destyear,int destNumHours)
{
	//double* srcfracsRepLastWeek = new double[8760 + 24 * 7];
	QDate srcDate(srcyear, 1, 1);
	QDate destDate(destyear, 1, 1);
	int lagdays = 0;
	int hoursinweek = 24 * 7;
	double* lastweek = new double[hoursinweek];
	double* plastweek = srcfracs + (8760 - hoursinweek - 1);
	for (size_t ihour = 0; ihour < hoursinweek; ihour++)
	{
		lastweek[ihour] = *plastweek++;
	}
	while (destDate.dayOfWeek() != srcDate.dayOfWeek())
	{
		srcDate = srcDate.addDays(1);
		lagdays++;
	}
	int laghours = lagdays * 24;
	

	for (size_t ihour = laghours; ihour < 8760; ihour++)
	{
		destfracs[ihour - laghours] = srcfracs[ihour];
	}

	int remainder = 8760 - laghours;

	for (size_t ihour = 0; ihour < laghours; ihour++)
	{
		destfracs[remainder + ihour] = lastweek[ihour];
	}
	if (laghours > 0)
	{
		double sum = 0;
		for (int i = 0; i < 12; i++)
		{
			sum = 0;
			for (int ioffset = -3; ioffset <= 3; ioffset++)
			{
				if (ioffset == 0)
					continue;
				sum += destfracs[remainder + i + ioffset];
			}
			destfracs[remainder + i] = sum / 6.0;
		}

		sum = 0;

		if (destNumHours == 8784)
		{
			padLastDay(destfracs);
		}
		for (size_t ihour = 0; ihour < destNumHours; ihour++)
		{
			sum += destfracs[ihour];
		}
		if (sum != 1)
		{
			if (sum == 0)
			{
				double frac = 1.0 / destNumHours;
				for (size_t ihour = 0; ihour < destNumHours; ihour++)
				{
					destfracs[ihour] = frac;
				}
			}
			else
			{
				for (size_t ihour = 0; ihour < destNumHours; ihour++)
				{
					destfracs[ihour] = destfracs[ihour] / sum;
				}
			}
		}
	}



}
double * TimestructTool::timeshift(double * srcfracs, int srcyear, int destyear)
{
	int destNumHours = numOfHours(destyear);
	double* destfracs = new double[destNumHours];
	timeshift(srcfracs, srcyear, destfracs, destyear, destNumHours);
	return destfracs;
}
void TimestructTool::timeshiftOverwrite(double* srcfracs, int srcyear, int destyear)
{
	int destNumHours = numOfHours(destyear);
	double* destfracs = new double[destNumHours];
	timeshift(srcfracs, srcyear, destfracs, destyear, destNumHours);
	/*memcpy(srcfracs, destfracs, destNumHours * sizeof(double));*/
	double* pdestfracs = destfracs;
	for (size_t ihour = 0; ihour < destNumHours; ihour++)
	{
		*srcfracs++ = *pdestfracs++;
	}
	delete[] destfracs;
}
void TimestructTool::txt2binary(int numhours, std::string txtdir, std::string binaryfilename,bool normalize)
{
	if (numhours < 8760)
	{
		printf("%d\n", numhours);
	}
	QDir input_dir(txtdir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	std::vector<QFileInfo> txtfiles;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		//if (!fileInfo.fileName().endsWith(".txt"))
		txtfiles.push_back(list.at(i));
	}
	
	int numstructs = txtfiles.size();
	char* idbuf = new char[numstructs * 20];
	memset(idbuf, 0, (size_t)numstructs * 20);
	double* fractionbuf = new double[numstructs * numhours];
	double* pfractionbuf = fractionbuf;
	char* pidbuf = idbuf;
	numstructs = 0;
	for (int i = 0; i < txtfiles.size(); i++) {
		QFileInfo& fileInfo = txtfiles[i];

		numstructs++;
		std::string tsid;
		std::vector<double> fracs = readText(fileInfo.absoluteFilePath().toLocal8Bit().data(), tsid, true);
		if (fracs.size() < numhours)
		{
			for (int hourofday = 0; hourofday < 24; hourofday++)
			{
				fracs.push_back(fracs[8764 - 24 + hourofday]);
			}
		}
		memcpy(pfractionbuf, &fracs[0], sizeof(double)*numhours);
		pfractionbuf += numhours;
		//std::string tsid = fileInfo.baseName().toLocal8Bit().data();
		//numstructs++;
		//double* fracbegin = pfractionbuf;
		//std::ifstream ifs(fileInfo.absoluteFilePath().toLocal8Bit().data());
		//std::string line;
		//std::string newid;
		//std::vector<double> newset = readText()
		//while (std::getline(ifs, line))
		//{
		//	std::getline(ifs, line);
		//	*pfractionbuf++ = atof(line.data());
		//}
		//ifs.close();

		memcpy(pidbuf, tsid.data(), tsid.size());
		pidbuf += 20;
		if (i % 100 == 0)
		{
			printf("%d/%d\n", i, txtfiles.size());
		}
		//if (tsid == "560375237")
		//	break;
	}
	//FILE *fileout = fopen(binaryfilename.data(), "wb");
	//fwrite(&numstructs, sizeof(int), 1, fileout);
	//fwrite(idbuf, numstructs * 20, 1, fileout);
	//fwrite(fractionbuf, numstructs * sizeof(double) * 8760, 1, fileout);
	//fclose(fileout);

	std::ofstream fileout(binaryfilename.data(), std::ios::binary);
	fileout.write((char*)&numstructs, sizeof(int));
	fileout.write(idbuf, numstructs * 20);
	fileout.write((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	fileout.close();

	/*filein.read((char*)&numstructs, sizeof(int) * 1);
	char* idbuf = new char[numstructs * 20];
	filein.read(idbuf, numstructs * 20);*/

	delete[] idbuf;
	delete[] fractionbuf;
}
void TimestructTool::mergeBinary(int numhours, std::string bindir, std::string binaryfilename, bool normalize)
{
	if (numhours < 8760)
	{
		printf("%d\n", numhours);
	}
	QDir input_dir(bindir.data());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks | QDir::NoDotAndDotDot);
	input_dir.setSorting(QDir::Name);
	QFileInfoList list = input_dir.entryInfoList();
	std::vector<QFileInfo> binfiles;
	for (int i = 0; i < list.size(); ++i) {
		QFileInfo fileInfo = list.at(i);
		//if (!fileInfo.fileName().endsWith(".txt"))
		binfiles.push_back(list.at(i));
	}

	int numstructs = binfiles.size();
	char* idbuf = new char[numstructs * 20];
	memset(idbuf, 0, (size_t)numstructs * 20);
	double* fractionbuf = new double[numstructs * numhours];
	double* pfractionbuf = fractionbuf;
	char* pidbuf = idbuf;
	numstructs = 0;
	for (int i = 0; i < binfiles.size(); i++) {
		QFileInfo& fileInfo = binfiles[i];

		numstructs++;
		std::string tsid;
		std::vector<double> fracs = readBin(fileInfo.absoluteFilePath().toLocal8Bit().data(), tsid, true);
		if (fracs.size() < numhours)
		{
			for (int hourofday = 0; hourofday < 24; hourofday++)
			{
				fracs.push_back(fracs[8764 - 24 + hourofday]);
			}
		}
		memcpy(pfractionbuf, &fracs[0], sizeof(double)*numhours);
		pfractionbuf += numhours;


		memcpy(pidbuf, tsid.data(), tsid.size());
		pidbuf += 20;
		if (i % 100 == 0)
		{
			printf("%d/%d\n", i, binfiles.size());
		}
	
	}
	//FILE *fileout = fopen(binaryfilename.data(), "wb");
	//fwrite(&numstructs, sizeof(int), 1, fileout);
	//fwrite(idbuf, numstructs * 20, 1, fileout);
	//fwrite(fractionbuf, numstructs * sizeof(double) * 8760, 1, fileout);
	//fclose(fileout);

	std::ofstream fileout(binaryfilename.data(), std::ios::binary);
	fileout.write((char*)&numstructs, sizeof(int));
	fileout.write(idbuf, numstructs * 20);
	fileout.write((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	fileout.close();

	/*filein.read((char*)&numstructs, sizeof(int) * 1);
	char* idbuf = new char[numstructs * 20];
	filein.read(idbuf, numstructs * 20);*/

	delete[] idbuf;
	delete[] fractionbuf;
}
void TimestructTool::timeshiftbinary(std::string infile, std::string outfile, int srcyear,int destyear,bool normalize)
{



	int destNumHours = numOfHours(destyear);
	std::ifstream filein(infile.data(), std::ios::binary);
	filein.seekg(0, filein.end);
	size_t fileSize = filein.tellg();
	filein.seekg(0, filein.beg);
	int numstructs = 0;
	filein.read((char*)&numstructs, sizeof(int) * 1);
	int numhours = (fileSize - (numstructs * 20) - sizeof(int)) / (sizeof(double)*numstructs);
	char* idbuf = new char[numstructs * 20];
	double* fractionbuf = new double[numstructs * numhours];
	filein.read(idbuf, numstructs * 20);
	filein.read((char*)&fractionbuf[0], (size_t)numstructs * sizeof(double) * numhours);
	filein.close();


	double* destfractionbuf = new double[numstructs * destNumHours];
	char* pidbuf = idbuf;
	double* pfractionbuf = fractionbuf;
	double* pdestfractionbuf = destfractionbuf;
	double* destfracs = new double[destNumHours];
	for (size_t i = 0; i < numstructs; i++)
	{
		std::string tsid = pidbuf;
		//timestructMap[tsid] = pfractionbuf;
		//timeshiftOverwrite(pfractionbuf, srcyear, destyear);
		timeshift(pfractionbuf, srcyear, destfracs, destyear, destNumHours);
		double* pdestfracs = destfracs;
		for (size_t ihour = 0; ihour < destNumHours; ihour++)
		{
			*pdestfractionbuf++ = *pdestfracs++;
		}
		pfractionbuf += numhours;
		//std::memcpy(pfractionbuf, destfracs, 8760 * 8);
		//pidbuf += 20;
		//pfractionbuf += 8760;
	}
	delete[] destfracs;
	//FILE *fileout = fopen(outfile.data(), "wb");
	//fwrite(&numstructs, sizeof(int), 1, fileout);
	//fwrite(idbuf, numstructs * 20, 1, fileout);
	//fwrite(fractionbuf, (size_t)numstructs * sizeof(double) * 8760, 1, fileout);
	//fclose(fileout);
	std::ofstream fileout(outfile.data(), std::ios::binary);
	fileout.write((char*)&numstructs, sizeof(int));
	fileout.write(idbuf, numstructs * 20);
	fileout.write((char*)&destfractionbuf[0], (size_t)numstructs * sizeof(double) * destNumHours);
	fileout.close();
	delete[] idbuf;
	delete[] fractionbuf;
	delete[] destfractionbuf;


}
double* TimestructTool::readBinary(std::string binaryfilename, std::string id,int& numhours)
{
	std::ifstream filein(binaryfilename.data(), std::ios::binary);
	int numstructs = 0;
	filein.read((char*)&numstructs, sizeof(int)*1);
	filein.seekg(0, filein.end);
	size_t fileSize = filein.tellg();
	filein.seekg(0, filein.beg);
	numhours = (fileSize - (numstructs * 20) - sizeof(int)) / (sizeof(double)*numstructs);
	char* idbuf = new char[numstructs * 20];
	filein.read(idbuf, numstructs * 20);
	char* pidbuf = idbuf;
	int numofid = -1;
	for (size_t i = 0; i < numstructs; i++)
	{
		std::string tsid = pidbuf;
		if (tsid == id)
		{
			numofid = i;
			break;
		}
		pidbuf += 20;
	}
	if(numofid == -1)
		return nullptr;
	//numofid = 70000;
	size_t offset = 4 + numstructs * 20 + numofid * numhours * sizeof(double);
	//4296380415
	//4907009684
	//4206209684
	printf("%d,%ld\n", numofid, numofid * numhours * sizeof(double));
	filein.seekg(offset, std::ios::beg);
	double* fractionbuf = new double[numhours];
	filein.read((char*)(&fractionbuf[0]), sizeof(double) * numhours);
	filein.close();
	return fractionbuf;
}
void TimestructTool::selectFromBinary2Text(std::string binaryfilename, std::string id, std::string txtfilename)
{
	int numhours;
	double* buf = readBinary(binaryfilename, id, numhours);
	if (!buf)
		return;
	toText(txtfilename, buf, numhours);
	delete[] buf;
}

void TimestructTool::padLastDay(double*& fracs)
{
	for (size_t ihour = 0; ihour < 24; ihour++)
	{
		fracs[8760 + ihour] = fracs[8760 + ihour - 24];
	}
}
