#include "PeMSStation.h"



PeMSStation::PeMSStation()
{
	m_Fractions = NULL;
}


PeMSStation::~PeMSStation()
{
	if (m_Fractions)
		delete[] m_Fractions;
}
void PeMSStation::toBin(std::string filename)
{
	std::ofstream ofs(filename.data(), std::ios::binary);
	ofs.write((char*)m_Fractions, sizeof(double)*m_NumHours);
	ofs.close();

	//FILE *f;
	//f = fopen(filename.data(), "wb");
	//fwrite(&m_Fractions, sizeof(double), m_NumHours, f);
	//fclose(f);

}
void PeMSStation::toTxt(std::string filename)
{
	std::ofstream ofs(filename.data());

	for (size_t i = 0; i < m_NumHours; i++)
	{
		ofs << m_Fractions[i] << std::endl;
	}
	ofs.close();
}

void PeMSStation::fromStream(std::ifstream& ifs)
{
	ifs.read(((char*)&m_ID), sizeof(int));
	ifs.read(((char*)&m_NumHours), sizeof(int));
	ifs.read(((char*)&m_X), sizeof(double));
	ifs.read(((char*)&m_Y), sizeof(double));
	if (m_Fractions)
		delete[] m_Fractions;	
	int sized = sizeof(double);
	m_Fractions = (double*)(new char[sized*m_NumHours]);
	char* fractions = (char*)m_Fractions;
	ifs.read(fractions, sized*m_NumHours);

	//double sumField = 0;
	//for (int i = 0; i< m_NumHours; i++)
	//{
	//	sumField += m_Fractions[i];
	//}
	//printf("%f\n", sumField);
}
