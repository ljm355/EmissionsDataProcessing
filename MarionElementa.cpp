#include "MarionElementa.h"
#include "gdal_priv.h"
#include "TemporalGridding.h"
//band1 : LAT
//band2 : LON
//band3 : err_post
//band4 : err_prior
//band5 : post
//band6 : prior
//The emissions are in tkC per km ^ 2, aggregated for the period September 1, 2012 to April 30, 2013
MarionElementa::MarionElementa()
{

}


MarionElementa::~MarionElementa()
{

}
void MarionElementa::getSpatialRefInfo(std::string reffile,std::string& proj, double* adftransform)
{
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(reffile.data(), GA_ReadOnly);
	poDataset->GetGeoTransform(adftransform);
	proj = poDataset->GetProjectionRef();
	GDALClose(poDataset);
}


float* MarionElementa::readAttribute(std::string filename,std::string atrname)
{
	int file;
	int err;
	nc_open(filename.data(), NC_WRITE | NC_SHARE, &file);
	//nc_open(_filename.data(), NC_WRITE | NC_SHARE, &_file);
	NetCDFVar var;// = _vars[varIndex];
	var.name = atrname;
	err = nc_inq_varid(file, var.name.data(), &var.id);
	int xsize = 87;
	int ysize = 87;
	std::vector<size_t> startArr;
	startArr.push_back(0);
	startArr.push_back(0);
	std::vector<size_t> latlondim;
	latlondim.push_back(87);
	latlondim.push_back(87);
	size_t* start = &startArr[0];
	size_t* count = &latlondim[0];
	float* data = new float[xsize * ysize];
	err = nc_get_vara_float(file, var.id, start, count, data);

	float* datacpy = new float[xsize * ysize];
	memcpy(datacpy, data, xsize * ysize * sizeof(float));

	float* pdata = data;
	for (int irow = ysize-1; irow >= 0; irow--)
	{
		for (int icol = 0; icol < xsize; icol++)
		{
			*pdata = datacpy[irow * xsize + icol];
			*pdata++;
		}
	}
	delete[] datacpy;
	nc_close(file);
	return data;
}

void MarionElementa::crop(std::string filename, std::string maskfile, int mask, std::string outfile)
{
	GDALDataset  *poInflux = (GDALDataset *)GDALOpen(filename.data(), GA_ReadOnly);
	GDALDataset  *poMask = (GDALDataset *)GDALOpen(maskfile.data(), GA_ReadOnly);

	int xsize = poInflux->GetRasterXSize();
	int ysize = poInflux->GetRasterYSize();
	int ncells = xsize * ysize;
	float* prior = new float[ncells];
	float* post = new float[ncells];
	double adftransform[6];
	poMask->GetGeoTransform(adftransform);
	std::string proj = poMask->GetProjectionRef();
	poInflux->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, xsize, ysize,prior, xsize, ysize, GDT_Float32,0, 0);
	poInflux->GetRasterBand(2)->RasterIO(GF_Read, 0, 0, xsize, ysize, post, xsize, ysize, GDT_Float32, 0, 0);

	float* masks = new float[ncells];
	poMask->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, xsize, ysize, masks, xsize, ysize, GDT_Float32, 0, 0);
	
	float* pprior = prior;
	float* ppost = post;
	float* pmasks = masks;
	for (int i = 0; i < ncells; i++)
	{
		int pmask = (int)(*pmasks);
		if (pmask != mask)
		{
			*pprior = 0;
			*ppost = 0;
		}
		pprior++;
		ppost++;
		pmasks++;
	}

	GDALClose(poInflux);
	GDALClose(poMask);

	GDALAllRegister();
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset* pDataset = poDriver->Create(outfile.data(), xsize, ysize, 2, GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(adftransform);
	pDataset->SetProjection(proj.data());

	pDataset->GetRasterBand(1)->SetNoDataValue(0);
	pDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, xsize, ysize, prior, xsize, ysize, GDT_Float32, 0, 0);


	pDataset->GetRasterBand(2)->SetNoDataValue(0);
	pDataset->GetRasterBand(2)->RasterIO(GF_Write, 0, 0, xsize, ysize, post, xsize, ysize, GDT_Float32, 0, 0);
	GDALClose(pDataset);
	delete[] prior;
	delete[] post;
	delete[] masks;
}

void MarionElementa::dif(std::string filename, std::string outfile,std::string outcsv)
{
	float* prior = readdata(filename, 1);
	float* post = readdata(filename, 2);
	GDALDataset  *poDataset = (GDALDataset *)GDALOpen(filename.data(), GA_ReadOnly);
	int xsize = poDataset->GetRasterXSize();
	int ysize = poDataset->GetRasterYSize();
	int ncells = xsize * ysize;
	double adftransform[6];
	poDataset->GetGeoTransform(adftransform);
	std::string proj = poDataset->GetProjectionRef();
	GDALClose(poDataset);

	float* absdif = new float[ncells];
	float* perdif = new float[ncells];

	memset(absdif, 0, ncells * sizeof(float));
	memset(perdif, 0, ncells * sizeof(float));
	std::ofstream ofs;
	ofs.open(outcsv.data());
	ofs << "prior,post,absdif,perdif" << std::endl;
	for (size_t i = 0; i < ncells; i++)
	{
		if (post[i] > 0 && prior[i] > 0)
		{
			absdif[i] = post[i] - prior[i];
			perdif[i] = absdif[i] / prior[i] * 100;
			if (perdif[i] > 100)
				perdif[i] = 100;
			if (perdif[i] < -100)
				perdif[i] = -100;
			ofs << prior[i] << "," << post[i] << "," << absdif[i] << "," << perdif[i] << std::endl;
		}

	}
	ofs.close();
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset* pDataset = poDriver->Create(outfile.data(), xsize, ysize, 4, GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(adftransform);
	pDataset->SetProjection(proj.data());
	pDataset->GetRasterBand(1)->SetNoDataValue(0);
	pDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, xsize, ysize, prior, xsize, ysize, GDT_Float32, 0, 0);


	pDataset->GetRasterBand(2)->SetNoDataValue(0);
	pDataset->GetRasterBand(2)->RasterIO(GF_Write, 0, 0, xsize, ysize, post, xsize, ysize, GDT_Float32, 0, 0);


	pDataset->GetRasterBand(3)->SetNoDataValue(0);
	pDataset->GetRasterBand(3)->RasterIO(GF_Write, 0, 0, xsize, ysize, absdif, xsize, ysize, GDT_Float32, 0, 0);


	pDataset->GetRasterBand(4)->SetNoDataValue(0);
	pDataset->GetRasterBand(4)->RasterIO(GF_Write, 0, 0, xsize, ysize, perdif, xsize, ysize, GDT_Float32, 0, 0);
	GDALClose(pDataset);




	delete[] prior;
	delete[] post;
	delete[] absdif;
	delete[] perdif;


}

void MarionElementa::nc2tif(std::string ncfile, std::string tiffile, std::string reffile)
{
	std::string proj;
	double adftransform[6];
	getSpatialRefInfo(reffile, proj, adftransform);

	int ncols = 87;
	int nrows = 87;
	GDALAllRegister();
	const char *pszFormat = "GTiff";
	char **papszOptions = NULL;
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset* pDataset = poDriver->Create(tiffile.data(), ncols, nrows, 2, GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(adftransform);
	float* prior = readAttribute(ncfile,"prior");
	for (size_t i = 0; i < ncols * nrows; i++)
	{
		float val = prior[i];
		if (val <= 0 || val != val || isinf(val))
			prior[i] = 0;
	}

	float* post = readAttribute(ncfile, "post");
	for (size_t i = 0; i < ncols * nrows; i++)
	{
		float val = post[i];
		if (val <= 0 || val != val || isinf(val))
			post[i] = 0;
	}


	pDataset->GetRasterBand(1)->SetNoDataValue(0);
	pDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, ncols, nrows, prior, ncols, nrows, GDT_Float32, 0, 0);


	pDataset->GetRasterBand(2)->SetNoDataValue(0);
	pDataset->GetRasterBand(2)->RasterIO(GF_Write, 0, 0, ncols, nrows, post, ncols, nrows, GDT_Float32, 0, 0);

	GDALClose((GDALDatasetH)pDataset);

}

float * MarionElementa::readdata(std::string filename, int band)
{
	GDALDataset  *poDataset= (GDALDataset *)GDALOpen(filename.data(), GA_ReadOnly);
	int xsize = poDataset->GetRasterXSize();
	int ysize = poDataset->GetRasterYSize();
	int ncells = xsize * ysize;
	float* data = new float[ncells];
	poDataset->GetRasterBand(band)->RasterIO(GF_Read, 0, 0, xsize, ysize, data, xsize, ysize, GDT_Float32, 0, 0);
	GDALClose(poDataset);
	return data;
}


