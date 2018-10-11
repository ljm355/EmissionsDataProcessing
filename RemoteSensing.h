#pragma once
#include <string>
#include "ShapeFile.h"
//USGS Landsat 8 TOA Reflectance(Orthorectified) with Fmask
//
//Data availability(time)
//Apr 11, 2013 - Oct 15, 2016
//Provider
//USGS
//Tags
//landsat, usgs, l8, lc8, oli, tirs, cloud, shadow, snow, water, fmask
//ImageCollection BlockGroupID
//LANDSAT / LC8_L1T_TOA_FMASK
//Landsat 8 calibrated top - of - atmosphere reflectance with an Fmask quality band.Orthorectified scenes only.
//
//This collection appends the following cloud / shadow / snow / water mask band to the USGS Landsat 8 TOA Reflectance(Orthorectified) collection:
//
//fmask: cloud mask
//	0 = clear
//	1 = water
//	2 = shadow
//	3 = snow
//	4 = cloud
//	The cloud mask band is based on the Fmask algorithm, described in :
//
//	   Zhu, Z. and Woodcock, C.E., 2012. Object - based cloud and cloud shadow detection in Landsat imagery.Remote Sensing of Environment, 118, pp.83 - 94.
//		   Zhu, Z., Wang, S. and Woodcock, C.E., 2015. Improvement and expansion of the Fmask algorithm : cloud, cloud shadow, and snow detection for Landsats 4 - 7, 8, and Sentinel 2 images.Remote Sensing of Environment, 159, pp.269 - 277.
//		   GitHub. (2016).prs021 / fmask.[online] Available at : https ://github.com/prs021/fmask [Accessed 12 Nov. 2015].
class RemoteSensing
{
public:
	RemoteSensing();
	~RemoteSensing();
	std::string getProj(std::string filename);
	OGRSpatialReference create(std::string wkt);
	static OGRSpatialReference getWGS84();
	OGREnvelope readBound(std::string filename, OGRSpatialReference* targetSRS);
	OGREnvelope transformBound(OGREnvelope srcBound, OGRSpatialReference* sourceSRS, OGRSpatialReference* targetSRS);
	void writeBound(OGREnvelope bound, OGRSpatialReference* srs, std::string outshapefile);
	void extractUrbanBound(OGREnvelope targetBound, const char* rasterDir, const char* outdir);
	std::vector<std::string> constructShapeFilelist(std::string dir, std::string * filearr, int numoffiles);
	void gridBaltimore();
	void processBaltimore();
	void test();
	void processLandsat(std::string indir, std::string outdir);

};

