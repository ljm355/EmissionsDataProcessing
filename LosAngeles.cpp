#include "LosAngeles.h"
#include "qfile.h"


LosAngeles::LosAngeles()
{

}


LosAngeles::~LosAngeles()
{

}

void LosAngeles::mergeNonPoint()
{
	std::string indir = "B:/LA_Version2/gridPrep_SHP_master/";
	std::string outdir = "B:/LA_Version2/gridPrep_SHP_master/Vulcan/Spatial/";
	std::string folders[5] = { "Los_Angeles","Orange","Riverside","San_Bernardino","Ventura" };
	std::string fips[5] = { "06037","06059","06065 ","06071","06111" };
	//std::string fields[9] = { "yearbuilt","tfs","bt","ca_ng","ca_p","ca_c","timestruct","neeui_p","neeui_ng" };
	std::vector<std::string> fields = Utils::buildVector("", new std::string[9]{ "yearbuilt","tfs","bt","ca_ng","ca_p","ca_c","timestruct","neeui_p","neeui_ng" }, 9);
	std::vector<std::string> sectors = Utils::buildVector("",new std::string[3]{ "ResNonPoint","ComNonPoint","IndNonPoint" }, 3);
	for (size_t sector = 0; sector < 3; sector++)
	{
		std::string targetfile = outdir + sectors[sector] + ".shp";
		if (!QFile(targetfile.data()).exists())
		{
			ShapeFile::copy(indir + folders[0] + "/" + sectors[sector] + ".shp", targetfile, fields);
			ShapeFile targetshp(targetfile, 1);
			OGRFeature* fea = NULL;
			int fipsfield = targetshp.getOrCreateField("fips", OGRFieldType::OFTString);
			while ((fea = targetshp.poLayer->GetNextFeature()) != NULL)
			{
				fea->SetField(fipsfield, fips[0].data());
				targetshp.poLayer->SetFeature(fea);
				OGRFeature::DestroyFeature(fea);
			}
			targetshp.close();
		}
		OGRFeature* fea = NULL;
		for (size_t folder = 1; folder < 5; folder++)
		{
			ShapeFile targetshp(targetfile, 1);
			int fipsfield = targetshp.getOrCreateField("fips", OGRFieldType::OFTString);
			ShapeFile sourceshp(indir + folders[folder] + "/" + sectors[sector] + ".shp", 0);
			std::vector<int> fieldindices;
			for (size_t field = 0; field < fields.size(); field++)
			{
				fieldindices.push_back(sourceshp.poLayer->GetLayerDefn()->GetFieldIndex(fields[field].data()));
				printf("%d,%s\n", fieldindices[field], fields[field].data());
			}
			OGRFeature* fea = NULL;
			sourceshp.poLayer->ResetReading();
			while ((fea = sourceshp.poLayer->GetNextFeature()) != NULL)
			{
				OGRFeature* newfea = OGRFeature::CreateFeature(targetshp.poLayer->GetLayerDefn());
				for (size_t field = 0; field < fields.size(); field++)
				{
					newfea->SetField(field, fea->GetRawFieldRef(fieldindices[field]));
				}
				newfea->SetField(fipsfield, fips[folder].data());
				newfea->SetGeometry(fea->GetGeometryRef());
				targetshp.poLayer->CreateFeature(newfea);
				OGRFeature::DestroyFeature(newfea);
				OGRFeature::DestroyFeature(fea);
			}
			sourceshp.close();
			targetshp.close();
		}
	
	}
}


void LosAngeles::ChangeDataType()
{
	std::string indir = "B:/LA_Version2/gridPrep_SHP_master/Vulcan/Spatial/src/";
	std::string outdir = "B:/LA_Version2/gridPrep_SHP_master/Vulcan/Spatial/";
	std::vector<std::string> fields = Utils::buildVector("", new std::string[9]{ "yearbuilt","tfs","bt","ca_ng","ca_p","ca_c","timestruct","neeui_p","neeui_ng" }, 9);
	std::vector<std::string> sectors = Utils::buildVector("", new std::string[3]{ "ResNonPoint","ComNonPoint","IndNonPoint" }, 3);
	for (size_t sector = 0; sector < 3; sector++)
	{
		std::string targetfile = outdir + sectors[sector] + ".shp";
		ShapeFile targetshp;
		ShapeFile sourceshp(indir  + sectors[sector] + ".shp", 0);
		targetshp.create(targetfile, sourceshp.poLayer->GetSpatialRef(),NULL, sourceshp.poLayer->GetGeomType());
		std::vector<int> fieldindices;
		for (size_t field = 0; field < fields.size(); field++)
		{
			OGRFieldDefn* fd = sourceshp.poLayer->GetLayerDefn()->GetFieldDefn(field);
			std::string fieldname = fd->GetNameRef();
			OGRFieldDefn* newfd = fd;
			if (field < 3)
			{
				newfd = new OGRFieldDefn(fd->GetNameRef(), OGRFieldType::OFTInteger);
			}
			targetshp.poLayer->CreateField(newfd);
		}
		OGRFeature* fea = NULL;
		sourceshp.poLayer->ResetReading();
		while ((fea = sourceshp.poLayer->GetNextFeature()) != NULL)
		{
			OGRFeature* newfea = OGRFeature::CreateFeature(targetshp.poLayer->GetLayerDefn());
			for (size_t field = 0; field < fields.size(); field++)
			{
				if (field == 0)
				{
					newfea->SetField(field, fea->GetFieldAsInteger(field));
				}
				else if (field ==  1)
				{
					newfea->SetField(field, (int)fea->GetFieldAsDouble(field));
				}
				else if (field == 2)
				{
					newfea->SetField(field, (int)fea->GetFieldAsDouble(field));
				}
				else
				{
					newfea->SetField(field, fea->GetRawFieldRef(field));
				}

			}
			newfea->SetGeometry(fea->GetGeometryRef());
			targetshp.poLayer->CreateFeature(newfea);
			OGRFeature::DestroyFeature(newfea);
			OGRFeature::DestroyFeature(fea);
		}
		sourceshp.close();
		targetshp.close();

	}
}
