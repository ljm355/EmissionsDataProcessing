# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# ReprojectMarion2.py
# Created on: 2016-04-11 15:32:47.00000
#   (generated by ArcGIS/ModelBuilder)
# Description: 
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy


# Local variables:
inshp = "INPUT"
outshp = "OUTPUT"

# Process: Project
arcpy.Project_management(inshp, outshp, "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]", "GCS_North_American_1983_To_WGS84", "PROJCS['NAD_1983_StatePlane_Indiana_East_FIPS_1301_USFeet',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['False_Easting',328083.3333333333],PARAMETER['False_Northing',820208.3333333331],PARAMETER['Central_Meridian',-85.66666666666667],PARAMETER['Scale_Factor',0.9999666666666667],PARAMETER['Latitude_Of_Origin',37.5],UNIT['Foot_US',0.3048006096012192]]", "NO_PRESERVE_SHAPE", "")

