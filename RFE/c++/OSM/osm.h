#ifndef OSM_H
#define OSM_H

//#include <QCoreApplication>
//#include <sqlite3.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <iostream>
#include <ogr_geometry.h>
#include <list>
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogrsf_frmts.h"
#include "gdal_alg.h"


class OSM
{

public:
    OSM(std::string PathToFile, std::string NameOfInputFile, std::string PathTmpFiles, std::string OSM_DB_Name);
    ~OSM();
    void createOSMMask();

private:

    std::string TiffFileName;
    std::string NameOfInputFile;
    std::string PathTmpFiles;
    std::string OSM_DB_Name;

    void rasterizeShp(const std::string& tmpPath_bufferImg, const std::string& tmpPath_bufferVec, std::string& tmpMaskImage);
    void bufferClip(const std::string& tmpPath_clip, const std::string& tmpPath_buffer);
    void clipDB(const std::string& AbsPathDB, const std::string& tmpPath_Rectangle,
                const std::string& tmpPath_clip);
    void createTmpLayer(OGRPolygon& MyRing, const std::string& tmpPath_Rectangle);
    void createBoundingBox(std::string& AbsPath, OGRLinearRing& MyRing);
    //void loadCompleteOSMToMem(OGRPolygon& MyRing);
    void computeAngle(double& angle, OGRPoint& Point1, OGRPoint& Point2);
    void computeFeatureClassAndBufferSize(std::string& FeatureClass, double& BufferInMeters);
    void computeBufferValue(const double& lat, double& BufferInMeters, double& BufferInDeg);
    inline double deg2rad (double degrees);

};

#endif // OSM_H
