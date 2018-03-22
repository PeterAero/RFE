#include "osm.h"


inline double OSM::deg2rad (double degrees) {
    return degrees * M_PI / 180.;
}

void OSM::computeBufferValue(const double& lat, double& BufferInMeters, double& BufferInDeg){

    // transformation from floating point angle to radians
    double rad = deg2rad(lat);
    // distance of 1 degree latitude at equator in m
    static const double distAtEquatore = 69.172 * 1609.34;
    // distance of 1 degree latitude at current latitude in km
    double currentDist = rad * distAtEquatore;
    BufferInDeg = currentDist/BufferInMeters;
    BufferInDeg = 1./BufferInDeg;

}


void OSM::computeFeatureClassAndBufferSize(std::string& FeatureClass, double& BufferInMeters){

    static const std::string bridleway      = "bridleway";
    static const std::string cycleway       = "cycleway";
    static const std::string footway        = "footway";
    static const std::string living_street  = "living_street";
    static const std::string motorway       = "motorway";
    static const std::string motorway_link  = "motorway_link";
    static const std::string path           = "path";
    static const std::string pedestrian     = "pedestrian";
    static const std::string primary        = "primary";
    static const std::string primary_link   = "primary_link";
    static const std::string residential    = "residential";
    static const std::string secondary      = "secondary";
    static const std::string secondary_link = "secondary_link";
    static const std::string service        = "service";
    static const std::string steps          = "steps";
    static const std::string tertiary       = "tertiary";
    static const std::string tertiary_link  = "tertiary_link";
    static const std::string track          = "track";
    static const std::string track_grade1   = "track_grade1";
    static const std::string track_grade2   = "track_grade2";
    static const std::string track_grade3   = "track_grade3";
    static const std::string track_grade4   = "track_grade4";
    static const std::string track_grade5   = "track_grade5";
    static const std::string trunck         = "trunck";
    static const std::string trunck_link    = "trunck_link";
    static const std::string unclassified   = "unclassified";
    static const std::string unknown        = "unknown";

    if(FeatureClass.compare(bridleway) == 0){
        BufferInMeters = 3.;
    }else if(FeatureClass.compare(cycleway) == 0){
        BufferInMeters = 4.;
    }else if(FeatureClass.compare(footway) == 0){
        BufferInMeters = 2.;
    }else if(FeatureClass.compare(living_street) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(motorway) == 0){
        BufferInMeters = 30.;
    }else if(FeatureClass.compare(motorway_link) == 0){
        BufferInMeters = 40.;
    }else if(FeatureClass.compare(path) == 0){
        BufferInMeters = 4.;
    }else if(FeatureClass.compare(primary) == 0){
        BufferInMeters = 15.;
    }else if(FeatureClass.compare(primary_link) == 0){
        BufferInMeters = 15.;
    }else if(FeatureClass.compare(residential) == 0){
        BufferInMeters = 15.;
    }else if(FeatureClass.compare(secondary) == 0){
        BufferInMeters = 10.;
    }else if(FeatureClass.compare(secondary_link) == 0){
        BufferInMeters = 10.;
    }else if(FeatureClass.compare(service) == 0){
        BufferInMeters = 8.;
    }else if(FeatureClass.compare(steps) == 0){
        BufferInMeters = 3.;
        std::cout << "steps" << std::endl;
    }else if(FeatureClass.compare(tertiary) == 0){
        BufferInMeters = 7.5;
    }else if(FeatureClass.compare(tertiary_link) == 0){
        BufferInMeters = 7.5;
    }else if(FeatureClass.compare(track) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(track_grade1) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(track_grade2) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(track_grade3) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(track_grade4) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(track_grade5) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(trunck) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(trunck_link) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(unclassified) == 0){
        BufferInMeters = 5.;
    }else if(FeatureClass.compare(unknown) == 0){
        BufferInMeters = 5.;
    }else{
        BufferInMeters = 5.;
    }

}

void OSM::computeAngle(double& angle, OGRPoint& Point1, OGRPoint& Point2){

    angle = std::atan2(Point1.getY() - Point2.getY(), Point1.getX() - Point2.getX());

}
/*
void OSM::loadCompleteOSMToMem(OGRPolygon& MyRing){

    GDALDataset * poDS;
    poDS = (GDALDataset*) GDALOpenEx( "/home/peter/Downloads/hessen/gis.osm_roads_free_1.shp", GDAL_OF_VECTOR, NULL, NULL, NULL );

    if( poDS == NULL )
    {
        printf( "Open failed.\n" );
        exit( 1 );
    }

    std::cout << poDS->GetLayerCount() << std::endl;

    OGRGeometry * NewGeometry;
    NewGeometry = &MyRing;

    OGRLayer * poLayer;
    poLayer = poDS->GetLayer(0);
    NewGeometry->transformTo(poLayer->GetSpatialRef());

    poLayer->SetSpatialFilter(NewGeometry);

    OGRFeature *poFeature;
    poLayer->ResetReading();

    // * * * * * * * * * * * * * * * * print out db content * * * * * * * * * * * * * * * *
    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {

        OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
        int iField;
        for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ )
        {
            OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );
            if( poFieldDefn->GetType() == OFTInteger )
                printf( "%d,", poFeature->GetFieldAsInteger( iField ) );
            else if( poFieldDefn->GetType() == OFTInteger64 )
                printf( CPL_FRMT_GIB ",", poFeature->GetFieldAsInteger64( iField ) );
            else if( poFieldDefn->GetType() == OFTReal )
                printf( "%.3f,", poFeature->GetFieldAsDouble(iField) );
            else if( poFieldDefn->GetType() == OFTString )
                printf( "%s,", poFeature->GetFieldAsString(iField) );
            else
                printf( "%s,", poFeature->GetFieldAsString(iField) );
        }
        OGRGeometry *poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL
                && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )
        {
            OGRPoint *poPoint = (OGRPoint *) poGeometry;
            printf( "%.3f,%3.f\n", poPoint->getX(), poPoint->getY() );
        }
        else
        {
            printf( "no point geometry\n" );
        }
        OGRFeature::DestroyFeature( poFeature );
    }

    GDALClose( poDS );

}
*/

void OSM::createBoundingBox(std::string& AbsPath, OGRLinearRing& MyRing){

    /* function to create the bounding box of a given image
     * as OGRLinearRing object in geographic coordinates
     *
     * */

    // * * * * * get corner coordinates of image in image coordinate system * * * * * * *
    GDALDataset  * srcDataset;
    srcDataset = (GDALDataset *) GDALOpen( AbsPath.c_str(), GA_ReadOnly );
    double geoParam[6];
    double x[4];
    double y[4];

    srcDataset->GetGeoTransform(geoParam);
    double Cols = srcDataset->GetRasterXSize();
    double Rows = srcDataset->GetRasterYSize();

    x[0] = geoParam[0];
    x[1] = geoParam[0] + Cols * geoParam[1];
    x[2] = geoParam[0] + Cols * geoParam[1];
    x[3] = geoParam[0];
    y[0] = geoParam[3];
    y[1] = geoParam[3];
    y[2] = geoParam[3] + Rows * geoParam[5];
    y[3] = geoParam[3] + Rows * geoParam[5];

    // * * * * * * * * * * * * * * create SRS from input image * * * * * * * * * * * * * *
    OGRSpatialReference oSourceSRS;
    const char * poWKT = srcDataset->GetProjectionRef();
    char * poWKT_tmp = (char *) poWKT;
    oSourceSRS.importFromWkt(&poWKT_tmp);

    // * * * * * * * create transformation object from input image to polygon SRS * * * * * * *
    OGRCoordinateTransformation *poCT;
    poCT = OGRCreateCoordinateTransformation( &oSourceSRS, MyRing.getSpatialReference() );
    poCT->Transform(4, x, y);

    OGRPoint MyPoint1(x[0], y[0]);
    OGRPoint MyPoint2(x[1], y[1]);
    OGRPoint MyPoint3(x[2], y[2]);
    OGRPoint MyPoint4(x[3], y[3]);

    MyRing.addPoint(&MyPoint1);
    MyRing.addPoint(&MyPoint2);
    MyRing.addPoint(&MyPoint3);
    MyRing.addPoint(&MyPoint4);
    MyRing.addPoint(&MyPoint1);

    GDALClose(srcDataset);

}

void OSM::createTmpLayer(OGRPolygon& MyRing, const std::string& tmpPath_Rectangle){

    /* function which creates a shapefile consisting of a single rectangle,
     * usually the foodprint of an aerial/satellite imagery
     *
     * */

    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ){
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }

    GDALDataset *poDS;
    poDS = poDriver->Create( tmpPath_Rectangle.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ){
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer( "point_out", MyRing.getSpatialReference(), wkbPolygon, NULL );

    if( poLayer == NULL ){
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    OGRFieldDefn oField( "Name", OFTString );
    oField.SetWidth(32);
    if( poLayer->CreateField( &oField ) != OGRERR_NONE ){
        printf( "Creating Name field failed.\n" );
        exit( 1 );
    }

    OGRFeature *poFeature;
    poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
    poFeature->SetField( "Name", "whatever" );

    poFeature->SetGeometry( &MyRing );
    if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE ){
        printf( "Failed to create feature in shapefile.\n" );
        exit( 1 );
    }
    OGRFeature::DestroyFeature( poFeature );

    GDALClose( poDS );
}

void OSM::clipDB(const std::string& AbsPathDB, const std::string& tmpPath_Rectangle,
            const std::string& tmpPath_clip){

    /* function to clip the data of a given vector file (e.g. OSM DB) with respect to
     * a second vector file (e.g. footprint of an image)
     *
     *
     * */

    // * * * * * * * * * * * * * * * * * * load OSM DB * * * * * * * * * * * * * * * * * *
    GDALDataset * poDS1;
    poDS1 = (GDALDataset*) GDALOpenEx( AbsPathDB.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );

    if( poDS1 == NULL ){
        printf( "Open of DB failed.\n" );
        exit( 1 );
    }
    OGRLayer * poDBLayer;
    poDBLayer = poDS1->GetLayer(0);

    // * * * * * * * * * * * * * * * * * * load Mask Layer * * * * * * * * * * * * * * * * * *
    GDALDataset * poDS2;
    poDS2 = (GDALDataset*) GDALOpenEx( tmpPath_Rectangle.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );
    if( poDS2 == NULL ){
        printf( "Open of Rectangle file failed.\n" );
        exit( 1 );
    }
    OGRLayer * poMaskLayer;
    poMaskLayer = poDS2->GetLayer(0);

    // * * * * * * * * * * * * * * * * * * create Output Layer * * * * * * * * * * * * * * * * * *
    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ){
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }

    GDALDataset *poDS;
    poDS = poDriver->Create( tmpPath_clip.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ){
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer( "point_out", poMaskLayer->GetSpatialRef(), wkbLineString, NULL );
    if( poLayer == NULL ){
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    poDBLayer->Clip(poMaskLayer, poLayer);

    GDALClose( poDS );
    GDALClose( poDBLayer );
    GDALClose( poMaskLayer );
}

void OSM::bufferClip(const std::string& tmpPath_clip, const std::string& tmpPath_buffer){

    /* function which buffers a given shp-file holding LineStrings
     *
     *
     *
     *
     *
     *
     * */

    GDALDataset * poDS1;
    poDS1 = (GDALDataset*) GDALOpenEx( tmpPath_clip.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );

    if( poDS1 == NULL ){
        printf( "Open of clipped file failed.\n" );
        exit( 1 );
    }
    OGRLayer * poDBLayer;
    poDBLayer = poDS1->GetLayer(0);


    // * * * * * * * * * * * * * * * * * * create Output Layer * * * * * * * * * * * * * * * * * *
    const char *pszDriverName = "ESRI Shapefile";
    GDALDriver *poDriver;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName );
    if( poDriver == NULL ){
        printf( "%s driver not available.\n", pszDriverName );
        exit( 1 );
    }

    GDALDataset *poDS;
    poDS = poDriver->Create( tmpPath_buffer.c_str(), 0, 0, 0, GDT_Unknown, NULL );
    if( poDS == NULL ){
        printf( "Creation of output file failed.\n" );
        exit( 1 );
    }

    OGRLayer *poLayer;
    poLayer = poDS->CreateLayer( "point_out", poDBLayer->GetSpatialRef(), wkbPolygon, NULL );
    if( poLayer == NULL ){
        printf( "Layer creation failed.\n" );
        exit( 1 );
    }

    OGRFieldDefn oField( "Angle", OFTReal );
    oField.SetPrecision(15);
    oField.SetWidth(32);

    if( poLayer->CreateField( &oField ) != OGRERR_NONE ){
        printf( "Creating Angle field failed.\n" );
        exit( 1 );
    }

    OGRFeature *poFeature;
    poDBLayer->ResetReading();

    double BufferValueInMeters = 0.;
    double BufferValueInDeg = 0.;
    double Angle = 0.;
    std::string currentFeatureClass = "";
    // * * * * * * * * * * * * * * * * print out db content * * * * * * * * * * * * * * * *
    while( (poFeature = poDBLayer->GetNextFeature()) != NULL )
    {

        OGRGeometry * poGeometry;
        poGeometry = poFeature->GetGeometryRef();
        const char * pszFName = "fclass";
        currentFeatureClass = poFeature->GetFieldAsString(pszFName);

        computeFeatureClassAndBufferSize(currentFeatureClass, BufferValueInMeters);

        OGRLineString * poLineString;
        poLineString = (OGRLineString *) poGeometry;

        // std::cout << poLineString->getGeometryType() << std::endl;
        // std::cout << poLineString->getGeometryName() << std::endl;


        if(poLineString->getGeometryType() == 2){
            for(int i = 0; i < poLineString->getNumPoints() - 1; i++){
                OGRPoint MyPoint1, MyPoint2;
                poLineString->getPoint(i, &MyPoint1);
                poLineString->getPoint(i + 1, &MyPoint2);
                computeAngle(Angle, MyPoint1, MyPoint2);

                computeBufferValue(MyPoint1.getY(), BufferValueInMeters, BufferValueInDeg);

                OGRLineString LineString2;
                OGRLineString * poLineString2;
                LineString2.addPoint(MyPoint1.getX(), MyPoint1.getY());
                LineString2.addPoint(MyPoint2.getX(), MyPoint2.getY());
                poLineString2 = &LineString2;

                OGRGeometry * poBufferGeometry;
                poBufferGeometry = (OGRGeometry *) poLineString2;
                poBufferGeometry = poBufferGeometry->Buffer(BufferValueInDeg, 30);
                OGRLinearRing * poBufferPolygon = (OGRLinearRing *) poBufferGeometry;
                OGRFeature *poFeature;
                poFeature = OGRFeature::CreateFeature( poLayer->GetLayerDefn() );
                poFeature->SetField( "Angle", Angle );
                poBufferPolygon->closeRings();
                poFeature->SetGeometry( poBufferPolygon );

                if( poLayer->CreateFeature( poFeature ) != OGRERR_NONE )
                {
                    printf( "Failed to create feature in shapefile.\n" );
                    exit( 1 );
                }

                OGRFeature::DestroyFeature( poFeature );

            }
        }
        delete(poLineString);

    }

    GDALClose( poDS );

}

void OSM::rasterizeShp(const std::string& tmpPath_bufferImg, const std::string& tmpPath_bufferVec, std::string& tmpMaskImage){

    /* function to rasterize a given shape file
     * as input a shape file with polygons is expected, where each polygon holds an attribute "Angle" as floating point number
     * the attribute value "Angle" is then burned into the output tiff image of type float32 which shares the Geodetic Datum and
     * GSD  of the given input image.
     *
     *
     *
     * */

    //* * * * * * * * * * * get spatial extent and geoinfo of mask image * * * * * * * * * * *
    GDALDataset  *poDataset;
    poDataset = (GDALDataset *) GDALOpen( tmpMaskImage.c_str(), GA_ReadOnly );
    double adfGeoTransform[6];
    poDataset->GetGeoTransform( adfGeoTransform );

    //* * * * * * * * * * * create output file for rasterized buffer image * * * * * * * * * * *
    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    char **papszMetadata;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    papszMetadata = poDriver->GetMetadata();

    GDALDataset *poDstDS;
    char **papszOptions = NULL;
    poDstDS = poDriver->Create( tmpPath_bufferImg.c_str(), poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), 1, GDT_Float32,
                                papszOptions );
    poDstDS->SetGeoTransform( adfGeoTransform );
    poDstDS->SetProjection(poDataset->GetProjectionRef());

    GDALRasterBand *poBand;
    float * ImgData;
    ImgData = (float *) CPLMalloc(sizeof(float) * poDataset->GetRasterXSize() * poDataset->GetRasterYSize());

    for (int i = 0; i < poDataset->GetRasterXSize() * poDataset->GetRasterYSize(); i++){
        ImgData[i] = 0.0;
    }

    poBand = poDstDS->GetRasterBand(1);
    poBand->RasterIO( GF_Write, 0, 0, poDataset->GetRasterXSize(), poDataset->GetRasterYSize(),
                      ImgData, poDataset->GetRasterXSize(), poDataset->GetRasterYSize(), GDT_Float32, 0, 0 );

    /* close Dataset from which we got geoinfo and raster size */
    GDALClose( (GDALDatasetH) poDataset );

    //* * * * * * * * * * * get geometries to rasterize from vector file * * * * * * * * * * *

    GDALDataset * poDS1;
    poDS1 = (GDALDataset*) GDALOpenEx( tmpPath_bufferVec.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );

    if( poDS1 == NULL ){
        printf( "Open of clipped file failed.\n" );
        exit( 1 );
    }
    OGRLayer * poDBLayer;
    poDBLayer = poDS1->GetLayer(0);

    char** options = nullptr;

    //options = CSLSetNameValue(options, "ALL_TOUCHED", "TRUE");
    options = CSLSetNameValue(options, "ATTRIBUTE", "Angle");

    int panBandList = 1;
    GDALRasterizeLayers(poDstDS, 1, &panBandList , 1, (OGRLayerH*)&poDBLayer,  NULL, NULL, NULL, options, NULL, NULL);

    GDALClose( (GDALDatasetH) poDBLayer );
    GDALClose( (GDALDatasetH) poDstDS );

}


void OSM::createOSMMask(){

    std::cout << this->TiffFileName << std::endl;

    const std::string tmpPath_Rectangle  = this->PathTmpFiles + "/vec/" + this->NameOfInputFile + "_rectangle.shp";
    const std::string tmpPath_clip       = this->PathTmpFiles + "/vec/" + this->NameOfInputFile + "_clip.shp";
    const std::string tmpPath_buffer     = this->PathTmpFiles + "/vec/" + this->NameOfInputFile + "_buffer.shp";
    const std::string AbsPathBufferImage = this->PathTmpFiles + "/OSM_Mask/" + this->NameOfInputFile + ".tif";

    OGRLinearRing MyRing;
    OGRSpatialReference SRS;
    SRS.SetWellKnownGeogCS("WGS84");
    MyRing.assignSpatialReference(&SRS);
    createBoundingBox(this->TiffFileName, MyRing);
    OGRPolygon MyPolygon;
    MyPolygon.assignSpatialReference(&SRS);
    MyPolygon.addRing(&MyRing);
    createTmpLayer(MyPolygon, tmpPath_Rectangle);
    clipDB(this->OSM_DB_Name, tmpPath_Rectangle, tmpPath_clip);

    bufferClip(tmpPath_clip, tmpPath_buffer);
    rasterizeShp(AbsPathBufferImage, tmpPath_buffer, this->TiffFileName);

}

OSM::OSM(std::string PathToFile, std::string NameOfInputFile, std::string PathTmpFiles, std::string OSM_DB_Name){

    this->TiffFileName = PathToFile;
    this->NameOfInputFile = NameOfInputFile;
    this->PathTmpFiles = PathTmpFiles;
    this->OSM_DB_Name = OSM_DB_Name;
}

OSM::~OSM(){
}
