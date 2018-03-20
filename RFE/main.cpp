#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "datamanager.h"
#include "imageproc.h"
#include "parallelprocessing.h"
#include "osm.h"
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>


using namespace std;




int main(int argc, char *argv[])
{
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * Initialization - create Car-Free Mosaic * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    std::cout << std::endl;
    std::cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
    std::cout << "*              Road Feature Extractor v. 1.0                      *" << std::endl;
    std::cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;

    std::cout << "Should a car-free mosaic be generated (Y/n)?: "<< std::endl;
    std::string CarFreeMosaic_YesNo;
    std::getline(cin, CarFreeMosaic_YesNo);

    std::cout << "Please specify tile size: "<< std::endl;
    std::string TileSize;
    std::getline(cin, TileSize);
    int TileSize_Int = std::stoi(TileSize);

    std::cout << "Please specify absolute path to input imagery: "<< std::endl;
    std::string PathInputImages;
    std::getline(cin, PathInputImages);


    std::cout << "Please specify absolute path to directory for temporary files: "<< std::endl;
    std::string PathTmpFiles;
    std::getline(cin, PathTmpFiles);

    std::cout << "Please specify absolute path to OpenStreetMap database: "<< std::endl;
    std::string PathOSM_DB;
    std::getline(cin, PathOSM_DB);

    std::cout << "Please specify the number of CPUs to use: "<< std::endl;
    std::string NumberOfCPUs;
    std::getline(cin, NumberOfCPUs);
    int NumberOfCPUs_Int = std::stoi(NumberOfCPUs);

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    // TODO be aware of absolute pathes and fixed variables typesetted here
    PathInputImages = "/home/fisc_p0/Desktop/data/InputImg";
    PathTmpFiles    = "/home/fisc_p0/Desktop/data/.tmp";
    PathOSM_DB      = "/home/fisc_p0/Desktop/data/GermanyRoads.sqlite";
    NumberOfCPUs_Int = 8;
    const float SpatialResolution = 0.2;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    const boost::filesystem::path PathOfTiff(PathInputImages);
    const string ImageFileExtension(".tif");

    dataManager MyMosaicManager(PathOfTiff, ImageFileExtension);
    if (CarFreeMosaic_YesNo == "Y"){
        MyMosaicManager.createComposite(PathOfTiff);
    }else{
        MyMosaicManager.tmpMosaicFile = boost::filesystem::path(PathTmpFiles + "/CarFreeMosaic/MyMosaic.tif");
    }



    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * *  Initialization - create Tiles * * * * * * *  * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    std::cout << "Start Tiling of Input Image ... "<< std::endl;
    std::cout << "Tile Size: " << TileSize << " px" << std::endl;
    int BorderSize = 0;
    float InnerRadius = 0.5;
    float OuterRadius = 2.0;
    BorderSize = int(std::roundf(OuterRadius / SpatialResolution));

    MyMosaicManager.tileData(TileSize_Int, BorderSize);
    std::cout << "Tiling finished!"<< std::endl;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * Image Processing - compute Features Tiles * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    std::cout << "Start Image Processing of Input Image ... "<< std::endl;

    boost::thread_group threadgroup;
    std::vector<imageProc *> SmartObjects;

    for(auto Tile : MyMosaicManager.TilesTifs){
        SmartObjects.push_back(new imageProc(Tile.string(), PathTmpFiles + "/DEV/" + Tile.filename().c_str()));
    }

    int KernelMaskSize = 0.;
    KernelMaskSize = int(std::roundf(OuterRadius/SpatialResolution));

    for(size_t i = 0; i < SmartObjects.size(); i++){
        SmartObjects[i]->initKernel(KernelMaskSize, KernelMaskSize);
        SmartObjects[i]->setDonutKernel(int(std::roundf(InnerRadius/SpatialResolution)),
                              int(std::roundf(OuterRadius/SpatialResolution)));
    }

    int threadCounter = 0;
    size_t iterCounter = 0;
    while(iterCounter < SmartObjects.size()){

        if (threadCounter < NumberOfCPUs_Int){

            threadgroup.add_thread(new boost::thread (&imageProc::filterImg, SmartObjects[iterCounter],
                                                      InnerRadius, OuterRadius, true));
            threadCounter++;
            iterCounter++;

            if(iterCounter == SmartObjects.size()){
                threadgroup.join_all();
            }

        }else{
            threadgroup.join_all();
            threadCounter = 0;
        }


    }

    for(size_t i = 0; i < SmartObjects.size(); i++){
        delete SmartObjects[i];
    }
    SmartObjects.clear();

    std::cout << "Image Processing finished!"<< std::endl;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * Image Processing - create OSM Tiles * * * * * * * * * * * * * */
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


    std::cout << "Start OpenStreetMap Processor ... "<< std::endl;


    std::vector<OSM *> SmartOSMObjects;
    for(auto Tile : MyMosaicManager.TilesTifs){
        SmartOSMObjects.push_back(new OSM(Tile.string(), Tile.filename().replace_extension().c_str(), PathTmpFiles, PathOSM_DB));
    }

    threadCounter = 0;
    iterCounter = 0;
    while(iterCounter < SmartOSMObjects.size()){
        if (threadCounter < NumberOfCPUs_Int){
            threadgroup.add_thread(new boost::thread (&OSM::createOSMMask, SmartOSMObjects[iterCounter]));
            threadCounter++;
            iterCounter++;
            if(iterCounter == SmartOSMObjects.size()){
                threadgroup.join_all();
            }
        }else{
            threadgroup.join_all();
            threadCounter = 0;
        }
    }


    for(size_t i = 0; i < SmartOSMObjects.size(); i++){
        delete SmartOSMObjects[i];
    }

    SmartOSMObjects.clear();

    std::cout << "OpenStreetMap Processor finished! "<< std::endl;


    return 0;
}
