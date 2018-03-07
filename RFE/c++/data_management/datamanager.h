#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_port.h"
#include <stdarg.h>
#include <stddef.h>
#include "cpl_string.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

class dataManager
{
/*
 * Core Class for everything related to data handling, main topics are
 *
 * 1. Parsing the input directory, generating a list INPUT_LIST of filenames related to overlapping GTifs
 * 2. Generating a Car-Free-Mosaic MOSAIC.gtif using all files of INPUT_LIST
 * 3. Tiling the MOSAIC.gtif into a bigger number of quadratic tiles
 *    a) preventing us to run into memory issues
 *    b) enabling the parallelized image processing
 * 4. Combining several Output Tiles of the Image Processing Module to one single Tile GTif
 *
 */

public:

    // Constructor if a directory + file Extension is passed, creates Car-Free-Mosaic
    dataManager(const boost::filesystem::path& PathToInputTiffs, const string& FileExtension);

    // Directory / Filenames of temporary data, e.g. Mosaic, Tiles
    boost::filesystem::path tmpDir;
    boost::filesystem::path tmpMosaicDir;
    boost::filesystem::path tmpMosaicFile;

    boost::filesystem::path tmpTilesNoBorderDir;

    // Vector of pathes to all Input GTifs for Mosaicing, orthorectified 3K data for instance
    vector<boost::filesystem::path> InputTifs;

    // Vector of pathes + filenames to all Tiles for further image processing
    vector<boost::filesystem::path> TilesTifs;

    void createComposite(const boost::filesystem::path &PathToInputTiffs);

    void tileData(const int &TileSize, const int &BorderSize);

private:

    int getTilePosition(int StartRow, int StartColumn, int EndRow, int EndColumn,
                                      int NbrColumns, int NbrRows);

    void readTileData(int &StartRow, int &StartColumn, int &EndRow, int &EndColumn,
                      GDALDataset * poDataset, const int &Border);

    void writeTileData(int &StartRow, int &StartColumn, int &NbrRows, int &NbrColumns,
                       float * TileData,
                       double adfGeoTransform[6],
                       GDALDataset * srcDataset, std::string Prefix);

    void ParseDirForTifs(const boost::filesystem::path& PathToInputTiffs, const string& FileExtension);

    void computePixelValue(std::vector<std::string> &VecGeoTiffs, double &Current_R, double &Current_H,
                           float &PixelValueRed, float &PixelValueGreen, float &PixelValueBlue);
};

#endif // DATAMANAGER_H
