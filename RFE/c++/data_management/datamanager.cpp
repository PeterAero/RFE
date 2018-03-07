#include "datamanager.h"



int dataManager::getTilePosition(int StartRow, int StartColumn, int EndRow, int EndColumn,
                                  int NbrColumns, int NbrRows){

    /* Find out whether a tile is neighbouring to a border and if yes which border
     *
     *
     *
     *
     */

    enum BorderType {NoBorder = 0, Left = 1, Right = 2, Up = 3, Down = 4,
                     UpLeft = 5, UpRight = 6, LowLeft = 7, LowRight = 8};

    BorderType MyBorder;

    // if no NoData-Space around the Tile
    if(StartRow != 0 && StartColumn != 0 && EndRow != NbrRows && EndColumn != NbrColumns){
        MyBorder = NoBorder;
    // if upper left corner Tile
    }else if(StartRow == 0 && StartColumn == 0){
        MyBorder = UpLeft;
    // if upper tile, no left or right corner
    }else if(StartRow == 0 && StartColumn != 0 && EndColumn != NbrColumns){
        MyBorder = Up;
    // if upper right corner
    }else if(StartRow == 0 && StartColumn != 0 && EndColumn == NbrColumns){
        MyBorder = UpRight;
    // if left border
    }else if(StartRow != 0 && StartColumn == 0 && EndColumn != NbrColumns && EndRow != NbrRows){
        MyBorder = Left;
    // if right boarder
    }else if(StartRow != 0 && StartColumn != 0 && EndColumn == NbrColumns && EndRow != NbrRows){
        MyBorder = Right;
    // if lower border
    }else if(StartRow != 0 && EndRow == NbrRows && StartColumn != 0 && EndColumn != NbrColumns){
        MyBorder = Down;
        // if lower left corner
     }else if(StartColumn == 0 && EndRow == NbrRows){
        MyBorder = LowLeft;
    // if lower right corner
    }else if(StartRow != 0 && StartColumn != 0 && EndColumn == NbrColumns && EndRow == NbrRows){
        MyBorder = LowRight;
    }

    return MyBorder;

}


void dataManager::tileData(const int &TileSize, const int &BorderSize){

    GDALAllRegister();
    GDALDataset  * srcDataset;
    srcDataset = (GDALDataset *) GDALOpen( this->tmpMosaicFile.c_str(), GA_ReadOnly );

    // iterate over Image
    int StartRow = 0, EndRow = 0, StartColumn = 0, EndColumn = 0;

    int NbrColumns = srcDataset->GetRasterXSize();
    int NbrRows = srcDataset->GetRasterYSize();

    for(int i = 0; i < NbrRows; i += TileSize){
        for(int j = 0; j < NbrColumns; j += TileSize){

            StartRow = i;
            StartColumn = j;

            if(NbrColumns < StartColumn + TileSize){
                EndColumn = NbrColumns;
            }else{
                EndColumn = StartColumn + TileSize;
            }

            if(NbrRows < StartRow + TileSize){
                EndRow = NbrRows;
            }else{
                EndRow = StartRow + TileSize;
            }

            readTileData(StartRow, StartColumn, EndRow, EndColumn, srcDataset, BorderSize);

        }

    }

    GDALClose( (GDALDatasetH) srcDataset );

}



void dataManager::readTileData(int &StartRow, int &StartColumn, int &EndRow, int &EndColumn,
                               GDALDataset * poDataset, const int &BorderSize){

    float * TileData;
    GDALRasterBand * poBand;
    poBand = poDataset->GetRasterBand( 1 );

    int NbrRowsBorderTile = (EndRow - StartRow) + 2 * BorderSize;
    int NbrColumnsBorderTile = (EndColumn - StartColumn) + 2 * BorderSize;

    int NbrColumns = poDataset->GetRasterXSize();
    int NbrRows = poDataset->GetRasterYSize();


    TileData = (float *) CPLMalloc(sizeof(float)*NbrRowsBorderTile*NbrColumnsBorderTile);
    for(int i = 0; i < NbrRowsBorderTile*NbrColumnsBorderTile; i++){
        TileData[i] = 0.;
    }

    if(BorderSize != 0){

        int BorderPosition = this->getTilePosition(StartRow, StartColumn, EndRow, EndColumn, NbrColumns, NbrRows);

        if(BorderPosition == 0){

            poBand->RasterIO( GF_Read, StartColumn - BorderSize, StartRow - BorderSize, NbrColumnsBorderTile, NbrRowsBorderTile,
                              TileData, NbrColumnsBorderTile, NbrRowsBorderTile, GDT_Float32,
                              0, 0 );

        }else{

            float * tmpTileData;
            int Index = 0;
            int tmpRowSize = 0, tmpColumnSize = 0;
            int tmpStartRow = StartRow, tmpStartColumn = StartColumn;
            int tmpIterRowStart = 0, tmpIterColStart = 0;
            int tmpIterRowStop = NbrRowsBorderTile, tmpIterColStop = NbrColumnsBorderTile;

            if(BorderPosition == 1 || BorderPosition == 2){
                tmpRowSize = (EndRow - StartRow) + 2 * BorderSize;
                tmpColumnSize = (EndColumn - StartColumn) + 1 * BorderSize;
            }else if (BorderPosition == 3 || BorderPosition == 4){
                tmpRowSize = (EndRow - StartRow) + 1 * BorderSize;
                tmpColumnSize = (EndColumn - StartColumn) + 2 * BorderSize;
            }else{
                tmpRowSize = (EndRow - StartRow) + 1 * BorderSize;
                tmpColumnSize = (EndColumn - StartColumn) + 1 * BorderSize;
            }

            switch (BorderPosition){
            /*
             * 1 = Left
             * 2 = Right
             * 3 = Top
             * 4 = Down
             * 5 = UpLeft
             * 6 = UpRight
             * 7 = DownLeft
             * 8 = DownRight
            */
            case 1:
                tmpStartRow -=BorderSize;
                tmpIterColStart += BorderSize;
                break;
            case 2:
                tmpStartColumn -= BorderSize;
                tmpStartRow -= BorderSize;
                tmpIterColStop -= BorderSize;
                break;
            case 3:
                tmpStartColumn -= BorderSize;
                tmpIterRowStart += BorderSize;
                break;
            case 4:
                tmpStartColumn -= BorderSize;
                tmpStartRow -= BorderSize;
                break;
            case 5:
                tmpStartColumn = StartColumn;
                tmpIterRowStart += BorderSize;
                tmpIterColStart += BorderSize;
                break;
            case 6:
                tmpStartColumn -= BorderSize;
                tmpIterRowStart += BorderSize;
                tmpIterColStop -= BorderSize;
                break;
            case 7:
                tmpStartRow -= BorderSize;
                tmpIterColStart += BorderSize;
                break;
            case 8:
                tmpStartColumn -= BorderSize;
                tmpStartRow -= BorderSize;
                tmpIterColStop -= BorderSize;
                break;
            }
            tmpTileData = (float *) CPLMalloc(sizeof(float) * tmpRowSize * tmpColumnSize);
            poBand->RasterIO( GF_Read, tmpStartColumn, tmpStartRow, tmpColumnSize, tmpRowSize,
                              tmpTileData, tmpColumnSize, tmpRowSize, GDT_Float32,
                              0, 0 );


            for(int i = tmpIterRowStart; i < tmpIterRowStop; i++){

                for(int j = tmpIterColStart; j < tmpIterColStop; j++){

                    Index = j + i * NbrColumnsBorderTile;
                    TileData[Index] = tmpTileData[j - tmpIterColStart + ((i - tmpIterRowStart)*tmpColumnSize)];

                }

            }

        }

    }else{
        poBand->RasterIO( GF_Read, StartColumn, StartRow, NbrColumnsBorderTile, NbrRowsBorderTile,
                          TileData, NbrColumnsBorderTile, NbrRowsBorderTile, GDT_Float32,
                          0, 0 );
    }

    double adfGeoTransform[6];
    poDataset->GetGeoTransform( adfGeoTransform );
    adfGeoTransform[0] += (StartColumn - BorderSize) * adfGeoTransform[1] ;
    adfGeoTransform[3] += (StartRow - BorderSize) * adfGeoTransform[5];

    if (BorderSize == 0){
        writeTileData(StartRow, StartColumn, NbrRowsBorderTile, NbrColumnsBorderTile, TileData, adfGeoTransform, poDataset, "NoBorder");
    }else{
        writeTileData(StartRow, StartColumn, NbrRowsBorderTile, NbrColumnsBorderTile, TileData, adfGeoTransform, poDataset, "Border");
    }

}


void dataManager::writeTileData(int &StartRow, int &StartColumn, int &NbrRows, int &NbrColumns,
                   float * TileData,
                   double adfGeoTransform[6],
                   GDALDataset * srcDataset,
                   std::string prefix){

    const char * pszFormat = "GTiff";
    GDALDriver * poDriver;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == NULL )
        exit( 1 );

    GDALDataset * poDstDataset;
    char **papszOptions = NULL;

    std::string FileName = this->tmpTilesNoBorderDir.string() + "/";
    std::string FileSuffix (".tif");
    std::stringstream s1, s2;
    std::string StartRowStr, StartColumnStr;

    if (StartRow == 0){
        StartRowStr = "000";
    }else{
        s1 << StartRow;
        StartRowStr = s1.str();
    }
    if (StartColumn == 0){
        StartColumnStr = "000";
    }else{
        s2 << StartColumn;
        StartColumnStr = s2.str();
    }

    FileName += prefix + StartRowStr + StartColumnStr + FileSuffix;

    this->TilesTifs.push_back(FileName);
    std::cout << "Create Tile: " << FileName << std::endl;

    poDstDataset = poDriver->Create( FileName.c_str(), NbrColumns, NbrRows, 1, GDT_Byte,
                                papszOptions );

    GDALRasterBand * poDstBand;
    poDstDataset->SetGeoTransform( adfGeoTransform );
    poDstDataset->SetProjection(srcDataset->GetProjectionRef());

    poDstBand = poDstDataset->GetRasterBand(1);

    poDstBand->RasterIO( GF_Write, 0, 0, NbrColumns, NbrRows,
                      TileData, NbrColumns, NbrRows, GDT_Float32, 0, 0 );

    GDALClose( (GDALDatasetH) poDstDataset );

}


void dataManager::ParseDirForTifs(const boost::filesystem::path& root, const string& ext)
{
    if(!boost::filesystem::exists(root) || !boost::filesystem::is_directory(root)) return;

    boost::filesystem::recursive_directory_iterator it(root);
    boost::filesystem::recursive_directory_iterator endit;


    while(it != endit)
    {
        if(boost::filesystem::is_regular_file(*it) && it->path().extension() == ext){

            this->InputTifs.push_back(it->path().filename());

        }
        ++it;

    }


}

void dataManager::computePixelValue(std::vector<std::string> &VecGeoTiffs, double &Current_R, double &Current_H,
                        float &PixelValueRed, float &PixelValueGreen, float &PixelValueBlue){

    double adfGeoTransform[6];
    int Column = 0;
    int Row = 0;
    float * TileData1;
    TileData1 = (float *) CPLMalloc(sizeof(float) * 1 * 1);
    float * TileData2;
    TileData2 = (float *) CPLMalloc(sizeof(float) * 1 * 1);
    float * TileData3;
    TileData3 = (float *) CPLMalloc(sizeof(float) * 1 * 1);

    boost::numeric::ublas::matrix<float> PixelValues (VecGeoTiffs.size(), 5);

    for(unsigned int i = 0; i < VecGeoTiffs.size(); i++){
        for(int j = 0; j < 5; j++){
            PixelValues(i,j) = 0.;
        }
    }

    for(unsigned int i = 0; i < VecGeoTiffs.size(); i++){

        GDALDataset  * srcDataset;

        srcDataset = (GDALDataset *) GDALOpen( VecGeoTiffs[i].c_str(), GA_ReadOnly );

        srcDataset->GetGeoTransform( adfGeoTransform );

        Column = std::roundl((Current_R - adfGeoTransform[0]) / adfGeoTransform[1]);
        Row = std::roundl((adfGeoTransform[3] - Current_H) / adfGeoTransform[1]);

        GDALRasterBand * poBand1;
        GDALRasterBand * poBand2;
        GDALRasterBand * poBand3;


        poBand1 = srcDataset->GetRasterBand( 1 );
        poBand2 = srcDataset->GetRasterBand( 2 );
        poBand3 = srcDataset->GetRasterBand( 3 );

        if(Row < srcDataset->GetRasterYSize() && Column < srcDataset->GetRasterXSize()){
            poBand1->RasterIO( GF_Read, Column, Row, 1, 1,
                              TileData1, 1, 1, GDT_Float32,
                              0, 0 );
            poBand2->RasterIO( GF_Read, Column, Row, 1, 1,
                              TileData2, 1, 1, GDT_Float32,
                              0, 0 );
            poBand3->RasterIO( GF_Read, Column, Row, 1, 1,
                              TileData3, 1, 1, GDT_Float32,
                              0, 0 );

            PixelValues(i, 0) = *TileData1;
            PixelValues(i, 1) = *TileData2;
            PixelValues(i, 2) = *TileData3;
            PixelValues(i, 3) = (PixelValues(i, 0) + PixelValues(i, 1) + PixelValues(i, 2)) / 3.;
            for(int j = 0; j < 3; j++){
                PixelValues(i, 4) += std::pow(PixelValues(i, j) - PixelValues(i, 3), 2);
            }
            PixelValues(i, 4) = std::sqrt(PixelValues(i, 4));
        }else{
            *TileData1 = 0.;
            *TileData2 = 0.;
            *TileData3 = 0.;
        }
        GDALClose(srcDataset);


    }

    CPLFree(TileData1);
    CPLFree(TileData2);
    CPLFree(TileData3);

    if( VecGeoTiffs.size()  == 1){
        PixelValueRed = PixelValues(0, 0);
        PixelValueGreen = PixelValues(0, 1);
        PixelValueBlue = PixelValues(0, 2);
    }else{
        float MaxMean = 0., MaxStd = 0.;
        for(unsigned int i = 0; i < VecGeoTiffs.size(); i++){
            if (PixelValues(i, 3) > MaxMean){
                MaxMean = PixelValues(i, 3);
            }
            if(PixelValues(i, 4) > MaxStd){
                MaxStd = PixelValues(i, 4);
            }
        }

        float MinDistToMean = 99999.;
        float tmpFit = 0.;
        if(VecGeoTiffs.size() > 2){
            for(unsigned int i = 0; i < VecGeoTiffs.size(); i++){
                if(PixelValues(i, 3) != MaxMean){
                    tmpFit = MaxMean - PixelValues(i, 3);
                    if(tmpFit < MinDistToMean){
                        if(PixelValues(i, 4) != MaxStd){
                            if(PixelValues(i,0) != 0. && PixelValues(i, 1) != 0. && PixelValues(i,2) != 0.){
                                MinDistToMean = tmpFit;
                                PixelValueRed = PixelValues(i, 0);
                                PixelValueGreen = PixelValues(i, 1);
                                PixelValueBlue = PixelValues(i, 2);
                            }
                        }
                    }
                }
            }
        }else if(VecGeoTiffs.size()  == 2) {
            if(PixelValues(0, 4) != MaxStd){
                if(PixelValues(0, 0) != 0. && PixelValues(0, 1) != 0. && PixelValues(0, 2) != 0.){
                    PixelValueRed = PixelValues(0, 0);
                    PixelValueGreen = PixelValues(0, 1);
                    PixelValueBlue = PixelValues(0, 2);
                }
            }
            if(PixelValueRed == 0. && PixelValueGreen == 0. && PixelValueBlue == 0.){
                PixelValueRed = PixelValues(1, 0);
                PixelValueGreen = PixelValues(1, 1);
                PixelValueBlue = PixelValues(1, 2);
            }
        }else{
            PixelValueRed = 0;
            PixelValueGreen = 0;
            PixelValueBlue = 0;
        }
    }

}

void dataManager::createComposite(const boost::filesystem::path& root){

    // first -> get Extent
    double UL_R = 9999999., UL_H = 0., LR_R = 0., LR_H = 9999999.;

    GDALAllRegister();

    double adfGeoTransform[6];

    std::string AbsPath;

    for(unsigned int i = 0; i < this->InputTifs.size(); i++){

        GDALDataset  * srcDataset;

        AbsPath = root.string() + "/" + this->InputTifs[i].string();

        srcDataset = (GDALDataset *) GDALOpen( AbsPath.c_str(), GA_ReadOnly );

        srcDataset->GetGeoTransform( adfGeoTransform );

        if(adfGeoTransform [0] < UL_R){
            UL_R = adfGeoTransform [0];
        }

        if(adfGeoTransform [3] > UL_H){
            UL_H = adfGeoTransform [3];
        }

        if(adfGeoTransform [0] + adfGeoTransform [1] * srcDataset->GetRasterXSize() > LR_R){
            LR_R = adfGeoTransform [0] + adfGeoTransform [1] * srcDataset->GetRasterXSize();
        }

        if(adfGeoTransform [3] - std::fabs(adfGeoTransform [5]) * srcDataset->GetRasterYSize() < LR_H){
            LR_H = adfGeoTransform [3] - std::fabs(adfGeoTransform [5]) * srcDataset->GetRasterYSize();
        }

        GDALClose(srcDataset);

    }

    // second -> create Output Image

    GDALDataset  * srcDataset;
    AbsPath = root.string() + "/" + this->InputTifs[0].string();
    srcDataset = (GDALDataset *) GDALOpen( AbsPath.c_str(), GA_ReadOnly );
    srcDataset->GetGeoTransform( adfGeoTransform );

    adfGeoTransform[0] = UL_R;
    adfGeoTransform[3] = UL_H;
    int NbrColumns = 0, NbrRows = 0;
    NbrColumns = std::abs(std::lroundf((UL_R - LR_R) / adfGeoTransform [1]));
    NbrRows = std::abs(std::lroundf((UL_H - LR_H) / adfGeoTransform [5]));

    const char * pszFormat = "GTiff";
    GDALDriver * poDriver;
    GDALDataset * poDstDataset;
    char **papszOptions = NULL;

    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == NULL )
        exit( 1 );

    poDstDataset = poDriver->Create( this->tmpMosaicFile.c_str(), NbrColumns, NbrRows, 3, GDT_Float32,
                                papszOptions );


    GDALRasterBand * poDstBand;
    poDstDataset->SetGeoTransform( adfGeoTransform );
    poDstDataset->SetProjection(srcDataset->GetProjectionRef());


    poDstBand = poDstDataset->GetRasterBand(1);
    float * MosaikDataRed;
    MosaikDataRed = (float *) CPLMalloc(sizeof(float)*NbrRows*NbrColumns);
    float * MosaikDataGreen;
    MosaikDataGreen = (float *) CPLMalloc(sizeof(float)*NbrRows*NbrColumns);
    float * MosaikDataBlue;
    MosaikDataBlue = (float *) CPLMalloc(sizeof(float)*NbrRows*NbrColumns);

    for(int i = 0; i < NbrRows*NbrColumns; i++){
        MosaikDataRed[i] = 0.;
        MosaikDataGreen[i] = 0.;
        MosaikDataBlue[i] = 0.;
    }

    poDstBand->RasterIO( GF_Write, 0, 0, NbrColumns, NbrRows,
                      MosaikDataRed, NbrColumns, NbrRows, GDT_Float32, 0, 0 );

    GDALClose(srcDataset);
    GDALClose(poDstDataset);

    // third -> write Output Image

    double Current_R = adfGeoTransform[0];
    double Current_H = adfGeoTransform[3];

    std::vector<std::string> MatchingTiffs;
    double tmpGeoTransform[6];

    for(int i = 0; i < NbrRows; i++){
        std::cout << i << std::endl;
        Current_R = adfGeoTransform[0];
        for(int j = 0; j < NbrColumns; j++){

            for(std::vector<boost::filesystem::path>::iterator it = this->InputTifs.begin(); it != this->InputTifs.end(); ++it) {
                // super-nasty hack, (*it).string() equivalent to it->string() ...
                AbsPath = root.string() + "/" + (*it).string();

                srcDataset = (GDALDataset *) GDALOpen( AbsPath.c_str(), GA_ReadOnly );

                srcDataset->GetGeoTransform( tmpGeoTransform );

                double tmpImg_UL_H = 0., tmpImg_UL_R = 0.;
                double tmpImg_LR_H = 0., tmpImg_LR_R = 0.;
                tmpImg_UL_H = tmpGeoTransform [3];
                tmpImg_UL_R = tmpGeoTransform [0];
                tmpImg_LR_H = tmpGeoTransform [3] + tmpGeoTransform [5] * srcDataset->GetRasterYSize();
                tmpImg_LR_R = tmpGeoTransform [0] + tmpGeoTransform [1] * srcDataset->GetRasterXSize();
                GDALClose(srcDataset);
                if(Current_H <= tmpImg_UL_H && Current_H >= tmpImg_LR_H){
                    if(Current_R >= tmpImg_UL_R && Current_R <= tmpImg_LR_R){
                        MatchingTiffs.push_back(AbsPath);
                    }
                }
            }


            float PixelValueRed = 0., PixelValueGreen = 0., PixelValueBlue = 0.;
            this->computePixelValue(MatchingTiffs, Current_R, Current_H, PixelValueRed, PixelValueGreen, PixelValueBlue);
            MatchingTiffs.clear();
            MosaikDataRed[j + i * NbrColumns] = PixelValueRed;
            MosaikDataGreen[j + i * NbrColumns] = PixelValueGreen;
            MosaikDataBlue[j + i * NbrColumns] = PixelValueBlue;

            Current_R += adfGeoTransform[1];
        }
        Current_H -= adfGeoTransform[1];
    }

    GDALDataset  * mosaicDataset;
    mosaicDataset = (GDALDataset *) GDALOpen( this->tmpMosaicFile.c_str(), GA_Update );
    GDALRasterBand * poMosaicBandRed;
    poMosaicBandRed = mosaicDataset->GetRasterBand(1);
    GDALRasterBand * poMosaicBandGreen;
    poMosaicBandGreen = mosaicDataset->GetRasterBand(2);
    GDALRasterBand * poMosaicBandBlue;
    poMosaicBandBlue = mosaicDataset->GetRasterBand(3);

    poMosaicBandRed->RasterIO( GF_Write, 0, 0, NbrColumns, NbrRows,
                      MosaikDataRed, NbrColumns, NbrRows, GDT_Float32, 0, 0 );
    poMosaicBandGreen->RasterIO( GF_Write, 0, 0, NbrColumns, NbrRows,
                      MosaikDataGreen, NbrColumns, NbrRows, GDT_Float32, 0, 0 );
    poMosaicBandBlue->RasterIO( GF_Write, 0, 0, NbrColumns, NbrRows,
                      MosaikDataBlue, NbrColumns, NbrRows, GDT_Float32, 0, 0 );

    GDALClose(mosaicDataset);
    GDALClose(poMosaicBandRed);
    GDALClose(poMosaicBandGreen);
    GDALClose(poMosaicBandBlue);

}


dataManager::dataManager(const boost::filesystem::path& PathToInputTiffs, const string& FileExtension)
{
    /* Constructor of Data Manager, following tasks:
     * 1. Parse Input directory with given file extension (e.g. .tif), create vector of image files
     * 2. Create all mandatory temporary directories
     *    .tmp/
     *    .tmp/TilesNoBorderDir
     *    .tmp/CarFreeMosaic
     * 3. specify nomenclature of temporary files, e.g.
     *    a) name of car-free mosaic (e.g. MyMosaik.tif)
     *    c) name of Tiles with/without border
     *
     *
     *
     */

    // parse Input Directory for files with given Extension
    this->ParseDirForTifs(PathToInputTiffs, FileExtension);

    // create .tmp Directory for Files created during mosaicing and tiling
    std::vector<boost::filesystem::path> tmpString;
    for(auto& part : PathToInputTiffs){
            tmpString.push_back(part);
    }

    for(unsigned int i = 0; i < tmpString.size() - 1; i++){
        this->tmpDir /= tmpString[i];
    }

    this->tmpDir /= ".tmp";
    if(boost::filesystem::create_directory(this->tmpDir))
    {
        std::cout << "Directory Created: " << this->tmpDir << std::endl;
    }


    this->tmpMosaicDir = this->tmpDir;
    this->tmpMosaicDir /= "CarFreeMosaic";
    if(boost::filesystem::create_directory(this->tmpMosaicDir))
    {
        std::cout << "Directory Created: " << this->tmpMosaicDir << std::endl;
    }
    this->tmpMosaicFile = this->tmpMosaicDir;
    this->tmpMosaicFile /= "MyMosaic.tif";


    this->tmpTilesNoBorderDir = this->tmpDir;
    this->tmpTilesNoBorderDir /= "TilesNoBorderDir";
    if(boost::filesystem::create_directory(this->tmpTilesNoBorderDir))
    {
        std::cout << "Directory Created: " << this->tmpTilesNoBorderDir << std::endl;
    }


}
