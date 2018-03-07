#include "imageproc.h"
#include <math.h>

using namespace std;
using namespace boost::numeric::ublas;


void imageProc::initMatrix(int size1, int size2){

    this->FilterKernel.resize(2 * size1 + 1, 2 * size2 + 1);

    for (unsigned int i = 0; i < this->FilterKernel.size1(); i++){
        for (unsigned int j = 0; j < this->FilterKernel.size2(); j++){
            this->FilterKernel(i,j) = 0;
        }
    }

}
void imageProc::printCurrentKernel(){

    for(unsigned int i = 0; i < this->FilterKernel.size1(); i++){

        for(unsigned int j = 0; j < this->FilterKernel.size2(); j++){

            std::cout << this->FilterKernel(i,j) << " ";

        }

        std::cout << std::endl;

    }

}

template <typename T> void imageProc::createDonutMask(T innerRadius, T outerRadius){

    float x_center = outerRadius;
    float y_center = outerRadius;

    float dist = 0;

    for(unsigned int i = 0; i < this->FilterKernel.size1(); i++){

        for(unsigned int j = 0; j < this->FilterKernel.size2(); j++){

            dist = sqrt(pow(x_center - i, 2.) + pow(y_center - j, 2.));

            if (dist <= outerRadius &&  dist >= innerRadius){
                this->FilterKernel(i,j) = 1.;
                this->NbrWeights++;
            }

        }

    }

}


template <typename T> void imageProc::createCircleMask(T Radius){

    std::cout << this->FilterKernel << std::endl;
    int x_center = round(Radius);//round(this->NbrRowsFilterKernel/2.);
    int y_center = round(Radius);//round(this->NbrColumnsFilterKernel/2.);

    double dist = 0;

    for(unsigned int i = 0; i < this->FilterKernel.size1(); i++){

        for(unsigned int j = 0; j < this->FilterKernel.size2(); j++){

            dist = pow(x_center - i, 2.) + pow(y_center - j, 2.);

            if (dist <= pow(Radius, 2.)){
                this->FilterKernel(i,j) = 1.;
                this->NbrWeights++;
            }

        }

    }

}

void imageProc::saveImg(GDALDataset * srcDataset, float * ImgData, int& BorderSize){

    const char * pszFormat = "GTiff";
    GDALDriver * poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    GDALDataset * poDstDataset;
    char **papszOptions = NULL;
    int NbrColumns = srcDataset->GetRasterXSize() - 2*BorderSize;
    int NbrRows = srcDataset->GetRasterYSize() - 2*BorderSize;

    poDstDataset = poDriver->Create( this->OutputFileName.c_str(), NbrColumns, NbrRows, 1, GDT_Float32,
                                papszOptions );


    GDALRasterBand * poDstBand;
    double adfGeoTransform[6];
    srcDataset->GetGeoTransform( adfGeoTransform );
    adfGeoTransform[0] += BorderSize * adfGeoTransform[1];
    adfGeoTransform[3] += BorderSize * adfGeoTransform[5];
    poDstDataset->SetGeoTransform( adfGeoTransform );
    poDstDataset->SetProjection(srcDataset->GetProjectionRef());

    poDstBand = poDstDataset->GetRasterBand(1);

    poDstBand->RasterIO( GF_Write, 0, 0, NbrColumns, NbrRows,
                      ImgData, NbrColumns, NbrRows, GDT_Float32, 0, 0 );

    GDALClose(poDstDataset);
    std::cout << "Stored DEV-Image at: "<< this->OutputFileName << std::endl;

}

void imageProc::computeTPI(const boost::numeric::ublas::matrix<float>& InputData,
                           boost::numeric::ublas::matrix<float>& OutputData,
                           int& KernelMaskSize){

    boost::numeric::ublas::matrix<float> tmpArray;
    tmpArray.resize(this->FilterKernel.size1(), this->FilterKernel.size2());

    boost::numeric::ublas::matrix<float> elemProd;
    elemProd.resize(this->FilterKernel.size1(), this->FilterKernel.size2());

    float tmpMean = 0.;
    float tmpStd = 0.;

    for (size_t i = 0; i <= InputData.size1() - this->FilterKernel.size1(); ++ i){

        for (size_t j = 0; j <= InputData.size2() - this->FilterKernel.size2(); ++ j){


            tmpArray = boost::numeric::ublas::subrange(InputData,
                                                       i, i + this->FilterKernel.size1(),
                                                       j, j + this->FilterKernel.size2());
            elemProd = boost::numeric::ublas::element_prod(tmpArray, this->FilterKernel);

            for (size_t m = 0; m < elemProd.size1(); ++ m){

                for (size_t n = 0; n < elemProd.size2(); ++ n){

                    tmpMean +=elemProd(m,n);

                }

            }

            tmpMean /= this->NbrWeights;

            for (size_t m = 0; m < elemProd.size1(); ++ m){

                for (size_t n = 0; n < elemProd.size2(); ++ n){

                    tmpStd += std::pow(elemProd(m,n) - tmpMean, 2.);

                }

            }

            tmpStd = std::sqrt(tmpStd);
            if (tmpStd != 0.){
                // Versatz berücksichtigen, (i,j) jeweils + halbe Kantenstärke
                OutputData(i, j) = (InputData(i + KernelMaskSize , j + KernelMaskSize) - tmpMean)/tmpStd;

            }else{
                OutputData(i, j) = 0.;
            }

        }
    }
}

void imageProc::filterImg(const double& InnerRadius,
                          const double& OuterRadius,
                          const bool& DonutYesNo){

    int Index = 0;
    float * InputData;
    GDALDataset * SrcDataset = (GDALDataset *) GDALOpen( this->FileName.c_str(), GA_ReadOnly );
    GDALRasterBand * poBand;
    poBand = SrcDataset->GetRasterBand( 1 );
    int NbrColumns = SrcDataset->GetRasterXSize();
    int NbrRows = SrcDataset->GetRasterYSize();

    InputData = (float *) CPLMalloc(sizeof(float) * NbrRows * NbrColumns);

    poBand->RasterIO( GF_Read, 0, 0, NbrColumns, NbrRows,
                      InputData, NbrColumns, NbrRows, GDT_Float32,
                      0, 0 );

    boost::numeric::ublas::matrix<float> InputDataBoost;
    InputDataBoost.resize(NbrRows, NbrColumns);

    for(size_t i = 0; i < InputDataBoost.size1(); i++){
        for(size_t j = 0; j < InputDataBoost.size2(); j++){
            Index = j + i*NbrColumns;
            InputDataBoost(i,j) = InputData[Index];
        }
    }

    /* * * * * * * * * creation of kernel mask for filtering * * * * * * * * */
    double adfGeoTransform[6];

    SrcDataset->GetGeoTransform( adfGeoTransform );
    int KernelMaskSize = 0.;
    KernelMaskSize = int(std::roundf(OuterRadius/adfGeoTransform[1]));

    this->initMatrix(KernelMaskSize, KernelMaskSize);

    if(DonutYesNo == true){
        this->createDonutMask(int(std::roundf(InnerRadius/adfGeoTransform[1])),
                              int(std::roundf(OuterRadius/adfGeoTransform[1])));
    }else{
        this->createCircleMask(int(std::roundf(OuterRadius/adfGeoTransform[1])));
    }

    /* * * * * * * * *  Create and initialize output Array with zeros * * * * * * * * */
    float * OutputData;
    OutputData = (float *) CPLMalloc(sizeof(float) * (NbrRows-2*KernelMaskSize) * (NbrColumns-2*KernelMaskSize));
    for(int i = 0; i < (NbrRows - 2 * KernelMaskSize) * (NbrColumns - 2 * KernelMaskSize); i++){
        OutputData[i] = 0.;
    }

    boost::numeric::ublas::matrix<float> OutputDataBoost;
    OutputDataBoost.resize(NbrRows - 2*KernelMaskSize, NbrColumns - 2*KernelMaskSize);
    for(size_t i = 0; i < OutputDataBoost.size1(); i++){
        for(size_t j = 0; j < OutputDataBoost.size2(); j++){
            OutputDataBoost(i,j) = 0.;
        }
    }


    /* * * * * * * * * * * * * * computing tpi Image * * * * * * * * * * * * */
    this->computeTPI(InputDataBoost, OutputDataBoost, KernelMaskSize);

    /* * * * * * * * * * * * * * save tpi Image * * * * * * * * * * * * */
    for(size_t i = 0; i < OutputDataBoost.size1(); i++){
        for(size_t j = 0; j < OutputDataBoost.size2(); j++){
            Index = j + i*OutputDataBoost.size2();
            OutputData[Index] = OutputDataBoost(i, j);
        }
    }

    //this->saveImg(SrcDataset, OutputData, int(std::roundf(OuterRadius/adfGeoTransform[1])) );
    this->saveImg(SrcDataset, OutputData, KernelMaskSize);
    GDALClose(SrcDataset);

}


imageProc::imageProc(std::string FileName, std::string OutputFileName)
{
    GDALAllRegister();

    this->FileName = FileName;
    this->OutputFileName = OutputFileName;
    this->NbrWeights = 0;

}

imageProc::~imageProc(){

}
