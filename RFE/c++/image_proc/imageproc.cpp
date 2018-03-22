#include "imageproc.h"
#include <math.h>

using namespace std;
using namespace boost::numeric::ublas;


void imageProc::initKernel(int size1, int size2){

    this->BorderSize = size1;

    this->FilterKernel.resize(2 * size1 + 1, 2 * size2 + 1);

    for (size_t i = 0; i < this->FilterKernel.size1(); i++){
        for (size_t j = 0; j < this->FilterKernel.size2(); j++){
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

void imageProc::setDonutKernel(float innerRadius, float outerRadius){

    float x_center = outerRadius;
    float y_center = outerRadius;

    float dist = 0;

    for(size_t i = 0; i < this->FilterKernel.size1(); i++){

        for(size_t j = 0; j < this->FilterKernel.size2(); j++){

            dist = std::sqrt(std::pow(x_center - i, 2.) + std::pow(y_center - j, 2.));

            if (dist <= outerRadius &&  dist >= innerRadius){
                this->FilterKernel(i,j) = 1.;
                this->NbrWeights++;
            }

        }

    }

}

void imageProc::setCircleKernel(float Radius){

    std::cout << this->FilterKernel << std::endl;
    int x_center = std::round(Radius);
    int y_center = std::round(Radius);

    double dist = 0;

    for(size_t i = 0; i < this->FilterKernel.size1(); i++){

        for(size_t j = 0; j < this->FilterKernel.size2(); j++){

            dist = std::pow(x_center - i, 2.) + std::pow(y_center - j, 2.);

            if (dist <= std::pow(Radius, 2.)){
                this->FilterKernel(i,j) = 1.;
                this->NbrWeights++;
            }

        }

    }

}

float imageProc::sigma_prefactor(float bandwidth){
    float b = bandwidth;
    return 1.0/ M_PI * std::sqrt(std::log(2.)/2.) * (std::pow( 2., b+1.) / std::pow(2., b-1.));
}

void imageProc::setGaborKernel(float frequency, float theta, float bandwidth, float sigma_x, float sigma_y,
                               float n_stds, float offset){

    // sigma_x and sigma_y influencing the spread width of the function response
    if(sigma_x == 0.){
        sigma_x = sigma_prefactor(bandwidth) / frequency;
    }

    if(sigma_y == 0.){
        sigma_y = sigma_prefactor(bandwidth) / frequency;
    }

    float value1 = 0.;
    float value2 = 0.;
    value1 = std::abs(n_stds * sigma_x * std::cos(theta));
    value2 = std::abs(n_stds * sigma_y * std::sin(theta));
    int x0 = std::ceil(std::max(value1, value2));

    value1 = std::abs(n_stds * sigma_y * std::cos(theta));
    value2 = std::abs(n_stds * sigma_x * std::sin(theta));
    int y0 = std::ceil(std::max(value1, value2));

    boost::numeric::ublas::matrix<float> y(this->FilterKernel.size1(), this->FilterKernel.size2());
    boost::numeric::ublas::matrix<float> x(this->FilterKernel.size1(), this->FilterKernel.size2());

    float PixelValue = - float(std::floor(this->FilterKernel.size1()/2.));

    for(int i = 0; i < y.size1(); ++i){
        for(int j = 0; j < y.size2(); ++j){
            y(i,j) = PixelValue;
            x(j,i) = PixelValue;
        }
        PixelValue += 1.;
    }

    boost::numeric::ublas::matrix<float> roty(this->FilterKernel.size1(), this->FilterKernel.size2());
    boost::numeric::ublas::matrix<float> rotx(this->FilterKernel.size1(), this->FilterKernel.size2());

    for(int i = 0; i < y.size1(); ++i){
        for(int j = 0; j < y.size2(); ++j){
            rotx(i,j) = x(i,j) * std::cos(theta) + y(i,j) * std::sin(theta);
            roty(i,j) = -x(i,j) * std::sin(theta) + y(i,j) * std::cos(theta);
        }
    }
    float MyImag = 0.;

    for(int i = 0; i < this->FilterKernel.size1(); ++i){
        for(int j = 0; j < this->FilterKernel.size2(); ++j){
            this->FilterKernel(i, j) = std::exp(-0.5 * (std::pow(rotx(i, j), 2) / std::pow(sigma_x, 2)
                                                      + std::pow(roty(i, j), 2) / std::pow(sigma_y, 2)));
            this->FilterKernel(i, j) /= 2 * M_PI * sigma_x * sigma_y;
            MyImag = (2 * M_PI * frequency * rotx(i,j) + offset);
            this->FilterKernel(i, j) *= std::cos(MyImag) + 0.*std::sin(MyImag);

        }
    }

    this->NbrWeights = this->FilterKernel.size1() * this->FilterKernel.size2();

}

void imageProc::saveImg(float * ImgData){

    GDALDataset * poSrcDataset = (GDALDataset *) GDALOpen( this->InputFileName.c_str(), GA_ReadOnly );

    const char * pszFormat = "GTiff";
    GDALDriver * poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    GDALDataset * poDstDataset;
    char **papszOptions = NULL;
    int NbrColumns = poSrcDataset->GetRasterXSize() - 2*this->BorderSize;
    int NbrRows = poSrcDataset->GetRasterYSize() - 2*this->BorderSize;

    poDstDataset = poDriver->Create( this->OutputFileName.c_str(), NbrColumns, NbrRows, 1, GDT_Float32,
                                papszOptions );


    GDALRasterBand * poDstBand;
    double adfGeoTransform[6];
    poSrcDataset->GetGeoTransform( adfGeoTransform );
    adfGeoTransform[0] += this->BorderSize * adfGeoTransform[1];
    adfGeoTransform[3] += this->BorderSize * adfGeoTransform[5];
    poDstDataset->SetGeoTransform( adfGeoTransform );
    poDstDataset->SetProjection(poSrcDataset->GetProjectionRef());

    poDstBand = poDstDataset->GetRasterBand(1);

    poDstBand->RasterIO( GF_Write, 0, 0, NbrColumns, NbrRows,
                      ImgData, NbrColumns, NbrRows, GDT_Float32, 0, 0 );

    GDALClose(poDstDataset);
    GDALClose(poSrcDataset);
    std::cout << "Stored Image at: "<< this->OutputFileName << std::endl;

}

void imageProc::computeTPI(const boost::numeric::ublas::matrix<float>& InputData,
                           const boost::numeric::ublas::matrix<float>& MaskData,
                           boost::numeric::ublas::matrix<float>& OutputData){

    boost::numeric::ublas::matrix<float> tmpArray;
    tmpArray.resize(this->FilterKernel.size1(), this->FilterKernel.size2());

    boost::numeric::ublas::matrix<float> elemProd;
    elemProd.resize(this->FilterKernel.size1(), this->FilterKernel.size2());

    float tmpMean = 0.;
    float tmpStd = 0.;

    for (size_t i = 0; i <= InputData.size1() - this->FilterKernel.size1(); ++ i){

        for (size_t j = 0; j <= InputData.size2() - this->FilterKernel.size2(); ++ j){

            if(MaskData(i + this->BorderSize , j + this->BorderSize) != 0){

                tmpMean = 0.;
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
                    OutputData(i, j) = (InputData(i + this->BorderSize , j + this->BorderSize) - tmpMean)/tmpStd;

                }else{
                    OutputData(i, j) = 0.;
                }

            }

        }
    }
}

void imageProc::computeGabor(const boost::numeric::ublas::matrix<float>& InputData,
                             const boost::numeric::ublas::matrix<float>& MaskData,
                           boost::numeric::ublas::matrix<float>& OutputData){

    boost::numeric::ublas::matrix<float> tmpArray;
    tmpArray.resize(this->FilterKernel.size1(), this->FilterKernel.size2());

    boost::numeric::ublas::matrix<float> elemProd;
    elemProd.resize(this->FilterKernel.size1(), this->FilterKernel.size2());

    float tmpMean = 0.;

    for (size_t i = 0; i <= InputData.size1() - this->FilterKernel.size1(); ++ i){

        for (size_t j = 0; j <= InputData.size2() - this->FilterKernel.size2(); ++ j){

            if(MaskData(i,j) != 0){
                this->setGaborKernel(0.25, MaskData(i,j), 30.0*M_PI/180.0, 0.58/0.25, 0.58/0.25, 3.0, 0.);
                tmpMean = 0.;
                tmpArray = boost::numeric::ublas::subrange(InputData,
                                                           i, i + this->FilterKernel.size1(),
                                                           j, j + this->FilterKernel.size2());
                elemProd = boost::numeric::ublas::element_prod(tmpArray, this->FilterKernel);

                for (size_t m = 0; m < elemProd.size1(); ++ m){

                    for (size_t n = 0; n < elemProd.size2(); ++ n){

                        tmpMean +=elemProd(m,n);

                    }

                }
                OutputData(i,j) = tmpMean/this->NbrWeights;
            }else{
                OutputData(i,j) = 0;
            }
        }
    }
}

boost::numeric::ublas::matrix<float> imageProc::getInputData(std::string InputFile){

    int Index = 0;
    float * InputData;

    GDALDataset * SrcDataset = (GDALDataset *) GDALOpen( InputFile.c_str(), GA_ReadOnly );
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

    GDALClose(SrcDataset);

    return InputDataBoost;
}

boost::numeric::ublas::matrix<float> imageProc::initOutputData(int NbrRows, int NbrColumns, int KernelMaskSize){

    boost::numeric::ublas::matrix<float> OutputDataBoost;
    OutputDataBoost.resize(NbrRows - 2*KernelMaskSize, NbrColumns - 2*KernelMaskSize);
    for(size_t i = 0; i < OutputDataBoost.size1(); i++){
        for(size_t j = 0; j < OutputDataBoost.size2(); j++){
            OutputDataBoost(i,j) = 0.;
        }
    }
    return OutputDataBoost;

}

void imageProc::filterImg(bool GaborYesNo){

    boost::numeric::ublas::matrix<float> InputDataBoost = this->getInputData(this->InputFileName);
    boost::numeric::ublas::matrix<float> MaskDataBoost = this->getInputData(this->MaskFileName);

    /* * * * * * * * *  Create and initialize output Array with zeros * * * * * * * * */

    int Index = 0;
    int NbrRows = InputDataBoost.size1();
    int NbrColumns = InputDataBoost.size2();

    float * OutputData;
    OutputData = (float *) CPLMalloc(sizeof(float) * (NbrRows-2*this->BorderSize) * (NbrColumns-2*this->BorderSize));
    for(int i = 0; i < (NbrRows - 2 * this->BorderSize) * (NbrColumns - 2 * this->BorderSize); i++){
        OutputData[i] = 0.;
    }

    boost::numeric::ublas::matrix<float> OutputDataBoost;
    OutputDataBoost = this->initOutputData(NbrRows, NbrColumns, this->BorderSize);

    /* * * * * * * * * * * * * * computing tpi Image * * * * * * * * * * * * */
    if(GaborYesNo == false){
        this->computeTPI(InputDataBoost, MaskDataBoost, OutputDataBoost);
    }else{
        this->computeGabor(InputDataBoost, MaskDataBoost, OutputDataBoost);
    }


    /* * * * * * * * * * * * * * save tpi Image * * * * * * * * * * * * */
    for(size_t i = 0; i < OutputDataBoost.size1(); i++){
        for(size_t j = 0; j < OutputDataBoost.size2(); j++){
            Index = j + i*OutputDataBoost.size2();
            OutputData[Index] = OutputDataBoost(i, j);
        }
    }

    this->saveImg(OutputData);


}


imageProc::imageProc(std::string InputFileName, std::string MaskFileName, std::string OutputFileName)
{
    GDALAllRegister();

    this->InputFileName = InputFileName;
    this->MaskFileName = MaskFileName;
    this->OutputFileName = OutputFileName;
    this->NbrWeights = 0;

}

imageProc::~imageProc(){

}
