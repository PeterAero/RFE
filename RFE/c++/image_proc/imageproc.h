#ifndef IMAGEPROC_H
#define IMAGEPROC_H
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <iostream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>



class imageProc
{
public:
    imageProc(std::string InputFileName, std::string MaskFileName, std::string OutputFileName);
    ~imageProc();

    void filterImg(bool GaborYesNo);

    void initKernel(int size1, int size2);
    void setCircleKernel(float Radius);
    void setDonutKernel(float innerRadius, float outerRadius);
    void setGaborKernel(float frequency, float theta, float bandwidth, float sigma_x, float sigma_y,
                                   float n_stds, float offset);


private:
    boost::numeric::ublas::matrix<float> getInputData(std::string InputFile);
    boost::numeric::ublas::matrix<float> initOutputData(int NbrRows, int NbrColumns, int KernelMaskSize);

    float sigma_prefactor(float bandwidth);

    void saveImg(float * ImgData);
    void printCurrentKernel();

    void computeTPI(const boost::numeric::ublas::matrix<float>& InputData,
                    const boost::numeric::ublas::matrix<float>& MaskData,
                    boost::numeric::ublas::matrix<float>& OutputData);

    void computeGabor(const boost::numeric::ublas::matrix<float>& InputData,
                    const boost::numeric::ublas::matrix<float>& MaskData,
                    boost::numeric::ublas::matrix<float>& OutputData);

    float NbrWeights;
    int BorderSize;
    boost::numeric::ublas::matrix<float> FilterKernel;
    std::string InputFileName;
    std::string MaskFileName;
    std::string OutputFileName;

};

#endif // IMAGEPROC_H
