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
    imageProc(std::string FileName, std::string OutputFileName);
    ~imageProc();

    void filterImg(const double& InnerRadius, const double& OuterRadius, const bool& DonutYesNo);
    void setCircleKernel(float Radius);
    void setDonutKernel(float innerRadius, float outerRadius);
    void setGaborKernel(float frequency, float theta, float bandwidth, float sigma_x, float sigma_y,
                                   float n_stds, float offset);
    void initKernel(int size1, int size2);
    std::string FileName;

private:
    boost::numeric::ublas::matrix<float> getInputData();
    boost::numeric::ublas::matrix<float> initOutputData(int NbrRows, int NbrColumns, int KernelMaskSize);

    float sigma_prefactor(float bandwidth);

    void saveImg(float * ImgData, int& BorderSize);
    void printCurrentKernel();

    void computeTPI(const boost::numeric::ublas::matrix<float>& InputData,
                    boost::numeric::ublas::matrix<float>& OutputData,
                    int& KernelMaskSize);


    float NbrWeights;
    int BorderSize;
    boost::numeric::ublas::matrix<float> FilterKernel;

    std::string OutputFileName;

};

#endif // IMAGEPROC_H
