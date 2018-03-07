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
    template <typename T> void createCircleMask(T Radius);
    template <typename T> void createDonutMask(T innerRadius, T outerRadius);
    std::string FileName;

private:
    void saveImg(GDALDataset * poDataset, float * ImgData, int& BorderSize);
    void printCurrentKernel();

    void computeTPI(const boost::numeric::ublas::matrix<float>& InputData,
                    boost::numeric::ublas::matrix<float>& OutputData,
                    int& KernelMaskSize);

    void initMatrix(int size1, int size2);
    float NbrWeights;
    boost::numeric::ublas::matrix<float> FilterKernel;

    std::string OutputFileName;

};

#endif // IMAGEPROC_H
