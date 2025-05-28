#pragma once
#include <cmath>
#include <string>
#include <limits>
#include <optional>
#include <stdexcept>
#include "helper/roots.h"
#include <iostream>

/**
 * Very senstive to floating point errors and it can be the case that the merit function minimum is closest to a point on the complex plane even when you only care about real answers.
 * 
 * 
* Estimates the weight W as Y=WX using weighted total least squares, where Y and X are lists of measurements recursively added; it also weights both Y and X separately.
*
* It can only work with Single Input Single Output data. 
* 
* @warning This method can fail if the variance is significantly larger in magnitude than the measurement value due to floating-point precision limitations.
*
* Gregory L. Plett,
* Recursive approximate weighted total least squares estimation of battery cell total capacity,
* Journal of Power Sources,
* Volume 196, Issue 4,
* 2011,
* Pages 2319-2331,
* ISSN 0378-7753,
* https://doi.org/10.1016/j.jpowsour.2010.09.048
*/
class DualVarianceWeightedTotalLeastSquares {
    public:
        /**
         * @brief Constructor for DualVarianceWeightedTotalLeastSquares
         * 
         * @param nominalValue Initial estimate of the weight
         * @param forgettingFactor Factor to reduce influence of older measurements (0 < f <= 1)
         * @param initialXVariance Variance of a hypothetical (imaginary) measurement of x when x = 1 and y = nominalValue.
         *                  This represents the initial uncertainty in the relationship between x and y.
         * @param initialYVariance Variance of a hypothetical (imaginary) measurement of y when x = 1 and y = nominalValue.
         *                  This represents the initial uncertainty in the relationship between x and y.
         * @param varianceRatio The average relative uncertainty between x and y values is used to improve convergence; it only needs to be an order of magnitude value.  
         *                      By default this class uses the ratio of the first x and y variance found.
         */
        DualVarianceWeightedTotalLeastSquares(
            double nominalValue=0.0, double forgettingFactor=1.0, 
            double initialXVariance=100.0, double initialYVariance=100.0,
            double varianceRatio=-1
        );

        /**
         * @brief Update with a new measurement
         * 
         * @param x measurement for first variable
         * @param y measurement for second variable
         * @param xVariance Variance (uncertainty) of the x measurement (must be more than 0)
         * @param yVariance Variance (uncertainty) of the y measurement (must be more than 0)
         */
        void update(double x, double y, double xVariance, double yVariance);

         /**
         * @brief Get the current variance of the weight estimate
         * 
         * @return Estimated variance of the weight
         */
        double getVariance();

        /**
         * @brief Get the current estimate
         * 
         * @return Current estimate
         */
        double getEstimate();

    private:
        double forgettingFactor;
        double c1;
        double c2;
        double c3;
        double c4;
        double c5;
        double c6;
        double varianceRatio;
        bool hasVarianceRatio;

        /**
         * @brief Get the value of the merit function at a certain estimate
         * 
         * @return value of the merit function
         */
        double getEstimateMerit(double estimate); 
        
        /**
         * @brief Get the current estimate without correction for varianceRatio
         * 
         * @return Current estimate
         */
        double getEstimateUncorrected();
        
};