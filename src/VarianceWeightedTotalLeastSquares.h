#pragma once
#include <cmath>
#include <string>
#include <stdexcept>

/*
Estmates the weight W as Y=WX by doing weighted total least sqears, where Y and X are a list of mesurements recusivly.
Can only work with SISO data. 

Can fail if the varaince is proptonaly larger then the mesurement value because of floating point presion problems.

Gregory L. Plett,
Recursive approximate weighted total least squares estimation of battery cell total capacity,
Journal of Power Sources,
Volume 196, Issue 4,
2011,
Pages 2319-2331,
ISSN 0378-7753,
https://doi.org/10.1016/j.jpowsour.2010.09.048.
*/

class VarianceWeightedTotalLeastSquares {
    public:
        /**
         * @brief Constructor for VarianceWeightedTotalLeastSquares
         * 
         * @param nominalValue Initial estimate of the capacity
         * @param varianceRatio Ratio of output measurement variance to input measurement variance of x over y
         * @param forgettingFactor Factor to reduce influence of older measurements (0 < f <= 1)
         * @param yVariance Variance of a hypothetical (imaginary) measurement of y when x = 1 and y = nominalValue.
         *                  This represents the initial uncertainty in the relationship between x and y.
         */
        VarianceWeightedTotalLeastSquares(
            double nominalValue=0.0, double varainceRatio=1.0,
            double forgettingFactor=1.0, double initialVariance=1.0
        );

        /**
         * @brief Update with a new measurement
         * 
         * @param x mesurement for first variabile
         * @param y mesurement for second variabile
         * @param yVariance Variance (uncertainty) of the y measurement (must be more then 0)
         */
        void update(double x, double y, double yVariance);

        
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
        double varianceRatioSquared; // because it is allways used as sqeared
        double c1;
        double c2;
        double c3;
        
};