#include "VarianceWeightedTotalLeastSquares.h"


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


/**
 * @brief Constructor for VarianceWeightedTotalLeastSquares
 * 
 * @param nominalValue Initial estimate of the capacity
 * @param varianceRatio Ratio of output measurement variance to input measurement variance of x over y
 * @param forgettingFactor Factor to reduce influence of older measurements (0 < f <= 1)
 * @param yVariance Variance of a hypothetical (imaginary) measurement of y when x = 1 and y = nominalValue.
 *                  This represents the initial uncertainty in the relationship between x and y.
 */
VarianceWeightedTotalLeastSquares::VarianceWeightedTotalLeastSquares (
    double nominalValue, double varianceRatio,
    double forgettingFactor, double yVariance
) {
    if (forgettingFactor > 1 || forgettingFactor <= 0) {
        throw std::invalid_argument( "Forgetting Factor must be in the range 0 to 1 (exluding zero) got " + std::to_string(forgettingFactor) );
    }

    this->forgettingFactor = forgettingFactor;

    if (varianceRatio <= 0) {
        throw std::invalid_argument( "Variance Ratio must grater then 0 got " + std::to_string(varianceRatio) );
    }

    this->varianceRatioSquared = varianceRatio * varianceRatio;


    if (yVariance <= 0) {
        throw std::invalid_argument( "Initial Variance must grater then 0 got " + std::to_string(varianceRatio) );
    }

    // You can't get this yVariance
    this->c1 = 1 / yVariance;
    this->c2 = nominalValue / yVariance;
    this->c3 = (nominalValue * nominalValue) / yVariance;
}

/**
 * @brief Constructor for VarianceWeightedTotalLeastSquares
 * 
 * @param x mesurement for first variabile
 * @param y mesurement for second variabile
 * @param yVariance Variance (uncertainty) of the x measurement (must be more then 0)
 */
void VarianceWeightedTotalLeastSquares::update(double x, double y, double yVariance) {
    // Don't check input because it would massivly slow down this.
    this->c1 = this->forgettingFactor * this->c1 + x * x / yVariance;
    this->c2 = this->forgettingFactor * this->c2 + x * y / yVariance;
    this->c3 = this->forgettingFactor * this->c3 + y * y / yVariance;
    return;
}

/**
 * @brief Get the current capacity estimate
 * 
 * @return Current estimate of battery capacity
 */
double VarianceWeightedTotalLeastSquares:: getEstimate() {
    if (this->varianceRatioSquared == 0 || this->c2 == 0) {
        return 0.0;
    }
    
    double top_left = -this->c1 + this->varianceRatioSquared * this->c3;
    double top_right_inner = (this->c1 - this->varianceRatioSquared * this->c3);
    double top_right = std::sqrt(top_right_inner * top_right_inner + 4 * this->varianceRatioSquared * this->c2 * this->c2);

    return (top_left + top_right) / (2 * this->varianceRatioSquared * this->c2);
}

/**
 * @brief Get the current variance of the weight estimate
 * 
 * @return Estimated variance of the weight
 */
double VarianceWeightedTotalLeastSquares:: getVariance() {
    // TODO rewrite this based on page 6/2034 of http://mocha-java.uccs.edu/dossier/RESEARCH/2011jps-.pdf
    double estimate = this->getEstimate();
    
    double bottom = (estimate * estimate * this->varianceRatioSquared + 1);

    double top = (-4.0 * this->varianceRatioSquared * this->varianceRatioSquared * this->c2) * estimate * estimate * estimate
           + 6.0 * this->varianceRatioSquared * this->varianceRatioSquared * this->c3 * estimate * estimate
           + (-6.0 * this->c1 + 12.0 * this->c2) * this->varianceRatioSquared * estimate
           + 2.0 * (this->c1 - this->varianceRatioSquared * this->c3);

    double hessian = top / (bottom * bottom * bottom);
    
    return 2 / hessian;
}
