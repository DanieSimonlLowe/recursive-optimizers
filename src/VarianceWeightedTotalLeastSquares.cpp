#include "VarianceWeightedTotalLeastSquares.h"


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
        throw std::invalid_argument( "Initial Variance must grater then 0 got " + std::to_string(yVariance) );
    }

    // You can't get this yVariance
    this->c1 = 1 / yVariance;
    this->c2 = nominalValue / yVariance;
    this->c3 = (nominalValue * nominalValue) / yVariance;
}


void VarianceWeightedTotalLeastSquares::update(double x, double y, double yVariance) {
    // Don't check input because it would massivly slow down this.
    this->c1 = this->forgettingFactor * this->c1 + x * x / yVariance;
    this->c2 = this->forgettingFactor * this->c2 + x * y / yVariance;
    this->c3 = this->forgettingFactor * this->c3 + y * y / yVariance;
    return;
}


double VarianceWeightedTotalLeastSquares:: getEstimate() {
    if (this->varianceRatioSquared == 0 || this->c2 == 0) {
        return 0.0;
    }
    
    double top_left = -this->c1 + this->varianceRatioSquared * this->c3;
    double top_right_inner = (this->c1 - this->varianceRatioSquared * this->c3);
    double top_right = std::sqrt(top_right_inner * top_right_inner + 4 * this->varianceRatioSquared * this->c2 * this->c2);

    return (top_left + top_right) / (2 * this->varianceRatioSquared * this->c2);
}


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
