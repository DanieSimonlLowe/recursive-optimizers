#include "DualVarianceWeightedTotalLeastSquares.h"




DualVarianceWeightedTotalLeastSquares::DualVarianceWeightedTotalLeastSquares(double nominalValue, double forgettingFactor, 
            double initialXVariance, double initialYVariance, double varianceRatio) {
    if (forgettingFactor > 1 || forgettingFactor <= 0) {
        throw std::invalid_argument( "Forgetting Factor must be in the range 0 to 1 (exluding zero) got " + std::to_string(forgettingFactor) );
    }

    
    this->forgettingFactor = forgettingFactor;

    if (initialXVariance <= 0) {
        throw std::invalid_argument( "Initial X Variance must grater then 0 got " + std::to_string(initialXVariance) );
    }

    if (initialYVariance <= 0) {
        throw std::invalid_argument( "Initial Y Variance must grater then 0 got " + std::to_string(initialYVariance) );
    }

    if (varianceRatio != -1) {
        if (varianceRatio <= 0) {
            throw std::invalid_argument( "Variance Ratio must grater then 0 got " + std::to_string(varianceRatio) );
        }
    
        this->varianceRatio = varianceRatio;
        this->hasVarianceRatio = true;

        nominalValue = this->varianceRatio * nominalValue;
    }

    this->c1 = 1 / initialYVariance;
    this->c2 = nominalValue / initialYVariance;
    this->c3 = nominalValue * nominalValue / initialYVariance;

    this->c4 = 1 / initialXVariance;
    this->c5 = nominalValue / initialXVariance;
    this->c6 = nominalValue * nominalValue / initialXVariance;
    
    
}



void DualVarianceWeightedTotalLeastSquares::update(double x, double y, double xVariance, double yVariance) {
    if (!this->hasVarianceRatio) {
        // Asumes a value 
        this->varianceRatio = std::sqrt(xVariance) / std::sqrt(yVariance);
        this->hasVarianceRatio = true;

        this->c1 /= this->varianceRatio * this->varianceRatio;
        this->c2 /= this->varianceRatio;
        // c3 is not affected because they cancle out
        // c4 is not affected because it has no y factor
        this->c5 *= this->varianceRatio;
        this->c6 *= this->varianceRatio * this->varianceRatio;
    }
    

    double correctedY = y * this->varianceRatio;
    double yBottom = yVariance * this->varianceRatio * this->varianceRatio;

    this->c1 = this->forgettingFactor * this->c1 + x * x / yBottom;
    this->c2 = this->forgettingFactor * this->c2 + x * correctedY / yBottom;
    this->c3 = this->forgettingFactor * this->c3 + correctedY * correctedY / yBottom;

    this->c4 = this->forgettingFactor * this->c4 + x * x /  xVariance;
    this->c5 = this->forgettingFactor * this->c5 + x * correctedY /  xVariance;
    this->c6 = this->forgettingFactor * this->c6 + correctedY * correctedY /  xVariance;


}


double DualVarianceWeightedTotalLeastSquares::getEstimateMerit(double estimate) {
    
    double estimateSq = estimate * estimate;
    double top = this->c4 * estimateSq * estimateSq + 
                -2 * this->c5 * estimateSq * estimate +
                (this->c1 + this->c6) * estimateSq +
                -2 * this->c2 * estimate + 
                this->c3;
    
    double bottom = estimateSq + 1;

    return top / (bottom * bottom);
}



double DualVarianceWeightedTotalLeastSquares::getEstimateUncorrected() {
    double a = this->c5;
    double b = 2 * this->c4 - this->c1 - this->c6;
    double c = 3 * this->c2 - 3 * this->c5;
    double d = this->c1 - 2 * this->c3 + this->c6;
    double e = -this->c2;
    
    std::vector<double> roots = calculate_real_roots(a,b,c,d,e);

    int bestRootPos = -1;
    double bestMerit;
    for (int i = 0; i < roots.size(); ++i) {
        double end = a * roots[i] * roots[i] * roots[i] * roots[i] + 
                     b * roots[i] * roots[i] * roots[i] +
                     c * roots[i] * roots[i] +
                     d * roots[i] + e;
        
        if (bestRootPos == -1) {
            bestRootPos = i;
            bestMerit = this->getEstimateMerit(roots[i]);
        } else {
            double merit = this->getEstimateMerit(roots[i]);
            if (merit <= bestMerit) {
                bestRootPos = i;
                bestMerit = merit;
            }
        }
    }

    if (bestRootPos == -1) {
        throw std::domain_error("All roots are complex.");
    }


    return roots[bestRootPos];
}



double DualVarianceWeightedTotalLeastSquares::getVariance() {
    double estimate = this->getEstimateUncorrected();
    double estimateSq = estimate * estimate;

    double top = -2 * this->c5 * estimateSq * estimateSq * estimate + 
                (3 * this->c3 - 6 * this->c4 + 3 * this->c6) * estimateSq * estimateSq +
                (-12 * this->c2 + 16 * this->c5) * estimateSq * estimate +
                (-8 * this->c1 + 10 * this->c3 + 6 * this->c4 - 8 * this->c6) * estimateSq +
                (12 * this->c2 - 6 * this->c5) * estimate +
                this->c1 - 2 * this->c3 + this->c6;

    double bottom = estimateSq + 1;

    double hessian = 2 * top / (bottom * bottom * bottom * bottom);
    //hessian = hessian / (this->varianceRatio * this->varianceRatio); // Correcting the hassian by the varianceRatio
    return 2.0 * this->varianceRatio * this->varianceRatio / hessian;
}

double DualVarianceWeightedTotalLeastSquares::getEstimate() {
    return this->getEstimateUncorrected() / this->varianceRatio;
}