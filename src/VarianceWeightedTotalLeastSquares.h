#pragma once
#include <cmath>
#include <string>
#include <stdexcept>

class VarianceWeightedTotalLeastSquares {
    public:
        VarianceWeightedTotalLeastSquares(
            double nominalValue=0.0, double varainceRatio=1.0,
            double forgettingFactor=1.0, double initialVariance=1.0
        );

        void update(double x, double y, double yVariance);

        double getVariance();


        double getEstimate();

    private:
        double forgettingFactor;
        double varianceRatioSquared; // because it is allways used as sqeared
        double c1;
        double c2;
        double c3;
        
};