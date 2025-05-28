#include <gtest/gtest.h>
#include <DualVarianceWeightedTotalLeastSquares.h>
#include <tuple>  

TEST(DVWTLSUnitTest, InitialEstimateIsZero) {
    DualVarianceWeightedTotalLeastSquares estimator = DualVarianceWeightedTotalLeastSquares();
    EXPECT_NEAR(estimator.getEstimate(), 0.0, 1e-8);
}

TEST(DVWTLSUnitTest, LowForgettingFactor) {
    EXPECT_THROW(DualVarianceWeightedTotalLeastSquares(0.0,0.0), std::invalid_argument);
}

TEST(DVWTLSUnitTest, HighForgettingFactor) {
    EXPECT_THROW(DualVarianceWeightedTotalLeastSquares(0.0,1.001), std::invalid_argument);
}

TEST(DVWTLSUnitTest, InvalidInitialXVariance) {
    EXPECT_THROW(DualVarianceWeightedTotalLeastSquares(0.0,1.0,0.0), std::invalid_argument);
}

TEST(DVWTLSUnitTest, InvalidInitialYVariance) {
    EXPECT_THROW(DualVarianceWeightedTotalLeastSquares(0.0,1.0,1.0,0.0), std::invalid_argument);
}

TEST(DVWTLSUnitTest, InvalidVarianceRatio) {
    EXPECT_THROW(DualVarianceWeightedTotalLeastSquares(0.0,1.0,1.0,1.0,0.0), std::invalid_argument);
}

TEST(DVWTLSUnitTest, ValidVarianceRatio) {
    EXPECT_NO_THROW(DualVarianceWeightedTotalLeastSquares(0.0,1.0,1.0,1.0,1e6));
}

TEST(DVWTLSUnitTest, ValidCreation) {
    EXPECT_NO_THROW(DualVarianceWeightedTotalLeastSquares());
}

TEST(DVWTLSUnitTest, ValidCreationWithInitalEstimate) {
    EXPECT_NO_THROW(DualVarianceWeightedTotalLeastSquares(11.0));
}

TEST(DVWTLSUnitTest, UpdateToOne) {
    DualVarianceWeightedTotalLeastSquares estimator = DualVarianceWeightedTotalLeastSquares();
    
    estimator.update(1,1,1e-4,1e-4);
    EXPECT_NEAR(estimator.getEstimate(), 1.0, 1e-4);
}

TEST(DVWTLSUnitTest, UpdateToTwo) {
    DualVarianceWeightedTotalLeastSquares estimator = DualVarianceWeightedTotalLeastSquares();   
    estimator.update(1,2,1e-4,1e-4);
    EXPECT_NEAR(estimator.getEstimate(), 2.0, 1e-4);
}

