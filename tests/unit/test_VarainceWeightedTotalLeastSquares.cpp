#include <gtest/gtest.h>
#include <VarianceWeightedTotalLeastSquares.h>
#include <tuple>  

TEST(VWTLSUnitTest, InitialEstimateIsZero) {
    VarianceWeightedTotalLeastSquares estimator = VarianceWeightedTotalLeastSquares();
    EXPECT_NEAR(estimator.getEstimate(), 0.0, 1e-8);
}

TEST(VWTLSUnitTest, InitialVarianceIsOne) {
    VarianceWeightedTotalLeastSquares estimator = VarianceWeightedTotalLeastSquares();
    EXPECT_NEAR(estimator.getVariance(), 1.0, 1e-8);
}

TEST(VWTLSUnitTest, InvalidVarianceRatio) {
    EXPECT_THROW(VarianceWeightedTotalLeastSquares(0.0,0.0), std::invalid_argument);
}

TEST(VWTLSUnitTest, LowForgettingFactor) {
    EXPECT_THROW(VarianceWeightedTotalLeastSquares(0.0,1.0,0.0), std::invalid_argument);
}

TEST(VWTLSUnitTest, HighForgettingFactor) {
    EXPECT_THROW(VarianceWeightedTotalLeastSquares(0.0,1.0,1.001), std::invalid_argument);
}

TEST(VWTLSUnitTest, InvalidInitialVariance) {
    EXPECT_THROW(VarianceWeightedTotalLeastSquares(0.0,1.0,1.0,0.0), std::invalid_argument);
}

class VWTLSEstimateParamTest : public ::testing::TestWithParam<std::tuple<double, double, double>> {
    protected:
        double initialEstimate;
        double variance;
        double varianceRatio;
    
        void SetUp() override {
            std::tie(initialEstimate, variance, varianceRatio) = GetParam();
        }
};
    
TEST_P(VWTLSEstimateParamTest, InitialEstimateTest) {
    VarianceWeightedTotalLeastSquares estimator(initialEstimate, varianceRatio, 1.0, variance);
    EXPECT_NEAR(estimator.getEstimate(), initialEstimate, 1e-8);
}

TEST_P(VWTLSEstimateParamTest, PositiveVarianceTest) {
    VarianceWeightedTotalLeastSquares estimator(initialEstimate, varianceRatio, 1.0, variance);
    EXPECT_GT(estimator.getVariance(), 0);
}


INSTANTIATE_TEST_SUITE_P(
    VWTLSEstimateParamTests,
    VWTLSEstimateParamTest,
    ::testing::Values(
        std::make_tuple(0.0, 1.0,1.0),
        std::make_tuple(5.0, 1.0,2.0),
        std::make_tuple(-3.0, 1.0,0.5),
        std::make_tuple(2.5, 1.0,1.54),
        std::make_tuple(1e3, 1.0,104.2),
        std::make_tuple(1e-4, 1e-4,23.1),
        std::make_tuple(205, 203.0,2.05),
        std::make_tuple(2.5e5, 1e4,0.001),
        std::make_tuple(2.5, 1e-4, 102),
        std::make_tuple(7.25786, 2.15085, 1.0),
        std::make_tuple(9.02267, 4.07965, 1.0)
    )
);


class VWTLSVarianceParamTest : public ::testing::TestWithParam<std::tuple<double, double, double, double>> {
    protected:
        double initialEstimate;
        double variance1;
        double variance2;
        double varianceRatio;
       
        void SetUp() override {
            std::tie(initialEstimate, variance1, variance2, varianceRatio) = GetParam();
        }
};

TEST_P(VWTLSVarianceParamTest, InitialVarianceTest) {
    VarianceWeightedTotalLeastSquares estimator1(initialEstimate, varianceRatio, 1.0, variance1);
    VarianceWeightedTotalLeastSquares estimator2(initialEstimate, varianceRatio, 1.0, variance2);
    if (variance1 > variance2) {
        EXPECT_GT(estimator1.getVariance(), estimator2.getVariance());
    } else {
        EXPECT_GT(estimator2.getVariance(), estimator1.getVariance());
    }
    
}

INSTANTIATE_TEST_SUITE_P(
    VWTLSVarianceParamTests,
    VWTLSVarianceParamTest,
    ::testing::Values(
        std::make_tuple(0.0, 1.0, 1.1, 1.0),
        std::make_tuple(5.0, 1.0, 1.1, 2.0),
        std::make_tuple(-3.0, 1.0, 0.2, 0.5),
        std::make_tuple(2.5, 1.0, 2.0, 1.654),
        std::make_tuple(1e3, 1.0, 10, 1234),
        std::make_tuple(1e-4, 1e-4, 1e-5, 0.512),
        std::make_tuple(205, 203.0, 204.0, 12.21),
        std::make_tuple(2.5e5, 1e4, 1e3, 1e3),
        std::make_tuple(2.5, 1e-4, 1e-3, 1.42),
        std::make_tuple(9.27450, 8.18914, 3.13302, 0.00123),
        std::make_tuple(2.78746, 1.70555, 0.66330, 450.541)
    )
);