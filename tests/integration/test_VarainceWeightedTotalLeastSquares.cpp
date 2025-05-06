#include <gtest/gtest.h>
#include <VarianceWeightedTotalLeastSquares.h>
#include <tuple>
#include <random>


class VWTLSParamTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, long unsigned int>> {
    protected:

        double estimate;

        double varianceRatio;
        
        double min_x;
        double max_x;
        double max_y_noise;

        long unsigned int seed;

    
        void SetUp() override {
            std::tie(estimate, varianceRatio, min_x, max_x, max_y_noise, seed) = GetParam();
        }
};
    
TEST_P(VWTLSParamTest, VWTLSEstimateIntegrationTest) {
    std::mt19937 gen{seed};

    std::normal_distribution<double> std_generator(0, 1);
    std::uniform_real_distribution<> x_generator(min_x, max_x);
    std::uniform_real_distribution<> yStd_generator(0, max_y_noise);

    VarianceWeightedTotalLeastSquares estimator(0.0, varianceRatio);

    for (int i = 0; i < 500; i++) {
        
        double real_x = x_generator(gen);
        double real_y = real_x * estimate;

        double std_y = yStd_generator(gen);
        double std_x = std_y * varianceRatio;

        double x = real_x + std_x * std_generator(gen);
        double y = real_y + std_y * std_generator(gen);

        estimator.update(x, y, std_y * std_y);
    }

    EXPECT_NEAR(estimator.getEstimate(), estimate, 5e-3 * abs(estimate));
}

TEST_P(VWTLSParamTest, VWTLSVarianceIntegrationTest) {
    std::mt19937 gen{seed};

    std::normal_distribution<double> std_generator(0, 1);
    std::uniform_real_distribution<> x_generator(min_x, max_x);
    std::uniform_real_distribution<> yStd_generator(0, max_y_noise);

    VarianceWeightedTotalLeastSquares estimator(0.0, varianceRatio);

    for (int i = 0; i < 500; i++) {
        
        double real_x = x_generator(gen);
        double real_y = real_x * estimate;

        double std_y = yStd_generator(gen);
        double std_x = std_y * varianceRatio;

        double x = real_x + std_x * std_generator(gen);
        double y = real_y + std_y * std_generator(gen);

        estimator.update(x, y, std_y * std_y);

        if (i > 10) {
            EXPECT_NEAR(estimator.getEstimate(), estimate, std::sqrt(estimator.getVariance()) * 3);
        }
    }
}


INSTANTIATE_TEST_SUITE_P(
    VWTLSParamTests,
    VWTLSParamTest,
    ::testing::Values(
        std::make_tuple(2.0, 1.0, -1.0, 1.0, 0.1, 1),
        std::make_tuple(8.0, 0.4, -10.0, 19.0, 0.1, 2),
        std::make_tuple(58.0, 1.4, 3.0, 1e4, 0.5, 3),
        std::make_tuple(-1.0, 4.0, 13.0, 14.0, 0.15, 4),
        std::make_tuple(-1e4, 1e-2, -2.0, 0.0, 5e2, 5),

        std::make_tuple(2.0, 1.0, -1.0, 1.0, 0.1, 10),
        std::make_tuple(8.0, 0.4, -10.0, 19.0, 0.1, 20),
        std::make_tuple(58.0, 1.4, 3.0, 1e4, 0.5, 30),
        std::make_tuple(-1.0, 4.0, 13.0, 14.0, 0.15, 40),
        std::make_tuple(-1e4, 1e-2, -2.0, 0.0, 5e2, 50),

        std::make_tuple(2.0, 1.0, -1.0, 1.0, 0.1, 100),
        std::make_tuple(8.0, 0.4, -10.0, 19.0, 0.1, 200),
        std::make_tuple(58.0, 1.4, 3.0, 1e4, 0.5, 300),
        std::make_tuple(-1.0, 4.0, 13.0, 14.0, 0.15, 400),
        std::make_tuple(-1e4, 1e-2, -2.0, 0.0, 5e2, 500)
    )
);


class VWTLSFFParamTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, double, double, long unsigned int>> {
    protected:

        double estimate1;
        double estimate2;


        double varianceRatio;
        double forgettingFactor;

        double min_x;
        double max_x;
        double max_y_noise;

        long unsigned int seed;

    
        void SetUp() override {
            std::tie(estimate1, estimate2, varianceRatio, forgettingFactor, min_x, max_x, max_y_noise, seed) = GetParam();
        }
};

TEST_P(VWTLSFFParamTest, VWTLSForgettingFactorIntegrationTest) {
    std::mt19937 gen{seed};

    std::normal_distribution<double> std_generator(0, 1);
    std::uniform_real_distribution<> x_generator(min_x, max_x);
    std::uniform_real_distribution<> yStd_generator(0, max_y_noise);

    VarianceWeightedTotalLeastSquares estimator(0.0, varianceRatio, forgettingFactor);

    for (int i = 0; i < 1000; i++) {
        double estimate = (i < 500) ? estimate1 : estimate2;
        double real_x = x_generator(gen);
        double real_y = real_x * estimate;

        double std_y = yStd_generator(gen);
        double std_x = std_y * varianceRatio;

        double x = real_x + std_x * std_generator(gen);
        double y = real_y + std_y * std_generator(gen);

        estimator.update(x, y, std_y * std_y);

    }

    EXPECT_NEAR(estimator.getEstimate(), estimate2, 1e-2 * abs(estimate2));
}

TEST_P(VWTLSFFParamTest, VWTLSForgettingFactorVarianceIntegrationTest) {
    std::mt19937 gen{seed};

    std::normal_distribution<double> std_generator(0, 1);
    std::uniform_real_distribution<> x_generator(min_x, max_x);
    std::uniform_real_distribution<> yStd_generator(0, max_y_noise);

    VarianceWeightedTotalLeastSquares estimator(0.0, varianceRatio, forgettingFactor);

    for (int i = 0; i < 1000; i++) {
        double estimate = (i < 500) ? estimate1 : estimate2;
        double real_x = x_generator(gen);
        double real_y = real_x * estimate;

        double std_y = yStd_generator(gen);
        double std_x = std_y * varianceRatio;

        double x = real_x + std_x * std_generator(gen);
        double y = real_y + std_y * std_generator(gen);

        estimator.update(x, y, std_y * std_y);

    }
    EXPECT_NEAR(estimator.getEstimate(), estimate2, std::sqrt(estimator.getVariance()) * 3);
}


INSTANTIATE_TEST_SUITE_P(
    VWTLSFFParamTests,
    VWTLSFFParamTest,
    ::testing::Values(
        std::make_tuple(2.0, -2.0, 1.0, 0.97, -1.0, 1.0, 0.1, 1),
        std::make_tuple(8.0, 1.0, 0.4, 0.98, -10.0, 19.0, 1, 2),
        std::make_tuple(58.0, 103.2, 1.4, 0.97, 3.0, 1e4, 0.5, 3),
        std::make_tuple(-2.0, -0.1, 4.0, 0.975, 13.0, 14.0, 0.5, 4),
        std::make_tuple(-1e4, 52.0, 1e-2, 0.96, -2.0, 0.0, 5, 5),

        std::make_tuple(2.0, -2.0, 1.0, 0.97, -1.0, 1.0, 0.1, 10),
        std::make_tuple(8.0, 1.0, 0.4, 0.98, -10.0, 19.0, 1, 20),
        std::make_tuple(58.0, 103.2, 1.4, 0.97, 3.0, 1e4, 0.5, 30),
        std::make_tuple(-2.0, -0.1, 4.0, 0.975, 13.0, 14.0, 0.5, 40),
        std::make_tuple(-1e4, 52.0, 1e-2, 0.96, -2.0, 0.0, 5, 50),

        std::make_tuple(2.0, -2.0, 1.0, 0.97, -1.0, 1.0, 0.1, 100),
        std::make_tuple(8.0, 1.0, 0.4, 0.98, -10.0, 19.0, 1, 200),
        std::make_tuple(58.0, 103.2, 1.4, 0.97, 3.0, 1e4, 0.5, 300),
        std::make_tuple(-2.0, -0.1, 4.0, 0.975, 13.0, 14.0, 0.5, 400),
        std::make_tuple(-1e4, 52.0, 1e-2, 0.96, -2.0, 0.0, 5, 500)
    )
);