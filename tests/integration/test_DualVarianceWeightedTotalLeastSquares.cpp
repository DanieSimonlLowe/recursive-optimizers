#include <gtest/gtest.h>
#include <DualVarianceWeightedTotalLeastSquares.h>
#include <tuple>
#include <random>


class DVWTLSParamTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double, long unsigned int>> {
    protected:

        double estimate;
        
        double min_x;
        double max_x;

        double max_x_noise;
        double max_y_noise;

        long unsigned int seed;

    
        void SetUp() override {
            std::tie(estimate, min_x, max_x, max_x_noise, max_y_noise, seed) = GetParam();
        }
};
    
TEST_P(DVWTLSParamTest, DVWTLSEstimateIntegrationTest) {
    std::mt19937 gen{seed};

    std::normal_distribution<double> std_generator(0, 1);
    std::uniform_real_distribution<> x_generator(min_x, max_x);
    std::uniform_real_distribution<> xStd_generator(1e-5, max_x_noise);
    std::uniform_real_distribution<> yStd_generator(1e-5, max_y_noise);

    DualVarianceWeightedTotalLeastSquares estimator(0.0,1.0, 100, 100, 1.0);

    for (int i = 0; i < 500; i++) {
        
        double real_x = x_generator(gen);
        double real_y = real_x * estimate;

        double std_y = yStd_generator(gen);
        double std_x = xStd_generator(gen);

        double x = real_x + std_x * std_generator(gen);
        double y = real_y + std_y * std_generator(gen);

        estimator.update(x, y, std_x*std_x, std_y*std_y);
    }

    EXPECT_NEAR(estimator.getEstimate(), estimate, 1e-8);
}

INSTANTIATE_TEST_SUITE_P(
    DVWTLSParamTests,
    DVWTLSParamTest,
    ::testing::Values(
        std::make_tuple(1.0, -1.0, 1.0, 1.0, 1.0, 1),
        std::make_tuple(5.0, -100.0, 0.0, 100.0, 1.0, 2),
        std::make_tuple(200.0, 0.025, 2, 1.0, 100.0, 3),
        std::make_tuple(0.0, -5, 1, 0.5, 0.1, 4),
        std::make_tuple(0.01, 100, 110, 5, 0.1, 5),
        std::make_tuple(-3.26, 20, 50, 3, 1, 6),

        std::make_tuple(1.0, -1.0, 1.0, 1.0, 1.0, 10),
        std::make_tuple(5.0, -100.0, 0.0, 100.0, 1.0, 20),
        std::make_tuple(200.0, 0.025, 2, 1.0, 100.0, 30),
        std::make_tuple(0.0, -5, 1, 0.5, 0.1, 40),
        std::make_tuple(0.01, 100, 110, 5, 0.1, 50),
        std::make_tuple(-3.26, 20, 50, 3, 1, 60),

        std::make_tuple(1.0, -1.0, 1.0, 1.0, 1.0, 100),
        std::make_tuple(5.0, -100.0, 0.0, 100.0, 1.0, 200),
        std::make_tuple(200.0, 0.025, 2, 1.0, 100.0, 300),
        std::make_tuple(0.0, -5, 1, 0.5, 0.1, 400),
        std::make_tuple(0.01, 100, 110, 5, 0.1, 500),
        std::make_tuple(-3.26, 20, 50, 3, 1, 600)
    )
);

