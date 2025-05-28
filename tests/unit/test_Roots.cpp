#include <gtest/gtest.h>
#include <helper/roots.h>
#include <tuple>  
#include <random>

bool ContainsCloseValue(const std::vector<double>& vec, double target, double tolerance) {
    for (double value : vec) {
        if (std::fabs(value - target) <= tolerance) {
            return true;
        }
    }
    return false;
}

void expect_double_in(const std::vector<double>& values, double test, double error) {
    bool found = ContainsCloseValue(values, test, error);
    std::ostringstream oss;
    oss << "Expected vector to contain a value close to " << test
        << " within error " << error
        << ", but was: { ";
    for (size_t i = 0; i < values.size(); ++i) {
        oss << values[i];
        if (i + 1 < values.size()) oss << ", ";
    }
    oss << " }";
    EXPECT_TRUE(found) << oss.str();
}


class Poly2DegreeParamTest : public ::testing::TestWithParam<std::tuple<double, double, double,
                                                                        double, double, unsigned int>> {
    protected:
        double val1;
        double val2;
        double val3;
        double out_1;
        double out_2;
        unsigned int out_size;
        

    
        void SetUp() override {
            std::tie(val1,val2,val3, out_1,out_2, out_size) = GetParam();
        }
};

TEST_P(Poly2DegreeParamTest, Poly2DegreeParamSizeTest) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3);
    EXPECT_EQ(out.size(),out_size);
}

TEST_P(Poly2DegreeParamTest, Poly2DegreeParamOut1) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3);
    if (out_size >= 1) {
        expect_double_in(out,out_1,1e-8);
    }
}

TEST_P(Poly2DegreeParamTest, Poly2DegreeParamOut2) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3);
    if (out_size >= 2) {
        expect_double_in(out,out_2,1e-8);
    }
}

TEST_P(Poly2DegreeParamTest, Poly2DegreeParamValidValues) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3);
    for (int i = 0 ; i < out_size; i++) {
        
        double func = out[i] * out[i] * val1 + out[i] * val2 + val3;
        EXPECT_NEAR(func,0,1e-8);
    }
}


INSTANTIATE_TEST_SUITE_P(
    Poly2DegreeParamTests,
    Poly2DegreeParamTest,
    ::testing::Values(
        std::make_tuple(1.0,1.0,1.0, 
                        0.0,0.0,0),
        std::make_tuple(1.0,2.0,1.0, 
                        -1.0,0.0,1),
        std::make_tuple(1.0,5.0,4.0, 
                        -1.0,-4.0,2),
        std::make_tuple(1.0,5.0,4.0, 
                        -1.0,-4.0,2),
        std::make_tuple(0.0,5.0,4.0, 
                        -0.8,0.0,1),
        std::make_tuple(2.0,1.0,-1.0, 
                        -1.0,0.5,2),
        std::make_tuple(-2.0,0.0,8.0, 
                        -2.0,2.0,2)
    )
);


class Poly3DegreeParamTest : public ::testing::TestWithParam<std::tuple<double, double, double, double,
                                                                        double, double, double, unsigned int>> {
    protected:
        double val1;
        double val2;
        double val3;
        double val4;
        double out_1;
        double out_2;
        double out_3;

        unsigned int out_size;
        

    
        void SetUp() override {
            std::tie(val1,val2,val3,val4, out_1,out_2,out_3, out_size) = GetParam();
        }
};

TEST_P(Poly3DegreeParamTest, Poly3DegreeParamSizeTest) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4);
    EXPECT_EQ(out.size(),out_size);
}

TEST_P(Poly3DegreeParamTest, Poly3DegreeParamOut1) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4);
    if (out_size >= 1) {
        expect_double_in(out,out_1,1e-8);
    }
}

TEST_P(Poly3DegreeParamTest, Poly3DegreeParamOut2) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4);
    if (out_size >= 2) {
        expect_double_in(out,out_2,1e-8);
    }
}

TEST_P(Poly3DegreeParamTest, Poly3DegreeParamOut3) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4);
    if (out_size >= 3) {
        expect_double_in(out,out_3,1e-8);
    }
}

TEST_P(Poly3DegreeParamTest, Poly3DegreeParamValidValues) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4);
    for (int i = 0 ; i < out_size; i++) {
        
        double func = out[i] * out[i] * out[i] * val1 + out[i] * out[i] * val2 + out[i] * val3 + val4;
        EXPECT_NEAR(func,0,1e-8);
    }
}


INSTANTIATE_TEST_SUITE_P(
    Poly3DegreeParamTests,
    Poly3DegreeParamTest,
    ::testing::Values(
        std::make_tuple(1.0,1.0,1.0,1.0, 
                        -1.0,0.0,0.0, 1),
        std::make_tuple(1.0,1.0,0.0,0.0, 
                        -1.0,0.0,0.0, 2),
        std::make_tuple(1.0,-1.0,-1.0,1.0, 
                        -1.0,1.0,0.0, 2),
        std::make_tuple(2.0,9.0,3.0,-4.0, 
                        -1.0,-4.0,0.5, 3),
        std::make_tuple(4.0,-3.0,20.0,-15.0, 
                        0.75,0.0,0.0, 1),
        std::make_tuple(0.0,1.0,1.0,1.0, 
                        0.0,0.0,0.0, 0),
        std::make_tuple(1.0,-6.0,11.0,-6.0, 
                        1.0,2.0,3.0, 3),
        std::make_tuple(1.0,0.0,-3.0,2.0, 
                        1.0,-2.0,0.0, 2)
    )
);




class Poly4DegreeParamTest : public ::testing::TestWithParam<std::tuple<double, double, double, double, double,
                                                                        double, double, double, double, unsigned int>> {
    protected:
        double val1;
        double val2;
        double val3;
        double val4;
        double val5;
        double out_1;
        double out_2;
        double out_3;
        double out_4;

        unsigned int out_size;
        

    
        void SetUp() override {
            std::tie(val1,val2,val3,val4, val5, out_1,out_2,out_3, out_4, out_size) = GetParam();
        }
};

TEST_P(Poly4DegreeParamTest, Poly4DegreeParamSizeTest) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4,val5);
    EXPECT_EQ(out.size(),out_size);
}

TEST_P(Poly4DegreeParamTest, Poly4DegreeParamOut1) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4,val5);
    if (out_size >= 1) {
        expect_double_in(out,out_1,1e-8);
    }
}

TEST_P(Poly4DegreeParamTest, Poly4DegreeParamOut2) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4,val5);
    if (out_size >= 2) {
        expect_double_in(out,out_2,1e-8);
    }
}

TEST_P(Poly4DegreeParamTest, Poly4DegreeParamOut3) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4,val5);
    if (out_size >= 3) {
        expect_double_in(out,out_3,1e-8);
    }
}

TEST_P(Poly4DegreeParamTest, Poly4DegreeParamOut4) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4,val5);
    if (out_size >= 4) {
        expect_double_in(out,out_4,1e-8);
    }
}


TEST_P(Poly4DegreeParamTest, Poly4DegreeParamValidValues) {
    std::vector<double> out = calculate_real_roots(val1,val2,val3,val4,val5);
    for (int i = 0 ; i < out.size(); i++) {
        double x = out[i];
        double x2 = x * x;
        
        double func = x2 * x2 * val1 + x2 * x * val2 + x2 * val3 + x * val4 + val5;
        EXPECT_NEAR(func,0, 1e-7 * std::max(std::max(val1,val2),std::max(std::max(val3,val4),val5)) );
    }
}

INSTANTIATE_TEST_SUITE_P(
    Poly4DegreeParamTests,
    Poly4DegreeParamTest,
    ::testing::Values(
        std::make_tuple(1.0,1.0,1.0,1.0,1.0, 
                        0.0,0.0,0.0,0.0, 0),
        std::make_tuple(1.0,10.0,35.0,50.0,24.0,
                        -1.0,-2.0,-3.0,-4.0, 4),
        std::make_tuple(1.0,-6.0,17.0,-24.0,12.0, 
                        1.0,2.0,0.0,0.0, 2),
        std::make_tuple(1.0,0.0,0.0,0.0,-1.0, 
                        1.0,-1.0,0.0,0.0, 2),
        std::make_tuple(1.0,0.0,0.0,0.0,1.0, 
                        0.0,0.0,0.0,0.0, 0),
        std::make_tuple(2.0,-8.0,-10.0,0.0,0.0, 
                        -1.0,5.0,0.0,0.0, 4),
        std::make_tuple(1.0e10,10.0e10,35.0e10,50.0e10,24.0e10,
                        -1.0,-2.0,-3.0,-4.0, 4),
        std::make_tuple(-4.81942e+07,-5.94044e+10,-5.8153e+11,-1.20916e+12,1.93892e+11,
                        0.14944716860361495, -3.2717463095442554, -6.729130492137074, -1222.753323635988, 4),
        std::make_tuple(4.81799e+06,-6.04703e+11,1.82146e+10,6.05613e+11,-6.07634e+09,
                        125509.36279156385, 1.0109525048610521, 0.010031546719519668, -0.9908545010184753, 4)
    )
);


// TEST(Poly4DegreeRandomTest, Poly4DegreeRandomTest) {
//     std::mt19937 gen(0);
//     std::uniform_real_distribution<double> dis(1e-8, 1.0);
//     //std::uniform_real_distribution<double> dis2(0.0, 10.0);
//     int counter = 0;
    
//     for (int i = 0; i < 999999; i++) {
//         double val1 = dis(gen);
//         double val2 = dis(gen);
//         double val3 = dis(gen);
//         double val4 = dis(gen);
//         double val5 = dis(gen);

//         std::vector<double> out = calculate_real_roots(val1,val2,val3,val4,val5);
//         for (int i = 0 ; i < out.size(); i++) {
//             double x = out[i];
//             double x2 = x * x;
            
//             double func = x2 * x2 * val1 + x2 * x * val2 + x2 * val3 + x * val4 + val5;
//             if (fabs(func) > 1e-5) {
//                 counter++;
//             }
//         }
//     }
//     EXPECT_LE(counter,250);
// }