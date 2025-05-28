#include "roots.h"

// ignore complex and imaganery roots

namespace {

constexpr double sqrtConstExpr(double x, double curr) {
    double prev = -1;
    while (curr != prev) {
        prev = curr;
        curr = 0.5 * (curr + x / curr);
    }
    return curr;
}


double constexpr cosConstExpr(double x, double curr) {
    double prev = -1;
    while (curr != prev) {
        prev = curr;
        curr = 0.5 * (curr + x / curr);
    }
    return curr;
}

double approximate_2_cos_arccos_over_3_plus_4pi_over_3_helper(double x) {
    // https://www.wolframalpha.com/input?i=taylor+approximation+of+2+*+cos%28arccos%28x%29+%2F+3%2B4pi%2F3%29+at+x+%3D+-1+of+order+3    double x_diff = x + 1;
    double x_diff = x + 1;
    double x_diff_Sq = x_diff * x_diff;
    double x_diff_SqRoot = sqrt(x_diff);

    constexpr double c1 = -sqrtConstExpr(2.0/3.0,0.816);
    constexpr double c2 = -1/9;
    constexpr double c3 = -5.0/(54.0 * sqrtConstExpr(6,2.449));
    constexpr double c4 = -4.0/243.0;
    constexpr double c5 = -(77.0)/(3888.0 * sqrtConstExpr(6,2.449));
    constexpr double c6 = -(28.0)/6561.0;
    return 1 + c1 * x_diff_SqRoot + c2 * x_diff + c3 * x_diff_SqRoot * x_diff 
                 + c4 *  x_diff_Sq + c5 * x_diff_Sq * x_diff_SqRoot + c6 * x_diff_Sq * x_diff;
}

double approximate_2_cos_arccos_over_3_plus_4pi_over_3(double x) {
    if (x < -0.818) {
        return approximate_2_cos_arccos_over_3_plus_4pi_over_3_helper(x);
    }
    if (x > 0.818) {
        return -approximate_2_cos_arccos_over_3_plus_4pi_over_3_helper(-x);
    }

    // https://www.wolframalpha.com/input?i=pade+approximation+of+2+*+cos%28arccos%28x%29+%2F+3%2B4pi%2F3%29+at+x+%3D+0+of+order+%5B8%2F8%5D

    constexpr double t1 = 4544.0 / 98415.0;
    constexpr double t2 = -4768.0 / 10935.0;
    constexpr double t3 = 412.0 / 405.0;
    constexpr double t4 = -2.0 / 3.0;

    constexpr double b1 = 4864.0 / 2657205.0;
    constexpr double b2 = -800.0 / 6561.0;
    constexpr double b3 = 1016.0 / 1215.0;
    constexpr double b4 = -226.0 / 135.0;

    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x2*x2;
    double x5 = x2*x3;
    double x6 = x4*x2;
    double x7 = x4*x3;
    double x8 = x4*x4;

    return (t1 * x7 + t2 * x5 + t3 * x3 + t4 * x) / (b1 * x8 + b2 * x6 + b3 * x4 + b4 * x2 + 1);

}


double approximate_2_cos_arccos_over_3(double x) {
    if (x < -0.7681) {
        // Because pade fails for this case
        // Don't use more terms because of expensive sqrt and they don't give much more accercy.
        // https://www.wolframalpha.com/input?i=taylor+approximation+of+2+*+cos%28arccos%28x%29+%2F+3%29+at+x+%3D+-1+of+order+4        double x_diff = x + 1;
        double x_diff = x + 1; 
        double x_diff_Sq = x_diff * x_diff;
        double x_diff_SqRoot = sqrt(x_diff);

        constexpr double c1 = sqrtConstExpr(2.0/3.0,0.816);
        constexpr double c2 = -1.0/9.0;
        constexpr double c3 = 5.0/(54.0*sqrtConstExpr(6,2.449));
        constexpr double c4 = -4.0/243.0;
        constexpr double c5 = (77.0)/(3888.0 * sqrtConstExpr(6,2.449));
        constexpr double c6 = -(28.0)/6561.0;
        constexpr double c7 = (2431.0)/419904.0;
        constexpr double c8 = -80/59049;

        return 1 + c1 * x_diff_SqRoot + c2 * x_diff + c3 * x_diff * x_diff_SqRoot
                 + c4 * x_diff_Sq + c5 * x_diff_Sq * x_diff_SqRoot + c6 * x_diff_Sq * x_diff 
                 + c7 * x_diff_Sq * x_diff * x_diff_Sq + c8 * x_diff_Sq * x_diff_Sq;
    }

    // Approximates 2*cos(arccos(x)/3)) using a [6/6] Padé approximation.
    // https://www.wolframalpha.com/input?i=pade+approximation+of+2+*+cos%28arccos%28x%29+%2F+3%29+at+x+%3D+0+of+order+%5B6%2F6%5D    double x = imaginary / real;
    double x2 = x*x;
    double x3 = x2 * x;
    double x4 = x2*x2;
    double x5 = x4*x;
    double x6 = x3 * x3;

    constexpr double sqrt3 = sqrtConstExpr(3.0,1.732);
    constexpr double t1 = 6367150827790091.0 / (1500694954217744832.0*sqrt3);
    constexpr double t2 = (21315389368883117.0/(250115825702957472.0));
    constexpr double t3 = (9617895791423501.0/(6947661825082152.0*sqrt3));
    constexpr double t4 = (1807789764256883.0/(578971818756846.0));
    constexpr double t5 = (432592647843845.0/(42886801389396.0 * sqrt3));
    constexpr double t6 = (110360394453383.0/(21443400694698.0));

    constexpr double b1 = 1599678636998003.0/4502084862653234496.0;
    constexpr double b2 = 3425084203314289.0/(83371941900985824.0 * sqrt3);
    constexpr double b3 = 6169664756291261.0/20842985475246456.0;
    constexpr double b4 = 459206458924015.0/(192990606252282.0 * sqrt3);
    constexpr double b5 = 370932051927533.0/128660404168188.0;
    constexpr double b6 = 34404198073939.0/(7147800231566.0 * sqrt3);

    return (t1*x6 + t2 * x5 + t3*x4 + t4 * x3 + t5 * x2 + t6 * x + sqrt3) / 
            (b1*x6 + b2 * x5 + b3*x4 + b4 * x3 + b5 * x2 + b6 * x + 1.0);
}

double calculate_real_root_helper(double b,double c,double d) {
    std::vector<double> v = calculate_real_roots(1.0,b,c,d);
    return *max_element(v.begin(), v.end());
}

inline double safe_sqrt(double value) {
    return sqrt((value <= 0) ? 0 : value);
}

#define MIN_ZERO -1e-11

unsigned int factorial(unsigned int n)
{
    double out = 1.0;
    while (n > 1) {
        out *= n;
        n -= 1;
    }
    return out;
}

inline double newton_step(double x, double a, double b, double c, double d, double e) {
    // https://arxiv.org/pdf/2003.00372v1
    double xSq = x * x;

    double function = xSq * xSq * a + xSq * x * b + xSq * c + x * d + e;
    double dir1 = xSq * x * a * 4.0 + xSq * b * 3.0 + x * c * 2.0 + d;
    double dir2 = xSq * a * 12.0 + x * b * 6.0 + c * 2.0;
    double dir3 = x * a * 24.0 + b * 6.0;
    double dir4 = a * 24.0;

    // max of p(j)(x) / j!
    double A = std::max(
        std::max(fabs(function),fabs(dir1)),
        std::max(std::max(fabs(dir2)/2.0, fabs(dir3)/6.0), fabs(dir4)/24.0)
    );

    // as all are real there complex conjagate are equal
    if (fabs(dir1) > 1e-8) {
        return x - function * dir1 / (9.0 * A * A);
    } else {
        unsigned int k = 1;
        double dirk = function;
        if (fabs(dir2) > 1e-8) {
            k = 2;
            dirk = dir2;
        } else if (fabs(dir3) > 1e-8) {
            k = 3;
            dirk = dir3;
        } else if (fabs(dir4) > 1e-8) {
            k = 4;
            dirk = dir4;
        }

        double u = function * dirk / factorial(k);
        //double c = fabs(2 * pow(u,k-1));
        double C = 2.0*u/(6*A*A); // (did some calulation and took inot account the sign of UK)
        if (u < 0 and k % 2 == 0) {// pow(u,k-1) < 0 
            // anlge = 0
            return x + C / 3.0;
        } else {
            if (k == 1) {
                return x - C / 3.0;
            } if (k == 2) {
                return x; // ignore imagnery movement
            } else if (k == 3) {
                return x + 0.5 * C / 3.0;
            } else {
                return x + sqrtConstExpr(0.5,0.7071) * C / 3.0;
            } 
        }
    }
    
}

inline double newton_step(double x, double a, double b, double c, double d) {
    // https://arxiv.org/pdf/2003.00372v1
    double xSq = x * x;

    double function = xSq * x * a + xSq * b + x * c + d;
    double dir1 = 3.0 * xSq * x * a + 2.0 * x * b + c;
    double dir2 = 6.0 * x * a + 2.0 * b;
    double dir3 = 6.0 * a;

    // max of p(j)(x) / j!
    double A = std::max(
        std::max(fabs(function),fabs(dir1)),
        std::max(fabs(dir2)/2.0, fabs(dir3)/6.0)
    );

    // as all are real there complex conjagate are equal
    return x - function * dir1 / (9.0 * A * A);
}


}



std::vector<double> calculate_real_roots(double a,double b,double c,double d,double e) {
    if (a == 0) {
        return calculate_real_roots(b,c,d,e);
    }

    // https://quarticequations.com/Quartic2.pdf use modifyed NBS method
    double A3 = b/a;
    double A2 = c/a;
    double A1 = d/a;
    double A0 = e/a;
    
    
    double u = calculate_real_root_helper(
        -A2,
        A1*A3-4.0*A0,
        4.0*A0*A2 - A1*A1 - A0*A3*A3
    );
    // the following sqrts should all be defined but because of float prescion error it can fail, so use safe_sqrt
    double psub = safe_sqrt(A3*A3/4.0 + u - A2);
    double p1 = A3/2.0 - psub;
    double p2 = A3/2.0 + psub;


    double qsign = (A1 - A3*u/2.0) > 0 ? 1 : -1;
    double qsub = safe_sqrt(u*u/4.0 - A0);
    double q1 = u/2.0 + qsign * qsub;
    double q2 = u/2.0 - qsign * qsub;

    double inner1 = p1*p1/4.0 - q1;
    double inner2 = p2*p2/4.0 - q2;

    std::vector<double> roots;

    if (inner1 >= MIN_ZERO ) {
        double root = safe_sqrt(inner1);
        roots.push_back(-p1/2.0 + root);
        roots.push_back(-p1/2.0 - root);
    }

    if (inner2 >= MIN_ZERO ) {
        double root = safe_sqrt(inner2);
        roots.push_back(-p2/2.0 + root);
        roots.push_back(-p2/2.0 - root);
    }

    for (std::vector<double>::iterator it = roots.begin(); it != roots.end(); ++it) {
        for (int i = 0; i < 10; i++) {
        
            double x = *it;
            double x2 = x*x;
            double x3 = x2*x;
            double x4 = x2*x2;

            double func = a * x4 + x3 * b + x2 * c + x * d + e;
            if (fabs(func) < 1e-8) {
                break;
            }

            
            *it = newton_step(x,a,b,c,d,e);
        }
    }
    
    return roots;
}


std::vector<double> calculate_real_roots(double a,double b,double c,double d) {
    if (abs(a) <= 1e-8) {
        return calculate_real_roots(b,c,d);
    }

    // https://proofwiki.org/wiki/Cardano%27s_Formula
    double Q = (3.0*a*c - b*b) / (9.0*a*a);
    double R = (9.0*a*b*c - 27.0*a*a*d - 2.0*b*b*b) / (54.0*a*a*a);

    std::vector<double> roots;
    double D = Q*Q*Q + R*R;
    if (D > 0) {
        double inner = sqrt(D);
        double S = cbrt(R+inner);
        double T = cbrt(R-inner);

        
        roots = std::vector<double>({
            S + T - b / (3.0*a)
        });
    } else if (D >= MIN_ZERO) {
        double S = cbrt(R);
        roots = std::vector<double>({
            2.0 * S - b / (3.0*a),
            -S - b / (3.0*a)
        });
    } else {
        // https://proofwiki.org/wiki/Cardano%27s_Formula/Trigonometric_Form
        double sqQ = safe_sqrt(-Q);
        double ratio = R / safe_sqrt(-(Q*Q*Q));
        double part2 = -b / (3.0*a);
        
        roots = std::vector<double>({
            sqQ * approximate_2_cos_arccos_over_3(ratio) + part2,
            -sqQ * approximate_2_cos_arccos_over_3(-ratio) + part2, // cos(arccos(x)/3 + 2 * pi / 3) = -cos(arccos(-x)/3)
            sqQ * approximate_2_cos_arccos_over_3_plus_4pi_over_3(ratio) + part2
        });
    }

    for (std::vector<double>::iterator it = roots.begin(); it != roots.end(); ++it) {
        for (int i = 0; i < 5; i++) {
            *it = newton_step(*it,a,b,c,d);
        }
    }

    return roots;
}


std::vector<double> calculate_real_roots(double a,double b,double c) {
    if (a == 0) {
        return std::vector<double>({-c / b});
    }
    double inner = b*b-4.0*a*c;
    if (inner > 0) {
        inner = sqrt(inner);
        return std::vector<double>({
            (-b+inner)/(2.0*a),
            (-b-inner)/(2.0*a)
        });
    } else if (inner >= MIN_ZERO) {
        return std::vector<double>({
            -b/(2.0*a)
        });
    } else {
        return std::vector<double>();
    }
}