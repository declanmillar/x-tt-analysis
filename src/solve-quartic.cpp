#include <math.h>
#include <complex>
#include <vector>
#include <iostream>
#include "solve-quartic.hpp"

const double PI = 3.14159265358979323846;

double sign(double a) {
     return (a < 0) ? -1 : (a > 0) ? 1 : 0;
}

std::vector<double> solveQuartic(double h0, double h1, double h2, double h3, double h4) {
    if (false)
        std::cout << "Solve Quartic\n";


    const double nh1 = h1 / h0;
    const double nh2 = h2 / h0;
    const double nh3 = h3 / h0;
    const double nh4 = h4 / h0;


    const double k1 = nh2 - 3 * nh1 * nh1 / 8.;
    const double k2 = nh3 + nh1 * nh1 * nh1 / 8. - nh1 * nh2 / 2.;
    const double k3 = nh4 - 3 * pow(nh1, 4) / 256. + nh1 * nh1 * nh2 / 16. - nh1 * nh3 / 4.;

    if (false) {
        std::cout << "e " << k1 << "\n";
        std::cout << "f " << k2 << "\n";
        std::cout << "g " << k3 << "\n";
    }

    std::vector<double> sols = solveCubic(2. * k1, k1 * k1 -4. * k3, -k2 * k2);


    double t1 = -1.;

    for (std::vector<double>::const_iterator it = sols.begin(); it != sols.end(); ++it) {
        if (*it > 0.)
            t1 = sqrt(*it);
    }



    double t2 = (k1 + t1 * t1 - k2 / t1) / 2.;

    if (false)
        std::cout << "first quadratic\n";

    std::vector<double> s1 = solveQuadratic(1., t1, t2);

    std::vector<double> pvxs;

    for (std::vector<double>::const_iterator it = s1.begin(); it != s1.end(); ++it) {

        if (false)
            std::cout << "pvx " << *it - nh1 / 4. << "\n";

        pvxs.push_back(*it - nh1 / 4.);
    }


    if (false)
        std::cout << "second quadratic\n";

    std::vector<double> s2 = solveQuadratic(1., -t1, k3 / t2);

    for (std::vector<double>::const_iterator it = s2.begin(); it != s2.end(); ++it) {
        if (false)
            std::cout << "pvx " << *it - nh1 / 4. << "\n";

        pvxs.push_back(*it - nh1 / 4.);
    }

    return pvxs;
}

std::vector<double> solveCubic(double s1, double s2, double s3) {
    std::vector<double> sols;

    const double q = (s1 * s1 - 3 * s2) / 9.;
    const double r = (2 * s1 * s1 * s1 - 9 * s1 * s2 + 27 * s3) / 54.;

    if (false) {
        std::cout << "Solve Cubic\n";
        std::cout << "q " << q << "\n";
        std::cout << "r " << r << "\n";
    }

    if (r*r < q*q*q) {
        const double theta = acos(r / sqrt(q*q*q) );
        const double z1_1 = -2 * sqrt(q) * cos(theta / 3.) - s1 / 3.;
        const double z1_2 = -2 * sqrt(q) * cos((theta + 2 * M_PI) / 3.) - s1 / 3.;
        const double z1_3 = -2 * sqrt(q) * cos((theta - 2 * M_PI) / 3.) - s1 / 3.;

        if (false) {
            std::cout << "r^2 < q^3 -> 3 solutions\n";
            std::cout << "z1_1 " << z1_1 << "\n";
            std::cout << "z1_2 " << z1_2 << "\n";
            std::cout << "z1_3 " << z1_3 << "\n";
        }

        sols.push_back(z1_1);
        sols.push_back(z1_2);
        sols.push_back(z1_3);
    } else {
        if (false)
            std::cout << "r^2 >= q^3 -> 1 solution" << "\n";

        const double radicant = -r + sqrt( r * r - q * q * q);
        const double powthrd = pow( fabs(radicant), 1./3.);
        const double a = sign(radicant) * powthrd;

        const double b = q / a;

        const double z1_1 = a + b - s1 /3.;


        sols.push_back(z1_1);
    }

    return sols;
}

std::vector<double> solveQuadratic(double a, double b, double c) {
    std::vector<double> sols;

    if (false)
        std::cout << "Solve Quadratic\n";

    const double b2m4ac = b * b - 4. * a * c;

    if (b2m4ac > 0.) {
        const double sol1 = (-b + sqrt(b2m4ac)) / 2. * a;
        const double sol2 = (-b - sqrt(b2m4ac)) / 2. * a;

        if (false) {
            std::cout << "Quadratic sol1 : " << sol1 << "\n";
            std::cout << "Quadratic sol2 : " << sol2 << "\n";
        }

        sols.push_back(sol1);
        sols.push_back(sol2);
    } else if (false)
        std::cout << "Quadratic has no solutions\n";

    return sols;
}

//
// bool solveQuadratic(double &a, double &b, double &c, double &root)
// {
// 	if(a == 0.0 || std::abs(a/b) < 1.0e-6)
// 	{
// 		if(std::abs(b) < 1.0e-4)
// 			return false;
// 		else
// 		{
// 			root = -c/b;
// 			return true;
// 		}
// 	}
//
// 	double discriminant = b * b - 4.0 * a * c;
// 	if(discriminant >= 0.0)
// 	{
// 		discriminant = sqrt(discriminant);
// 		root = (b + discriminant) * -0.5 / a;
// 		return true;
// 	}
//
// 	return false;
// }
//
// bool solveQuadraticOther(double a, double b, double c, double &root)
// {
// 	if(a == 0.0 || std::abs(a/b) < 1.0e-4)
// 	{
// 		if(std::abs(b) < 1.0e-4)
// 			return false;
// 		else
// 		{
// 			root = -c/b;
// 			return true;
// 		}
// 	}
//
// 	double discriminant = b * b - 4.0 * a * c;
// 	if(discriminant >= 0.0)
// 	{
// 		discriminant = sqrt(discriminant);
// 		root = (b - discriminant) * -0.5 / a;
// 		return true;
// 	}
//
// 	return false;
// }
//
// bool solveCubic(double &a, double &b, double &c, double &d, double &root)
// {
//     if(a == 0.0 || std::abs(a/b) < 1.0e-6)
//         return solveQuadratic(b, c, d, root);
//
//     double B = b/a, C = c/a, D = d/a;
//
//     double Q = (B*B - C*3.0)/9.0, QQQ = Q*Q*Q;
//     double R = (2.0*B*B*B - 9.0*B*C + 27.0*D)/54.0, RR = R*R;
//
//     // 3 real roots
//     if(RR<QQQ)
//     {
//         /* This sqrt and division is safe, since RR >= 0, so QQQ > RR,    */
//         /* so QQQ > 0.  The acos is also safe, since RR/QQQ < 1, and      */
//         /* thus R/sqrt(QQQ) < 1.                                     */
//         double theta = acos(R/sqrt(QQQ));
//         /* This sqrt is safe, since QQQ >= 0, and thus Q >= 0             */
//         double r1, r2, r3;
//         r1 = r2 = r3 = -2.0*sqrt(Q);
//         r1 *= cos(theta/3.0);
//         r2 *= cos((theta+2*PI)/3.0);
//         r3 *= cos((theta-2*PI)/3.0);
//
//         r1 -= B/3.0;
//         r2 -= B/3.0;
//         r3 -= B/3.0;
//
//         root = 1000000.0;
//
//         if(r1 >= 0.0) root = r1;
//         if(r2 >= 0.0 && r2 < root) root = r2;
//         if(r3 >= 0.0 && r3 < root) root = r3;
//
//         return true;
//     }
//     // 1 real root
//     else
//     {
//         double A2 = -pow(std::abs(R)+sqrt(RR-QQQ),1.0/3.0);
//         if (A2!=0.0) {
//             if (R<0.0) A2 = -A2;
//             root = A2 + Q/A2;
//         }
//         root -= B/3.0;
//         return true;
//     }
// }
//
// bool solveQuartic(double a, double b, double c, double d, double e, double &root)
// {
//     // I switched to this method, and it seems to be more numerically stable.
//     // http://www.gamedev.n...topic_id=451048
//
//     // When a or (a and b) are magnitudes of order smaller than C,D,E
//     // just ignore them entirely. This seems to happen because of numerical
//     // inaccuracies of the line-circle algorithm. I wanted a robust solver,
//     // so I put the fix here instead of there.
//     if(a == 0.0 || std::abs(a/b) < 1.0e-5 || std::abs(a/c) < 1.0e-5 || std::abs(a/d) < 1.0e-5)
//         return solveCubic(b, c, d, e, root);
//
//     double B = b/a, C = c/a, D = d/a, E = e/a;
//     double BB = B*B;
//     double I = -3.0*BB*0.125 + C;
//     double J = BB*B*0.125 - B*C*0.5 + D;
//     double K = -3*BB*BB/256.0 + C*BB/16.0 - B*D*0.25 + E;
//
//     double z;
//     bool foundRoot2 = false, foundRoot3 = false, foundRoot4 = false, foundRoot5 = false;
//     double a0 = 1.0;
//     double a1 = I + I;
//     double a2 = I*I - 4*K;
//     double a3 =  -(J*J);
//     if(solveCubic(a0, a1, a2, a3, z))
//     {
//         double value = z*z*z + z*z*(I+I) + z*(I*I - 4*K) - J*J;
//
//         double p = sqrt(z);
//         double r = -p;
//         double q = (I + z - J/p)*0.5;
//         double s = (I + z + J/p)*0.5;
//
//         bool foundRoot = false, foundARoot;
//         double aRoot;
//         double b0 = 1.0;
//         double b1 = p;
//         double b2 = q;
//         foundRoot = solveQuadratic(b0, b1, b2, root);
//         root -= B/4.0;
//
//         double c0 = 1.0;
//         double c1 = r;
//         double c2 = s;
//         foundARoot = solveQuadratic(c0, c1, c2, aRoot);
//         aRoot -= B/4.0;
//         if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0)
//             || root < 0.0)) || (!foundRoot && foundARoot))
//         {
//             root = aRoot;
//             foundRoot = true;
//         }
//
//         foundARoot = solveQuadraticOther(1.0, p, q, aRoot);
//         aRoot -= B/4.0;
//         if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0)
//             || root < 0.0)) || (!foundRoot && foundARoot))
//         {
//             root = aRoot;
//             foundRoot = true;
//         }
//
//         foundARoot = solveQuadraticOther(1.0, r, s, aRoot);
//         aRoot -= B/4.0;
//         if((foundRoot && foundARoot && ((aRoot < root && aRoot >= 0.0)
//             || root < 0.0)) || (!foundRoot && foundARoot))
//         {
//             root = aRoot;
//             foundRoot = true;
//         }
//         return foundRoot;
//     }
//     return false;
// }
