#include "util.hpp"
#include <cmath>
#include <cstdlib>

namespace anyprog {

// https://github.com/janantala/GPS-distance
// http://www.movable-type.co.uk/scripts/latlong-vincenty.html
double to_radians(double degree)
{
    double r = degree * M_PI / 180;
    return r;
}

double gps_distance(double lat1, double lon1, double lat2, double lon2, double u)
{
    double a = 6378137, b = 6356752.314245, f = 1 / 298.257223563;
    double L = to_radians(lon2 - lon1);

    double U1 = atan((1 - f) * tan(to_radians(lat1)));
    double U2 = atan((1 - f) * tan(to_radians(lat2)));
    double sinU1 = sin(U1), cosU1 = cos(U1);
    double sinU2 = sin(U2), cosU2 = cos(U2);
    double cosSqAlpha;
    double sinSigma;
    double cos2SigmaM;
    double cosSigma;
    double sigma;

    double lambda = L, lambdaP, iterLimit = 100;
    do {
        double sinLambda = sin(lambda), cosLambda = cos(lambda);
        sinSigma = sqrt((cosU2 * sinLambda)
                * (cosU2 * sinLambda)
            + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda)
                * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
        if (sinSigma == 0) {
            return 0;
        }

        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
        sigma = atan2(sinSigma, cosSigma);
        double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1 - sinAlpha * sinAlpha;
        cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

        double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));

    } while (abs(lambda - lambdaP) > 1e-12 && --iterLimit > 0);

    if (iterLimit == 0) {
        return 0;
    }

    double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
    double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
    double deltaSigma = B * sinSigma
        * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));

    double s = b * A * (sigma - deltaSigma);

    return s / u;
}
}