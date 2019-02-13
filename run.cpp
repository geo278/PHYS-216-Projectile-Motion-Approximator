#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;


double accelDrag(double v, double z) {
    double mass = 10;
    return -0.043 * pow(2.71828182846, (-1 * z / 8000)) * v * fabs(v) / mass;
}
double accelGravity(double z) {
    return -9.8 * pow( (1 + z / 6370000), -2);
}
double getXratio(double azimuth, double altitude) {
    return cos(altitude) * sin(azimuth);
}
double getYratio(double azimuth, double altitude) {
    return cos(altitude) * cos(azimuth);
}
double getZratio(double altitude) {
    return sin(altitude);
}

int main() {
    double azimuth;
    double altitude;

    double x = 0;
    double y = 0;
    double z = 30000;

    double v = 100;
    double vx = getXratio(azimuth, altitude) * v;
    double vy = getYratio(azimuth, altitude) * v;
    double vz = getZratio(altitude) * v;

    double a;
    double ax = 0;
    double ay = 0;
    double az = accelGravity(z);

    double t = 0;
    double dt = 0.1; // step size

    ofstream data;
    data.open("data.csv");

    double xp;
    double xc;
    double yp;
    double yc;
    double zp;
    double zc;

    double vxp;
    double vxc;
    double vyp;
    double vyc;
    double vzp;
    double vzc;

    double axp;
    double axc;
    double ayp;
    double ayc;
    double azp;
    double azc;
    
    while (z > 0) {
        azp = (accelGravity(z) + accelDrag(vz,z));
        zp = z + vz * dt;
        vzp = v + azp * dt;
    
        // corrector:
        azc = (accelGravity(zp) + accelDrag(vzp,zp));
        zc = z + vzp * dt;
        vzc = v + azc * dt;

        // Average solutions for final z value:
        z = 0.5 * (zp + zc);
        v = 0.5 * (vzp + vzc);
        a = 0.5 * (azp + azc);

        // To eliminate velocity-dependent forces, 
        // take the midpoint of the initial and final states using the predictor-corrector to give a final update to the velocity. 
        // z0 = z * 1;
        // v0 = v * 1;
        // zm = 0.5 * (z + z0);
        // vm = 0.5 * (v + v0);
        // a = (accelGravity(zm) + accelDrag(vm,zm));
        // v = v0 + a * dt;

        t += dt;

        if (z < 1) {
            cout<<" fall time: "<<t;
            cout<<" fall position: "<<z;
            cout<<" fall speed: "<<v;
            cout<<" fall accel: "<<a;
        }

        data << t << "," << z << "," << v << "," << a << endl;
    }

    double range = sqrt (pow(x, 2) + pow(y, 2));

    return 1;
}
