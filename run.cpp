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
double getTotalFromComponents(double x, double y, double z) {
    return sqrt (pow(x, 2) + pow(y, 2) + pow(z, 2));
}

int main() {
    double azimuth = 0.400 * 360;
    double altitude = 0.400 * 70;

    double range = 0;
    double x = 0;
    double y = 0;
    double z = 4;

    double v = 100;
    double vx = getXratio(azimuth, altitude) * v;
    double vy = getYratio(azimuth, altitude) * v;
    double vz = getZratio(altitude) * v;

    double a;
    double ax = 0;
    double ay = 0;
    double az = accelGravity(z);

    double t = 0;
    double dt = 0.01; // step size

    ofstream data;
    data.open("data.csv");

    double xp;
    double xc;
    double yp;
    double yc;
    double zp;
    double zc;

    double vp;
    double vc;
    double vxp;
    double vxc;
    double vyp;
    double vyc;
    double vzp;
    double vzc;

    double ap;
    double ac;
    double axp;
    double axc;
    double ayp;
    double ayc;
    double azp;
    double azc;
    
    while (z > 0) {
        axp = accelDrag(vx,z);
        ayp = accelDrag(vy,z);
        azp = accelGravity(z) + accelDrag(vz,z);
        ap = getTotalFromComponents(axp, ayp, azp);

        xp = x + vx * dt;
        yp = y + vy * dt;
        zp = z + vz * dt;

        vp = v + ap * dt;
        vxp = vx + axp * dt;
        vyp = vy + ayp * dt;
        vzp = vz + azp * dt;


        // corrector:
        axc = (accelDrag(vxp,zp));
        ayc = (accelDrag(vyp,zp));
        azc = (accelGravity(zp) + accelDrag(vzp,zp));
        ac = getTotalFromComponents(axc, ayc, azc);
        
        xc = x + vzp * dt;
        yc = y + vzp * dt;
        zc = z + vzp * dt;

        
        vc = v + azc * dt;
        vxc = vx + axc * dt;
        vyc = vy + ayc * dt;
        vzc = vz + azc * dt;


        // Average solutions for final z value:
        x = 0.5 * (xp + xc);
        y = 0.5 * (yp + yc);
        z = 0.5 * (zp + zc);

        v = 0.5 * (vp + vc);
        vx = 0.5 * (vxp + vxc);
        vy = 0.5 * (vyp + vyc);
        vz = 0.5 * (vzp + vzc);

        a = 0.5 * (ap + ac);
        ax = 0.5 * (axp + axc);
        ay = 0.5 * (ayp + ayc);
        az = 0.5 * (azp + azc);
        
        // outdated: 
        // To eliminate velocity-dependent forces, 
        // take the midpoint of the initial and final states using the predictor-corrector to give a final update to the velocity. 
        // z0 = z * 1;
        // v0 = v * 1;
        // zm = 0.5 * (z + z0);
        // vm = 0.5 * (v + v0);
        // a = (accelGravity(zm) + accelDrag(vm,zm));
        // v = v0 + a * dt;

        range = getTotalFromComponents(x, y, 0);
        t += dt;

        if (z < 1) {
            cout<<"  time: "<<t;
            cout<<"  range: "<<range;
            cout<<"  height: "<<z;
            cout<<"  speed: "<<v;
            cout<<"  accel: "<<a;
        }

        data << t <<","<< range <<","<< z <<","<< v <<","<< a << endl;
    }
    return 1;
}
