#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;


double accelDrag(double v, double z, double mass) {
    return -0.043 * pow(2.71828182846, (-1 * z / 8000)) * fabs(v) / mass;
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
    // TAN 54 DEGREES IS 1.37, MY RATIO IS DEFINITLY WRONG

    // an azimuth of 90◦ is east, 180◦ is south, and 270◦ is west, 0 is north
    double azimuth = 144;
    double altitude = 28;

    //double DragCoefficient1 = 0.043;
    double mass = 10;

    // the x direction is aligned along east,y along north, and z along height.
    double range = 0;
    double x = 0;
    double y = 0;
    double z = 4;

    double v = 100;
    double vx = getXratio(azimuth, altitude) * v;
    double vy = getYratio(azimuth, altitude) * v;
    double vz = getZratio(altitude) * v;

    double a;
    double ax;
    double ay;
    double az;

    double t = 0;
    double dt = 0.01; // step size

    ofstream data;
    data.open("data.csv");

    double xp = 0;
    double xc = 0;
    double yp = 0;
    double yc = 0;
    double zp = 4;
    double zc = 4;

    double vp = v;
    double vc = v;
    double vxp = vx;
    double vxc = vx;
    double vyp = vy;
    double vyc = vy;
    double vzp = vz;
    double vzc = vz;

    double ap;
    double ac;
    double axp;
    double axc;
    double ayp;
    double ayc;
    double azp;
    double azc;
    
    while (z > 0) {
//  while r[2]>0 or it==0:
//         # first take trial step
//         a = acceleration(r,v,mass)
//         for i in range(3):
//             rp[i]=r[i]+v[i]*dt
//             vp[i]=v[i]+a[i]*dt
    
//         # now add a corrector, which is just an estimate using the acceleration at the end of the trial step
//         ap = acceleration(rp,vp,mass)
//         for i in range(3):
//             rc[i] = r[i]+vp[i]*dt
//             vc[i] = v[i]+ap[i]*dt

//         ac = acceleration(rc,vc,mass) # this step is mainly for record keeping
//         # now average the solutions for final z estimate

//         # this is just a shortcut to ensure the new array is not just a pointer to the old array
//         r0=r*1
//         v0=v*1

//         r = 0.5*(rp+rc)
//         v = 0.5*(vp+vc)
//         a = 0.5*(ap+ac)

//         t+=dt

        double dragTerm = accelDrag(v, z, mass);
        ax = dragTerm*vx;
        ay = dragTerm*vy;
        az = dragTerm*vz + accelGravity(z);
        a = getTotalFromComponents(ax, ay, az);

        xp = x + vx * dt;
        yp = y + vy * dt;
        zp = z + vz * dt;

        vxp = vx + ax * dt;
        vyp = vy + ay * dt;
        vzp = vz + az * dt;
        vc = getTotalFromComponents(vxp, vyp, vzp);

        // corrector:

        dragTerm = accelDrag(vp, zp, mass);
        axp = dragTerm*vxp;
        ayp = dragTerm*vyp;
        azp = dragTerm*vzp + accelGravity(zp);
        ap = getTotalFromComponents(axp, ayp, azp);
        
        xc = x + vxp * dt;
        yc = y + vyp * dt;
        zc = z + vzp * dt;

        vxc = vx + axp * dt;
        vyc = vy + ayp * dt;
        vzc = vz + azp * dt;
        vc = getTotalFromComponents(vxc, vyc, vzc);
//
        dragTerm = accelDrag(vc, zc, mass);
        axc = dragTerm*vxc;
        ayc = dragTerm*vyc;
        azc = dragTerm*vzc + accelGravity(zc);
        ac = getTotalFromComponents(axc, ayc, azc);

        // Average solutions for final values:
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
            cout<<"  vx: "<<vx;
            cout<<"  vy: "<<vy;
            cout<<"  vz: "<<vz;
            cout<<"  height: "<<z;
            cout<<"  speed: "<<v;
            cout<<"  accel: "<<a;
            cout<<"  x: "<<x;
            cout<<"  y: "<<y;
        }

        data << range <<","<< z <<","<< t <<","<< az <<","<< v << endl;
    }
    return 1;
}
