#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;
#define PI 3.14159265

double getTotalFromComponents(double * vector) {
    return sqrt (pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}
double accelGravity(double * a) {
    return -9.8 * pow( (1 + a[2] / 6370000), -2);
}
double * accel(double * v, double * r, double mass, double b) {
        static double a[3] = {0, 0, 0};
        double Vabs=getTotalFromComponents(v);
        double multiple = -b * pow(2.71828182846, (-1 * r[2] / 8000)) * Vabs / mass;
        a[0] = multiple * v[0];
        a[1] = multiple * v[1];
        a[2] = multiple * v[2] + accelGravity(r);
        return a;
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

    // double * rotateX(double * v, double angle) {
    //     static double s[3] = {0, 0, 0};
    //     s[0] = v[0];
    //     s[1] = cos(angle) * v[1] - sin(angle) * v[2];
    //     s[2] = sin(angle) * v[1] + cos(angle) * v[2];
    //     return s;
    // }
    // double * rotateZ(double * v, double angle) {
    //     static double s[3] = {0, 0, 0};
    //     s[0] = cos(angle) * v[0] - sin(angle) * v[1];
    //     s[1] = sin(angle) * v[0] + cos(angle) * v[1];
    //     s[2] = v[2];
    //     return s;
    // }

int main() {
    // an azimuth of 90◦ is east, 180◦ is south, and 270◦ is west, 0 is north
    double azimuth = 144*PI/180;
    double altitude = 28*PI/180;

    double mass = 10;
    double b = 0.043;
    // the x direction is aligned along east,y along north, and z along height.
    double range = 0;

    double * r = new double[3];
    r[0] = 0;
    r[1] = 0;
    r[2] = 4;

    double * v = new double[3];
    double Vmag = 100;
    v[0] = getXratio(azimuth, altitude) * Vmag;
    v[1] = getYratio(azimuth, altitude) * Vmag;
    v[2] = getZratio(altitude) * Vmag;

    double * a = new double[3];
    a[0] = 0;
    a[1] = 0;
    a[2] = 0;

    double t = 0;
    double dt = 0.01; // step size

    ofstream data;
    data.open("data.csv");

    double * rp = new double[3];
    double * rc = new double[3];

    double * vp = new double[3];
    double * vc = new double[3];

    double * ap = new double[3];
    double * ac = new double[3];

    while (r[2] > 0) {
        a = accel(v, r, mass, b);
        for (int i = 0; i < 3; i++) {
            rp[i] = r[i] + v[i] * dt;
            vp[i] = v[i] + a[i] * dt;
        }

        // corrector:
        ap = accel(vp, rp, mass, b);
        for (int i = 0; i < 3; i++) {
            rc[i] = r[i] + vp[i] * dt;
            vc[i] = v[i] + ap[i] * dt;
        }

        ac = accel(vc, rc, mass, b);

        // Average solutions for final values:
        for (int i = 0; i < 3; i++) {
            r[i] = 0.5 * (rp[i] + rc[i]);
            v[i] = 0.5 * (vp[i] + vc[i]);
            a[i] = 0.5 * (ap[i] + ac[i]);
        }

        range = sqrt(r[0]*r[0] + r[1]*r[1]);
        t += dt;

        double zzz[3] = {2, 3, 6};
        if (r[2] < 1) {
            cout<<"  time: "<<t;
            cout<<"  range: "<<range;
            cout<<"  vx: "<<v[0];
            cout<<"  vy: "<<v[1];
            cout<<"  vz: "<<v[2];
            cout<<"  velocity: "<<getTotalFromComponents(v);
            cout<<"  accel: "<<getTotalFromComponents(a);
            cout<<"  x: "<<r[0];
            cout<<"  y: "<<r[1];
            cout<<"  height: "<<r[2];
        }

        data << range <<","<< r[2] <<","<< t <<","<< getTotalFromComponents(a) <<","<< getTotalFromComponents(v) << endl;

        // delete [] r;
        // r = NULL;
        // delete [] v;
        // v = NULL;
        // delete [] a;
        // a = NULL;
        // delete [] rp;
        // rp = NULL;
        // delete [] rc;
        // rc = NULL;
        // delete [] vp;
        // vp = NULL;
        // delete [] vc;
        // vc = NULL;
        // delete [] ap;
        // ap = NULL;
        // delete [] ac;
        // ac = NULL;
    }
    return 1;
}
