#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;
#define PI 3.14159265

double getTotalFromComponents(double vector[3]) {
    return sqrt (pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}
// Cross product of two vectors. 
double * crossProduct(double a[3], double b[3]) { 
    static double crossProduct[3] = {0, 0, 0};
    crossProduct[0] = a[1] * b[2] - a[2] * b[1]; 
    crossProduct[1] = a[0] * b[2] - a[2] * b[0]; 
    crossProduct[2] = a[0] * b[1] - a[1] * b[0]; 
    return crossProduct;
}
// multiplys scalar a by vector b
double * scalarVectorMultiply(double a, double b[3]) { 
    static double result[3] = {0, 0, 0};
    result[0] = a * b[0]; 
    result[1] = a * b[1]; 
    result[2] = a * b[2]; 
    return result;
}
double * addVectors(double a[3], double b[3]) { 
    static double result[3] = {0, 0, 0};
    result[0] = a[0] + b[0]; 
    result[1] = a[1] + b[1]; 
    result[2] = a[2] + b[2]; 
    return result;
}
double * subVectors(double a[3], double b[3]) { 
    static double result[3] = {0, 0, 0};
    result[0] = a[0] - b[0]; 
    result[1] = a[1] - b[1]; 
    result[2] = a[2] - b[2]; 
    return result;
}

double accelGravity(double r[3]) {
    return -9.8 * pow( (1 + r[2] / 6370000), -2);
}
void accel(double a[3], double omega[3], double v[3], double r[3], double mass, double b) {
    double Vabs = getTotalFromComponents(v);
    double dragValues = -b * pow(2.71828182846, (-1 * r[2] / 8000)) * Vabs / mass;
    a[0] = dragValues * v[0];
    a[1] = dragValues * v[1];
    a[2] = dragValues * v[2] + accelGravity(r);
    // corrections for Earth's rotation: a' = a - 2 * omega X v' - omega X (omega X r')
    a = subVectors(a, addVectors( scalarVectorMultiply(2, crossProduct(omega, v)), crossProduct(omega, crossProduct(omega, r)) ) );
    //return a;
}

// the following helpers will break down the xyz components of any vector using asimuth and altitude
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
    // an azimuth of 90◦ is east, 180◦ is south, and 270◦ is west, 0 is north
    double azimuth = 144 * PI/180;
    double altitude = 28 * PI/180;
    double mass = 10;
    double b = 0.043; // drag coefficient
    double range = 0;
    double lambda = 49 * PI/180; // lattitude
    double omegaMag = 0.0000729;
    double * omega = new double[3]; // angular acceleration, corrected for the east-north-up coordinate system below
    omega[0] = 0;
    omega[1] = cos(lambda) * omegaMag;
    omega[2] = sin(lambda) * omegaMag;

    // the x direction is aligned along east,y along north, and z along height.
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
    a[0] = 0, a[1] = 0, a[2] = 0;

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
        accel(a, omega, v, r, mass, b);
        for (int i = 0; i < 3; i++) {
            rp[i] = r[i] + v[i] * dt;
            vp[i] = v[i] + a[i] * dt;
        }

        // corrector:
        accel(ap, omega, vp, rp, mass, b);
        for (int i = 0; i < 3; i++) {
            rc[i] = r[i] + vp[i] * dt;
            vc[i] = v[i] + ap[i] * dt;
        }

        accel(ac, omega, vc, rc, mass, b);

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
    }

    delete [] r;
    r = NULL;
    delete [] rp;
    rp = NULL;
    delete [] rc;
    rc = NULL;
    delete [] v;
    v = NULL;
    delete [] vp;
    vp = NULL;
    delete [] vc;
    vc = NULL;
    delete [] a;
    a = NULL;
    delete [] ap;
    ap = NULL;
    delete [] ac;
    ac = NULL;
    
    return 1;
}
