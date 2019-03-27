#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;
#define PI 3.14159265

double getTotalFromComponents(double vector[3]) {
    return sqrt (pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}
// Cross product of two vectors. 
double * cross(double a[3], double b[3]) { 
    static double cross[3];
    cross[0] = a[1] * b[2] - a[2] * b[1]; 
    cross[1] = -1 * (a[0] * b[2] - a[2] * b[0]); 
    cross[2] = a[0] * b[1] - a[1] * b[0]; 
    return cross;
}
// multiplys scalar s by vector v
double * svMult(double s, double v[3]) { 
    static double result[3];
    result[0] = s * v[0]; 
    result[1] = s * v[1]; 
    result[2] = s * v[2]; 
    return result;
}
double * addVectors(double a[3], double b[3]) { 
    static double result[3];
    result[0] = a[0] + b[0]; 
    result[1] = a[1] + b[1]; 
    result[2] = a[2] + b[2]; 
    return result;
}
double * subVectors(double a[3], double b[3]) { 
    static double result[3];
    result[0] = a[0] - b[0]; 
    result[1] = a[1] - b[1]; 
    result[2] = a[2] - b[2]; 
    return result;
}

double accelGravity(double r[3]) {
    return -9.8 * pow( (1 + r[2] / 6370000), -2);
}
void accel(double a[3], double omega[3], double lambda, double v[3], double r[3], double mass, double massFuel, double burnRate, double b) {
    double Vabs = getTotalFromComponents(v);
    double dragValues = -b * pow(2.71828182846, (-1 * r[2] / 8000)) * Vabs / mass;
    a[0] = dragValues * v[0];
    a[1] = dragValues * v[1];
    a[2] = dragValues * v[2] + accelGravity(r);
    // rocket is propelled in direction of v, a = mDot * V / mass :
    if (massFuel > 0) {
        double thrustMag = burnRate * 2000 / mass;
        a[0] += thrustMag * v[0] / Vabs;
        a[1] += thrustMag * v[1] / Vabs;
        a[2] += thrustMag * v[2] / Vabs;
    }
    // corrections for Earth's rotation: a' = a - 2 * omega X v' - omega X (omega X r') :
    a[0] = subVectors(a, addVectors(svMult(2, cross(omega, v)), cross(omega, cross(omega, r))))[0];
    a[1] = subVectors(a, addVectors(svMult(2, cross(omega, v)), cross(omega, cross(omega, r))))[1];
    a[2] = subVectors(a, addVectors(svMult(2, cross(omega, v)), cross(omega, cross(omega, r))))[2];
    // Non-interial term: A=ω^2 R cosλsinλ en − ω^2 R (cosλ)^2 eup
    double omegaSquared = getTotalFromComponents(omega) * getTotalFromComponents(omega);
    a[1] -= omegaSquared * 6371000 * cos(lambda) * sin(lambda);
    a[2] -= omegaSquared * 6371000 * pow(cos(lambda), 2);
}

// the following helpers will break down the xyz components of initial velocity using asimuth and altitude
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
    double altitude = 80 * PI/180;
    double azimuth = 180 * PI/180;
    double mass = 10; // kg, half is chemical propellant
    double massFuel = 5;
    double burnRate = 0.25; // kg/s
    double b = 0.0025; // drag coefficient 0.043
    double range = 0;
    double lambda = 49 * PI/180; // lattitude
    double omegaMag = 0.0000729; // 0.0000729
    double * omega = new double[3]; // angular acceleration, corrected for the east-north-up coordinate system below
    omega[0] = 0;
    omega[1] = cos(lambda) * omegaMag;
    omega[2] = sin(lambda) * omegaMag;

    // the x direction is aligned along east,y along north, and z along height.
    double * r = new double[3];
    r[0] = 0;
    r[1] = 0;
    r[2] = 4; // initial launch height (meters)

    double * v = new double[3];
    double Vmag = 10; // launch speed m/s
    v[0] = getXratio(azimuth, altitude) * Vmag;
    v[1] = getYratio(azimuth, altitude) * Vmag;
    v[2] = getZratio(altitude) * Vmag;

    double * a = new double[3];
    a[0] = 0, a[1] = 0, a[2] = 0;

    ofstream data;
    data.open("data.csv");

    double * rp = new double[3];
    double * rc = new double[3];
    double * vp = new double[3];
    double * vc = new double[3];
    double * ap = new double[3];
    double * ac = new double[3];

    double t = 0; // seconds
    double dt = 0.001; // step size

    while (r[2] > 0) {
        accel(a, omega, lambda, v, r, mass, massFuel, burnRate, b);
        for (int i = 0; i < 3; i++) {
            rp[i] = r[i] + v[i] * dt;
            vp[i] = v[i] + a[i] * dt;
        }

        // corrector:
        accel(ap, omega, lambda, vp, rp, mass, massFuel, burnRate, b);
        for (int i = 0; i < 3; i++) {
            rc[i] = r[i] + vp[i] * dt;
            vc[i] = v[i] + ap[i] * dt;
        }

        accel(ac, omega, lambda, vc, rc, mass, massFuel, burnRate, b);

        // Average solutions for final values:
        for (int i = 0; i < 3; i++) {
            r[i] = 0.5 * (rp[i] + rc[i]);
            v[i] = 0.5 * (vp[i] + vc[i]);
            a[i] = 0.5 * (ap[i] + ac[i]);
        }

        range = sqrt(r[0]*r[0] + r[1]*r[1]);
        if (massFuel > 0) {
            mass -= burnRate * dt;
            massFuel -= burnRate * dt;
        }
        t += dt;

        if (r[2] < 0.1) {
            cout<<"  time: "<<t;
            cout<<"  range: "<<range;
            cout<<"  vx: "<<v[0];
            cout<<"  vy: "<<v[1];
            cout<<"  vz: "<<v[2];
            cout<<"  velocity: "<<getTotalFromComponents(v);
            cout<<"  accel: "<<getTotalFromComponents(a);
            cout<<"  x: "<<r[0];
            cout<<"  y: "<<r[1];
            cout<<"  z: "<<r[2];
            cout<<"  //// ";
        }

        data 
        << r[0] <<","<< r[1] 
        <<","<< getTotalFromComponents(v) <<","<< t 
        <<","<< getTotalFromComponents(a) <<","<< t 
        <<","<< r[2] <<","<< t 
        <<","<< r[2] <<","<< range 
        << endl;
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
