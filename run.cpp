#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;


double accelDrag(double v, double z) {
    return -0.5 * pow(2.71828182846, (-1 * z / 8000)) * v * fabs(v) / 100;
}
double accelGravity(double z) {
    return -9.8 * pow( (1 + z / 6370000), -2);
}

int main() {   
    double z = 30000;

    double v = 0;
    double t = 0;
    double dt = 0.1; // step size
    double a = -9.8;

    ofstream data;
    data.open("data.csv");

    double ap;
    double ac;
    double zp;
    double zc;
    double vp;
    double vc;

    while (z > 0) {
        // euler step:
        // a = ( accelGravity(z) + accelDrag(v,z) );
        // z += v * dt;
        // v += a * dt;
        // t += dt;

        ap = (accelGravity(z) + accelDrag(v,z));
        zp = z + v * dt;
        vp = v + ap * dt;
    
        // corrector:
        ac = (accelGravity(zp) + accelDrag(vp,zp));
        zc = z + vp * dt;
        vc = v + ac * dt;

        // Average solutions for final z value:
        z = 0.5 * (zp + zc);
        v = 0.5 * (vp + vc);
        a = 0.5 * (ap + ac);

        // To eliminate velocity-dependent forces, 
        // take the midpoint of the initial and final states using the predictor-corrector, 
        // and use that to give a final update to the velocity. 
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



    return 1;
}
