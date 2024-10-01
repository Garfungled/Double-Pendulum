// Stricly for creating lyapunov data quickly for a triple pendulum, adjust if needed
// Possibly make all l1, l2, m1, m2 = 1.0 hard coded ? improve memory ? speed ?

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define g 9.80665L
#define pi 3.141592653589793238462643383279502884L

typedef struct {
    // properties = [theta1, theta2, theta3, omega1, omega2, omega3]
    long double properties[6];
} Solution;

typedef struct {
    long double m1; long double m2; long double m3;
    long double l1; long double l2; long double l3;
    Solution sol;
} Pendulum;

Pendulum pendulumInit(long double theta1, long double theta2, long double theta3, 
                      long double m1, long double m2, long double m3, 
                      long double l1, long double l2, long double l3) {
    Pendulum retP;
    retP.m1 = m1; retP.m2 = m2; retP.m3 = m3;
    retP.l1 = l1; retP.l2 = l2; retP.l3 = l3;
    retP.sol.properties[0] = theta1; retP.sol.properties[1] = theta2; retP.sol.properties[2] = theta3;
    retP.sol.properties[3] = 0.0L;   retP.sol.properties[4] = 0.0L;   retP.sol.properties[5] = 0.0L;

    return retP;
}

Solution solutionInit(long double theta1, long double theta2, long double theta3, 
                      long double omega1, long double omega2, long double omega3) {
    Solution retS;
    retS.properties[0] = theta1; retS.properties[1] = theta2; retS.properties[2] = theta3; 
    retS.properties[3] = omega1; retS.properties[4] = omega2; retS.properties[5] = omega3;
    return retS;
}

Solution addSol(Solution sol1, Solution sol2) {
    Solution addedSol = sol1;
    addedSol.properties[0] += sol2.properties[0]; addedSol.properties[1] += sol2.properties[1];
    addedSol.properties[2] += sol2.properties[2]; addedSol.properties[3] += sol2.properties[3];
    addedSol.properties[4] += sol2.properties[4]; addedSol.properties[5] += sol2.properties[5];

    return addedSol;
}

Solution multSol(Solution sol, long double k) {
    Solution newsol = sol;
    newsol.properties[0] *= k; newsol.properties[1] *= k;
    newsol.properties[2] *= k; newsol.properties[3] *= k;
    newsol.properties[4] *= k; newsol.properties[5] *= k;
    return newsol;
}

Solution subSol(Solution sol1, Solution sol2){
    return addSol(sol1, multSol(sol2, -1));
}

long double normSol(Solution sol) {
    return sqrt(pow(sol.properties[0], 2) + pow(sol.properties[1], 2) + pow(sol.properties[2], 2) + pow(sol.properties[3], 2) + pow(sol.properties[4], 2) + pow(sol.properties[5], 2));
}


// X = B * inverse A, where X is the list of alphas(second derivative of theta)
Solution F(Pendulum self, long double t, Solution y) {
    long double theta1 = y.properties[0]; long double theta2 = y.properties[1]; long double theta3 = y.properties[2];
    long double omega1 = y.properties[3]; long double omega2 = y.properties[4]; long double omega3 = y.properties[5];
    long double m1 = self.m1; long double m2 = self.m2; long double m3 = self.m3;
    long double l1 = self.l1; long double l2 = self.l2; long double l3 = self.l3;

    // Future person, make this call a function
    long double a11 = l1*(m1 + m2 + m3); 
    long double a12 = l2*cos(theta1 - theta2)*(m2 + m3); 
    long double a13 = l3*m3*cos(theta1 - theta3);

    long double a21 = l1*cos(theta2 - theta1)*(m2 + m3); 
    long double a22 = l2*(m2 + m3); 
    long double a23 = l3*m3*cos(theta2 - theta3);

    long double a31 = l1*m3*cos(theta3 - theta1); 
    long double a32 = l2*m3*cos(theta3 - theta2); 
    long double a33 = l3*m3;

    long double b1 = -g*sin(theta1)*(m1 + m2 + m3) - (l2*pow(omega2, 2)*sin(theta1 - theta2)*(m2 + m3) + l3*pow(omega3, 2)*sin(theta1 - theta3)*m3);
    long double b2 = -g*sin(theta2)*(m2 + m3) - (l1*pow(omega1, 2)*sin(theta2 - theta1)*(m2 + m3) + l3*pow(omega3, 2)*sin(theta2 - theta3)*m3);
    long double b3 = -g*sin(theta3)*m3 - (l1*pow(omega1, 2)*sin(theta3 - theta1)*m3 + l2*pow(omega2, 2)*sin(theta3 - theta2)*m3);

    long double detA = a13*(a21*a32 - a22*a31) + a12*(a23*a31 - a21*a33) + a11*(a22*a33 - a23*a32);

    long double alpha1 = (b1*(a22*a33 - a23*a32) + b2*(a23*a31 - a21*a33) + b3*(a21*a32 - a22*a31))/detA;
    long double alpha2 = (b1*(a13*a32 - a12*a33) + b2*(a11*a33 - a13*a31) + b3*(a12*a31 - a11*a32))/detA;
    long double alpha3 = (b1*(a12*a23 - a13*a22) + b2*(a13*a21 - a11*a23) + b3*(a11*a22 - a12*a21))/detA;

    Solution newsol;
    newsol.properties[0] = omega1; newsol.properties[1] = omega2; newsol.properties[2] = omega3;
    newsol.properties[3] = alpha1; newsol.properties[4] = alpha2; newsol.properties[5] = alpha3;

    return newsol;
}

Solution RK4(Pendulum self, long double t, Solution y, long double dt) {
    Solution k1 = F(self, t, y);
    Solution k2 = F(self, t + 0.5L*dt, addSol(y, multSol(k1, 0.5L*dt)));
    Solution k3 = F(self, t + 0.5L*dt, addSol(y, multSol(k2, 0.5L*dt)));
    Solution k4 = F(self, t + dt, addSol(y, multSol(k3, dt)));

    return multSol(addSol(k1, addSol(multSol(k2, 2.0L), addSol(multSol(k3, 2.0L), k4))), dt/6.0L);
}

Pendulum step(Pendulum self, long double t, long double dt) {
    self.sol = addSol(self.sol, RK4(self, t, self.sol, dt));
    return self;
}

Solution updateDelta(long double deltaF, long double deltaI, Pendulum delta, Pendulum x) {
    return addSol(x.sol, multSol(subSol(delta.sol, x.sol), deltaI/deltaF));
}

long double partialLyapunov(long double deltaF, long double deltaI) {
    return log(fabs(deltaF/deltaI));
}

long double lyapunov(long double theta1, long double theta2, long double theta3) {
    long double t = 0.0L;
    long double dt = 1.0L/100.0L;
    long double totalTime = 25.0L;
    int averages = 20.0L;

    long double lyapunovSum = 0.0L;
    Pendulum x = pendulumInit(theta1, theta2, theta3, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L);
    Pendulum delta = pendulumInit(theta1 + 0.001L, theta2, theta3, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L, 1.0L);

    long double deltaI = normSol(subSol(delta.sol, x.sol));

    int i = 0;
    for (t = 0.0L; t < totalTime; t += dt) {
        i++;
        x = step(x, t, dt);
        delta = step(delta, t, dt);

        if(i % averages == 0) {
            long double deltaF = normSol(subSol(delta.sol, x.sol));

            lyapunovSum += partialLyapunov(deltaF, deltaI);
            delta.sol = updateDelta(deltaF, deltaI, delta, x);
            deltaI = normSol(subSol(delta.sol, x.sol));
        }
    }

    return lyapunovSum/totalTime;
}

int main() {
    // Constants
    int dimension = 1000;
    long double pixelStep = (2.0L * pi)/dimension;
    long double theta3step = 0.015625;
    long double dtheta3 = theta3step * 2.0L * pi;
    long double maxIteration = 2.0L*pi + 2.0L;
    long double maxDigitLength = (int)floor(log10(maxIteration));

    int i = 0;
    for(long double theta3 = 0; theta3 <= 2.0L*pi; theta3 += dtheta3) {
        FILE *fpt;
        char buffer[40];
        snprintf(buffer, 40, "data/1000x1000/data_%d.csv", i);

        fpt = fopen(buffer, "w+");

        int row;
        int col;
        for(row = 0; row < dimension; row++) {
            for (col = 0; col < dimension; col++) {
                fprintf(fpt, "%Lf", lyapunov(pixelStep * row, pixelStep * col, theta3));
                fprintf(fpt, (col == dimension - 1)? "" : ",");
            }
            fprintf(fpt, "\n");
        }
        fclose(fpt);
        i++;
    }

    return 0;
}