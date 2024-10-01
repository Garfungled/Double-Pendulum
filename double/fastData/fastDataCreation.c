// Stricly for creating lyapunov data quickly for a double pendulum, adjust if needed

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define g 9.80665L
#define pi 3.141592653589793238462643383279502884L

typedef struct {
    long double properties[4];
} Solution;

typedef struct {
    long double m1; long double m2;
    long double l1; long double l2;
    Solution sol;
} Pendulum;

Pendulum pendulumInit(long double theta1, long double theta2, long double m1, long double m2, long double l1, long double l2) {
    Pendulum retP;
    retP.m1 = m1; retP.m2 = m2;
    retP.l1 = l1; retP.l2 = l2;
    retP.sol.properties[0] = theta1; retP.sol.properties[1] = theta2;
    retP.sol.properties[2] = 0.0L; retP.sol.properties[3] = 0.0L;

    return retP;
}

Solution addSol(Solution sol1, Solution sol2) {
    Solution addedSol = sol1;
    addedSol.properties[0] += sol2.properties[0]; addedSol.properties[1] += sol2.properties[1];
    addedSol.properties[2] += sol2.properties[2]; addedSol.properties[3] += sol2.properties[3];

    return addedSol;
}

Solution multSol(Solution sol, long double k) {
    Solution newsol = sol;
    newsol.properties[0] *= k; newsol.properties[1] *= k;
    newsol.properties[2] *= k; newsol.properties[3] *= k;
    return newsol;
}

Solution subSol(Solution sol1, Solution sol2){
    return addSol(sol1, multSol(sol2, -1));
}

long double normSol(Solution sol) {
    return sqrt(pow(sol.properties[0], 2) + pow(sol.properties[1], 2) + pow(sol.properties[2], 2) + pow(sol.properties[3], 2));
}

Solution F(Pendulum self, long double t, Solution y) {
    long double theta1 = y.properties[0]; long double theta2 = y.properties[1];
    long double omega1 = y.properties[2]; long double omega2 = y.properties[3];
    long double m1 = self.m1; long double m2 = self.m2;
    long double l1 = self.l1; long double l2 = self.l2;

    long double omega1dot = (-g * (2.0L * m1 + m2) * sin(theta1) - m2 * g * sin(theta1 - 2.0L * theta2) - 2.0L * sin(theta1 - theta2) * m2 * (pow(omega2, 2.0L) * l2 + pow(omega1, 2.0L) * l1 * cos(theta1 - theta2))) / (l1 * (2.0L * m1 + m2 - m2 * cos(2.0L * theta1 - 2.0L * theta2)));
    long double omega2dot = (2.0L * sin(theta1 - theta2) * (pow(omega1, 2.0L) * l1 * (m1 + m2) + g * (m1 + m2) * cos(theta1) + pow(omega2, 2.0L) * l2 * m2 * cos(theta1 - theta2))) / (l2 * (2.0L * m1 + m2 - m2 * cos(2.0L * theta1 - 2.0L * theta2)));

    Solution newsol;
    newsol.properties[0] = omega1; newsol.properties[1] = omega2;
    newsol.properties[2] = omega1dot; newsol.properties[3] = omega2dot;

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

long double lyapunov(long double theta1, long double theta2) {
    long double t = 0.0L;
    long double dt = 1.0L/100.0L;
    long double totalTime = 25.0L;
    int averages = 20.0L;

    long double lyapunovSum = 0.0L;
    Pendulum x = pendulumInit(theta1, theta2, 1.0L, 1.0L, 1.0L, 1.0L);
    Pendulum delta = pendulumInit(theta1 + 0.001L, theta2, 1.0L, 1.0L, 1.0L, 1.0L);

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
    int dimension = 100;

    FILE *fpt;
    char buffer[(int)floor(log10(dimension)) + 40];
    snprintf(buffer, 2*(int)floor(log10(dimension)) + 40, "data/data_%dx%d.csv", dimension, dimension);

    fpt = fopen(buffer, "w+");

    int row;
    int col;
    for(row = 0; row < dimension; row++) {
        for (col = 0; col < dimension; col++) {
            fprintf(fpt, "%Lf", lyapunov((2.0 * pi)/dimension * row, (2.0 * pi)/dimension * col));
            fprintf(fpt, (col == dimension - 1)? "" : ",");
        }
        fprintf(fpt, "\n");
    }
    fclose(fpt);

    return 0;
}