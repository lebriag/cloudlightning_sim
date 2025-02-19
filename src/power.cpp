#include <inputs.h>
#include <power.h>

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

power::power()
  : alloc(0),
    typeCpu(0),
    typeAcc(0),
    accelerator(0),
    cpuPmin(0.0),
    cpuPmax(0.0),
    cpuC(0.0),
    numOfPoints(0),
    cpubins(nullptr),
    cpuP(nullptr),
    accPmin(0.0),
    accPmax(0.0),
    accC(0.0),
    a(nullptr),
    b(nullptr),
    c(nullptr),
    d(nullptr)
{
}

power::power(const powinputs& t)
{
  int i;
  double *ta, *tb, *tc, *h, *s;
  alloc = 1;
  typeCpu = t.typeCpu;
  cpuPmin = t.cpuPmin;
  cpuPmax = t.cpuPmax;
  cpuC = t.cpuC;

  if (typeCpu < 0) {
    cpubins = nullptr;
    cpuP = nullptr;
    a = nullptr;
    b = nullptr;
    c = nullptr;
    d = nullptr;
    numOfPoints = 0;
  } else if (typeCpu > 0) {
    numOfPoints = t.numOfPoints;
    cpubins = new double[numOfPoints];
    cpuP = new double[numOfPoints];

    for (i = 0; i < numOfPoints; i++) {
      cpubins[i] = t.cpubins[i];
      cpuP[i] = t.cpuP[i];
    }

    if (typeCpu == 2) {
      a = new double[numOfPoints - 1];
      b = new double[numOfPoints - 1];
      c = new double[numOfPoints];
      d = new double[numOfPoints - 1];
      ta = new double[numOfPoints];
      tb = new double[numOfPoints];
      tc = new double[numOfPoints];
      s = new double[numOfPoints - 1];
      h = new double[numOfPoints - 1];
      for (i = 0; i < numOfPoints - 1; i++) {
        a[i] = cpuP[i];
      }
      for (i = 0; i < numOfPoints - 1; i++) {
        h[i] = cpubins[i + 1] - cpubins[i];
      }

      // Eliminate exceeding terms to form tridiagonal
      s[0] = 0;
      s[numOfPoints - 1] = 0;
      for (i = 1; i < numOfPoints - 1; i++)
        s[i] = 3.0 * (cpuP[i + 1] - cpuP[i]) / h[i] - 3.0 * (cpuP[i] - cpuP[i - 1]) / h[i - 1];
      s[0] = -h[0] * s[1];
      s[numOfPoints - 1] = -h[numOfPoints - 2] * s[numOfPoints - 2];

      tb[0] = h[1] * h[1] - h[0] * h[0];
      for (i = 1; i < numOfPoints - 1; i++)
        tb[i] = 2.0 * (h[i - 1] + h[i]);
      tb[numOfPoints - 1] = h[numOfPoints - 3] * h[numOfPoints - 3] - h[numOfPoints - 2] * h[numOfPoints - 2];

      tc[0] = -2.0 * h[0] * h[0] - 3.0 * h[0] * h[1] - h[1] * h[1];
      for (i = 1; i < numOfPoints - 1; i++)
        tc[i] = h[i];
      for (i = 1; i < numOfPoints - 1; i++)
        ta[i] = h[i - 1];
      ta[numOfPoints - 1] = -2.0 * h[numOfPoints - 2] * h[numOfPoints - 2] -
                            3.0 * h[numOfPoints - 3] * h[numOfPoints - 2] - h[numOfPoints - 3] * h[numOfPoints - 3];

      // Thomas algorithm
      tc[0] = tc[0] / tb[0];
      for (i = 1; i < numOfPoints - 1; i++) {
        tc[i] = tc[i] / (tb[i] - ta[i] * tc[i - 1]);
      }

      s[0] = s[0] / tb[0];
      for (i = 1; i < numOfPoints; i++) {
        s[i] = (s[i] - ta[i] * s[i - 1]) / (tb[i] - ta[i] * tc[i]);
      }

      c[numOfPoints - 1] = s[numOfPoints - 1];
      for (i = numOfPoints - 2; i >= 0; i--) {
        c[i] = s[i] - tc[i] * c[i + 1];
      }
      for (i = 0; i < numOfPoints - 1; i++) {
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
      }
      for (i = 0; i < numOfPoints - 1; i++) {
        b[i] = (cpuP[i + 1] - cpuP[i]) / h[i] - c[i] * h[i] - d[i] * h[i] * h[i];
      }

      delete[] h;
      delete[] ta;
      delete[] tb;
      delete[] tc;
      delete[] s;
    }
  }
  accelerator = t.accelerator;
  typeAcc = t.typeAcc;
  accPmin = t.accPmin;
  accPmax = t.accPmax;
  accC = t.accC;
}

power::power(const power& t)
{
  if (t.galloc()) {
    alloc = 1;
    typeCpu = t.getTypeCpu();
    typeAcc = t.getTypeAccelerator();
    accelerator = t.getAccelerator();
    cpuPmin = t.gcpuPmin();
    cpuPmax = t.gcpuPmax();
    cpuC = t.gcpuC();
    numOfPoints = t.gnumOfPoints();
    a = nullptr;
    b = nullptr;
    c = nullptr;
    d = nullptr;
    cpubins = nullptr;
    cpuP = nullptr;
    if (numOfPoints > 0) {
      cpubins = new double[numOfPoints];
      cpuP = new double[numOfPoints];

      for (int i = 0; i < numOfPoints; i++) {
        cpubins[i] = *(t.gcpubins() + i);
        cpuP[i] = *(t.gcpuP() + i);
      }

      if (typeCpu == 2) {
        a = new double[numOfPoints - 1];
        b = new double[numOfPoints - 1];
        c = new double[numOfPoints];
        d = new double[numOfPoints - 1];

        for (int i = 0; i < numOfPoints - 1; i++) {
          a[i] = *(t.ga() + i);
          b[i] = *(t.gb() + i);
          c[i] = *(t.gc() + i);
          d[i] = *(t.gd() + i);
        }
        c[numOfPoints - 1] = *(t.gc() + numOfPoints - 1);
      }
    }
    accPmin = t.gaccPmin();
    accPmax = t.gaccPmax();
    accC = t.gaccC();
  }
}

power::~power()
{
  if (alloc) {
    alloc = 0;
    typeAcc = 0;
    accelerator = 0;
    cpuPmin = 0.0;
    cpuPmax = 0.0;
    cpuC = 0.0;
    if (numOfPoints > 0) {
      numOfPoints = 0;
      delete[] cpubins;
      delete[] cpuP;
      if (typeCpu == 2) {
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
      }
    }
    typeCpu = 0;
    accPmin = 0.0;
    accPmax = 0.0;
    accC = 0.0;
  }
}

power& power::operator=(const power& t)
{
  int i;
  if (this != &t) {
    if (alloc) {
      alloc = 0;
      typeAcc = 0;
      accelerator = 0;
      cpuPmin = 0.0;
      cpuPmax = 0.0;
      cpuC = 0.0;
      if (numOfPoints > 0) {
        numOfPoints = 0;
        delete[] cpubins;
        delete[] cpuP;
        if (typeCpu == 2) {
          delete[] a;
          delete[] b;
          delete[] c;
          delete[] d;
        }
      }
      typeCpu = 0;
      accPmin = 0.0;
      accPmax = 0.0;
      accC = 0.0;
    }
    alloc = t.galloc();
    if (alloc) {
      alloc = t.galloc();
      typeCpu = t.getTypeCpu();
      typeAcc = t.getTypeAccelerator();
      accelerator = t.getAccelerator();
      cpuPmin = t.gcpuPmin();
      cpuPmax = t.gcpuPmax();
      cpuC = t.gcpuC();
      numOfPoints = t.gnumOfPoints();
      a = nullptr;
      b = nullptr;
      c = nullptr;
      d = nullptr;
      cpubins = nullptr;
      cpuP = nullptr;
      if (numOfPoints > 0) {
        cpubins = new double[numOfPoints];
        cpuP = new double[numOfPoints];
        for (i = 0; i < numOfPoints; i++) {
          cpubins[i] = *(t.gcpubins() + i);
          cpuP[i] = *(t.gcpuP() + i);
        }
        if (typeCpu == 2) {
          a = new double[numOfPoints - 1];
          b = new double[numOfPoints - 1];
          c = new double[numOfPoints];
          d = new double[numOfPoints - 1];
          for (i = 0; i < numOfPoints - 1; i++) {
            a[i] = *(t.ga() + i);
            b[i] = *(t.gb() + i);
            c[i] = *(t.gc() + i);
            d[i] = *(t.gd() + i);
          }
          c[numOfPoints - 1] = *(t.gc() + numOfPoints - 1);
        }
      }
      accPmin = t.gaccPmin();
      accPmax = t.gaccPmax();
      accC = t.gaccC();
    }
  }
  return *this;
}

double power::modelCPU(double& u)
{
  double pcons = 0.0, Pmid;
  int i;
  switch (typeCpu) {
    // Global models
    case -1:
      pcons = cpuPmin + (cpuPmax - cpuPmin) * u;
      break;
    case -2:
      pcons = cpuPmin + (cpuPmax - cpuPmin) * u * u;
      break;
    case -3:
      pcons = cpuPmin + (cpuPmax - cpuPmin) * u * u * u;
      break;
    case -4:
      Pmid = cpuPmin + (cpuPmax - cpuPmin) / 2;
      pcons =
        (4.0 / 3.0 * Pmid - cpuPmin / 6.0 - cpuPmax / 3.0) +
        (4.0 / 3.0 * Pmid - 2.0 * cpuPmin / 3.0 - cpuPmax / 3.0) * u +
        (2.0 * cpuPmax + 2.0 * cpuPmin - 4.0 * Pmid) * u * u +
        (4.0 / 3.0 * Pmid - 7.0 / 6.0 * cpuPmin - cpuPmax / 3.0) * (2.0 * u - 1.0) * (2.0 * u - 1.0) * (2.0 * u - 1.0);
      break;
    case -5:
      Pmid = 5.0 * cpuPmax / 9.0;
      pcons =
        (4.0 / 3.0 * Pmid - cpuPmin / 6.0 - cpuPmax / 3.0) +
        (4.0 / 3.0 * Pmid - 2.0 * cpuPmin / 3.0 - cpuPmax / 3.0) * u +
        (2.0 * cpuPmax + 2.0 * cpuPmin - 4.0 * Pmid) * u * u +
        (4.0 / 3.0 * Pmid - 7.0 / 6.0 * cpuPmin - cpuPmax / 3.0) * (2.0 * u - 1.0) * (2.0 * u - 1.0) * (2.0 * u - 1.0);
      break;
    // Piecewise models
    case 1:
      if (u < cpubins[0])
        pcons = cpuP[0] + (cpuP[1] - cpuP[0]) * (u - cpubins[0]) / (cpubins[1] - cpubins[0]);
      else if (u > cpubins[numOfPoints - 1])
        pcons = cpuP[numOfPoints - 2] +
                (cpuP[numOfPoints - 1] - cpuP[numOfPoints - 2]) * (u - cpubins[numOfPoints - 2]) /
                  (cpubins[numOfPoints - 1] - cpubins[numOfPoints - 2]);
      else {
        for (i = 1; i < numOfPoints; i++) {
          if (u <= cpubins[i]) {
            break;
          }
        }
        pcons = cpuP[i - 1] + (cpuP[i] - cpuP[i - 1]) * (u - cpubins[i - 1]) / (cpubins[i] - cpubins[i - 1]);
      }
      break;
    case 2:
      if (u < cpubins[0]) {
        pcons = a[0] + b[0] * (u - cpubins[0]) + c[0] * (u - cpubins[0]) * (u - cpubins[0]) +
                d[0] * (u - cpubins[0]) * (u - cpubins[0]) * (u - cpubins[0]);
      } else if (u > cpubins[numOfPoints - 1]) {
        pcons = a[numOfPoints - 2] + b[numOfPoints - 2] * (u - cpubins[numOfPoints - 2]) +
                c[numOfPoints - 2] * (u - cpubins[numOfPoints - 2]) * (u - cpubins[numOfPoints - 2]) +
                d[numOfPoints - 2] * (u - cpubins[numOfPoints - 2]) * (u - cpubins[numOfPoints - 2]) *
                  (u - cpubins[numOfPoints - 2]);
      } else {
        for (i = 1; i < numOfPoints; i++) {
          if (u <= cpubins[i]) {
            break;
          }
        }
        i--;
        pcons = a[i] + b[i] * (u - cpubins[i]) + c[i] * (u - cpubins[i]) * (u - cpubins[i]) +
                d[i] * (u - cpubins[i]) * (u - cpubins[i]) * (u - cpubins[i]);
      }
      break;
    case 3:
      int ii = floor(u * 10);
      pcons = cpuP[ii] + (cpuP[ii + 1] - cpuP[ii]) * (u - 0.1 * ii) / (0.1 * (ii + 1) - 0.1 * ii);
      break;
  }
  return pcons;
}

double power::modelACC(double& rho, int& numAcc)
{
  double pcons = 0.0;
  switch (typeAcc) {
    case -1:
      pcons = accPmin * numAcc + (rho) * (accPmax - accPmin) * numAcc;
      break;
  }
  return pcons;
}

double power::consumption(double& u, double& rho, int& active, int& numAcc)
{
  if (active) {
    return (modelCPU(u) + modelACC(rho, numAcc)) * (1.0e-9) / 3600;
  } else {
    return (cpuC + numAcc * accC) * (1.0e-9) / 3600;
  }
}

int power::galloc() const { return alloc; }
int power::getTypeCpu() const { return typeCpu; }
int power::getTypeAccelerator() const { return typeAcc; }
int power::getAccelerator() const { return accelerator; }
double power::gcpuPmin() const { return cpuPmin; }
double power::gcpuPmax() const { return cpuPmax; }
double power::gcpuC() const { return cpuC; }
int power::gnumOfPoints() const { return numOfPoints; }
double* power::gcpubins() const { return cpubins; }
double* power::gcpuP() const { return cpuP; }
double power::gaccPmin() const { return accPmin; }
double power::gaccPmax() const { return accPmax; }
double power::gaccC() const { return accC; }
double* power::ga() const { return a; }
double* power::gb() const { return b; }
double* power::gc() const { return c; }
double* power::gd() const { return d; }
void power::print() const
{
  if (alloc) {
    cout << "         CPU Power Consumption model: " << typeCpu << endl;
    if (typeCpu < 0) {
      cout << "            CPU Idle Power Consumption: " << cpuPmin << " Watts" << endl;
      cout << "            CPU Max Power Consumption: " << cpuPmax << " Watts" << endl;
    } else if (typeCpu > 0) {
      cout << "            CPU Number of Points for Interpolation:  " << numOfPoints << endl;
      cout << "            CPU Utilization Bins: ";
      for (int j = 0; j < numOfPoints; j++) {
        cout << cpubins[j] << " ";
      }
      cout << endl;
      cout << "            CPU Power Consumption: ";
      for (int j = 0; j < numOfPoints; j++) {
        cout << cpuP[j] << " ";
      }
      cout << endl;
    }
    cout << "            CPU Sleep Power Consumption: " << cpuC << " Watts" << endl;
    cout << "         ACC Power Consumption model: " << typeAcc << endl;
    cout << "            ACC Idle Power Consumption: " << accPmin << " Watts" << endl;
    cout << "            ACC Max Power Consumption: " << accPmax << " Watts" << endl;
    cout << "            ACC Sleep Power Consumption: " << accC << " Watts" << endl;
  }
}
