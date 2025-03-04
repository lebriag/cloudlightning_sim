#ifndef POWER_H
#define POWER_H

class powinputs;

class power
{
 private:
  int alloc;
  int typeCpu;
  int typeAcc;
  int accelerator;
  double cpuPmin, cpuPmax;
  double cpuC;
  int numOfPoints;
  double* cpubins;
  double* cpuP;
  double accPmin, accPmax;
  double accC;
  double *a, *b, *c, *d;

 public:
  power();

  ~power();

  power(const powinputs& t);

  power(const power& t);

  power& operator=(const power& t);

  /// Defines the CPU power models of the simulator. There are 7 different models included, 5 global and 2 piecewise
  double modelCPU(double& u);

  /// Defines the accelerator power models of the simulator. A single composite accelerator power model is included
  double modelACC(double& rho, int& numAcc);

  /// Returns the total consumption of the system
  double consumption(double& u, double& rho, int& active, int& numAcc);

  int galloc() const;
  int getTypeCpu() const;
  int getTypeAccelerator() const;
  int getAccelerator() const;
  int gnumOfPoints() const;

  double gcpuPmin() const;
  double gcpuPmax() const;
  double gcpuC() const;

  double* gcpubins() const;
  double* gcpuP() const;
  double gaccPmin() const;
  double gaccPmax() const;
  double gaccC() const;
  double* ga() const;
  double* gb() const;
  double* gc() const;
  double* gd() const;
  void print() const;
};

#endif
