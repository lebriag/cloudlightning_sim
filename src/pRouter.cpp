#include <mpi.h>
#include <omp.h>
#include <pRouter.h>
#include <pSwitch.h>
#include <stat.h>
#include <task.h>
#include <cstdlib>
#include <iostream>
#include <list>

class netw;
class resource;

using std::cout;
using std::endl;

pRouter::pRouter()
  : alloc(0),
    numberOfpSwitches(0),
    numberOfFunctions(0),
    pollIntervalpRouter(0.0),
    pSwitches(nullptr),
    Fs(nullptr),
    Ws(nullptr),
    availableProcessors(nullptr),
    totalProcessors(nullptr),
    availableMemory(nullptr),
    totalMemory(nullptr),
    availableAccelerators(nullptr),
    totalAccelerators(nullptr),
    availableStorage(nullptr),
    totalStorage(nullptr),
    sPMSA(nullptr),
    SIs(nullptr),
    SI(0.0),
    C(0.0),
    P(0.0),
    Pi(0.0)
{
}

pRouter::pRouter(const int& start, const int& end, const int& type, list<pSwitch>** LpSwitches,
                 const double& L_pollIntervalpRouter, const double& L_C, const double& L_P, const double& L_Pi,
                 const int& L_numberOfFunctions, const double* L_Ws)
{
  alloc = 1;
  numberOfpSwitches = end - start;
  numberOfFunctions = L_numberOfFunctions;
  pollIntervalpRouter = L_pollIntervalpRouter;
  list<pSwitch>::iterator it = LpSwitches[type]->begin();
  pSwitches = new list<pSwitch*>[1];
  for (int i = 0; i < start; i++) {
    it++;
  }
  for (int i = start; i < end; i++) {
    pSwitches->push_back(&(*it));
    it++;
  }
  Ws = new double[numberOfFunctions];
  Fs = new double[numberOfFunctions];
  for (int i = 0; i < numberOfFunctions; i++) {
    Ws[i] = L_Ws[i];
    Fs[i] = 0.0;
  }
  availableProcessors = new double[numberOfpSwitches];
  totalProcessors = new double[numberOfpSwitches];
  availableMemory = new double[numberOfpSwitches];
  totalMemory = new double[numberOfpSwitches];
  availableAccelerators = new double[numberOfpSwitches];
  totalAccelerators = new double[numberOfpSwitches];
  availableStorage = new double[numberOfpSwitches];
  totalStorage = new double[numberOfpSwitches];
  SIs = new double[numberOfpSwitches];
  SI = 0.0;
  sPMSA = new double[8];
  C = L_C;
  P = L_P;
  Pi = L_Pi;
  updateStateInfo(0.0);
}

pRouter::pRouter(const pRouter& t)
{
  if (t.galloc()) {
    alloc = 1;
    numberOfpSwitches = t.getNumberOfpSwitches();
    numberOfFunctions = t.getNumberOfFunctions();
    pollIntervalpRouter = t.gpollIntervalpRouter();
    pSwitches = new list<pSwitch*>[1];

    for (list<pSwitch*>::iterator it = t.gpSwitches()->begin(); it != t.gpSwitches()->end(); it++) {
      pSwitches->push_back(*it);
    }

    Ws = new double[numberOfFunctions];
    Fs = new double[numberOfFunctions];
    for (int i = 0; i < numberOfFunctions; i++) {
      Ws[i] = t.gWs()[i];
      Fs[i] = t.gFs()[i];
    }

    availableProcessors = new double[numberOfpSwitches];
    totalProcessors = new double[numberOfpSwitches];
    availableMemory = new double[numberOfpSwitches];
    totalMemory = new double[numberOfpSwitches];
    availableAccelerators = new double[numberOfpSwitches];
    totalAccelerators = new double[numberOfpSwitches];
    availableStorage = new double[numberOfpSwitches];
    totalStorage = new double[numberOfpSwitches];
    SIs = new double[numberOfpSwitches];

    for (int i = 0; i < numberOfpSwitches; i++) {
      availableProcessors[i] = t.getAvailableProcessors()[i];
      totalProcessors[i] = t.getTotalProcessors()[i];
      availableMemory[i] = t.getAvailableMemory()[i];
      totalMemory[i] = t.getTotalMemory()[i];
      availableStorage[i] = t.getAvailableStorage()[i];
      totalStorage[i] = t.getTotalStorage()[i];
      availableAccelerators[i] = t.getAvailableAccelerators()[i];
      totalAccelerators[i] = t.getTotalAccelerators()[i];
      SIs[i] = t.gSIs()[i];
    }

    SI = t.gSI();
    sPMSA = new double[8];

    for (int i = 0; i < 8; i++) {
      sPMSA[i] = t.gsPMSA()[i];
    }

    C = t.gC();
    P = t.gP();
    Pi = t.gPi();
  }
}

pRouter& pRouter::operator=(const pRouter& t)
{
  if (this != &t) {
    if (alloc) {
      alloc = 0;
      numberOfpSwitches = 0;
      numberOfFunctions = 0;
      pollIntervalpRouter = 0.0;
      pSwitches->clear();
      delete[] pSwitches;
      pSwitches = nullptr;
      delete[] Ws;
      delete[] Fs;
      Ws = nullptr;
      Fs = nullptr;
      delete[] availableProcessors;
      delete[] totalProcessors;
      delete[] availableMemory;
      delete[] totalMemory;
      delete[] availableStorage;
      delete[] totalStorage;
      delete[] availableAccelerators;
      delete[] totalAccelerators;
      delete[] sPMSA;
      delete[] SIs;
      availableProcessors = nullptr;
      totalProcessors = nullptr;
      availableMemory = nullptr;
      totalMemory = nullptr;
      availableAccelerators = nullptr;
      totalAccelerators = nullptr;
      availableStorage = nullptr;
      totalStorage = nullptr;
      sPMSA = nullptr;
      SIs = nullptr;
      SI = 0.0;
      C = 0.0;
      P = 0.0;
      Pi = 0.0;
    }
    alloc = t.galloc();
    if (alloc) {
      numberOfpSwitches = t.getNumberOfpSwitches();
      numberOfFunctions = t.getNumberOfFunctions();
      pollIntervalpRouter = t.gpollIntervalpRouter();
      pSwitches = new list<pSwitch*>[1];

      for (list<pSwitch*>::iterator it = t.gpSwitches()->begin(); it != t.gpSwitches()->end(); it++) {
        pSwitches->push_back(*it);
      }

      Ws = new double[numberOfFunctions];
      Fs = new double[numberOfFunctions];

      for (int i = 0; i < numberOfFunctions; i++) {
        Ws[i] = t.gWs()[i];
        Fs[i] = t.gFs()[i];
      }

      availableProcessors = new double[numberOfpSwitches];
      totalProcessors = new double[numberOfpSwitches];
      availableMemory = new double[numberOfpSwitches];
      totalMemory = new double[numberOfpSwitches];
      availableAccelerators = new double[numberOfpSwitches];
      totalAccelerators = new double[numberOfpSwitches];
      availableStorage = new double[numberOfpSwitches];
      totalStorage = new double[numberOfpSwitches];
      SIs = new double[numberOfpSwitches];

      for (int i = 0; i < numberOfpSwitches; i++) {
        availableProcessors[i] = t.getAvailableProcessors()[i];
        totalProcessors[i] = t.getTotalProcessors()[i];
        availableMemory[i] = t.getAvailableMemory()[i];
        totalMemory[i] = t.getTotalMemory()[i];
        availableStorage[i] = t.getAvailableStorage()[i];
        totalStorage[i] = t.getTotalStorage()[i];
        availableAccelerators[i] = t.getAvailableAccelerators()[i];
        totalAccelerators[i] = t.getTotalAccelerators()[i];
        SIs[i] = t.gSIs()[i];
      }

      SI = t.gSI();
      sPMSA = new double[8];

      for (int i = 0; i < 8; i++) {
        sPMSA[i] = t.gsPMSA()[i];
      }

      C = t.gC();
      P = t.gP();
      Pi = t.gPi();
    }
  }
  return *this;
}

pRouter::~pRouter()
{
  if (alloc) {
    alloc = 0;
    numberOfpSwitches = 0;
    numberOfFunctions = 0;
    pollIntervalpRouter = 0.0;
    pSwitches->clear();
    delete[] pSwitches;
    pSwitches = nullptr;
    delete[] Ws;
    delete[] Fs;
    Ws = nullptr;
    Fs = nullptr;
    delete[] availableProcessors;
    delete[] totalProcessors;
    delete[] availableMemory;
    delete[] totalMemory;
    delete[] availableStorage;
    delete[] totalStorage;
    delete[] availableAccelerators;
    delete[] totalAccelerators;
    delete[] sPMSA;
    delete[] SIs;
    availableProcessors = nullptr;
    totalProcessors = nullptr;
    availableMemory = nullptr;
    totalMemory = nullptr;
    availableAccelerators = nullptr;
    totalAccelerators = nullptr;
    availableStorage = nullptr;
    totalStorage = nullptr;
    sPMSA = nullptr;
    SIs = nullptr;
    SI = 0.0;
    C = 0.0;
    P = 0.0;
    Pi = 0.0;
  }
}

void pRouter::computeFs()
{
  if (alloc) {
    for (int i = 0; i < numberOfFunctions; i++) {
      Fs[i] = 0.0;
    }
    list<pSwitch*>::iterator it = pSwitches->begin();
    for (int i = 0; i < numberOfpSwitches; i++) {
      for (int j = 0; j < numberOfFunctions; j++) {
        Fs[j] += (*it)->gFs()[j];
      }
      it++;
    }
    for (int j = 0; j < numberOfFunctions; j++) {
      Fs[j] /= numberOfpSwitches;
    }
  }
}

void pRouter::computeSI()
{
  if (alloc) {
    SI = (1e-4) * (((double)rand()) / RAND_MAX);
    for (int i = 0; i < numberOfFunctions; i++) {
      SI += Ws[i] * Fs[i];
    }
  }
}

int pRouter::probe(const double& Proc, const double& Mem, const double& Sto, const int& Acc)
{
  if (Proc <= sPMSA[0] && Mem <= sPMSA[2] && Sto <= sPMSA[4] && Acc <= (int)sPMSA[6]) {
    return 1;
  } else {
    return -1;
  }
}

void pRouter::deploy(resource** resources, netw* network, stat* stats, task& task_)
{
  if (alloc) {
    int L_numberOfVMs = task_.getNumberOfVMs();
    double* L_reqPMNS = task_.greqPMNS();
    double maxSI = 0.0, ssum = 0.0, *maxSIs;
    int L_avAcc = task_.gavAcc()[0];
    int omp_thr = atoi(getenv("OMP_NUM_THREADS"));
    int i, tid;
    int choice = -1, *choices;

    list<pSwitch*>::iterator itf;
    choices = new int[omp_thr];
    maxSIs = new double[omp_thr];

    // Initialize the choices and maxSIs arrays
    for (i = 0; i < omp_thr; i++) {
      choices[i] = -1;
      maxSIs[i] = 0.0;
    }

    // The total requested units = number of task VMs * task units
    double reqProc = L_numberOfVMs * L_reqPMNS[0];
    double reqMem = L_numberOfVMs * L_reqPMNS[1];
    double reqSto = L_numberOfVMs * L_reqPMNS[3];
    int reqAcc = L_numberOfVMs * L_avAcc;

// Partition the list of pSwithces per thread. Then, for every pSwitch, determine their maximum SI. Set the
// index number of the switch with the maximum SI in the array choices[thread_number]
#pragma omp parallel default(shared) private(i, tid) num_threads(omp_thr)
    {
      tid = omp_get_thread_num();
#pragma omp for
      for (i = 0; i < numberOfpSwitches; i++) {
        if (maxSIs[tid] < SIs[i] && reqProc <= availableProcessors[i] && reqMem <= availableMemory[i] &&
            reqSto <= availableStorage[i] && reqAcc <= availableAccelerators[i]) {
          maxSIs[tid] = SIs[i];
          choices[tid] = i;
        }
      }
    }

    choice = choices[0];
    maxSI = maxSIs[0];

    // From the list of pSwitches with the maximum SI per thread, selected the highest for the task deployment
    for (i = 1; i < omp_thr; i++) {
      if (maxSI < maxSIs[i]) {
        maxSI = maxSIs[i];
        choice = choices[i];
      }
    }
    delete[] choices;
    delete[] maxSIs;

    // If no suitable pSwitch was found, reject the task
    if (choice == -1) {
      stats[task_.getAvailableImplementations()[0]].rejectedTasks++;
      return;
    }

    // Set the pSwitch iterator to the position of the selected pSwitch
    itf = pSwitches->begin();
    for (i = 0; i < choice; i++) {
      itf++;
    }

    // Reduce the task's units from the list of available units
    availableProcessors[choice] -= reqProc;
    availableMemory[choice] -= reqMem;
    availableStorage[choice] -= reqSto;
    availableAccelerators[choice] -= reqAcc;
    sPMSA[0] -= reqProc;
    sPMSA[2] -= reqMem;
    sPMSA[4] -= reqSto;
    sPMSA[6] -= (double)reqAcc;

    // Re-compute the assessment functions
    if (reqAcc > 0) {
      for (i = 0; i < 4; i++) {
        ssum += Ws[i] * deassessmentFunctions(-(double)reqAcc, sPMSA[7], -reqMem, sPMSA[3], i);
      }
    } else {
      for (i = 0; i < 4; i++) {
        ssum += Ws[i] * deassessmentFunctions(-reqProc, sPMSA[1], -reqMem, sPMSA[3], i);
      }
    }
    SI += ssum;
    ssum = 0.0;
    if (sPMSA[7] > 0) {
      for (i = 0; i < 4; i++) {
        ssum += Ws[i] * deassessmentFunctions(-(double)reqAcc, (double)totalAccelerators[choice], -reqMem,
                                              totalMemory[choice], i);
      }
    } else {
      for (i = 0; i < 4; i++) {
        ssum += Ws[i] * deassessmentFunctions(-reqProc, totalProcessors[choice], -reqMem, totalMemory[choice], i);
      }
    }
    // Re-compute the suitability index
    SIs[choice] += ssum;

    // Call the pSwitch::deploy method recursively
    (*itf)->deploy(resources, network, stats, task_);
  }
}

void pRouter::updateStateInfo(const double& tstep)
{
  int i;
  if (alloc) {
    // On every time interval
    if (((int)tstep % (int)pollIntervalpRouter) == 0) {
      int len = (int)pSwitches->size();

      for (i = 0; i < 8; i++) {
        sPMSA[i] = 0.0;
      }

      i = 0;
      for (list<pSwitch*>::iterator it = pSwitches->begin(); it != pSwitches->end(); it++) {
        availableProcessors[i] = ((*it)->gsPMSA())[0];
        totalProcessors[i] = ((*it)->gsPMSA())[1];
        availableMemory[i] = ((*it)->gsPMSA())[2];
        totalMemory[i] = ((*it)->gsPMSA())[3];
        availableStorage[i] = ((*it)->gsPMSA())[4];
        totalStorage[i] = ((*it)->gsPMSA())[5];
        availableAccelerators[i] = ((*it)->gsPMSA())[6];
        totalAccelerators[i] = ((*it)->gsPMSA())[7];
        SIs[i] = ((*it)->gSI());
        i++;
      }

      for (i = 0; i < len; i++) {
        sPMSA[0] += availableProcessors[i];
        sPMSA[1] += totalProcessors[i];
        sPMSA[2] += availableMemory[i];
        sPMSA[3] += totalMemory[i];
        sPMSA[4] += availableStorage[i];
        sPMSA[5] += totalStorage[i];
        sPMSA[6] += availableAccelerators[i];
        sPMSA[7] += totalAccelerators[i];
      }
      computeFs();
      computeSI();
    }
  }
}

// Derivatives used for Taylor expansions (principal linear part)
double pRouter::deassessmentFunctions(const double& dNu, const double& totNu, const double& dNmem,
                                      const double& totalMemory, const int& choice)
{
  switch (choice) {
    case 0:
      return (dNu * C / totNu);
      break;
    case 1:
      return (dNmem / totalMemory);
      break;
    case 2:
      return (dNu * Pi * P * totNu) / ((P * (totNu - dNu) + Pi * dNu) * (P * (totNu - dNu) + Pi * dNu));
      break;
    case 3:
      return (0.2 * dNu / (totNu));
      break;
    default:
      return 0.0;
  }
}

int pRouter::galloc() const { return alloc; }
int pRouter::getNumberOfpSwitches() const { return numberOfpSwitches; }
int pRouter::getNumberOfFunctions() const { return numberOfFunctions; }
double pRouter::gpollIntervalpRouter() const { return pollIntervalpRouter; }
list<pSwitch*>* pRouter::gpSwitches() const { return pSwitches; }
double* pRouter::getAvailableProcessors() const { return availableProcessors; }
double* pRouter::getTotalProcessors() const { return totalProcessors; }
double* pRouter::getAvailableMemory() const { return availableMemory; }
double* pRouter::getTotalMemory() const { return totalMemory; }
double* pRouter::getAvailableAccelerators() const { return availableAccelerators; }
double* pRouter::getTotalAccelerators() const { return totalAccelerators; }
double* pRouter::getAvailableStorage() const { return availableStorage; }
double* pRouter::getTotalStorage() const { return totalStorage; }
double* pRouter::gsPMSA() const { return sPMSA; }
double* pRouter::gWs() const { return Ws; }
double* pRouter::gFs() const { return Fs; }
double pRouter::gSI() const { return SI; }
double* pRouter::gSIs() const { return SIs; }
double pRouter::gC() const { return C; }
double pRouter::gP() const { return P; }
double pRouter::gPi() const { return Pi; }
void pRouter::print() const
{
  if (alloc) {
    cout << "pRouter Stats" << endl;
    cout << "Available Processing Units: " << sPMSA[0] << endl;
    cout << "Total Processing Units: " << sPMSA[1] << endl;
    cout << "Available Memory: " << sPMSA[2] << endl;
    cout << "Total Memory: " << sPMSA[3] << endl;
    cout << "Available Storage: " << sPMSA[4] << endl;
    cout << "Total Storage: " << sPMSA[5] << endl;
    cout << "Available Accelerators: " << sPMSA[6] << endl;
    cout << "Total Accelerators: " << sPMSA[7] << endl;
  }
}
