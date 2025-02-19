#include <task.h>
#include <cstdlib>  // for rand, RAND_MAX
#include <iostream> // for cout, endl

using std::cout;
using std::endl;

task::task()
  : type(-1),
    numberOfAvailableImplementations(0),
    availableImplementations(nullptr),
    requestedInstructions(0.0),
    numberOfVMs(0),
    reqPMNS(nullptr),
    typeactPMN(nullptr),
    minmaxactPMN(nullptr),
    rhoAcc(nullptr),
    avAcc(nullptr),
    alloc(0),
    resourceIDs(nullptr),
    cUtilPMNr(nullptr)
{
}

task::task(const task& t)
{
  int i, j;
  alloc = t.galloc();
  if (alloc) {
    type = t.getType();
    numberOfAvailableImplementations = t.getNumberOfAvailableImplementations();
    availableImplementations = new int[numberOfAvailableImplementations];
    for (i = 0; i < numberOfAvailableImplementations; i++) {
      availableImplementations[i] = t.getAvailableImplementations()[i];
    }
    requestedInstructions = t.grequestedInstructions();
    numberOfVMs = t.getNumberOfVMs();
    reqPMNS = new double[4];
    for (i = 0; i < 4; i++) {
      reqPMNS[i] = t.greqPMNS()[i];
    }
    typeactPMN = new int[3];
    for (i = 0; i < 3; i++) {
      typeactPMN[i] = t.getTypeactPMN()[i];
    }
    minmaxactPMN = new double*[3];
    for (i = 0; i < 3; i++) {
      minmaxactPMN[i] = new double[2];
    }
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 2; j++) {
        minmaxactPMN[i][j] = t.gminmaxactPMN()[i][j];
      }
    }
    rhoAcc = new double[numberOfAvailableImplementations];
    avAcc = new int[numberOfAvailableImplementations];
    for (i = 0; i < numberOfAvailableImplementations; i++) {
      rhoAcc[i] = t.grhoAcc()[i];
      avAcc[i] = t.gavAcc()[i];
    }
    if (t.gresourceIDs() != nullptr) {
      resourceIDs = new int[numberOfVMs];
      for (i = 0; i < numberOfVMs; i++)
        resourceIDs[i] = t.gresourceIDs()[i];
    } else
      resourceIDs = nullptr;

    cUtilPMNr = new double[4];
    for (i = 0; i < 4; i++) {
      cUtilPMNr[i] = t.gcUtilPMNr()[i];
    }
  }
}

task::task(const int& L_type, const int& L_numberOfAvailableImplementations, const int* L_availableImplementations,
           const double& L_requestedInstructions, const int& L_numberOfVMs, const double& L_reqP, const double& L_reqM,
           const double& L_reqN, const double& L_reqS, const int& L_typeactP, const int& L_typeactM,
           const int& L_typeactN, const double* L_minmaxactP, const double* L_minmaxactM, const double* L_minmaxactN,
           const int* L_avAcc, const double* L_rhoAcc)
{
  alloc = 1;
  type = L_type;
  numberOfAvailableImplementations = L_numberOfAvailableImplementations;
  availableImplementations = new int[numberOfAvailableImplementations];

  for (int i = 0; i < numberOfAvailableImplementations; i++) {
    availableImplementations[i] = L_availableImplementations[i];
  }

  requestedInstructions = L_requestedInstructions;
  numberOfVMs = L_numberOfVMs;
  reqPMNS = new double[4];
  reqPMNS[0] = L_reqP;
  reqPMNS[1] = L_reqM;
  reqPMNS[2] = L_reqN;
  reqPMNS[3] = L_reqS;
  typeactPMN = new int[3];
  typeactPMN[0] = L_typeactP;
  typeactPMN[1] = L_typeactM;
  typeactPMN[2] = L_typeactN;
  minmaxactPMN = new double*[3];

  for (int i = 0; i < 3; i++) {
    minmaxactPMN[i] = new double[2];
  }

  for (int i = 0; i < 2; i++) {
    minmaxactPMN[0][i] = L_minmaxactP[i];
  }

  for (int i = 0; i < 2; i++) {
    minmaxactPMN[1][i] = L_minmaxactM[i];
  }

  for (int i = 0; i < 2; i++) {
    minmaxactPMN[2][i] = L_minmaxactN[i];
  }

  rhoAcc = new double[numberOfAvailableImplementations];
  avAcc = new int[numberOfAvailableImplementations];

  for (int i = 0; i < numberOfAvailableImplementations; i++) {
    avAcc[i] = L_avAcc[i];
    rhoAcc[i] = L_rhoAcc[i];
  }

  resourceIDs = nullptr;

  cUtilPMNr = new double[4];
  for (int i = 0; i < 4; i++) {
    cUtilPMNr[i] = 0.0;
  }
  // print();
}

task& task::operator=(const task& t)
{
  int i, j;
  if (this != &t) {
    if (alloc) {
      numberOfAvailableImplementations = 0;
      type = -1;
      delete[] availableImplementations;
      availableImplementations = nullptr;
      requestedInstructions = 0.0;
      numberOfVMs = 0;
      delete[] reqPMNS;
      reqPMNS = nullptr;
      delete[] typeactPMN;
      typeactPMN = nullptr;
      for (i = 0; i < 3; i++) {
        delete[] minmaxactPMN[i];
      }
      delete[] minmaxactPMN;
      minmaxactPMN = nullptr;
      delete[] rhoAcc;
      delete[] avAcc;
      rhoAcc = nullptr;
      avAcc = nullptr;
      alloc = 0;
      if (resourceIDs != nullptr) {
        delete[] resourceIDs;
      }
      delete[] cUtilPMNr;
      cUtilPMNr = nullptr;
    }
    alloc = t.galloc();
    if (alloc) {
      type = t.getType();
      numberOfAvailableImplementations = t.getNumberOfAvailableImplementations();
      availableImplementations = new int[numberOfAvailableImplementations];
      for (i = 0; i < numberOfAvailableImplementations; i++) {
        availableImplementations[i] = t.getAvailableImplementations()[i];
      }
      requestedInstructions = t.grequestedInstructions();
      numberOfVMs = t.getNumberOfVMs();
      reqPMNS = new double[4];
      for (i = 0; i < 4; i++) {
        reqPMNS[i] = t.greqPMNS()[i];
      }
      typeactPMN = new int[3];
      for (i = 0; i < 3; i++) {
        typeactPMN[i] = t.getTypeactPMN()[i];
      }
      minmaxactPMN = new double*[3];
      for (i = 0; i < 3; i++) {
        minmaxactPMN[i] = new double[2];
      }
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 2; j++) {
          minmaxactPMN[i][j] = t.gminmaxactPMN()[i][j];
        }
      }
      rhoAcc = new double[numberOfAvailableImplementations];
      avAcc = new int[numberOfAvailableImplementations];
      for (i = 0; i < numberOfAvailableImplementations; i++) {
        rhoAcc[i] = t.grhoAcc()[i];
        avAcc[i] = t.gavAcc()[i];
      }
      if (t.gresourceIDs() != nullptr) {
        resourceIDs = new int[numberOfVMs];
        for (i = 0; i < numberOfVMs; i++)
          resourceIDs[i] = t.gresourceIDs()[i];
      } else {
        resourceIDs = nullptr;
      }
      cUtilPMNr = new double[4];
      for (i = 0; i < 4; i++) {
        cUtilPMNr[i] = t.gcUtilPMNr()[i];
      }
    }
  }
  return *this;
}

task::~task()
{
  if (alloc) {
    type = -1;
    numberOfAvailableImplementations = 0;
    delete[] availableImplementations;
    availableImplementations = nullptr;
    requestedInstructions = 0.0;
    numberOfVMs = 0;
    delete[] reqPMNS;
    reqPMNS = nullptr;
    delete[] typeactPMN;
    typeactPMN = nullptr;

    for (int i = 0; i < 3; i++) {
      delete[] minmaxactPMN[i];
    }

    delete[] minmaxactPMN;
    minmaxactPMN = nullptr;
    delete[] avAcc;
    delete[] rhoAcc;
    rhoAcc = nullptr;
    avAcc = nullptr;
    alloc = 0;

    if (resourceIDs != nullptr) {
      delete[] resourceIDs;
    }

    resourceIDs = nullptr;
    delete[] cUtilPMNr;
    cUtilPMNr = nullptr;
  }
}

void task::reduceIns(const double& amount) { requestedInstructions -= amount; }
double task::getactP()
{
  double r = ((double)rand() / RAND_MAX);
  switch (typeactPMN[0]) {
    case 1:
      return (minmaxactPMN[0][0] + (minmaxactPMN[0][1] - minmaxactPMN[0][0]) * r);
      break;
    case 2:
      return 0;
      break;
    default:
      return 0;
      break;
  }
}

double task::getactM()
{
  double r = ((double)rand() / RAND_MAX);
  switch (typeactPMN[1]) {
    case 1:
      return (minmaxactPMN[1][0] + (minmaxactPMN[1][1] - minmaxactPMN[1][0]) * r);
      break;
    case 2:
      return 0;
      break;
    default:
      return 0;
      break;
  }
}

double task::getactN()
{
  double r = ((double)rand() / RAND_MAX);
  switch (typeactPMN[2]) {
    case 1:
      return (minmaxactPMN[2][0] + (minmaxactPMN[2][1] - minmaxactPMN[2][0]) * r);
      break;
    case 2:
      return 0;
      break;
    default:
      return 0;
      break;
  }
}

int task::getType() const { return type; }
int task::getNumberOfAvailableImplementations() const { return numberOfAvailableImplementations; }
int* task::getAvailableImplementations() const { return availableImplementations; }
double task::grequestedInstructions() const { return requestedInstructions; }
int task::getNumberOfVMs() const { return numberOfVMs; }
double* task::greqPMNS() const { return reqPMNS; }
int* task::getTypeactPMN() const { return typeactPMN; }
double** task::gminmaxactPMN() const { return minmaxactPMN; }
int* task::gavAcc() const { return avAcc; }
double* task::grhoAcc() const { return rhoAcc; }
int task::galloc() const { return alloc; }
int* task::gresourceIDs() const { return resourceIDs; }
double* task::gcUtilPMNr() const { return cUtilPMNr; }
void task::compcUtilPMNr()
{
  if (alloc) {
    cUtilPMNr[0] = getactP() * reqPMNS[0];
    cUtilPMNr[1] = getactM() * reqPMNS[1];
    cUtilPMNr[2] = getactN() * reqPMNS[2];
    cUtilPMNr[3] = rhoAcc[0] * ((double)avAcc[0]);
  }
}

void task::attachResources(const int* IDs)
{
  int i;
  if (alloc) {
    if (resourceIDs == nullptr) {
      resourceIDs = new int[numberOfVMs];
      for (i = 0; i < numberOfVMs; i++) {
        resourceIDs[i] = IDs[i];
      }
    }
  }
}

void task::detachResources()
{
  if (alloc) {
    if (resourceIDs != nullptr) {
      delete[] resourceIDs;
      resourceIDs = nullptr;
    }
  }
}

void task::remapType(const int* type, const int& num)
{
  if (alloc) {
    for (int i = 0; i < num; i++)
      availableImplementations[i] = type[i];
  }
}

void task::reduceImpl(const int* type)
{
  if (alloc) {
    numberOfAvailableImplementations = 1;
    delete[] availableImplementations;
    availableImplementations = new int[1];
    availableImplementations[0] = *type;
    int L_avAcc = avAcc[*type];
    delete[] avAcc;
    avAcc = new int[1];
    avAcc[0] = L_avAcc;
    double L_rhoAcc = rhoAcc[*type];
    delete[] rhoAcc;
    rhoAcc = new double[1];
    rhoAcc[0] = L_rhoAcc;
  }
}

void task::print() const
{
  if (alloc) {
    cout << "-----------------------------------------------" << endl;
    cout << "Task type (oil, genomics etc): " << type << endl;
    cout << "Number of available implementations: " << numberOfAvailableImplementations << endl;
    cout << "Available Implementations: ";
    for (int i = 0; i < numberOfAvailableImplementations; i++) {
      cout << availableImplementations[i] << " ";
    }
    cout << endl;
    cout << "Number of Instructions: " << requestedInstructions << endl;
    cout << "Number of VMs: " << numberOfVMs << endl;
    cout << "vCPUs per VM: " << reqPMNS[0] << endl;
    cout << "Memory per VM: " << reqPMNS[1] << " GBytes" << endl;
    cout << "Storage per VM: " << reqPMNS[3] << " TBytes" << endl;
    cout << "Network per App: " << reqPMNS[2] << " Gbps" << endl;
    cout << "Type of Actual Utilization (Proc,Mem,Sto): ";
    for (int i = 0; i < 3; i++) {
      cout << typeactPMN[i] << " ";
    }
    cout << endl;
    cout << "Minimum - Maximum actual utilization Processors: " << minmaxactPMN[0][0] << " " << minmaxactPMN[0][1]
         << endl;
    cout << "Minimum - Maximum actual utilization Memory: " << minmaxactPMN[1][0] << " " << minmaxactPMN[1][1] << endl;
    cout << "Minimum - Maximum actual utilization Network: " << minmaxactPMN[2][0] << " " << minmaxactPMN[2][1] << endl;
    cout << "Accelerator support per Implementation: ";
    for (int i = 0; i < numberOfAvailableImplementations; i++) {
      cout << avAcc[i] << " ";
    }
    cout << endl;
    cout << "Actual accelerator usage: ";
    for (int i = 0; i < numberOfAvailableImplementations; i++) {
      cout << rhoAcc[i] << " ";
    }
    cout << endl;
    cout << "-----------------------------------------------" << endl;
  }
}
