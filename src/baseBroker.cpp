#include <baseBroker.h>
#include <cell.h>

baseBroker::baseBroker()
  : alloc(0),
    numberOfTypes(0),
    types(nullptr),
    numberOfResourcesPerType(nullptr),
    availableNetwork(0.0),
    totalNetwork(0.0)
{
}

baseBroker::~baseBroker() {}

void baseBroker::init(const cell* clCell)
{
  alloc = 1;
  numberOfTypes = clCell->getNumberOfTypes();
  types = new int[numberOfTypes];
  numberOfResourcesPerType = new int[numberOfTypes];

  for (int i = 0; i < numberOfTypes; i++) {
    types[i] = clCell->getTypes()[i];
    numberOfResourcesPerType[i] = clCell->getNumberOfResourcesPerType()[i];
  }
}

int baseBroker::galloc() const { return alloc; }
int baseBroker::getNumberOfTypes() const { return numberOfTypes; }
int* baseBroker::getTypes() const { return types; }
int* baseBroker::getNumberOfResourcesPerType() const { return numberOfResourcesPerType; }
double baseBroker::getTotalNetwork() const { return totalNetwork; }
double baseBroker::getAvailableNetwork() const { return availableNetwork; }
