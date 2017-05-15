#include <iostream>
#include "printme.hh"

printme::printme()
{
  fName = "me";
}

void printme::printname()
{
  printf("%s\n", fName);
}

void printme::setname(char *nm)
{
  fName = nm; 
}
