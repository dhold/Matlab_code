#ifndef PhiParameters_h
#define PhiParameters_h
#include <string>
#include <iostream>
#include "ComplexMatrix.h"

#define PARAMPOSITIVE 0
#define PARAMNONNEGATIVE 1
#define PARAMNONE 2
#define PARAMBOOLEAN 3

using namespace std;

const std::string NUMS = "1234567890";
const std::string ARRAY = ".Ee-+,";
const std::string CARRAY = ".Ee-+,()";

// vector
bool checkComplexEntryString(const std::string &);
bool checkEntryString(const std::string &);

// scalar 
bool checkIntString(const std::string &);
bool checkFloatString(const std::string &);
bool checkBoolString(const std::string &);
bool checkComplexString(const std::string &);



bool getVal(const string &, const string &, Float &, const int);
bool getVal(const string &, const string &, int &, const int);
bool getVal(const string &, const string &, bool &);
bool getVal(const string &, const string &, Complex &);

//
// Reads string line of Float or Int values
//
bool getLineOfVal(const string &, int *, const int, const int );
bool getLineOfVal(const string &, Float *, const int, const int );
bool getLineOfVal(const string &, Complex *, const int, const int);
 
#endif
