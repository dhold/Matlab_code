#include "PhiParameters.h"


template <class T>
bool readVal(const string & line, const string & varname, T & var,const int range) {
  stringstream strStream;
  strStream.str(line);
  strStream >> var;
  switch (range) {
    case (PARAMPOSITIVE):
      if ( var <= 0 ){
        cerr << "getVal() Error: " << varname << " must be positive.\n";
        return false;
      }
      break;
    case (PARAMNONNEGATIVE):
      if ( var < 0 ) {
        cerr << "getVal() Error: " << varname << " must be non-negative.\n";
        return false;
      }
      break;
    case (PARAMBOOLEAN):
      if ( var != 0 && var != 1 ) {
        cerr << "getVal() Error: " << varname << " must be 0 (for 'no' or 'false') or 1 (for 'yes' or 'true').\n";
        return false;
      }
    default:
      break;
  }
  return true;
}


//Read integer from line
bool getVal(const string & line, const string & varname, int & var,const int range) {
  //cout << "Getting (int): " << line << endl;
  bool isok = true;
  isok = checkIntString(line);
  if (!isok) {
    cerr << "getVal() Error reading " << varname << "\n";
    return false;
  }
  isok = readVal<int>(line, varname, var, range);
  return isok;
}

//Read Float from line
bool getVal(const string & line, const string & varname, Float & var,const int range) {
  bool isok = true;
  isok = checkFloatString(line);
  if (!isok) {
    cerr << "getVal() Error reading " << varname << "\n";
    return false;
  }
  isok = readVal<Float>(line, varname, var, range);
  return isok;
}

//Read bool from line
bool getVal(const string & line, const string & varname, bool & var) {
  bool isok = true;
  isok = checkBoolString(line);
  if (!isok) {
    cerr << "getVal() Error reading " << varname << ".\nBoolean parameters must be 0 (for no/false) or 1 (for yes/true).\n";
    return false;
  }
  int v;
  isok = readVal<int>(line, varname, v, PARAMBOOLEAN);
  if (isok) 
    var = ( v == 1 ) ? true : false;
  return isok;
}

//Read Complex from line
bool getVal(const string & line, const string & varname, Complex & var) {
  bool isok = true;
  isok = checkComplexString(line);
  if (!isok) {
    cerr << "getVal() Error reading " << varname << ".\n";
    return false;
  }

  int ilb = line.find("(",0);
  int i = line.find(",",ilb+1);
  int irb = line.find(")",i+1);

  Float bufre = 0.;
  Float bufim = 0.;

  isok = readVal<Float>(line.substr(ilb+1,i-ilb-1), varname, bufre, PARAMNONE);
  if (!isok)
    return false;
  isok = readVal<Float>(line.substr(i+1,irb-i-1), varname, bufim, PARAMNONE);
  if (!isok)
    return false;
  var = Complex(bufre,bufim);  
  return isok;
}


bool getLineOfVal(const string & line, int * vec, const int numVals, const int jskip ) {
  if (!checkEntryString(line))
    return false;

  int i = 0;
  int j = 0;
  int n = 0;
  int numComma = 0;
  if ((line.compare(line.size()-1,1,",") == 0)) {
      cerr << "getLineOfVal() Error: should end with , in " << line << endl;
      return false;
  }

  i = line.find(",",0);
  while(i != string::npos) {
    numComma++;
    i = line.find(",",i+1);
  }

  if (numComma != numVals-1) {
    cerr << "getLineOfVal() Error: incorrect number of ,s in " << line << endl;
    return false;
  }


  int iprev = 0;
  stringstream ss;
  int buf;
  j = 0;
  i = line.find(",",0);
  while ( i != string::npos ) {
      if (i>iprev) {
        ss.str(line.substr(iprev,i-iprev));
        ss >> buf;
        //debug
        //cout << "j= " << j << " buf= " << buf << endl;
        //debug
        vec[j] = buf;
        ss.clear();
        iprev = i+1;
        i = line.find(",",iprev);
        n += 1;
        j += jskip;
      }
      else {
        cerr << "getLineOfVal() Error: no value to unpack at element " << n << " in " << line << endl;
        return false;
      }
  }
  ss.str(line.substr(iprev,line.length()-iprev));
  ss >> buf;
  vec[j] = buf;
  ss.clear();
  n += 1;

  if (n != numVals) {
    cerr << "getLineOfVal() Error: Expected to read " << numVals << " but got " << n << endl;
    return false;
  }

  return true;
};


//
// Reads in a line of non-complex values
//
bool getLineOfVal(const string & line, Float * vec, const int numVals, const int jskip ) {
  if (!checkEntryString(line))
      return false;

  int i = 0;
  int j = 0;
  int n = 0;
  int numComma = 0;
  if ((line.compare(line.size()-1,1,",") == 0)) {
      cerr << "getLineOfVal() Error: should end with , in " << line << endl;
      return false;
  }

  i = line.find(",",0);
  while(i != string::npos) {
    numComma++;
    i = line.find(",",i+1);
  }

  if (numComma != numVals-1) {
    cerr << "getLineOfVal() Error: incorrect number of ,s in " << line << endl;
    return false;
  }


  int iprev = 0;
  stringstream ss;
  Float buf;
  j = 0;
  i = line.find(",",0);
  while ( i != string::npos ) {
      if (i>iprev) {
        ss.str(line.substr(iprev,i-iprev));
        ss >> buf;
        //debug
        //cout << "j= " << j << " buf= " << buf << endl;
        //debug
        vec[j] = buf;
        ss.clear();
        iprev = i+1;
        i = line.find(",",iprev);
        n += 1;
        j += jskip;
      }
      else {
        cerr << "getLineOfVal() Error: no value to unpack at element " << n << " in " << line << endl;
        return false;
      }
  }
  ss.str(line.substr(iprev,line.length()-iprev));
  ss >> buf;
  vec[j] = buf;
  ss.clear();
  n += 1;

  if (n != numVals) {
    cerr << "getLineOfVal() Error: Expected to read " << numVals << " but got " << n << endl;
    return false;
  }

  return true;
};


bool getLineOfVal(const string & line, Complex * vec, const int numVals, const int jskip ) {
  int i = 0;
  int j = 0;
  int n = 0;
  int numLeftBrac = 0;
  int numRightBrac = 0;
  int numComma = 0;
  if ((line.compare(line.size()-1,1,",") == 0)) {
      cerr << "getLineOfComplex() Error: should end with , in " << line << endl;
      return false;
  }

  i = line.find(",",0);
  while(i != string::npos) {
    numComma++;
    i = line.find(",",i+1);
  }

  if (numComma != 2*numVals - 1 && numComma != numVals-1) {
    cerr << "getLineOfComplex() Error: incorrect number of ,s in " << line << endl;
    return false;
  }

  i = line.find("(",0);
  while(i != string::npos) {
    numLeftBrac++;
    i = line.find("(",i+1);
  }

  if (numLeftBrac != 0 && numLeftBrac != numVals){
    cerr << "getLineOfComplex() Error: incorrect number of ('s in " << line << endl;
    return false;
  }

  i = line.find(")",0);
  while(i != string::npos) {
    numRightBrac++;
    i = line.find(")",i+1);
  }

  if (numRightBrac != 0 && numRightBrac != numVals) {
    cerr << "getLineOfComplex() Error: incorrect number of )'s in " << line << endl;
    return false;
  }

  if (numRightBrac == 0) {
    if (!checkEntryString(line))
        return false;
  }
  else {
    if (!checkComplexEntryString(line))
        return false;
  }

  if (numRightBrac == 0 && numLeftBrac == 0 && numComma == numVals-1){
    int iprev = 0;
    stringstream ss;
    Float buf;
    j = 0;
    i = line.find(",",0);
    while ( i != string::npos ) {
        if (i>iprev) {
          ss.str(line.substr(iprev,i-iprev));
          ss >> buf;
          //debug
          //cout << "j= " << j << " buf= " << buf << endl;
          //debug
          vec[j] = Complex(buf);
          ss.clear();
          iprev = i+1;
          i = line.find(",",iprev);
          n += 1;
          j += jskip;
        }
        else {
          cerr << "getLineOfComplex() Error: no value to unpack at element " << n << " in " << line << endl;
          return false;
        }
    }
    ss.str(line.substr(iprev,line.length()-iprev));
    ss >> buf;
    vec[j] = Complex(buf);
    ss.clear();
    n += 1;
  }
  else if (numRightBrac == numVals && numLeftBrac == numVals && numComma == 2*numVals-1) {
    stringstream ss;
    Float buf_re;
    Float buf_im;
    j = 0;
    int ilb = 0; //line.find("(",0);
    i = line.find(",",ilb+1);
    int irb = line.find(",",i+1);
    while (ilb != string::npos && i != string::npos) {
        ilb = line.find("(",ilb);
        irb = line.rfind(")",irb);
        //cout << "ilb->i: " << line.substr(ilb+1,i-ilb-1) << " i->irb: " << line.substr(i+1,irb-i-1) << endl;
        if ( i-ilb > 1 &&  irb-i > 1) {
          ss.str(line.substr(ilb+1,i-ilb-1));
          ss >> buf_re;
          ss.clear();
          ss.str(line.substr(i+1,irb-i-1));
          ss >> buf_im;
          ss.clear();
          vec[j] = Complex(buf_re,buf_im);
          ilb = line.find(",",irb);
          i   = line.find(",",ilb+1);
          irb = line.find(",",i+1);
          n += 1;
          j += jskip;
        }
        else {
          cerr << "getLineOfComplex() Error: no value to unpack at element " << n << " in " << line << endl;
          return false;
        }
    } 
  }
  else {
      cerr << "getLineOfComplex() Error: real and complex input should not be mixed in " << line << endl;
      return false;
  }

  if (n != numVals) {
    cerr << "getLineOfComplexError: Expected to read " << numVals << " entries but only got " << n << " in " << line << endl;
    return false;
  }

  return true;

} 

bool checkBoolString(const std::string &s) {
    bool isnum = true;
    int i = 0;
    bool whitespace = false;
    //skip initial whitespace
    while (s.at(i) == ' ')
      ++i;
    if (s.at(i) != '0' && s.at(i) != '1')
      isnum = false;
    ++i;
    while (isnum && i < s.size()) {
      if (s.at(i) != ' ')
        isnum = false;
      ++i;
    }
    if (!isnum) {
      stringstream ss;
      for (int j = 0; j < i-1; ++j)
        ss << " ";
      ss << "^";
      cerr << "checkBoolString() Error at position " << i-1 << " in: \n" << s << "\n" << ss.str() << "\n";
    }
    return isnum;
}

bool checkIntString(const std::string &s) {
    bool isnum = true;
    int i = 0;
    bool whitespace = false;
    //skip initial whitespace
    while (s.at(i) == ' ')
      ++i;
    while (isnum && i < s.size()) {
      if (NUMS.find(s.at(i)) == string::npos || (NUMS.find(s.at(i)) == string::npos && whitespace))
        isnum = false;
      ++i;
      while ( i < s.size() && s.at(i) == ' ') {
         ++i;
         whitespace = true;
      }
    }
    if (!isnum) {
      stringstream ss;
      for (int j = 0; j < i-1; ++j)
        ss << " ";
      ss << "^";
      cerr << "checkIntString() Error at position " << i-1 << " in: \n" << s << "\n" << ss.str() << "\n";
    }
    return isnum;
}

bool checkFloatString(const std::string &s) {
    bool isnum = true;
    bool period = false;
    bool e = false;
    bool whitespace = false;
    int numi;
    int modi;
    int m = -1;
    int i = 0;
    //skip initial whitespace
    while (s.at(i) == ' ')
      ++i;
    while (isnum && i < s.size()) {
      ++m;
      numi = NUMS.find(s.at(i));
      modi = ARRAY.find(s.at(i));
      if (numi == string::npos && modi == string::npos)
         isnum = false;
      else if ( numi != string::npos  ){
         if (whitespace) 
            isnum = false;
      }
      else if ( modi != string::npos ){
         switch (modi) {
           case 0: // '.'
             if (period || e || whitespace) 
               isnum = false;
             period = true;
             break;
           case 1: // 'E'
           case 2: // 'e'
             if (m == 0 || e || whitespace)
               isnum = false;
             e = true;
             break;
           case 3: // '-'
           case 4: // '+'
             if ((m != 0 && s.at(i-1) != 'e' && s.at(i-1) != 'E' ) || whitespace )
               isnum = false;
             break;
           case 5: // ','
           default:
             isnum = false;
             break;
         }
      }
      ++i;
      while ( i < s.size() && s.at(i) == ' ') {
         ++i;
         whitespace = true;
      }
    }
    if (!isnum) {
      stringstream ss;
      for (int j = 0; j < i-1; ++j)
        ss << " ";
      ss << "^";
      cerr << "checkFloatString() Error at position " << i-1 << " in: \n" << s << "\n" << ss.str() << "\n";
    }
    return isnum;
};


bool checkComplexString(const std::string &s) {
    bool isnum = true;
    bool period = false;
    bool e = false;
    bool comma = false;
    bool whitespace = false;
    bool leftbrace = false;
    bool rightbrace = false;
    int numi;
    int modi;
    int m = -1;
    int i = 0;
    //skip initial whitespace
    while (s.at(i) == ' ')
      ++i;
    while (isnum && i < s.size()) {
      numi = NUMS.find(s.at(i));
      modi = CARRAY.find(s.at(i));
      if (numi == string::npos && modi == string::npos)
         isnum = false;
      else if ( numi != string::npos  ){
         ++m;
         if (whitespace || !leftbrace) {
            isnum = false;
        }
        whitespace = false;
      }
      else if ( modi != string::npos ){
         switch (modi) {
           case 0: // '.'
             if (period || e || whitespace || !leftbrace) {
               isnum = false;
               //cout << "period error\n";
             }
             period = true;
             break;
           case 1: // 'E'
           case 2: // 'e'
             if (m == 0 || e || whitespace || !leftbrace) {
               isnum = false;
               //cout << "E error\n";
             }
             e = true;
             m = 0;
             break;
           case 3: // '-'
           case 4: // '+'
             if ((m != 0 ) 
                 || whitespace || !leftbrace ) {
               //cout << "+/- error (" << m << ") \n";
               isnum = false;
             }
             break;
           case 5: // ','
             if (m == 0 || i == s.size()-1 || comma) {
               isnum = false;
               //cout << ", error\n";
             }
             period = false;
             e = false;
             comma = true;
             whitespace = false;
             m = 0;
             break;
           case 6: // '('
             if (( m > 0 && !comma ) || leftbrace) 
               isnum = false;
             leftbrace = true;
             m = 0;
             break;
           case 7: // ')'
             if (!leftbrace || !comma || m == 0)
                isnum = false;
             leftbrace = false;
             break;
           default:
             isnum = false;
             break;
         }
         whitespace = false;
      }
      ++i;
      while ( isnum && i < s.size() && s.at(i) == ' ') {
         ++i;
         if (modi != 6 && modi != 7 && modi != 5)
           whitespace = true;
      }
    }
    if (!isnum) {
      stringstream ss;
      for (int j = 0; j < i-1; ++j)
        ss << " ";
      ss << "^";
      cerr << "checkComplexString() Error at position " << i-1 << " in: \n" << s << "\n" << ss.str() << "\n";
    }
    return isnum;
}


bool checkEntryString(const std::string &s) {
    bool isnum = true;
    bool period = false;
    bool e = false;
    bool comma = false;
    bool whitespace = false;
    int numi;
    int modi;
    int m = -1;
    int i = 0;
    //skip initial whitespace
    while (s.at(i) == ' ')
      ++i;
    while (isnum && i < s.size()) {
      ++m;
      numi = NUMS.find(s.at(i));
      modi = ARRAY.find(s.at(i));
      if (numi == string::npos && modi == string::npos)
         isnum = false;
      else if ( numi != string::npos  ){
         comma = false;
         if (whitespace) 
            isnum = false;
      }
      else if ( modi != string::npos ){
         if (modi != 5)
           comma = false;
         switch (modi) {
           case 0: // '.'
             if (period || e || whitespace) 
               isnum = false;
             period = true;
             break;
           case 1: // 'E'
           case 2: // 'e'
             if (m == 0 || e || whitespace)
               isnum = false;
             e = true;
             break;
           case 3: // '-'
           case 4: // '+'
             if ((m != 0 && s.at(i-1) != 'e' && s.at(i-1) != 'E' ) || whitespace )
               isnum = false;
             break;
           case 5: // ','
             if (m == 0 || i == s.size()-1 || comma)
               isnum = false;
             period = false;
             e = false;
             comma = true;
             whitespace = false;
             m = -1;
             break;
           default:
             isnum = false;
             break;
         }
      }
      ++i;
      while ( i < s.size() && s.at(i) == ' ') {
         ++i;
         if (modi != 5)
           whitespace = true;
      }
    }
    if (!isnum) {
      stringstream ss;
      for (int j = 0; j < i-1; ++j)
        ss << " ";
      ss << "^";
      cerr << "checkEntryString() Error at position " << i-1 << " in: \n" << s << "\n" << ss.str() << "\n";
    }
    return isnum;
};


bool checkComplexEntryString(const std::string &s) {
    bool isnum = true;
    bool period = false;
    bool e = false;
    bool comma = false;
    bool whitespace = false;
    bool leftbrace = false;
    bool rightbrace = false;
    int numi;
    int modi;
    int m = -1;
    int i = 0;
    //skip initial whitespace
    while (s.at(i) == ' ')
      ++i;
    while (isnum && i < s.size()) {
      numi = NUMS.find(s.at(i));
      modi = CARRAY.find(s.at(i));
      if (numi == string::npos && modi == string::npos)
         isnum = false;
      else if ( numi != string::npos  ){
         ++m;
         if (whitespace || !leftbrace) {
            isnum = false;
            //cout << "num error\n";
        }
         comma = false;
         whitespace = false;
      }
      else if ( modi != string::npos ){
         switch (modi) {
           case 0: // '.'
             if (period || e || whitespace || !leftbrace) {
               isnum = false;
               //cout << "period error\n";
             }
             period = true;
             break;
           case 1: // 'E'
           case 2: // 'e'
             if (m == 0 || e || whitespace || !leftbrace) {
               isnum = false;
               //cout << "E error\n";
             }
             e = true;
             m = 0;
             break;
           case 3: // '-'
           case 4: // '+'
             if ((m != 0 ) 
                 || whitespace || !leftbrace ) {
               //cout << "+/- error (" << m << ") \n";
               isnum = false;
             }
             break;
           case 5: // ','
             if (m == 0 || i == s.size()-1 || comma) {
               isnum = false;
               //cout << ", error\n";
             }
             period = false;
             e = false;
             comma = true;
             whitespace = false;
             m = 0;
             break;
           case 6: // '('
             if (( m > 0 && !comma ) || leftbrace) 
               isnum = false;
             leftbrace = true;
             m = 0;
             break;
           case 7: // ')'
             if (!leftbrace || m == 0)
                isnum = false;
             leftbrace = false;
             break;
           default:
             isnum = false;
             break;
         }
         if (modi != 5)
           comma = false;
         whitespace = false;
      }
      ++i;
      while ( isnum && i < s.size() && s.at(i) == ' ') {
         ++i;
         if (modi != 6 && modi != 7 && modi != 5)
           whitespace = true;
      }
    }
    if (!isnum) {
      stringstream ss;
      for (int j = 0; j < i-1; ++j)
        ss << " ";
      ss << "^";
      cerr << "checkComplexEntryString() Error at position " << i-1 << " in: \n" << s << "\n" << ss.str() << "\n";
    }
    return isnum;
};
