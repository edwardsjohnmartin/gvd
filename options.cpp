/*******************************************************
 ** Generalized Voronoi Diagram Project               **
 ** Copyright (c) 2015 John Martin Edwards            **
 ** Scientific Computing and Imaging Institute        **
 ** 72 S Central Campus Drive, Room 3750              **
 ** Salt Lake City, UT 84112                          **
 **                                                   **
 ** For information about this project contact        **
 ** John Edwards at                                   **
 **    edwardsjohnmartin@gmail.com                    **
 ** or visit                                          **
 **    sci.utah.edu/~jedwards/research/gvd/index.html **
 *******************************************************/

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>

#include "./options.h"

namespace oct {

using namespace std;

// const float Constants::kWidthf = Constants::kWidth;
// const double Constants::kWidthd = Constants::kWidth;

bool OctreeOptions::ProcessArg(int& i, char** argv) {
  OctreeOptions& o = *this;
  int orig_i = i;
  if (strcmp(argv[i], "-l") == 0) {
    ++i;
    o.max_level = atoi(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "-a") == 0) {
    ++i;
    o.ambiguous_max_level = atoi(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "--no-buffer") == 0) {
    o.make_buffer = false;
    ++i;
  } else if (strcmp(argv[i], "-h") == 0) {
    o.help = true;
    ++i;
  } else if (strcmp(argv[i], "--gpu") == 0) {
    o.gpu = true;
    ++i;
  } else if (strcmp(argv[i], "--cpu") == 0) {
    o.gpu = false;
    ++i;
  } else if (strcmp(argv[i], "--opencl-log") == 0) {
    o.opencl_log = true;
    ++i;
  } else if (strcmp(argv[i], "--cell-of-interest") == 0) {
    ++i;
    o.cell_of_interest = atoi(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "--level-of-interest") == 0) {
    ++i;
    o.level_of_interest = atoi(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "--test") == 0) {
    ++i;
    o.test = atoi(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "--bb-scale") == 0) {
    ++i;
    o.bb_scale = atof(argv[i]);
    ++i;
  } else if (strcmp(argv[i], "--rs") == 0) {
    o.restricted_surface = true;
    ++i;
  } else if (strcmp(argv[i], "--center") == 0) {
    ++i;
    o.center = atoi(argv[i]);
    ++i;
  } 
  return i != orig_i;
}

// trim from start
string& ltrim(string &s) {
  s.erase(s.begin(), find_if(s.begin(),
                             s.end(), not1(ptr_fun<int, int>(isspace))));
  return s;
}

// trim from end
string& rtrim(string &s) {
  s.erase(find_if(s.rbegin(), s.rend(),
                  not1(ptr_fun<int, int>(isspace))).base(), s.end());
  return s;
}

// trim from both ends
string& trim(string &s) {
  return ltrim(rtrim(s));
}

void OctreeOptions::ReadOptionsFile() {
  ifstream in("gvd.config");
  if (!in) return;
  while (!in.eof()) {
    string key;
    in >> key;
    string value;
    getline(in, value);
    
    if (!key.empty() && key[0] != '#')
      key2value[key] = trim(value);
  }
  in.close();
}

string OctreeOptions::Value(
    const string& key, const string& default_value) const {
  if (key2value.find(key) == key2value.end())
    return default_value;
  const string value = key2value.find(key)->second;
  return value;
}

bool OctreeOptions::BoolValue(
    const string& key, const bool default_value) const {
  const string value = Value(key, default_value?"true":"false");
  if (value == "0" || value == "false" || value == "False" || value == "FALSE")
    return false;
  return true;
}

int OctreeOptions::IntValue(
    const string& key, const int default_value) const {
  stringstream ss;
  ss << default_value;
  const string value = Value(key, ss.str());
  return atoi(value.c_str());
}

}
