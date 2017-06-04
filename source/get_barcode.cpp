#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <time.h>
#include <omp.h>
#include <boost\filesystem.hpp>
#include <boost\property_tree\ptree.hpp>
#include <boost\property_tree\json_parser.hpp>
#include <boost\foreach.hpp>  
#include <zlib.h>
#include "kseq.h"

using namespace std;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

KSEQ_INIT(gzFile, gzread);

unsigned int code_map(char ch)
{
	if (ch == 'A') return 0;
	if (ch == 'G') return 1;
	if (ch == 'C') return 2;
	if (ch == 'T') return 3;
	if (ch == 'N') return 0*(rand()%3);
}
unsigned int encode(string str)
{
	unsigned int code = 0;
	for (int i = 0; str[i]; ++i)
	{
		code <<= 2;
		code |= code_map(str[i]);
	}
	return code;
}

char* str_AGCT = "AGCT";
string decode(unsigned int code, int len)
{
	string ret;
	ret.resize(len);
	for (int i = 0; i < len; ++i)
	{
		ret[len - i - 1] = str_AGCT[code & 3];
		code >>= 2;
	}
	return ret;
}

int main(int argc, char *argv[])
{
	return 0;
}