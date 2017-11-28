#include <iostream>
#include "fasthamming.hpp"

using namespace std;
int main()
{
	nacgt<uint32_t> *a = new nacgt<uint32_t>("ACGG", 4);
	nacgt<uint32_t> *b = new nacgt<uint32_t>("ACGC", 4);

	hamming<uint32_t> hammingEval=hamming<uint32_t>(4);
	cout << hammingEval.distance(a, b, 4) << endl;

	return 0;
}