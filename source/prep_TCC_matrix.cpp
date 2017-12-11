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
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>  
#include <zlib.h>
#include "kseq.h"

using namespace std;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

double TCCdistance(int p, int q, int lastp, int lastq, int* gene, double* TCCdata, double countSum)
{
	double distance = countSum;
	while (1)
	{
		while (p < lastp && gene[p] < gene[q])
			++p;
		if (p >= lastp)
			return distance;
		if (gene[p] == gene[q])
			distance -= 2*min(TCCdata[p++], TCCdata[q++]);
		while (q < lastq && gene[q] < gene[p])
			++q;
		if (q >= lastq)
			return distance;
	}
	return distance;
}

int main(int argc, char* argv[])
{

	if(argc != 2)
	{
		cout << "arg cnt" << endl;
		return -1;
	}
	if (!fs::is_regular_file(argv[1]))
	{
		cout << "json" << endl;
		return -1;
	}
	pt::ptree root;
	pt::read_json(argv[1], root);
	
	int D_MIN = root.get<int>("dmin");
	int BARCODE_LENGTH = root.get<int>("BARCODE_LENGTH");
	int NUM_THREADS = root.get<int>("NUM_THREADS");

	string SOURCE_DIR = root.get<string>("SOURCE_DIR");
	string SAVE_DIR = root.get<string>("SAVE_DIR");
	string BASE_DIR = root.get<string>("BASE_DIR");
	string OUTPUT_DIR = root.get<string>("OUTPUT_DIR");
	string TCC_output = root.get<string>("kallisto.TCC_output");
	string kallisto_binary = root.get<string>("kallisto.binary");
	string kallisto_index = root.get<string>("kallisto.index");

	vector<int> WINDOWS;
	pt::ptree child_win = root.get_child("WINDOW");
	BOOST_FOREACH(pt::ptree::value_type &v, child_win)
		WINDOWS.push_back(v.second.get_value<int>());

	string ecfile_dir = TCC_output + "matrix.ec";
	string tsvfile_dir = TCC_output + "matrix.tsv";

	cout << "Loading TCCs..." << endl;

	FILE *fp;
	fp = fopen(tsvfile_dir.c_str(), "rb");
	int a, b, x;
	vector<int> rows, cols, data;
	while (fscanf(fp, "%d%d%d", &a, &b, &x) != EOF)
	{
		cols.push_back(a);
		rows.push_back(b);
		data.push_back(x);
	}
	fclose(fp);

	unordered_map<int, int> map_rows;
	vector<int> uni_rows(rows), uni_cols(cols);
	sort(uni_rows.begin(), uni_rows.end());
	uni_rows.erase(unique(uni_rows.begin(), uni_rows.end()), uni_rows.end());
	for (int i = 0; i < uni_rows.size(); ++i)
		map_rows[uni_rows[i]] = i;

	int* gene = new int[cols.size()];
	for (int i=0;i<cols.size();++i)
		gene[i] = cols[i];

	int NUM_OF_CELLS = uni_rows.size();
	cout << "NUM_OF_CELLS = " << NUM_OF_CELLS << endl;
	double* rows_sum = new double[NUM_OF_CELLS];
	memset(rows_sum, 0, uni_rows.size() * sizeof(double));
	for (int i = 0; i < rows.size(); ++i)
		rows_sum[map_rows[rows[i]]] += data[i];

	double* TCCsum = new double[NUM_OF_CELLS];
	memset(TCCsum, 0, NUM_OF_CELLS * sizeof(double));
	double* TCCdata = new double[rows.size()];
	for (int i = 0; i < rows.size(); ++i)
	{
		TCCdata[i] = data[i] / rows_sum[map_rows[rows[i]]];
		TCCsum[map_rows[rows[i]]] += TCCdata[i];
	}
	vector<int> TCCidx;
	TCCidx.push_back(0);
	for (int i = 1; i < rows.size(); ++i)
		if (rows[i] != rows[i-1])
			TCCidx.push_back(i);
	TCCidx.push_back(rows.size());
	cout << "NUM idx " << TCCidx.size() << endl;

	cout << "Calculating pairwise L1 distances..." << endl;
	double** dist = new double*[NUM_OF_CELLS];
	for (int i = 0; i < NUM_OF_CELLS; ++i)
	{
		dist[i] = new double[NUM_OF_CELLS];
		memset(dist[i], 0, NUM_OF_CELLS * sizeof(double));
	}

	int cnt = 0;
	int t0 = clock();
	time_t tt0 = time(NULL);

	for (int i = 0; i < NUM_OF_CELLS; ++i)
	{	
		//cout << endl << i << endl;
	#pragma omp parallel for num_threads(NUM_THREADS)
		for (int j = i + 1; j < NUM_OF_CELLS; ++j)
		{
			//cout << j << ' ';
			dist[i][j] = TCCdistance(TCCidx[i], TCCidx[j], TCCidx[i+1], TCCidx[j+1], gene, TCCdata, TCCsum[i] + TCCsum[j]);
			dist[j][i] = dist[i][j];
		}
	}
	
	cout << "time: " << (time(NULL) - tt0) << " s" << endl;

	string pwise_dist_file = SAVE_DIR + "pwise_dist_L1.txt";
	fp = fopen(pwise_dist_file.c_str(), "w");
	for (int i=0;i<NUM_OF_CELLS;++i)
	{
		for (int j=0;j<NUM_OF_CELLS;++j)
			fprintf(fp, "%lf ", dist[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 0;
}
