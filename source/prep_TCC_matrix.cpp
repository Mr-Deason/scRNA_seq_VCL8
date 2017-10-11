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
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>  
#include <zlib.h>
#include "kseq.h"

using namespace std;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

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

	int NUM_OF_CELLS = uni_rows.size();
	cout << "NUM_OF_CELLS = " << NUM_OF_CELLS << endl;
	double* rows_sum = new double[uni_rows.size()];
	memset(rows_sum, 0, uni_rows.size() * sizeof(double));
	for (int i = 0; i < rows.size(); ++i)
		rows_sum[map_rows[rows[i]]] += data[i];

	vector < pair<int, pair<int, double> > > TCC;
	for (int i = 0; i < rows.size(); ++i)
	{
		TCC.push_back(make_pair(rows[i], make_pair(cols[i], data[i] / rows_sum[map_rows[rows[i]]])));
	}
	sort(TCC.begin(), TCC.end());
	vector<int> TCCidx;
	TCCidx.push_back(0);
	for (int i = 1; i < TCC.size(); ++i)
		if (TCC[i].first != TCC[i - 1].first)
			TCCidx.push_back(i);
	TCCidx.push_back(TCC.size());
	cout << "NUM idx " << TCCidx.size() << endl;

	cout << "Calculating pairwise L1 distances..." << endl;
	double** dist = new double*[NUM_OF_CELLS];
	for (int i = 0; i < NUM_OF_CELLS; ++i)
	{
		dist[i] = new double[NUM_OF_CELLS];
		memset(dist[i], 0, NUM_OF_CELLS * sizeof(double));
	}

	cout << uni_cols.size() << endl;
	int cnt = 0;
	int t0 = clock();
	time_t tt0 = time(NULL);
//#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < NUM_OF_CELLS; ++i)
	{	
		//cout << endl << i << endl;
		for (int j = i + 1; j < NUM_OF_CELLS; ++j)
		{
			//cout << j << ' ';
			for (int p = TCCidx[i], q = TCCidx[j]; p < TCCidx[i + 1] || q < TCCidx[j + 1];)
			{
				//cout << p << ' ' << q << endl;
				if (q >= TCCidx[j + 1] || (p < TCCidx[i + 1] && TCC[p].second.first < TCC[q].second.first))
				{
					dist[i][j] += TCC[p].second.second;
					++p;
				}
				else if (p >= TCCidx[i + 1] || (q < TCCidx[j + 1] && TCC[p].second.first > TCC[q].second.first))
				{
					dist[i][j] += TCC[q].second.second;
					++q;
				}
				else
				{
					dist[i][j] += fabs(TCC[p].second.second - TCC[q].second.second);
					++p;
					++q;
				}
				++cnt;
			}
			dist[j][i] = dist[i][j];
		}
	}
	
	cout << "time: " << (time(NULL) - tt0) << " s" << endl;
	cout << "DONE" << endl;

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
