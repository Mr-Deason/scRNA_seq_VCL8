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
#include <assert.h>
#include <omp.h>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>  
#include <zlib.h>
#include <cblas.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "kseq.h"

using namespace std;
using namespace Eigen;
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

	unordered_map<int, int> map_rows, map_cols;
	vector<int> uni_rows(rows), uni_cols(cols);
	sort(uni_rows.begin(), uni_rows.end());
	sort(uni_cols.begin(), uni_cols.end());
	uni_rows.erase(unique(uni_rows.begin(), uni_rows.end()), uni_rows.end());
	uni_cols.erase(unique(uni_cols.begin(), uni_cols.end()), uni_cols.end());
	for (int i = 0; i < uni_rows.size(); ++i)
		map_rows[uni_rows[i]] = i;
	for (int i = 0; i < uni_cols.size(); ++i)
		map_cols[uni_cols[i]] = i;

	int NUM_OF_CELLS = uni_rows.size();
	cout << "NUM_OF_CELLS = " << NUM_OF_CELLS << endl;

	Eigen::initParallel();

	Matrix<double, Dynamic, Dynamic, RowMajor> TCCmat = Matrix<double, Dynamic, Dynamic, RowMajor>::Zero(uni_rows.size(), uni_cols.size());
	SparseMatrix<double,RowMajor> TCCspMat(uni_rows.size(), uni_cols.size());
	vector<Triplet<double> > tripletList;
	for (int i=0;i<rows.size();++i)
		tripletList.push_back(Triplet<double>(map_rows[rows[i]], map_cols[cols[i]], data[i]*1.0));
	TCCspMat.setFromTriplets(tripletList.begin(), tripletList.end());
	double row_sum;
	for (int i=0;i<TCCspMat.rows();++i)
	{
		row_sum = TCCspMat.row(i).sum();
		TCCspMat.row(i) /= row_sum;
	}

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

	Eigen::setNbThreads(NUM_THREADS);
//#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < NUM_OF_CELLS; ++i)
	{
		const auto &row = TCCspMat.row(i);
		//#pragma omp parallel for shared(row)
		#pragma omp parallel for num_threads(NUM_THREADS)
		for (int j = i+1; j < NUM_OF_CELLS; ++j)
		{
			dist[i][j] = (row-TCCspMat.row(j)).lpNorm<1>();
			dist[j][i] = dist[i][j];
		}
	}

	time_t tt1 = time(NULL);

	cout << "time: " << (tt1 - tt0) << " s" << endl;
	
	cout << "time: " << (time(NULL) - tt1) << " s" << endl;
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
