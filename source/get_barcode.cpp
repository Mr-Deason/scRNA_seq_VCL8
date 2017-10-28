#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cinttypes>
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
#include "fasthamming.hpp"

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
uint32_t code_map_n(char ch)
{
	if (ch == 'A') return 0;//000
	if (ch == 'G') return 1;//001
	if (ch == 'C') return 2;//010
	if (ch == 'T') return 3;//011
	if (ch == 'N') return 4;//100
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
uint64_t encode_n(string str)
{
	uint64_t code = 0;
	for (int i=0;str[i];++i)
	{
		code <<= 3;
		code |= code_map_n(str[i]);
	}
	return code;
}

char* str_AGCTN = "AGCTN";
string decode(uint64_t code, int len)
{
	string ret;
	ret.resize(len);
	for (int i = 0; i < len; ++i)
	{
		ret[len - i - 1] = str_AGCTN[code & 7];
		code >>= 3;
	}
	return ret;
}
void devode_cstr(uint64_t code, int len, const char* str)
{
	str[len] = 0;
	for (int i = 0; i < len; ++i)
	{
		str[len - i - 1] = str_AGCTN[code & 7];
		code >>= 3;
	}
}



#define  MMAX  6

typedef  double MAT[MMAX + 2][MMAX + 2];

void LUBKSB(MAT A, int N, int np, int *INDX, double *B) {

	float SUM;
	int I, II, J, LL;

	II = 0;

	for (I = 1; I <= N; I++) {
		LL = INDX[I];
		SUM = B[LL];
		B[LL] = B[I];
		if (II != 0)
			for (J = II; J<I; J++)
				SUM -= A[I][J] * B[J];
		else if (SUM != 0.0)
			II = I;
		B[I] = SUM;
	}

	for (I = N; I>0; I--) {
		SUM = B[I];
		if (I < N)
			for (J = I + 1; J <= N; J++)
				SUM -= A[I][J] * B[J];
		B[I] = SUM / A[I][I];
	}

}
void LUDCMP(MAT A, int N, int np, int *INDX, int *D, int *CODE) {

	double AMAX, DUM, SUM, TINY;
	double VV[1000];
	int   I, IMAX, J, K;

	TINY = (double)1e-12;

	*D = 1; *CODE = 0;

	for (I = 1; I <= N; I++) {
		AMAX = 0.0;
		for (J = 1; J <= N; J++)
			if (fabs(A[I][J]) > AMAX) AMAX = (double)fabs(A[I][J]);
		if (AMAX < TINY) {
			*CODE = 1;
			return;
		}
		VV[I] = (double)1.0 / AMAX;
	}

	for (J = 1; J <= N; J++) {
		for (I = 1; I<J; I++) {
			SUM = A[I][J];
			for (K = 1; K<I; K++)
				SUM -= A[I][K] * A[K][J];
			A[I][J] = SUM;
		}
		AMAX = 0.0;
		for (I = J; I <= N; I++) {
			SUM = A[I][J];
			for (K = 1; K<J; K++)
				SUM -= A[I][K] * A[K][J];
			A[I][J] = SUM;
			DUM = VV[I] * (double)fabs(SUM);
			if (DUM >= AMAX) {
				IMAX = I;
				AMAX = DUM;
			}
		}

		if (J != IMAX) {
			for (K = 1; K <= N; K++) {
				DUM = A[IMAX][K];
				A[IMAX][K] = A[J][K];
				A[J][K] = DUM;
			}
			*D = -(*D);
			VV[IMAX] = VV[J];
		}

		INDX[J] = IMAX;
		if ((double)fabs(A[J][J]) < TINY)  A[J][J] = TINY;

		if (J != N) {
			DUM = (double)1.0 / A[J][J];
			for (I = J + 1; I <= N; I++)  A[I][J] *= DUM;
		}
	} //j loop

} //LUDCMP()

int polyfit(const double* const dependentValues,
	const double* const independentValues,
	unsigned int        countOfElements,
	unsigned int        order,
	double*             coefficients)
{
	// Declarations...
	// ----------------------------------
	enum { maxOrder = 5 };

	double B[maxOrder + 1] = { 0.0f };
	double P[((maxOrder + 1) * 2) + 1] = { 0.0f };
	double A[(maxOrder + 1) * 2 * (maxOrder + 1)] = { 0.0f };

	double x, y, powx;

	unsigned int ii, jj, kk;

	// Verify initial conditions....
	// ----------------------------------

	// This method requires that the countOfElements > 
	// (order+1) 
	if (countOfElements <= order)
		return -1;

	// This method has imposed an arbitrary bound of
	// order <= maxOrder.  Increase maxOrder if necessary.
	if (order > maxOrder)
		return -1;

	// Begin Code...
	// ----------------------------------

	// Identify the column vector
	for (ii = 0; ii < countOfElements; ii++)
	{
		x = dependentValues[ii];
		y = independentValues[ii];
		powx = 1;

		for (jj = 0; jj < (order + 1); jj++)
		{
			B[jj] = B[jj] + (y * powx);
			powx = powx * x;
		}
	}

	// Initialize the PowX array
	P[0] = countOfElements;

	// Compute the sum of the Powers of X
	for (ii = 0; ii < countOfElements; ii++)
	{
		x = dependentValues[ii];
		powx = dependentValues[ii];

		for (jj = 1; jj < ((2 * (order + 1)) + 1); jj++)
		{
			P[jj] = P[jj] + powx;
			powx = powx * x;
		}
	}

	// Initialize the reduction matrix
	//
	for (ii = 0; ii < (order + 1); ii++)
	{
		for (jj = 0; jj < (order + 1); jj++)
		{
			A[(ii * (2 * (order + 1))) + jj] = P[ii + jj];
		}

		A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1;
	}

	// Move the Identity matrix portion of the redux matrix
	// to the left side (find the inverse of the left side
	// of the redux matrix
	for (ii = 0; ii < (order + 1); ii++)
	{
		x = A[(ii * (2 * (order + 1))) + ii];
		if (x != 0)
		{
			for (kk = 0; kk < (2 * (order + 1)); kk++)
			{
				A[(ii * (2 * (order + 1))) + kk] =
					A[(ii * (2 * (order + 1))) + kk] / x;
			}

			for (jj = 0; jj < (order + 1); jj++)
			{
				if ((jj - ii) != 0)
				{
					y = A[(jj * (2 * (order + 1))) + ii];
					for (kk = 0; kk < (2 * (order + 1)); kk++)
					{
						A[(jj * (2 * (order + 1))) + kk] =
							A[(jj * (2 * (order + 1))) + kk] -
							y * A[(ii * (2 * (order + 1))) + kk];
					}
				}
			}
		}
		else
		{
			// Cannot work with singular matrices
			return -1;
		}
	}

	// Calculate and Identify the coefficients
	for (ii = 0; ii < (order + 1); ii++)
	{
		for (jj = 0; jj < (order + 1); jj++)
		{
			x = 0;
			for (kk = 0; kk < (order + 1); kk++)
			{
				x = x + (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] *
					B[kk]);
			}
			coefficients[ii] = x;
		}
	}

	return 0;
}

void savgol_coeff(double *c, int np, int nl, int nr, int ld, int m) {
	/*-------------------------------------------------------------------------------------------
	USES lubksb,ludcmp given below.
	Returns in c(np), in wrap-around order (see reference) consistent with the argument respns
	in routine convlv, a set of Savitzky-Golay filter coefficients. nl is the number of leftward
	(past) data points used, while nr is the number of rightward (future) data points, making
	the total number of data points used nl +nr+1. ld is the order of the derivative desired
	(e.g., ld = 0 for smoothed function). m is the order of the smoothing polynomial, also
	equal to the highest conserved moment; usual values are m = 2 or m = 4.
	-------------------------------------------------------------------------------------------*/
	int d, icode, imj, ipj, i, j, k, kk, mm;
	int indx[MMAX + 2];
	double fac, sum;
	MAT   a;
	double b[MMAX + 2];

	if (np<nl + nr + 1 || nl<0 || nr<0 || ld>m || m>MMAX || nl + nr<m) {
		printf("\n Bad args in savgol.\n");
		return;
	}

	for (i = 1; i <= MMAX + 1; i++) {
		for (j = 1; j <= MMAX + 1; j++) a[i][j] = 0.0;
		b[i] = 0.0;
		indx[i] = 0;
	}

	for (ipj = 0; ipj <= 2 * m; ipj++) { //Set up the normal equations of the desired leastsquares fit.
		sum = 0.0;
		if (ipj == 0) sum = 1.0;
		for (k = 1; k <= nr; k++) sum += (double)pow(k, ipj);
		for (k = 1; k <= nl; k++) sum += (double)pow(-k, ipj);
		mm = ipj <=  2 * m - ipj ? ipj : 2 * m - ipj;
		imj = -mm;
		do {
			a[1 + (ipj + imj) / 2][1 + (ipj - imj) / 2] = sum;
			imj += 2;
		} while (imj <= mm);
	}

	LUDCMP(a, m + 1, MMAX + 1, indx, &d, &icode);    //Solve them: LU decomposition

	for (j = 1; j <= m + 1; j++) b[j] = 0.0;
	b[ld + 1] = 1.0;    //Right-hand side vector is unit vector, depending on which derivative we want.

	LUBKSB(a, m + 1, MMAX + 1, indx, b);      //Backsubstitute, giving one row of the inverse matrix.

	for (kk = 1; kk <= np; kk++)          //Zero the output array (it may be bigger than the number
		c[kk] = 0.0;                      //of coefficients.

	for (k = -nl; k <= nr; k++) {         //Each Savitzky-Golay coefficient is the dot product
		sum = b[1];                       //of powers of an integer with the inverse matrix row.
		fac = 1.0;
		for (mm = 1; mm <= m; mm++) {
			fac *= k;
			sum += b[mm + 1] * fac;
		}
		kk = ((np - k) % np) + 1;           //Store in wrap-around order}
		c[kk] = sum;
	}
}

void edge_fit(double *data, int window_start, int window_end, int interp_start, int interp_end, int polyorder, double *y)
{

	double *poly_coeff = (double*)malloc((polyorder + 1) * sizeof(double));
	memset(poly_coeff, 0, (polyorder + 1) * sizeof(double));

	int win_len = window_end - window_start;
	double *x = (double*)malloc(win_len * sizeof(double));
	memset(x, 0, win_len * sizeof(double));
	for (int i = 0; i < win_len; ++i)
		x[i] = i;
	polyfit(x, y + window_start, win_len, polyorder, poly_coeff);

	int interp_len = interp_end - interp_start;
	for (int i = 0; i < interp_len; ++i)
	{
		double tt = 1.0;
		data[i + interp_start] = 0.0;
		for (int j = 0; j <= polyorder; ++j)
		{
			data[i + interp_start] += tt*poly_coeff[j];
			tt *= x[i];
		}
	}

	free(poly_coeff);
	free(x);
}

double* savgol_filter(double *data, int n, int win_len, int polyorder)
{
	double* ret = (double*)malloc((n+2) * sizeof(double));
	memset(ret, 0, (n + 2) *sizeof(double));
	int nl = win_len / 2;
	int nr = win_len / 2;
	int m = polyorder;

	double *c = (double*)malloc((win_len + 2) * sizeof(double));
	memset(c, 0, (win_len + 2) * sizeof(double));
	savgol_coeff(c, nl + nr + 1, nl, nr, 0, m);

	//default mode - constant
	int *index = (int*)malloc((win_len + 2) * sizeof(int));
	memset(index, 0, (win_len + 2) * sizeof(int));
	index[1] = 0;
	int j = 3;
	for (int i = 2; i <= nl + 1; i++) {
		index[i] = i - j;
		j += 2;
	}
	j = 2;
	for (int i = nl + 2; i <= nl + nr + 1; i++) {
		index[i] = i - j;
		j += 2;
	}
	for (int i = 1; i <= n - nr; i++) {
		ret[i] = 0.0;
		for (j = 1; j <= nl + nr + 1; j++)
			if (i + index[j]>0)   //skip left points that do not exist.
				ret[i] += c[j] * data[i + index[j]-1];
	}
	free(index);
	
	//interp mode - using poly fit twice
	int half = win_len / 2;
	edge_fit(ret+1, 0, win_len, 0, half, m, data);
	edge_fit(ret+1, n - win_len, n, n - half, n, m, data);

	free(c);
	return ret+1;
}

bool brc_cmp(pair<int, unsigned int> a, pair<int, unsigned int> b)
{
	return a.first > b.first;
}
bool brc_cmp1(pair<int, unsigned int> a, pair<int, unsigned int> b)
{
	return a.first >= b.first;
}
bool brc_cmp2(pair<int, unsigned int> a, pair<int, unsigned int> b)
{
	return a.first == b.first ? (a.second > b.second) : (a.first > b.first);
}
bool brc_cmp3(pair<int, unsigned int> a, pair<int, unsigned int> b)
{
	return a.first == b.first ? (a.second >= b.second) : (a.first > b.first);
}
bool brc_cmp4(pair<int, unsigned int> a, pair<int, unsigned int> b)
{
	return a.first == b.first ? (a.second < b.second) : (a.first > b.first);
}
bool brc_cmp5(pair<int, unsigned int> a, pair<int, unsigned int> b)
{
	return a.first == b.first ? (a.second <= b.second) : (a.first > b.first);
}

int main(int argc, char *argv[])
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

	vector<string> files;
	srand(clock());
	
	fs::path fastq_dir(BASE_DIR.c_str());
	if (exists(fastq_dir) && fs::is_directory(fastq_dir))
	{
		for (fs::directory_entry ite : fs::directory_iterator(fastq_dir))
			if (fs::is_regular_file(ite))
				files.push_back(ite.path().string());
		sort(files.begin(), files.end());
	}

	vector<uint64_t> barcodes;
	vector<nacgt<uint64_t> > barcodes_nacgt;
	map<uint64_t, int> barcodes_nacgt_cnt;
	//all "read-I1_*" files
//#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < files.size() / 2; ++i)
	{
		cout << files[i] << endl;
		gzFile gfp = gzopen(files[i].c_str(), "r");
		kseq_t *seq1;
		seq1 = kseq_init(gfp);
		int l;
		while ((l = kseq_read(seq1) >= 0))
		{
			string bar = seq1->seq.s;
			nacgt<uint64_t> barcode(seq1->seq.s, BARCODE_LENGTH);
			barcodes_nacgt.push_back(barcode);
			barcodes_nacgt_cnt[encode_n(bar)]++;
			//printf("%s\n", bar.c_str());
			barcodes.push_back(encode_n(bar));
		}
		gzclose(gfp);
		kseq_destroy(seq1);
	}

	FILE* fp;
	/*
	map<unsigned int, int> barcodes_cnt;
	for (int i = 0; i < barcodes.size(); ++i)
	{
		barcodes_cnt[barcodes[i]]++;
	}
	*/

	vector<pair<int, uint64_t> > cnt_bar;
	for (map<uint64_t, int>::iterator ite = barcodes_nacgt_cnt.begin(); ite != barcodes_nacgt_cnt.end(); ++ite)
	{
		cnt_bar.push_back(make_pair(ite->second, ite->first));
	}

//---
/*
	vector<pair<int, unsigned int> > cnt_bar1;
	cnt_bar1.assign(cnt_bar.begin(), cnt_bar.end());
	sort(cnt_bar1.begin(), cnt_bar1.end(), brc_cmp1);
	fp = fopen("barcodes_cnt_cpp1.txt", "w");
	for (int i=0;i<cnt_bar.size();++i)
		fprintf(fp, "%d %s\n", cnt_bar1[i].first, decode(cnt_bar1[i].second, BARCODE_LENGTH).c_str());
	fclose(fp);
	
	vector<pair<int, unsigned int> > cnt_bar2;
	cnt_bar2.assign(cnt_bar.begin(), cnt_bar.end());
	sort(cnt_bar2.begin(), cnt_bar2.end(), brc_cmp2);
	fp = fopen("barcodes_cnt_cpp2.txt", "w");
	for (int i=0;i<cnt_bar.size();++i)
		fprintf(fp, "%d %s\n", cnt_bar2[i].first, decode(cnt_bar2[i].second, BARCODE_LENGTH).c_str());
	fclose(fp);
	vector<pair<int, unsigned int> > cnt_bar3;
	vector<pair<int, unsigned int> > cnt_bar4;
	vector<pair<int, unsigned int> > cnt_bar5;
	cnt_bar3.assign(cnt_bar.begin(), cnt_bar.end());
	cnt_bar4.assign(cnt_bar.begin(), cnt_bar.end());
	cnt_bar5.assign(cnt_bar.begin(), cnt_bar.end());
	sort(cnt_bar3.begin(), cnt_bar3.end(), brc_cmp3);
	sort(cnt_bar4.begin(), cnt_bar4.end(), brc_cmp4);
	sort(cnt_bar5.begin(), cnt_bar5.end(), brc_cmp5);
	fp = fopen("barcodes_cnt_cpp3.txt", "w");
	for (int i=0;i<cnt_bar.size();++i)
		fprintf(fp, "%d %s\n", cnt_bar3[i].first, decode(cnt_bar3[i].second, BARCODE_LENGTH).c_str());
	fclose(fp);
	fp = fopen("barcodes_cnt_cpp4.txt", "w");
	for (int i=0;i<cnt_bar.size();++i)
		fprintf(fp, "%d %s\n", cnt_bar4[i].first, decode(cnt_bar4[i].second, BARCODE_LENGTH).c_str());
	fclose(fp);
	fp = fopen("barcodes_cnt_cpp5.txt", "w");
	for (int i=0;i<cnt_bar.size();++i)
		fprintf(fp, "%d %s\n", cnt_bar5[i].first, decode(cnt_bar5[i].second, BARCODE_LENGTH).c_str());
	fclose(fp);
*/
//---


	sort(cnt_bar.begin(), cnt_bar.end(), brc_cmp5);
	for (int i = 0; i < 10 && i < cnt_bar.size(); ++i)
		cout << cnt_bar[i].first << ' ' << decode(cnt_bar[i].second, BARCODE_LENGTH) << endl;
	fp = fopen("barcodes_cnt_cpp5.txt", "w");
	for (int i=0;i<cnt_bar.size();++i)
		fprintf(fp, "%d %s\n", cnt_bar[i].first, decode(cnt_bar[i].second, BARCODE_LENGTH).c_str());
	fclose(fp);


	double *diff = (double*)malloc(cnt_bar.size() * sizeof(double));
	for (int i = 0; i < cnt_bar.size() - 1; ++i)
	{
		diff[i] = cnt_bar[i + 1].first - cnt_bar[i].first;
	}

	double *yhat = savgol_filter(diff, cnt_bar.size() - 1, 151, 1);
	
	
	int border = WINDOWS[0];
	for (int i = WINDOWS[0]; i < WINDOWS[1]; ++i)
		if (yhat[border] > yhat[i])
			border = i;	

	int num_barcodes = border;
	int num_reads = 0;
	uint64_t *codewords = new uint64_t[num_barcodes];
	for (int i = 0; i < num_barcodes; ++i)
	{
		num_reads += cnt_bar[i].first;
		codewords[i] = cnt_bar[i].second;
	}


	int **dist = new int*[num_barcodes];
	for (int i = 0; i < num_barcodes; ++i)
		dist[i] = new int[num_barcodes];
	
	char *buffa = (char*)malloc((BARCODE_LENGTH+1)*sizeof(char));
	char *buffb = (char*)malloc((BARCODE_LENGTH+1)*sizeof(char));

	hamming<uint64_t> hammingEval=hamming<uint64_t>(BARCODE_LENGTH);

	for (int i = 0; i < num_barcodes; ++i)
	{
		dist[i][i] = (BARCODE_LENGTH + 1)*16;
		for (int j = i + 1; j < num_barcodes; ++j)
		{
			devode_cstr(codewords[i], BARCODE_LENGTH, buffa);
			devode_cstr(codewords[j], BARCODE_LENGTH, buffb);
			dist[i][j] = hammingEval.slowDistance(buffa, buffb, BARCODE_LENGTH, 1);
			dist[j][i] = dist[i][j];
		}
	}

	vector<int> brc_to_correct;
	for (int i = 0; i < num_barcodes; ++i)
	{
		int dmin = (BARCODE_LENGTH + 1)*16;
		for (int j = 0; j < num_barcodes; ++j)
			dmin = min(dmin, dist[i][j]);
		if (dmin >= D_MIN*16)
			brc_to_correct.push_back(i);
	}
	cout << "save" << endl;
	fs::path save_dir(SAVE_DIR.c_str());
	fs::create_directory(save_dir);
	string bar_file = "barcodes.dat";
	fp = fopen((SAVE_DIR + bar_file).c_str(), "wb");
	for (int i = 0; i < barcodes.size(); ++i)
		fprintf(fp, "%"PRIu64"\n", barcodes[i]);
	fclose(fp);
	string codes_file = "codewords.dat";
	fp = fopen((SAVE_DIR + codes_file).c_str(), "wb");
	for (int i = 0; i < num_barcodes; ++i)
		fprintf(fp, "%"PRIu64"\n", codewords[i]);
	fclose(fp);
	string brc_crt_file = "brc_idx_to_correct.dat";
	fp = fopen((SAVE_DIR + brc_crt_file).c_str(), "wb");
	for (int i = 0; i < brc_to_correct.size(); ++i)
		fprintf(fp, "%d\n", brc_to_correct[i]);
	fclose(fp);


	cout << "Cell_barcodes_detected: " << num_barcodes << endl;
	cout << "NUM_OF_READS_in_CELL_BARCODES = " << num_reads << endl;
	cout << endl << "number of cell barcodes to error-correct: " << brc_to_correct.size() << endl;
	for (int i = 0; i < num_barcodes; i++)
		delete[] dist[i];
	delete[] dist;
	
	return 0;
}

