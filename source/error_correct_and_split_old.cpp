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
void decode_cstr(uint64_t code, int len, char* str)
{
	str[len] = 0;
	for (int i = 0; i < len; ++i)
	{
		str[len - i - 1] = str_AGCTN[code & 7];
		code >>= 3;
	}
}

vector<uint64_t> hamming_circle(uint64_t barcode, int len, int d)
{
	vector<uint64_t> cousins;
	if (d == 1)
	{
		for (int i = 0; i < len; ++i)
		{
			uint64_t b = (barcode >> (i*3)) & 7;
			for (int k = 1; k < 5; ++k)
			{
				uint64_t cousin = (barcode ^ (b << (i*3))) | ((b+k)%5 << (i*3));
				cousins.push_back(cousin);
			}
		}
	}
	else if (d == 2)
	{
		for (int i = 0; i < len; ++i)
		{
			uint64_t b1 = (barcode >> (i*3)) & 7;
			for (int j = i+1; j < len; ++j)
			{
				uint64_t b2 = (barcode >> (j*3)) & 7;
				for (int k = 1; k < 5; ++k)
					for (int l = 1; l < 5; ++l)
					{
						uint64_t cousin = (((barcode ^ (b1 << (i*3))) | ((b1+k)%5 << (i*3))) ^ (b2 << (j*3))) | ((b2+l)%5 << (j*3)) ;
						cousins.push_back(cousin);
					}
			}
		}
	}
	else
		cout << "invalid parameter!" << endl;
	return cousins;
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

	int t0 = clock();

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

	vector<uint64_t> barcodes;
	vector<uint64_t> codewords;
	vector<int> brc_idx_to_correct;

	uint64_t read;
	int idx;

	FILE *fp;
	string bar_file = "barcodes.dat";
	fp = fopen((SAVE_DIR + bar_file).c_str(), "rb");
	while (fscanf(fp, "%" SCNu64, &read) != EOF)
		barcodes.push_back(read);
	fclose(fp);
	string codes_file = "codewords.dat";
	fp = fopen((SAVE_DIR + codes_file).c_str(), "rb");
	while (fscanf(fp, "%" SCNu64, &read) != EOF)
		codewords.push_back(read);
	fclose(fp);
	string brc_crt_file = "brc_idx_to_correct.dat";
	fp = fopen((SAVE_DIR + brc_crt_file).c_str(), "rb");
	while (fscanf(fp, "%d", &idx) != EOF)
		brc_idx_to_correct.push_back(idx);
	fclose(fp);

	map<uint64_t, int> barcode_to_idx;
	for (int i = 0; i < codewords.size(); ++i)
		barcode_to_idx[codewords[i]] = i;

	set<uint64_t> brc_to_correct;
	for (int i = 0; i < brc_idx_to_correct.size(); ++i)
	{
		brc_to_correct.insert(codewords[brc_idx_to_correct[i]]);
	}

	vector<vector<int> > ret;
	ret.resize(codewords.size());
//#pragma omp parallel for num_threads(NUM_THREADS)
	for (int i = 0; i < barcodes.size(); ++i)
	{
		uint64_t barcode = barcodes[i];
		if (barcode_to_idx.find(barcode) != barcode_to_idx.end())
		{
			ret[barcode_to_idx[barcode]].push_back(i);
		}
		else
		{
			vector<uint64_t> cousins = hamming_circle(barcode, BARCODE_LENGTH, 1);
			for (int j=0;j<cousins.size();++j)
				if (brc_to_correct.find(cousins[j]) != brc_to_correct.end())
				{
					ret[barcode_to_idx[cousins[j]]].push_back(i);
					break;
				}
		}
	}

	int t1 = clock();
	int NUM_OF_READS_in_CELL_BARCODES = 0;
	for (int i = 0; i < codewords.size(); ++i)
		NUM_OF_READS_in_CELL_BARCODES += ret[i].size();
	cout << "NUM_OF_READS_in_CELL_BARCODES (after error-correct): " << NUM_OF_READS_in_CELL_BARCODES << endl;


	vector<string> files;
	fs::path fastq_dir(BASE_DIR.c_str());
	if (exists(fastq_dir) && fs::is_directory(fastq_dir))
	{
		for (fs::directory_entry ite : fs::directory_iterator(fastq_dir))
			if (fs::is_regular_file(ite))
				files.push_back(ite.path().string());
		sort(files.begin(), files.end());
	}

	vector<uint64_t> line_byte_idx;
	line_byte_idx.push_back(0);
	uint64_t byte_cnt = 0;

	//all "read-RA_*" files
	cout << "merge all reads..." << endl;
	string all_reads_file = "all_reads.fastq";
	fp = fopen((SAVE_DIR+all_reads_file).c_str(), "wb");
	char buff[1024];
	int buff_len = 1024;
	for (int i = files.size() / 2; i < files.size(); ++i)
	{
		cout << files[i] << endl;
		gzFile gfp = gzopen(files[i].c_str(), "r");
		while (gzgets(gfp, buff, buff_len))
		{
			int len = strlen(buff);
			fputs(buff, fp);
			line_byte_idx.push_back(byte_cnt + len);
			byte_cnt += len;
		}
		gzclose(gfp);
	}
	fclose(fp);

	int t2 = clock();
	/*
	string line_byte_idx_ = "line_byte_idx";
	fp = fopen((SAVE_DIR+line_byte_idx_).c_str(), "wb");
	for (int i=0;i<line_byte_idx.size();++i)
		fprintf(fp, "%ld\n", line_byte_idx[i]);
	fclose(fp);

	vector<uint64_t> line_byte_idx_o;
	line_byte_idx_o.push_back(0);
	byte_cnt = 0;

	//all "read-RA_*" files
	cout << "merge all reads..." << endl;
	fp = fopen((SAVE_DIR+all_reads_file).c_str(), "wb");
	for (int i = files.size() / 2; i < files.size(); ++i)
	{
		gzFile gfp = gzopen(files[i].c_str(), "r");
		int len;
		while (len = gzread(gfp, buff, buff_len))
		{
			fwrite(buff, sizeof(char), len, fp);
			for (int j = 0; j < len; ++j)
			{
				if (buff[j] == '\n')
					line_byte_idx_o.push_back(byte_cnt + j + 1);
			}
			byte_cnt += len + 1;
			if (len < buff_len)
				break;
		}
		gzclose(gfp);
	}
	fclose(fp);

	string line_byte_idx_ori = "line_byte_idx_ori";
	fp = fopen((SAVE_DIR+line_byte_idx_ori).c_str(), "wb");
	for (int i=0;i<line_byte_idx_o.size();++i)
		fprintf(fp, "%ld\n", line_byte_idx_o[i]);

	line_byte_idx = line_byte_idx_o;

	*/

	fp = fopen((SAVE_DIR+all_reads_file).c_str(), "r");
	string umi_read_file = "umi_read_list.txt";
	fs::path output_dir(OUTPUT_DIR.c_str());
	fs::create_directory(output_dir);
	FILE *fp_umi_list = fopen((OUTPUT_DIR+umi_read_file).c_str(), "wb");
	int flag = 0;
	char filename[40], fastq_file[100], fastq_gz_file[100], umi_file[100];
	for (int i = 0; i < codewords.size(); ++i)
	{
		sprintf(filename, "cell_%04d_%s", i, decode(codewords[i], BARCODE_LENGTH).c_str());
		//sprintf(fastq_file, "%s%s.fastq", OUTPUT_DIR.c_str(), filename);
		sprintf(fastq_gz_file, "%s%s.fastq.gz", OUTPUT_DIR.c_str(), filename);
		sprintf(umi_file, "%s%s.umi", OUTPUT_DIR.c_str(), filename);
		//FILE *ffq = fopen(fastq_file, "wb");
		FILE *fum = fopen(umi_file, "wb");
		gzFile gffq = gzopen(fastq_gz_file, "w");
		//cout << "writing " << filename << endl;

		fprintf(fp_umi_list, "%s\t%s\t%s\n", filename, umi_file, fastq_gz_file);
		for (int j = 0; j < ret[i].size(); ++j)
		{
			uint64_t r = ret[i][j];
			uint64_t line = r * 8;
			fseek(fp, line_byte_idx[line] , SEEK_SET);
			for (;line < r*8 +6;++line)
			{ 
				//int len = line_byte_idx[line + 1] - line_byte_idx[line];
				//fread(buff, sizeof(char), len, fp);
				fgets(buff, buff_len, fp);
				if (line < r * 8 + 4)
				{
					//fwrite(buff, sizeof(char), len, ffq);
					//fputs(buff, ffq);
					//gzwrite(gffq, buff, len);
					gzputs(gffq, buff);
				}
				if (line > r * 8 + 4)
					//fwrite(buff, sizeof(char), len, fum);
					fputs(buff, fum);
			}
		}
		//fclose(ffq);
		gzclose(gffq);
		fclose(fum);
	}
	fclose(fp);
	fclose(fp_umi_list);

	int t3 = clock();

	cout << "calc ret " <<  (t1 - t0) << endl;
	cout << "merge " << (t2 - t1) << endl;
	cout << "split " << (t3 - t2) << endl;
	return 0;
}

