#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <zlib.h>
#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>  

using namespace std;
namespace fs = boost::filesystem;
namespace pt = boost::property_tree;

int main(int argc, char *argv[])
{
	if (argc != 2)
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
	vector<string> filenames;
	fs::path fastq_dir((BASE_DIR + "fastq_files/").c_str());
	if (exists(fastq_dir) && fs::is_directory(fastq_dir))
	{
		for (fs::directory_entry ite : fs::directory_iterator(fastq_dir))
			if (fs::is_regular_file(ite) && ite.path().filename().string()[0] == 'r')
			{
				files.push_back(ite.path().string());
				filenames.push_back(ite.path().filename().string());
			}
		//sort(files.begin(), files.end());
	}


	int subset_len[] = {10000, 100000, 1000000, 5000000, 10000000};
	string sub_dir[] = {"subset_10k/", "subset_100k/", "subset_1m/", "subset_5m/", "subset_10m/"};
	for (int case = 0; case < 5; ++case)
	{

		FILE *fs = fopen("ids", "r");
		vector<int> idx;
		int id;
		for (int i=0;i<subset_len[case];++i)
		{
			fscanf(fs, "%d", &id);
			idx.push_back(id);
		}
		sort(idx.begin(), idx.end());
		fclose(fs);


		string gen_dir = BASE_DIR + sub_dir[case];
		fs::create_directory(gen_dir);

		char buff[256];
		int buff_len = 256;
		for (int i = 0; i < files.size(); ++i)
		{
			cout << files[i] << "..." << endl;
			gzFile gfp = gzopen(files[i].c_str(), "r");
			gzFile gf_out = gzopen((gen_dir + filenames[i]).c_str(), "w");
			int times = 4;
			if (filenames[i][5] == 'R')
				times = 8;
			for (int j = 0, line = 0; j < idx.size(); ++line)
			{
				gzgets(gfp, buff, buff_len);
				if (line / times == idx[j])
				{
					gzputs(gf_out, buff);
					if (line % times == times-1)
						++j;
				}
			}
			gzclose(gfp);
			gzclose(gf_out);
		}
	}
}
