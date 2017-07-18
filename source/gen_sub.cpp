#include <iostream>
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
	fs::path fastq_dir(BASE_DIR.c_str());
	if (exists(fastq_dir) && fs::is_directory(fastq_dir))
	{
		for (fs::directory_entry ite : fs::directory_iterator(fastq_dir))
			if (fs::is_regular_file(ite))
			{
				files.push_back(ite.path().string());
				filenames.push_back(ite.path().filename().string());
			}
		//sort(files.begin(), files.end());
	}

	string sub_dir = "subset/";
	string gen_dir = BASE_DIR + sub_dir;
	fs::create_directory(gen_dir);

	char buff[256];
	int buff_len = 256;
	int subset_len = 10000;
	for (int i = 2; i < files.size(); ++i)
	{
		cout << files[i] << "..." << endl;
		gzFile gfp = gzopen(files[i].c_str(), "r");
		gzFile gf_out = gzopen((gen_dir + filenames[i]).c_str(), "w");
		int len;
		int line_cnt = 0;
		int qt = 0;
		while (len = gzread(gfp, buff, buff_len))
		{
			for (int j = 0; (!qt) && j < len; ++j)
			{
				if (buff[j] == '\n')
					++line_cnt;
				if (line_cnt == subset_len * 4)
				{
					qt = 1;
					len = j+1;
				}
			}
			buff[len] = 0;
			gzwrite(gf_out, buff, len);
			if (qt)
				break;
		}
		gzclose(gfp);
		gzclose(gf_out);
	}
}
