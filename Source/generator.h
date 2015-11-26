#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

void generate_canonical_matrix(const int number_of_row,
		const int limit_number_generate) {

	srand(time(NULL));
	int number = rand() % 100 + 10;
	ostringstream convert;
	convert << number;

	string a = "../data/problem";
	string b = convert.str();
	string c = ".txt";
	string namefile = a + b + c;

	ofstream myfile(namefile.c_str());
	if (myfile.is_open()) {

		for (int i = 0; i < number_of_row * 2; ++i)
			myfile << rand() % limit_number_generate << " ";

		myfile << endl << endl;

		for (int i = 0; i < number_of_row; ++i)
			myfile << rand() % limit_number_generate << " ";

		myfile << endl << endl;

		for (int i = 0; i < number_of_row; ++i) {
			for (int j = 0; j < number_of_row * 2; ++j) {
				if (j < number_of_row) {
					if (j == i)
						myfile << 1 << " ";
					else
						myfile << 0 << " ";
				} else
					myfile << rand() % limit_number_generate << " ";

			}
			myfile << endl;
		}

		myfile.close();
	} else {
		cerr << "Unable to save file!" << endl;
		exit(0);
	}

}
