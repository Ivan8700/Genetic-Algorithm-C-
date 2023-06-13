#include <iostream>
#include <fstream>
#include <time.h>
#include <chrono>
#include "RunAlgorithm.h"
using namespace std::chrono;
using std::cout; using std::cin; using std::endl;
const static int fromMicroToSec = 1000000;

int main()
{
	int read_input, print;
	cout << "Would you like to print everything or only the final solution ? \n 0 - everything , 1 - only final : ";
	cin >> print;
	cout << "from where would you like to get the input ?\n 0 - read from .txt file, 1 - manual input, 2 - randomization : ";
	cin >> read_input;
	/*output file to desktop*/
	//std::ofstream out_file; out_file.open("C:/Users/Ivan/Desktop/output.txt");  
	//std::streambuf* cout_buffer = cout.rdbuf();  
	//cout.rdbuf(out_file.rdbuf()); 
	auto start = high_resolution_clock::now();
	RunAlgorithm alg(read_input, print, 100, 2, 5, 2);
	//Default Properties of the algorithm
	//RunAlgorithm alg(0, 1, 100, 2, 10, 2); //(option_to_read_input, print or not, population size, mutation size, amount of different genes at the end, function)
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << endl << " Time taken by function : " << float(duration.count()) / fromMicroToSec << " seconds" << endl;
	/*Check if the file is closed*/
	//out_file.close(); // 
	//cout.rdbuf(cout_buffer); //
}
