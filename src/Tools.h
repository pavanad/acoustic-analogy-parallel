

#ifndef TOOLS
#define TOOLS

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <complex>
#include <map>


class Tools
{

	public:

		Tools();
		virtual ~Tools();

		void SetParameters(int argc, char const **argv);

		void GetDataRealPressure(std::vector<double> &vec);

		void GetDataPressureImaginary(std::vector<double> &vec);

		void GetDataNACA(std::vector<double> &VecX, std::vector<double> &VecY);		

		void GetPressComplex(std::vector<double> &VecPress, std::vector<std::complex<double> > &VecPressComplex,
							 const std::vector<double> &VecRealPress, const std::vector<double> &VecImPress);

		int GetObservers()
		{
			return NObservers;
		}

		int  GetThreadCounts()
		{
			return thread_counts;
		}

		bool ExecSerial()
		{
			return exec_serial;
		}

		bool ExecParallel()
		{
			return exec_parallel;
		}

		bool PrintOutput()
		{
			return print_output;
		}

		std::string Filename()
		{
			return filename;
		}

		void SaveResults(std::string inputFilename, std::map<double, double> &PatTeta, std::vector<double> &VecX, std::vector<double> &VecY, std::vector<double> &VecPress);


	private:

		int NObservers;
		int thread_counts;

		bool exec_serial;
		bool exec_parallel;

		bool print_output;
		std::string filename;


};


#endif