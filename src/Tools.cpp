
#include "Tools.h"

/*------------------------------------------------------------------------------------------------------------------------------------------------*/
Tools::Tools()
{
	// initialize parameters

	this->NObservers = 360;
	this->thread_counts = 1;

	this->exec_serial = false;
	this->exec_parallel = true;
	this->print_output = true;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
Tools::~Tools()
{
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void Tools::SetParameters(int argc, char const **argv)
{
	if (argc > 0)
	{
		for (int i = 0; i < argc; ++i)
		{
			std::string command = argv[i];
			if (command == "-t")
			{
				this->thread_counts = atoi(argv[i+1]);
			}
			else if (command == "-o")
			{
				this->NObservers = atoi(argv[i+1]);	
			}
			else if (command == "-serial")
			{
				this->exec_serial = true;	
			}
			else if (command == "-parallel")
			{
				this->exec_parallel = true;	
			}
			else if (command == "-output")
			{
				this->print_output = true;
				this->filename = argv[i+1];
			}
		}
	}
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void Tools::GetDataRealPressure(std::vector<double> &vec)
{

#if defined(_WIN32) || defined(_WIN64)
    std::ifstream fileRealPress("input\\real_pressure.in");
#else        
    std::ifstream fileRealPress("input//real_pressure.in");
#endif

#ifdef DEBUG
    if (fileRealPress.is_open())
    {
		std::cout << "##### FILE IS OPEN #####"<< std::endl;
    }
	else 
    {
		std::cout << "##### FILE IS NOT OPEN #####"<< std::endl;
		return;
	}
#endif

	std::string line;
	fileRealPress >> line;
    int Np = atoi( line.c_str() );

    for(int i = 0; i < Np; i++)
    {
        fileRealPress >> line;
        double RealPress = atof( line.c_str() );
        vec.push_back(RealPress);
    }

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void Tools::GetDataPressureImaginary(std::vector<double> &vec)
{

#if defined(_WIN32) || defined(_WIN64)
    std::ifstream fileImPress("input\\pressure_imaginary.in");
#else        
    std::ifstream fileImPress("input//pressure_imaginary.in");
#endif

#ifdef DEBUG
    if (fileImPress.is_open())
    {
		std::cout << "##### FILE IS OPEN #####"<< std::endl;
    }
	else 
    {
		std::cout << "##### FILE IS NOT OPEN #####"<< std::endl;
		return;
	}
#endif

	std::string line;
	fileImPress >> line;
    int Np = atoi( line.c_str() );

    for(int i = 0; i < Np; i++)
    {
        fileImPress >> line;
        double ImPress = atof( line.c_str() );
        vec.push_back(ImPress);
    }

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void Tools::GetDataNACA(std::vector<double> &VecX, std::vector<double> &VecY)
{

#if defined(_WIN32) || defined(_WIN64)
    std::ifstream fileNaca0012("input\\NACA0012.dat");
#else
    std::ifstream fileNaca0012("input//NACA0012.dat");
#endif

#ifdef DEBUG
    if (fileNaca0012.is_open())
    {
		std::cout << "##### FILE IS OPEN #####"<< std::endl;
    }
	else 
    {
		std::cout << "##### FILE IS NOT OPEN ##### 333"<< std::endl;
		return;
	}
#endif

	std::string line;
	fileNaca0012 >> line;
    int Np = atoi( line.c_str() );

    for(int i = 0; i < Np; i++)
    {
        fileNaca0012 >> line;
        double X = atof(line.c_str());
        fileNaca0012 >> line;
        double Y = atof(line.c_str());

        VecX.push_back(X);
        VecY.push_back(Y);
    }

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void Tools::GetPressComplex(std::vector<double> &VecPress, std::vector<std::complex<double> > &VecPressComplex,
							const std::vector<double> &VecRealPress, const std::vector<double> &VecImPress)
{
	for(unsigned int i = 0; i < VecImPress.size(); i++)
    {
        double Press = sqrt(VecImPress[i]*VecImPress[i] + VecRealPress[i]*VecRealPress[i]);
        VecPress.push_back(Press);

        std::complex<double> ComplexValue(VecRealPress[i], VecImPress[i]);
        VecPressComplex.push_back(ComplexValue);
    }
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void Tools::SaveResults(std::string inputFilename, std::map<double, double> &PatTeta, std::vector<double> &VecX,
						std::vector<double> &VecY, std::vector<double> &VecPress)
{
	// Pressao At Teta
    
    std::string name = "output//PressaoAtTeta_" + inputFilename + ".nb";
	std::ofstream fileMath(name.c_str());

    fileMath.setf(std::ios::fixed);
    fileMath.precision(15);

    fileMath << "ListA = {";

    std::map<double, double>::iterator it;
    it = PatTeta.begin();

    fileMath << "{" << it->first << ", " << it->second <<"}";

    for(it = PatTeta.begin(); it != PatTeta.end(); it++)
    {
        if(it == PatTeta.begin()) continue;
        fileMath << ",{" << it->first << ", " << it->second <<"}";
    }

    fileMath << "};" << std::endl;
    fileMath << "ListPolarPlot[" << "ListA, AspectRatio->Automatic]";
    fileMath.close();


    // NACA version Mathematica

    std::ofstream fileMathNACA("output//NACA0012.nb");

    fileMathNACA << "ListA = {";
    fileMathNACA << "{" << VecX[0] << ", " << VecY[0] <<"}";

    for(unsigned int i = 1; i < VecX.size(); i ++)
    {
        fileMathNACA << ",{" << VecX[i] << ", " << VecY[i] <<"}";
    }
    fileMathNACA << "};" << std::endl;
    fileMathNACA << "ListPlot[" << "ListA, AspectRatio->Automatic]";
    fileMathNACA.close();


    // Pressao At Point

    name = "output//PressAtPoint_" + inputFilename + ".nb";
    std::ofstream fileMathPress(name.c_str());

    fileMathPress << "ListA = {";
    fileMathPress << "{" << 1 << ", " << VecPress[0] <<"}";

    for(unsigned int i = 1; i < VecPress.size(); i ++)
    {
        fileMathPress << ",{" << i+1 << ", " << VecPress[i] <<"}";
    }
    fileMathPress << "};" << std::endl;
    fileMathPress << "ListPlot[" << "ListA, AxesLabel -> {Ponto,p'}]";
    fileMathPress.close();

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
