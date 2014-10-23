
/**
 * Author: Adilson Pavan
 * RA:     159018
 *
 * Description: Numerical implementation of the acoustic analogy of Ffowcs Williams & Hawkings.
 *
 * Unicamp MO644/MC900 (Parallel Programming) 2014s1.
 *
 * Compile: g++ -pedantic -Wall -lpthread -lm -fopenmp -DUse_Complex -o surf *.cpp
 *
 * Run: Use the parameters below
 *      -t          = number of threads (default=1)
 *      -o          = number of observers (default=360)
 *      -serial     = activates the serial version (default=false)
 *      -parallel   = activates the parallel version (default=true)
 *      -output     = activates output on standard discipline MO644/MC900 (default=true)
 *
 *      example: ./surf -t 4 -o 360 -serial -parallel -output input01
 */


#include <omp.h>
#include <math.h>
#include <iomanip> 


#include "Tools.h" 
#include "TSurfIntegral.h"

int main(int argc, char const **argv)
{

    // receives parameters
    Tools tools;
    tools.SetParameters(argc,argv);


    // get data input files of pressure
    std::vector<double> VecRealPress, VecImPress;
    tools.GetDataRealPressure(VecRealPress);
    tools.GetDataPressureImaginary(VecImPress);

    if(VecRealPress.empty() && VecImPress.empty()) return 0;


    // get data pressure complex
    std::vector<double> VecPress;
    std::vector<std::complex<double> > VecPressComplex;
    tools.GetPressComplex(VecPress,VecPressComplex,VecRealPress,VecImPress);
    
    if (VecPress.empty() && VecPressComplex.empty()) return 0;
    

    // get NACA0012 airfoil points
    std::vector<double> VecX, VecY;
    tools.GetDataNACA(VecX,VecY);

    if (VecX.empty() && VecY.empty()) return 0;
    

    /* execute simulation */

    std::cout << "\nStarting simulation\n" << std::endl;
    std::cout << "Threads: " << tools.GetThreadCounts() << std::endl;
    std::cout << "Observers: " << tools.GetObservers() << std::endl;

    std::map<double, double> PatTeta;
    double elapsed_serial, elapsed_parallel;

    TSurfIntegral SurfIntegral;
    SurfIntegral.SetValues(VecX,VecY,VecPressComplex);

    if (tools.ExecSerial())
    {
        elapsed_serial = SurfIntegral.PressCalcForAllObservers(PatTeta,tools.GetObservers());
        std::cout << "Elapsed serial: " << elapsed_serial << std::endl;
    }

    if (tools.ExecParallel())
    {
        elapsed_parallel = SurfIntegral.PressCalcForAllObserversParallel(PatTeta,tools.GetObservers(),tools.GetThreadCounts());
        std::cout << "Elapsed parallel: " << elapsed_parallel << std::endl;
    }

    double speedup, efficiency;
    if (tools.ExecSerial() && tools.ExecParallel())
    {
        speedup = elapsed_serial / elapsed_parallel;
        efficiency = speedup / tools.GetThreadCounts();

        std::cout << "Speedup: " << speedup << std::endl; 
        std::cout << "Efficiency: " << efficiency << std::endl; 
    }

    std::cout << std::endl;

    if (tools.PrintOutput())
    {
        std::fstream fileOutput;
        fileOutput.open("output//output.dat",std::fstream::out | std::fstream::app);
        fileOutput << std::endl;
        fileOutput << tools.Filename() << ".dat --------------" << std::endl;
        fileOutput << std::endl;
        fileOutput << "Speedup: " << speedup << std::endl; 
        fileOutput.close();

        tools.SaveResults(tools.Filename(),PatTeta,VecX,VecY,VecPress);
    }
    
    /* end of the simulation */


    return 0;
}
