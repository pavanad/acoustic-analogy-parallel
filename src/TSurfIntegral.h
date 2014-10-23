
#ifndef TSURFINTEGRAL
#define TSURFINTEGRAL

#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <map>


class TGaussIntegrateRule
{

  public:
      TGaussIntegrateRule();
      virtual ~TGaussIntegrateRule();

      void GetPoints(int NRule,
                     std::vector<double> &points,
                     std::vector<double> &weigths);
};


class TSurfIntegral 
{

  public:

      TSurfIntegral();
      ~TSurfIntegral();

      double PressCalcForAllObservers(std::map<double, double> &PatTeta, int NObservers);
      double PressCalcForAllObserversParallel(std::map<double, double> &PatTeta, int NObservers, int thread_counts);

      void SetValues(const std::vector<double> &X,
                     const std::vector<double> &Y,
                     const std::vector<double> &P);

      void SetValues(const std::vector<double> &X,
                     const std::vector<double> &Y,
                     const std::vector<std::complex<double> > &P);

      TGaussIntegrateRule fGaussIntRule;

      void GetNormalVec(double &Xi,
                        double &Yi,
                        double &Xj,
                        double &Yj,
                        double &Nx,
                        double &Ny);

  #ifdef Use_Complex
      void IntegrateOneElem(std::complex<double> &Pm,
                            double &Xi,
                            double &Yi,
                            double &Xj,
                            double &Yj,
                            double &Xo,
                            double &Yo,
                            std::complex<double> &IntValue);

      void IntegrateAllElements(double &Xo,
                                double &Yo,
                                std::complex<double> &IntValue);
  #else
      void IntegrateOneElem(double &Pm,
                            double &Xi,
                            double &Yi,
                            double &Xj,
                            double &Yj,
                            double &Xo,
                            double &Yo,
                            double &IntValue);

      void IntegrateAllElements(double &Xo,
                                double &Yo,
                                double &IntValue);
  #endif


  #ifdef Use_Complex
      void GetFX(double &Xs,
                 double &Ys,
                 double &Xo,
                 double &Yo,
                 std::complex<double> &Fx);

      void GetFY(double &Xs,
                 double &Ys,
                 double &Xo,
                 double &Yo,
                 std::complex<double> &Fy);

      inline void GetFvalues(double &Xs,
                      double &Ys,
                      double &Xo,
                      double &Yo,
                      std::complex<double> &Fx,
                      std::complex<double> &Fy);

  #else

      void GetFX(double &Xs,
                 double &Ys,
                 double &Xo,
                 double &Yo,
                 double &Fx);

      void GetFY(double &Xs,
                 double &Ys,
                 double &Xo,
                 double &Yo,
                 double &Fy);

      inline void GetFvalues(double &Xs,
                      double &Ys,
                      double &Xo,
                      double &Yo,
                      double &Fx,
                      double &Fy);
  #endif

      // função de Bessel tipo 1
      void BesselFuncJ(int &n, double &x, double &result );

      //função de Bessel tipo 2,ordem 0
      void BesselFuncY(int &n, double &x, double &result );

      //função de Hankel do segundo tipo, ordem 0
      void HankelFuncT2(int &n, double &x, std::complex<double> &result );

      int Fatorial(const int &k);

      double HFunc(const int &m);

      std::vector<double> fX;
      std::vector<double> fY;
      std::vector<double> fP;

      std::vector<std::complex<double> > fPcomplex;

};


#endif
