
#include <omp.h>
#include "TSurfIntegral.h"

/*------------------------------------------------------------------------------------------------------------------------------------------------*/
TGaussIntegrateRule::TGaussIntegrateRule()
{
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
TGaussIntegrateRule::~TGaussIntegrateRule()
{
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void TGaussIntegrateRule::GetPoints(int NRule,
                                    std::vector<double> &points,
                                    std::vector<double> &weigths)
{

    points.clear();
    weigths.clear();

    switch (NRule) 
    {

        case 2:
            points.resize(2);
            weigths.resize(2);

            points[0] = -0.5773502692;
            points[1] = 0.5773502692;

            weigths[0] = 1.;
            weigths[1] = 1.;

            break;

        case 3:
            points.resize(3);
            weigths.resize(3);

            points[0] = -0.774596669241483;
            points[1] = 0.0;
            points[2] = 0.774596669241483;

            weigths[0] = 0.555555555555555;
            weigths[1] = 0.8888888888888888;
            weigths[2] = 0.555555555555555;

            break;

        case 4:
            points.resize(4);
            weigths.resize(4);

            points[0] = -0.861136311594053;
            points[1] = -0.339981043584856;
            points[2] = 0.339981043584856;
            points[3] = 0.861136311594053;

            weigths[0] = 0.3478548451;
            weigths[1] = 0.6521451549;
            weigths[2] = 0.6521451549;
            weigths[3] = 0.3478548451;

            break;

        case 5:
            points.resize(5);
            weigths.resize(5);

            points[0] = -0.9061798459;
            points[1] = -0.5384693101;
            points[2] = 0.0;
            points[3] = 0.5384693101;
            points[4] = 0.9061798459;

            weigths[0] = 0.2369268850;
            weigths[1] = 0.4786286705;
            weigths[2] = 0.568888888888888;
            weigths[3] = 0.4786286705;
            weigths[4] = 0.2369268850;

            break;

        default:
            break;
    }

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
TSurfIntegral::TSurfIntegral()
{

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
TSurfIntegral::~TSurfIntegral()
{

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
double TSurfIntegral::PressCalcForAllObservers(std::map<double, double> &PatTeta, int NObservers)
{
    double start, end, elapsed;
    start = omp_get_wtime();

    int i;
    PatTeta.clear();

    const double mPI = 4.*atan(1.);
    double Radius = 5.;

    //inicia com y = 0 e x = 5, segunido no sentido anti-horario

    double DTeta = 2.*mPI / (NObservers);

    for(i = 0; i < NObservers; i++)
    {

        double Teta = i * DTeta;

        double Xo = Radius * cos(Teta);
        double Yo = Radius * sin(Teta);

#ifdef Use_Complex
    std::complex<double> Press;
#else
    double Press;
#endif 

        this->IntegrateAllElements(Xo, Yo, Press);

        double AbsPress;
#ifdef Use_Complex
        AbsPress = std::abs(Press);
#else
        AbsPress = Press;
#endif

        PatTeta.insert(std::make_pair(Teta, AbsPress));
    }

    end = omp_get_wtime();
    elapsed = end - start;

    return elapsed;    
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
double TSurfIntegral::PressCalcForAllObserversParallel(std::map<double, double> &PatTeta, int NObservers, int thread_counts)
{
    double start, end, elapsed;
    start = omp_get_wtime();

    int i;
    PatTeta.clear();

    const double mPI = 4.*atan(1.);
    double Radius = 5.;

    //inicia com y = 0 e x = 5, segunido no sentido anti-horario

    double DTeta = 2.*mPI / (NObservers);

    #pragma omp parallel for num_threads(thread_counts) default(none) shared(NObservers,Radius,PatTeta,DTeta) private(i)
    for(i = 0; i < NObservers; i++)
    {

        double Teta = i * DTeta;

        double Xo = Radius * cos(Teta);
        double Yo = Radius * sin(Teta);

#ifdef Use_Complex
    std::complex<double> Press;
#else
    double Press;
#endif 

        this->IntegrateAllElements(Xo, Yo, Press);

        double AbsPress;
#ifdef Use_Complex
        AbsPress = std::abs(Press);
#else
        AbsPress = Press;
#endif

        PatTeta.insert(std::make_pair(Teta, AbsPress));
    }

    end = omp_get_wtime();
    elapsed = end - start;

    return elapsed;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void TSurfIntegral::SetValues(const std::vector<double> &X,
                              const std::vector<double> &Y,
                              const std::vector<double> &P)
{
    this->fX = X;
    this->fY = Y;
    this->fP = P;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void TSurfIntegral::SetValues(const std::vector<double> &X,
                   const std::vector<double> &Y,
                   const std::vector<std::complex<double> > &P)
{
    this->fX = X;
    this->fY = Y;
    this->fPcomplex = P;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
#ifdef Use_Complex
void TSurfIntegral::IntegrateAllElements(double &Xo,
                                         double &Yo,
                                         std::complex<double> &IntValue)
{

    std::complex<double> NullValue(0.,0.);
    IntValue = NullValue;

#else
void TSurfIntegral::IntegrateAllElements(double &Xo,
                                         double &Yo,
                                         double &IntValue)
{
    IntValue = 0.;
#endif // Use_Complex

    int i;
    int NPoints = this->fX.size();

    for(i = 0; i < NPoints; i++) 
    {
        double Xi, Yi, Xj, Yj;

        Xi = this->fX[i];
        Yi = this->fY[i];

        Xj = this->fX[i+1];
        Yj = this->fY[i+1];

#ifdef Use_Complex
        std::complex<double> Pi, Pj, Pm;
        Pi = this->fPcomplex[i];
        Pj = this->fPcomplex[i+1];
        std::complex<double> IntValueOneElem(0.,0.);
#else
        double Pi, Pj, Pm;
        Pi = this->fP[i];
        Pj = this->fP[i+1];
        double IntValueOneElem = 0.;
#endif

        Pm =(Pi+Pj)/2.;

        this->IntegrateOneElem(Pm, Xi, Yi, Xj, Yj, Xo, Yo, IntValueOneElem);
        
        IntValue += IntValueOneElem;
    }

    IntValue *= -1.;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void TSurfIntegral::GetNormalVec(double &Xi,
                      double &Yi,
                      double &Xj,
                      double &Yj,
                      double &Nx,
                      double &Ny)
{

        double Rx = Xi - Xj;
        double Ry = Yi - Yj;

        Ny = (-1.) * 1. / sqrt(1. + Ry*Ry/(Rx*Rx));
        Nx = (-1.) * sqrt(1. - Ny*Ny);
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
#ifdef Use_Complex
void TSurfIntegral::IntegrateOneElem(std::complex<double> &Pm,
                          double &Xi,
                          double &Yi,
                          double &Xj,
                          double &Yj,
                          double &Xo,
                          double &Yo,
                          std::complex<double> &IntValue)
{

    std::complex<double> NullValue(0.,0.);
    IntValue = NullValue;

#else
void TSurfIntegral::IntegrateOneElem(double &Pm,
                                     double &Xi,
                                     double &Yi,
                                     double &Xj,
                                     double &Yj,
                                     double &Xo,
                                     double &Yo,
                                     double &IntValue)
{

    IntValue = 0.;
#endif
    int NRule = 5, i; //de 2 a 5

    std::vector<double> points;
    std::vector<double> weigths;

    this->fGaussIntRule.GetPoints(NRule, points, weigths);

    double Nx, Ny;
    this->GetNormalVec(Xi, Yi, Xj, Yj, Nx, Ny);

    for(i = 0; i < NRule; i++)
    {
        double Z = points[i];
        double w = weigths[i];
        double X = 0.5*(Xj-Xi)*Z + 0.5*(Xj+Xi);
        double Y = 0.5*(Yj-Yi)*Z + 0.5*(Yj+Yi);

#ifdef Use_Complex
        std::complex<double> Fx, Fy;
#else
        double Fx, Fy;
#endif // Use_Complex

        this->GetFvalues(X, Y, Xo, Yo, Fx, Fy);

        IntValue += Pm*(Nx*Fx+Ny*Fy)*w;
    }

    double dx = Xj-Xi;
    double dy = Yj-Yi;

    double L = sqrt(dx*dx + dy*dy);

    IntValue *= L/2.;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
#ifdef Use_Complex
void TSurfIntegral::GetFvalues(double &Xs,
                    double &Ys,
                    double &Xo,
                    double &Yo,
                    std::complex<double> &Fx,
                    std::complex<double> &Fy)
{
#else
void TSurfIntegral::GetFvalues(double &Xs,
                    double &Ys,
                    double &Xo,
                    double &Yo,
                    double &Fx,
                    double &Fy)
{

#endif

    this->GetFX(Xs, Ys, Xo, Yo, Fx);
    this->GetFY(Xs, Ys, Xo, Yo, Fy);
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
#ifdef Use_Complex
void TSurfIntegral::GetFX(double &Xs,
               double &Ys,
               double &Xo,
               double &Yo,
               std::complex<double> &Fx)
{
#else
void TSurfIntegral::GetFX(double &Xs,
               double &Ys,
               double &Xo,
               double &Yo,
               double &Fx)
{

#endif

    double Mach = 0.5;
    double k = 5.5;

    double ArgHankel = k * sqrt( (Xo-Xs)*(Xo-Xs) + (1-Mach*Mach)*(Yo-Ys)*(Yo-Ys) ) / (1-Mach*Mach);
    std::complex<double> Hankel20;
    int n0 = 0;
    this->HankelFuncT2(n0, ArgHankel, Hankel20);
    std::complex<double> Hankel21;
    int n1 = 1;
    this->HankelFuncT2(n1, ArgHankel, Hankel21);

    std::complex<double> i(0,1);

    std::complex<double> FirstPart, SecondPart;

    FirstPart = - 0.25*exp(i*k*Mach*(Xo-Xs)/(1 - Mach*Mach))*k*Mach*Hankel20 / pow( (1 - Mach*Mach), 3./2.);

    SecondPart = - 0.25*i*exp(i*k*Mach*(Xo-Xs)/(1 - Mach*Mach))*k*(Xo-Xs)*Hankel21;
    SecondPart = SecondPart / pow( (1-Mach*Mach), 3./2. );
    SecondPart = SecondPart / sqrt( (Xo-Xs)*(Xo-Xs) + (1 - Mach*Mach)*(Yo-Ys)*(Yo-Ys) );

#ifdef Use_Complex
    Fx = FirstPart + SecondPart;
#else
    Fx = std::abs(FirstPart) + std::abs(SecondPart);
#endif // Use_Complex

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
#ifdef Use_Complex
void TSurfIntegral::GetFY(double &Xs,
               double &Ys,
               double &Xo,
               double &Yo,
               std::complex<double> &Fy)
{
#else
void TSurfIntegral::GetFY(double &Xs,
               double &Ys,
               double &Xo,
               double &Yo,
               double &Fy)
{

#endif

    double k = 5.5;
    double Mach = 0.5;

    double ArgHankel = k * sqrt( (Xo-Xs)*(Xo-Xs) + (1-Mach*Mach)*(Yo-Ys)*(Yo-Ys) ) / (1-Mach*Mach);

    std::complex<double> Hankel21;
    int n1 =1;
    this->HankelFuncT2(n1, ArgHankel, Hankel21);

    std::complex<double> func;

    std::complex<double> i(0,1);

    func = - 0.25*i*exp(i*k*Mach*(Xo-Xs)/(1 - Mach*Mach))*k*(Yo-Ys)*Hankel21;
    func = func / pow( (1-Mach*Mach), 1./2. );
    func = func / sqrt( (Xo-Xs)*(Xo-Xs) + (1 - Mach*Mach)*(Yo-Ys)*(Yo-Ys) );

#ifdef Use_Complex
    Fy = func;
#else
    Fy = std::abs(func);
#endif // Use_Complex

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void TSurfIntegral::BesselFuncJ(int &n, double &x, double &result )
{

    /*const int NMax = 10;
    double func = 0.;*/

   /* for(int k = 0; k < NMax; k++){

        double value = pow(-1,k) * pow(x, 2*k) / ( pow(2., 2*k+n) * this->Fatorial(k) * this->Fatorial(k+n) );
        func += value;

    }*/
/*
    for(int s = 0; s < NMax; s++){

        double Prod = 1.;

        for(int i = 1; i < s; i++){
            Prod *= (0.5*x)*(0.5*x) / (i*(n+i));
        }

        func += pow(-1, s) * Prod;

    }


    result = pow(0.5*x, n) * func / this->Fatorial(n);*/

#if defined(_WIN32) || defined(_WIN64)
    result = _jn(n, x);
#else
    result = jn(n, x);
#endif

}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
int TSurfIntegral::Fatorial(const int &k)
{
    int result = 1;
    int number = k;
    for(int i = 0; i < k; i++)
    {
        result *= number;
        number--;
    }

    return result;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void TSurfIntegral::BesselFuncY(int &n, double &x, double &result )
{
/*
    const double mPI = 4.*atan(1.);
    const double gama = 0.577215664901532;

    double Jn;
    this->BesselFuncJ(n,x,Jn);

    double LnTerm = log(x/2.) + gama;

    const int mMax = 10;
    //primeiro somatorio
    double sum1 = 0.;
    for(int m = 0; m < mMax; m++){
        double Hm = this->HFunc(m) + this->HFunc(m+n);
        sum1 += (pow(-1,m-1) * Hm) / (pow(2, 2*m+n) * this->Fatorial(m) * this->Fatorial(m+n) ) * pow(x, 2*m);

    }
    sum1 *= pow(x,n);

    //segundo somatorio
    double sum2 = 0.;
    for(int m = 0; m < n-1; m++){

        sum2 += this->Fatorial(n-m-1)/(pow(2,2*m-n)*this->Fatorial(m)) * pow(x, 2*m);
    }
    sum2 *= pow(x,-n) / n;

    result = (2./mPI) * (Jn*LnTerm ) + sum1 - sum2;*/

#if defined(_WIN32) || defined(_WIN64)
    result = _yn(n,x);
#else
    result = yn(n,x);
#endif
    
 }
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
double TSurfIntegral::HFunc(const int &m)
{
    if(m == 0) return 0.0;

    double result = 1.;
    int number = m;

    for(int i = 1; i < m; i++)
    {
        result += 1. / number;
        number--;
    }

    return result;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/
void TSurfIntegral::HankelFuncT2(int &n, double &x, std::complex<double> &result )
{
    double Jn;
    this->BesselFuncJ(n, x, Jn);

    double Yn;
    this->BesselFuncY(n, x, Yn);

    std::complex<double> HankelValue(Jn, -Yn);

    result = HankelValue;
}
/*------------------------------------------------------------------------------------------------------------------------------------------------*/