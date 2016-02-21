#ifndef __Integral1D_h__wanghan__
#define __Integral1D_h__wanghan__

#include <iostream>
#include <vector>
#include <list>
#include <math.h>

#define INIT_N 1

struct Integral1DInfo 
{
public:
  typedef int Method;
  static Method Simpson;
  static Method Gauss3;
  static Method Gauss4;
}
    ;

template <typename UnitaryFunction, typename ValueType >
class IntegralBox1DGauss3;
template <typename UnitaryFunction, typename ValueType >
class IntegralBox1DGauss4;


template <typename UnitaryFunction, typename ValueType >
class Integral1D
{
  ValueType hx;
  IntegralBox1DGauss3 <UnitaryFunction,ValueType >  box_Gauss3;
  IntegralBox1DGauss4 <UnitaryFunction,ValueType >  box_Gauss4;
public:
  ValueType cal_int(const Integral1DInfo::Method & method,
		    const UnitaryFunction & f,
                    const ValueType & xlo, const ValueType & xup,
                    const ValueType & tol = 1e-6,
                    const int & init_nx = INIT_N);
}
  ;

template <typename UnitaryFunction, typename ValueType >
class IntegralBox1DGauss3
{
public:
  IntegralBox1DGauss3 () {};
  void reinit (const UnitaryFunction & f, 
	       const ValueType & xlo, const ValueType & xup,
	       const ValueType & tol = 1e-6);
  ValueType & coarseInt () {return coarse_int;}
  const ValueType & coarseInt () const {return coarse_int;}
  
  void cal_int (ValueType & value);
private:
  ValueType hx;
  ValueType xlo, xup;
  ValueType tol;
  const UnitaryFunction * f;
  
  ValueType coarse_int;
};

template <typename UnitaryFunction, typename ValueType >
class IntegralBox1DGauss4
{
public:
  IntegralBox1DGauss4 () {};
  void reinit (const UnitaryFunction & f, 
	       const ValueType & xlo, const ValueType & xup,
	       const ValueType & tol = 1e-6);
  ValueType & coarseInt () {return coarse_int;}
  const ValueType & coarseInt () const {return coarse_int;}
  
  void cal_int (ValueType & value);
private:
  ValueType hx;
  ValueType xlo, xup;
  ValueType tol;
  const UnitaryFunction * f;
  
  ValueType coarse_int;
};


#endif
