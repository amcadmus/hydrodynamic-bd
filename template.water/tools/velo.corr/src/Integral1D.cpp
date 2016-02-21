#include "Integral1D.h"

Integral1DInfo::Method Integral1DInfo::Simpson	= 0;
Integral1DInfo::Method Integral1DInfo::Gauss3	= 1;
Integral1DInfo::Method Integral1DInfo::Gauss4	= 2;

static double global_hx;

template <typename UnitaryFunction, typename ValueType >
void IntegralBox1DGauss3<UnitaryFunction, ValueType>::
reinit (const UnitaryFunction & f_, 
	const ValueType & xlo_, const ValueType & xup_,
	const ValueType & tol_)
{
  hx = (xup_ - xlo_);
  xlo = xlo_;
  xup = xup_;
  tol = tol_;
  f = &f_;
  
  ValueType sqrt3o5 = sqrt (0.6);
  ValueType x0 = 0.5 * (1 + sqrt3o5) * xlo + 0.5 * (1 - sqrt3o5) * xup;
  ValueType x1 = 0.5 * (xlo + xup);
  ValueType x2 = 0.5 * (1 - sqrt3o5) * xlo + 0.5 * (1 + sqrt3o5) * xup;
  coarse_int = 1./18. * hx * (5 * f_(x0) + 8 * f_(x1) + 5 * f_(x2));
}

template <typename UnitaryFunction, typename ValueType >
void IntegralBox1DGauss3<UnitaryFunction, ValueType>
::cal_int (ValueType & value)
{
  IntegralBox1DGauss3 box0;
  IntegralBox1DGauss3 box1;
  ValueType hhx = 0.5 * (xup - xlo);

  box0.reinit (*f, xlo, xlo+hhx, tol*0.5);
  box1.reinit (*f, xlo+hhx, xup, tol*0.5);
  
  value = box0.coarseInt () + box1.coarseInt();
  
  if (fabs(coarse_int - value) < tol * 63) {
    return;
  }
  else if (hhx*2 < 1e-13 * global_hx){
    // // std::cout << "lala"<< std::endl;
    return;
  }
  else {
    ValueType tmp;
    value = 0;
    box0.cal_int (tmp);
    value += tmp;
    box1.cal_int (tmp);
    value += tmp;
    return;
  }
}



template <typename UnitaryFunction, typename ValueType >
void IntegralBox1DGauss4<UnitaryFunction, ValueType>::
reinit (const UnitaryFunction & f_, 
	const ValueType & xlo_, const ValueType & xup_,
	const ValueType & tol_)
{
  hx = (xup_ - xlo_);
  xlo = xlo_;
  xup = xup_;
  tol = tol_;
  f = &f_;

  ValueType value0 = sqrt ((3 - 2 * sqrt(1.2)) / 7.);
  ValueType value1 = sqrt ((3 + 2 * sqrt(1.2)) / 7.);
  
  ValueType x0 = 0.5 * (xup - xlo) * (-value1) + 0.5 * (xup + xlo);
  ValueType x1 = 0.5 * (xup - xlo) * (-value0) + 0.5 * (xup + xlo);
  ValueType x2 = 0.5 * (xup - xlo) * (+value0) + 0.5 * (xup + xlo);
  ValueType x3 = 0.5 * (xup - xlo) * (+value1) + 0.5 * (xup + xlo);

  double weigh0 = (18 + sqrt(30.)) / 36.;
  double weigh1 = (18 - sqrt(30.)) / 36.;

  coarse_int = 0.5 * hx * (weigh1 * (f_(x0) + f_(x3)) + 
			   weigh0 * (f_(x1) + f_(x2)) );
  
}

template <typename UnitaryFunction, typename ValueType >
void IntegralBox1DGauss4<UnitaryFunction, ValueType>
::cal_int (ValueType & value)
{
  IntegralBox1DGauss4 box0;
  IntegralBox1DGauss4 box1;
  ValueType hhx = 0.5 * (xup - xlo);

  box0.reinit (*f, xlo, xlo+hhx, tol*0.5);
  box1.reinit (*f, xlo+hhx, xup, tol*0.5);
  
  value = box0.coarseInt () + box1.coarseInt();
  
  if (fabs(coarse_int - value) < tol * 255) {
    return;
  }
  else if (hhx*2 < 1e-13 * global_hx){
    // // std::cout << "lala"<< std::endl;
    return;
  }
  else {
    ValueType tmp;
    value = 0;
    box0.cal_int (tmp);
    value += tmp;
    box1.cal_int (tmp);
    value += tmp;
    return;
  }
}


template <typename BinaryFunction, typename ValueType >
ValueType Simpson_cal_int (const BinaryFunction & f,
			   const ValueType & c, const ValueType & d, 
			   const ValueType & tol, 
			   const int & ny)
{
  ValueType integration = 0;
    
  std::list<ValueType> r (ny+1);
  std::list<ValueType> fr (ny+1);
  std::list<int> tag (ny+1, 0);
  ValueType ddy = (d - c)/ny;
    
  typename std::list<ValueType>::iterator i;
  int count = 0;
  for (i = r.begin(); i != r.end(); i++){
    *i = (count ++) * ddy + c;
  }

  typename std::list<ValueType>::iterator this_r = r.begin();
  typename std::list<ValueType>::iterator next_r = ++ r.begin();
  typename std::list<ValueType>::iterator this_fr = fr.begin();
  typename std::list<ValueType>::iterator next_fr = ++ fr.begin();
  typename std::list<int>::iterator this_tag = tag.begin();
  typename std::list<int>::iterator next_tag = ++ tag.begin();
    
  int goon = 0;
  while (next_r != r.end()){
    ValueType dx = * next_r - * this_r;
    if (fabs(dx) < 1e-13){
      goon = 1;
    }
    ValueType my_tol = tol * dx / fabs(d-c);
    ValueType mid_r = 0.5 * (* next_r + * this_r);
    ValueType mid_fr;
	
    if (* this_tag == 0){
      * this_fr = f(* this_r);
      * this_tag = 1;
    }
    if (* next_tag == 0){
      * next_fr = f(* next_r);
      * next_tag = 1;
    }
    //////////////////////////////
    // here maybe optimazed !!
    //////////////////////////////
    mid_fr = f(mid_r);
	
    ValueType s1 = 1./6. * (*this_fr + 4*mid_fr + *next_fr) * dx;
    ValueType s21 = 1./6. * (*this_fr + 4 * f(0.5*(mid_r+*this_r)) + mid_fr) * 0.5 * dx;
    ValueType s22 = 1./6. * (*next_fr + 4 * f(0.5*(mid_r+*next_r)) + mid_fr) * 0.5 * dx;

	
    if (fabs (s1 - s21 - s22) < 31 * my_tol || goon == 1){
      integration += s21 + s22;
      this_r ++;
      next_r ++;
      this_fr ++;
      next_fr ++;
      this_tag ++;
      next_tag ++;
      if (goon == 1){
	goon = 0;
      }
    }
    else {
      next_r = r.insert (next_r, mid_r);
      next_fr = fr.insert (next_fr, mid_fr);
      next_tag = tag.insert (next_tag, 1);
    }
  }
    
  return integration;
}



template <typename UnitaryFunction, typename ValueType >
ValueType Integral1D<UnitaryFunction, ValueType >::
cal_int (const Integral1DInfo::Method & method,
	 const UnitaryFunction & f,
	 const ValueType & xlo, const ValueType & xup,
	 const ValueType & tol,
	 const int & init_nx)
{
  
  hx = (xup - xlo) / ValueType (init_nx);
  global_hx = hx;
  ValueType nxi = ValueType(1)/ValueType(init_nx);
  ValueType value = 0;

  if (method == Integral1DInfo::Simpson){
    value = Simpson_cal_int (f, xlo, xup, tol, init_nx);
  }
  else if (method == Integral1DInfo::Gauss3){
    for (int i = 0; i < init_nx; i ++){
      ValueType x = xlo + i * hx;
      box_Gauss3.reinit (f, x, x+hx, tol*nxi);
      ValueType tmp;
      box_Gauss3.cal_int (tmp);
      value += tmp;
    }
  }
  else if (method == Integral1DInfo::Gauss4){
    for (int i = 0; i < init_nx; i ++){
      ValueType x = xlo + i * hx;
      box_Gauss4.reinit (f, x, x+hx, tol*nxi);
      ValueType tmp;
      box_Gauss4.cal_int (tmp);
      value += tmp;
    }
  }
  
  return value;
}

// typedef double (*MFP) (double );      

// double global_x;
// double global_u;
// unsigned global_m;


// double intPow (const double & x, const unsigned & alpha_)
// {
//   double sum = 1;
//   unsigned alpha = alpha_;
//   for (; alpha > 0; --alpha){
//     sum *= x;
//   }
//   return sum;
// }

// double value (double t) {
//   return (cos(2*M_PI*global_x/t) - 1) * intPow (t, 2*global_m) / t / t
//       / intPow (global_u*t + 2*M_PI, 2*global_m);
// }


// int main(int argc, char * argv[])
// {
//   global_x = 32;
//   global_u = 0.5*M_PI;
//   global_m = 3;
  
//   Integral1D <MFP, double > inte;
//   std::cout.precision (16);
//   std::cout << inte.cal_int (Integral1DInfo::Gauss4, value, 0, 1, 1e-12, 2*64) << std::endl;

//   return 0;
  
// }



#include "ModeOverlapping.h"

// typedef double (*F1) ( double);

template class Integral1D<ModeOverlapping,double >;
