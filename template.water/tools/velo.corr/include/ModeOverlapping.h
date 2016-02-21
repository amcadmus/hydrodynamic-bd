#pragma once

#include <cmath>

struct ModeOverlapping 
{
  ModeOverlapping () : kk(0), xi(0) {};
  void reinit (const double & kk_,
	       const double & xi_ ) {kk = kk_; xi = xi_;}
  double operator () (const double & xx) const {return sin(kk * xx) * sin(xi *xx);}
private:
  double kk, xi;
};
