#ifndef __CBRT_H__
#define __CBRT_H__
/**
 * Vectorized evaluation of cubic roots.
 */

namespace CBRT
{
  /**
   * Pade approximatin for obtaining the intial guess
   * Cost 5*  4+  1/
   */
  template <typename Number>
  constexpr inline Number
  inital_guess(const Number &t)
  {
    const Number t2 = t * t;
    // Evaluate the numerator: 5 + 35t + 14t^2
    Number numerator = Number(5) + Number(35) * t + Number(14) * t2;
    // Evaluate the denominator: 14 + 35t + 5t^2
    Number denominator = Number(14) + Number(35) * t + Number(5) * t2;
    // Return the result of the division
    return numerator / denominator;
  }

  /**
   * Halley's method for cubic root.Cubic convergence.
   * Cost 4*  2+  1/
   */
  template <typename Number>
  constexpr inline Number
  halley(const Number &a, int iterations = 2)
  {
    Number       x    = inital_guess(a);
    const Number two  = Number(2.);
    const Number twoA = two * a;
    for (int i = 0; i < iterations; ++i)
      {
        Number x3 = x * x * x;
        x         = x * ((x3 + twoA) / (two * x3 + a)); // Iteration formula
      }
    return x;
  }

  /**
   * Newtons's method for cubic root. Quadratic convergence.
   */
  template <typename Number>
  constexpr inline Number
  newton(Number a, int iterations = 8)
  {
    Number x     = inital_guess(a);
    Number three = 3.;
    for (int i = 0; i < iterations; ++i)
      x = x - (x * x * x - a) / (three * x * x); // Newton's iteration formula
    return x;
  }

} // namespace CBRT
#endif // __CBRT_H__