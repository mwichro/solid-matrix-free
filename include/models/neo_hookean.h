/*************************************************************
* AceGen    8.002 Linux (11 Oct 23)                          *
*           Co. J. Korelc  2020           28 Nov 23 19:01:29 *
**************************************************************
User     : Limited evaluation version
Notebook : local_vector
Evaluation time                 : 0 s     Mode  : Optimal
Number of formulae              : 2       Method: Automatic
Subroutine                      : Residual size: 64
Total size of Mathematica  code : 64 subexpressions
Total size of C code            : 365 bytes */

/******************* S U B R O U T I N E *********************/


#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

using namespace dealii;

template <typename Number>
inline constexpr Number
Power(Number x, unsigned int n)
{
  Number result = 1;
  while (n > 0)
    {
      if (n % 2 == 1)
        {
          result *= x;
        }
      x *= x;
      n /= 2;
    }
  return result;
}

template <int DIMENSION>
struct SolidModel
{
  const static unsigned int dim      = DIMENSION;
  const static unsigned int n_cached = dim == 3 ? 27 : 9;
  template <typename Number>
  static inline void
  residual([[maybe_unused]] const Tensor<1, dim, Number> &uIn,
           const Tensor<2, dim, Number> &                 graduIn,
           Tensor<1, dim, Number> &                       valueOut,
           Tensor<2, dim, Number> &                       gradientOut,
           const Number &                                 mu,
           const Number &                                 lambda);

  template <typename Number>
  static inline void
  tangent(const Tensor<2, dim, Number> &graduIn,
          const Tensor<2, dim, Number> &gradduIn,
          Tensor<2, dim, Number> &      gradientOut,
          const Number &                mu,
          const Number &                lambda);

  template <typename Number>
  static inline void
  cache(const Tensor<2, dim, Number> &graduIn,
        Tensor<1, dim, Number> &      valueOut,
        Tensor<2, dim, Number> &      gradientOut,
        ArrayView<Number> &           cacheOut,
        const Number &                mu,
        const Number &                lambda);

  template <typename Number>
  static inline void
  cache_deformed(const Tensor<2, dim, Number> &graduIn,
                 Tensor<2, dim, Number> &      gradientOut,
                 ArrayView<Number> &           cacheOut,
                 const Number &                mu,
                 const Number &                lambda);
};



// =====
// dim=2
// =====


template <>
template <typename Number>
inline void
SolidModel<2>::residual([[maybe_unused]] const Tensor<1, dim, Number> &uIn,
                        const Tensor<2, dim, Number> &                 graduIn,
                        Tensor<1, dim, Number> &                       valueOut,
                        Tensor<2, dim, Number> &gradientOut,
                        const Number &          mu,
                        const Number &          lambda)
{
  Number acegen_scratch__10__, acegen_scratch__11__, acegen_scratch__12__,
    acegen_scratch__13__, acegen_scratch__14__, acegen_scratch__15__,
    acegen_scratch__19__, acegen_scratch__20__, acegen_scratch__21__,
    acegen_scratch__22__, acegen_scratch__23__, acegen_scratch__29__,
    acegen_scratch__30__, acegen_scratch__31__, acegen_scratch__32__,
    acegen_scratch__34__, acegen_scratch__40__, acegen_scratch__41__,
    acegen_scratch__9__;
  acegen_scratch__41__ = mu / 2e0;
  acegen_scratch__9__  = 1e0 + graduIn[0][0];
  acegen_scratch__19__ = 2e0 * acegen_scratch__9__;
  acegen_scratch__10__ = graduIn[0][1];
  acegen_scratch__21__ = 2e0 * acegen_scratch__10__;
  acegen_scratch__11__ = graduIn[1][0];
  acegen_scratch__20__ = 2e0 * acegen_scratch__11__;
  acegen_scratch__12__ = 1e0 + graduIn[1][1];
  acegen_scratch__22__ = 2e0 * acegen_scratch__12__;
  acegen_scratch__13__ = (acegen_scratch__11__ * acegen_scratch__11__) +
                         (acegen_scratch__9__ * acegen_scratch__9__);
  acegen_scratch__14__ = acegen_scratch__11__ * acegen_scratch__12__ +
                         acegen_scratch__10__ * acegen_scratch__9__;
  acegen_scratch__15__ = (acegen_scratch__10__ * acegen_scratch__10__) +
                         (acegen_scratch__12__ * acegen_scratch__12__);
  acegen_scratch__23__ = sqrt(-(acegen_scratch__14__ * acegen_scratch__14__) +
                              acegen_scratch__13__ * acegen_scratch__15__);
  acegen_scratch__40__ =
    1e0 / (2e0 * (acegen_scratch__23__ * acegen_scratch__23__));
  acegen_scratch__32__ = (-(acegen_scratch__14__ * acegen_scratch__20__) +
                          acegen_scratch__13__ * acegen_scratch__22__) *
                         acegen_scratch__40__;
  acegen_scratch__31__ = (acegen_scratch__15__ * acegen_scratch__20__ -
                          acegen_scratch__14__ * acegen_scratch__22__) *
                         acegen_scratch__40__;
  acegen_scratch__30__ = (-(acegen_scratch__14__ * acegen_scratch__19__) +
                          acegen_scratch__13__ * acegen_scratch__21__) *
                         acegen_scratch__40__;
  acegen_scratch__29__ = (acegen_scratch__15__ * acegen_scratch__19__ -
                          acegen_scratch__14__ * acegen_scratch__21__) *
                         acegen_scratch__40__;
  acegen_scratch__34__ = 2e0 * lambda * log(acegen_scratch__23__);
  valueOut[0]          = 0e0;
  valueOut[1]          = 0e0;
  gradientOut[0][0] =
    acegen_scratch__29__ * acegen_scratch__34__ +
    (acegen_scratch__19__ - 2e0 * acegen_scratch__29__) * acegen_scratch__41__;
  gradientOut[0][1] =
    acegen_scratch__30__ * acegen_scratch__34__ +
    (acegen_scratch__21__ - 2e0 * acegen_scratch__30__) * acegen_scratch__41__;
  gradientOut[1][0] =
    acegen_scratch__31__ * acegen_scratch__34__ +
    (acegen_scratch__20__ - 2e0 * acegen_scratch__31__) * acegen_scratch__41__;
  gradientOut[1][1] =
    acegen_scratch__32__ * acegen_scratch__34__ +
    (acegen_scratch__22__ - 2e0 * acegen_scratch__32__) * acegen_scratch__41__;
}


template <>
template <typename Number>
inline void
SolidModel<2>::tangent(const Tensor<2, dim, Number> &graduIn,
                       const Tensor<2, dim, Number> &gradduIn,
                       Tensor<2, dim, Number> &      gradientOut,
                       const Number &                mu,
                       const Number &                lambda)
{
  Number acegen_scratch__10__, acegen_scratch__100__, acegen_scratch__102__,
    acegen_scratch__103__, acegen_scratch__104__, acegen_scratch__107__,
    acegen_scratch__11__, acegen_scratch__113__, acegen_scratch__12__,
    acegen_scratch__122__, acegen_scratch__123__, acegen_scratch__124__,
    acegen_scratch__125__, acegen_scratch__126__, acegen_scratch__135__,
    acegen_scratch__136__, acegen_scratch__137__, acegen_scratch__138__,
    acegen_scratch__139__, acegen_scratch__140__, acegen_scratch__141__,
    acegen_scratch__15__, acegen_scratch__16__, acegen_scratch__17__,
    acegen_scratch__18__, acegen_scratch__19__, acegen_scratch__20__,
    acegen_scratch__21__, acegen_scratch__25__, acegen_scratch__26__,
    acegen_scratch__27__, acegen_scratch__28__, acegen_scratch__40__,
    acegen_scratch__44__, acegen_scratch__45__, acegen_scratch__46__,
    acegen_scratch__47__, acegen_scratch__53__, acegen_scratch__59__,
    acegen_scratch__65__, acegen_scratch__70__, acegen_scratch__9__,
    acegen_scratch__97__, acegen_scratch__98__;
  acegen_scratch__9__   = gradduIn[0][0];
  acegen_scratch__97__  = 2e0 * acegen_scratch__9__;
  acegen_scratch__10__  = gradduIn[0][1];
  acegen_scratch__98__  = 2e0 * acegen_scratch__10__;
  acegen_scratch__11__  = gradduIn[1][0];
  acegen_scratch__100__ = 2e0 * acegen_scratch__11__;
  acegen_scratch__12__  = gradduIn[1][1];
  acegen_scratch__102__ = 2e0 * acegen_scratch__12__;
  acegen_scratch__140__ = mu / 2e0;
  acegen_scratch__139__ = 2e0 * lambda;
  acegen_scratch__15__  = 1e0 + graduIn[0][0];
  acegen_scratch__25__  = 2e0 * acegen_scratch__15__;
  acegen_scratch__16__  = graduIn[0][1];
  acegen_scratch__27__  = 2e0 * acegen_scratch__16__;
  acegen_scratch__17__  = graduIn[1][0];
  acegen_scratch__103__ = acegen_scratch__100__ * acegen_scratch__17__ +
                          acegen_scratch__15__ * acegen_scratch__97__;
  acegen_scratch__26__  = 2e0 * acegen_scratch__17__;
  acegen_scratch__18__  = 1e0 + graduIn[1][1];
  acegen_scratch__104__ = acegen_scratch__10__ * acegen_scratch__15__ +
                          acegen_scratch__12__ * acegen_scratch__17__ +
                          acegen_scratch__11__ * acegen_scratch__18__ +
                          acegen_scratch__16__ * acegen_scratch__9__;
  acegen_scratch__107__ = acegen_scratch__102__ * acegen_scratch__18__ +
                          acegen_scratch__16__ * acegen_scratch__98__;
  acegen_scratch__28__ = 2e0 * acegen_scratch__18__;
  acegen_scratch__19__ = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__17__ * acegen_scratch__17__);
  acegen_scratch__20__ = acegen_scratch__15__ * acegen_scratch__16__ +
                         acegen_scratch__17__ * acegen_scratch__18__;
  acegen_scratch__65__ = -(acegen_scratch__20__ * acegen_scratch__25__) +
                         acegen_scratch__19__ * acegen_scratch__27__;
  acegen_scratch__53__ = -(acegen_scratch__20__ * acegen_scratch__26__) +
                         acegen_scratch__19__ * acegen_scratch__28__;
  acegen_scratch__21__ = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__135__ = -(acegen_scratch__20__ * acegen_scratch__20__) +
                          acegen_scratch__19__ * acegen_scratch__21__;
  acegen_scratch__137__ = sqrt(acegen_scratch__135__);
  acegen_scratch__136__ = 1e0 / (2e0 * acegen_scratch__137__);
  acegen_scratch__122__ = acegen_scratch__136__ *
                          (acegen_scratch__107__ * acegen_scratch__19__ -
                           2e0 * acegen_scratch__104__ * acegen_scratch__20__ +
                           acegen_scratch__103__ * acegen_scratch__21__);
  acegen_scratch__113__ =
    -0.5e0 * acegen_scratch__122__ / acegen_scratch__135__;
  acegen_scratch__70__ = acegen_scratch__21__ * acegen_scratch__25__ -
                         acegen_scratch__20__ * acegen_scratch__27__;
  acegen_scratch__59__ = acegen_scratch__21__ * acegen_scratch__26__ -
                         acegen_scratch__20__ * acegen_scratch__28__;
  acegen_scratch__47__ = acegen_scratch__136__ * acegen_scratch__53__;
  acegen_scratch__46__ = acegen_scratch__136__ * acegen_scratch__59__;
  acegen_scratch__45__ = acegen_scratch__136__ * acegen_scratch__65__;
  acegen_scratch__44__ = acegen_scratch__136__ * acegen_scratch__70__;
  acegen_scratch__138__ =
    -(acegen_scratch__122__ / (acegen_scratch__137__ * acegen_scratch__137__));
  acegen_scratch__141__ = (acegen_scratch__122__ * acegen_scratch__139__) /
                          (acegen_scratch__137__ * acegen_scratch__137__);
  acegen_scratch__126__ =
    acegen_scratch__138__ * acegen_scratch__44__ +
    (acegen_scratch__113__ * acegen_scratch__70__ +
     acegen_scratch__136__ * (acegen_scratch__107__ * acegen_scratch__25__ -
                              acegen_scratch__104__ * acegen_scratch__27__ +
                              acegen_scratch__21__ * acegen_scratch__97__ -
                              acegen_scratch__20__ * acegen_scratch__98__)) /
      acegen_scratch__137__;
  acegen_scratch__125__ =
    acegen_scratch__138__ * acegen_scratch__45__ +
    (acegen_scratch__113__ * acegen_scratch__65__ +
     acegen_scratch__136__ * (-(acegen_scratch__104__ * acegen_scratch__25__) +
                              acegen_scratch__103__ * acegen_scratch__27__ -
                              acegen_scratch__20__ * acegen_scratch__97__ +
                              acegen_scratch__19__ * acegen_scratch__98__)) /
      acegen_scratch__137__;
  acegen_scratch__124__ =
    acegen_scratch__138__ * acegen_scratch__46__ +
    (acegen_scratch__136__ * (-(acegen_scratch__102__ * acegen_scratch__20__) +
                              acegen_scratch__100__ * acegen_scratch__21__ +
                              acegen_scratch__107__ * acegen_scratch__26__ -
                              acegen_scratch__104__ * acegen_scratch__28__) +
     acegen_scratch__113__ * acegen_scratch__59__) /
      acegen_scratch__137__;
  acegen_scratch__123__ =
    acegen_scratch__138__ * acegen_scratch__47__ +
    (acegen_scratch__136__ * (acegen_scratch__102__ * acegen_scratch__19__ -
                              acegen_scratch__100__ * acegen_scratch__20__ -
                              acegen_scratch__104__ * acegen_scratch__26__ +
                              acegen_scratch__103__ * acegen_scratch__28__) +
     acegen_scratch__113__ * acegen_scratch__53__) /
      acegen_scratch__137__;
  acegen_scratch__40__ = acegen_scratch__139__ * log(acegen_scratch__137__);
  gradientOut[0][0]    = acegen_scratch__126__ * acegen_scratch__40__ +
                      acegen_scratch__141__ * acegen_scratch__44__ +
                      acegen_scratch__140__ *
                        (-2e0 * acegen_scratch__126__ + acegen_scratch__97__);
  gradientOut[0][1] = acegen_scratch__125__ * acegen_scratch__40__ +
                      acegen_scratch__141__ * acegen_scratch__45__ +
                      acegen_scratch__140__ *
                        (-2e0 * acegen_scratch__125__ + acegen_scratch__98__);
  gradientOut[1][0] = (acegen_scratch__100__ - 2e0 * acegen_scratch__124__) *
                        acegen_scratch__140__ +
                      acegen_scratch__124__ * acegen_scratch__40__ +
                      acegen_scratch__141__ * acegen_scratch__46__;
  gradientOut[1][1] = (acegen_scratch__102__ - 2e0 * acegen_scratch__123__) *
                        acegen_scratch__140__ +
                      acegen_scratch__123__ * acegen_scratch__40__ +
                      acegen_scratch__141__ * acegen_scratch__47__;
}



template <>
template <typename Number>
inline void
SolidModel<2>::cache(const Tensor<2, dim, Number> &graduIn,
                     Tensor<1, dim, Number> &      valueOut,
                     Tensor<2, dim, Number> &      gradientOut,
                     ArrayView<Number> &           cacheOut,
                     const Number &                mu,
                     const Number &                lambda)
{
  Number acegen_scratch__10__, acegen_scratch__100__, acegen_scratch__101__,
    acegen_scratch__11__, acegen_scratch__12__, acegen_scratch__16__,
    acegen_scratch__17__, acegen_scratch__19__, acegen_scratch__29__,
    acegen_scratch__30__, acegen_scratch__31__, acegen_scratch__33__,
    acegen_scratch__36__, acegen_scratch__37__, acegen_scratch__39__,
    acegen_scratch__45__, acegen_scratch__46__, acegen_scratch__47__,
    acegen_scratch__49__, acegen_scratch__51__, acegen_scratch__65__,
    acegen_scratch__66__, acegen_scratch__7__, acegen_scratch__9__,
    acegen_scratch__91__, acegen_scratch__92__, acegen_scratch__93__,
    acegen_scratch__94__, acegen_scratch__95__, acegen_scratch__96__,
    acegen_scratch__97__, acegen_scratch__98__, acegen_scratch__99__;
  acegen_scratch__7__  = mu;
  acegen_scratch__97__ = 2e0 * lambda;
  acegen_scratch__9__  = 1e0 + graduIn[0][0];
  acegen_scratch__10__ = graduIn[0][1];
  acegen_scratch__11__ = graduIn[1][0];
  acegen_scratch__12__ = 1e0 + graduIn[1][1];
  acegen_scratch__16__ = (acegen_scratch__11__ * acegen_scratch__11__) +
                         (acegen_scratch__9__ * acegen_scratch__9__);
  acegen_scratch__17__ = (acegen_scratch__10__ * acegen_scratch__10__) +
                         (acegen_scratch__12__ * acegen_scratch__12__);
  acegen_scratch__19__ = acegen_scratch__11__ * acegen_scratch__12__ +
                         acegen_scratch__10__ * acegen_scratch__9__;
  acegen_scratch__96__ = -0.7071067811865476e0 * acegen_scratch__19__;
  acegen_scratch__92__ = sqrt(acegen_scratch__16__ * acegen_scratch__17__ -
                              (acegen_scratch__19__ * acegen_scratch__19__));
  acegen_scratch__91__ = 1e0 / (2e0 * acegen_scratch__92__);
  acegen_scratch__47__ = acegen_scratch__96__ / acegen_scratch__92__;
  acegen_scratch__46__ = acegen_scratch__16__ * acegen_scratch__91__;
  acegen_scratch__45__ = acegen_scratch__17__ * acegen_scratch__91__;
  acegen_scratch__49__ = 1e0 / (acegen_scratch__92__ * acegen_scratch__92__);
  acegen_scratch__101__ =
    -((acegen_scratch__47__ * acegen_scratch__47__) * acegen_scratch__49__);
  acegen_scratch__95__  = acegen_scratch__45__ * acegen_scratch__49__;
  acegen_scratch__94__  = -2e0 * acegen_scratch__49__;
  acegen_scratch__51__  = -(acegen_scratch__47__ * acegen_scratch__49__);
  acegen_scratch__93__  = 2e0 * acegen_scratch__51__;
  acegen_scratch__31__  = acegen_scratch__47__ / acegen_scratch__92__;
  acegen_scratch__30__  = acegen_scratch__46__ / acegen_scratch__92__;
  acegen_scratch__99__  = 2e0 * acegen_scratch__30__;
  acegen_scratch__66__  = acegen_scratch__30__ * acegen_scratch__97__;
  acegen_scratch__29__  = acegen_scratch__45__ / acegen_scratch__92__;
  acegen_scratch__98__  = 2e0 * acegen_scratch__29__;
  acegen_scratch__65__  = acegen_scratch__29__ * acegen_scratch__97__;
  acegen_scratch__33__  = acegen_scratch__97__ * log(acegen_scratch__92__);
  acegen_scratch__100__ = acegen_scratch__33__ - acegen_scratch__7__;
  acegen_scratch__36__  = acegen_scratch__7__ * (1e0 - acegen_scratch__98__) +
                         acegen_scratch__33__ * acegen_scratch__98__;
  acegen_scratch__37__ = acegen_scratch__7__ * (1e0 - acegen_scratch__99__) +
                         acegen_scratch__33__ * acegen_scratch__99__;
  acegen_scratch__39__ =
    0.1414213562373095e1 * acegen_scratch__100__ * acegen_scratch__31__;
  valueOut[0]       = 0e0;
  valueOut[1]       = 1000e0;
  gradientOut[0][0] = acegen_scratch__10__ * acegen_scratch__39__ +
                      acegen_scratch__36__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__10__ * acegen_scratch__37__ +
                      acegen_scratch__39__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__11__ * acegen_scratch__36__ +
                      acegen_scratch__12__ * acegen_scratch__39__;
  gradientOut[1][1] = acegen_scratch__12__ * acegen_scratch__37__ +
                      acegen_scratch__11__ * acegen_scratch__39__;
  cacheOut[0] = 4e0 * (acegen_scratch__29__ * acegen_scratch__65__ +
                       acegen_scratch__100__ *
                         (acegen_scratch__45__ * acegen_scratch__45__) *
                         acegen_scratch__94__);
  cacheOut[1] =
    4e0 * (acegen_scratch__30__ * acegen_scratch__65__ +
           acegen_scratch__100__ *
             (-(acegen_scratch__46__ * acegen_scratch__95__) -
              (-acegen_scratch__91__ +
               (acegen_scratch__16__ * acegen_scratch__95__) / 2e0) /
                acegen_scratch__92__));
  cacheOut[2] =
    4e0 * (acegen_scratch__31__ * acegen_scratch__65__ +
           acegen_scratch__100__ * acegen_scratch__45__ * acegen_scratch__93__);
  cacheOut[3] = 4e0 * (acegen_scratch__30__ * acegen_scratch__66__ +
                       acegen_scratch__100__ *
                         (acegen_scratch__46__ * acegen_scratch__46__) *
                         acegen_scratch__94__);
  cacheOut[4] =
    4e0 * (acegen_scratch__31__ * acegen_scratch__66__ +
           acegen_scratch__100__ * acegen_scratch__46__ * acegen_scratch__93__);
  cacheOut[5] =
    4e0 *
    (acegen_scratch__100__ *
       (acegen_scratch__101__ +
        (-acegen_scratch__91__ + acegen_scratch__51__ * acegen_scratch__96__) /
          acegen_scratch__92__) -
     acegen_scratch__101__ * acegen_scratch__97__);
  cacheOut[6] = acegen_scratch__36__;
  cacheOut[7] = acegen_scratch__39__;
  cacheOut[8] = acegen_scratch__37__;
}



template <>
template <typename Number>
inline void
SolidModel<2>::cache_deformed(const Tensor<2, dim, Number> &graduIn,
                              Tensor<2, dim, Number> &      gradientOut,
                              ArrayView<Number> &           cacheOut,
                              const Number &                mu,
                              const Number &                lambda)
{
  AssertThrow(false, ExcNotImplemented());
}


// =====
// dim=3
// =====


template <>
template <typename Number>
inline void
SolidModel<3>::residual([[maybe_unused]] const Tensor<1, dim, Number> &uIn,
                        const Tensor<2, dim, Number> &                 graduIn,
                        Tensor<1, dim, Number> &                       valueOut,
                        Tensor<2, dim, Number> &gradientOut,
                        const Number &          mu,
                        const Number &          lambda)
{
  Number acegen_scratch__13__, acegen_scratch__15__, acegen_scratch__16__,
    acegen_scratch__17__, acegen_scratch__18__, acegen_scratch__19__,
    acegen_scratch__20__, acegen_scratch__21__, acegen_scratch__22__,
    acegen_scratch__23__, acegen_scratch__24__, acegen_scratch__25__,
    acegen_scratch__26__, acegen_scratch__27__, acegen_scratch__28__,
    acegen_scratch__29__, acegen_scratch__34__, acegen_scratch__35__,
    acegen_scratch__37__, acegen_scratch__38__, acegen_scratch__39__,
    acegen_scratch__40__, acegen_scratch__41__, acegen_scratch__43__,
    acegen_scratch__44__, acegen_scratch__58__, acegen_scratch__59__,
    acegen_scratch__60__, acegen_scratch__61__, acegen_scratch__62__,
    acegen_scratch__63__;
  acegen_scratch__13__ = mu;
  acegen_scratch__40__ = acegen_scratch__13__ / 2e0;
  acegen_scratch__15__ = 1e0 + graduIn[0][0];
  acegen_scratch__16__ = graduIn[0][1];
  acegen_scratch__17__ = graduIn[0][2];
  acegen_scratch__18__ = graduIn[1][0];
  acegen_scratch__19__ = 1e0 + graduIn[1][1];
  acegen_scratch__20__ = graduIn[1][2];
  acegen_scratch__21__ = graduIn[2][0];
  acegen_scratch__22__ = graduIn[2][1];
  acegen_scratch__23__ = 1e0 + graduIn[2][2];
  acegen_scratch__24__ = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__18__ * acegen_scratch__18__) +
                         (acegen_scratch__21__ * acegen_scratch__21__);
  acegen_scratch__25__ = acegen_scratch__15__ * acegen_scratch__16__ +
                         acegen_scratch__18__ * acegen_scratch__19__ +
                         acegen_scratch__21__ * acegen_scratch__22__;
  acegen_scratch__34__ = (acegen_scratch__25__ * acegen_scratch__25__);
  acegen_scratch__26__ = acegen_scratch__15__ * acegen_scratch__17__ +
                         acegen_scratch__18__ * acegen_scratch__20__ +
                         acegen_scratch__21__ * acegen_scratch__23__;
  acegen_scratch__58__ = 2e0 * acegen_scratch__26__;
  acegen_scratch__41__ = (acegen_scratch__26__ * acegen_scratch__26__);
  acegen_scratch__38__ = acegen_scratch__25__ * acegen_scratch__58__;
  acegen_scratch__27__ = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__19__ * acegen_scratch__19__) +
                         (acegen_scratch__22__ * acegen_scratch__22__);
  acegen_scratch__28__ = acegen_scratch__16__ * acegen_scratch__17__ +
                         acegen_scratch__19__ * acegen_scratch__20__ +
                         acegen_scratch__22__ * acegen_scratch__23__;
  acegen_scratch__59__ = 2e0 * acegen_scratch__28__;
  acegen_scratch__29__ = (acegen_scratch__17__ * acegen_scratch__17__) +
                         (acegen_scratch__20__ * acegen_scratch__20__) +
                         (acegen_scratch__23__ * acegen_scratch__23__);
  acegen_scratch__60__ = -(acegen_scratch__28__ * acegen_scratch__28__) +
                         acegen_scratch__27__ * acegen_scratch__29__;
  acegen_scratch__35__ = -(acegen_scratch__29__ * acegen_scratch__34__) +
                         acegen_scratch__28__ * acegen_scratch__38__ -
                         acegen_scratch__27__ * acegen_scratch__41__ +
                         acegen_scratch__24__ * acegen_scratch__60__;
  acegen_scratch__37__ =
    -0.5e0 *
    (acegen_scratch__13__ - 2e0 * lambda * log(sqrt(acegen_scratch__35__))) /
    acegen_scratch__35__;
  acegen_scratch__63__ =
    2e0 *
    ((acegen_scratch__24__ * acegen_scratch__27__ - acegen_scratch__34__) *
       acegen_scratch__37__ +
     acegen_scratch__40__);
  acegen_scratch__39__ =
    acegen_scratch__37__ *
    (acegen_scratch__38__ - acegen_scratch__24__ * acegen_scratch__59__);
  acegen_scratch__62__ =
    2e0 * (acegen_scratch__40__ +
           acegen_scratch__37__ * (acegen_scratch__24__ * acegen_scratch__29__ -
                                   acegen_scratch__41__));
  acegen_scratch__43__ =
    acegen_scratch__37__ * (-(acegen_scratch__27__ * acegen_scratch__58__) +
                            acegen_scratch__25__ * acegen_scratch__59__);
  acegen_scratch__44__ = 2e0 *
                         (acegen_scratch__26__ * acegen_scratch__28__ -
                          acegen_scratch__25__ * acegen_scratch__29__) *
                         acegen_scratch__37__;
  acegen_scratch__61__ =
    2e0 * (acegen_scratch__40__ + acegen_scratch__37__ * acegen_scratch__60__);
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  valueOut[2]       = 0e0;
  gradientOut[0][0] = acegen_scratch__17__ * acegen_scratch__43__ +
                      acegen_scratch__16__ * acegen_scratch__44__ +
                      acegen_scratch__15__ * acegen_scratch__61__;
  gradientOut[0][1] = acegen_scratch__17__ * acegen_scratch__39__ +
                      acegen_scratch__15__ * acegen_scratch__44__ +
                      acegen_scratch__16__ * acegen_scratch__62__;
  gradientOut[0][2] = acegen_scratch__16__ * acegen_scratch__39__ +
                      acegen_scratch__15__ * acegen_scratch__43__ +
                      acegen_scratch__17__ * acegen_scratch__63__;
  gradientOut[1][0] = acegen_scratch__20__ * acegen_scratch__43__ +
                      acegen_scratch__19__ * acegen_scratch__44__ +
                      acegen_scratch__18__ * acegen_scratch__61__;
  gradientOut[1][1] = acegen_scratch__20__ * acegen_scratch__39__ +
                      acegen_scratch__18__ * acegen_scratch__44__ +
                      acegen_scratch__19__ * acegen_scratch__62__;
  gradientOut[1][2] = acegen_scratch__19__ * acegen_scratch__39__ +
                      acegen_scratch__18__ * acegen_scratch__43__ +
                      acegen_scratch__20__ * acegen_scratch__63__;
  gradientOut[2][0] = acegen_scratch__23__ * acegen_scratch__43__ +
                      acegen_scratch__22__ * acegen_scratch__44__ +
                      acegen_scratch__21__ * acegen_scratch__61__;
  gradientOut[2][1] = acegen_scratch__23__ * acegen_scratch__39__ +
                      acegen_scratch__21__ * acegen_scratch__44__ +
                      acegen_scratch__22__ * acegen_scratch__62__;
  gradientOut[2][2] = acegen_scratch__22__ * acegen_scratch__39__ +
                      acegen_scratch__21__ * acegen_scratch__43__ +
                      acegen_scratch__23__ * acegen_scratch__63__;
}



template <>
template <typename Number>
inline void
SolidModel<3>::tangent(const Tensor<2, dim, Number> &graduIn,
                       const Tensor<2, dim, Number> &gradduIn,
                       Tensor<2, dim, Number> &      gradientOut,
                       const Number &                mu,
                       const Number &                lambda)
{
  Number acegen_scratch__10__, acegen_scratch__100__, acegen_scratch__101__,
    acegen_scratch__102__, acegen_scratch__103__, acegen_scratch__105__,
    acegen_scratch__107__, acegen_scratch__109__, acegen_scratch__11__,
    acegen_scratch__110__, acegen_scratch__111__, acegen_scratch__12__,
    acegen_scratch__124__, acegen_scratch__126__, acegen_scratch__127__,
    acegen_scratch__128__, acegen_scratch__129__, acegen_scratch__13__,
    acegen_scratch__130__, acegen_scratch__131__, acegen_scratch__133__,
    acegen_scratch__14__, acegen_scratch__15__, acegen_scratch__16__,
    acegen_scratch__17__, acegen_scratch__18__, acegen_scratch__19__,
    acegen_scratch__20__, acegen_scratch__21__, acegen_scratch__22__,
    acegen_scratch__23__, acegen_scratch__24__, acegen_scratch__25__,
    acegen_scratch__26__, acegen_scratch__27__, acegen_scratch__28__,
    acegen_scratch__29__, acegen_scratch__36__, acegen_scratch__37__,
    acegen_scratch__38__, acegen_scratch__42__, acegen_scratch__43__,
    acegen_scratch__44__, acegen_scratch__47__, acegen_scratch__48__,
    acegen_scratch__49__, acegen_scratch__50__, acegen_scratch__51__,
    acegen_scratch__52__, acegen_scratch__53__, acegen_scratch__54__,
    acegen_scratch__55__, acegen_scratch__56__, acegen_scratch__58__,
    acegen_scratch__61__, acegen_scratch__64__, acegen_scratch__65__,
    acegen_scratch__66__, acegen_scratch__77__, acegen_scratch__78__,
    acegen_scratch__79__, acegen_scratch__80__, acegen_scratch__81__,
    acegen_scratch__82__, acegen_scratch__87__, acegen_scratch__88__,
    acegen_scratch__90__, acegen_scratch__95__, acegen_scratch__96__,
    acegen_scratch__97__, acegen_scratch__98__, acegen_scratch__99__;
  acegen_scratch__10__ = gradduIn[0][0];
  acegen_scratch__11__ = gradduIn[0][1];
  acegen_scratch__12__ = gradduIn[0][2];
  acegen_scratch__13__ = gradduIn[1][0];
  acegen_scratch__14__ = gradduIn[1][1];
  acegen_scratch__15__ = gradduIn[1][2];
  acegen_scratch__16__ = gradduIn[2][0];
  acegen_scratch__17__ = gradduIn[2][1];
  acegen_scratch__18__ = gradduIn[2][2];
  acegen_scratch__19__ = mu;
  acegen_scratch__51__ = acegen_scratch__19__ / 2e0;
  acegen_scratch__20__ = lambda;
  acegen_scratch__21__ = 1e0 + graduIn[0][0];
  acegen_scratch__22__ = graduIn[0][1];
  acegen_scratch__23__ = graduIn[0][2];
  acegen_scratch__24__ = graduIn[1][0];
  acegen_scratch__25__ = 1e0 + graduIn[1][1];
  acegen_scratch__26__ = graduIn[1][2];
  acegen_scratch__27__ = graduIn[2][0];
  acegen_scratch__77__ = 2e0 * (acegen_scratch__10__ * acegen_scratch__21__ +
                                acegen_scratch__13__ * acegen_scratch__24__ +
                                acegen_scratch__16__ * acegen_scratch__27__);
  acegen_scratch__28__ = graduIn[2][1];
  acegen_scratch__80__ = 2e0 * (acegen_scratch__11__ * acegen_scratch__22__ +
                                acegen_scratch__14__ * acegen_scratch__25__ +
                                acegen_scratch__17__ * acegen_scratch__28__);
  acegen_scratch__78__ = acegen_scratch__11__ * acegen_scratch__21__ +
                         acegen_scratch__10__ * acegen_scratch__22__ +
                         acegen_scratch__14__ * acegen_scratch__24__ +
                         acegen_scratch__13__ * acegen_scratch__25__ +
                         acegen_scratch__17__ * acegen_scratch__27__ +
                         acegen_scratch__16__ * acegen_scratch__28__;
  acegen_scratch__29__ = 1e0 + graduIn[2][2];
  acegen_scratch__82__ = 2e0 * (acegen_scratch__12__ * acegen_scratch__23__ +
                                acegen_scratch__15__ * acegen_scratch__26__ +
                                acegen_scratch__18__ * acegen_scratch__29__);
  acegen_scratch__81__ = acegen_scratch__12__ * acegen_scratch__22__ +
                         acegen_scratch__11__ * acegen_scratch__23__ +
                         acegen_scratch__15__ * acegen_scratch__25__ +
                         acegen_scratch__14__ * acegen_scratch__26__ +
                         acegen_scratch__18__ * acegen_scratch__28__ +
                         acegen_scratch__17__ * acegen_scratch__29__;
  acegen_scratch__79__ = acegen_scratch__12__ * acegen_scratch__21__ +
                         acegen_scratch__10__ * acegen_scratch__23__ +
                         acegen_scratch__15__ * acegen_scratch__24__ +
                         acegen_scratch__13__ * acegen_scratch__26__ +
                         acegen_scratch__18__ * acegen_scratch__27__ +
                         acegen_scratch__16__ * acegen_scratch__29__;
  acegen_scratch__36__ = (acegen_scratch__21__ * acegen_scratch__21__) +
                         (acegen_scratch__24__ * acegen_scratch__24__) +
                         (acegen_scratch__27__ * acegen_scratch__27__);
  acegen_scratch__37__ = (acegen_scratch__22__ * acegen_scratch__22__) +
                         (acegen_scratch__25__ * acegen_scratch__25__) +
                         (acegen_scratch__28__ * acegen_scratch__28__);
  acegen_scratch__38__ = (acegen_scratch__23__ * acegen_scratch__23__) +
                         (acegen_scratch__26__ * acegen_scratch__26__) +
                         (acegen_scratch__29__ * acegen_scratch__29__);
  acegen_scratch__42__ = acegen_scratch__21__ * acegen_scratch__22__ +
                         acegen_scratch__24__ * acegen_scratch__25__ +
                         acegen_scratch__27__ * acegen_scratch__28__;
  acegen_scratch__126__ = 2e0 * acegen_scratch__42__;
  acegen_scratch__131__ = acegen_scratch__37__ * acegen_scratch__77__ -
                          acegen_scratch__126__ * acegen_scratch__78__ +
                          acegen_scratch__36__ * acegen_scratch__80__;
  acegen_scratch__55__ = (acegen_scratch__42__ * acegen_scratch__42__);
  acegen_scratch__101__ =
    acegen_scratch__36__ * acegen_scratch__37__ - acegen_scratch__55__;
  acegen_scratch__43__ = acegen_scratch__21__ * acegen_scratch__23__ +
                         acegen_scratch__24__ * acegen_scratch__26__ +
                         acegen_scratch__27__ * acegen_scratch__29__;
  acegen_scratch__124__ = 2e0 * acegen_scratch__43__;
  acegen_scratch__88__  = acegen_scratch__124__ * acegen_scratch__79__;
  acegen_scratch__87__  = 2e0 * (acegen_scratch__43__ * acegen_scratch__78__ +
                                acegen_scratch__42__ * acegen_scratch__79__);
  acegen_scratch__58__  = acegen_scratch__124__ * acegen_scratch__42__;
  acegen_scratch__53__  = (acegen_scratch__43__ * acegen_scratch__43__);
  acegen_scratch__99__ =
    acegen_scratch__36__ * acegen_scratch__38__ - acegen_scratch__53__;
  acegen_scratch__44__ = acegen_scratch__22__ * acegen_scratch__23__ +
                         acegen_scratch__25__ * acegen_scratch__26__ +
                         acegen_scratch__28__ * acegen_scratch__29__;
  acegen_scratch__127__ = 2e0 * acegen_scratch__44__;
  acegen_scratch__107__ = -(acegen_scratch__126__ * acegen_scratch__38__) +
                          acegen_scratch__124__ * acegen_scratch__44__;
  acegen_scratch__105__ = -(acegen_scratch__124__ * acegen_scratch__37__) +
                          acegen_scratch__126__ * acegen_scratch__44__;
  acegen_scratch__103__ =
    -(acegen_scratch__127__ * acegen_scratch__36__) + acegen_scratch__58__;
  acegen_scratch__90__ = acegen_scratch__127__ * acegen_scratch__81__;
  acegen_scratch__48__ = (acegen_scratch__44__ * acegen_scratch__44__);
  acegen_scratch__97__ =
    acegen_scratch__37__ * acegen_scratch__38__ - acegen_scratch__48__;
  acegen_scratch__49__ = -(acegen_scratch__37__ * acegen_scratch__53__) -
                         acegen_scratch__38__ * acegen_scratch__55__ +
                         acegen_scratch__44__ * acegen_scratch__58__ +
                         acegen_scratch__36__ * acegen_scratch__97__;
  acegen_scratch__133__ = (acegen_scratch__131__ * acegen_scratch__38__ -
                           acegen_scratch__48__ * acegen_scratch__77__ -
                           acegen_scratch__53__ * acegen_scratch__80__ +
                           acegen_scratch__58__ * acegen_scratch__81__ +
                           acegen_scratch__101__ * acegen_scratch__82__ +
                           acegen_scratch__44__ * acegen_scratch__87__ -
                           acegen_scratch__37__ * acegen_scratch__88__ -
                           acegen_scratch__36__ * acegen_scratch__90__) /
                          (acegen_scratch__49__ * acegen_scratch__49__);
  acegen_scratch__129__ = acegen_scratch__133__ * acegen_scratch__20__;
  acegen_scratch__47__ =
    -acegen_scratch__19__ + acegen_scratch__20__ * log(acegen_scratch__49__);
  acegen_scratch__128__ = -(acegen_scratch__133__ * acegen_scratch__47__);
  acegen_scratch__96__  = (acegen_scratch__128__ + acegen_scratch__129__) / 2e0;
  acegen_scratch__95__ =
    0.7071067811865476e0 * (acegen_scratch__128__ + acegen_scratch__129__);
  acegen_scratch__61__ =
    (0.7071067811865476e0 * acegen_scratch__47__) / acegen_scratch__49__;
  acegen_scratch__130__ = 2e0 * acegen_scratch__61__;
  acegen_scratch__109__ =
    0.7071067811865476e0 *
    (acegen_scratch__130__ * (-(acegen_scratch__38__ * acegen_scratch__78__) +
                              acegen_scratch__44__ * acegen_scratch__79__ +
                              acegen_scratch__43__ * acegen_scratch__81__ -
                              acegen_scratch__42__ * acegen_scratch__82__) +
     acegen_scratch__107__ * acegen_scratch__95__);
  acegen_scratch__110__ =
    0.7071067811865476e0 *
    (acegen_scratch__130__ * (acegen_scratch__44__ * acegen_scratch__78__ -
                              acegen_scratch__37__ * acegen_scratch__79__ -
                              acegen_scratch__43__ * acegen_scratch__80__ +
                              acegen_scratch__42__ * acegen_scratch__81__) +
     acegen_scratch__105__ * acegen_scratch__95__);
  acegen_scratch__111__ =
    0.7071067811865476e0 *
    (acegen_scratch__61__ * (-(acegen_scratch__127__ * acegen_scratch__77__) -
                             2e0 * acegen_scratch__36__ * acegen_scratch__81__ +
                             acegen_scratch__87__) +
     acegen_scratch__103__ * acegen_scratch__95__);
  acegen_scratch__52__  = 0.7071067811865476e0 * acegen_scratch__61__;
  acegen_scratch__102__ = 2e0 * (acegen_scratch__131__ * acegen_scratch__52__ +
                                 acegen_scratch__101__ * acegen_scratch__96__);
  acegen_scratch__100__ =
    2e0 * (acegen_scratch__52__ * (acegen_scratch__38__ * acegen_scratch__77__ +
                                   acegen_scratch__36__ * acegen_scratch__82__ -
                                   acegen_scratch__88__) +
           acegen_scratch__96__ * acegen_scratch__99__);
  acegen_scratch__98__ =
    2e0 * (acegen_scratch__52__ * (acegen_scratch__38__ * acegen_scratch__80__ +
                                   acegen_scratch__37__ * acegen_scratch__82__ -
                                   acegen_scratch__90__) +
           acegen_scratch__96__ * acegen_scratch__97__);
  acegen_scratch__50__ =
    2e0 * (acegen_scratch__51__ + acegen_scratch__52__ * acegen_scratch__97__);
  acegen_scratch__54__ =
    2e0 * (acegen_scratch__51__ + acegen_scratch__52__ * acegen_scratch__99__);
  acegen_scratch__56__ =
    2e0 * (acegen_scratch__51__ + acegen_scratch__101__ * acegen_scratch__52__);
  acegen_scratch__64__ = acegen_scratch__107__ * acegen_scratch__52__;
  acegen_scratch__65__ = acegen_scratch__105__ * acegen_scratch__52__;
  acegen_scratch__66__ = 1e0 * acegen_scratch__103__ * acegen_scratch__52__;
  gradientOut[0][0]    = acegen_scratch__109__ * acegen_scratch__22__ +
                      acegen_scratch__110__ * acegen_scratch__23__ +
                      acegen_scratch__10__ * acegen_scratch__50__ +
                      acegen_scratch__11__ * acegen_scratch__64__ +
                      acegen_scratch__12__ * acegen_scratch__65__ +
                      acegen_scratch__21__ * acegen_scratch__98__;
  gradientOut[0][1] = acegen_scratch__109__ * acegen_scratch__21__ +
                      acegen_scratch__100__ * acegen_scratch__22__ +
                      acegen_scratch__111__ * acegen_scratch__23__ +
                      acegen_scratch__11__ * acegen_scratch__54__ +
                      acegen_scratch__10__ * acegen_scratch__64__ +
                      acegen_scratch__12__ * acegen_scratch__66__;
  gradientOut[0][2] = acegen_scratch__110__ * acegen_scratch__21__ +
                      acegen_scratch__111__ * acegen_scratch__22__ +
                      acegen_scratch__102__ * acegen_scratch__23__ +
                      acegen_scratch__12__ * acegen_scratch__56__ +
                      acegen_scratch__10__ * acegen_scratch__65__ +
                      acegen_scratch__11__ * acegen_scratch__66__;
  gradientOut[1][0] = acegen_scratch__109__ * acegen_scratch__25__ +
                      acegen_scratch__110__ * acegen_scratch__26__ +
                      acegen_scratch__13__ * acegen_scratch__50__ +
                      acegen_scratch__14__ * acegen_scratch__64__ +
                      acegen_scratch__15__ * acegen_scratch__65__ +
                      acegen_scratch__24__ * acegen_scratch__98__;
  gradientOut[1][1] = acegen_scratch__109__ * acegen_scratch__24__ +
                      acegen_scratch__100__ * acegen_scratch__25__ +
                      acegen_scratch__111__ * acegen_scratch__26__ +
                      acegen_scratch__14__ * acegen_scratch__54__ +
                      acegen_scratch__13__ * acegen_scratch__64__ +
                      acegen_scratch__15__ * acegen_scratch__66__;
  gradientOut[1][2] = acegen_scratch__110__ * acegen_scratch__24__ +
                      acegen_scratch__111__ * acegen_scratch__25__ +
                      acegen_scratch__102__ * acegen_scratch__26__ +
                      acegen_scratch__15__ * acegen_scratch__56__ +
                      acegen_scratch__13__ * acegen_scratch__65__ +
                      acegen_scratch__14__ * acegen_scratch__66__;
  gradientOut[2][0] = acegen_scratch__109__ * acegen_scratch__28__ +
                      acegen_scratch__110__ * acegen_scratch__29__ +
                      acegen_scratch__16__ * acegen_scratch__50__ +
                      acegen_scratch__17__ * acegen_scratch__64__ +
                      acegen_scratch__18__ * acegen_scratch__65__ +
                      acegen_scratch__27__ * acegen_scratch__98__;
  gradientOut[2][1] = acegen_scratch__109__ * acegen_scratch__27__ +
                      acegen_scratch__100__ * acegen_scratch__28__ +
                      acegen_scratch__111__ * acegen_scratch__29__ +
                      acegen_scratch__17__ * acegen_scratch__54__ +
                      acegen_scratch__16__ * acegen_scratch__64__ +
                      acegen_scratch__18__ * acegen_scratch__66__;
  gradientOut[2][2] = acegen_scratch__110__ * acegen_scratch__27__ +
                      acegen_scratch__111__ * acegen_scratch__28__ +
                      acegen_scratch__102__ * acegen_scratch__29__ +
                      acegen_scratch__18__ * acegen_scratch__56__ +
                      acegen_scratch__16__ * acegen_scratch__65__ +
                      acegen_scratch__17__ * acegen_scratch__66__;
}



template <>
template <typename Number>
inline void
SolidModel<3>::cache(const Tensor<2, dim, Number> &graduIn,
                     Tensor<1, dim, Number> &      valueOut,
                     Tensor<2, dim, Number> &      gradientOut,
                     ArrayView<Number> &           cacheOut,
                     const Number &                mu,
                     const Number &                lambda)
{
  Number acegen_scratch__106__, acegen_scratch__113__, acegen_scratch__114__,
    acegen_scratch__116__, acegen_scratch__117__, acegen_scratch__118__,
    acegen_scratch__119__, acegen_scratch__120__, acegen_scratch__13__,
    acegen_scratch__136__, acegen_scratch__14__, acegen_scratch__140__,
    acegen_scratch__141__, acegen_scratch__144__, acegen_scratch__15__,
    acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__170__,
    acegen_scratch__171__, acegen_scratch__172__, acegen_scratch__173__,
    acegen_scratch__174__, acegen_scratch__175__, acegen_scratch__176__,
    acegen_scratch__178__, acegen_scratch__179__, acegen_scratch__18__,
    acegen_scratch__180__, acegen_scratch__181__, acegen_scratch__182__,
    acegen_scratch__183__, acegen_scratch__184__, acegen_scratch__185__,
    acegen_scratch__186__, acegen_scratch__187__, acegen_scratch__188__,
    acegen_scratch__19__, acegen_scratch__20__, acegen_scratch__21__,
    acegen_scratch__22__, acegen_scratch__23__, acegen_scratch__30__,
    acegen_scratch__31__, acegen_scratch__32__, acegen_scratch__33__,
    acegen_scratch__34__, acegen_scratch__35__, acegen_scratch__42__,
    acegen_scratch__44__, acegen_scratch__45__, acegen_scratch__46__,
    acegen_scratch__47__, acegen_scratch__48__, acegen_scratch__49__,
    acegen_scratch__50__, acegen_scratch__51__, acegen_scratch__53__,
    acegen_scratch__56__, acegen_scratch__59__, acegen_scratch__60__,
    acegen_scratch__61__, acegen_scratch__77__, acegen_scratch__78__,
    acegen_scratch__79__, acegen_scratch__80__, acegen_scratch__81__,
    acegen_scratch__82__, acegen_scratch__84__, acegen_scratch__87__;
  acegen_scratch__13__  = mu;
  acegen_scratch__46__  = acegen_scratch__13__ / 2e0;
  acegen_scratch__14__  = lambda;
  acegen_scratch__15__  = 1e0 + graduIn[0][0];
  acegen_scratch__16__  = graduIn[0][1];
  acegen_scratch__17__  = graduIn[0][2];
  acegen_scratch__18__  = graduIn[1][0];
  acegen_scratch__19__  = 1e0 + graduIn[1][1];
  acegen_scratch__20__  = graduIn[1][2];
  acegen_scratch__21__  = graduIn[2][0];
  acegen_scratch__22__  = graduIn[2][1];
  acegen_scratch__170__ = acegen_scratch__15__ * acegen_scratch__16__ +
                          acegen_scratch__18__ * acegen_scratch__19__ +
                          acegen_scratch__21__ * acegen_scratch__22__;
  acegen_scratch__23__  = 1e0 + graduIn[2][2];
  acegen_scratch__172__ = acegen_scratch__16__ * acegen_scratch__17__ +
                          acegen_scratch__19__ * acegen_scratch__20__ +
                          acegen_scratch__22__ * acegen_scratch__23__;
  acegen_scratch__171__ = acegen_scratch__15__ * acegen_scratch__17__ +
                          acegen_scratch__18__ * acegen_scratch__20__ +
                          acegen_scratch__21__ * acegen_scratch__23__;
  acegen_scratch__30__ = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__18__ * acegen_scratch__18__) +
                         (acegen_scratch__21__ * acegen_scratch__21__);
  acegen_scratch__31__ = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__19__ * acegen_scratch__19__) +
                         (acegen_scratch__22__ * acegen_scratch__22__);
  acegen_scratch__32__ = (acegen_scratch__17__ * acegen_scratch__17__) +
                         (acegen_scratch__20__ * acegen_scratch__20__) +
                         (acegen_scratch__23__ * acegen_scratch__23__);
  acegen_scratch__33__  = 0.1414213562373095e1 * acegen_scratch__172__;
  acegen_scratch__34__  = 0.1414213562373095e1 * acegen_scratch__171__;
  acegen_scratch__35__  = 0.1414213562373095e1 * acegen_scratch__170__;
  acegen_scratch__173__ = 2e0 * acegen_scratch__170__;
  acegen_scratch__50__  = (acegen_scratch__170__ * acegen_scratch__170__);
  acegen_scratch__79__ =
    acegen_scratch__30__ * acegen_scratch__31__ - acegen_scratch__50__;
  acegen_scratch__187__ = 4e0 * acegen_scratch__79__;
  acegen_scratch__174__ = 2e0 * acegen_scratch__171__;
  acegen_scratch__53__  = acegen_scratch__171__ * acegen_scratch__173__;
  acegen_scratch__80__  = -(acegen_scratch__30__ * acegen_scratch__33__) +
                         0.7071067811865476e0 * acegen_scratch__53__;
  acegen_scratch__48__ = (acegen_scratch__171__ * acegen_scratch__171__);
  acegen_scratch__78__ =
    acegen_scratch__30__ * acegen_scratch__32__ - acegen_scratch__48__;
  acegen_scratch__186__ = 4e0 * acegen_scratch__78__;
  acegen_scratch__144__ = acegen_scratch__172__ * acegen_scratch__174__ -
                          acegen_scratch__173__ * acegen_scratch__32__;
  acegen_scratch__141__ = acegen_scratch__172__ * acegen_scratch__173__ -
                          acegen_scratch__174__ * acegen_scratch__31__;
  acegen_scratch__136__ =
    -2e0 * acegen_scratch__172__ * acegen_scratch__30__ + acegen_scratch__53__;
  acegen_scratch__82__ = acegen_scratch__172__ * acegen_scratch__34__ -
                         acegen_scratch__32__ * acegen_scratch__35__;
  acegen_scratch__81__ = -(acegen_scratch__31__ * acegen_scratch__34__) +
                         acegen_scratch__172__ * acegen_scratch__35__;
  acegen_scratch__77__ = -(acegen_scratch__172__ * acegen_scratch__172__) +
                         acegen_scratch__31__ * acegen_scratch__32__;
  acegen_scratch__44__ = -(acegen_scratch__31__ * acegen_scratch__48__) -
                         acegen_scratch__32__ * acegen_scratch__50__ +
                         acegen_scratch__172__ * acegen_scratch__53__ +
                         acegen_scratch__30__ * acegen_scratch__77__;
  acegen_scratch__188__ = 0.7071067811865476e0 / acegen_scratch__44__;
  acegen_scratch__175__ = sqrt(acegen_scratch__44__);
  acegen_scratch__176__ =
    acegen_scratch__14__ / (acegen_scratch__175__ * acegen_scratch__175__);
  acegen_scratch__180__ =
    (acegen_scratch__176__ * acegen_scratch__82__) / acegen_scratch__44__;
  acegen_scratch__182__ =
    (acegen_scratch__176__ * acegen_scratch__81__) / acegen_scratch__44__;
  acegen_scratch__106__ = acegen_scratch__176__ * acegen_scratch__80__;
  acegen_scratch__84__  = 1e0 / (acegen_scratch__44__ * acegen_scratch__44__);
  acegen_scratch__87__  = -(acegen_scratch__80__ * acegen_scratch__84__);
  acegen_scratch__42__  = -acegen_scratch__13__ + 2e0 * acegen_scratch__14__ *
                                                   log(acegen_scratch__175__);
  acegen_scratch__181__ =
    -(acegen_scratch__42__ * acegen_scratch__81__ * acegen_scratch__84__);
  acegen_scratch__179__ =
    -(acegen_scratch__42__ * acegen_scratch__82__ * acegen_scratch__84__);
  acegen_scratch__178__ = (acegen_scratch__176__ / acegen_scratch__44__ -
                           acegen_scratch__42__ * acegen_scratch__84__) /
                          2e0;
  acegen_scratch__120__ = (acegen_scratch__179__ + acegen_scratch__180__) / 2e0;
  acegen_scratch__185__ = 4e0 * acegen_scratch__120__;
  acegen_scratch__119__ = (acegen_scratch__181__ + acegen_scratch__182__) / 2e0;
  acegen_scratch__184__ = 4e0 * acegen_scratch__119__;
  acegen_scratch__118__ = (acegen_scratch__106__ / acegen_scratch__44__ +
                           acegen_scratch__42__ * acegen_scratch__87__) /
                          2e0;
  acegen_scratch__117__ = acegen_scratch__178__ * acegen_scratch__79__;
  acegen_scratch__116__ = acegen_scratch__178__ * acegen_scratch__78__;
  acegen_scratch__114__ =
    0.7071067811865476e0 * (acegen_scratch__179__ + acegen_scratch__180__);
  acegen_scratch__113__ =
    0.7071067811865476e0 * (acegen_scratch__181__ + acegen_scratch__182__);
  acegen_scratch__56__  = acegen_scratch__188__ * acegen_scratch__42__;
  acegen_scratch__140__ = -0.1414213562373095e1 * acegen_scratch__56__;
  acegen_scratch__47__  = 0.7071067811865476e0 * acegen_scratch__56__;
  acegen_scratch__183__ = 1e0 * acegen_scratch__47__;
  acegen_scratch__45__ =
    2e0 * (acegen_scratch__46__ + acegen_scratch__47__ * acegen_scratch__77__);
  acegen_scratch__49__ =
    2e0 * (acegen_scratch__46__ + acegen_scratch__47__ * acegen_scratch__78__);
  acegen_scratch__51__ =
    2e0 * (acegen_scratch__46__ + acegen_scratch__47__ * acegen_scratch__79__);
  acegen_scratch__59__ = acegen_scratch__144__ * acegen_scratch__183__;
  acegen_scratch__60__ = acegen_scratch__141__ * acegen_scratch__183__;
  acegen_scratch__61__ = 1e0 * acegen_scratch__136__ * acegen_scratch__183__;
  valueOut[0]          = 0e0;
  valueOut[1]          = 1000e0;
  valueOut[2]          = 0e0;
  gradientOut[0][0]    = acegen_scratch__15__ * acegen_scratch__45__ +
                      acegen_scratch__16__ * acegen_scratch__59__ +
                      acegen_scratch__17__ * acegen_scratch__60__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__49__ +
                      acegen_scratch__15__ * acegen_scratch__59__ +
                      acegen_scratch__17__ * acegen_scratch__61__;
  gradientOut[0][2] = acegen_scratch__17__ * acegen_scratch__51__ +
                      acegen_scratch__15__ * acegen_scratch__60__ +
                      acegen_scratch__16__ * acegen_scratch__61__;
  gradientOut[1][0] = acegen_scratch__18__ * acegen_scratch__45__ +
                      acegen_scratch__19__ * acegen_scratch__59__ +
                      acegen_scratch__20__ * acegen_scratch__60__;
  gradientOut[1][1] = acegen_scratch__19__ * acegen_scratch__49__ +
                      acegen_scratch__18__ * acegen_scratch__59__ +
                      acegen_scratch__20__ * acegen_scratch__61__;
  gradientOut[1][2] = acegen_scratch__20__ * acegen_scratch__51__ +
                      acegen_scratch__18__ * acegen_scratch__60__ +
                      acegen_scratch__19__ * acegen_scratch__61__;
  gradientOut[2][0] = acegen_scratch__21__ * acegen_scratch__45__ +
                      acegen_scratch__22__ * acegen_scratch__59__ +
                      acegen_scratch__23__ * acegen_scratch__60__;
  gradientOut[2][1] = acegen_scratch__22__ * acegen_scratch__49__ +
                      acegen_scratch__21__ * acegen_scratch__59__ +
                      acegen_scratch__23__ * acegen_scratch__61__;
  gradientOut[2][2] = acegen_scratch__23__ * acegen_scratch__51__ +
                      acegen_scratch__21__ * acegen_scratch__60__ +
                      acegen_scratch__22__ * acegen_scratch__61__;
  cacheOut[0] =
    4e0 * acegen_scratch__178__ * (acegen_scratch__77__ * acegen_scratch__77__);
  cacheOut[1]  = 4e0 * (acegen_scratch__32__ * acegen_scratch__47__ +
                       acegen_scratch__116__ * acegen_scratch__77__);
  cacheOut[2]  = 4e0 * (acegen_scratch__31__ * acegen_scratch__47__ +
                       acegen_scratch__117__ * acegen_scratch__77__);
  cacheOut[3]  = 4e0 * (-(acegen_scratch__33__ * acegen_scratch__47__) +
                       acegen_scratch__118__ * acegen_scratch__77__);
  cacheOut[4]  = acegen_scratch__184__ * acegen_scratch__77__;
  cacheOut[5]  = acegen_scratch__185__ * acegen_scratch__77__;
  cacheOut[6]  = acegen_scratch__116__ * acegen_scratch__186__;
  cacheOut[7]  = 4e0 * (acegen_scratch__30__ * acegen_scratch__47__ +
                       acegen_scratch__117__ * acegen_scratch__78__);
  cacheOut[8]  = acegen_scratch__118__ * acegen_scratch__186__;
  cacheOut[9]  = 4e0 * (-(acegen_scratch__34__ * acegen_scratch__47__) +
                       acegen_scratch__119__ * acegen_scratch__78__);
  cacheOut[10] = acegen_scratch__185__ * acegen_scratch__78__;
  cacheOut[11] = acegen_scratch__117__ * acegen_scratch__187__;
  cacheOut[12] = acegen_scratch__118__ * acegen_scratch__187__;
  cacheOut[13] = acegen_scratch__184__ * acegen_scratch__79__;
  cacheOut[14] = 4e0 * (-(acegen_scratch__35__ * acegen_scratch__47__) +
                        acegen_scratch__120__ * acegen_scratch__79__);
  cacheOut[15] =
    2e0 *
    (acegen_scratch__140__ * acegen_scratch__30__ +
     acegen_scratch__136__ *
       (acegen_scratch__106__ * acegen_scratch__188__ +
        acegen_scratch__44__ * acegen_scratch__56__ * acegen_scratch__87__));
  cacheOut[16] = 2e0 * (acegen_scratch__113__ * acegen_scratch__136__ +
                        acegen_scratch__35__ * acegen_scratch__56__);
  cacheOut[17] = 2e0 * (acegen_scratch__114__ * acegen_scratch__136__ +
                        acegen_scratch__34__ * acegen_scratch__56__);
  cacheOut[18] = 2e0 * (acegen_scratch__113__ * acegen_scratch__141__ +
                        acegen_scratch__140__ * acegen_scratch__31__);
  cacheOut[19] = 2e0 * (acegen_scratch__114__ * acegen_scratch__141__ -
                        acegen_scratch__140__ * acegen_scratch__172__);
  cacheOut[20] = 2e0 * (acegen_scratch__114__ * acegen_scratch__144__ +
                        acegen_scratch__140__ * acegen_scratch__32__);
  cacheOut[21] = acegen_scratch__45__;
  cacheOut[22] = acegen_scratch__59__;
  cacheOut[23] = acegen_scratch__60__;
  cacheOut[24] = acegen_scratch__49__;
  cacheOut[25] = acegen_scratch__61__;
  cacheOut[26] = acegen_scratch__51__;
}



template <>
template <typename Number>
inline void
SolidModel<3>::cache_deformed(const Tensor<2, dim, Number> &graduIn,
                              Tensor<2, dim, Number> &      gradientOut,
                              ArrayView<Number> &           cacheOut,
                              const Number &                mu,
                              const Number &                lambda)
{
  Number acegen_scratch__10__, acegen_scratch__103__, acegen_scratch__11__,
    acegen_scratch__110__, acegen_scratch__111__, acegen_scratch__113__,
    acegen_scratch__114__, acegen_scratch__115__, acegen_scratch__116__,
    acegen_scratch__117__, acegen_scratch__12__, acegen_scratch__13__,
    acegen_scratch__133__, acegen_scratch__137__, acegen_scratch__138__,
    acegen_scratch__14__, acegen_scratch__141__, acegen_scratch__143__,
    acegen_scratch__144__, acegen_scratch__145__, acegen_scratch__146__,
    acegen_scratch__147__, acegen_scratch__148__, acegen_scratch__149__,
    acegen_scratch__15__, acegen_scratch__150__, acegen_scratch__151__,
    acegen_scratch__152__, acegen_scratch__153__, acegen_scratch__154__,
    acegen_scratch__155__, acegen_scratch__156__, acegen_scratch__157__,
    acegen_scratch__158__, acegen_scratch__159__, acegen_scratch__16__,
    acegen_scratch__160__, acegen_scratch__161__, acegen_scratch__162__,
    acegen_scratch__163__, acegen_scratch__164__, acegen_scratch__165__,
    acegen_scratch__167__, acegen_scratch__17__, acegen_scratch__171__,
    acegen_scratch__172__, acegen_scratch__174__, acegen_scratch__175__,
    acegen_scratch__176__, acegen_scratch__178__, acegen_scratch__179__,
    acegen_scratch__18__, acegen_scratch__180__, acegen_scratch__184__,
    acegen_scratch__19__, acegen_scratch__191__, acegen_scratch__193__,
    acegen_scratch__194__, acegen_scratch__195__, acegen_scratch__196__,
    acegen_scratch__197__, acegen_scratch__198__, acegen_scratch__20__,
    acegen_scratch__27__, acegen_scratch__273__, acegen_scratch__276__,
    acegen_scratch__279__, acegen_scratch__28__, acegen_scratch__280__,
    acegen_scratch__282__, acegen_scratch__283__, acegen_scratch__284__,
    acegen_scratch__29__, acegen_scratch__30__, acegen_scratch__308__,
    acegen_scratch__31__, acegen_scratch__310__, acegen_scratch__32__,
    acegen_scratch__39__, acegen_scratch__390__, acegen_scratch__41__,
    acegen_scratch__42__, acegen_scratch__43__, acegen_scratch__44__,
    acegen_scratch__45__, acegen_scratch__46__, acegen_scratch__47__,
    acegen_scratch__48__, acegen_scratch__483__, acegen_scratch__484__,
    acegen_scratch__485__, acegen_scratch__486__, acegen_scratch__487__,
    acegen_scratch__488__, acegen_scratch__489__, acegen_scratch__491__,
    acegen_scratch__492__, acegen_scratch__493__, acegen_scratch__494__,
    acegen_scratch__495__, acegen_scratch__496__, acegen_scratch__497__,
    acegen_scratch__498__, acegen_scratch__50__, acegen_scratch__500__,
    acegen_scratch__502__, acegen_scratch__506__, acegen_scratch__508__,
    acegen_scratch__511__, acegen_scratch__512__, acegen_scratch__516__,
    acegen_scratch__518__, acegen_scratch__520__, acegen_scratch__521__,
    acegen_scratch__524__, acegen_scratch__527__, acegen_scratch__528__,
    acegen_scratch__53__, acegen_scratch__530__, acegen_scratch__532__,
    acegen_scratch__533__, acegen_scratch__535__, acegen_scratch__537__,
    acegen_scratch__541__, acegen_scratch__546__, acegen_scratch__547__,
    acegen_scratch__548__, acegen_scratch__549__, acegen_scratch__56__,
    acegen_scratch__57__, acegen_scratch__572__, acegen_scratch__573__,
    acegen_scratch__576__, acegen_scratch__58__, acegen_scratch__581__,
    acegen_scratch__59__, acegen_scratch__592__, acegen_scratch__593__,
    acegen_scratch__598__, acegen_scratch__599__, acegen_scratch__60__,
    acegen_scratch__600__, acegen_scratch__601__, acegen_scratch__61__,
    acegen_scratch__62__, acegen_scratch__63__, acegen_scratch__64__,
    acegen_scratch__65__, acegen_scratch__66__, acegen_scratch__67__,
    acegen_scratch__74__, acegen_scratch__75__, acegen_scratch__76__,
    acegen_scratch__77__, acegen_scratch__78__, acegen_scratch__79__,
    acegen_scratch__81__, acegen_scratch__84__;
  acegen_scratch__10__  = (mu);
  acegen_scratch__43__  = acegen_scratch__10__ / 2e0;
  acegen_scratch__11__  = (lambda);
  acegen_scratch__12__  = 1e0 + graduIn[0][0];
  acegen_scratch__164__ = (acegen_scratch__12__ * acegen_scratch__12__);
  acegen_scratch__541__ = 2e0 * acegen_scratch__164__;
  acegen_scratch__13__  = graduIn[0][1];
  acegen_scratch__518__ = acegen_scratch__12__ * acegen_scratch__13__;
  acegen_scratch__165__ = (acegen_scratch__13__ * acegen_scratch__13__);
  acegen_scratch__14__  = graduIn[0][2];
  acegen_scratch__533__ = acegen_scratch__13__ * acegen_scratch__14__;
  acegen_scratch__191__ = 2e0 * acegen_scratch__12__ * acegen_scratch__14__;
  acegen_scratch__598__ = acegen_scratch__191__ / 2e0;
  acegen_scratch__167__ = (acegen_scratch__14__ * acegen_scratch__14__);
  acegen_scratch__524__ = acegen_scratch__12__ * acegen_scratch__167__;
  acegen_scratch__15__  = graduIn[1][0];
  acegen_scratch__178__ = (acegen_scratch__15__ * acegen_scratch__15__);
  acegen_scratch__16__  = 1e0 + graduIn[1][1];
  acegen_scratch__547__ = acegen_scratch__15__ * acegen_scratch__16__;
  acegen_scratch__179__ = (acegen_scratch__16__ * acegen_scratch__16__);
  acegen_scratch__17__  = graduIn[1][2];
  acegen_scratch__572__ = acegen_scratch__17__ / 2e0;
  acegen_scratch__521__ = acegen_scratch__16__ * acegen_scratch__17__;
  acegen_scratch__284__ = 2e0 * acegen_scratch__15__ * acegen_scratch__17__;
  acegen_scratch__576__ = acegen_scratch__284__ / 2e0;
  acegen_scratch__546__ = 0.7071067811865476e0 * acegen_scratch__284__;
  acegen_scratch__180__ = (acegen_scratch__17__ * acegen_scratch__17__);
  acegen_scratch__18__  = graduIn[2][0];
  acegen_scratch__498__ = acegen_scratch__15__ * acegen_scratch__18__;
  acegen_scratch__390__ = acegen_scratch__18__ * acegen_scratch__546__;
  acegen_scratch__193__ = (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__19__  = graduIn[2][1];
  acegen_scratch__573__ = acegen_scratch__16__ * acegen_scratch__19__;
  acegen_scratch__520__ = acegen_scratch__18__ * acegen_scratch__19__;
  acegen_scratch__483__ =
    acegen_scratch__518__ + acegen_scratch__520__ + acegen_scratch__547__;
  acegen_scratch__194__ = (acegen_scratch__19__ * acegen_scratch__19__);
  acegen_scratch__20__  = 1e0 + graduIn[2][2];
  acegen_scratch__535__ = acegen_scratch__19__ * acegen_scratch__20__;
  acegen_scratch__512__ = acegen_scratch__17__ * acegen_scratch__20__;
  acegen_scratch__485__ =
    acegen_scratch__521__ + acegen_scratch__533__ + acegen_scratch__535__;
  acegen_scratch__283__ = 2e0 * acegen_scratch__18__ * acegen_scratch__20__;
  acegen_scratch__500__ = acegen_scratch__283__ / 2e0;
  acegen_scratch__484__ =
    acegen_scratch__500__ + acegen_scratch__576__ + acegen_scratch__598__;
  acegen_scratch__195__ = (acegen_scratch__20__ * acegen_scratch__20__);
  acegen_scratch__27__ =
    acegen_scratch__164__ + acegen_scratch__178__ + acegen_scratch__193__;
  acegen_scratch__28__ =
    acegen_scratch__165__ + acegen_scratch__179__ + acegen_scratch__194__;
  acegen_scratch__29__ =
    acegen_scratch__167__ + acegen_scratch__180__ + acegen_scratch__195__;
  acegen_scratch__30__  = 0.1414213562373095e1 * acegen_scratch__485__;
  acegen_scratch__31__  = 0.1414213562373095e1 * acegen_scratch__484__;
  acegen_scratch__32__  = 0.1414213562373095e1 * acegen_scratch__483__;
  acegen_scratch__486__ = 2e0 * acegen_scratch__483__;
  acegen_scratch__47__  = (acegen_scratch__483__ * acegen_scratch__483__);
  acegen_scratch__76__ =
    acegen_scratch__27__ * acegen_scratch__28__ - acegen_scratch__47__;
  acegen_scratch__511__ = 4e0 * acegen_scratch__76__;
  acegen_scratch__487__ = 2e0 * acegen_scratch__484__;
  acegen_scratch__50__  = acegen_scratch__484__ * acegen_scratch__486__;
  acegen_scratch__77__  = -(acegen_scratch__27__ * acegen_scratch__30__) +
                         0.7071067811865476e0 * acegen_scratch__50__;
  acegen_scratch__45__ = (acegen_scratch__484__ * acegen_scratch__484__);
  acegen_scratch__75__ =
    acegen_scratch__27__ * acegen_scratch__29__ - acegen_scratch__45__;
  acegen_scratch__506__ = 4e0 * acegen_scratch__75__;
  acegen_scratch__532__ = -2e0 * acegen_scratch__485__;
  acegen_scratch__141__ = -(acegen_scratch__29__ * acegen_scratch__486__) +
                          acegen_scratch__485__ * acegen_scratch__487__;
  acegen_scratch__138__ = acegen_scratch__485__ * acegen_scratch__486__ -
                          acegen_scratch__28__ * acegen_scratch__487__;
  acegen_scratch__133__ =
    acegen_scratch__50__ + acegen_scratch__27__ * acegen_scratch__532__;
  acegen_scratch__79__ = -(acegen_scratch__29__ * acegen_scratch__32__) +
                         acegen_scratch__31__ * acegen_scratch__485__;
  acegen_scratch__78__ = -(acegen_scratch__28__ * acegen_scratch__31__) +
                         acegen_scratch__32__ * acegen_scratch__485__;
  acegen_scratch__74__ = acegen_scratch__28__ * acegen_scratch__29__ -
                         (acegen_scratch__485__ * acegen_scratch__485__);
  acegen_scratch__502__ = 4e0 * acegen_scratch__74__;
  acegen_scratch__41__  = -(acegen_scratch__28__ * acegen_scratch__45__) -
                         acegen_scratch__29__ * acegen_scratch__47__ +
                         acegen_scratch__485__ * acegen_scratch__50__ +
                         acegen_scratch__27__ * acegen_scratch__74__;
  acegen_scratch__516__ = 0.7071067811865476e0 / acegen_scratch__41__;
  acegen_scratch__488__ = sqrt(acegen_scratch__41__);
  acegen_scratch__489__ =
    acegen_scratch__11__ / (acegen_scratch__488__ * acegen_scratch__488__);
  acegen_scratch__493__ =
    (acegen_scratch__489__ * acegen_scratch__79__) / acegen_scratch__41__;
  acegen_scratch__495__ =
    (acegen_scratch__489__ * acegen_scratch__78__) / acegen_scratch__41__;
  acegen_scratch__103__ = acegen_scratch__489__ * acegen_scratch__77__;
  acegen_scratch__81__  = 1e0 / (acegen_scratch__41__ * acegen_scratch__41__);
  acegen_scratch__84__  = -(acegen_scratch__77__ * acegen_scratch__81__);
  acegen_scratch__39__  = -acegen_scratch__10__ + 2e0 * acegen_scratch__11__ *
                                                   log(acegen_scratch__488__);
  acegen_scratch__601__ = (-2e0 * acegen_scratch__39__) / acegen_scratch__488__;
  acegen_scratch__497__ = acegen_scratch__489__ / acegen_scratch__41__ -
                          acegen_scratch__39__ * acegen_scratch__81__;
  acegen_scratch__494__ =
    -(acegen_scratch__39__ * acegen_scratch__78__ * acegen_scratch__81__);
  acegen_scratch__492__ =
    -(acegen_scratch__39__ * acegen_scratch__79__ * acegen_scratch__81__);
  acegen_scratch__491__ = acegen_scratch__497__ / 2e0;
  acegen_scratch__117__ = (acegen_scratch__492__ + acegen_scratch__493__) / 2e0;
  acegen_scratch__116__ = (acegen_scratch__494__ + acegen_scratch__495__) / 2e0;
  acegen_scratch__115__ = (acegen_scratch__103__ / acegen_scratch__41__ +
                           acegen_scratch__39__ * acegen_scratch__84__) /
                          2e0;
  acegen_scratch__114__ = acegen_scratch__491__ * acegen_scratch__76__;
  acegen_scratch__113__ = acegen_scratch__491__ * acegen_scratch__75__;
  acegen_scratch__111__ =
    0.7071067811865476e0 * (acegen_scratch__492__ + acegen_scratch__493__);
  acegen_scratch__110__ =
    0.7071067811865476e0 * (acegen_scratch__494__ + acegen_scratch__495__);
  acegen_scratch__53__  = acegen_scratch__39__ * acegen_scratch__516__;
  acegen_scratch__137__ = -0.1414213562373095e1 * acegen_scratch__53__;
  acegen_scratch__44__  = 0.7071067811865476e0 * acegen_scratch__53__;
  acegen_scratch__508__ = -4e0 * acegen_scratch__44__;
  acegen_scratch__496__ = 1e0 * acegen_scratch__44__;
  acegen_scratch__42__ =
    2e0 * (acegen_scratch__43__ + acegen_scratch__44__ * acegen_scratch__74__);
  acegen_scratch__46__ =
    2e0 * (acegen_scratch__43__ + acegen_scratch__44__ * acegen_scratch__75__);
  acegen_scratch__48__ =
    2e0 * (acegen_scratch__43__ + acegen_scratch__44__ * acegen_scratch__76__);
  acegen_scratch__56__ = acegen_scratch__141__ * acegen_scratch__496__;
  acegen_scratch__57__ = acegen_scratch__138__ * acegen_scratch__496__;
  acegen_scratch__58__ = 1e0 * acegen_scratch__133__ * acegen_scratch__496__;
  acegen_scratch__59__ = acegen_scratch__12__ * acegen_scratch__42__ +
                         acegen_scratch__13__ * acegen_scratch__56__ +
                         acegen_scratch__14__ * acegen_scratch__57__;
  acegen_scratch__60__ = acegen_scratch__13__ * acegen_scratch__46__ +
                         acegen_scratch__12__ * acegen_scratch__56__ +
                         acegen_scratch__14__ * acegen_scratch__58__;
  acegen_scratch__61__ = acegen_scratch__14__ * acegen_scratch__48__ +
                         acegen_scratch__12__ * acegen_scratch__57__ +
                         acegen_scratch__13__ * acegen_scratch__58__;
  acegen_scratch__62__ = acegen_scratch__15__ * acegen_scratch__42__ +
                         acegen_scratch__16__ * acegen_scratch__56__ +
                         acegen_scratch__17__ * acegen_scratch__57__;
  acegen_scratch__63__ = acegen_scratch__16__ * acegen_scratch__46__ +
                         acegen_scratch__15__ * acegen_scratch__56__ +
                         acegen_scratch__17__ * acegen_scratch__58__;
  acegen_scratch__64__ = acegen_scratch__17__ * acegen_scratch__48__ +
                         acegen_scratch__15__ * acegen_scratch__57__ +
                         acegen_scratch__16__ * acegen_scratch__58__;
  acegen_scratch__65__ = acegen_scratch__18__ * acegen_scratch__42__ +
                         acegen_scratch__19__ * acegen_scratch__56__ +
                         acegen_scratch__20__ * acegen_scratch__57__;
  acegen_scratch__66__ = acegen_scratch__19__ * acegen_scratch__46__ +
                         acegen_scratch__18__ * acegen_scratch__56__ +
                         acegen_scratch__20__ * acegen_scratch__58__;
  acegen_scratch__67__ = acegen_scratch__20__ * acegen_scratch__48__ +
                         acegen_scratch__18__ * acegen_scratch__57__ +
                         acegen_scratch__19__ * acegen_scratch__58__;
  acegen_scratch__143__ =
    2e0 * acegen_scratch__497__ * (acegen_scratch__74__ * acegen_scratch__74__);
  acegen_scratch__592__ = acegen_scratch__143__ * acegen_scratch__193__;
  acegen_scratch__144__ = 4e0 * (acegen_scratch__29__ * acegen_scratch__44__ +
                                 acegen_scratch__113__ * acegen_scratch__74__);
  acegen_scratch__145__ = 4e0 * (acegen_scratch__28__ * acegen_scratch__44__ +
                                 acegen_scratch__114__ * acegen_scratch__74__);
  acegen_scratch__146__ = acegen_scratch__115__ * acegen_scratch__502__ +
                          acegen_scratch__30__ * acegen_scratch__508__;
  acegen_scratch__147__ = acegen_scratch__116__ * acegen_scratch__502__;
  acegen_scratch__148__ = acegen_scratch__117__ * acegen_scratch__502__;
  acegen_scratch__149__ = acegen_scratch__113__ * acegen_scratch__506__;
  acegen_scratch__593__ = acegen_scratch__149__ * acegen_scratch__194__;
  acegen_scratch__150__ = 4e0 * (acegen_scratch__27__ * acegen_scratch__44__ +
                                 acegen_scratch__114__ * acegen_scratch__75__);
  acegen_scratch__151__ = acegen_scratch__115__ * acegen_scratch__506__;
  acegen_scratch__152__ = acegen_scratch__116__ * acegen_scratch__506__ +
                          acegen_scratch__31__ * acegen_scratch__508__;
  acegen_scratch__153__ = acegen_scratch__117__ * acegen_scratch__506__;
  acegen_scratch__154__ = acegen_scratch__114__ * acegen_scratch__511__;
  acegen_scratch__155__ = acegen_scratch__115__ * acegen_scratch__511__;
  acegen_scratch__156__ = acegen_scratch__116__ * acegen_scratch__511__;
  acegen_scratch__157__ = acegen_scratch__32__ * acegen_scratch__508__ +
                          acegen_scratch__117__ * acegen_scratch__511__;
  acegen_scratch__158__ =
    2e0 *
    (acegen_scratch__137__ * acegen_scratch__27__ +
     acegen_scratch__133__ *
       (acegen_scratch__103__ * acegen_scratch__516__ +
        acegen_scratch__41__ * acegen_scratch__53__ * acegen_scratch__84__));
  acegen_scratch__159__ = 2e0 * (acegen_scratch__110__ * acegen_scratch__133__ +
                                 acegen_scratch__32__ * acegen_scratch__53__);
  acegen_scratch__549__ = acegen_scratch__159__ * acegen_scratch__20__;
  acegen_scratch__310__ = acegen_scratch__159__ * acegen_scratch__521__;
  acegen_scratch__160__ = 2e0 * (acegen_scratch__111__ * acegen_scratch__133__ +
                                 acegen_scratch__31__ * acegen_scratch__53__);
  acegen_scratch__548__ = acegen_scratch__12__ * acegen_scratch__160__;
  acegen_scratch__528__ = acegen_scratch__15__ * acegen_scratch__160__;
  acegen_scratch__527__ = acegen_scratch__160__ * acegen_scratch__17__;
  acegen_scratch__161__ = 2e0 * (acegen_scratch__110__ * acegen_scratch__138__ +
                                 acegen_scratch__137__ * acegen_scratch__28__);
  acegen_scratch__530__ = acegen_scratch__161__ / 2e0;
  acegen_scratch__162__ = 2e0 * acegen_scratch__111__ * acegen_scratch__138__ +
                          acegen_scratch__137__ * acegen_scratch__532__;
  acegen_scratch__537__ = acegen_scratch__16__ * acegen_scratch__162__;
  acegen_scratch__308__ = acegen_scratch__15__ * acegen_scratch__537__;
  acegen_scratch__163__ = 2e0 * (acegen_scratch__111__ * acegen_scratch__141__ +
                                 acegen_scratch__137__ * acegen_scratch__29__);
  acegen_scratch__171__ = Power(acegen_scratch__12__, 3);
  acegen_scratch__172__ = 0.282842712474619e1 * acegen_scratch__13__;
  acegen_scratch__174__ = Power(acegen_scratch__13__, 3);
  acegen_scratch__175__ = 0.282842712474619e1 * acegen_scratch__12__;
  acegen_scratch__176__ = Power(acegen_scratch__14__, 3);
  acegen_scratch__599__ =
    (Power(acegen_scratch__12__, 4) * acegen_scratch__143__ +
     Power(acegen_scratch__13__, 4) * acegen_scratch__149__ +
     Power(acegen_scratch__14__, 4) * acegen_scratch__154__ +
     0.282842712474619e1 * acegen_scratch__14__ *
       (acegen_scratch__12__ * acegen_scratch__152__ * acegen_scratch__165__ +
        acegen_scratch__147__ * acegen_scratch__171__ +
        acegen_scratch__151__ * acegen_scratch__174__) +
     acegen_scratch__175__ * (acegen_scratch__153__ * acegen_scratch__174__ +
                              acegen_scratch__156__ * acegen_scratch__176__) +
     4e0 * acegen_scratch__13__ * acegen_scratch__159__ *
       acegen_scratch__524__ +
     acegen_scratch__172__ * (acegen_scratch__148__ * acegen_scratch__171__ +
                              acegen_scratch__155__ * acegen_scratch__176__ +
                              acegen_scratch__157__ * acegen_scratch__524__) +
     (0.282842712474619e1 * acegen_scratch__146__ +
      4e0 * acegen_scratch__162__) *
       acegen_scratch__164__ * acegen_scratch__533__ +
     (acegen_scratch__145__ + acegen_scratch__161__) * acegen_scratch__167__ *
       acegen_scratch__541__ +
     acegen_scratch__165__ *
       (2e0 * (acegen_scratch__150__ + acegen_scratch__158__) *
          acegen_scratch__167__ +
        (acegen_scratch__144__ + acegen_scratch__163__) *
          acegen_scratch__541__ +
        4e0 * acegen_scratch__14__ * acegen_scratch__548__)) /
    acegen_scratch__488__;
  acegen_scratch__184__ = 0.1414213562373095e1 * acegen_scratch__547__;
  acegen_scratch__196__ = 0.1414213562373095e1 * acegen_scratch__535__;
  acegen_scratch__197__ = 0.7071067811865476e0 * acegen_scratch__283__;
  acegen_scratch__198__ = 0.1414213562373095e1 * acegen_scratch__520__;
  acegen_scratch__273__ = 0.282842712474619e1 * acegen_scratch__16__;
  acegen_scratch__276__ = 0.282842712474619e1 * acegen_scratch__15__;
  acegen_scratch__581__ = acegen_scratch__16__ * acegen_scratch__276__;
  acegen_scratch__279__ = acegen_scratch__193__ * acegen_scratch__572__;
  acegen_scratch__280__ = acegen_scratch__194__ * acegen_scratch__572__;
  acegen_scratch__282__ = 2e0 * acegen_scratch__573__;
  acegen_scratch__600__ =
    (acegen_scratch__184__ * (acegen_scratch__148__ * acegen_scratch__193__ +
                              acegen_scratch__153__ * acegen_scratch__194__) +
     acegen_scratch__180__ * (acegen_scratch__145__ * acegen_scratch__193__ +
                              acegen_scratch__150__ * acegen_scratch__194__ +
                              acegen_scratch__154__ * acegen_scratch__195__ +
                              acegen_scratch__155__ * acegen_scratch__196__ +
                              acegen_scratch__156__ * acegen_scratch__197__ +
                              acegen_scratch__157__ * acegen_scratch__198__) +
     (acegen_scratch__146__ * acegen_scratch__273__ +
      acegen_scratch__147__ * acegen_scratch__276__) *
       acegen_scratch__279__ +
     (acegen_scratch__151__ * acegen_scratch__273__ +
      acegen_scratch__152__ * acegen_scratch__276__) *
       acegen_scratch__280__ +
     acegen_scratch__282__ * (acegen_scratch__163__ * acegen_scratch__498__ +
                              acegen_scratch__158__ * acegen_scratch__512__ +
                              acegen_scratch__18__ * acegen_scratch__527__ +
                              acegen_scratch__20__ * acegen_scratch__528__) +
     acegen_scratch__283__ * (acegen_scratch__308__ + acegen_scratch__310__ +
                              acegen_scratch__284__ * acegen_scratch__530__) +
     acegen_scratch__19__ *
       (0.1414213562373095e1 * acegen_scratch__162__ * acegen_scratch__390__ +
        acegen_scratch__284__ * acegen_scratch__549__) +
     acegen_scratch__195__ *
       (acegen_scratch__145__ * acegen_scratch__178__ +
        acegen_scratch__150__ * acegen_scratch__179__ +
        (acegen_scratch__155__ * acegen_scratch__273__ +
         acegen_scratch__156__ * acegen_scratch__276__) *
          acegen_scratch__572__ +
        (acegen_scratch__157__ * acegen_scratch__581__) / 2e0) +
     acegen_scratch__178__ *
       (acegen_scratch__144__ * acegen_scratch__194__ +
        acegen_scratch__146__ * acegen_scratch__196__ +
        acegen_scratch__147__ * acegen_scratch__197__ +
        acegen_scratch__148__ * acegen_scratch__198__ + acegen_scratch__592__) +
     acegen_scratch__179__ * (acegen_scratch__144__ * acegen_scratch__193__ +
                              acegen_scratch__151__ * acegen_scratch__196__ +
                              acegen_scratch__152__ * acegen_scratch__197__ +
                              acegen_scratch__153__ * acegen_scratch__198__ +
                              acegen_scratch__593__)) /
    acegen_scratch__488__;
  gradientOut[0][0] = acegen_scratch__59__;
  gradientOut[0][1] = acegen_scratch__60__;
  gradientOut[0][2] = acegen_scratch__61__;
  gradientOut[1][0] = acegen_scratch__62__;
  gradientOut[1][1] = acegen_scratch__63__;
  gradientOut[1][2] = acegen_scratch__64__;
  gradientOut[2][0] = acegen_scratch__65__;
  gradientOut[2][1] = acegen_scratch__66__;
  gradientOut[2][2] = acegen_scratch__67__;
  cacheOut[0]       = acegen_scratch__599__;
  cacheOut[1]       = acegen_scratch__600__;
  cacheOut[2]       = acegen_scratch__600__;
  cacheOut[3]       = 0e0;
  cacheOut[4]       = 0e0;
  cacheOut[5]       = 0e0;
  cacheOut[6]       = acegen_scratch__599__;
  cacheOut[7]       = acegen_scratch__600__;
  cacheOut[8]       = 0e0;
  cacheOut[9]       = 0e0;
  cacheOut[10]      = 0e0;
  cacheOut[11]      = acegen_scratch__599__;
  cacheOut[12]      = 0e0;
  cacheOut[13]      = 0e0;
  cacheOut[14]      = 0e0;
  cacheOut[15]      = acegen_scratch__601__;
  cacheOut[16]      = 0e0;
  cacheOut[17]      = 0e0;
  cacheOut[18]      = acegen_scratch__601__;
  cacheOut[19]      = 0e0;
  cacheOut[20]      = acegen_scratch__601__;
  cacheOut[21]      = (acegen_scratch__12__ * acegen_scratch__59__ +
                  acegen_scratch__13__ * acegen_scratch__60__ +
                  acegen_scratch__14__ * acegen_scratch__61__) /
                 acegen_scratch__488__;
  cacheOut[22] = (acegen_scratch__15__ * acegen_scratch__59__ +
                  acegen_scratch__16__ * acegen_scratch__60__ +
                  acegen_scratch__17__ * acegen_scratch__61__) /
                 acegen_scratch__488__;
  cacheOut[23] = (acegen_scratch__18__ * acegen_scratch__59__ +
                  acegen_scratch__19__ * acegen_scratch__60__ +
                  acegen_scratch__20__ * acegen_scratch__61__) /
                 acegen_scratch__488__;
  cacheOut[24] = (acegen_scratch__15__ * acegen_scratch__62__ +
                  acegen_scratch__16__ * acegen_scratch__63__ +
                  acegen_scratch__17__ * acegen_scratch__64__) /
                 acegen_scratch__488__;
  cacheOut[25] = (acegen_scratch__18__ * acegen_scratch__62__ +
                  acegen_scratch__19__ * acegen_scratch__63__ +
                  acegen_scratch__20__ * acegen_scratch__64__) /
                 acegen_scratch__488__;
  cacheOut[26] = (acegen_scratch__18__ * acegen_scratch__65__ +
                  acegen_scratch__19__ * acegen_scratch__66__ +
                  acegen_scratch__20__ * acegen_scratch__67__) /
                 acegen_scratch__488__;
}
