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
    acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__19__,
    acegen_scratch__24__, acegen_scratch__33__, acegen_scratch__36__,
    acegen_scratch__37__, acegen_scratch__39__, acegen_scratch__47__,
    acegen_scratch__48__, acegen_scratch__49__, acegen_scratch__50__,
    acegen_scratch__7__, acegen_scratch__9__;
  acegen_scratch__7__  = mu;
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
  acegen_scratch__24__ = sqrt(acegen_scratch__16__ * acegen_scratch__17__ -
                              (acegen_scratch__19__ * acegen_scratch__19__));
  acegen_scratch__47__ = 1e0 / (acegen_scratch__24__ * acegen_scratch__24__);
  acegen_scratch__48__ = acegen_scratch__47__ / 2e0;
  acegen_scratch__50__ = 2e0 * acegen_scratch__16__ * acegen_scratch__48__;
  acegen_scratch__49__ = 2e0 * acegen_scratch__17__ * acegen_scratch__48__;
  acegen_scratch__33__ = 2e0 * lambda * log(acegen_scratch__24__);
  acegen_scratch__36__ = acegen_scratch__33__ * acegen_scratch__49__ +
                         (1e0 - acegen_scratch__49__) * acegen_scratch__7__;
  acegen_scratch__37__ = acegen_scratch__33__ * acegen_scratch__50__ +
                         (1e0 - acegen_scratch__50__) * acegen_scratch__7__;
  acegen_scratch__39__ = -1e0 * acegen_scratch__19__ * acegen_scratch__47__ *
                         (acegen_scratch__33__ - acegen_scratch__7__);
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  gradientOut[0][0] = acegen_scratch__10__ * acegen_scratch__39__ +
                      acegen_scratch__36__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__10__ * acegen_scratch__37__ +
                      acegen_scratch__39__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__11__ * acegen_scratch__36__ +
                      acegen_scratch__12__ * acegen_scratch__39__;
  gradientOut[1][1] = acegen_scratch__12__ * acegen_scratch__37__ +
                      acegen_scratch__11__ * acegen_scratch__39__;
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
  Number acegen_scratch__10__, acegen_scratch__11__, acegen_scratch__12__,
    acegen_scratch__13__, acegen_scratch__139__, acegen_scratch__140__,
    acegen_scratch__141__, acegen_scratch__144__, acegen_scratch__146__,
    acegen_scratch__147__, acegen_scratch__15__, acegen_scratch__155__,
    acegen_scratch__159__, acegen_scratch__16__, acegen_scratch__160__,
    acegen_scratch__162__, acegen_scratch__17__, acegen_scratch__171__,
    acegen_scratch__172__, acegen_scratch__173__, acegen_scratch__174__,
    acegen_scratch__177__, acegen_scratch__178__, acegen_scratch__179__,
    acegen_scratch__18__, acegen_scratch__180__, acegen_scratch__181__,
    acegen_scratch__182__, acegen_scratch__183__, acegen_scratch__22__,
    acegen_scratch__23__, acegen_scratch__25__, acegen_scratch__34__,
    acegen_scratch__35__, acegen_scratch__36__, acegen_scratch__37__,
    acegen_scratch__39__, acegen_scratch__42__, acegen_scratch__43__,
    acegen_scratch__45__, acegen_scratch__9__;
  acegen_scratch__9__   = gradduIn[0][0];
  acegen_scratch__10__  = gradduIn[0][1];
  acegen_scratch__11__  = gradduIn[1][0];
  acegen_scratch__12__  = gradduIn[1][1];
  acegen_scratch__13__  = mu;
  acegen_scratch__177__ = 2e0 * lambda;
  acegen_scratch__15__  = 1e0 + graduIn[0][0];
  acegen_scratch__16__  = graduIn[0][1];
  acegen_scratch__17__  = graduIn[1][0];
  acegen_scratch__139__ = 2e0 * (acegen_scratch__11__ * acegen_scratch__17__ +
                                 acegen_scratch__15__ * acegen_scratch__9__);
  acegen_scratch__18__  = 1e0 + graduIn[1][1];
  acegen_scratch__140__ = acegen_scratch__10__ * acegen_scratch__15__ +
                          acegen_scratch__12__ * acegen_scratch__17__ +
                          acegen_scratch__11__ * acegen_scratch__18__ +
                          acegen_scratch__16__ * acegen_scratch__9__;
  acegen_scratch__141__ = 2e0 * (acegen_scratch__10__ * acegen_scratch__16__ +
                                 acegen_scratch__12__ * acegen_scratch__18__);
  acegen_scratch__22__  = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__17__ * acegen_scratch__17__);
  acegen_scratch__23__ = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__25__ = acegen_scratch__15__ * acegen_scratch__16__ +
                         acegen_scratch__17__ * acegen_scratch__18__;
  acegen_scratch__144__ = acegen_scratch__22__ * acegen_scratch__23__ -
                          (acegen_scratch__25__ * acegen_scratch__25__);
  acegen_scratch__171__ = sqrt(acegen_scratch__144__);
  acegen_scratch__173__ = 1e0 / (2e0 * acegen_scratch__171__);
  acegen_scratch__181__ = acegen_scratch__173__ / acegen_scratch__171__;
  acegen_scratch__146__ = acegen_scratch__173__ *
                          (acegen_scratch__141__ * acegen_scratch__22__ +
                           acegen_scratch__139__ * acegen_scratch__23__ -
                           2e0 * acegen_scratch__140__ * acegen_scratch__25__);
  acegen_scratch__174__ = 1e0 / (acegen_scratch__171__ * acegen_scratch__171__);
  acegen_scratch__172__ = -0.7071067811865476e0 / acegen_scratch__171__;
  acegen_scratch__147__ = -(acegen_scratch__146__ * acegen_scratch__174__);
  acegen_scratch__183__ = 2e0 * acegen_scratch__147__ * acegen_scratch__171__;
  acegen_scratch__34__  = acegen_scratch__172__ * acegen_scratch__25__;
  acegen_scratch__155__ =
    (acegen_scratch__146__ * acegen_scratch__177__) / acegen_scratch__171__;
  acegen_scratch__37__  = acegen_scratch__34__ / acegen_scratch__171__;
  acegen_scratch__36__  = acegen_scratch__181__ * acegen_scratch__22__;
  acegen_scratch__180__ = 2e0 * acegen_scratch__36__;
  acegen_scratch__35__  = acegen_scratch__181__ * acegen_scratch__23__;
  acegen_scratch__179__ = 2e0 * acegen_scratch__35__;
  acegen_scratch__39__  = acegen_scratch__177__ * log(acegen_scratch__171__);
  acegen_scratch__178__ = -acegen_scratch__13__ + acegen_scratch__39__;
  acegen_scratch__182__ = (acegen_scratch__174__ * acegen_scratch__178__) / 2e0;
  acegen_scratch__159__ =
    2e0 *
    (acegen_scratch__182__ *
       (acegen_scratch__141__ + acegen_scratch__183__ * acegen_scratch__23__) +
     acegen_scratch__155__ * acegen_scratch__35__);
  acegen_scratch__160__ =
    2e0 *
    (acegen_scratch__182__ *
       (acegen_scratch__139__ + acegen_scratch__183__ * acegen_scratch__22__) +
     acegen_scratch__155__ * acegen_scratch__36__);
  acegen_scratch__162__ =
    0.1414213562373095e1 *
    (acegen_scratch__178__ *
       ((acegen_scratch__140__ * acegen_scratch__172__ -
         0.7071067811865476e0 * acegen_scratch__147__ * acegen_scratch__25__) /
          acegen_scratch__171__ -
        (acegen_scratch__146__ * acegen_scratch__34__) /
          acegen_scratch__144__) +
     acegen_scratch__155__ * acegen_scratch__37__);
  acegen_scratch__42__ = acegen_scratch__13__ * (1e0 - acegen_scratch__179__) +
                         acegen_scratch__179__ * acegen_scratch__39__;
  acegen_scratch__43__ = acegen_scratch__13__ * (1e0 - acegen_scratch__180__) +
                         acegen_scratch__180__ * acegen_scratch__39__;
  acegen_scratch__45__ =
    0.1414213562373095e1 * acegen_scratch__178__ * acegen_scratch__37__;
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  gradientOut[0][0] = acegen_scratch__15__ * acegen_scratch__159__ +
                      acegen_scratch__16__ * acegen_scratch__162__ +
                      acegen_scratch__10__ * acegen_scratch__45__ +
                      acegen_scratch__42__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__160__ +
                      acegen_scratch__15__ * acegen_scratch__162__ +
                      acegen_scratch__10__ * acegen_scratch__43__ +
                      acegen_scratch__45__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__159__ * acegen_scratch__17__ +
                      acegen_scratch__162__ * acegen_scratch__18__ +
                      acegen_scratch__11__ * acegen_scratch__42__ +
                      acegen_scratch__12__ * acegen_scratch__45__;
  gradientOut[1][1] = acegen_scratch__162__ * acegen_scratch__17__ +
                      acegen_scratch__160__ * acegen_scratch__18__ +
                      acegen_scratch__12__ * acegen_scratch__43__ +
                      acegen_scratch__11__ * acegen_scratch__45__;
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
  valueOut[1]       = 0e0;
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
    acegen_scratch__23__, acegen_scratch__30__, acegen_scratch__31__,
    acegen_scratch__32__, acegen_scratch__36__, acegen_scratch__37__,
    acegen_scratch__38__, acegen_scratch__44__, acegen_scratch__45__,
    acegen_scratch__46__, acegen_scratch__47__, acegen_scratch__48__,
    acegen_scratch__49__, acegen_scratch__50__, acegen_scratch__51__,
    acegen_scratch__53__, acegen_scratch__56__, acegen_scratch__59__,
    acegen_scratch__60__, acegen_scratch__61__, acegen_scratch__76__,
    acegen_scratch__77__, acegen_scratch__78__;
  acegen_scratch__13__ = mu;
  acegen_scratch__46__ = acegen_scratch__13__ / 2e0;
  acegen_scratch__15__ = 1e0 + graduIn[0][0];
  acegen_scratch__16__ = graduIn[0][1];
  acegen_scratch__17__ = graduIn[0][2];
  acegen_scratch__18__ = graduIn[1][0];
  acegen_scratch__19__ = 1e0 + graduIn[1][1];
  acegen_scratch__20__ = graduIn[1][2];
  acegen_scratch__21__ = graduIn[2][0];
  acegen_scratch__22__ = graduIn[2][1];
  acegen_scratch__23__ = 1e0 + graduIn[2][2];
  acegen_scratch__30__ = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__18__ * acegen_scratch__18__) +
                         (acegen_scratch__21__ * acegen_scratch__21__);
  acegen_scratch__31__ = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__19__ * acegen_scratch__19__) +
                         (acegen_scratch__22__ * acegen_scratch__22__);
  acegen_scratch__32__ = (acegen_scratch__17__ * acegen_scratch__17__) +
                         (acegen_scratch__20__ * acegen_scratch__20__) +
                         (acegen_scratch__23__ * acegen_scratch__23__);
  acegen_scratch__36__ = acegen_scratch__15__ * acegen_scratch__16__ +
                         acegen_scratch__18__ * acegen_scratch__19__ +
                         acegen_scratch__21__ * acegen_scratch__22__;
  acegen_scratch__50__ = (acegen_scratch__36__ * acegen_scratch__36__);
  acegen_scratch__37__ = acegen_scratch__15__ * acegen_scratch__17__ +
                         acegen_scratch__18__ * acegen_scratch__20__ +
                         acegen_scratch__21__ * acegen_scratch__23__;
  acegen_scratch__53__ = 2e0 * acegen_scratch__36__ * acegen_scratch__37__;
  acegen_scratch__48__ = (acegen_scratch__37__ * acegen_scratch__37__);
  acegen_scratch__38__ = acegen_scratch__16__ * acegen_scratch__17__ +
                         acegen_scratch__19__ * acegen_scratch__20__ +
                         acegen_scratch__22__ * acegen_scratch__23__;
  acegen_scratch__78__ = -1e0 * acegen_scratch__38__;
  acegen_scratch__76__ = acegen_scratch__31__ * acegen_scratch__32__ -
                         (acegen_scratch__38__ * acegen_scratch__38__);
  acegen_scratch__44__ = -(acegen_scratch__31__ * acegen_scratch__48__) -
                         acegen_scratch__32__ * acegen_scratch__50__ +
                         acegen_scratch__38__ * acegen_scratch__53__ +
                         acegen_scratch__30__ * acegen_scratch__76__;
  acegen_scratch__56__ =
    (0.7071067811865476e0 *
     (-acegen_scratch__13__ + 2e0 * lambda * log(sqrt(acegen_scratch__44__)))) /
    acegen_scratch__44__;
  acegen_scratch__77__ = -0.1414213562373095e1 * acegen_scratch__56__;
  acegen_scratch__47__ = 0.7071067811865476e0 * acegen_scratch__56__;
  acegen_scratch__45__ =
    2e0 * (acegen_scratch__46__ + acegen_scratch__47__ * acegen_scratch__76__);
  acegen_scratch__49__ =
    2e0 * (acegen_scratch__46__ +
           acegen_scratch__47__ * (acegen_scratch__30__ * acegen_scratch__32__ -
                                   acegen_scratch__48__));
  acegen_scratch__51__ =
    2e0 * (acegen_scratch__46__ +
           acegen_scratch__47__ * (acegen_scratch__30__ * acegen_scratch__31__ -
                                   acegen_scratch__50__));
  acegen_scratch__59__ =
    acegen_scratch__77__ * (acegen_scratch__32__ * acegen_scratch__36__ +
                            acegen_scratch__37__ * acegen_scratch__78__);
  acegen_scratch__60__ =
    acegen_scratch__77__ * (acegen_scratch__31__ * acegen_scratch__37__ +
                            acegen_scratch__36__ * acegen_scratch__78__);
  acegen_scratch__61__ =
    1e0 * acegen_scratch__47__ *
    (-2e0 * acegen_scratch__30__ * acegen_scratch__38__ + acegen_scratch__53__);
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  valueOut[2]       = 0e0;
  gradientOut[0][0] = acegen_scratch__15__ * acegen_scratch__45__ +
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
  Number acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__19__, acegen_scratch__196__, acegen_scratch__20__,
    acegen_scratch__206__, acegen_scratch__21__, acegen_scratch__216__,
    acegen_scratch__22__, acegen_scratch__226__, acegen_scratch__23__,
    acegen_scratch__24__, acegen_scratch__243__, acegen_scratch__25__,
    acegen_scratch__26__, acegen_scratch__262__, acegen_scratch__27__,
    acegen_scratch__28__, acegen_scratch__29__, acegen_scratch__30__,
    acegen_scratch__31__, acegen_scratch__32__, acegen_scratch__33__,
    acegen_scratch__34__, acegen_scratch__35__, acegen_scratch__356__,
    acegen_scratch__357__, acegen_scratch__358__, acegen_scratch__359__,
    acegen_scratch__360__, acegen_scratch__361__, acegen_scratch__367__,
    acegen_scratch__368__, acegen_scratch__374__, acegen_scratch__376__,
    acegen_scratch__381__, acegen_scratch__382__, acegen_scratch__383__,
    acegen_scratch__384__, acegen_scratch__385__, acegen_scratch__389__,
    acegen_scratch__390__, acegen_scratch__391__, acegen_scratch__405__,
    acegen_scratch__407__, acegen_scratch__408__, acegen_scratch__409__,
    acegen_scratch__410__, acegen_scratch__411__, acegen_scratch__412__,
    acegen_scratch__42__, acegen_scratch__43__, acegen_scratch__44__,
    acegen_scratch__48__, acegen_scratch__49__, acegen_scratch__50__,
    acegen_scratch__51__, acegen_scratch__54__, acegen_scratch__55__,
    acegen_scratch__56__, acegen_scratch__57__, acegen_scratch__58__,
    acegen_scratch__59__, acegen_scratch__60__, acegen_scratch__61__,
    acegen_scratch__62__, acegen_scratch__63__, acegen_scratch__65__,
    acegen_scratch__68__, acegen_scratch__71__, acegen_scratch__72__,
    acegen_scratch__73__;
  acegen_scratch__16__  = gradduIn[0][0];
  acegen_scratch__17__  = gradduIn[0][1];
  acegen_scratch__18__  = gradduIn[0][2];
  acegen_scratch__19__  = gradduIn[1][0];
  acegen_scratch__20__  = gradduIn[1][1];
  acegen_scratch__21__  = gradduIn[1][2];
  acegen_scratch__22__  = gradduIn[2][0];
  acegen_scratch__23__  = gradduIn[2][1];
  acegen_scratch__24__  = gradduIn[2][2];
  acegen_scratch__25__  = mu;
  acegen_scratch__58__  = acegen_scratch__25__ / 2e0;
  acegen_scratch__26__  = lambda;
  acegen_scratch__27__  = 1e0 + graduIn[0][0];
  acegen_scratch__28__  = graduIn[0][1];
  acegen_scratch__29__  = graduIn[0][2];
  acegen_scratch__30__  = graduIn[1][0];
  acegen_scratch__31__  = 1e0 + graduIn[1][1];
  acegen_scratch__32__  = graduIn[1][2];
  acegen_scratch__33__  = graduIn[2][0];
  acegen_scratch__356__ = 2e0 * (acegen_scratch__16__ * acegen_scratch__27__ +
                                 acegen_scratch__19__ * acegen_scratch__30__ +
                                 acegen_scratch__22__ * acegen_scratch__33__);
  acegen_scratch__34__  = graduIn[2][1];
  acegen_scratch__357__ = acegen_scratch__17__ * acegen_scratch__27__ +
                          acegen_scratch__16__ * acegen_scratch__28__ +
                          acegen_scratch__20__ * acegen_scratch__30__ +
                          acegen_scratch__19__ * acegen_scratch__31__ +
                          acegen_scratch__23__ * acegen_scratch__33__ +
                          acegen_scratch__22__ * acegen_scratch__34__;
  acegen_scratch__359__ = 2e0 * (acegen_scratch__17__ * acegen_scratch__28__ +
                                 acegen_scratch__20__ * acegen_scratch__31__ +
                                 acegen_scratch__23__ * acegen_scratch__34__);
  acegen_scratch__35__  = 1e0 + graduIn[2][2];
  acegen_scratch__360__ = acegen_scratch__18__ * acegen_scratch__28__ +
                          acegen_scratch__17__ * acegen_scratch__29__ +
                          acegen_scratch__21__ * acegen_scratch__31__ +
                          acegen_scratch__20__ * acegen_scratch__32__ +
                          acegen_scratch__24__ * acegen_scratch__34__ +
                          acegen_scratch__23__ * acegen_scratch__35__;
  acegen_scratch__358__ = acegen_scratch__18__ * acegen_scratch__27__ +
                          acegen_scratch__16__ * acegen_scratch__29__ +
                          acegen_scratch__21__ * acegen_scratch__30__ +
                          acegen_scratch__19__ * acegen_scratch__32__ +
                          acegen_scratch__24__ * acegen_scratch__33__ +
                          acegen_scratch__22__ * acegen_scratch__35__;
  acegen_scratch__361__ = 2e0 * (acegen_scratch__18__ * acegen_scratch__29__ +
                                 acegen_scratch__21__ * acegen_scratch__32__ +
                                 acegen_scratch__24__ * acegen_scratch__35__);
  acegen_scratch__42__  = (acegen_scratch__27__ * acegen_scratch__27__) +
                         (acegen_scratch__30__ * acegen_scratch__30__) +
                         (acegen_scratch__33__ * acegen_scratch__33__);
  acegen_scratch__43__ = (acegen_scratch__28__ * acegen_scratch__28__) +
                         (acegen_scratch__31__ * acegen_scratch__31__) +
                         (acegen_scratch__34__ * acegen_scratch__34__);
  acegen_scratch__44__ = (acegen_scratch__29__ * acegen_scratch__29__) +
                         (acegen_scratch__32__ * acegen_scratch__32__) +
                         (acegen_scratch__35__ * acegen_scratch__35__);
  acegen_scratch__48__ = acegen_scratch__27__ * acegen_scratch__28__ +
                         acegen_scratch__30__ * acegen_scratch__31__ +
                         acegen_scratch__33__ * acegen_scratch__34__;
  acegen_scratch__407__ = 2e0 * acegen_scratch__48__;
  acegen_scratch__412__ = -(acegen_scratch__357__ * acegen_scratch__407__) +
                          acegen_scratch__359__ * acegen_scratch__42__ +
                          acegen_scratch__356__ * acegen_scratch__43__;
  acegen_scratch__62__ = (acegen_scratch__48__ * acegen_scratch__48__);
  acegen_scratch__216__ =
    acegen_scratch__42__ * acegen_scratch__43__ - acegen_scratch__62__;
  acegen_scratch__49__ = acegen_scratch__27__ * acegen_scratch__29__ +
                         acegen_scratch__30__ * acegen_scratch__32__ +
                         acegen_scratch__33__ * acegen_scratch__35__;
  acegen_scratch__405__ = 2e0 * acegen_scratch__49__;
  acegen_scratch__368__ = acegen_scratch__358__ * acegen_scratch__405__;
  acegen_scratch__367__ = 2e0 * (acegen_scratch__358__ * acegen_scratch__48__ +
                                 acegen_scratch__357__ * acegen_scratch__49__);
  acegen_scratch__65__  = acegen_scratch__405__ * acegen_scratch__48__;
  acegen_scratch__60__  = (acegen_scratch__49__ * acegen_scratch__49__);
  acegen_scratch__206__ =
    acegen_scratch__42__ * acegen_scratch__44__ - acegen_scratch__60__;
  acegen_scratch__50__ = acegen_scratch__28__ * acegen_scratch__29__ +
                         acegen_scratch__31__ * acegen_scratch__32__ +
                         acegen_scratch__34__ * acegen_scratch__35__;
  acegen_scratch__408__ = 2e0 * acegen_scratch__50__;
  acegen_scratch__374__ = acegen_scratch__360__ * acegen_scratch__408__;
  acegen_scratch__262__ = -(acegen_scratch__407__ * acegen_scratch__44__) +
                          acegen_scratch__405__ * acegen_scratch__50__;
  acegen_scratch__243__ = -(acegen_scratch__405__ * acegen_scratch__43__) +
                          acegen_scratch__407__ * acegen_scratch__50__;
  acegen_scratch__226__ =
    -(acegen_scratch__408__ * acegen_scratch__42__) + acegen_scratch__65__;
  acegen_scratch__55__  = (acegen_scratch__50__ * acegen_scratch__50__);
  acegen_scratch__376__ = acegen_scratch__216__ * acegen_scratch__361__ -
                          acegen_scratch__374__ * acegen_scratch__42__ -
                          acegen_scratch__368__ * acegen_scratch__43__ +
                          acegen_scratch__412__ * acegen_scratch__44__ +
                          acegen_scratch__367__ * acegen_scratch__50__ -
                          acegen_scratch__356__ * acegen_scratch__55__ -
                          acegen_scratch__359__ * acegen_scratch__60__ +
                          acegen_scratch__360__ * acegen_scratch__65__;
  acegen_scratch__196__ =
    acegen_scratch__43__ * acegen_scratch__44__ - acegen_scratch__55__;
  acegen_scratch__56__ = acegen_scratch__196__ * acegen_scratch__42__ -
                         acegen_scratch__43__ * acegen_scratch__60__ -
                         acegen_scratch__44__ * acegen_scratch__62__ +
                         acegen_scratch__50__ * acegen_scratch__65__;
  acegen_scratch__51__ = sqrt(acegen_scratch__56__);
  acegen_scratch__410__ =
    (acegen_scratch__26__ * acegen_scratch__376__) /
    ((acegen_scratch__51__ * acegen_scratch__51__) * acegen_scratch__56__);
  acegen_scratch__54__ = -acegen_scratch__25__ +
                         2e0 * acegen_scratch__26__ * log(acegen_scratch__51__);
  acegen_scratch__409__ = -((acegen_scratch__376__ * acegen_scratch__54__) /
                            (acegen_scratch__56__ * acegen_scratch__56__));
  acegen_scratch__382__ = (acegen_scratch__409__ + acegen_scratch__410__) / 2e0;
  acegen_scratch__381__ =
    0.7071067811865476e0 * (acegen_scratch__409__ + acegen_scratch__410__);
  acegen_scratch__68__ =
    (0.7071067811865476e0 * acegen_scratch__54__) / acegen_scratch__56__;
  acegen_scratch__411__ = 2e0 * acegen_scratch__68__;
  acegen_scratch__389__ =
    0.7071067811865476e0 *
    (acegen_scratch__262__ * acegen_scratch__381__ +
     acegen_scratch__411__ * (-(acegen_scratch__357__ * acegen_scratch__44__) -
                              acegen_scratch__361__ * acegen_scratch__48__ +
                              acegen_scratch__360__ * acegen_scratch__49__ +
                              acegen_scratch__358__ * acegen_scratch__50__));
  acegen_scratch__390__ =
    0.7071067811865476e0 *
    (acegen_scratch__243__ * acegen_scratch__381__ +
     acegen_scratch__411__ * (-(acegen_scratch__358__ * acegen_scratch__43__) +
                              acegen_scratch__360__ * acegen_scratch__48__ -
                              acegen_scratch__359__ * acegen_scratch__49__ +
                              acegen_scratch__357__ * acegen_scratch__50__));
  acegen_scratch__391__ =
    0.7071067811865476e0 *
    (acegen_scratch__226__ * acegen_scratch__381__ +
     (acegen_scratch__367__ - acegen_scratch__356__ * acegen_scratch__408__ -
      2e0 * acegen_scratch__360__ * acegen_scratch__42__) *
       acegen_scratch__68__);
  acegen_scratch__59__  = 0.7071067811865476e0 * acegen_scratch__68__;
  acegen_scratch__385__ = 2e0 * (acegen_scratch__216__ * acegen_scratch__382__ +
                                 acegen_scratch__412__ * acegen_scratch__59__);
  acegen_scratch__384__ =
    2e0 *
    (acegen_scratch__206__ * acegen_scratch__382__ +
     (-acegen_scratch__368__ + acegen_scratch__361__ * acegen_scratch__42__ +
      acegen_scratch__356__ * acegen_scratch__44__) *
       acegen_scratch__59__);
  acegen_scratch__383__ =
    2e0 *
    (acegen_scratch__196__ * acegen_scratch__382__ +
     (-acegen_scratch__374__ + acegen_scratch__361__ * acegen_scratch__43__ +
      acegen_scratch__359__ * acegen_scratch__44__) *
       acegen_scratch__59__);
  acegen_scratch__57__ =
    2e0 * (acegen_scratch__58__ + acegen_scratch__196__ * acegen_scratch__59__);
  acegen_scratch__61__ =
    2e0 * (acegen_scratch__58__ + acegen_scratch__206__ * acegen_scratch__59__);
  acegen_scratch__63__ =
    2e0 * (acegen_scratch__58__ + acegen_scratch__216__ * acegen_scratch__59__);
  acegen_scratch__71__ = acegen_scratch__262__ * acegen_scratch__59__;
  acegen_scratch__72__ = acegen_scratch__243__ * acegen_scratch__59__;
  acegen_scratch__73__ = 1e0 * acegen_scratch__226__ * acegen_scratch__59__;
  valueOut[0]          = 0e0;
  valueOut[1]          = 0e0;
  valueOut[2]          = 0e0;
  gradientOut[0][0]    = acegen_scratch__27__ * acegen_scratch__383__ +
                      acegen_scratch__28__ * acegen_scratch__389__ +
                      acegen_scratch__29__ * acegen_scratch__390__ +
                      acegen_scratch__16__ * acegen_scratch__57__ +
                      acegen_scratch__17__ * acegen_scratch__71__ +
                      acegen_scratch__18__ * acegen_scratch__72__;
  gradientOut[0][1] = acegen_scratch__28__ * acegen_scratch__384__ +
                      acegen_scratch__27__ * acegen_scratch__389__ +
                      acegen_scratch__29__ * acegen_scratch__391__ +
                      acegen_scratch__17__ * acegen_scratch__61__ +
                      acegen_scratch__16__ * acegen_scratch__71__ +
                      acegen_scratch__18__ * acegen_scratch__73__;
  gradientOut[0][2] = acegen_scratch__29__ * acegen_scratch__385__ +
                      acegen_scratch__27__ * acegen_scratch__390__ +
                      acegen_scratch__28__ * acegen_scratch__391__ +
                      acegen_scratch__18__ * acegen_scratch__63__ +
                      acegen_scratch__16__ * acegen_scratch__72__ +
                      acegen_scratch__17__ * acegen_scratch__73__;
  gradientOut[1][0] = acegen_scratch__30__ * acegen_scratch__383__ +
                      acegen_scratch__31__ * acegen_scratch__389__ +
                      acegen_scratch__32__ * acegen_scratch__390__ +
                      acegen_scratch__19__ * acegen_scratch__57__ +
                      acegen_scratch__20__ * acegen_scratch__71__ +
                      acegen_scratch__21__ * acegen_scratch__72__;
  gradientOut[1][1] = acegen_scratch__31__ * acegen_scratch__384__ +
                      acegen_scratch__30__ * acegen_scratch__389__ +
                      acegen_scratch__32__ * acegen_scratch__391__ +
                      acegen_scratch__20__ * acegen_scratch__61__ +
                      acegen_scratch__19__ * acegen_scratch__71__ +
                      acegen_scratch__21__ * acegen_scratch__73__;
  gradientOut[1][2] = acegen_scratch__32__ * acegen_scratch__385__ +
                      acegen_scratch__30__ * acegen_scratch__390__ +
                      acegen_scratch__31__ * acegen_scratch__391__ +
                      acegen_scratch__21__ * acegen_scratch__63__ +
                      acegen_scratch__19__ * acegen_scratch__72__ +
                      acegen_scratch__20__ * acegen_scratch__73__;
  gradientOut[2][0] = acegen_scratch__33__ * acegen_scratch__383__ +
                      acegen_scratch__34__ * acegen_scratch__389__ +
                      acegen_scratch__35__ * acegen_scratch__390__ +
                      acegen_scratch__22__ * acegen_scratch__57__ +
                      acegen_scratch__23__ * acegen_scratch__71__ +
                      acegen_scratch__24__ * acegen_scratch__72__;
  gradientOut[2][1] = acegen_scratch__34__ * acegen_scratch__384__ +
                      acegen_scratch__33__ * acegen_scratch__389__ +
                      acegen_scratch__35__ * acegen_scratch__391__ +
                      acegen_scratch__23__ * acegen_scratch__61__ +
                      acegen_scratch__22__ * acegen_scratch__71__ +
                      acegen_scratch__24__ * acegen_scratch__73__;
  gradientOut[2][2] = acegen_scratch__35__ * acegen_scratch__385__ +
                      acegen_scratch__33__ * acegen_scratch__390__ +
                      acegen_scratch__34__ * acegen_scratch__391__ +
                      acegen_scratch__24__ * acegen_scratch__63__ +
                      acegen_scratch__22__ * acegen_scratch__72__ +
                      acegen_scratch__23__ * acegen_scratch__73__;
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
  valueOut[1]          = 0e0;
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
