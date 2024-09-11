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

#include "cbrt.h"

using namespace dealii;

template <typename Number, std::size_t width>
inline ::dealii::VectorizedArray<Number, width>
Power(const ::dealii::VectorizedArray<Number, width> &                    x,
      const typename ::dealii::VectorizedArray<Number, width>::value_type p)
{
  ::dealii::VectorizedArray<Number, width> out;
  for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size();
       ++i)
    out[i] = std::pow(x[i], p);
  return out;
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
    acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__20__, acegen_scratch__28__, acegen_scratch__32__,
    acegen_scratch__44__, acegen_scratch__45__, acegen_scratch__47__,
    acegen_scratch__54__, acegen_scratch__55__, acegen_scratch__56__,
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
  acegen_scratch__18__ =
    0.1414213562373095e1 * (acegen_scratch__11__ * acegen_scratch__12__ +
                            acegen_scratch__10__ * acegen_scratch__9__);
  acegen_scratch__20__ = acegen_scratch__16__ * acegen_scratch__17__ -
                         0.5e0 * (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__56__ =
    ((2e0 / 3e0) * acegen_scratch__7__ + lambda) / acegen_scratch__20__;
  acegen_scratch__55__ =
    ((-1e0 + acegen_scratch__20__) * acegen_scratch__56__) / 2e0;
  acegen_scratch__54__ =
    1e0 / (3e0 * Power(acegen_scratch__20__, 0.13333333333333333e1));
  acegen_scratch__28__ = -(acegen_scratch__17__ * acegen_scratch__54__);
  acegen_scratch__32__ =
    1e0 / Power(acegen_scratch__20__, 0.3333333333333333e0) +
    acegen_scratch__16__ * acegen_scratch__28__;
  acegen_scratch__44__ =
    acegen_scratch__17__ * acegen_scratch__55__ +
    (acegen_scratch__17__ * acegen_scratch__28__ + acegen_scratch__32__) *
      acegen_scratch__7__;
  acegen_scratch__45__ =
    acegen_scratch__16__ * acegen_scratch__55__ +
    (acegen_scratch__32__ -
     (acegen_scratch__16__ * acegen_scratch__16__) * acegen_scratch__54__) *
      acegen_scratch__7__;
  acegen_scratch__47__ =
    acegen_scratch__18__ *
    ((0.3535533905932738e0 - 0.3535533905932738e0 * acegen_scratch__20__) *
       acegen_scratch__56__ +
     0.7071067811865476e0 * (acegen_scratch__16__ + acegen_scratch__17__) *
       acegen_scratch__54__ * acegen_scratch__7__);
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  gradientOut[0][0] = acegen_scratch__10__ * acegen_scratch__47__ +
                      acegen_scratch__44__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__10__ * acegen_scratch__45__ +
                      acegen_scratch__47__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__11__ * acegen_scratch__44__ +
                      acegen_scratch__12__ * acegen_scratch__47__;
  gradientOut[1][1] = acegen_scratch__12__ * acegen_scratch__45__ +
                      acegen_scratch__11__ * acegen_scratch__47__;
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
    acegen_scratch__13__, acegen_scratch__15__, acegen_scratch__16__,
    acegen_scratch__17__, acegen_scratch__178__, acegen_scratch__179__,
    acegen_scratch__18__, acegen_scratch__180__, acegen_scratch__181__,
    acegen_scratch__182__, acegen_scratch__183__, acegen_scratch__185__,
    acegen_scratch__193__, acegen_scratch__195__, acegen_scratch__200__,
    acegen_scratch__201__, acegen_scratch__203__, acegen_scratch__210__,
    acegen_scratch__211__, acegen_scratch__212__, acegen_scratch__213__,
    acegen_scratch__214__, acegen_scratch__215__, acegen_scratch__216__,
    acegen_scratch__217__, acegen_scratch__218__, acegen_scratch__219__,
    acegen_scratch__22__, acegen_scratch__221__, acegen_scratch__222__,
    acegen_scratch__223__, acegen_scratch__23__, acegen_scratch__24__,
    acegen_scratch__26__, acegen_scratch__27__, acegen_scratch__35__,
    acegen_scratch__36__, acegen_scratch__37__, acegen_scratch__38__,
    acegen_scratch__39__, acegen_scratch__47__, acegen_scratch__51__,
    acegen_scratch__52__, acegen_scratch__54__, acegen_scratch__9__;
  acegen_scratch__9__   = gradduIn[0][0];
  acegen_scratch__10__  = gradduIn[0][1];
  acegen_scratch__11__  = gradduIn[1][0];
  acegen_scratch__12__  = gradduIn[1][1];
  acegen_scratch__13__  = (mu);
  acegen_scratch__47__  = (2e0 / 3e0) * acegen_scratch__13__ + (lambda);
  acegen_scratch__15__  = 1e0 + graduIn[0][0];
  acegen_scratch__16__  = graduIn[0][1];
  acegen_scratch__17__  = graduIn[1][0];
  acegen_scratch__178__ = 2e0 * (acegen_scratch__11__ * acegen_scratch__17__ +
                                 acegen_scratch__15__ * acegen_scratch__9__);
  acegen_scratch__18__  = 1e0 + graduIn[1][1];
  acegen_scratch__210__ = acegen_scratch__15__ * acegen_scratch__16__ +
                          acegen_scratch__17__ * acegen_scratch__18__;
  acegen_scratch__179__ = acegen_scratch__10__ * acegen_scratch__15__ +
                          acegen_scratch__12__ * acegen_scratch__17__ +
                          acegen_scratch__11__ * acegen_scratch__18__ +
                          acegen_scratch__16__ * acegen_scratch__9__;
  acegen_scratch__181__ = 0.1414213562373095e1 * acegen_scratch__179__;
  acegen_scratch__180__ = 2e0 * (acegen_scratch__10__ * acegen_scratch__16__ +
                                 acegen_scratch__12__ * acegen_scratch__18__);
  acegen_scratch__22__  = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__17__ * acegen_scratch__17__);
  acegen_scratch__23__ = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__213__ =
    -0.7071067811865476e0 * (acegen_scratch__22__ + acegen_scratch__23__);
  acegen_scratch__24__  = 0.1414213562373095e1 * acegen_scratch__210__;
  acegen_scratch__182__ = -2e0 * acegen_scratch__179__ * acegen_scratch__210__ +
                          acegen_scratch__180__ * acegen_scratch__22__ +
                          acegen_scratch__178__ * acegen_scratch__23__;
  acegen_scratch__26__ = -(acegen_scratch__210__ * acegen_scratch__210__) +
                         acegen_scratch__22__ * acegen_scratch__23__;
  acegen_scratch__218__ = -1e0 + acegen_scratch__26__;
  acegen_scratch__217__ = acegen_scratch__47__ / (2e0 * acegen_scratch__26__);
  acegen_scratch__222__ = acegen_scratch__217__ * acegen_scratch__218__;
  acegen_scratch__214__ =
    acegen_scratch__182__ / (acegen_scratch__26__ * acegen_scratch__26__);
  acegen_scratch__219__ = acegen_scratch__214__ * acegen_scratch__26__;
  acegen_scratch__27__  = CBRT::halley(acegen_scratch__26__);
  acegen_scratch__211__ = acegen_scratch__27__ / (3e0 * acegen_scratch__26__);
  acegen_scratch__183__ = acegen_scratch__182__ * acegen_scratch__211__;
  acegen_scratch__185__ = (-2e0 * acegen_scratch__183__) / acegen_scratch__26__;
  acegen_scratch__212__ =
    (2e0 / 9e0) * acegen_scratch__214__ * acegen_scratch__27__;
  acegen_scratch__193__ = acegen_scratch__180__ * acegen_scratch__211__ -
                          acegen_scratch__212__ * acegen_scratch__23__;
  acegen_scratch__39__  = 1e0 / (acegen_scratch__27__ * acegen_scratch__27__);
  acegen_scratch__223__ = acegen_scratch__213__ * acegen_scratch__39__;
  acegen_scratch__215__ = acegen_scratch__178__ * acegen_scratch__39__;
  acegen_scratch__216__ =
    acegen_scratch__215__ + acegen_scratch__185__ * acegen_scratch__22__;
  acegen_scratch__37__ = -(acegen_scratch__211__ * acegen_scratch__24__);
  acegen_scratch__203__ =
    acegen_scratch__13__ *
      (acegen_scratch__223__ *
         (-(acegen_scratch__181__ * acegen_scratch__211__) +
          acegen_scratch__212__ * acegen_scratch__24__) +
       acegen_scratch__37__ *
         (acegen_scratch__185__ * acegen_scratch__213__ -
          0.7071067811865476e0 *
            (acegen_scratch__178__ + acegen_scratch__180__) *
            acegen_scratch__39__)) +
    (-0.3535533905932738e0 * acegen_scratch__214__ * acegen_scratch__24__ +
     acegen_scratch__181__ *
       (-0.3535533905932738e0 + 0.3535533905932738e0 / acegen_scratch__26__)) *
      acegen_scratch__47__;
  acegen_scratch__36__  = acegen_scratch__211__ * acegen_scratch__22__;
  acegen_scratch__35__  = acegen_scratch__211__ * acegen_scratch__23__;
  acegen_scratch__221__ = -(acegen_scratch__35__ * acegen_scratch__39__);
  acegen_scratch__195__ =
    -(acegen_scratch__216__ * acegen_scratch__35__) -
    (acegen_scratch__183__ + acegen_scratch__193__ * acegen_scratch__22__) *
      acegen_scratch__39__;
  acegen_scratch__201__ =
    acegen_scratch__217__ * (acegen_scratch__178__ * acegen_scratch__218__ +
                             acegen_scratch__219__ * acegen_scratch__22__) +
    acegen_scratch__13__ *
      (acegen_scratch__195__ - acegen_scratch__215__ * acegen_scratch__36__ -
       acegen_scratch__216__ * acegen_scratch__36__ +
       acegen_scratch__212__ * (acegen_scratch__22__ * acegen_scratch__22__) *
         acegen_scratch__39__);
  acegen_scratch__200__ =
    acegen_scratch__217__ * (acegen_scratch__180__ * acegen_scratch__218__ +
                             acegen_scratch__219__ * acegen_scratch__23__) -
    acegen_scratch__13__ *
      (-acegen_scratch__195__ + 2e0 * acegen_scratch__193__ *
                                  acegen_scratch__23__ * acegen_scratch__39__);
  acegen_scratch__38__ =
    acegen_scratch__22__ * acegen_scratch__221__ + 1e0 / acegen_scratch__27__;
  acegen_scratch__51__ =
    acegen_scratch__222__ * acegen_scratch__23__ +
    acegen_scratch__13__ *
      (acegen_scratch__221__ * acegen_scratch__23__ + acegen_scratch__38__);
  acegen_scratch__52__ =
    acegen_scratch__22__ * acegen_scratch__222__ +
    acegen_scratch__13__ *
      (acegen_scratch__38__ -
       acegen_scratch__22__ * acegen_scratch__36__ * acegen_scratch__39__);
  acegen_scratch__54__ =
    2e0 * acegen_scratch__217__ * acegen_scratch__24__ *
      (0.3535533905932738e0 - 0.3535533905932738e0 * acegen_scratch__26__) +
    acegen_scratch__13__ * acegen_scratch__223__ * acegen_scratch__37__;
  gradientOut[0][0] = acegen_scratch__15__ * acegen_scratch__200__ +
                      acegen_scratch__16__ * acegen_scratch__203__ +
                      acegen_scratch__10__ * acegen_scratch__54__ +
                      acegen_scratch__51__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__201__ +
                      acegen_scratch__15__ * acegen_scratch__203__ +
                      acegen_scratch__10__ * acegen_scratch__52__ +
                      acegen_scratch__54__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__17__ * acegen_scratch__200__ +
                      acegen_scratch__18__ * acegen_scratch__203__ +
                      acegen_scratch__11__ * acegen_scratch__51__ +
                      acegen_scratch__12__ * acegen_scratch__54__;
  gradientOut[1][1] = acegen_scratch__18__ * acegen_scratch__201__ +
                      acegen_scratch__17__ * acegen_scratch__203__ +
                      acegen_scratch__12__ * acegen_scratch__52__ +
                      acegen_scratch__11__ * acegen_scratch__54__;
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
  Number acegen_scratch__10__, acegen_scratch__104__, acegen_scratch__105__,
    acegen_scratch__106__, acegen_scratch__107__, acegen_scratch__108__,
    acegen_scratch__109__, acegen_scratch__11__, acegen_scratch__110__,
    acegen_scratch__111__, acegen_scratch__112__, acegen_scratch__113__,
    acegen_scratch__114__, acegen_scratch__12__, acegen_scratch__16__,
    acegen_scratch__17__, acegen_scratch__18__, acegen_scratch__20__,
    acegen_scratch__28__, acegen_scratch__30__, acegen_scratch__31__,
    acegen_scratch__32__, acegen_scratch__40__, acegen_scratch__44__,
    acegen_scratch__45__, acegen_scratch__47__, acegen_scratch__54__,
    acegen_scratch__58__, acegen_scratch__61__, acegen_scratch__69__,
    acegen_scratch__7__, acegen_scratch__71__, acegen_scratch__73__,
    acegen_scratch__74__, acegen_scratch__76__, acegen_scratch__80__,
    acegen_scratch__83__, acegen_scratch__85__, acegen_scratch__9__;
  acegen_scratch__7__   = mu;
  acegen_scratch__107__ = acegen_scratch__7__ / 2e0;
  acegen_scratch__40__  = (2e0 / 3e0) * acegen_scratch__7__ + lambda;
  acegen_scratch__111__ = acegen_scratch__40__ / 2e0;
  acegen_scratch__9__   = 1e0 + graduIn[0][0];
  acegen_scratch__10__  = graduIn[0][1];
  acegen_scratch__11__  = graduIn[1][0];
  acegen_scratch__12__  = 1e0 + graduIn[1][1];
  acegen_scratch__16__  = (acegen_scratch__11__ * acegen_scratch__11__) +
                         (acegen_scratch__9__ * acegen_scratch__9__);
  acegen_scratch__114__ = (acegen_scratch__16__ * acegen_scratch__16__);
  acegen_scratch__17__  = (acegen_scratch__10__ * acegen_scratch__10__) +
                         (acegen_scratch__12__ * acegen_scratch__12__);
  acegen_scratch__112__ = acegen_scratch__16__ * acegen_scratch__17__;
  acegen_scratch__109__ = (acegen_scratch__17__ * acegen_scratch__17__);
  acegen_scratch__18__ =
    0.1414213562373095e1 * (acegen_scratch__11__ * acegen_scratch__12__ +
                            acegen_scratch__10__ * acegen_scratch__9__);
  acegen_scratch__105__ = (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__110__ = 0.5e0 * acegen_scratch__105__;
  acegen_scratch__20__  = -acegen_scratch__110__ + acegen_scratch__112__;
  acegen_scratch__83__  = 1e0 / (2e0 * acegen_scratch__20__);
  acegen_scratch__106__ =
    (-1e0 + acegen_scratch__20__) * acegen_scratch__40__ * acegen_scratch__83__;
  acegen_scratch__71__ = 1e0 / (acegen_scratch__20__ * acegen_scratch__20__);
  acegen_scratch__85__ = (acegen_scratch__40__ * acegen_scratch__71__) / 4e0;
  acegen_scratch__69__ = -(acegen_scratch__18__ * acegen_scratch__85__);
  acegen_scratch__54__ =
    1e0 / Power(acegen_scratch__20__, 0.23333333333333334e1);
  acegen_scratch__104__ = (-4e0 / 3e0) * acegen_scratch__54__;
  acegen_scratch__113__ = (acegen_scratch__104__ * acegen_scratch__16__) / 3e0;
  acegen_scratch__58__  = acegen_scratch__113__ * acegen_scratch__18__;
  acegen_scratch__108__ = (acegen_scratch__104__ * acegen_scratch__17__) / 3e0;
  acegen_scratch__73__ =
    1e0 / (3e0 * Power(acegen_scratch__20__, 0.13333333333333333e1));
  acegen_scratch__74__ =
    -(acegen_scratch__108__ * acegen_scratch__16__) - acegen_scratch__73__;
  acegen_scratch__31__ = acegen_scratch__18__ * acegen_scratch__73__;
  acegen_scratch__61__ =
    acegen_scratch__31__ + acegen_scratch__17__ * acegen_scratch__58__;
  acegen_scratch__30__ = -(acegen_scratch__16__ * acegen_scratch__73__);
  acegen_scratch__76__ =
    acegen_scratch__30__ + acegen_scratch__16__ * acegen_scratch__74__;
  acegen_scratch__28__ = -(acegen_scratch__17__ * acegen_scratch__73__);
  acegen_scratch__80__ =
    acegen_scratch__28__ + acegen_scratch__17__ * acegen_scratch__74__;
  acegen_scratch__32__ =
    1e0 / Power(acegen_scratch__20__, 0.3333333333333333e0) +
    acegen_scratch__16__ * acegen_scratch__28__;
  acegen_scratch__44__ =
    acegen_scratch__106__ * acegen_scratch__17__ +
    (acegen_scratch__17__ * acegen_scratch__28__ + acegen_scratch__32__) *
      acegen_scratch__7__;
  acegen_scratch__45__ =
    acegen_scratch__106__ * acegen_scratch__16__ +
    (acegen_scratch__16__ * acegen_scratch__30__ + acegen_scratch__32__) *
      acegen_scratch__7__;
  acegen_scratch__47__ =
    (acegen_scratch__18__ *
     (0.3535533905932738e0 - 0.3535533905932738e0 * acegen_scratch__20__) *
     acegen_scratch__40__) /
      acegen_scratch__20__ +
    0.7071067811865476e0 * (acegen_scratch__16__ + acegen_scratch__17__) *
      acegen_scratch__31__ * acegen_scratch__7__;
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  gradientOut[0][0] = acegen_scratch__10__ * acegen_scratch__47__ +
                      acegen_scratch__44__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__10__ * acegen_scratch__45__ +
                      acegen_scratch__47__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__11__ * acegen_scratch__44__ +
                      acegen_scratch__12__ * acegen_scratch__47__;
  gradientOut[1][1] = acegen_scratch__12__ * acegen_scratch__45__ +
                      acegen_scratch__11__ * acegen_scratch__47__;
  cacheOut[0] = 4e0 * (acegen_scratch__107__ *
                         (-(acegen_scratch__108__ * acegen_scratch__109__) +
                          acegen_scratch__80__) +
                       acegen_scratch__109__ * acegen_scratch__85__);
  cacheOut[1] =
    4e0 *
    (acegen_scratch__107__ * (acegen_scratch__76__ + acegen_scratch__80__) +
     acegen_scratch__111__ *
       (0.5e0 + (acegen_scratch__112__ * acegen_scratch__71__) / 2e0 -
        acegen_scratch__83__));
  cacheOut[2] =
    4e0 *
    (acegen_scratch__107__ *
       (acegen_scratch__108__ * acegen_scratch__17__ * acegen_scratch__18__ +
        acegen_scratch__61__) +
     acegen_scratch__17__ * acegen_scratch__69__);
  cacheOut[3] = 4e0 * (acegen_scratch__107__ *
                         (-(acegen_scratch__113__ * acegen_scratch__114__) +
                          acegen_scratch__76__) +
                       acegen_scratch__114__ * acegen_scratch__85__);
  cacheOut[4] =
    4e0 *
    (acegen_scratch__107__ *
       (acegen_scratch__16__ * acegen_scratch__58__ + acegen_scratch__61__) +
     acegen_scratch__16__ * acegen_scratch__69__);
  cacheOut[5] =
    4e0 *
    (acegen_scratch__107__ * (acegen_scratch__16__ + acegen_scratch__17__) *
       ((4e0 / 9e0) * acegen_scratch__105__ * acegen_scratch__54__ +
        acegen_scratch__73__) +
     acegen_scratch__111__ *
       (-0.5e0 + acegen_scratch__110__ * acegen_scratch__71__ +
        acegen_scratch__83__));
  cacheOut[6] = acegen_scratch__44__;
  cacheOut[7] = acegen_scratch__47__;
  cacheOut[8] = acegen_scratch__45__;
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
    acegen_scratch__38__, acegen_scratch__39__, acegen_scratch__50__,
    acegen_scratch__52__, acegen_scratch__53__, acegen_scratch__54__,
    acegen_scratch__55__, acegen_scratch__56__, acegen_scratch__57__,
    acegen_scratch__58__, acegen_scratch__62__, acegen_scratch__63__,
    acegen_scratch__64__, acegen_scratch__79__, acegen_scratch__80__,
    acegen_scratch__81__;
  acegen_scratch__13__ = (mu);
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
  acegen_scratch__56__ = (acegen_scratch__36__ * acegen_scratch__36__);
  acegen_scratch__37__ = acegen_scratch__15__ * acegen_scratch__17__ +
                         acegen_scratch__18__ * acegen_scratch__20__ +
                         acegen_scratch__21__ * acegen_scratch__23__;
  acegen_scratch__58__ = 2e0 * acegen_scratch__36__ * acegen_scratch__37__;
  acegen_scratch__54__ = (acegen_scratch__37__ * acegen_scratch__37__);
  acegen_scratch__38__ = acegen_scratch__16__ * acegen_scratch__17__ +
                         acegen_scratch__19__ * acegen_scratch__20__ +
                         acegen_scratch__22__ * acegen_scratch__23__;
  acegen_scratch__81__ = -1e0 * acegen_scratch__38__;
  acegen_scratch__79__ = acegen_scratch__31__ * acegen_scratch__32__ -
                         (acegen_scratch__38__ * acegen_scratch__38__);
  acegen_scratch__39__ = -(acegen_scratch__31__ * acegen_scratch__54__) -
                         acegen_scratch__32__ * acegen_scratch__56__ +
                         acegen_scratch__38__ * acegen_scratch__58__ +
                         acegen_scratch__30__ * acegen_scratch__79__;
  acegen_scratch__53__ =
    acegen_scratch__13__ /
    (2e0 * Power(acegen_scratch__39__, 0.3333333333333333e0));
  acegen_scratch__50__ =
    (-1e0 / 6e0) *
      (acegen_scratch__13__ *
       (acegen_scratch__30__ + acegen_scratch__31__ + acegen_scratch__32__)) /
      Power(acegen_scratch__39__, 0.13333333333333333e1) +
    ((-1e0 + acegen_scratch__39__) *
     (2e0 * acegen_scratch__13__ + 3e0 * (lambda))) /
      (12e0 * acegen_scratch__39__);
  acegen_scratch__80__ = -2e0 * acegen_scratch__50__;
  acegen_scratch__52__ =
    2e0 * (acegen_scratch__53__ + acegen_scratch__50__ * acegen_scratch__79__);
  acegen_scratch__55__ =
    2e0 * (acegen_scratch__53__ +
           acegen_scratch__50__ * (acegen_scratch__30__ * acegen_scratch__32__ -
                                   acegen_scratch__54__));
  acegen_scratch__57__ =
    2e0 * (acegen_scratch__53__ +
           acegen_scratch__50__ * (acegen_scratch__30__ * acegen_scratch__31__ -
                                   acegen_scratch__56__));
  acegen_scratch__62__ =
    acegen_scratch__80__ * (acegen_scratch__32__ * acegen_scratch__36__ +
                            acegen_scratch__37__ * acegen_scratch__81__);
  acegen_scratch__63__ =
    acegen_scratch__80__ * (acegen_scratch__31__ * acegen_scratch__37__ +
                            acegen_scratch__36__ * acegen_scratch__81__);
  acegen_scratch__64__ =
    1e0 * acegen_scratch__50__ *
    (-2e0 * acegen_scratch__30__ * acegen_scratch__38__ + acegen_scratch__58__);
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  valueOut[2]       = 0e0;
  gradientOut[0][0] = acegen_scratch__15__ * acegen_scratch__52__ +
                      acegen_scratch__16__ * acegen_scratch__62__ +
                      acegen_scratch__17__ * acegen_scratch__63__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__55__ +
                      acegen_scratch__15__ * acegen_scratch__62__ +
                      acegen_scratch__17__ * acegen_scratch__64__;
  gradientOut[0][2] = acegen_scratch__17__ * acegen_scratch__57__ +
                      acegen_scratch__15__ * acegen_scratch__63__ +
                      acegen_scratch__16__ * acegen_scratch__64__;
  gradientOut[1][0] = acegen_scratch__18__ * acegen_scratch__52__ +
                      acegen_scratch__19__ * acegen_scratch__62__ +
                      acegen_scratch__20__ * acegen_scratch__63__;
  gradientOut[1][1] = acegen_scratch__19__ * acegen_scratch__55__ +
                      acegen_scratch__18__ * acegen_scratch__62__ +
                      acegen_scratch__20__ * acegen_scratch__64__;
  gradientOut[1][2] = acegen_scratch__20__ * acegen_scratch__57__ +
                      acegen_scratch__18__ * acegen_scratch__63__ +
                      acegen_scratch__19__ * acegen_scratch__64__;
  gradientOut[2][0] = acegen_scratch__21__ * acegen_scratch__52__ +
                      acegen_scratch__22__ * acegen_scratch__62__ +
                      acegen_scratch__23__ * acegen_scratch__63__;
  gradientOut[2][1] = acegen_scratch__22__ * acegen_scratch__55__ +
                      acegen_scratch__21__ * acegen_scratch__62__ +
                      acegen_scratch__23__ * acegen_scratch__64__;
  gradientOut[2][2] = acegen_scratch__23__ * acegen_scratch__57__ +
                      acegen_scratch__21__ * acegen_scratch__63__ +
                      acegen_scratch__22__ * acegen_scratch__64__;
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
    acegen_scratch__181__, acegen_scratch__19__, acegen_scratch__191__,
    acegen_scratch__20__, acegen_scratch__201__, acegen_scratch__21__,
    acegen_scratch__211__, acegen_scratch__22__, acegen_scratch__228__,
    acegen_scratch__23__, acegen_scratch__24__, acegen_scratch__247__,
    acegen_scratch__25__, acegen_scratch__27__, acegen_scratch__28__,
    acegen_scratch__29__, acegen_scratch__30__, acegen_scratch__31__,
    acegen_scratch__32__, acegen_scratch__33__, acegen_scratch__34__,
    acegen_scratch__341__, acegen_scratch__342__, acegen_scratch__343__,
    acegen_scratch__344__, acegen_scratch__345__, acegen_scratch__346__,
    acegen_scratch__35__, acegen_scratch__352__, acegen_scratch__353__,
    acegen_scratch__358__, acegen_scratch__360__, acegen_scratch__361__,
    acegen_scratch__364__, acegen_scratch__366__, acegen_scratch__367__,
    acegen_scratch__368__, acegen_scratch__369__, acegen_scratch__370__,
    acegen_scratch__374__, acegen_scratch__375__, acegen_scratch__376__,
    acegen_scratch__390__, acegen_scratch__392__, acegen_scratch__393__,
    acegen_scratch__394__, acegen_scratch__395__, acegen_scratch__396__,
    acegen_scratch__397__, acegen_scratch__398__, acegen_scratch__399__,
    acegen_scratch__403__, acegen_scratch__42__, acegen_scratch__43__,
    acegen_scratch__44__, acegen_scratch__48__, acegen_scratch__49__,
    acegen_scratch__50__, acegen_scratch__51__, acegen_scratch__52__,
    acegen_scratch__65__, acegen_scratch__66__, acegen_scratch__67__,
    acegen_scratch__68__, acegen_scratch__69__, acegen_scratch__70__,
    acegen_scratch__71__, acegen_scratch__72__, acegen_scratch__73__,
    acegen_scratch__74__, acegen_scratch__78__, acegen_scratch__79__,
    acegen_scratch__80__;
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
  acegen_scratch__396__ = acegen_scratch__25__ / 2e0;
  acegen_scratch__397__ = 3e0 * ((2e0 / 3e0) * acegen_scratch__25__ + (lambda));
  acegen_scratch__27__  = 1e0 + graduIn[0][0];
  acegen_scratch__28__  = graduIn[0][1];
  acegen_scratch__29__  = graduIn[0][2];
  acegen_scratch__30__  = graduIn[1][0];
  acegen_scratch__31__  = 1e0 + graduIn[1][1];
  acegen_scratch__32__  = graduIn[1][2];
  acegen_scratch__33__  = graduIn[2][0];
  acegen_scratch__341__ = 2e0 * (acegen_scratch__16__ * acegen_scratch__27__ +
                                 acegen_scratch__19__ * acegen_scratch__30__ +
                                 acegen_scratch__22__ * acegen_scratch__33__);
  acegen_scratch__34__  = graduIn[2][1];
  acegen_scratch__342__ = acegen_scratch__17__ * acegen_scratch__27__ +
                          acegen_scratch__16__ * acegen_scratch__28__ +
                          acegen_scratch__20__ * acegen_scratch__30__ +
                          acegen_scratch__19__ * acegen_scratch__31__ +
                          acegen_scratch__23__ * acegen_scratch__33__ +
                          acegen_scratch__22__ * acegen_scratch__34__;
  acegen_scratch__344__ = 2e0 * (acegen_scratch__17__ * acegen_scratch__28__ +
                                 acegen_scratch__20__ * acegen_scratch__31__ +
                                 acegen_scratch__23__ * acegen_scratch__34__);
  acegen_scratch__35__  = 1e0 + graduIn[2][2];
  acegen_scratch__345__ = acegen_scratch__18__ * acegen_scratch__28__ +
                          acegen_scratch__17__ * acegen_scratch__29__ +
                          acegen_scratch__21__ * acegen_scratch__31__ +
                          acegen_scratch__20__ * acegen_scratch__32__ +
                          acegen_scratch__24__ * acegen_scratch__34__ +
                          acegen_scratch__23__ * acegen_scratch__35__;
  acegen_scratch__343__ = acegen_scratch__18__ * acegen_scratch__27__ +
                          acegen_scratch__16__ * acegen_scratch__29__ +
                          acegen_scratch__21__ * acegen_scratch__30__ +
                          acegen_scratch__19__ * acegen_scratch__32__ +
                          acegen_scratch__24__ * acegen_scratch__33__ +
                          acegen_scratch__22__ * acegen_scratch__35__;
  acegen_scratch__346__ = 2e0 * (acegen_scratch__18__ * acegen_scratch__29__ +
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
  acegen_scratch__395__ =
    acegen_scratch__42__ + acegen_scratch__43__ + acegen_scratch__44__;
  acegen_scratch__48__ = acegen_scratch__27__ * acegen_scratch__28__ +
                         acegen_scratch__30__ * acegen_scratch__31__ +
                         acegen_scratch__33__ * acegen_scratch__34__;
  acegen_scratch__392__ = 2e0 * acegen_scratch__48__;
  acegen_scratch__403__ = -(acegen_scratch__342__ * acegen_scratch__392__) +
                          acegen_scratch__344__ * acegen_scratch__42__ +
                          acegen_scratch__341__ * acegen_scratch__43__;
  acegen_scratch__72__ = (acegen_scratch__48__ * acegen_scratch__48__);
  acegen_scratch__201__ =
    acegen_scratch__42__ * acegen_scratch__43__ - acegen_scratch__72__;
  acegen_scratch__49__ = acegen_scratch__27__ * acegen_scratch__29__ +
                         acegen_scratch__30__ * acegen_scratch__32__ +
                         acegen_scratch__33__ * acegen_scratch__35__;
  acegen_scratch__390__ = 2e0 * acegen_scratch__49__;
  acegen_scratch__353__ = acegen_scratch__343__ * acegen_scratch__390__;
  acegen_scratch__352__ = 2e0 * (acegen_scratch__343__ * acegen_scratch__48__ +
                                 acegen_scratch__342__ * acegen_scratch__49__);
  acegen_scratch__74__  = acegen_scratch__390__ * acegen_scratch__48__;
  acegen_scratch__70__  = (acegen_scratch__49__ * acegen_scratch__49__);
  acegen_scratch__191__ =
    acegen_scratch__42__ * acegen_scratch__44__ - acegen_scratch__70__;
  acegen_scratch__50__ = acegen_scratch__28__ * acegen_scratch__29__ +
                         acegen_scratch__31__ * acegen_scratch__32__ +
                         acegen_scratch__34__ * acegen_scratch__35__;
  acegen_scratch__393__ = 2e0 * acegen_scratch__50__;
  acegen_scratch__358__ = acegen_scratch__345__ * acegen_scratch__393__;
  acegen_scratch__247__ = -(acegen_scratch__392__ * acegen_scratch__44__) +
                          acegen_scratch__390__ * acegen_scratch__50__;
  acegen_scratch__228__ = -(acegen_scratch__390__ * acegen_scratch__43__) +
                          acegen_scratch__392__ * acegen_scratch__50__;
  acegen_scratch__211__ =
    -(acegen_scratch__393__ * acegen_scratch__42__) + acegen_scratch__74__;
  acegen_scratch__67__  = (acegen_scratch__50__ * acegen_scratch__50__);
  acegen_scratch__360__ = acegen_scratch__201__ * acegen_scratch__346__ -
                          acegen_scratch__358__ * acegen_scratch__42__ -
                          acegen_scratch__353__ * acegen_scratch__43__ +
                          acegen_scratch__403__ * acegen_scratch__44__ +
                          acegen_scratch__352__ * acegen_scratch__50__ -
                          acegen_scratch__341__ * acegen_scratch__67__ -
                          acegen_scratch__344__ * acegen_scratch__70__ +
                          acegen_scratch__345__ * acegen_scratch__74__;
  acegen_scratch__181__ =
    acegen_scratch__43__ * acegen_scratch__44__ - acegen_scratch__67__;
  acegen_scratch__51__ = acegen_scratch__181__ * acegen_scratch__42__ -
                         acegen_scratch__43__ * acegen_scratch__70__ -
                         acegen_scratch__44__ * acegen_scratch__72__ +
                         acegen_scratch__50__ * acegen_scratch__74__;
  acegen_scratch__394__ = -2e0 / acegen_scratch__51__;
  acegen_scratch__52__  = CBRT::halley(acegen_scratch__51__);
  acegen_scratch__361__ = (acegen_scratch__360__ * acegen_scratch__52__) /
                          (3e0 * acegen_scratch__51__);
  acegen_scratch__65__  = 1e0 / (acegen_scratch__52__ * acegen_scratch__52__);
  acegen_scratch__364__ = -(acegen_scratch__395__ * acegen_scratch__65__);
  acegen_scratch__398__ =
    -2e0 * acegen_scratch__25__ * acegen_scratch__364__ * acegen_scratch__52__;
  acegen_scratch__367__ =
    ((acegen_scratch__360__ * (acegen_scratch__397__ + acegen_scratch__398__)) /
       (acegen_scratch__51__ * acegen_scratch__51__) +
     acegen_scratch__25__ * acegen_scratch__394__ *
       (-(acegen_scratch__361__ * acegen_scratch__364__) +
        acegen_scratch__52__ * (acegen_scratch__361__ * acegen_scratch__394__ *
                                  acegen_scratch__395__ +
                                (acegen_scratch__341__ + acegen_scratch__344__ +
                                 acegen_scratch__346__) *
                                  acegen_scratch__65__))) /
    12e0;
  acegen_scratch__366__ =
    -(acegen_scratch__361__ * acegen_scratch__396__ * acegen_scratch__65__);
  acegen_scratch__69__ = acegen_scratch__396__ / acegen_scratch__52__;
  acegen_scratch__66__ =
    (-acegen_scratch__398__ +
     acegen_scratch__397__ * (-1e0 + acegen_scratch__51__)) /
    (12e0 * acegen_scratch__51__);
  acegen_scratch__399__ = 2e0 * acegen_scratch__66__;
  acegen_scratch__374__ =
    acegen_scratch__247__ * acegen_scratch__367__ +
    acegen_scratch__399__ * (-(acegen_scratch__342__ * acegen_scratch__44__) -
                             acegen_scratch__346__ * acegen_scratch__48__ +
                             acegen_scratch__345__ * acegen_scratch__49__ +
                             acegen_scratch__343__ * acegen_scratch__50__);
  acegen_scratch__375__ =
    acegen_scratch__228__ * acegen_scratch__367__ +
    acegen_scratch__399__ * (-(acegen_scratch__343__ * acegen_scratch__43__) +
                             acegen_scratch__345__ * acegen_scratch__48__ -
                             acegen_scratch__344__ * acegen_scratch__49__ +
                             acegen_scratch__342__ * acegen_scratch__50__);
  acegen_scratch__376__ =
    1e0 * acegen_scratch__211__ * acegen_scratch__367__ +
    (acegen_scratch__352__ - acegen_scratch__341__ * acegen_scratch__393__ -
     2e0 * acegen_scratch__345__ * acegen_scratch__42__) *
      acegen_scratch__66__;
  acegen_scratch__370__ = 2e0 * (acegen_scratch__366__ +
                                 acegen_scratch__201__ * acegen_scratch__367__ +
                                 acegen_scratch__403__ * acegen_scratch__66__);
  acegen_scratch__369__ =
    2e0 *
    (acegen_scratch__366__ + acegen_scratch__191__ * acegen_scratch__367__ +
     (-acegen_scratch__353__ + acegen_scratch__346__ * acegen_scratch__42__ +
      acegen_scratch__341__ * acegen_scratch__44__) *
       acegen_scratch__66__);
  acegen_scratch__368__ =
    2e0 *
    (acegen_scratch__366__ + acegen_scratch__181__ * acegen_scratch__367__ +
     (-acegen_scratch__358__ + acegen_scratch__346__ * acegen_scratch__43__ +
      acegen_scratch__344__ * acegen_scratch__44__) *
       acegen_scratch__66__);
  acegen_scratch__68__ =
    2e0 * (acegen_scratch__181__ * acegen_scratch__66__ + acegen_scratch__69__);
  acegen_scratch__71__ =
    2e0 * (acegen_scratch__191__ * acegen_scratch__66__ + acegen_scratch__69__);
  acegen_scratch__73__ =
    2e0 * (acegen_scratch__201__ * acegen_scratch__66__ + acegen_scratch__69__);
  acegen_scratch__78__ = acegen_scratch__247__ * acegen_scratch__66__;
  acegen_scratch__79__ = acegen_scratch__228__ * acegen_scratch__66__;
  acegen_scratch__80__ = acegen_scratch__211__ * acegen_scratch__66__;
  gradientOut[0][0]    = acegen_scratch__27__ * acegen_scratch__368__ +
                      acegen_scratch__28__ * acegen_scratch__374__ +
                      acegen_scratch__29__ * acegen_scratch__375__ +
                      acegen_scratch__16__ * acegen_scratch__68__ +
                      acegen_scratch__17__ * acegen_scratch__78__ +
                      acegen_scratch__18__ * acegen_scratch__79__;
  gradientOut[0][1] = acegen_scratch__28__ * acegen_scratch__369__ +
                      acegen_scratch__27__ * acegen_scratch__374__ +
                      acegen_scratch__29__ * acegen_scratch__376__ +
                      acegen_scratch__17__ * acegen_scratch__71__ +
                      acegen_scratch__16__ * acegen_scratch__78__ +
                      acegen_scratch__18__ * acegen_scratch__80__;
  gradientOut[0][2] = acegen_scratch__29__ * acegen_scratch__370__ +
                      acegen_scratch__27__ * acegen_scratch__375__ +
                      acegen_scratch__28__ * acegen_scratch__376__ +
                      acegen_scratch__18__ * acegen_scratch__73__ +
                      acegen_scratch__16__ * acegen_scratch__79__ +
                      acegen_scratch__17__ * acegen_scratch__80__;
  gradientOut[1][0] = acegen_scratch__30__ * acegen_scratch__368__ +
                      acegen_scratch__31__ * acegen_scratch__374__ +
                      acegen_scratch__32__ * acegen_scratch__375__ +
                      acegen_scratch__19__ * acegen_scratch__68__ +
                      acegen_scratch__20__ * acegen_scratch__78__ +
                      acegen_scratch__21__ * acegen_scratch__79__;
  gradientOut[1][1] = acegen_scratch__31__ * acegen_scratch__369__ +
                      acegen_scratch__30__ * acegen_scratch__374__ +
                      acegen_scratch__32__ * acegen_scratch__376__ +
                      acegen_scratch__20__ * acegen_scratch__71__ +
                      acegen_scratch__19__ * acegen_scratch__78__ +
                      acegen_scratch__21__ * acegen_scratch__80__;
  gradientOut[1][2] = acegen_scratch__32__ * acegen_scratch__370__ +
                      acegen_scratch__30__ * acegen_scratch__375__ +
                      acegen_scratch__31__ * acegen_scratch__376__ +
                      acegen_scratch__21__ * acegen_scratch__73__ +
                      acegen_scratch__19__ * acegen_scratch__79__ +
                      acegen_scratch__20__ * acegen_scratch__80__;
  gradientOut[2][0] = acegen_scratch__33__ * acegen_scratch__368__ +
                      acegen_scratch__34__ * acegen_scratch__374__ +
                      acegen_scratch__35__ * acegen_scratch__375__ +
                      acegen_scratch__22__ * acegen_scratch__68__ +
                      acegen_scratch__23__ * acegen_scratch__78__ +
                      acegen_scratch__24__ * acegen_scratch__79__;
  gradientOut[2][1] = acegen_scratch__34__ * acegen_scratch__369__ +
                      acegen_scratch__33__ * acegen_scratch__374__ +
                      acegen_scratch__35__ * acegen_scratch__376__ +
                      acegen_scratch__23__ * acegen_scratch__71__ +
                      acegen_scratch__22__ * acegen_scratch__78__ +
                      acegen_scratch__24__ * acegen_scratch__80__;
  gradientOut[2][2] = acegen_scratch__35__ * acegen_scratch__370__ +
                      acegen_scratch__33__ * acegen_scratch__375__ +
                      acegen_scratch__34__ * acegen_scratch__376__ +
                      acegen_scratch__24__ * acegen_scratch__73__ +
                      acegen_scratch__22__ * acegen_scratch__79__ +
                      acegen_scratch__23__ * acegen_scratch__80__;
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
  Number acegen_scratch__101__, acegen_scratch__102__, acegen_scratch__103__,
    acegen_scratch__104__, acegen_scratch__105__, acegen_scratch__121__,
    acegen_scratch__125__, acegen_scratch__128__, acegen_scratch__13__,
    acegen_scratch__15__, acegen_scratch__154__, acegen_scratch__155__,
    acegen_scratch__156__, acegen_scratch__157__, acegen_scratch__158__,
    acegen_scratch__159__, acegen_scratch__16__, acegen_scratch__160__,
    acegen_scratch__161__, acegen_scratch__162__, acegen_scratch__163__,
    acegen_scratch__164__, acegen_scratch__165__, acegen_scratch__166__,
    acegen_scratch__167__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__19__, acegen_scratch__20__, acegen_scratch__21__,
    acegen_scratch__22__, acegen_scratch__23__, acegen_scratch__30__,
    acegen_scratch__31__, acegen_scratch__32__, acegen_scratch__33__,
    acegen_scratch__34__, acegen_scratch__35__, acegen_scratch__39__,
    acegen_scratch__49__, acegen_scratch__50__, acegen_scratch__52__,
    acegen_scratch__53__, acegen_scratch__54__, acegen_scratch__55__,
    acegen_scratch__56__, acegen_scratch__57__, acegen_scratch__58__,
    acegen_scratch__62__, acegen_scratch__63__, acegen_scratch__64__,
    acegen_scratch__77__, acegen_scratch__78__, acegen_scratch__79__,
    acegen_scratch__80__, acegen_scratch__81__, acegen_scratch__82__,
    acegen_scratch__91__, acegen_scratch__92__, acegen_scratch__93__,
    acegen_scratch__94__, acegen_scratch__95__, acegen_scratch__96__,
    acegen_scratch__97__;
  acegen_scratch__13__  = (mu);
  acegen_scratch__49__  = (2e0 / 3e0) * acegen_scratch__13__ + (lambda);
  acegen_scratch__15__  = 1e0 + graduIn[0][0];
  acegen_scratch__16__  = graduIn[0][1];
  acegen_scratch__17__  = graduIn[0][2];
  acegen_scratch__18__  = graduIn[1][0];
  acegen_scratch__19__  = 1e0 + graduIn[1][1];
  acegen_scratch__20__  = graduIn[1][2];
  acegen_scratch__21__  = graduIn[2][0];
  acegen_scratch__22__  = graduIn[2][1];
  acegen_scratch__154__ = acegen_scratch__15__ * acegen_scratch__16__ +
                          acegen_scratch__18__ * acegen_scratch__19__ +
                          acegen_scratch__21__ * acegen_scratch__22__;
  acegen_scratch__23__  = 1e0 + graduIn[2][2];
  acegen_scratch__156__ = acegen_scratch__16__ * acegen_scratch__17__ +
                          acegen_scratch__19__ * acegen_scratch__20__ +
                          acegen_scratch__22__ * acegen_scratch__23__;
  acegen_scratch__155__ = acegen_scratch__15__ * acegen_scratch__17__ +
                          acegen_scratch__18__ * acegen_scratch__20__ +
                          acegen_scratch__21__ * acegen_scratch__23__;
  acegen_scratch__30__ = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__18__ * acegen_scratch__18__) +
                         (acegen_scratch__21__ * acegen_scratch__21__);
  acegen_scratch__163__ = -2e0 * acegen_scratch__30__;
  acegen_scratch__31__  = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__19__ * acegen_scratch__19__) +
                         (acegen_scratch__22__ * acegen_scratch__22__);
  acegen_scratch__32__ = (acegen_scratch__17__ * acegen_scratch__17__) +
                         (acegen_scratch__20__ * acegen_scratch__20__) +
                         (acegen_scratch__23__ * acegen_scratch__23__);
  acegen_scratch__97__ =
    acegen_scratch__30__ + acegen_scratch__31__ + acegen_scratch__32__;
  acegen_scratch__33__  = 0.1414213562373095e1 * acegen_scratch__156__;
  acegen_scratch__34__  = 0.1414213562373095e1 * acegen_scratch__155__;
  acegen_scratch__35__  = 0.1414213562373095e1 * acegen_scratch__154__;
  acegen_scratch__157__ = 2e0 * acegen_scratch__154__;
  acegen_scratch__56__  = (acegen_scratch__154__ * acegen_scratch__154__);
  acegen_scratch__79__ =
    acegen_scratch__30__ * acegen_scratch__31__ - acegen_scratch__56__;
  acegen_scratch__158__ = 2e0 * acegen_scratch__155__;
  acegen_scratch__58__  = acegen_scratch__155__ * acegen_scratch__157__;
  acegen_scratch__80__  = -(acegen_scratch__30__ * acegen_scratch__33__) +
                         0.7071067811865476e0 * acegen_scratch__58__;
  acegen_scratch__54__ = (acegen_scratch__155__ * acegen_scratch__155__);
  acegen_scratch__78__ =
    acegen_scratch__30__ * acegen_scratch__32__ - acegen_scratch__54__;
  acegen_scratch__128__ = acegen_scratch__156__ * acegen_scratch__158__ -
                          acegen_scratch__157__ * acegen_scratch__32__;
  acegen_scratch__125__ = acegen_scratch__156__ * acegen_scratch__157__ -
                          acegen_scratch__158__ * acegen_scratch__31__;
  acegen_scratch__167__ = 0.1414213562373095e1 * acegen_scratch__125__;
  acegen_scratch__121__ =
    acegen_scratch__156__ * acegen_scratch__163__ + acegen_scratch__58__;
  acegen_scratch__166__ = 0.1414213562373095e1 * acegen_scratch__121__;
  acegen_scratch__82__  = acegen_scratch__156__ * acegen_scratch__34__ -
                         acegen_scratch__32__ * acegen_scratch__35__;
  acegen_scratch__81__ = -(acegen_scratch__31__ * acegen_scratch__34__) +
                         acegen_scratch__156__ * acegen_scratch__35__;
  acegen_scratch__77__ = -(acegen_scratch__156__ * acegen_scratch__156__) +
                         acegen_scratch__31__ * acegen_scratch__32__;
  acegen_scratch__39__ = -(acegen_scratch__31__ * acegen_scratch__54__) -
                         acegen_scratch__32__ * acegen_scratch__56__ +
                         acegen_scratch__156__ * acegen_scratch__58__ +
                         acegen_scratch__30__ * acegen_scratch__77__;
  acegen_scratch__159__ =
    acegen_scratch__49__ /
      (4e0 * (acegen_scratch__39__ * acegen_scratch__39__)) +
    ((2e0 / 9e0) * acegen_scratch__13__ * acegen_scratch__97__) /
      Power(acegen_scratch__39__, 0.23333333333333334e1);
  acegen_scratch__105__ = acegen_scratch__159__ * acegen_scratch__82__;
  acegen_scratch__104__ = acegen_scratch__159__ * acegen_scratch__81__;
  acegen_scratch__103__ = acegen_scratch__159__ * acegen_scratch__80__;
  acegen_scratch__96__  = (-1e0 / 6e0) * acegen_scratch__13__ /
                         Power(acegen_scratch__39__, 0.13333333333333333e1);
  acegen_scratch__102__ =
    acegen_scratch__159__ * acegen_scratch__79__ + acegen_scratch__96__;
  acegen_scratch__101__ =
    acegen_scratch__159__ * acegen_scratch__78__ + acegen_scratch__96__;
  acegen_scratch__95__ = acegen_scratch__82__ * acegen_scratch__96__;
  acegen_scratch__94__ = acegen_scratch__81__ * acegen_scratch__96__;
  acegen_scratch__93__ = acegen_scratch__80__ * acegen_scratch__96__;
  acegen_scratch__92__ = acegen_scratch__79__ * acegen_scratch__96__;
  acegen_scratch__91__ = acegen_scratch__78__ * acegen_scratch__96__;
  acegen_scratch__53__ =
    acegen_scratch__13__ /
    (2e0 * Power(acegen_scratch__39__, 0.3333333333333333e0));
  acegen_scratch__50__ =
    (0.25e0 - 1e0 / (4e0 * acegen_scratch__39__)) * acegen_scratch__49__ +
    acegen_scratch__96__ * acegen_scratch__97__;
  acegen_scratch__165__ = acegen_scratch__35__ * acegen_scratch__50__;
  acegen_scratch__164__ = acegen_scratch__34__ * acegen_scratch__50__;
  acegen_scratch__162__ = acegen_scratch__31__ * acegen_scratch__50__;
  acegen_scratch__161__ = acegen_scratch__32__ * acegen_scratch__50__;
  acegen_scratch__160__ = 1e0 * acegen_scratch__50__;
  acegen_scratch__52__ =
    2e0 * (acegen_scratch__53__ + acegen_scratch__50__ * acegen_scratch__77__);
  acegen_scratch__55__ =
    2e0 * (acegen_scratch__53__ + acegen_scratch__50__ * acegen_scratch__78__);
  acegen_scratch__57__ =
    2e0 * (acegen_scratch__53__ + acegen_scratch__50__ * acegen_scratch__79__);
  acegen_scratch__62__ = acegen_scratch__128__ * acegen_scratch__160__;
  acegen_scratch__63__ = acegen_scratch__125__ * acegen_scratch__160__;
  acegen_scratch__64__ = 1e0 * acegen_scratch__121__ * acegen_scratch__160__;
  valueOut[0]          = 0e0;
  valueOut[1]          = 0e0;
  valueOut[2]          = 0e0;
  gradientOut[0][0]    = acegen_scratch__15__ * acegen_scratch__52__ +
                      acegen_scratch__16__ * acegen_scratch__62__ +
                      acegen_scratch__17__ * acegen_scratch__63__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__55__ +
                      acegen_scratch__15__ * acegen_scratch__62__ +
                      acegen_scratch__17__ * acegen_scratch__64__;
  gradientOut[0][2] = acegen_scratch__17__ * acegen_scratch__57__ +
                      acegen_scratch__15__ * acegen_scratch__63__ +
                      acegen_scratch__16__ * acegen_scratch__64__;
  gradientOut[1][0] = acegen_scratch__18__ * acegen_scratch__52__ +
                      acegen_scratch__19__ * acegen_scratch__62__ +
                      acegen_scratch__20__ * acegen_scratch__63__;
  gradientOut[1][1] = acegen_scratch__19__ * acegen_scratch__55__ +
                      acegen_scratch__18__ * acegen_scratch__62__ +
                      acegen_scratch__20__ * acegen_scratch__64__;
  gradientOut[1][2] = acegen_scratch__20__ * acegen_scratch__57__ +
                      acegen_scratch__18__ * acegen_scratch__63__ +
                      acegen_scratch__19__ * acegen_scratch__64__;
  gradientOut[2][0] = acegen_scratch__21__ * acegen_scratch__52__ +
                      acegen_scratch__22__ * acegen_scratch__62__ +
                      acegen_scratch__23__ * acegen_scratch__63__;
  gradientOut[2][1] = acegen_scratch__22__ * acegen_scratch__55__ +
                      acegen_scratch__21__ * acegen_scratch__62__ +
                      acegen_scratch__23__ * acegen_scratch__64__;
  gradientOut[2][2] = acegen_scratch__23__ * acegen_scratch__57__ +
                      acegen_scratch__21__ * acegen_scratch__63__ +
                      acegen_scratch__22__ * acegen_scratch__64__;
  cacheOut[0] =
    4e0 * acegen_scratch__77__ *
    (acegen_scratch__159__ * acegen_scratch__77__ + 2e0 * acegen_scratch__96__);
  cacheOut[1] =
    4e0 * (acegen_scratch__161__ +
           acegen_scratch__101__ * acegen_scratch__77__ + acegen_scratch__91__);
  cacheOut[2] =
    4e0 * (acegen_scratch__162__ +
           acegen_scratch__102__ * acegen_scratch__77__ + acegen_scratch__92__);
  cacheOut[3] =
    4e0 * (-(acegen_scratch__33__ * acegen_scratch__50__) +
           acegen_scratch__103__ * acegen_scratch__77__ + acegen_scratch__93__);
  cacheOut[4] =
    4e0 * (acegen_scratch__104__ * acegen_scratch__77__ + acegen_scratch__94__);
  cacheOut[5] =
    4e0 * (acegen_scratch__105__ * acegen_scratch__77__ + acegen_scratch__95__);
  cacheOut[6] =
    4e0 * (acegen_scratch__101__ * acegen_scratch__78__ + acegen_scratch__91__);
  cacheOut[7] =
    4e0 * (acegen_scratch__30__ * acegen_scratch__50__ +
           acegen_scratch__102__ * acegen_scratch__78__ + acegen_scratch__92__);
  cacheOut[8] =
    4e0 * (acegen_scratch__103__ * acegen_scratch__78__ + acegen_scratch__93__);
  cacheOut[9] =
    4e0 * (-acegen_scratch__164__ +
           acegen_scratch__104__ * acegen_scratch__78__ + acegen_scratch__94__);
  cacheOut[10] =
    4e0 * (acegen_scratch__105__ * acegen_scratch__78__ + acegen_scratch__95__);
  cacheOut[11] =
    4e0 * (acegen_scratch__102__ * acegen_scratch__79__ + acegen_scratch__92__);
  cacheOut[12] =
    4e0 * (acegen_scratch__103__ * acegen_scratch__79__ + acegen_scratch__93__);
  cacheOut[13] =
    4e0 * (acegen_scratch__104__ * acegen_scratch__79__ + acegen_scratch__94__);
  cacheOut[14] =
    4e0 * (-acegen_scratch__165__ +
           acegen_scratch__105__ * acegen_scratch__79__ + acegen_scratch__95__);
  cacheOut[15] = 2e0 * (acegen_scratch__103__ * acegen_scratch__166__ +
                        acegen_scratch__163__ * acegen_scratch__50__);
  cacheOut[16] = 2e0 * (0.1414213562373095e1 * acegen_scratch__165__ +
                        1e0 * acegen_scratch__104__ * acegen_scratch__166__);
  cacheOut[17] = 2e0 * (0.1414213562373095e1 * acegen_scratch__164__ +
                        acegen_scratch__105__ * acegen_scratch__166__);
  cacheOut[18] = 2e0 * (-2e0 * acegen_scratch__162__ +
                        acegen_scratch__104__ * acegen_scratch__167__);
  cacheOut[19] = 2e0 * (acegen_scratch__105__ * acegen_scratch__167__ +
                        2e0 * acegen_scratch__156__ * acegen_scratch__50__);
  cacheOut[20] = 2e0 * (0.1414213562373095e1 * acegen_scratch__105__ *
                          acegen_scratch__128__ -
                        2e0 * acegen_scratch__161__);
  cacheOut[21] = acegen_scratch__52__;
  cacheOut[22] = acegen_scratch__62__;
  cacheOut[23] = acegen_scratch__63__;
  cacheOut[24] = acegen_scratch__55__;
  cacheOut[25] = acegen_scratch__64__;
  cacheOut[26] = acegen_scratch__57__;
}
