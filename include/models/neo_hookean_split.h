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
    acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__20__, acegen_scratch__21__, acegen_scratch__32__,
    acegen_scratch__33__, acegen_scratch__41__, acegen_scratch__42__,
    acegen_scratch__44__, acegen_scratch__46__, acegen_scratch__48__,
    acegen_scratch__49__, acegen_scratch__51__, acegen_scratch__58__,
    acegen_scratch__59__, acegen_scratch__60__, acegen_scratch__7__,
    acegen_scratch__9__;
  acegen_scratch__7__  = (mu);
  acegen_scratch__41__ = (2e0 / 3e0) * acegen_scratch__7__ + (lambda);
  acegen_scratch__60__ = acegen_scratch__41__ / 2e0;
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
  acegen_scratch__58__ = 1e0 / (2e0 * acegen_scratch__20__);
  acegen_scratch__46__ = acegen_scratch__18__ * acegen_scratch__58__;
  acegen_scratch__44__ = acegen_scratch__16__ * acegen_scratch__58__;
  acegen_scratch__42__ = acegen_scratch__17__ * acegen_scratch__58__;
  acegen_scratch__21__ = sqrt(acegen_scratch__20__);
  acegen_scratch__33__ = 1e0 / (acegen_scratch__21__ * acegen_scratch__21__);
  acegen_scratch__59__ =
    -(acegen_scratch__21__ * acegen_scratch__33__ * acegen_scratch__42__);
  acegen_scratch__32__ =
    1e0 / acegen_scratch__21__ + acegen_scratch__16__ * acegen_scratch__59__;
  acegen_scratch__48__ =
    (acegen_scratch__17__ - 2e0 * acegen_scratch__42__) * acegen_scratch__60__ +
    (acegen_scratch__32__ + acegen_scratch__17__ * acegen_scratch__59__) *
      acegen_scratch__7__;
  acegen_scratch__49__ =
    (acegen_scratch__16__ - 2e0 * acegen_scratch__44__) * acegen_scratch__60__ +
    (acegen_scratch__32__ - acegen_scratch__16__ * acegen_scratch__21__ *
                              acegen_scratch__33__ * acegen_scratch__44__) *
      acegen_scratch__7__;
  acegen_scratch__51__ =
    acegen_scratch__41__ * (-0.3535533905932738e0 * acegen_scratch__18__ +
                            0.7071067811865476e0 * acegen_scratch__46__) +
    0.7071067811865476e0 * (acegen_scratch__16__ + acegen_scratch__17__) *
      acegen_scratch__21__ * acegen_scratch__33__ * acegen_scratch__46__ *
      acegen_scratch__7__;
  valueOut[0]       = 0e0;
  valueOut[1]       = 0e0;
  gradientOut[0][0] = acegen_scratch__10__ * acegen_scratch__51__ +
                      acegen_scratch__48__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__10__ * acegen_scratch__49__ +
                      acegen_scratch__51__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__11__ * acegen_scratch__48__ +
                      acegen_scratch__12__ * acegen_scratch__51__;
  gradientOut[1][1] = acegen_scratch__12__ * acegen_scratch__49__ +
                      acegen_scratch__11__ * acegen_scratch__51__;
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
    acegen_scratch__17__, acegen_scratch__18__, acegen_scratch__185__,
    acegen_scratch__186__, acegen_scratch__187__, acegen_scratch__188__,
    acegen_scratch__189__, acegen_scratch__193__, acegen_scratch__194__,
    acegen_scratch__202__, acegen_scratch__204__, acegen_scratch__209__,
    acegen_scratch__210__, acegen_scratch__212__, acegen_scratch__219__,
    acegen_scratch__22__, acegen_scratch__220__, acegen_scratch__221__,
    acegen_scratch__222__, acegen_scratch__223__, acegen_scratch__224__,
    acegen_scratch__225__, acegen_scratch__226__, acegen_scratch__227__,
    acegen_scratch__228__, acegen_scratch__229__, acegen_scratch__23__,
    acegen_scratch__230__, acegen_scratch__24__, acegen_scratch__26__,
    acegen_scratch__27__, acegen_scratch__35__, acegen_scratch__36__,
    acegen_scratch__37__, acegen_scratch__38__, acegen_scratch__39__,
    acegen_scratch__47__, acegen_scratch__48__, acegen_scratch__50__,
    acegen_scratch__52__, acegen_scratch__54__, acegen_scratch__55__,
    acegen_scratch__57__, acegen_scratch__9__;
  acegen_scratch__9__   = gradduIn[0][0];
  acegen_scratch__10__  = gradduIn[0][1];
  acegen_scratch__11__  = gradduIn[1][0];
  acegen_scratch__12__  = gradduIn[1][1];
  acegen_scratch__13__  = (mu);
  acegen_scratch__47__  = (2e0 / 3e0) * acegen_scratch__13__ + (lambda);
  acegen_scratch__226__ = acegen_scratch__47__ / 2e0;
  acegen_scratch__15__  = 1e0 + graduIn[0][0];
  acegen_scratch__16__  = graduIn[0][1];
  acegen_scratch__17__  = graduIn[1][0];
  acegen_scratch__185__ = 2e0 * (acegen_scratch__11__ * acegen_scratch__17__ +
                                 acegen_scratch__15__ * acegen_scratch__9__);
  acegen_scratch__18__  = 1e0 + graduIn[1][1];
  acegen_scratch__219__ = acegen_scratch__15__ * acegen_scratch__16__ +
                          acegen_scratch__17__ * acegen_scratch__18__;
  acegen_scratch__186__ = acegen_scratch__10__ * acegen_scratch__15__ +
                          acegen_scratch__12__ * acegen_scratch__17__ +
                          acegen_scratch__11__ * acegen_scratch__18__ +
                          acegen_scratch__16__ * acegen_scratch__9__;
  acegen_scratch__188__ = 0.1414213562373095e1 * acegen_scratch__186__;
  acegen_scratch__187__ = 2e0 * (acegen_scratch__10__ * acegen_scratch__16__ +
                                 acegen_scratch__12__ * acegen_scratch__18__);
  acegen_scratch__22__  = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__17__ * acegen_scratch__17__);
  acegen_scratch__23__ = (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__223__ = acegen_scratch__22__ + acegen_scratch__23__;
  acegen_scratch__24__  = 0.1414213562373095e1 * acegen_scratch__219__;
  acegen_scratch__189__ = -2e0 * acegen_scratch__186__ * acegen_scratch__219__ +
                          acegen_scratch__187__ * acegen_scratch__22__ +
                          acegen_scratch__185__ * acegen_scratch__23__;
  acegen_scratch__26__ = -(acegen_scratch__219__ * acegen_scratch__219__) +
                         acegen_scratch__22__ * acegen_scratch__23__;
  acegen_scratch__221__ = 1e0 / (2e0 * acegen_scratch__26__);
  acegen_scratch__220__ =
    -(acegen_scratch__189__ / (acegen_scratch__26__ * acegen_scratch__26__));
  acegen_scratch__227__ = acegen_scratch__220__ * acegen_scratch__23__ +
                          acegen_scratch__187__ / acegen_scratch__26__;
  acegen_scratch__224__ = acegen_scratch__22__ * acegen_scratch__220__ +
                          acegen_scratch__185__ / acegen_scratch__26__;
  acegen_scratch__222__ = acegen_scratch__220__ * acegen_scratch__24__ +
                          acegen_scratch__188__ / acegen_scratch__26__;
  acegen_scratch__52__  = acegen_scratch__221__ * acegen_scratch__24__;
  acegen_scratch__50__  = acegen_scratch__22__ * acegen_scratch__221__;
  acegen_scratch__48__  = acegen_scratch__221__ * acegen_scratch__23__;
  acegen_scratch__27__  = sqrt(acegen_scratch__26__);
  acegen_scratch__230__ = acegen_scratch__27__ / 2e0;
  acegen_scratch__193__ =
    acegen_scratch__189__ * acegen_scratch__221__ * acegen_scratch__27__;
  acegen_scratch__202__ = acegen_scratch__227__ * acegen_scratch__230__ +
                          acegen_scratch__193__ * acegen_scratch__48__;
  acegen_scratch__194__ =
    (-2e0 * acegen_scratch__193__) / Power(acegen_scratch__27__, 3);
  acegen_scratch__39__  = 1e0 / (acegen_scratch__27__ * acegen_scratch__27__);
  acegen_scratch__229__ = acegen_scratch__23__ * acegen_scratch__39__;
  acegen_scratch__228__ = acegen_scratch__22__ * acegen_scratch__39__;
  acegen_scratch__225__ = acegen_scratch__194__ * acegen_scratch__22__ +
                          acegen_scratch__185__ * acegen_scratch__39__;
  acegen_scratch__37__ = -(acegen_scratch__27__ * acegen_scratch__52__);
  acegen_scratch__212__ =
    -0.3535533905932738e0 * (acegen_scratch__188__ - acegen_scratch__222__) *
      acegen_scratch__47__ -
    0.7071067811865476e0 * acegen_scratch__13__ *
      (acegen_scratch__37__ * (acegen_scratch__194__ * acegen_scratch__223__ +
                               (acegen_scratch__185__ + acegen_scratch__187__) *
                                 acegen_scratch__39__) +
       acegen_scratch__223__ * acegen_scratch__39__ *
         (-(acegen_scratch__222__ * acegen_scratch__230__) -
          acegen_scratch__193__ * acegen_scratch__52__));
  acegen_scratch__36__ = acegen_scratch__27__ * acegen_scratch__50__;
  acegen_scratch__35__ = acegen_scratch__27__ * acegen_scratch__48__;
  acegen_scratch__204__ =
    -(acegen_scratch__225__ * acegen_scratch__35__) -
    (acegen_scratch__193__ + acegen_scratch__202__ * acegen_scratch__22__) *
      acegen_scratch__39__;
  acegen_scratch__210__ =
    (acegen_scratch__185__ - acegen_scratch__224__) * acegen_scratch__226__ -
    acegen_scratch__13__ *
      (-acegen_scratch__204__ + acegen_scratch__225__ * acegen_scratch__36__ +
       acegen_scratch__228__ * (acegen_scratch__224__ * acegen_scratch__230__ +
                                acegen_scratch__193__ * acegen_scratch__50__));
  acegen_scratch__209__ =
    acegen_scratch__226__ * (acegen_scratch__187__ - acegen_scratch__227__) -
    acegen_scratch__13__ *
      (-acegen_scratch__204__ + acegen_scratch__202__ * acegen_scratch__229__ +
       acegen_scratch__35__ * (acegen_scratch__194__ * acegen_scratch__23__ +
                               acegen_scratch__187__ * acegen_scratch__39__));
  acegen_scratch__38__ =
    1e0 / acegen_scratch__27__ - acegen_scratch__228__ * acegen_scratch__35__;
  acegen_scratch__54__ =
    acegen_scratch__13__ *
      (-(acegen_scratch__229__ * acegen_scratch__35__) + acegen_scratch__38__) +
    acegen_scratch__226__ * (acegen_scratch__23__ - 2e0 * acegen_scratch__48__);
  acegen_scratch__55__ =
    acegen_scratch__13__ *
      (-(acegen_scratch__228__ * acegen_scratch__36__) + acegen_scratch__38__) +
    acegen_scratch__226__ * (acegen_scratch__22__ - 2e0 * acegen_scratch__50__);
  acegen_scratch__57__ =
    -0.7071067811865476e0 * acegen_scratch__13__ *
      (acegen_scratch__22__ + acegen_scratch__23__) * acegen_scratch__37__ *
      acegen_scratch__39__ +
    acegen_scratch__47__ * (-0.3535533905932738e0 * acegen_scratch__24__ +
                            0.7071067811865476e0 * acegen_scratch__52__);
  gradientOut[0][0] = acegen_scratch__15__ * acegen_scratch__209__ +
                      acegen_scratch__16__ * acegen_scratch__212__ +
                      acegen_scratch__10__ * acegen_scratch__57__ +
                      acegen_scratch__54__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__210__ +
                      acegen_scratch__15__ * acegen_scratch__212__ +
                      acegen_scratch__10__ * acegen_scratch__55__ +
                      acegen_scratch__57__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__17__ * acegen_scratch__209__ +
                      acegen_scratch__18__ * acegen_scratch__212__ +
                      acegen_scratch__11__ * acegen_scratch__54__ +
                      acegen_scratch__12__ * acegen_scratch__57__;
  gradientOut[1][1] = acegen_scratch__18__ * acegen_scratch__210__ +
                      acegen_scratch__17__ * acegen_scratch__212__ +
                      acegen_scratch__12__ * acegen_scratch__55__ +
                      acegen_scratch__11__ * acegen_scratch__57__;
}



template <>
template <typename Number>
inline void
SolidModel<2>::cache(const Tensor<2, dim, Number> &           graduIn,
                     [[maybe_unused]] Tensor<1, dim, Number> &valueOut,
                     [[maybe_unused]] Tensor<2, dim, Number> &gradientOut,
                     ArrayView<Number> &                      cacheOut,
                     const Number &                           mu,
                     const Number &                           lambda)
{
  Number acegen_scratch__10__, acegen_scratch__11__, acegen_scratch__12__,
    acegen_scratch__120__, acegen_scratch__121__, acegen_scratch__122__,
    acegen_scratch__123__, acegen_scratch__124__, acegen_scratch__125__,
    acegen_scratch__126__, acegen_scratch__127__, acegen_scratch__128__,
    acegen_scratch__129__, acegen_scratch__130__, acegen_scratch__131__,
    acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__20__, acegen_scratch__21__, acegen_scratch__30__,
    acegen_scratch__32__, acegen_scratch__33__, acegen_scratch__41__,
    acegen_scratch__42__, acegen_scratch__44__, acegen_scratch__46__,
    acegen_scratch__48__, acegen_scratch__49__, acegen_scratch__51__,
    acegen_scratch__60__, acegen_scratch__61__, acegen_scratch__62__,
    acegen_scratch__66__, acegen_scratch__68__, acegen_scratch__69__,
    acegen_scratch__7__, acegen_scratch__70__, acegen_scratch__71__,
    acegen_scratch__77__, acegen_scratch__78__, acegen_scratch__88__,
    acegen_scratch__89__, acegen_scratch__9__, acegen_scratch__91__,
    acegen_scratch__92__, acegen_scratch__93__, acegen_scratch__97__,
    acegen_scratch__98__;
  acegen_scratch__7__   = (mu);
  acegen_scratch__127__ = acegen_scratch__7__ / 2e0;
  acegen_scratch__41__  = (2e0 / 3e0) * acegen_scratch__7__ + (lambda);
  acegen_scratch__124__ = acegen_scratch__41__ / 2e0;
  acegen_scratch__9__   = 1e0 + graduIn[0][0];
  acegen_scratch__10__  = graduIn[0][1];
  acegen_scratch__11__  = graduIn[1][0];
  acegen_scratch__12__  = 1e0 + graduIn[1][1];
  acegen_scratch__16__  = (acegen_scratch__11__ * acegen_scratch__11__) +
                         (acegen_scratch__9__ * acegen_scratch__9__);
  acegen_scratch__17__ = (acegen_scratch__10__ * acegen_scratch__10__) +
                         (acegen_scratch__12__ * acegen_scratch__12__);
  acegen_scratch__130__ = acegen_scratch__16__ * acegen_scratch__17__;
  acegen_scratch__128__ = -0.15e1 * acegen_scratch__17__;
  acegen_scratch__18__ =
    0.1414213562373095e1 * (acegen_scratch__11__ * acegen_scratch__12__ +
                            acegen_scratch__10__ * acegen_scratch__9__);
  acegen_scratch__120__ = 0.5e0 * (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__20__  = -acegen_scratch__120__ + acegen_scratch__130__;
  acegen_scratch__62__  = 1e0 / (2e0 * acegen_scratch__20__);
  acegen_scratch__60__  = 1e0 / (acegen_scratch__20__ * acegen_scratch__20__);
  acegen_scratch__129__ = acegen_scratch__60__ / 2e0;
  acegen_scratch__126__ = (acegen_scratch__41__ * acegen_scratch__60__) / 4e0;
  acegen_scratch__61__ =
    acegen_scratch__120__ * acegen_scratch__60__ + acegen_scratch__62__;
  acegen_scratch__131__ =
    -(acegen_scratch__124__ * acegen_scratch__129__ * acegen_scratch__18__);
  acegen_scratch__46__  = acegen_scratch__18__ * acegen_scratch__62__;
  acegen_scratch__44__  = acegen_scratch__16__ * acegen_scratch__62__;
  acegen_scratch__42__  = acegen_scratch__17__ * acegen_scratch__62__;
  acegen_scratch__21__  = sqrt(acegen_scratch__20__);
  acegen_scratch__70__  = 1e0 / Power(acegen_scratch__21__, 3);
  acegen_scratch__121__ = -2e0 * acegen_scratch__70__;
  acegen_scratch__68__  = -(acegen_scratch__21__ * acegen_scratch__46__);
  acegen_scratch__88__  = -(acegen_scratch__46__ * acegen_scratch__68__);
  acegen_scratch__71__ =
    -((acegen_scratch__16__ * acegen_scratch__21__ * acegen_scratch__70__) /
      acegen_scratch__20__);
  acegen_scratch__93__ = -0.15e1 * acegen_scratch__16__ * acegen_scratch__71__;
  acegen_scratch__66__ = acegen_scratch__21__ * acegen_scratch__42__;
  acegen_scratch__89__ =
    acegen_scratch__44__ * acegen_scratch__66__ - acegen_scratch__88__;
  acegen_scratch__69__ = acegen_scratch__121__ * acegen_scratch__66__;
  acegen_scratch__33__ = 1e0 / (acegen_scratch__21__ * acegen_scratch__21__);
  acegen_scratch__122__ =
    -(acegen_scratch__33__ *
      (-(acegen_scratch__21__ * acegen_scratch__61__) + acegen_scratch__88__));
  acegen_scratch__97__  = acegen_scratch__122__ * acegen_scratch__17__;
  acegen_scratch__91__  = acegen_scratch__122__ * acegen_scratch__16__;
  acegen_scratch__125__ = acegen_scratch__33__ * acegen_scratch__68__;
  acegen_scratch__77__  = acegen_scratch__128__ * acegen_scratch__68__;
  acegen_scratch__78__ =
    -acegen_scratch__125__ + acegen_scratch__71__ * acegen_scratch__77__;
  acegen_scratch__30__ = acegen_scratch__21__ * acegen_scratch__44__;
  acegen_scratch__92__ =
    acegen_scratch__33__ *
      (-acegen_scratch__30__ - acegen_scratch__16__ * acegen_scratch__89__) +
    acegen_scratch__91__;
  acegen_scratch__123__ = -(acegen_scratch__33__ * acegen_scratch__66__);
  acegen_scratch__98__ =
    acegen_scratch__33__ *
      (-acegen_scratch__66__ - acegen_scratch__17__ * acegen_scratch__89__) +
    acegen_scratch__97__;
  acegen_scratch__32__ =
    acegen_scratch__123__ * acegen_scratch__16__ + 1e0 / acegen_scratch__21__;
  acegen_scratch__48__ =
    acegen_scratch__124__ *
      (acegen_scratch__17__ - 2e0 * acegen_scratch__42__) +
    (acegen_scratch__123__ * acegen_scratch__17__ + acegen_scratch__32__) *
      acegen_scratch__7__;
  acegen_scratch__49__ =
    acegen_scratch__124__ *
      (acegen_scratch__16__ - 2e0 * acegen_scratch__44__) +
    (acegen_scratch__32__ -
     acegen_scratch__16__ * acegen_scratch__30__ * acegen_scratch__33__) *
      acegen_scratch__7__;
  acegen_scratch__51__ =
    acegen_scratch__41__ * (-0.3535533905932738e0 * acegen_scratch__18__ +
                            0.7071067811865476e0 * acegen_scratch__46__) -
    0.7071067811865476e0 * acegen_scratch__125__ *
      (acegen_scratch__16__ + acegen_scratch__17__) * acegen_scratch__7__;
  gradientOut[0][0] = acegen_scratch__10__ * acegen_scratch__51__ +
                      acegen_scratch__48__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__10__ * acegen_scratch__49__ +
                      acegen_scratch__51__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__11__ * acegen_scratch__48__ +
                      acegen_scratch__12__ * acegen_scratch__51__;
  gradientOut[1][1] = acegen_scratch__12__ * acegen_scratch__49__ +
                      acegen_scratch__11__ * acegen_scratch__51__;
  cacheOut[0] =
    4e0 *
    (acegen_scratch__126__ * (acegen_scratch__17__ * acegen_scratch__17__) +
     acegen_scratch__127__ *
       (acegen_scratch__128__ * acegen_scratch__66__ * acegen_scratch__69__ +
        acegen_scratch__98__));
  cacheOut[1] =
    4e0 *
    (acegen_scratch__124__ *
       (0.5e0 + acegen_scratch__129__ * acegen_scratch__130__ -
        acegen_scratch__62__) +
     acegen_scratch__127__ * (acegen_scratch__92__ + acegen_scratch__98__));
  cacheOut[2] =
    4e0 *
    (acegen_scratch__131__ * acegen_scratch__17__ +
     acegen_scratch__127__ *
       (acegen_scratch__69__ * acegen_scratch__77__ + acegen_scratch__78__));
  cacheOut[3] =
    4e0 *
    (acegen_scratch__126__ * (acegen_scratch__16__ * acegen_scratch__16__) +
     acegen_scratch__127__ *
       (acegen_scratch__92__ + acegen_scratch__30__ * acegen_scratch__93__));
  cacheOut[4] =
    4e0 *
    (acegen_scratch__131__ * acegen_scratch__16__ +
     acegen_scratch__127__ *
       (acegen_scratch__78__ + acegen_scratch__68__ * acegen_scratch__93__));
  cacheOut[5] = 4e0 * (acegen_scratch__124__ * (-0.5e0 + acegen_scratch__61__) +
                       acegen_scratch__127__ *
                         (acegen_scratch__121__ *
                            (-acegen_scratch__16__ - acegen_scratch__17__) *
                            (acegen_scratch__68__ * acegen_scratch__68__) +
                          acegen_scratch__91__ + acegen_scratch__97__));
  cacheOut[6] = acegen_scratch__48__;
  cacheOut[7] = acegen_scratch__51__;
  cacheOut[8] = acegen_scratch__49__;
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
SolidModel<3>::cache(const Tensor<2, dim, Number> &           graduIn,
                     [[maybe_unused]] Tensor<1, dim, Number> &valueOut,
                     [[maybe_unused]] Tensor<2, dim, Number> &gradientOut,
                     ArrayView<Number> &                      cacheOut,
                     const Number &                           mu,
                     const Number &                           lambda)
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



template <>
template <typename Number>
inline void
SolidModel<2>::cache_deformed(const Tensor<2, dim, Number> &graduIn,
                              Tensor<2, dim, Number> &      gradientOut,
                              ArrayView<Number> &           cacheOut,
                              const Number &                mu,
                              const Number &                lambda)
{
  Number acegen_scratch__10__, acegen_scratch__109__, acegen_scratch__110__,
    acegen_scratch__111__, acegen_scratch__112__, acegen_scratch__113__,
    acegen_scratch__114__, acegen_scratch__115__, acegen_scratch__116__,
    acegen_scratch__119__, acegen_scratch__120__, acegen_scratch__121__,
    acegen_scratch__122__, acegen_scratch__123__, acegen_scratch__127__,
    acegen_scratch__14__, acegen_scratch__142__, acegen_scratch__143__,
    acegen_scratch__15__, acegen_scratch__153__, acegen_scratch__154__,
    acegen_scratch__155__, acegen_scratch__156__, acegen_scratch__157__,
    acegen_scratch__158__, acegen_scratch__159__, acegen_scratch__16__,
    acegen_scratch__160__, acegen_scratch__161__, acegen_scratch__163__,
    acegen_scratch__164__, acegen_scratch__166__, acegen_scratch__167__,
    acegen_scratch__169__, acegen_scratch__170__, acegen_scratch__171__,
    acegen_scratch__172__, acegen_scratch__173__, acegen_scratch__18__,
    acegen_scratch__19__, acegen_scratch__28__, acegen_scratch__30__,
    acegen_scratch__31__, acegen_scratch__39__, acegen_scratch__40__,
    acegen_scratch__42__, acegen_scratch__44__, acegen_scratch__46__,
    acegen_scratch__47__, acegen_scratch__49__, acegen_scratch__5__,
    acegen_scratch__50__, acegen_scratch__51__, acegen_scratch__52__,
    acegen_scratch__53__, acegen_scratch__58__, acegen_scratch__59__,
    acegen_scratch__60__, acegen_scratch__64__, acegen_scratch__66__,
    acegen_scratch__67__, acegen_scratch__68__, acegen_scratch__69__,
    acegen_scratch__7__, acegen_scratch__70__, acegen_scratch__75__,
    acegen_scratch__76__, acegen_scratch__8__, acegen_scratch__86__,
    acegen_scratch__87__, acegen_scratch__89__, acegen_scratch__9__,
    acegen_scratch__90__, acegen_scratch__91__, acegen_scratch__95__,
    acegen_scratch__96__;
  acegen_scratch__5__   = (mu);
  acegen_scratch__161__ = 2e0 * acegen_scratch__5__;
  acegen_scratch__39__  = (2e0 / 3e0) * acegen_scratch__5__ + (lambda);
  acegen_scratch__157__ = acegen_scratch__39__ / 2e0;
  acegen_scratch__7__   = 1e0 + graduIn[0][0];
  acegen_scratch__115__ = (acegen_scratch__7__ * acegen_scratch__7__);
  acegen_scratch__8__   = graduIn[0][1];
  acegen_scratch__116__ = (acegen_scratch__8__ * acegen_scratch__8__);
  acegen_scratch__9__   = graduIn[1][0];
  acegen_scratch__166__ = acegen_scratch__7__ * acegen_scratch__9__;
  acegen_scratch__120__ = (acegen_scratch__9__ * acegen_scratch__9__);
  acegen_scratch__170__ = acegen_scratch__116__ * acegen_scratch__120__;
  acegen_scratch__10__  = 1e0 + graduIn[1][1];
  acegen_scratch__167__ = acegen_scratch__10__ * acegen_scratch__9__;
  acegen_scratch__143__ = 2e0 * acegen_scratch__10__ * acegen_scratch__8__;
  acegen_scratch__121__ = (acegen_scratch__10__ * acegen_scratch__10__);
  acegen_scratch__171__ = acegen_scratch__115__ * acegen_scratch__121__;
  acegen_scratch__14__  = acegen_scratch__115__ + acegen_scratch__120__;
  acegen_scratch__15__  = acegen_scratch__116__ + acegen_scratch__121__;
  acegen_scratch__159__ = acegen_scratch__14__ * acegen_scratch__15__;
  acegen_scratch__16__ =
    0.1414213562373095e1 *
    (acegen_scratch__167__ + acegen_scratch__7__ * acegen_scratch__8__);
  acegen_scratch__153__ = 0.5e0 * (acegen_scratch__16__ * acegen_scratch__16__);
  acegen_scratch__18__  = -acegen_scratch__153__ + acegen_scratch__159__;
  acegen_scratch__119__ = 1e0 / sqrt(acegen_scratch__18__);
  acegen_scratch__60__  = 1e0 / (2e0 * acegen_scratch__18__);
  acegen_scratch__160__ = 2e0 * acegen_scratch__60__;
  acegen_scratch__58__  = 1e0 / (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__163__ = acegen_scratch__39__ * acegen_scratch__58__;
  acegen_scratch__59__ =
    acegen_scratch__153__ * acegen_scratch__58__ + acegen_scratch__60__;
  acegen_scratch__164__ =
    -0.5e0 *
    (acegen_scratch__16__ * acegen_scratch__39__ * acegen_scratch__58__);
  acegen_scratch__44__  = acegen_scratch__16__ * acegen_scratch__60__;
  acegen_scratch__42__  = acegen_scratch__14__ * acegen_scratch__60__;
  acegen_scratch__40__  = acegen_scratch__15__ * acegen_scratch__60__;
  acegen_scratch__19__  = sqrt(acegen_scratch__18__);
  acegen_scratch__68__  = 1e0 / Power(acegen_scratch__19__, 3);
  acegen_scratch__154__ = -2e0 * acegen_scratch__68__;
  acegen_scratch__66__  = -(acegen_scratch__19__ * acegen_scratch__44__);
  acegen_scratch__86__  = -(acegen_scratch__44__ * acegen_scratch__66__);
  acegen_scratch__70__  = acegen_scratch__154__ * acegen_scratch__66__;
  acegen_scratch__69__  = -(acegen_scratch__14__ * acegen_scratch__160__ *
                           acegen_scratch__19__ * acegen_scratch__68__);
  acegen_scratch__91__  = -0.15e1 * acegen_scratch__14__ * acegen_scratch__69__;
  acegen_scratch__64__  = acegen_scratch__19__ * acegen_scratch__40__;
  acegen_scratch__87__ =
    acegen_scratch__42__ * acegen_scratch__64__ - acegen_scratch__86__;
  acegen_scratch__67__ = acegen_scratch__154__ * acegen_scratch__64__;
  acegen_scratch__31__ = 1e0 / (acegen_scratch__19__ * acegen_scratch__19__);
  acegen_scratch__155__ =
    -(acegen_scratch__31__ *
      (-(acegen_scratch__19__ * acegen_scratch__59__) + acegen_scratch__86__));
  acegen_scratch__95__  = acegen_scratch__15__ * acegen_scratch__155__;
  acegen_scratch__89__  = acegen_scratch__14__ * acegen_scratch__155__;
  acegen_scratch__158__ = acegen_scratch__31__ * acegen_scratch__66__;
  acegen_scratch__75__  = -0.15e1 * acegen_scratch__15__ * acegen_scratch__66__;
  acegen_scratch__76__ =
    -acegen_scratch__158__ + acegen_scratch__69__ * acegen_scratch__75__;
  acegen_scratch__28__ = acegen_scratch__19__ * acegen_scratch__42__;
  acegen_scratch__90__ =
    acegen_scratch__31__ *
      (-acegen_scratch__28__ - acegen_scratch__14__ * acegen_scratch__87__) +
    acegen_scratch__89__;
  acegen_scratch__156__ = -(acegen_scratch__31__ * acegen_scratch__64__);
  acegen_scratch__96__ =
    acegen_scratch__31__ *
      (-acegen_scratch__64__ - acegen_scratch__15__ * acegen_scratch__87__) +
    acegen_scratch__95__;
  acegen_scratch__30__ =
    acegen_scratch__14__ * acegen_scratch__156__ + 1e0 / acegen_scratch__19__;
  acegen_scratch__46__ =
    acegen_scratch__157__ *
      (acegen_scratch__15__ - 2e0 * acegen_scratch__40__) +
    (acegen_scratch__15__ * acegen_scratch__156__ + acegen_scratch__30__) *
      acegen_scratch__5__;
  acegen_scratch__47__ =
    acegen_scratch__157__ *
      (acegen_scratch__14__ - 2e0 * acegen_scratch__42__) +
    (acegen_scratch__30__ -
     acegen_scratch__14__ * acegen_scratch__28__ * acegen_scratch__31__) *
      acegen_scratch__5__;
  acegen_scratch__49__ =
    acegen_scratch__39__ * (-0.3535533905932738e0 * acegen_scratch__16__ +
                            0.7071067811865476e0 * acegen_scratch__44__) -
    0.7071067811865476e0 * (acegen_scratch__14__ + acegen_scratch__15__) *
      acegen_scratch__158__ * acegen_scratch__5__;
  acegen_scratch__50__ = acegen_scratch__46__ * acegen_scratch__7__ +
                         acegen_scratch__49__ * acegen_scratch__8__;
  acegen_scratch__51__ = acegen_scratch__49__ * acegen_scratch__7__ +
                         acegen_scratch__47__ * acegen_scratch__8__;
  acegen_scratch__173__ = -(acegen_scratch__10__ * acegen_scratch__51__) -
                          acegen_scratch__50__ * acegen_scratch__9__;
  acegen_scratch__52__ = acegen_scratch__10__ * acegen_scratch__49__ +
                         acegen_scratch__46__ * acegen_scratch__9__;
  acegen_scratch__53__ = acegen_scratch__10__ * acegen_scratch__47__ +
                         acegen_scratch__49__ * acegen_scratch__9__;
  acegen_scratch__109__ =
    (acegen_scratch__15__ * acegen_scratch__15__) * acegen_scratch__163__ -
    3e0 * acegen_scratch__15__ * acegen_scratch__5__ * acegen_scratch__64__ *
      acegen_scratch__67__ +
    acegen_scratch__161__ * acegen_scratch__96__;
  acegen_scratch__110__ =
    acegen_scratch__39__ * (1e0 - acegen_scratch__160__ +
                            acegen_scratch__159__ * acegen_scratch__58__) +
    acegen_scratch__161__ * (acegen_scratch__90__ + acegen_scratch__96__);
  acegen_scratch__111__ =
    2e0 * (acegen_scratch__15__ * acegen_scratch__164__ +
           acegen_scratch__5__ * (acegen_scratch__67__ * acegen_scratch__75__ +
                                  acegen_scratch__76__));
  acegen_scratch__123__ =
    0.1414213562373095e1 * acegen_scratch__111__ * acegen_scratch__8__;
  acegen_scratch__112__ =
    (acegen_scratch__14__ * acegen_scratch__14__) * acegen_scratch__163__ +
    acegen_scratch__161__ *
      (acegen_scratch__90__ + acegen_scratch__28__ * acegen_scratch__91__);
  acegen_scratch__113__ =
    2e0 * (acegen_scratch__14__ * acegen_scratch__164__ +
           acegen_scratch__5__ * (acegen_scratch__76__ +
                                  acegen_scratch__66__ * acegen_scratch__91__));
  acegen_scratch__122__ =
    0.1414213562373095e1 * acegen_scratch__113__ * acegen_scratch__7__;
  acegen_scratch__114__ =
    acegen_scratch__39__ * (-1e0 + 2e0 * acegen_scratch__59__) +
    acegen_scratch__161__ *
      (acegen_scratch__70__ * (-(acegen_scratch__14__ * acegen_scratch__66__) +
                               (2e0 / 3e0) * acegen_scratch__75__) +
       acegen_scratch__89__ + acegen_scratch__95__);
  acegen_scratch__169__ = acegen_scratch__110__ + acegen_scratch__114__;
  acegen_scratch__127__ = acegen_scratch__114__ * acegen_scratch__166__;
  acegen_scratch__142__ =
    0.1414213562373095e1 *
      (acegen_scratch__111__ * acegen_scratch__115__ +
       acegen_scratch__113__ * acegen_scratch__116__) *
      acegen_scratch__167__ +
    acegen_scratch__120__ * (acegen_scratch__109__ * acegen_scratch__115__ +
                             acegen_scratch__123__ * acegen_scratch__7__) +
    acegen_scratch__121__ * (acegen_scratch__112__ * acegen_scratch__116__ +
                             acegen_scratch__122__ * acegen_scratch__8__);
  acegen_scratch__172__ =
    0.1414213562373095e1 * acegen_scratch__119__ * acegen_scratch__173__;
  gradientOut[0][0] = acegen_scratch__50__;
  gradientOut[0][1] = acegen_scratch__51__;
  gradientOut[1][0] = acegen_scratch__52__;
  gradientOut[1][1] = acegen_scratch__53__;
  cacheOut[0]       = acegen_scratch__119__ *
                (2e0 * acegen_scratch__115__ * acegen_scratch__116__ *
                   acegen_scratch__169__ +
                 2e0 * acegen_scratch__123__ * Power(acegen_scratch__7__, 3) +
                 acegen_scratch__109__ * Power(acegen_scratch__7__, 4) +
                 2e0 * acegen_scratch__122__ * Power(acegen_scratch__8__, 3) +
                 acegen_scratch__112__ * Power(acegen_scratch__8__, 4));
  cacheOut[1] =
    acegen_scratch__119__ *
    (acegen_scratch__142__ + acegen_scratch__127__ * acegen_scratch__143__ +
     acegen_scratch__110__ * (acegen_scratch__170__ + acegen_scratch__171__));
  cacheOut[2] = acegen_scratch__172__;
  cacheOut[3] = acegen_scratch__119__ *
                (Power(acegen_scratch__10__, 4) * acegen_scratch__112__ +
                 2e0 * acegen_scratch__120__ * acegen_scratch__121__ *
                   acegen_scratch__169__ +
                 0.282842712474619e1 * Power(acegen_scratch__10__, 3) *
                   acegen_scratch__113__ * acegen_scratch__9__ +
                 0.282842712474619e1 * acegen_scratch__10__ *
                   acegen_scratch__111__ * Power(acegen_scratch__9__, 3) +
                 acegen_scratch__109__ * Power(acegen_scratch__9__, 4));
  cacheOut[4] = acegen_scratch__172__;
  cacheOut[5] =
    2e0 * acegen_scratch__119__ *
    (acegen_scratch__142__ +
     acegen_scratch__143__ * (acegen_scratch__127__ / 2e0 +
                              acegen_scratch__110__ * acegen_scratch__166__) +
     (acegen_scratch__114__ * (acegen_scratch__170__ + acegen_scratch__171__)) /
       2e0);
  cacheOut[6] =
    acegen_scratch__119__ * (acegen_scratch__50__ * acegen_scratch__7__ +
                             acegen_scratch__51__ * acegen_scratch__8__);
  cacheOut[7] = -(acegen_scratch__119__ * acegen_scratch__173__);
  cacheOut[8] =
    acegen_scratch__119__ * (acegen_scratch__10__ * acegen_scratch__53__ +
                             acegen_scratch__52__ * acegen_scratch__9__);
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
  Number acegen_scratch__10__, acegen_scratch__104__, acegen_scratch__105__,
    acegen_scratch__107__, acegen_scratch__109__, acegen_scratch__110__,
    acegen_scratch__111__, acegen_scratch__112__, acegen_scratch__113__,
    acegen_scratch__114__, acegen_scratch__115__, acegen_scratch__12__,
    acegen_scratch__13__, acegen_scratch__131__, acegen_scratch__135__,
    acegen_scratch__138__, acegen_scratch__14__, acegen_scratch__140__,
    acegen_scratch__141__, acegen_scratch__142__, acegen_scratch__143__,
    acegen_scratch__144__, acegen_scratch__145__, acegen_scratch__146__,
    acegen_scratch__147__, acegen_scratch__148__, acegen_scratch__149__,
    acegen_scratch__15__, acegen_scratch__150__, acegen_scratch__151__,
    acegen_scratch__152__, acegen_scratch__153__, acegen_scratch__154__,
    acegen_scratch__155__, acegen_scratch__156__, acegen_scratch__157__,
    acegen_scratch__158__, acegen_scratch__159__, acegen_scratch__16__,
    acegen_scratch__160__, acegen_scratch__161__, acegen_scratch__162__,
    acegen_scratch__164__, acegen_scratch__168__, acegen_scratch__169__,
    acegen_scratch__17__, acegen_scratch__171__, acegen_scratch__172__,
    acegen_scratch__173__, acegen_scratch__175__, acegen_scratch__176__,
    acegen_scratch__177__, acegen_scratch__178__, acegen_scratch__179__,
    acegen_scratch__18__, acegen_scratch__180__, acegen_scratch__182__,
    acegen_scratch__183__, acegen_scratch__184__, acegen_scratch__186__,
    acegen_scratch__188__, acegen_scratch__189__, acegen_scratch__19__,
    acegen_scratch__191__, acegen_scratch__192__, acegen_scratch__193__,
    acegen_scratch__194__, acegen_scratch__195__, acegen_scratch__196__,
    acegen_scratch__197__, acegen_scratch__198__, acegen_scratch__199__,
    acegen_scratch__20__, acegen_scratch__203__, acegen_scratch__204__,
    acegen_scratch__209__, acegen_scratch__210__, acegen_scratch__211__,
    acegen_scratch__212__, acegen_scratch__213__, acegen_scratch__214__,
    acegen_scratch__215__, acegen_scratch__216__, acegen_scratch__217__,
    acegen_scratch__218__, acegen_scratch__219__, acegen_scratch__220__,
    acegen_scratch__221__, acegen_scratch__222__, acegen_scratch__223__,
    acegen_scratch__228__, acegen_scratch__236__, acegen_scratch__239__,
    acegen_scratch__245__, acegen_scratch__260__, acegen_scratch__269__,
    acegen_scratch__27__, acegen_scratch__270__, acegen_scratch__271__,
    acegen_scratch__273__, acegen_scratch__274__, acegen_scratch__275__,
    acegen_scratch__277__, acegen_scratch__278__, acegen_scratch__28__,
    acegen_scratch__280__, acegen_scratch__281__, acegen_scratch__282__,
    acegen_scratch__284__, acegen_scratch__287__, acegen_scratch__29__,
    acegen_scratch__291__, acegen_scratch__292__, acegen_scratch__294__,
    acegen_scratch__295__, acegen_scratch__30__, acegen_scratch__301__,
    acegen_scratch__303__, acegen_scratch__304__, acegen_scratch__306__,
    acegen_scratch__308__, acegen_scratch__31__, acegen_scratch__318__,
    acegen_scratch__319__, acegen_scratch__32__, acegen_scratch__321__,
    acegen_scratch__329__, acegen_scratch__338__, acegen_scratch__339__,
    acegen_scratch__340__, acegen_scratch__342__, acegen_scratch__343__,
    acegen_scratch__344__, acegen_scratch__350__, acegen_scratch__351__,
    acegen_scratch__355__, acegen_scratch__358__, acegen_scratch__36__,
    acegen_scratch__360__, acegen_scratch__363__, acegen_scratch__364__,
    acegen_scratch__365__, acegen_scratch__37__, acegen_scratch__377__,
    acegen_scratch__378__, acegen_scratch__380__, acegen_scratch__381__,
    acegen_scratch__382__, acegen_scratch__387__, acegen_scratch__388__,
    acegen_scratch__481__, acegen_scratch__482__, acegen_scratch__483__,
    acegen_scratch__484__, acegen_scratch__485__, acegen_scratch__486__,
    acegen_scratch__487__, acegen_scratch__488__, acegen_scratch__489__,
    acegen_scratch__49__, acegen_scratch__490__, acegen_scratch__491__,
    acegen_scratch__492__, acegen_scratch__493__, acegen_scratch__494__,
    acegen_scratch__495__, acegen_scratch__496__, acegen_scratch__497__,
    acegen_scratch__498__, acegen_scratch__499__, acegen_scratch__50__,
    acegen_scratch__500__, acegen_scratch__501__, acegen_scratch__502__,
    acegen_scratch__503__, acegen_scratch__504__, acegen_scratch__505__,
    acegen_scratch__506__, acegen_scratch__507__, acegen_scratch__508__,
    acegen_scratch__509__, acegen_scratch__51__, acegen_scratch__510__,
    acegen_scratch__511__, acegen_scratch__512__, acegen_scratch__513__,
    acegen_scratch__514__, acegen_scratch__515__, acegen_scratch__516__,
    acegen_scratch__517__, acegen_scratch__519__, acegen_scratch__520__,
    acegen_scratch__521__, acegen_scratch__522__, acegen_scratch__523__,
    acegen_scratch__524__, acegen_scratch__525__, acegen_scratch__526__,
    acegen_scratch__527__, acegen_scratch__528__, acegen_scratch__529__,
    acegen_scratch__53__, acegen_scratch__530__, acegen_scratch__532__,
    acegen_scratch__533__, acegen_scratch__534__, acegen_scratch__535__,
    acegen_scratch__536__, acegen_scratch__537__, acegen_scratch__538__,
    acegen_scratch__539__, acegen_scratch__54__, acegen_scratch__542__,
    acegen_scratch__543__, acegen_scratch__544__, acegen_scratch__545__,
    acegen_scratch__546__, acegen_scratch__547__, acegen_scratch__548__,
    acegen_scratch__549__, acegen_scratch__55__, acegen_scratch__550__,
    acegen_scratch__551__, acegen_scratch__552__, acegen_scratch__553__,
    acegen_scratch__554__, acegen_scratch__555__, acegen_scratch__556__,
    acegen_scratch__557__, acegen_scratch__558__, acegen_scratch__559__,
    acegen_scratch__56__, acegen_scratch__560__, acegen_scratch__561__,
    acegen_scratch__563__, acegen_scratch__564__, acegen_scratch__565__,
    acegen_scratch__567__, acegen_scratch__568__, acegen_scratch__57__,
    acegen_scratch__570__, acegen_scratch__571__, acegen_scratch__574__,
    acegen_scratch__576__, acegen_scratch__577__, acegen_scratch__578__,
    acegen_scratch__579__, acegen_scratch__58__, acegen_scratch__580__,
    acegen_scratch__581__, acegen_scratch__582__, acegen_scratch__583__,
    acegen_scratch__584__, acegen_scratch__585__, acegen_scratch__586__,
    acegen_scratch__587__, acegen_scratch__588__, acegen_scratch__589__,
    acegen_scratch__59__, acegen_scratch__590__, acegen_scratch__591__,
    acegen_scratch__592__, acegen_scratch__593__, acegen_scratch__594__,
    acegen_scratch__595__, acegen_scratch__596__, acegen_scratch__597__,
    acegen_scratch__598__, acegen_scratch__599__, acegen_scratch__600__,
    acegen_scratch__601__, acegen_scratch__602__, acegen_scratch__603__,
    acegen_scratch__604__, acegen_scratch__605__, acegen_scratch__607__,
    acegen_scratch__608__, acegen_scratch__609__, acegen_scratch__611__,
    acegen_scratch__612__, acegen_scratch__613__, acegen_scratch__614__,
    acegen_scratch__615__, acegen_scratch__616__, acegen_scratch__617__,
    acegen_scratch__618__, acegen_scratch__619__, acegen_scratch__620__,
    acegen_scratch__621__, acegen_scratch__623__, acegen_scratch__624__,
    acegen_scratch__625__, acegen_scratch__627__, acegen_scratch__628__,
    acegen_scratch__63__, acegen_scratch__64__, acegen_scratch__65__,
    acegen_scratch__66__, acegen_scratch__67__, acegen_scratch__68__,
    acegen_scratch__69__, acegen_scratch__70__, acegen_scratch__71__,
    acegen_scratch__72__, acegen_scratch__73__, acegen_scratch__74__,
    acegen_scratch__78__, acegen_scratch__79__, acegen_scratch__80__,
    acegen_scratch__81__, acegen_scratch__82__, acegen_scratch__83__,
    acegen_scratch__84__, acegen_scratch__85__, acegen_scratch__86__,
    acegen_scratch__91__, acegen_scratch__92__, acegen_scratch__93__,
    acegen_scratch__94__, acegen_scratch__95__, acegen_scratch__96__;
  acegen_scratch__10__  = (mu);
  acegen_scratch__493__ = acegen_scratch__10__ / 2e0;
  acegen_scratch__49__  = (2e0 / 3e0) * acegen_scratch__10__ + (lambda);
  acegen_scratch__12__  = 1e0 + graduIn[0][0];
  acegen_scratch__498__ = 0.7071067811865476e0 * acegen_scratch__12__;
  acegen_scratch__161__ = (acegen_scratch__12__ * acegen_scratch__12__);
  acegen_scratch__539__ = 2e0 * acegen_scratch__161__;
  acegen_scratch__13__  = graduIn[0][1];
  acegen_scratch__548__ = 2e0 * acegen_scratch__13__;
  acegen_scratch__513__ = acegen_scratch__12__ * acegen_scratch__13__;
  acegen_scratch__162__ = (acegen_scratch__13__ * acegen_scratch__13__);
  acegen_scratch__14__  = graduIn[0][2];
  acegen_scratch__530__ = acegen_scratch__13__ * acegen_scratch__14__;
  acegen_scratch__508__ = acegen_scratch__14__ / 2e0;
  acegen_scratch__507__ = 0.7071067811865476e0 * acegen_scratch__14__;
  acegen_scratch__576__ = 1e0 * acegen_scratch__507__;
  acegen_scratch__189__ = 2e0 * acegen_scratch__12__ * acegen_scratch__14__;
  acegen_scratch__596__ = acegen_scratch__189__ / 2e0;
  acegen_scratch__164__ = (acegen_scratch__14__ * acegen_scratch__14__);
  acegen_scratch__519__ = acegen_scratch__12__ * acegen_scratch__164__;
  acegen_scratch__15__  = graduIn[1][0];
  acegen_scratch__565__ = 0.7071067811865476e0 * acegen_scratch__15__;
  acegen_scratch__538__ = acegen_scratch__12__ * acegen_scratch__15__;
  acegen_scratch__377__ = acegen_scratch__15__ * acegen_scratch__507__;
  acegen_scratch__176__ = (acegen_scratch__15__ * acegen_scratch__15__);
  acegen_scratch__568__ = 2e0 * acegen_scratch__176__;
  acegen_scratch__16__  = 1e0 + graduIn[1][1];
  acegen_scratch__555__ = 0.7071067811865476e0 * acegen_scratch__16__;
  acegen_scratch__582__ = 1e0 * acegen_scratch__555__;
  acegen_scratch__545__ = acegen_scratch__15__ * acegen_scratch__16__;
  acegen_scratch__177__ = (acegen_scratch__16__ * acegen_scratch__16__);
  acegen_scratch__602__ = 2e0 * acegen_scratch__177__;
  acegen_scratch__17__  = graduIn[1][2];
  acegen_scratch__592__ = acegen_scratch__14__ * acegen_scratch__17__;
  acegen_scratch__581__ = acegen_scratch__17__ * acegen_scratch__176__;
  acegen_scratch__570__ = acegen_scratch__17__ / 2e0;
  acegen_scratch__554__ = 0.7071067811865476e0 * acegen_scratch__17__;
  acegen_scratch__516__ = acegen_scratch__16__ * acegen_scratch__17__;
  acegen_scratch__321__ = acegen_scratch__17__ * acegen_scratch__498__;
  acegen_scratch__282__ = 2e0 * acegen_scratch__15__ * acegen_scratch__17__;
  acegen_scratch__574__ = acegen_scratch__282__ / 2e0;
  acegen_scratch__544__ = 0.7071067811865476e0 * acegen_scratch__282__;
  acegen_scratch__178__ = (acegen_scratch__17__ * acegen_scratch__17__);
  acegen_scratch__18__  = graduIn[2][0];
  acegen_scratch__577__ = acegen_scratch__12__ * acegen_scratch__18__;
  acegen_scratch__561__ = acegen_scratch__13__ * acegen_scratch__18__;
  acegen_scratch__525__ = acegen_scratch__14__ * acegen_scratch__18__;
  acegen_scratch__495__ = acegen_scratch__15__ * acegen_scratch__18__;
  acegen_scratch__388__ = acegen_scratch__18__ * acegen_scratch__544__;
  acegen_scratch__191__ = (acegen_scratch__18__ * acegen_scratch__18__);
  acegen_scratch__583__ = 2e0 * acegen_scratch__191__;
  acegen_scratch__19__  = graduIn[2][1];
  acegen_scratch__620__ = acegen_scratch__17__ * acegen_scratch__19__;
  acegen_scratch__571__ = acegen_scratch__16__ * acegen_scratch__19__;
  acegen_scratch__515__ = acegen_scratch__18__ * acegen_scratch__19__;
  acegen_scratch__481__ =
    acegen_scratch__513__ + acegen_scratch__515__ + acegen_scratch__545__;
  acegen_scratch__192__ = (acegen_scratch__19__ * acegen_scratch__19__);
  acegen_scratch__20__  = 1e0 + graduIn[2][2];
  acegen_scratch__589__ = acegen_scratch__14__ * acegen_scratch__20__;
  acegen_scratch__563__ = acegen_scratch__13__ * acegen_scratch__20__;
  acegen_scratch__552__ = 0.7071067811865476e0 * acegen_scratch__20__;
  acegen_scratch__532__ = acegen_scratch__19__ * acegen_scratch__20__;
  acegen_scratch__522__ = acegen_scratch__12__ * acegen_scratch__20__;
  acegen_scratch__505__ = acegen_scratch__17__ * acegen_scratch__20__;
  acegen_scratch__483__ =
    acegen_scratch__516__ + acegen_scratch__530__ + acegen_scratch__532__;
  acegen_scratch__281__ = 2e0 * acegen_scratch__18__ * acegen_scratch__20__;
  acegen_scratch__595__ = acegen_scratch__281__ / 4e0;
  acegen_scratch__497__ = acegen_scratch__281__ / 2e0;
  acegen_scratch__482__ =
    acegen_scratch__497__ + acegen_scratch__574__ + acegen_scratch__596__;
  acegen_scratch__193__ = (acegen_scratch__20__ * acegen_scratch__20__);
  acegen_scratch__608__ = 2e0 * acegen_scratch__193__;
  acegen_scratch__27__ =
    acegen_scratch__161__ + acegen_scratch__176__ + acegen_scratch__191__;
  acegen_scratch__28__ =
    acegen_scratch__162__ + acegen_scratch__177__ + acegen_scratch__192__;
  acegen_scratch__29__ =
    acegen_scratch__164__ + acegen_scratch__178__ + acegen_scratch__193__;
  acegen_scratch__490__ =
    acegen_scratch__27__ + acegen_scratch__28__ + acegen_scratch__29__;
  acegen_scratch__30__  = 0.1414213562373095e1 * acegen_scratch__483__;
  acegen_scratch__31__  = 0.1414213562373095e1 * acegen_scratch__482__;
  acegen_scratch__32__  = 0.1414213562373095e1 * acegen_scratch__481__;
  acegen_scratch__484__ = 2e0 * acegen_scratch__481__;
  acegen_scratch__57__  = (acegen_scratch__481__ * acegen_scratch__481__);
  acegen_scratch__80__ =
    acegen_scratch__27__ * acegen_scratch__28__ - acegen_scratch__57__;
  acegen_scratch__485__ = 2e0 * acegen_scratch__482__;
  acegen_scratch__59__  = acegen_scratch__482__ * acegen_scratch__484__;
  acegen_scratch__81__  = -(acegen_scratch__27__ * acegen_scratch__30__) +
                         0.7071067811865476e0 * acegen_scratch__59__;
  acegen_scratch__55__ = (acegen_scratch__482__ * acegen_scratch__482__);
  acegen_scratch__79__ =
    acegen_scratch__27__ * acegen_scratch__29__ - acegen_scratch__55__;
  acegen_scratch__138__ = -(acegen_scratch__29__ * acegen_scratch__484__) +
                          acegen_scratch__483__ * acegen_scratch__485__;
  acegen_scratch__135__ = acegen_scratch__483__ * acegen_scratch__484__ -
                          acegen_scratch__28__ * acegen_scratch__485__;
  acegen_scratch__529__ = 0.282842712474619e1 * acegen_scratch__135__;
  acegen_scratch__131__ =
    -2e0 * acegen_scratch__27__ * acegen_scratch__483__ + acegen_scratch__59__;
  acegen_scratch__512__ = 0.282842712474619e1 * acegen_scratch__131__;
  acegen_scratch__83__  = -(acegen_scratch__29__ * acegen_scratch__32__) +
                         acegen_scratch__31__ * acegen_scratch__483__;
  acegen_scratch__82__ = -(acegen_scratch__28__ * acegen_scratch__31__) +
                         acegen_scratch__32__ * acegen_scratch__483__;
  acegen_scratch__78__ = acegen_scratch__28__ * acegen_scratch__29__ -
                         (acegen_scratch__483__ * acegen_scratch__483__);
  acegen_scratch__36__ = -(acegen_scratch__28__ * acegen_scratch__55__) -
                         acegen_scratch__29__ * acegen_scratch__57__ +
                         acegen_scratch__483__ * acegen_scratch__59__ +
                         acegen_scratch__27__ * acegen_scratch__78__;
  acegen_scratch__491__ = -2e0 / acegen_scratch__36__;
  acegen_scratch__175__ = 1e0 / sqrt(acegen_scratch__36__);
  acegen_scratch__615__ = 0.1414213562373095e1 * acegen_scratch__175__;
  acegen_scratch__104__ = 1e0 / (acegen_scratch__36__ * acegen_scratch__36__);
  acegen_scratch__107__ = (acegen_scratch__104__ * acegen_scratch__49__) / 4e0;
  acegen_scratch__37__  = CBRT::exact(acegen_scratch__36__);
  acegen_scratch__489__ = acegen_scratch__37__ / (3e0 * acegen_scratch__36__);
  acegen_scratch__486__ = acegen_scratch__10__ * acegen_scratch__37__;
  acegen_scratch__109__ = acegen_scratch__486__ / (6e0 * acegen_scratch__36__);
  acegen_scratch__50__  = 1e0 / (acegen_scratch__37__ * acegen_scratch__37__);
  acegen_scratch__105__ =
    -(acegen_scratch__486__ * acegen_scratch__490__ * acegen_scratch__50__);
  acegen_scratch__113__ =
    (-2e0 / 9e0) * acegen_scratch__104__ * acegen_scratch__105__;
  acegen_scratch__487__ = acegen_scratch__107__ + acegen_scratch__113__;
  acegen_scratch__115__ = acegen_scratch__487__ * acegen_scratch__83__;
  acegen_scratch__114__ = acegen_scratch__487__ * acegen_scratch__82__;
  acegen_scratch__112__ = acegen_scratch__487__ * acegen_scratch__81__;
  acegen_scratch__492__ = acegen_scratch__107__ + acegen_scratch__113__ / 2e0;
  acegen_scratch__91__  = -(acegen_scratch__493__ * acegen_scratch__50__);
  acegen_scratch__488__ = acegen_scratch__489__ * acegen_scratch__91__;
  acegen_scratch__96__  = acegen_scratch__488__ * acegen_scratch__83__;
  acegen_scratch__95__  = acegen_scratch__488__ * acegen_scratch__82__;
  acegen_scratch__94__  = acegen_scratch__488__ * acegen_scratch__81__;
  acegen_scratch__86__  = acegen_scratch__489__ * acegen_scratch__80__;
  acegen_scratch__111__ =
    acegen_scratch__492__ * acegen_scratch__80__ -
    acegen_scratch__109__ *
      (acegen_scratch__50__ +
       acegen_scratch__490__ * acegen_scratch__491__ * acegen_scratch__86__);
  acegen_scratch__93__ = acegen_scratch__86__ * acegen_scratch__91__;
  acegen_scratch__85__ = acegen_scratch__489__ * acegen_scratch__79__;
  acegen_scratch__110__ =
    acegen_scratch__492__ * acegen_scratch__79__ -
    acegen_scratch__109__ *
      (acegen_scratch__50__ +
       acegen_scratch__490__ * acegen_scratch__491__ * acegen_scratch__85__);
  acegen_scratch__92__ = acegen_scratch__85__ * acegen_scratch__91__;
  acegen_scratch__84__ = acegen_scratch__489__ * acegen_scratch__78__;
  acegen_scratch__54__ = acegen_scratch__493__ / acegen_scratch__37__;
  acegen_scratch__51__ =
    (2e0 * acegen_scratch__105__ +
     3e0 * (-1e0 + acegen_scratch__36__) * acegen_scratch__49__) /
    (12e0 * acegen_scratch__36__);
  acegen_scratch__536__ = 4e0 * acegen_scratch__51__;
  acegen_scratch__526__ = acegen_scratch__28__ * acegen_scratch__51__;
  acegen_scratch__521__ = acegen_scratch__31__ * acegen_scratch__51__;
  acegen_scratch__511__ = acegen_scratch__32__ * acegen_scratch__51__;
  acegen_scratch__509__ = acegen_scratch__27__ * acegen_scratch__51__;
  acegen_scratch__494__ = 1e0 * acegen_scratch__51__;
  acegen_scratch__53__ =
    2e0 * (acegen_scratch__54__ + acegen_scratch__51__ * acegen_scratch__78__);
  acegen_scratch__56__ =
    2e0 * (acegen_scratch__54__ + acegen_scratch__51__ * acegen_scratch__79__);
  acegen_scratch__58__ =
    2e0 * (acegen_scratch__54__ + acegen_scratch__51__ * acegen_scratch__80__);
  acegen_scratch__63__ = acegen_scratch__138__ * acegen_scratch__494__;
  acegen_scratch__64__ = acegen_scratch__135__ * acegen_scratch__494__;
  acegen_scratch__65__ = 1e0 * acegen_scratch__131__ * acegen_scratch__494__;
  acegen_scratch__66__ = acegen_scratch__12__ * acegen_scratch__53__ +
                         acegen_scratch__13__ * acegen_scratch__63__ +
                         acegen_scratch__14__ * acegen_scratch__64__;
  acegen_scratch__67__ = acegen_scratch__13__ * acegen_scratch__56__ +
                         acegen_scratch__12__ * acegen_scratch__63__ +
                         acegen_scratch__14__ * acegen_scratch__65__;
  acegen_scratch__68__ = acegen_scratch__14__ * acegen_scratch__58__ +
                         acegen_scratch__12__ * acegen_scratch__64__ +
                         acegen_scratch__13__ * acegen_scratch__65__;
  acegen_scratch__69__ = acegen_scratch__15__ * acegen_scratch__53__ +
                         acegen_scratch__16__ * acegen_scratch__63__ +
                         acegen_scratch__17__ * acegen_scratch__64__;
  acegen_scratch__70__ = acegen_scratch__16__ * acegen_scratch__56__ +
                         acegen_scratch__15__ * acegen_scratch__63__ +
                         acegen_scratch__17__ * acegen_scratch__65__;
  acegen_scratch__71__ = acegen_scratch__17__ * acegen_scratch__58__ +
                         acegen_scratch__15__ * acegen_scratch__64__ +
                         acegen_scratch__16__ * acegen_scratch__65__;
  acegen_scratch__72__ = acegen_scratch__18__ * acegen_scratch__53__ +
                         acegen_scratch__19__ * acegen_scratch__63__ +
                         acegen_scratch__20__ * acegen_scratch__64__;
  acegen_scratch__73__ = acegen_scratch__19__ * acegen_scratch__56__ +
                         acegen_scratch__18__ * acegen_scratch__63__ +
                         acegen_scratch__20__ * acegen_scratch__65__;
  acegen_scratch__74__ = acegen_scratch__20__ * acegen_scratch__58__ +
                         acegen_scratch__18__ * acegen_scratch__64__ +
                         acegen_scratch__19__ * acegen_scratch__65__;
  acegen_scratch__140__ =
    4e0 *
    (acegen_scratch__78__ *
       (acegen_scratch__492__ * acegen_scratch__78__ -
        acegen_scratch__109__ * (acegen_scratch__50__ +
                                 acegen_scratch__490__ * acegen_scratch__491__ *
                                   acegen_scratch__84__)) +
     acegen_scratch__84__ * acegen_scratch__91__);
  acegen_scratch__590__ = acegen_scratch__140__ * acegen_scratch__191__;
  acegen_scratch__141__ =
    4e0 * (acegen_scratch__29__ * acegen_scratch__51__ +
           acegen_scratch__110__ * acegen_scratch__78__ + acegen_scratch__92__);
  acegen_scratch__496__ = acegen_scratch__141__ * acegen_scratch__18__;
  acegen_scratch__142__ =
    4e0 * (acegen_scratch__526__ +
           acegen_scratch__111__ * acegen_scratch__78__ + acegen_scratch__93__);
  acegen_scratch__143__ =
    4e0 * (-(acegen_scratch__30__ * acegen_scratch__51__) +
           acegen_scratch__112__ * acegen_scratch__78__ + acegen_scratch__94__);
  acegen_scratch__144__ =
    4e0 * (acegen_scratch__114__ * acegen_scratch__78__ + acegen_scratch__95__);
  acegen_scratch__553__ = 0.7071067811865476e0 * acegen_scratch__144__;
  acegen_scratch__585__ = 1e0 * acegen_scratch__553__;
  acegen_scratch__499__ = 3e0 * acegen_scratch__144__;
  acegen_scratch__358__ = acegen_scratch__191__ * acegen_scratch__499__;
  acegen_scratch__145__ =
    4e0 * (acegen_scratch__115__ * acegen_scratch__78__ + acegen_scratch__96__);
  acegen_scratch__584__ = 0.7071067811865476e0 * acegen_scratch__145__;
  acegen_scratch__500__ = 3e0 * acegen_scratch__145__;
  acegen_scratch__360__ = acegen_scratch__191__ * acegen_scratch__500__;
  acegen_scratch__146__ =
    4e0 * (acegen_scratch__110__ * acegen_scratch__79__ + acegen_scratch__92__);
  acegen_scratch__591__ = acegen_scratch__146__ * acegen_scratch__192__;
  acegen_scratch__147__ =
    4e0 * (acegen_scratch__509__ +
           acegen_scratch__111__ * acegen_scratch__79__ + acegen_scratch__93__);
  acegen_scratch__501__ = acegen_scratch__14__ * acegen_scratch__147__;
  acegen_scratch__148__ =
    4e0 * (acegen_scratch__112__ * acegen_scratch__79__ + acegen_scratch__94__);
  acegen_scratch__502__ = 3e0 * acegen_scratch__148__;
  acegen_scratch__363__ = acegen_scratch__192__ * acegen_scratch__502__;
  acegen_scratch__149__ =
    4e0 * (-acegen_scratch__521__ +
           acegen_scratch__114__ * acegen_scratch__79__ + acegen_scratch__95__);
  acegen_scratch__600__ = acegen_scratch__149__ * acegen_scratch__162__;
  acegen_scratch__503__ = 0.7071067811865476e0 * acegen_scratch__149__;
  acegen_scratch__364__ = acegen_scratch__19__ * acegen_scratch__503__;
  acegen_scratch__150__ =
    4e0 * (acegen_scratch__115__ * acegen_scratch__79__ + acegen_scratch__96__);
  acegen_scratch__504__ = 3e0 * acegen_scratch__150__;
  acegen_scratch__365__ = acegen_scratch__192__ * acegen_scratch__504__;
  acegen_scratch__151__ =
    4e0 * (acegen_scratch__111__ * acegen_scratch__80__ + acegen_scratch__93__);
  acegen_scratch__152__ =
    4e0 * (acegen_scratch__112__ * acegen_scratch__80__ + acegen_scratch__94__);
  acegen_scratch__543__ = 0.7071067811865476e0 * acegen_scratch__152__;
  acegen_scratch__506__ = 3e0 * acegen_scratch__152__;
  acegen_scratch__304__ = acegen_scratch__152__ * acegen_scratch__505__;
  acegen_scratch__380__ = acegen_scratch__304__ * acegen_scratch__508__;
  acegen_scratch__153__ =
    4e0 * (acegen_scratch__114__ * acegen_scratch__80__ + acegen_scratch__95__);
  acegen_scratch__542__ = 0.7071067811865476e0 * acegen_scratch__153__;
  acegen_scratch__587__ = 1e0 * acegen_scratch__542__;
  acegen_scratch__329__ = acegen_scratch__153__ * acegen_scratch__507__;
  acegen_scratch__209__ = acegen_scratch__189__ * acegen_scratch__542__;
  acegen_scratch__154__ =
    4e0 * (-acegen_scratch__511__ +
           acegen_scratch__115__ * acegen_scratch__80__ + acegen_scratch__96__);
  acegen_scratch__155__ = -4e0 * acegen_scratch__509__ +
                          acegen_scratch__112__ * acegen_scratch__512__;
  acegen_scratch__599__ = acegen_scratch__147__ + acegen_scratch__155__;
  acegen_scratch__549__ = acegen_scratch__14__ * acegen_scratch__155__;
  acegen_scratch__510__ = acegen_scratch__155__ / 2e0;
  acegen_scratch__287__ = acegen_scratch__155__ * acegen_scratch__505__;
  acegen_scratch__156__ = 0.282842712474619e1 * acegen_scratch__511__ +
                          acegen_scratch__114__ * acegen_scratch__512__;
  acegen_scratch__559__ = acegen_scratch__14__ * acegen_scratch__156__;
  acegen_scratch__547__ = acegen_scratch__156__ * acegen_scratch__20__;
  acegen_scratch__520__ = acegen_scratch__156__ * acegen_scratch__19__;
  acegen_scratch__517__ = acegen_scratch__156__ * acegen_scratch__178__;
  acegen_scratch__514__ = acegen_scratch__156__ * acegen_scratch__193__;
  acegen_scratch__350__ = acegen_scratch__18__ * acegen_scratch__514__;
  acegen_scratch__308__ = acegen_scratch__156__ * acegen_scratch__516__;
  acegen_scratch__291__ = acegen_scratch__282__ * acegen_scratch__547__;
  acegen_scratch__213__ = acegen_scratch__189__ * acegen_scratch__520__;
  acegen_scratch__157__ = 1e0 * acegen_scratch__115__ * acegen_scratch__512__ +
                          0.282842712474619e1 * acegen_scratch__521__;
  acegen_scratch__546__ = acegen_scratch__12__ * acegen_scratch__157__;
  acegen_scratch__524__ = acegen_scratch__15__ * acegen_scratch__157__;
  acegen_scratch__523__ = acegen_scratch__157__ * acegen_scratch__17__;
  acegen_scratch__294__ = acegen_scratch__20__ * acegen_scratch__524__;
  acegen_scratch__292__ = acegen_scratch__18__ * acegen_scratch__523__;
  acegen_scratch__216__ = acegen_scratch__14__ * acegen_scratch__524__;
  acegen_scratch__215__ = acegen_scratch__157__ * acegen_scratch__525__;
  acegen_scratch__158__ = -4e0 * acegen_scratch__526__ +
                          acegen_scratch__114__ * acegen_scratch__529__;
  acegen_scratch__598__ = acegen_scratch__142__ + acegen_scratch__158__;
  acegen_scratch__528__ = acegen_scratch__158__ * acegen_scratch__189__;
  acegen_scratch__527__ = acegen_scratch__158__ / 2e0;
  acegen_scratch__220__ = acegen_scratch__15__ * acegen_scratch__528__;
  acegen_scratch__219__ = acegen_scratch__18__ * acegen_scratch__528__;
  acegen_scratch__159__ = acegen_scratch__115__ * acegen_scratch__529__ +
                          acegen_scratch__483__ * acegen_scratch__536__;
  acegen_scratch__594__ = acegen_scratch__159__ * acegen_scratch__18__;
  acegen_scratch__560__ = acegen_scratch__12__ * acegen_scratch__159__;
  acegen_scratch__535__ = acegen_scratch__159__ * acegen_scratch__189__;
  acegen_scratch__534__ = acegen_scratch__159__ * acegen_scratch__16__;
  acegen_scratch__533__ = acegen_scratch__159__ * acegen_scratch__191__;
  acegen_scratch__351__ = acegen_scratch__19__ * acegen_scratch__533__;
  acegen_scratch__306__ = acegen_scratch__15__ * acegen_scratch__534__;
  acegen_scratch__295__ =
    0.1414213562373095e1 * acegen_scratch__159__ * acegen_scratch__388__;
  acegen_scratch__221__ = acegen_scratch__19__ * acegen_scratch__535__;
  acegen_scratch__160__ =
    0.282842712474619e1 * acegen_scratch__115__ * acegen_scratch__138__ -
    acegen_scratch__29__ * acegen_scratch__536__;
  acegen_scratch__597__ = acegen_scratch__141__ + acegen_scratch__160__;
  acegen_scratch__537__ = acegen_scratch__160__ / 2e0;
  acegen_scratch__284__ = acegen_scratch__160__ * acegen_scratch__495__;
  acegen_scratch__627__ = (acegen_scratch__284__ + acegen_scratch__287__ +
                           acegen_scratch__292__ + acegen_scratch__294__) /
                          2e0;
  acegen_scratch__223__ = acegen_scratch__160__ * acegen_scratch__538__;
  acegen_scratch__222__ = acegen_scratch__160__ * acegen_scratch__577__;
  acegen_scratch__618__ = (acegen_scratch__215__ + acegen_scratch__222__) / 2e0;
  acegen_scratch__611__ = acegen_scratch__215__ + acegen_scratch__222__;
  acegen_scratch__168__ = Power(acegen_scratch__12__, 3);
  acegen_scratch__169__ = 0.282842712474619e1 * acegen_scratch__13__;
  acegen_scratch__171__ = Power(acegen_scratch__13__, 3);
  acegen_scratch__172__ = 0.282842712474619e1 * acegen_scratch__12__;
  acegen_scratch__173__ = Power(acegen_scratch__14__, 3);
  acegen_scratch__179__ = 0.1414213562373095e1 * acegen_scratch__516__;
  acegen_scratch__601__ = acegen_scratch__143__ * acegen_scratch__179__;
  acegen_scratch__180__ = acegen_scratch__176__ * acegen_scratch__508__;
  acegen_scratch__603__ = acegen_scratch__149__ * acegen_scratch__544__;
  acegen_scratch__182__ = 0.1414213562373095e1 * acegen_scratch__545__;
  acegen_scratch__183__ = acegen_scratch__177__ * acegen_scratch__508__;
  acegen_scratch__184__ = (acegen_scratch__13__ * acegen_scratch__172__) / 2e0;
  acegen_scratch__612__ = acegen_scratch__142__ * acegen_scratch__161__ +
                          acegen_scratch__154__ * acegen_scratch__184__;
  acegen_scratch__605__ =
    acegen_scratch__147__ * acegen_scratch__162__ + acegen_scratch__612__;
  acegen_scratch__556__ = acegen_scratch__145__ * acegen_scratch__184__;
  acegen_scratch__186__ = acegen_scratch__16__ * acegen_scratch__548__;
  acegen_scratch__567__ = acegen_scratch__186__ / 2e0;
  acegen_scratch__218__ = acegen_scratch__186__ * acegen_scratch__546__;
  acegen_scratch__211__ = acegen_scratch__186__ * acegen_scratch__549__;
  acegen_scratch__188__ = acegen_scratch__16__ * acegen_scratch__189__;
  acegen_scratch__194__ = 0.1414213562373095e1 * acegen_scratch__532__;
  acegen_scratch__607__ = acegen_scratch__143__ * acegen_scratch__194__;
  acegen_scratch__195__ = 0.7071067811865476e0 * acegen_scratch__281__;
  acegen_scratch__624__ = acegen_scratch__149__ * acegen_scratch__195__;
  acegen_scratch__196__ = 0.1414213562373095e1 * acegen_scratch__515__;
  acegen_scratch__609__ = acegen_scratch__154__ * acegen_scratch__196__;
  acegen_scratch__197__ = acegen_scratch__169__ * acegen_scratch__508__;
  acegen_scratch__198__ = acegen_scratch__19__ * acegen_scratch__548__;
  acegen_scratch__564__ = acegen_scratch__198__ / 2e0;
  acegen_scratch__217__ = acegen_scratch__198__ * acegen_scratch__546__;
  acegen_scratch__210__ = acegen_scratch__198__ * acegen_scratch__549__;
  acegen_scratch__199__ = acegen_scratch__13__ * acegen_scratch__281__;
  acegen_scratch__551__ = acegen_scratch__159__ * acegen_scratch__199__;
  acegen_scratch__550__ = acegen_scratch__156__ * acegen_scratch__199__;
  acegen_scratch__558__ = acegen_scratch__16__ * acegen_scratch__552__ +
                          acegen_scratch__19__ * acegen_scratch__554__;
  acegen_scratch__203__ = acegen_scratch__15__ * acegen_scratch__553__;
  acegen_scratch__204__ = acegen_scratch__18__ * acegen_scratch__554__;
  acegen_scratch__557__ = acegen_scratch__18__ * acegen_scratch__555__ +
                          acegen_scratch__19__ * acegen_scratch__565__;
  acegen_scratch__613__ =
    acegen_scratch__204__ + 1e0 * acegen_scratch__15__ * acegen_scratch__552__;
  acegen_scratch__625__ = acegen_scratch__147__ * acegen_scratch__505__ +
                          acegen_scratch__149__ * acegen_scratch__613__;
  acegen_scratch__212__ = acegen_scratch__17__ * acegen_scratch__561__;
  acegen_scratch__214__ = acegen_scratch__15__ * acegen_scratch__563__;
  acegen_scratch__228__ = 0.7071067811865476e0 * acegen_scratch__561__;
  acegen_scratch__236__ = 0.7071067811865476e0 * acegen_scratch__522__;
  acegen_scratch__578__ =
    acegen_scratch__236__ + 0.7071067811865476e0 * acegen_scratch__525__;
  acegen_scratch__239__ = 0.7071067811865476e0 * acegen_scratch__563__;
  acegen_scratch__245__ = acegen_scratch__13__ * acegen_scratch__565__;
  acegen_scratch__260__ = 1e0 * acegen_scratch__13__ * acegen_scratch__554__;
  acegen_scratch__269__ = 0.282842712474619e1 * acegen_scratch__17__;
  acegen_scratch__270__ = Power(acegen_scratch__15__, 3);
  acegen_scratch__271__ = 0.282842712474619e1 * acegen_scratch__16__;
  acegen_scratch__273__ = Power(acegen_scratch__16__, 3);
  acegen_scratch__274__ = 0.282842712474619e1 * acegen_scratch__15__;
  acegen_scratch__623__ = acegen_scratch__152__ * acegen_scratch__271__ +
                          acegen_scratch__153__ * acegen_scratch__274__;
  acegen_scratch__579__ = acegen_scratch__16__ * acegen_scratch__274__;
  acegen_scratch__604__ = acegen_scratch__142__ * acegen_scratch__176__ +
                          acegen_scratch__147__ * acegen_scratch__177__ +
                          (acegen_scratch__154__ * acegen_scratch__579__) / 2e0;
  acegen_scratch__275__ = Power(acegen_scratch__17__, 3);
  acegen_scratch__277__ = acegen_scratch__191__ * acegen_scratch__570__;
  acegen_scratch__278__ = acegen_scratch__192__ * acegen_scratch__570__;
  acegen_scratch__387__ =
    acegen_scratch__182__ * (acegen_scratch__145__ * acegen_scratch__191__ +
                             acegen_scratch__150__ * acegen_scratch__192__) +
    acegen_scratch__178__ * (acegen_scratch__151__ * acegen_scratch__193__ +
                             acegen_scratch__152__ * acegen_scratch__194__ +
                             acegen_scratch__153__ * acegen_scratch__195__) +
    acegen_scratch__144__ * acegen_scratch__274__ * acegen_scratch__277__ +
    acegen_scratch__148__ * acegen_scratch__271__ * acegen_scratch__278__ +
    acegen_scratch__176__ *
      (acegen_scratch__144__ * acegen_scratch__195__ +
       acegen_scratch__145__ * acegen_scratch__196__ + acegen_scratch__590__) +
    acegen_scratch__177__ *
      (acegen_scratch__148__ * acegen_scratch__194__ +
       acegen_scratch__150__ * acegen_scratch__196__ + acegen_scratch__591__) +
    acegen_scratch__193__ * acegen_scratch__570__ * acegen_scratch__623__;
  acegen_scratch__280__ = 2e0 * acegen_scratch__571__;
  acegen_scratch__588__ = acegen_scratch__280__ / 2e0;
  acegen_scratch__301__ = acegen_scratch__19__ * acegen_scratch__576__;
  acegen_scratch__616__ = acegen_scratch__239__ + acegen_scratch__301__;
  acegen_scratch__303__ = 1e0 * acegen_scratch__19__ * acegen_scratch__498__;
  acegen_scratch__580__ = acegen_scratch__228__ + acegen_scratch__303__;
  acegen_scratch__318__ = acegen_scratch__14__ * acegen_scratch__582__;
  acegen_scratch__619__ = acegen_scratch__260__ + acegen_scratch__318__;
  acegen_scratch__319__ = acegen_scratch__12__ * acegen_scratch__582__;
  acegen_scratch__593__ = (acegen_scratch__245__ + acegen_scratch__319__) / 3e0;
  acegen_scratch__338__ = 0.282842712474619e1 * acegen_scratch__20__;
  acegen_scratch__339__ = Power(acegen_scratch__18__, 3);
  acegen_scratch__340__ = 0.282842712474619e1 * acegen_scratch__19__;
  acegen_scratch__342__ = Power(acegen_scratch__19__, 3);
  acegen_scratch__586__ = 0.7071067811865476e0 * acegen_scratch__342__;
  acegen_scratch__343__ = 0.282842712474619e1 * acegen_scratch__18__;
  acegen_scratch__344__ = Power(acegen_scratch__20__, 3);
  acegen_scratch__355__ = acegen_scratch__191__ * acegen_scratch__589__;
  acegen_scratch__617__ =
    acegen_scratch__615__ *
    (acegen_scratch__154__ * acegen_scratch__281__ * acegen_scratch__301__ +
     (acegen_scratch__142__ + acegen_scratch__158__) * acegen_scratch__355__ +
     acegen_scratch__236__ * acegen_scratch__358__ +
     acegen_scratch__303__ * acegen_scratch__360__ +
     acegen_scratch__239__ * acegen_scratch__363__ +
     acegen_scratch__199__ *
       (acegen_scratch__157__ * acegen_scratch__19__ + acegen_scratch__364__) +
     acegen_scratch__228__ * acegen_scratch__365__ +
     acegen_scratch__13__ *
       (acegen_scratch__146__ * acegen_scratch__342__ + acegen_scratch__350__ +
        acegen_scratch__20__ * acegen_scratch__533__ +
        1e0 * acegen_scratch__344__ * acegen_scratch__543__ +
        acegen_scratch__339__ * acegen_scratch__584__) +
     acegen_scratch__14__ *
       (acegen_scratch__151__ * acegen_scratch__344__ + acegen_scratch__351__ +
        acegen_scratch__281__ * acegen_scratch__520__ +
        acegen_scratch__339__ * acegen_scratch__585__ +
        acegen_scratch__148__ * acegen_scratch__586__) +
     acegen_scratch__12__ *
       (acegen_scratch__140__ * acegen_scratch__339__ +
        acegen_scratch__192__ * acegen_scratch__496__ +
        acegen_scratch__19__ * (acegen_scratch__159__ * acegen_scratch__281__ +
                                acegen_scratch__514__) +
        acegen_scratch__150__ * acegen_scratch__586__ +
        acegen_scratch__344__ * acegen_scratch__587__) +
     acegen_scratch__191__ * acegen_scratch__564__ * acegen_scratch__597__ +
     acegen_scratch__193__ *
       (3e0 * acegen_scratch__18__ * acegen_scratch__329__ +
        acegen_scratch__301__ * acegen_scratch__506__ +
        acegen_scratch__154__ * acegen_scratch__580__ +
        acegen_scratch__577__ * acegen_scratch__598__ +
        acegen_scratch__564__ * acegen_scratch__599__) +
     acegen_scratch__192__ *
       (acegen_scratch__157__ * acegen_scratch__522__ +
        acegen_scratch__149__ * acegen_scratch__578__ +
        acegen_scratch__589__ * acegen_scratch__599__ + acegen_scratch__611__) +
     acegen_scratch__143__ * (acegen_scratch__281__ * acegen_scratch__303__ +
                              acegen_scratch__191__ * acegen_scratch__616__));
  acegen_scratch__378__ = acegen_scratch__538__ / 2e0;
  acegen_scratch__381__ =
    acegen_scratch__14__ * acegen_scratch__15__ * acegen_scratch__19__;
  acegen_scratch__382__ = acegen_scratch__12__ * acegen_scratch__620__;
  acegen_scratch__621__ =
    1e0 * acegen_scratch__615__ *
    (acegen_scratch__281__ * (acegen_scratch__12__ * acegen_scratch__203__ +
                              acegen_scratch__17__ * acegen_scratch__329__) +
     ((acegen_scratch__260__ + acegen_scratch__318__) * acegen_scratch__363__) /
       3e0 +
     (acegen_scratch__358__ * (acegen_scratch__321__ + acegen_scratch__377__)) /
       3e0 +
     (acegen_scratch__143__ * acegen_scratch__20__ * acegen_scratch__340__ +
      acegen_scratch__145__ * acegen_scratch__19__ * acegen_scratch__343__) *
       acegen_scratch__378__ +
     acegen_scratch__340__ * acegen_scratch__380__ +
     acegen_scratch__16__ * acegen_scratch__199__ * acegen_scratch__503__ +
     acegen_scratch__12__ * acegen_scratch__497__ *
       (acegen_scratch__158__ * acegen_scratch__17__ + acegen_scratch__534__) +
     acegen_scratch__192__ * (acegen_scratch__149__ * (acegen_scratch__321__ +
                                                       acegen_scratch__377__) +
                              acegen_scratch__17__ * acegen_scratch__501__ +
                              acegen_scratch__141__ * acegen_scratch__538__) +
     (acegen_scratch__381__ + acegen_scratch__382__) * acegen_scratch__547__ +
     (acegen_scratch__17__ * acegen_scratch__550__) / 2e0 +
     (acegen_scratch__15__ * acegen_scratch__551__) / 2e0 +
     acegen_scratch__20__ * (acegen_scratch__546__ + acegen_scratch__549__) *
       acegen_scratch__588__ +
     acegen_scratch__538__ * acegen_scratch__590__ +
     acegen_scratch__567__ * acegen_scratch__591__ +
     (acegen_scratch__360__ + acegen_scratch__365__) * acegen_scratch__593__ +
     (acegen_scratch__381__ + acegen_scratch__382__) * acegen_scratch__594__ +
     acegen_scratch__280__ *
       (acegen_scratch__150__ * acegen_scratch__228__ +
        acegen_scratch__148__ * acegen_scratch__239__ + acegen_scratch__618__) +
     acegen_scratch__191__ * (acegen_scratch__141__ * acegen_scratch__567__ +
                              acegen_scratch__142__ * acegen_scratch__592__ +
                              acegen_scratch__143__ * acegen_scratch__619__) +
     acegen_scratch__193__ * (acegen_scratch__154__ * (acegen_scratch__245__ +
                                                       acegen_scratch__319__) +
                              acegen_scratch__153__ * acegen_scratch__321__ +
                              acegen_scratch__15__ * acegen_scratch__329__ +
                              acegen_scratch__142__ * acegen_scratch__538__ +
                              acegen_scratch__147__ * acegen_scratch__567__ +
                              acegen_scratch__151__ * acegen_scratch__592__ +
                              acegen_scratch__152__ * acegen_scratch__619__) +
     acegen_scratch__508__ *
       ((acegen_scratch__15__ * acegen_scratch__158__ +
         acegen_scratch__156__ * acegen_scratch__16__) *
          acegen_scratch__281__ +
        acegen_scratch__154__ * acegen_scratch__343__ * acegen_scratch__620__) +
     acegen_scratch__198__ * acegen_scratch__627__);
  acegen_scratch__628__ =
    2e0 * acegen_scratch__175__ *
    (acegen_scratch__387__ +
     acegen_scratch__19__ *
       (acegen_scratch__291__ / 2e0 + acegen_scratch__295__ / 2e0 +
        acegen_scratch__143__ * acegen_scratch__388__) +
     acegen_scratch__178__ * (acegen_scratch__156__ * acegen_scratch__515__ +
                              acegen_scratch__191__ * acegen_scratch__527__) +
     acegen_scratch__176__ * (acegen_scratch__193__ * acegen_scratch__527__ +
                              acegen_scratch__159__ * acegen_scratch__532__) +
     acegen_scratch__16__ * (acegen_scratch__15__ * acegen_scratch__514__ +
                             acegen_scratch__17__ * acegen_scratch__533__) +
     acegen_scratch__192__ * (acegen_scratch__178__ * acegen_scratch__510__ +
                              acegen_scratch__176__ * acegen_scratch__537__) +
     acegen_scratch__177__ * (acegen_scratch__157__ * acegen_scratch__497__ +
                              acegen_scratch__193__ * acegen_scratch__510__ +
                              acegen_scratch__191__ * acegen_scratch__537__) +
     acegen_scratch__281__ *
       (acegen_scratch__306__ / 2e0 + acegen_scratch__308__ / 2e0 +
        1e0 *
          (acegen_scratch__143__ * acegen_scratch__15__ +
           acegen_scratch__154__ * acegen_scratch__17__) *
          acegen_scratch__582__) +
     acegen_scratch__282__ *
       ((acegen_scratch__157__ * acegen_scratch__192__) / 2e0 +
        acegen_scratch__142__ * acegen_scratch__497__ +
        acegen_scratch__154__ * acegen_scratch__19__ * acegen_scratch__552__ +
        acegen_scratch__158__ * acegen_scratch__595__) +
     acegen_scratch__280__ * (acegen_scratch__141__ * acegen_scratch__495__ +
                              acegen_scratch__625__ + acegen_scratch__627__));
  acegen_scratch__614__ =
    acegen_scratch__615__ *
    ((acegen_scratch__20__ *
      (acegen_scratch__211__ + acegen_scratch__218__ + acegen_scratch__220__)) /
       2e0 +
     acegen_scratch__197__ * acegen_scratch__304__ +
     acegen_scratch__15__ * (1e0 * acegen_scratch__143__ *
                               acegen_scratch__18__ * acegen_scratch__197__ +
                             acegen_scratch__221__ / 2e0 +
                             acegen_scratch__162__ * acegen_scratch__496__) +
     acegen_scratch__209__ * acegen_scratch__505__ +
     acegen_scratch__198__ *
       (acegen_scratch__216__ / 2e0 + acegen_scratch__223__ / 2e0 +
        (acegen_scratch__14__ * acegen_scratch__148__ +
         acegen_scratch__12__ * acegen_scratch__150__) *
          acegen_scratch__555__) +
     acegen_scratch__18__ * (acegen_scratch__189__ * acegen_scratch__203__ +
                             acegen_scratch__15__ * acegen_scratch__556__) +
     (acegen_scratch__212__ + acegen_scratch__214__) * acegen_scratch__559__ +
     (acegen_scratch__212__ + acegen_scratch__214__) * acegen_scratch__560__ +
     acegen_scratch__161__ * (acegen_scratch__20__ * acegen_scratch__203__ +
                              acegen_scratch__144__ * acegen_scratch__204__ +
                              acegen_scratch__140__ * acegen_scratch__495__ +
                              acegen_scratch__145__ * acegen_scratch__557__ +
                              acegen_scratch__143__ * acegen_scratch__558__ +
                              acegen_scratch__141__ * acegen_scratch__571__) +
     acegen_scratch__188__ *
       (acegen_scratch__364__ + acegen_scratch__547__ / 2e0 +
        acegen_scratch__594__ / 2e0) +
     acegen_scratch__17__ *
       (acegen_scratch__210__ / 2e0 + acegen_scratch__213__ / 2e0 +
        acegen_scratch__217__ / 2e0 + acegen_scratch__219__ / 2e0 +
        acegen_scratch__20__ * acegen_scratch__612__) +
     acegen_scratch__164__ * (acegen_scratch__142__ * acegen_scratch__495__ +
                              acegen_scratch__151__ * acegen_scratch__505__ +
                              acegen_scratch__154__ * acegen_scratch__557__ +
                              acegen_scratch__152__ * acegen_scratch__558__ +
                              acegen_scratch__147__ * acegen_scratch__571__ +
                              acegen_scratch__153__ * acegen_scratch__613__) +
     acegen_scratch__186__ * acegen_scratch__618__ +
     acegen_scratch__162__ *
       (acegen_scratch__150__ * acegen_scratch__557__ +
        acegen_scratch__148__ * acegen_scratch__558__ +
        acegen_scratch__146__ * acegen_scratch__571__ + acegen_scratch__625__));
  gradientOut[0][0] = acegen_scratch__66__;
  gradientOut[0][1] = acegen_scratch__67__;
  gradientOut[0][2] = acegen_scratch__68__;
  gradientOut[1][0] = acegen_scratch__69__;
  gradientOut[1][1] = acegen_scratch__70__;
  gradientOut[1][2] = acegen_scratch__71__;
  gradientOut[2][0] = acegen_scratch__72__;
  gradientOut[2][1] = acegen_scratch__73__;
  gradientOut[2][2] = acegen_scratch__74__;
  cacheOut[0] =
    acegen_scratch__175__ *
    (Power(acegen_scratch__12__, 4) * acegen_scratch__140__ +
     Power(acegen_scratch__13__, 4) * acegen_scratch__146__ +
     Power(acegen_scratch__14__, 4) * acegen_scratch__151__ +
     acegen_scratch__172__ * (acegen_scratch__150__ * acegen_scratch__171__ +
                              acegen_scratch__153__ * acegen_scratch__173__) +
     4e0 * acegen_scratch__13__ * acegen_scratch__156__ *
       acegen_scratch__519__ +
     acegen_scratch__169__ * (acegen_scratch__145__ * acegen_scratch__168__ +
                              acegen_scratch__152__ * acegen_scratch__173__ +
                              acegen_scratch__154__ * acegen_scratch__519__) +
     (0.282842712474619e1 * acegen_scratch__143__ +
      4e0 * acegen_scratch__159__) *
       acegen_scratch__161__ * acegen_scratch__530__ +
     acegen_scratch__164__ * acegen_scratch__539__ * acegen_scratch__598__ +
     acegen_scratch__162__ *
       (4e0 * acegen_scratch__14__ * acegen_scratch__546__ +
        acegen_scratch__539__ * acegen_scratch__597__ +
        2e0 * acegen_scratch__164__ * acegen_scratch__599__) +
     0.282842712474619e1 * acegen_scratch__14__ *
       (acegen_scratch__144__ * acegen_scratch__168__ +
        acegen_scratch__148__ * acegen_scratch__171__ +
        acegen_scratch__12__ * acegen_scratch__600__));
  cacheOut[1] =
    acegen_scratch__175__ *
    ((acegen_scratch__143__ * acegen_scratch__169__ +
      acegen_scratch__144__ * acegen_scratch__172__) *
       acegen_scratch__180__ +
     (acegen_scratch__148__ * acegen_scratch__169__ +
      acegen_scratch__149__ * acegen_scratch__172__) *
       acegen_scratch__183__ +
     (acegen_scratch__145__ * acegen_scratch__176__ +
      acegen_scratch__150__ * acegen_scratch__177__) *
       acegen_scratch__184__ +
     acegen_scratch__15__ * acegen_scratch__159__ * acegen_scratch__188__ +
     acegen_scratch__17__ *
       (acegen_scratch__156__ * acegen_scratch__188__ + acegen_scratch__211__ +
        acegen_scratch__218__ + acegen_scratch__220__) +
     acegen_scratch__186__ * (acegen_scratch__216__ + acegen_scratch__223__) +
     acegen_scratch__13__ * acegen_scratch__282__ *
       (acegen_scratch__559__ + acegen_scratch__560__) +
     acegen_scratch__161__ *
       (acegen_scratch__140__ * acegen_scratch__176__ +
        acegen_scratch__141__ * acegen_scratch__177__ +
        acegen_scratch__145__ * acegen_scratch__182__ +
        acegen_scratch__144__ * acegen_scratch__544__ + acegen_scratch__601__) +
     acegen_scratch__162__ *
       (acegen_scratch__141__ * acegen_scratch__176__ +
        acegen_scratch__146__ * acegen_scratch__177__ +
        acegen_scratch__148__ * acegen_scratch__179__ +
        acegen_scratch__150__ * acegen_scratch__182__ + acegen_scratch__603__) +
     acegen_scratch__164__ *
       (acegen_scratch__151__ * acegen_scratch__178__ +
        acegen_scratch__152__ * acegen_scratch__179__ +
        acegen_scratch__153__ * acegen_scratch__544__ + acegen_scratch__604__) +
     acegen_scratch__178__ * ((acegen_scratch__152__ * acegen_scratch__169__ +
                               acegen_scratch__153__ * acegen_scratch__172__) *
                                acegen_scratch__508__ +
                              acegen_scratch__605__));
  cacheOut[2] =
    acegen_scratch__175__ *
    ((acegen_scratch__148__ * acegen_scratch__162__ +
      acegen_scratch__152__ * acegen_scratch__164__) *
       acegen_scratch__194__ +
     (acegen_scratch__145__ * acegen_scratch__161__ +
      acegen_scratch__150__ * acegen_scratch__162__) *
       acegen_scratch__196__ +
     acegen_scratch__20__ * (acegen_scratch__210__ + acegen_scratch__213__ +
                             acegen_scratch__217__ + acegen_scratch__219__) +
     acegen_scratch__18__ * acegen_scratch__221__ +
     acegen_scratch__192__ *
       (acegen_scratch__141__ * acegen_scratch__161__ +
        acegen_scratch__146__ * acegen_scratch__162__ +
        acegen_scratch__147__ * acegen_scratch__164__ +
        acegen_scratch__150__ * acegen_scratch__184__ +
        acegen_scratch__148__ * acegen_scratch__197__ +
        1e0 * acegen_scratch__189__ * acegen_scratch__503__) +
     acegen_scratch__14__ * acegen_scratch__550__ +
     acegen_scratch__12__ * acegen_scratch__551__ +
     acegen_scratch__191__ *
       (acegen_scratch__140__ * acegen_scratch__161__ +
        acegen_scratch__141__ * acegen_scratch__162__ +
        acegen_scratch__142__ * acegen_scratch__164__ +
        acegen_scratch__143__ * acegen_scratch__197__ +
        acegen_scratch__189__ * acegen_scratch__553__ + acegen_scratch__556__) +
     acegen_scratch__195__ *
       (acegen_scratch__144__ * acegen_scratch__161__ +
        acegen_scratch__153__ * acegen_scratch__164__ + acegen_scratch__600__) +
     acegen_scratch__193__ * (acegen_scratch__151__ * acegen_scratch__164__ +
                              acegen_scratch__152__ * acegen_scratch__197__ +
                              acegen_scratch__209__ + acegen_scratch__605__) +
     acegen_scratch__161__ * acegen_scratch__607__ +
     acegen_scratch__164__ * acegen_scratch__609__ +
     acegen_scratch__198__ * acegen_scratch__611__);
  cacheOut[3] = acegen_scratch__614__;
  cacheOut[4] = acegen_scratch__617__;
  cacheOut[5] = acegen_scratch__621__;
  cacheOut[6] =
    acegen_scratch__175__ *
    (acegen_scratch__140__ * Power(acegen_scratch__15__, 4) +
     acegen_scratch__146__ * Power(acegen_scratch__16__, 4) +
     acegen_scratch__151__ * Power(acegen_scratch__17__, 4) +
     2e0 * acegen_scratch__154__ * acegen_scratch__178__ *
       acegen_scratch__182__ +
     acegen_scratch__270__ * (acegen_scratch__144__ * acegen_scratch__269__ +
                              acegen_scratch__145__ * acegen_scratch__271__) +
     acegen_scratch__273__ * (acegen_scratch__148__ * acegen_scratch__269__ +
                              acegen_scratch__150__ * acegen_scratch__274__) +
     4e0 * acegen_scratch__15__ *
       (acegen_scratch__16__ * acegen_scratch__517__ +
        acegen_scratch__177__ * acegen_scratch__523__) +
     4e0 * acegen_scratch__159__ * acegen_scratch__16__ *
       acegen_scratch__581__ +
     acegen_scratch__568__ *
       (acegen_scratch__177__ * acegen_scratch__597__ +
        acegen_scratch__178__ * acegen_scratch__598__ + acegen_scratch__601__) +
     acegen_scratch__602__ *
       (acegen_scratch__178__ * acegen_scratch__599__ + acegen_scratch__603__) +
     acegen_scratch__275__ * acegen_scratch__623__);
  cacheOut[7] =
    acegen_scratch__175__ *
    (acegen_scratch__143__ * acegen_scratch__271__ * acegen_scratch__277__ +
     acegen_scratch__149__ * acegen_scratch__274__ * acegen_scratch__278__ +
     acegen_scratch__280__ * (acegen_scratch__284__ + acegen_scratch__287__ +
                              acegen_scratch__292__ + acegen_scratch__294__) +
     acegen_scratch__19__ * (acegen_scratch__291__ + acegen_scratch__295__) +
     acegen_scratch__387__ +
     acegen_scratch__281__ * (acegen_scratch__306__ + acegen_scratch__308__ +
                              acegen_scratch__282__ * acegen_scratch__527__) +
     acegen_scratch__193__ * acegen_scratch__604__ +
     acegen_scratch__176__ *
       (acegen_scratch__141__ * acegen_scratch__192__ + acegen_scratch__607__) +
     acegen_scratch__178__ *
       (acegen_scratch__142__ * acegen_scratch__191__ +
        acegen_scratch__147__ * acegen_scratch__192__ + acegen_scratch__609__) +
     acegen_scratch__177__ *
       (acegen_scratch__141__ * acegen_scratch__191__ + acegen_scratch__624__));
  cacheOut[8]  = acegen_scratch__614__;
  cacheOut[9]  = acegen_scratch__617__;
  cacheOut[10] = acegen_scratch__621__;
  cacheOut[11] =
    acegen_scratch__175__ *
    (acegen_scratch__140__ * Power(acegen_scratch__18__, 4) +
     acegen_scratch__146__ * Power(acegen_scratch__19__, 4) +
     acegen_scratch__151__ * Power(acegen_scratch__20__, 4) +
     acegen_scratch__339__ * (acegen_scratch__144__ * acegen_scratch__338__ +
                              acegen_scratch__145__ * acegen_scratch__340__) +
     acegen_scratch__342__ * (acegen_scratch__148__ * acegen_scratch__338__ +
                              acegen_scratch__150__ * acegen_scratch__343__) +
     (acegen_scratch__152__ * acegen_scratch__340__ +
      acegen_scratch__153__ * acegen_scratch__343__) *
       acegen_scratch__344__ +
     4e0 * acegen_scratch__19__ * acegen_scratch__350__ +
     4e0 * acegen_scratch__20__ * acegen_scratch__351__ +
     acegen_scratch__583__ *
       (acegen_scratch__192__ * acegen_scratch__597__ +
        acegen_scratch__193__ * acegen_scratch__598__ + acegen_scratch__607__) +
     acegen_scratch__608__ *
       (acegen_scratch__192__ * acegen_scratch__599__ + acegen_scratch__609__) +
     2e0 * acegen_scratch__192__ *
       (acegen_scratch__157__ * acegen_scratch__281__ + acegen_scratch__624__));
  cacheOut[12] = acegen_scratch__614__;
  cacheOut[13] = acegen_scratch__617__;
  cacheOut[14] = acegen_scratch__621__;
  cacheOut[15] = acegen_scratch__628__;
  cacheOut[16] = 0e0;
  cacheOut[17] = 0e0;
  cacheOut[18] = acegen_scratch__628__;
  cacheOut[19] = 0e0;
  cacheOut[20] = acegen_scratch__628__;
  cacheOut[21] =
    acegen_scratch__175__ * (acegen_scratch__12__ * acegen_scratch__66__ +
                             acegen_scratch__13__ * acegen_scratch__67__ +
                             acegen_scratch__14__ * acegen_scratch__68__);
  cacheOut[22] =
    acegen_scratch__175__ * (acegen_scratch__15__ * acegen_scratch__66__ +
                             acegen_scratch__16__ * acegen_scratch__67__ +
                             acegen_scratch__17__ * acegen_scratch__68__);
  cacheOut[23] =
    acegen_scratch__175__ * (acegen_scratch__18__ * acegen_scratch__66__ +
                             acegen_scratch__19__ * acegen_scratch__67__ +
                             acegen_scratch__20__ * acegen_scratch__68__);
  cacheOut[24] =
    acegen_scratch__175__ * (acegen_scratch__15__ * acegen_scratch__69__ +
                             acegen_scratch__16__ * acegen_scratch__70__ +
                             acegen_scratch__17__ * acegen_scratch__71__);
  cacheOut[25] =
    acegen_scratch__175__ * (acegen_scratch__18__ * acegen_scratch__69__ +
                             acegen_scratch__19__ * acegen_scratch__70__ +
                             acegen_scratch__20__ * acegen_scratch__71__);
  cacheOut[26] =
    acegen_scratch__175__ * (acegen_scratch__18__ * acegen_scratch__72__ +
                             acegen_scratch__19__ * acegen_scratch__73__ +
                             acegen_scratch__20__ * acegen_scratch__74__);
}