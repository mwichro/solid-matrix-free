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
struct NeoHookean
{
  const static unsigned int dim = DIMENSION;
  template <typename Number>
  static inline void
  Residual_local([[maybe_unused]] const Tensor<1, dim, Number> &uIn,
                 const Tensor<2, dim, Number> &                 graduIn,
                 Tensor<1, dim, Number> &                       valueOut,
                 Tensor<2, dim, Number> &                       gradientOut);

  template <typename Number>
  static inline void
  Tangent_local(const Tensor<2, dim, Number> &graduIn,
                const Tensor<2, dim, Number> &gradduIn,
                Tensor<2, dim, Number> &      gradientOut);
};



// =====
// dim=3
// =====


template <>
template <typename Number>
inline void
NeoHookean<2>::Residual_local(
  [[maybe_unused]] const Tensor<1, dim, Number> &uIn,
  const Tensor<2, dim, Number> &                 graduIn,
  Tensor<1, dim, Number> &                       valueOut,
  Tensor<2, dim, Number> &                       gradientOut)
{}


template <>
template <typename Number>
inline void
NeoHookean<2>::Tangent_local(const Tensor<2, dim, Number> &graduIn,
                             const Tensor<2, dim, Number> &gradduIn,
                             Tensor<2, dim, Number> &      gradientOut)
{
  Number acegen_scratch__10__, acegen_scratch__100__, acegen_scratch__101__,
    acegen_scratch__104__, acegen_scratch__11__, acegen_scratch__110__,
    acegen_scratch__119__, acegen_scratch__12__, acegen_scratch__120__,
    acegen_scratch__121__, acegen_scratch__122__, acegen_scratch__123__,
    acegen_scratch__13__, acegen_scratch__131__, acegen_scratch__132__,
    acegen_scratch__133__, acegen_scratch__134__, acegen_scratch__135__,
    acegen_scratch__136__, acegen_scratch__14__, acegen_scratch__15__,
    acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__19__, acegen_scratch__23__, acegen_scratch__24__,
    acegen_scratch__25__, acegen_scratch__26__, acegen_scratch__41__,
    acegen_scratch__42__, acegen_scratch__43__, acegen_scratch__44__,
    acegen_scratch__50__, acegen_scratch__56__, acegen_scratch__62__,
    acegen_scratch__67__, acegen_scratch__9__, acegen_scratch__94__,
    acegen_scratch__95__, acegen_scratch__97__, acegen_scratch__99__;
  acegen_scratch__9__   = gradduIn[0][0];
  acegen_scratch__94__  = 2e0 * acegen_scratch__9__;
  acegen_scratch__10__  = gradduIn[0][1];
  acegen_scratch__95__  = 2e0 * acegen_scratch__10__;
  acegen_scratch__11__  = gradduIn[1][0];
  acegen_scratch__97__  = 2e0 * acegen_scratch__11__;
  acegen_scratch__12__  = gradduIn[1][1];
  acegen_scratch__99__  = 2e0 * acegen_scratch__12__;
  acegen_scratch__13__  = 1e0 + graduIn[0][0];
  acegen_scratch__23__  = 2e0 * acegen_scratch__13__;
  acegen_scratch__14__  = graduIn[0][1];
  acegen_scratch__25__  = 2e0 * acegen_scratch__14__;
  acegen_scratch__15__  = graduIn[1][0];
  acegen_scratch__100__ = acegen_scratch__13__ * acegen_scratch__94__ +
                          acegen_scratch__15__ * acegen_scratch__97__;
  acegen_scratch__24__  = 2e0 * acegen_scratch__15__;
  acegen_scratch__16__  = 1e0 + graduIn[1][1];
  acegen_scratch__101__ = acegen_scratch__10__ * acegen_scratch__13__ +
                          acegen_scratch__12__ * acegen_scratch__15__ +
                          acegen_scratch__11__ * acegen_scratch__16__ +
                          acegen_scratch__14__ * acegen_scratch__9__;
  acegen_scratch__104__ = acegen_scratch__14__ * acegen_scratch__95__ +
                          acegen_scratch__16__ * acegen_scratch__99__;
  acegen_scratch__26__ = 2e0 * acegen_scratch__16__;
  acegen_scratch__17__ = (acegen_scratch__13__ * acegen_scratch__13__) +
                         (acegen_scratch__15__ * acegen_scratch__15__);
  acegen_scratch__18__ = acegen_scratch__13__ * acegen_scratch__14__ +
                         acegen_scratch__15__ * acegen_scratch__16__;
  acegen_scratch__62__ = -(acegen_scratch__18__ * acegen_scratch__23__) +
                         acegen_scratch__17__ * acegen_scratch__25__;
  acegen_scratch__50__ = -(acegen_scratch__18__ * acegen_scratch__24__) +
                         acegen_scratch__17__ * acegen_scratch__26__;
  acegen_scratch__19__ = (acegen_scratch__14__ * acegen_scratch__14__) +
                         (acegen_scratch__16__ * acegen_scratch__16__);
  acegen_scratch__131__ = -(acegen_scratch__18__ * acegen_scratch__18__) +
                          acegen_scratch__17__ * acegen_scratch__19__;
  acegen_scratch__133__ = sqrt(acegen_scratch__131__);
  acegen_scratch__132__ = 1e0 / (2e0 * acegen_scratch__133__);
  acegen_scratch__119__ = acegen_scratch__132__ *
                          (acegen_scratch__104__ * acegen_scratch__17__ -
                           2e0 * acegen_scratch__101__ * acegen_scratch__18__ +
                           acegen_scratch__100__ * acegen_scratch__19__);
  acegen_scratch__110__ =
    -0.5e0 * acegen_scratch__119__ / acegen_scratch__131__;
  acegen_scratch__67__ = acegen_scratch__19__ * acegen_scratch__23__ -
                         acegen_scratch__18__ * acegen_scratch__25__;
  acegen_scratch__56__ = acegen_scratch__19__ * acegen_scratch__24__ -
                         acegen_scratch__18__ * acegen_scratch__26__;
  acegen_scratch__44__ = acegen_scratch__132__ * acegen_scratch__50__;
  acegen_scratch__43__ = acegen_scratch__132__ * acegen_scratch__56__;
  acegen_scratch__42__ = acegen_scratch__132__ * acegen_scratch__62__;
  acegen_scratch__41__ = acegen_scratch__132__ * acegen_scratch__67__;
  acegen_scratch__134__ =
    -(acegen_scratch__119__ / (acegen_scratch__133__ * acegen_scratch__133__));
  acegen_scratch__136__ = ((350000e0 / 81e0) * acegen_scratch__119__) /
                          (acegen_scratch__133__ * acegen_scratch__133__);
  acegen_scratch__123__ =
    acegen_scratch__134__ * acegen_scratch__41__ +
    (acegen_scratch__110__ * acegen_scratch__67__ +
     acegen_scratch__132__ * (acegen_scratch__104__ * acegen_scratch__23__ -
                              acegen_scratch__101__ * acegen_scratch__25__ +
                              acegen_scratch__19__ * acegen_scratch__94__ -
                              acegen_scratch__18__ * acegen_scratch__95__)) /
      acegen_scratch__133__;
  acegen_scratch__122__ =
    acegen_scratch__134__ * acegen_scratch__42__ +
    (acegen_scratch__110__ * acegen_scratch__62__ +
     acegen_scratch__132__ * (-(acegen_scratch__101__ * acegen_scratch__23__) +
                              acegen_scratch__100__ * acegen_scratch__25__ -
                              acegen_scratch__18__ * acegen_scratch__94__ +
                              acegen_scratch__17__ * acegen_scratch__95__)) /
      acegen_scratch__133__;
  acegen_scratch__121__ =
    acegen_scratch__134__ * acegen_scratch__43__ +
    (acegen_scratch__110__ * acegen_scratch__56__ +
     acegen_scratch__132__ * (acegen_scratch__104__ * acegen_scratch__24__ -
                              acegen_scratch__101__ * acegen_scratch__26__ +
                              acegen_scratch__19__ * acegen_scratch__97__ -
                              acegen_scratch__18__ * acegen_scratch__99__)) /
      acegen_scratch__133__;
  acegen_scratch__120__ =
    acegen_scratch__134__ * acegen_scratch__44__ +
    (acegen_scratch__110__ * acegen_scratch__50__ +
     acegen_scratch__132__ * (-(acegen_scratch__101__ * acegen_scratch__24__) +
                              acegen_scratch__100__ * acegen_scratch__26__ -
                              acegen_scratch__18__ * acegen_scratch__97__ +
                              acegen_scratch__17__ * acegen_scratch__99__)) /
      acegen_scratch__133__;
  acegen_scratch__135__ = (350000e0 / 81e0) * log(acegen_scratch__133__);
  gradientOut[0][0] =
    acegen_scratch__123__ * acegen_scratch__135__ +
    acegen_scratch__136__ * acegen_scratch__41__ +
    (12500e0 / 27e0) * (-2e0 * acegen_scratch__123__ + acegen_scratch__94__);
  gradientOut[0][1] =
    acegen_scratch__122__ * acegen_scratch__135__ +
    acegen_scratch__136__ * acegen_scratch__42__ +
    (12500e0 / 27e0) * (-2e0 * acegen_scratch__122__ + acegen_scratch__95__);
  gradientOut[1][0] =
    acegen_scratch__121__ * acegen_scratch__135__ +
    acegen_scratch__136__ * acegen_scratch__43__ +
    (12500e0 / 27e0) * (-2e0 * acegen_scratch__121__ + acegen_scratch__97__);
  gradientOut[1][1] =
    acegen_scratch__120__ * acegen_scratch__135__ +
    acegen_scratch__136__ * acegen_scratch__44__ +
    (12500e0 / 27e0) * (-2e0 * acegen_scratch__120__ + acegen_scratch__99__);
}



// =====
// dim=3
// =====


template <>
template <typename Number>
inline void
NeoHookean<3>::Residual_local(
  [[maybe_unused]] const Tensor<1, dim, Number> &uIn,
  const Tensor<2, dim, Number> &                 graduIn,
  Tensor<1, dim, Number> &                       valueOut,
  Tensor<2, dim, Number> &                       gradientOut)
{}



template <>
template <typename Number>
inline void
NeoHookean<3>::Tangent_local(const Tensor<2, dim, Number> &graduIn,
                             const Tensor<2, dim, Number> &gradduIn,
                             Tensor<2, dim, Number> &      gradientOut)
{
  Number acegen_scratch__143__, acegen_scratch__153__, acegen_scratch__16__,
    acegen_scratch__166__, acegen_scratch__17__, acegen_scratch__179__,
    acegen_scratch__18__, acegen_scratch__19__, acegen_scratch__195__,
    acegen_scratch__20__, acegen_scratch__21__, acegen_scratch__214__,
    acegen_scratch__22__, acegen_scratch__23__, acegen_scratch__24__,
    acegen_scratch__25__, acegen_scratch__254__, acegen_scratch__26__,
    acegen_scratch__261__, acegen_scratch__267__, acegen_scratch__27__,
    acegen_scratch__279__, acegen_scratch__28__, acegen_scratch__280__,
    acegen_scratch__282__, acegen_scratch__283__, acegen_scratch__284__,
    acegen_scratch__285__, acegen_scratch__287__, acegen_scratch__29__,
    acegen_scratch__290__, acegen_scratch__291__, acegen_scratch__30__,
    acegen_scratch__300__, acegen_scratch__301__, acegen_scratch__302__,
    acegen_scratch__303__, acegen_scratch__304__, acegen_scratch__305__,
    acegen_scratch__306__, acegen_scratch__31__, acegen_scratch__318__,
    acegen_scratch__319__, acegen_scratch__32__, acegen_scratch__320__,
    acegen_scratch__321__, acegen_scratch__322__, acegen_scratch__33__,
    acegen_scratch__34__, acegen_scratch__35__, acegen_scratch__36__,
    acegen_scratch__37__, acegen_scratch__38__, acegen_scratch__39__,
    acegen_scratch__40__, acegen_scratch__43__, acegen_scratch__44__,
    acegen_scratch__45__, acegen_scratch__47__, acegen_scratch__48__,
    acegen_scratch__49__, acegen_scratch__50__, acegen_scratch__52__,
    acegen_scratch__53__, acegen_scratch__54__, acegen_scratch__65__,
    acegen_scratch__66__, acegen_scratch__67__, acegen_scratch__83__,
    acegen_scratch__84__, acegen_scratch__85__, acegen_scratch__92__,
    acegen_scratch__93__, acegen_scratch__94__;
  acegen_scratch__16__  = gradduIn[0][0];
  acegen_scratch__17__  = gradduIn[0][1];
  acegen_scratch__18__  = gradduIn[0][2];
  acegen_scratch__19__  = gradduIn[1][0];
  acegen_scratch__20__  = gradduIn[1][1];
  acegen_scratch__21__  = gradduIn[1][2];
  acegen_scratch__22__  = gradduIn[2][0];
  acegen_scratch__23__  = gradduIn[2][1];
  acegen_scratch__24__  = gradduIn[2][2];
  acegen_scratch__25__  = 1e0 + graduIn[0][0];
  acegen_scratch__65__  = 2e0 * acegen_scratch__25__;
  acegen_scratch__26__  = graduIn[0][1];
  acegen_scratch__83__  = 2e0 * acegen_scratch__26__;
  acegen_scratch__27__  = graduIn[0][2];
  acegen_scratch__92__  = 2e0 * acegen_scratch__27__;
  acegen_scratch__28__  = graduIn[1][0];
  acegen_scratch__66__  = 2e0 * acegen_scratch__28__;
  acegen_scratch__29__  = 1e0 + graduIn[1][1];
  acegen_scratch__84__  = 2e0 * acegen_scratch__29__;
  acegen_scratch__30__  = graduIn[1][2];
  acegen_scratch__93__  = 2e0 * acegen_scratch__30__;
  acegen_scratch__31__  = graduIn[2][0];
  acegen_scratch__67__  = 2e0 * acegen_scratch__31__;
  acegen_scratch__279__ = acegen_scratch__16__ * acegen_scratch__65__ +
                          acegen_scratch__19__ * acegen_scratch__66__ +
                          acegen_scratch__22__ * acegen_scratch__67__;
  acegen_scratch__32__  = graduIn[2][1];
  acegen_scratch__280__ = acegen_scratch__17__ * acegen_scratch__25__ +
                          acegen_scratch__16__ * acegen_scratch__26__ +
                          acegen_scratch__20__ * acegen_scratch__28__ +
                          acegen_scratch__19__ * acegen_scratch__29__ +
                          acegen_scratch__23__ * acegen_scratch__31__ +
                          acegen_scratch__22__ * acegen_scratch__32__;
  acegen_scratch__85__  = 2e0 * acegen_scratch__32__;
  acegen_scratch__285__ = acegen_scratch__17__ * acegen_scratch__83__ +
                          acegen_scratch__20__ * acegen_scratch__84__ +
                          acegen_scratch__23__ * acegen_scratch__85__;
  acegen_scratch__33__  = 1e0 + graduIn[2][2];
  acegen_scratch__287__ = acegen_scratch__18__ * acegen_scratch__26__ +
                          acegen_scratch__17__ * acegen_scratch__27__ +
                          acegen_scratch__21__ * acegen_scratch__29__ +
                          acegen_scratch__20__ * acegen_scratch__30__ +
                          acegen_scratch__24__ * acegen_scratch__32__ +
                          acegen_scratch__23__ * acegen_scratch__33__;
  acegen_scratch__282__ = acegen_scratch__18__ * acegen_scratch__25__ +
                          acegen_scratch__16__ * acegen_scratch__27__ +
                          acegen_scratch__21__ * acegen_scratch__28__ +
                          acegen_scratch__19__ * acegen_scratch__30__ +
                          acegen_scratch__24__ * acegen_scratch__31__ +
                          acegen_scratch__22__ * acegen_scratch__33__;
  acegen_scratch__94__  = 2e0 * acegen_scratch__33__;
  acegen_scratch__291__ = acegen_scratch__18__ * acegen_scratch__92__ +
                          acegen_scratch__21__ * acegen_scratch__93__ +
                          acegen_scratch__24__ * acegen_scratch__94__;
  acegen_scratch__34__ = (acegen_scratch__25__ * acegen_scratch__25__) +
                         (acegen_scratch__28__ * acegen_scratch__28__) +
                         (acegen_scratch__31__ * acegen_scratch__31__);
  acegen_scratch__35__ = acegen_scratch__25__ * acegen_scratch__26__ +
                         acegen_scratch__28__ * acegen_scratch__29__ +
                         acegen_scratch__31__ * acegen_scratch__32__;
  acegen_scratch__320__ = 2e0 * acegen_scratch__35__;
  acegen_scratch__44__  = (acegen_scratch__35__ * acegen_scratch__35__);
  acegen_scratch__36__  = acegen_scratch__25__ * acegen_scratch__27__ +
                         acegen_scratch__28__ * acegen_scratch__30__ +
                         acegen_scratch__31__ * acegen_scratch__33__;
  acegen_scratch__318__ = 2e0 * acegen_scratch__36__;
  acegen_scratch__284__ = 2e0 * (acegen_scratch__282__ * acegen_scratch__35__ +
                                 acegen_scratch__280__ * acegen_scratch__36__);
  acegen_scratch__283__ = acegen_scratch__282__ * acegen_scratch__318__;
  acegen_scratch__50__  = (acegen_scratch__36__ * acegen_scratch__36__);
  acegen_scratch__48__  = acegen_scratch__318__ * acegen_scratch__35__;
  acegen_scratch__37__  = (acegen_scratch__26__ * acegen_scratch__26__) +
                         (acegen_scratch__29__ * acegen_scratch__29__) +
                         (acegen_scratch__32__ * acegen_scratch__32__);
  acegen_scratch__322__ = -(acegen_scratch__280__ * acegen_scratch__320__) +
                          acegen_scratch__285__ * acegen_scratch__34__ +
                          acegen_scratch__279__ * acegen_scratch__37__;
  acegen_scratch__143__ =
    acegen_scratch__34__ * acegen_scratch__37__ - acegen_scratch__44__;
  acegen_scratch__38__ = acegen_scratch__26__ * acegen_scratch__27__ +
                         acegen_scratch__29__ * acegen_scratch__30__ +
                         acegen_scratch__32__ * acegen_scratch__33__;
  acegen_scratch__319__ = 2e0 * acegen_scratch__38__;
  acegen_scratch__290__ = acegen_scratch__287__ * acegen_scratch__319__;
  acegen_scratch__179__ = acegen_scratch__319__ * acegen_scratch__35__ -
                          acegen_scratch__318__ * acegen_scratch__37__;
  acegen_scratch__153__ =
    -(acegen_scratch__319__ * acegen_scratch__34__) + acegen_scratch__48__;
  acegen_scratch__54__ = (acegen_scratch__38__ * acegen_scratch__38__);
  acegen_scratch__39__ = (acegen_scratch__27__ * acegen_scratch__27__) +
                         (acegen_scratch__30__ * acegen_scratch__30__) +
                         (acegen_scratch__33__ * acegen_scratch__33__);
  acegen_scratch__214__ =
    acegen_scratch__37__ * acegen_scratch__39__ - acegen_scratch__54__;
  acegen_scratch__195__ = acegen_scratch__319__ * acegen_scratch__36__ -
                          acegen_scratch__320__ * acegen_scratch__39__;
  acegen_scratch__166__ =
    acegen_scratch__34__ * acegen_scratch__39__ - acegen_scratch__50__;
  acegen_scratch__45__ = acegen_scratch__214__ * acegen_scratch__34__ -
                         acegen_scratch__39__ * acegen_scratch__44__ +
                         acegen_scratch__38__ * acegen_scratch__48__ -
                         acegen_scratch__37__ * acegen_scratch__50__;
  acegen_scratch__40__ = sqrt(acegen_scratch__45__);
  acegen_scratch__43__ =
    (25000e0 / 81e0) * (-3e0 + 14e0 * log(acegen_scratch__40__));
  acegen_scratch__300__ =
    ((-81e0 * acegen_scratch__43__ +
      (175000e0 * acegen_scratch__45__) /
        (acegen_scratch__40__ * acegen_scratch__40__)) *
     (acegen_scratch__143__ * acegen_scratch__291__ -
      acegen_scratch__290__ * acegen_scratch__34__ -
      acegen_scratch__283__ * acegen_scratch__37__ +
      acegen_scratch__284__ * acegen_scratch__38__ +
      acegen_scratch__322__ * acegen_scratch__39__ +
      acegen_scratch__287__ * acegen_scratch__48__ -
      acegen_scratch__285__ * acegen_scratch__50__ -
      acegen_scratch__279__ * acegen_scratch__54__)) /
    (162e0 * (acegen_scratch__45__ * acegen_scratch__45__));
  acegen_scratch__47__  = acegen_scratch__43__ / (2e0 * acegen_scratch__45__);
  acegen_scratch__321__ = 2e0 * acegen_scratch__47__;
  acegen_scratch__306__ =
    acegen_scratch__214__ * acegen_scratch__300__ +
    (-acegen_scratch__290__ + acegen_scratch__291__ * acegen_scratch__37__ +
     acegen_scratch__285__ * acegen_scratch__39__) *
      acegen_scratch__47__;
  acegen_scratch__305__ =
    acegen_scratch__195__ * acegen_scratch__300__ +
    acegen_scratch__321__ * (-(acegen_scratch__291__ * acegen_scratch__35__) +
                             acegen_scratch__287__ * acegen_scratch__36__ +
                             acegen_scratch__282__ * acegen_scratch__38__ -
                             acegen_scratch__280__ * acegen_scratch__39__);
  acegen_scratch__304__ =
    acegen_scratch__179__ * acegen_scratch__300__ +
    acegen_scratch__321__ * (acegen_scratch__287__ * acegen_scratch__35__ -
                             acegen_scratch__285__ * acegen_scratch__36__ -
                             acegen_scratch__282__ * acegen_scratch__37__ +
                             acegen_scratch__280__ * acegen_scratch__38__);
  acegen_scratch__303__ =
    acegen_scratch__166__ * acegen_scratch__300__ +
    (-acegen_scratch__283__ + acegen_scratch__291__ * acegen_scratch__34__ +
     acegen_scratch__279__ * acegen_scratch__39__) *
      acegen_scratch__47__;
  acegen_scratch__302__ =
    acegen_scratch__153__ * acegen_scratch__300__ +
    (acegen_scratch__284__ - acegen_scratch__279__ * acegen_scratch__319__ -
     2e0 * acegen_scratch__287__ * acegen_scratch__34__) *
      acegen_scratch__47__;
  acegen_scratch__301__ = acegen_scratch__143__ * acegen_scratch__300__ +
                          acegen_scratch__322__ * acegen_scratch__47__;
  acegen_scratch__267__ =
    (25000e0 / 27e0) + acegen_scratch__143__ * acegen_scratch__321__;
  acegen_scratch__49__ = acegen_scratch__153__ * acegen_scratch__47__;
  acegen_scratch__261__ =
    (25000e0 / 27e0) + acegen_scratch__166__ * acegen_scratch__321__;
  acegen_scratch__52__ = acegen_scratch__179__ * acegen_scratch__47__;
  acegen_scratch__53__ = acegen_scratch__195__ * acegen_scratch__47__;
  acegen_scratch__254__ =
    (25000e0 / 27e0) + acegen_scratch__214__ * acegen_scratch__321__;
  gradientOut[0][0] = acegen_scratch__16__ * acegen_scratch__254__ +
                      acegen_scratch__27__ * acegen_scratch__304__ +
                      acegen_scratch__26__ * acegen_scratch__305__ +
                      acegen_scratch__18__ * acegen_scratch__52__ +
                      acegen_scratch__17__ * acegen_scratch__53__ +
                      acegen_scratch__306__ * acegen_scratch__65__;
  gradientOut[0][1] = acegen_scratch__17__ * acegen_scratch__261__ +
                      acegen_scratch__27__ * acegen_scratch__302__ +
                      acegen_scratch__25__ * acegen_scratch__305__ +
                      acegen_scratch__18__ * acegen_scratch__49__ +
                      acegen_scratch__16__ * acegen_scratch__53__ +
                      acegen_scratch__303__ * acegen_scratch__83__;
  gradientOut[0][2] = acegen_scratch__18__ * acegen_scratch__267__ +
                      acegen_scratch__26__ * acegen_scratch__302__ +
                      acegen_scratch__25__ * acegen_scratch__304__ +
                      acegen_scratch__17__ * acegen_scratch__49__ +
                      acegen_scratch__16__ * acegen_scratch__52__ +
                      acegen_scratch__301__ * acegen_scratch__92__;
  gradientOut[1][0] = acegen_scratch__19__ * acegen_scratch__254__ +
                      acegen_scratch__30__ * acegen_scratch__304__ +
                      acegen_scratch__29__ * acegen_scratch__305__ +
                      acegen_scratch__21__ * acegen_scratch__52__ +
                      acegen_scratch__20__ * acegen_scratch__53__ +
                      acegen_scratch__306__ * acegen_scratch__66__;
  gradientOut[1][1] = acegen_scratch__20__ * acegen_scratch__261__ +
                      acegen_scratch__30__ * acegen_scratch__302__ +
                      acegen_scratch__28__ * acegen_scratch__305__ +
                      acegen_scratch__21__ * acegen_scratch__49__ +
                      acegen_scratch__19__ * acegen_scratch__53__ +
                      acegen_scratch__303__ * acegen_scratch__84__;
  gradientOut[1][2] = acegen_scratch__21__ * acegen_scratch__267__ +
                      acegen_scratch__29__ * acegen_scratch__302__ +
                      acegen_scratch__28__ * acegen_scratch__304__ +
                      acegen_scratch__20__ * acegen_scratch__49__ +
                      acegen_scratch__19__ * acegen_scratch__52__ +
                      acegen_scratch__301__ * acegen_scratch__93__;
  gradientOut[2][0] = acegen_scratch__22__ * acegen_scratch__254__ +
                      acegen_scratch__305__ * acegen_scratch__32__ +
                      acegen_scratch__304__ * acegen_scratch__33__ +
                      acegen_scratch__24__ * acegen_scratch__52__ +
                      acegen_scratch__23__ * acegen_scratch__53__ +
                      acegen_scratch__306__ * acegen_scratch__67__;
  gradientOut[2][1] = acegen_scratch__23__ * acegen_scratch__261__ +
                      acegen_scratch__305__ * acegen_scratch__31__ +
                      acegen_scratch__302__ * acegen_scratch__33__ +
                      acegen_scratch__24__ * acegen_scratch__49__ +
                      acegen_scratch__22__ * acegen_scratch__53__ +
                      acegen_scratch__303__ * acegen_scratch__85__;
  gradientOut[2][2] = acegen_scratch__24__ * acegen_scratch__267__ +
                      acegen_scratch__304__ * acegen_scratch__31__ +
                      acegen_scratch__302__ * acegen_scratch__32__ +
                      acegen_scratch__23__ * acegen_scratch__49__ +
                      acegen_scratch__22__ * acegen_scratch__52__ +
                      acegen_scratch__301__ * acegen_scratch__94__;
}