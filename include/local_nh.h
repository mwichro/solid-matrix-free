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
                Tensor<2, dim, Number> &      gradientOut,
                const Number &                mu,
                const Number &                lambda);
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
                             Tensor<2, dim, Number> &      gradientOut,
                             const Number &                mu,
                             const Number &                lambda)
{
  Number acegen_scratch__146__, acegen_scratch__156__, acegen_scratch__16__,
    acegen_scratch__169__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__182__, acegen_scratch__19__, acegen_scratch__198__,
    acegen_scratch__20__, acegen_scratch__21__, acegen_scratch__217__,
    acegen_scratch__22__, acegen_scratch__23__, acegen_scratch__24__,
    acegen_scratch__25__, acegen_scratch__257__, acegen_scratch__26__,
    acegen_scratch__264__, acegen_scratch__27__, acegen_scratch__270__,
    acegen_scratch__28__, acegen_scratch__282__, acegen_scratch__283__,
    acegen_scratch__285__, acegen_scratch__286__, acegen_scratch__287__,
    acegen_scratch__288__, acegen_scratch__29__, acegen_scratch__290__,
    acegen_scratch__293__, acegen_scratch__294__, acegen_scratch__30__,
    acegen_scratch__303__, acegen_scratch__304__, acegen_scratch__305__,
    acegen_scratch__306__, acegen_scratch__307__, acegen_scratch__308__,
    acegen_scratch__309__, acegen_scratch__31__, acegen_scratch__32__,
    acegen_scratch__321__, acegen_scratch__322__, acegen_scratch__323__,
    acegen_scratch__324__, acegen_scratch__325__, acegen_scratch__33__,
    acegen_scratch__34__, acegen_scratch__35__, acegen_scratch__36__,
    acegen_scratch__37__, acegen_scratch__38__, acegen_scratch__39__,
    acegen_scratch__40__, acegen_scratch__41__, acegen_scratch__42__,
    acegen_scratch__45__, acegen_scratch__46__, acegen_scratch__47__,
    acegen_scratch__49__, acegen_scratch__50__, acegen_scratch__51__,
    acegen_scratch__52__, acegen_scratch__53__, acegen_scratch__55__,
    acegen_scratch__56__, acegen_scratch__57__, acegen_scratch__68__,
    acegen_scratch__69__, acegen_scratch__70__, acegen_scratch__86__,
    acegen_scratch__87__, acegen_scratch__88__, acegen_scratch__95__,
    acegen_scratch__96__, acegen_scratch__97__;
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
  acegen_scratch__52__  = acegen_scratch__25__ / 2e0;
  acegen_scratch__26__  = lambda;
  acegen_scratch__27__  = 1e0 + graduIn[0][0];
  acegen_scratch__68__  = 2e0 * acegen_scratch__27__;
  acegen_scratch__28__  = graduIn[0][1];
  acegen_scratch__86__  = 2e0 * acegen_scratch__28__;
  acegen_scratch__29__  = graduIn[0][2];
  acegen_scratch__95__  = 2e0 * acegen_scratch__29__;
  acegen_scratch__30__  = graduIn[1][0];
  acegen_scratch__69__  = 2e0 * acegen_scratch__30__;
  acegen_scratch__31__  = 1e0 + graduIn[1][1];
  acegen_scratch__87__  = 2e0 * acegen_scratch__31__;
  acegen_scratch__32__  = graduIn[1][2];
  acegen_scratch__96__  = 2e0 * acegen_scratch__32__;
  acegen_scratch__33__  = graduIn[2][0];
  acegen_scratch__70__  = 2e0 * acegen_scratch__33__;
  acegen_scratch__282__ = acegen_scratch__16__ * acegen_scratch__68__ +
                          acegen_scratch__19__ * acegen_scratch__69__ +
                          acegen_scratch__22__ * acegen_scratch__70__;
  acegen_scratch__34__  = graduIn[2][1];
  acegen_scratch__283__ = acegen_scratch__17__ * acegen_scratch__27__ +
                          acegen_scratch__16__ * acegen_scratch__28__ +
                          acegen_scratch__20__ * acegen_scratch__30__ +
                          acegen_scratch__19__ * acegen_scratch__31__ +
                          acegen_scratch__23__ * acegen_scratch__33__ +
                          acegen_scratch__22__ * acegen_scratch__34__;
  acegen_scratch__88__  = 2e0 * acegen_scratch__34__;
  acegen_scratch__288__ = acegen_scratch__17__ * acegen_scratch__86__ +
                          acegen_scratch__20__ * acegen_scratch__87__ +
                          acegen_scratch__23__ * acegen_scratch__88__;
  acegen_scratch__35__  = 1e0 + graduIn[2][2];
  acegen_scratch__290__ = acegen_scratch__18__ * acegen_scratch__28__ +
                          acegen_scratch__17__ * acegen_scratch__29__ +
                          acegen_scratch__21__ * acegen_scratch__31__ +
                          acegen_scratch__20__ * acegen_scratch__32__ +
                          acegen_scratch__24__ * acegen_scratch__34__ +
                          acegen_scratch__23__ * acegen_scratch__35__;
  acegen_scratch__285__ = acegen_scratch__18__ * acegen_scratch__27__ +
                          acegen_scratch__16__ * acegen_scratch__29__ +
                          acegen_scratch__21__ * acegen_scratch__30__ +
                          acegen_scratch__19__ * acegen_scratch__32__ +
                          acegen_scratch__24__ * acegen_scratch__33__ +
                          acegen_scratch__22__ * acegen_scratch__35__;
  acegen_scratch__97__  = 2e0 * acegen_scratch__35__;
  acegen_scratch__294__ = acegen_scratch__18__ * acegen_scratch__95__ +
                          acegen_scratch__21__ * acegen_scratch__96__ +
                          acegen_scratch__24__ * acegen_scratch__97__;
  acegen_scratch__36__ = (acegen_scratch__27__ * acegen_scratch__27__) +
                         (acegen_scratch__30__ * acegen_scratch__30__) +
                         (acegen_scratch__33__ * acegen_scratch__33__);
  acegen_scratch__37__ = acegen_scratch__27__ * acegen_scratch__28__ +
                         acegen_scratch__30__ * acegen_scratch__31__ +
                         acegen_scratch__33__ * acegen_scratch__34__;
  acegen_scratch__323__ = 2e0 * acegen_scratch__37__;
  acegen_scratch__46__  = (acegen_scratch__37__ * acegen_scratch__37__);
  acegen_scratch__38__  = acegen_scratch__27__ * acegen_scratch__29__ +
                         acegen_scratch__30__ * acegen_scratch__32__ +
                         acegen_scratch__33__ * acegen_scratch__35__;
  acegen_scratch__321__ = 2e0 * acegen_scratch__38__;
  acegen_scratch__287__ = 2e0 * (acegen_scratch__285__ * acegen_scratch__37__ +
                                 acegen_scratch__283__ * acegen_scratch__38__);
  acegen_scratch__286__ = acegen_scratch__285__ * acegen_scratch__321__;
  acegen_scratch__53__  = (acegen_scratch__38__ * acegen_scratch__38__);
  acegen_scratch__50__  = acegen_scratch__321__ * acegen_scratch__37__;
  acegen_scratch__39__  = (acegen_scratch__28__ * acegen_scratch__28__) +
                         (acegen_scratch__31__ * acegen_scratch__31__) +
                         (acegen_scratch__34__ * acegen_scratch__34__);
  acegen_scratch__325__ = -(acegen_scratch__283__ * acegen_scratch__323__) +
                          acegen_scratch__288__ * acegen_scratch__36__ +
                          acegen_scratch__282__ * acegen_scratch__39__;
  acegen_scratch__146__ =
    acegen_scratch__36__ * acegen_scratch__39__ - acegen_scratch__46__;
  acegen_scratch__40__ = acegen_scratch__28__ * acegen_scratch__29__ +
                         acegen_scratch__31__ * acegen_scratch__32__ +
                         acegen_scratch__34__ * acegen_scratch__35__;
  acegen_scratch__322__ = 2e0 * acegen_scratch__40__;
  acegen_scratch__293__ = acegen_scratch__290__ * acegen_scratch__322__;
  acegen_scratch__182__ = acegen_scratch__322__ * acegen_scratch__37__ -
                          acegen_scratch__321__ * acegen_scratch__39__;
  acegen_scratch__156__ =
    -(acegen_scratch__322__ * acegen_scratch__36__) + acegen_scratch__50__;
  acegen_scratch__57__ = (acegen_scratch__40__ * acegen_scratch__40__);
  acegen_scratch__41__ = (acegen_scratch__29__ * acegen_scratch__29__) +
                         (acegen_scratch__32__ * acegen_scratch__32__) +
                         (acegen_scratch__35__ * acegen_scratch__35__);
  acegen_scratch__217__ =
    acegen_scratch__39__ * acegen_scratch__41__ - acegen_scratch__57__;
  acegen_scratch__198__ = acegen_scratch__322__ * acegen_scratch__38__ -
                          acegen_scratch__323__ * acegen_scratch__41__;
  acegen_scratch__169__ =
    acegen_scratch__36__ * acegen_scratch__41__ - acegen_scratch__53__;
  acegen_scratch__47__ = acegen_scratch__217__ * acegen_scratch__36__ -
                         acegen_scratch__41__ * acegen_scratch__46__ +
                         acegen_scratch__40__ * acegen_scratch__50__ -
                         acegen_scratch__39__ * acegen_scratch__53__;
  acegen_scratch__42__ = sqrt(acegen_scratch__47__);
  acegen_scratch__45__ = -acegen_scratch__25__ +
                         2e0 * acegen_scratch__26__ * log(acegen_scratch__42__);
  acegen_scratch__303__ =
    ((-acegen_scratch__45__ + (acegen_scratch__26__ * acegen_scratch__47__) /
                                (acegen_scratch__42__ * acegen_scratch__42__)) *
     (acegen_scratch__146__ * acegen_scratch__294__ -
      acegen_scratch__293__ * acegen_scratch__36__ -
      acegen_scratch__286__ * acegen_scratch__39__ +
      acegen_scratch__287__ * acegen_scratch__40__ +
      acegen_scratch__325__ * acegen_scratch__41__ +
      acegen_scratch__290__ * acegen_scratch__50__ -
      acegen_scratch__288__ * acegen_scratch__53__ -
      acegen_scratch__282__ * acegen_scratch__57__)) /
    (2e0 * (acegen_scratch__47__ * acegen_scratch__47__));
  acegen_scratch__49__  = acegen_scratch__45__ / (2e0 * acegen_scratch__47__);
  acegen_scratch__324__ = 2e0 * acegen_scratch__49__;
  acegen_scratch__309__ =
    acegen_scratch__217__ * acegen_scratch__303__ +
    (-acegen_scratch__293__ + acegen_scratch__294__ * acegen_scratch__39__ +
     acegen_scratch__288__ * acegen_scratch__41__) *
      acegen_scratch__49__;
  acegen_scratch__308__ =
    acegen_scratch__198__ * acegen_scratch__303__ +
    acegen_scratch__324__ * (-(acegen_scratch__294__ * acegen_scratch__37__) +
                             acegen_scratch__290__ * acegen_scratch__38__ +
                             acegen_scratch__285__ * acegen_scratch__40__ -
                             acegen_scratch__283__ * acegen_scratch__41__);
  acegen_scratch__307__ =
    acegen_scratch__182__ * acegen_scratch__303__ +
    acegen_scratch__324__ * (acegen_scratch__290__ * acegen_scratch__37__ -
                             acegen_scratch__288__ * acegen_scratch__38__ -
                             acegen_scratch__285__ * acegen_scratch__39__ +
                             acegen_scratch__283__ * acegen_scratch__40__);
  acegen_scratch__306__ =
    acegen_scratch__169__ * acegen_scratch__303__ +
    (-acegen_scratch__286__ + acegen_scratch__294__ * acegen_scratch__36__ +
     acegen_scratch__282__ * acegen_scratch__41__) *
      acegen_scratch__49__;
  acegen_scratch__305__ =
    acegen_scratch__156__ * acegen_scratch__303__ +
    (acegen_scratch__287__ - acegen_scratch__282__ * acegen_scratch__322__ -
     2e0 * acegen_scratch__290__ * acegen_scratch__36__) *
      acegen_scratch__49__;
  acegen_scratch__304__ = acegen_scratch__146__ * acegen_scratch__303__ +
                          acegen_scratch__325__ * acegen_scratch__49__;
  acegen_scratch__270__ =
    2e0 * (acegen_scratch__146__ * acegen_scratch__49__ + acegen_scratch__52__);
  acegen_scratch__51__ = acegen_scratch__156__ * acegen_scratch__49__;
  acegen_scratch__264__ =
    2e0 * (acegen_scratch__169__ * acegen_scratch__49__ + acegen_scratch__52__);
  acegen_scratch__55__ = acegen_scratch__182__ * acegen_scratch__49__;
  acegen_scratch__56__ = acegen_scratch__198__ * acegen_scratch__49__;
  acegen_scratch__257__ =
    2e0 * (acegen_scratch__217__ * acegen_scratch__49__ + acegen_scratch__52__);
  gradientOut[0][0] = acegen_scratch__16__ * acegen_scratch__257__ +
                      acegen_scratch__29__ * acegen_scratch__307__ +
                      acegen_scratch__28__ * acegen_scratch__308__ +
                      acegen_scratch__18__ * acegen_scratch__55__ +
                      acegen_scratch__17__ * acegen_scratch__56__ +
                      acegen_scratch__309__ * acegen_scratch__68__;
  gradientOut[0][1] = acegen_scratch__17__ * acegen_scratch__264__ +
                      acegen_scratch__29__ * acegen_scratch__305__ +
                      acegen_scratch__27__ * acegen_scratch__308__ +
                      acegen_scratch__18__ * acegen_scratch__51__ +
                      acegen_scratch__16__ * acegen_scratch__56__ +
                      acegen_scratch__306__ * acegen_scratch__86__;
  gradientOut[0][2] = acegen_scratch__18__ * acegen_scratch__270__ +
                      acegen_scratch__28__ * acegen_scratch__305__ +
                      acegen_scratch__27__ * acegen_scratch__307__ +
                      acegen_scratch__17__ * acegen_scratch__51__ +
                      acegen_scratch__16__ * acegen_scratch__55__ +
                      acegen_scratch__304__ * acegen_scratch__95__;
  gradientOut[1][0] = acegen_scratch__19__ * acegen_scratch__257__ +
                      acegen_scratch__308__ * acegen_scratch__31__ +
                      acegen_scratch__307__ * acegen_scratch__32__ +
                      acegen_scratch__21__ * acegen_scratch__55__ +
                      acegen_scratch__20__ * acegen_scratch__56__ +
                      acegen_scratch__309__ * acegen_scratch__69__;
  gradientOut[1][1] = acegen_scratch__20__ * acegen_scratch__264__ +
                      acegen_scratch__30__ * acegen_scratch__308__ +
                      acegen_scratch__305__ * acegen_scratch__32__ +
                      acegen_scratch__21__ * acegen_scratch__51__ +
                      acegen_scratch__19__ * acegen_scratch__56__ +
                      acegen_scratch__306__ * acegen_scratch__87__;
  gradientOut[1][2] = acegen_scratch__21__ * acegen_scratch__270__ +
                      acegen_scratch__30__ * acegen_scratch__307__ +
                      acegen_scratch__305__ * acegen_scratch__31__ +
                      acegen_scratch__20__ * acegen_scratch__51__ +
                      acegen_scratch__19__ * acegen_scratch__55__ +
                      acegen_scratch__304__ * acegen_scratch__96__;
  gradientOut[2][0] = acegen_scratch__22__ * acegen_scratch__257__ +
                      acegen_scratch__308__ * acegen_scratch__34__ +
                      acegen_scratch__307__ * acegen_scratch__35__ +
                      acegen_scratch__24__ * acegen_scratch__55__ +
                      acegen_scratch__23__ * acegen_scratch__56__ +
                      acegen_scratch__309__ * acegen_scratch__70__;
  gradientOut[2][1] = acegen_scratch__23__ * acegen_scratch__264__ +
                      acegen_scratch__308__ * acegen_scratch__33__ +
                      acegen_scratch__305__ * acegen_scratch__35__ +
                      acegen_scratch__24__ * acegen_scratch__51__ +
                      acegen_scratch__22__ * acegen_scratch__56__ +
                      acegen_scratch__306__ * acegen_scratch__88__;
  gradientOut[2][2] = acegen_scratch__24__ * acegen_scratch__270__ +
                      acegen_scratch__307__ * acegen_scratch__33__ +
                      acegen_scratch__305__ * acegen_scratch__34__ +
                      acegen_scratch__23__ * acegen_scratch__51__ +
                      acegen_scratch__22__ * acegen_scratch__55__ +
                      acegen_scratch__304__ * acegen_scratch__97__;
}