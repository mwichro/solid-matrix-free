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
struct NeoHookeanCached
{
  const static unsigned int dim      = DIMENSION;
  const static unsigned int n_cached = 27;


  template <typename Number>
  static inline void
  Tangent_local([[maybe_unused]] const Tensor<2, dim, Number> &graduIn,
                const Tensor<2, dim, Number> &                 gradduIn,
                const ArrayView<const Number> &                cacheIn,
                Tensor<2, dim, Number> &                       gradientOut);

  template <typename Number>
  static inline void
  cache_local(const Tensor<2, dim, Number> &graduIn,
              Tensor<1, dim, Number> &      valueOut,
              Tensor<2, dim, Number> &      gradientOut,
              ArrayView<Number> &           cacheOut,
              const Number &                mu,
              const Number &                lambda);
};


// =====
// dim=3
// =====


template <>
template <typename Number>
inline void
NeoHookeanCached<3>::cache_local(const Tensor<2, dim, Number> &graduIn,
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
NeoHookeanCached<3>::Tangent_local(
  [[maybe_unused]] const Tensor<2, dim, Number> &graduIn,
  const Tensor<2, dim, Number> &                 gradduIn,
  const ArrayView<const Number> &                cacheIn,
  Tensor<2, dim, Number> &                       gradientOut)
{
  Number acegen_scratch__107__, acegen_scratch__108__, acegen_scratch__109__,
    acegen_scratch__121__, acegen_scratch__122__, acegen_scratch__123__,
    acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__18__,
    acegen_scratch__19__, acegen_scratch__20__, acegen_scratch__21__,
    acegen_scratch__22__, acegen_scratch__23__, acegen_scratch__24__,
    acegen_scratch__26__, acegen_scratch__27__, acegen_scratch__28__,
    acegen_scratch__29__, acegen_scratch__30__, acegen_scratch__32__,
    acegen_scratch__33__, acegen_scratch__34__, acegen_scratch__35__,
    acegen_scratch__37__, acegen_scratch__38__, acegen_scratch__39__,
    acegen_scratch__41__, acegen_scratch__42__, acegen_scratch__44__,
    acegen_scratch__46__, acegen_scratch__47__, acegen_scratch__48__,
    acegen_scratch__49__, acegen_scratch__50__, acegen_scratch__51__,
    acegen_scratch__52__, acegen_scratch__53__, acegen_scratch__54__,
    acegen_scratch__55__, acegen_scratch__56__, acegen_scratch__57__,
    acegen_scratch__58__, acegen_scratch__59__, acegen_scratch__60__,
    acegen_scratch__61__, acegen_scratch__65__, acegen_scratch__69__,
    acegen_scratch__70__, acegen_scratch__71__, acegen_scratch__72__;
  acegen_scratch__16__ = gradduIn[0][0];
  acegen_scratch__17__ = gradduIn[0][1];
  acegen_scratch__18__ = gradduIn[0][2];
  acegen_scratch__19__ = gradduIn[1][0];
  acegen_scratch__20__ = gradduIn[1][1];
  acegen_scratch__21__ = gradduIn[1][2];
  acegen_scratch__22__ = gradduIn[2][0];
  acegen_scratch__23__ = gradduIn[2][1];
  acegen_scratch__24__ = gradduIn[2][2];
  acegen_scratch__26__ = cacheIn[1];
  acegen_scratch__27__ = cacheIn[2];
  acegen_scratch__28__ = cacheIn[3];
  acegen_scratch__29__ = cacheIn[4];
  acegen_scratch__30__ = cacheIn[5];
  acegen_scratch__32__ = cacheIn[7];
  acegen_scratch__33__ = cacheIn[8];
  acegen_scratch__34__ = cacheIn[9];
  acegen_scratch__35__ = cacheIn[10];
  acegen_scratch__37__ = cacheIn[12];
  acegen_scratch__38__ = cacheIn[13];
  acegen_scratch__39__ = cacheIn[14];
  acegen_scratch__41__ = cacheIn[16];
  acegen_scratch__42__ = cacheIn[17];
  acegen_scratch__44__ = cacheIn[19];
  acegen_scratch__46__ = cacheIn[21];
  acegen_scratch__47__ = cacheIn[22];
  acegen_scratch__48__ = cacheIn[23];
  acegen_scratch__49__ = cacheIn[24];
  acegen_scratch__50__ = cacheIn[25];
  acegen_scratch__51__ = cacheIn[26];
  acegen_scratch__52__ = 1e0 + graduIn[0][0];
  acegen_scratch__53__ = graduIn[0][1];
  acegen_scratch__54__ = graduIn[0][2];
  acegen_scratch__55__ = graduIn[1][0];
  acegen_scratch__56__ = 1e0 + graduIn[1][1];
  acegen_scratch__57__ = graduIn[1][2];
  acegen_scratch__58__ = graduIn[2][0];
  acegen_scratch__59__ = graduIn[2][1];
  acegen_scratch__60__ = 1e0 + graduIn[2][2];
  acegen_scratch__61__ = acegen_scratch__16__ * acegen_scratch__52__ +
                         acegen_scratch__19__ * acegen_scratch__55__ +
                         acegen_scratch__22__ * acegen_scratch__58__;
  acegen_scratch__65__ = acegen_scratch__17__ * acegen_scratch__53__ +
                         acegen_scratch__20__ * acegen_scratch__56__ +
                         acegen_scratch__23__ * acegen_scratch__59__;
  acegen_scratch__69__ = acegen_scratch__18__ * acegen_scratch__54__ +
                         acegen_scratch__21__ * acegen_scratch__57__ +
                         acegen_scratch__24__ * acegen_scratch__60__;
  acegen_scratch__70__ =
    0.7071067811865476e0 * (acegen_scratch__18__ * acegen_scratch__53__ +
                            acegen_scratch__17__ * acegen_scratch__54__ +
                            acegen_scratch__21__ * acegen_scratch__56__ +
                            acegen_scratch__20__ * acegen_scratch__57__ +
                            acegen_scratch__24__ * acegen_scratch__59__ +
                            acegen_scratch__23__ * acegen_scratch__60__);
  acegen_scratch__71__ =
    0.7071067811865476e0 * (acegen_scratch__18__ * acegen_scratch__52__ +
                            acegen_scratch__16__ * acegen_scratch__54__ +
                            acegen_scratch__21__ * acegen_scratch__55__ +
                            acegen_scratch__19__ * acegen_scratch__57__ +
                            acegen_scratch__24__ * acegen_scratch__58__ +
                            acegen_scratch__22__ * acegen_scratch__60__);
  acegen_scratch__72__ =
    0.7071067811865476e0 * (acegen_scratch__17__ * acegen_scratch__52__ +
                            acegen_scratch__16__ * acegen_scratch__53__ +
                            acegen_scratch__20__ * acegen_scratch__55__ +
                            acegen_scratch__19__ * acegen_scratch__56__ +
                            acegen_scratch__23__ * acegen_scratch__58__ +
                            acegen_scratch__22__ * acegen_scratch__59__);
  acegen_scratch__121__ =
    0.7071067811865476e0 * (acegen_scratch__30__ * acegen_scratch__61__ +
                            acegen_scratch__35__ * acegen_scratch__65__ +
                            acegen_scratch__39__ * acegen_scratch__69__ +
                            acegen_scratch__42__ * acegen_scratch__70__ +
                            acegen_scratch__44__ * acegen_scratch__71__ +
                            acegen_scratch__72__ * cacheIn[20]);
  acegen_scratch__122__ =
    0.7071067811865476e0 * (acegen_scratch__29__ * acegen_scratch__61__ +
                            acegen_scratch__34__ * acegen_scratch__65__ +
                            acegen_scratch__38__ * acegen_scratch__69__ +
                            acegen_scratch__41__ * acegen_scratch__70__ +
                            acegen_scratch__44__ * acegen_scratch__72__ +
                            acegen_scratch__71__ * cacheIn[18]);
  acegen_scratch__123__ =
    0.7071067811865476e0 * (acegen_scratch__28__ * acegen_scratch__61__ +
                            acegen_scratch__33__ * acegen_scratch__65__ +
                            acegen_scratch__37__ * acegen_scratch__69__ +
                            acegen_scratch__41__ * acegen_scratch__71__ +
                            acegen_scratch__42__ * acegen_scratch__72__ +
                            acegen_scratch__70__ * cacheIn[15]);
  acegen_scratch__107__ = acegen_scratch__27__ * acegen_scratch__61__ +
                          acegen_scratch__32__ * acegen_scratch__65__ +
                          acegen_scratch__37__ * acegen_scratch__70__ +
                          acegen_scratch__38__ * acegen_scratch__71__ +
                          acegen_scratch__39__ * acegen_scratch__72__ +
                          acegen_scratch__69__ * cacheIn[11];
  acegen_scratch__108__ = acegen_scratch__26__ * acegen_scratch__61__ +
                          acegen_scratch__32__ * acegen_scratch__69__ +
                          acegen_scratch__33__ * acegen_scratch__70__ +
                          acegen_scratch__34__ * acegen_scratch__71__ +
                          acegen_scratch__35__ * acegen_scratch__72__ +
                          acegen_scratch__65__ * cacheIn[6];
  acegen_scratch__109__ = acegen_scratch__26__ * acegen_scratch__65__ +
                          acegen_scratch__27__ * acegen_scratch__69__ +
                          acegen_scratch__28__ * acegen_scratch__70__ +
                          acegen_scratch__29__ * acegen_scratch__71__ +
                          acegen_scratch__30__ * acegen_scratch__72__ +
                          acegen_scratch__61__ * cacheIn[0];
  gradientOut[0][0] = acegen_scratch__16__ * acegen_scratch__46__ +
                      acegen_scratch__17__ * acegen_scratch__47__ +
                      acegen_scratch__18__ * acegen_scratch__48__ +
                      acegen_scratch__109__ * acegen_scratch__52__ +
                      acegen_scratch__121__ * acegen_scratch__53__ +
                      acegen_scratch__122__ * acegen_scratch__54__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__47__ +
                      acegen_scratch__17__ * acegen_scratch__49__ +
                      acegen_scratch__18__ * acegen_scratch__50__ +
                      acegen_scratch__121__ * acegen_scratch__52__ +
                      acegen_scratch__108__ * acegen_scratch__53__ +
                      acegen_scratch__123__ * acegen_scratch__54__;
  gradientOut[0][2] = acegen_scratch__16__ * acegen_scratch__48__ +
                      acegen_scratch__17__ * acegen_scratch__50__ +
                      acegen_scratch__18__ * acegen_scratch__51__ +
                      acegen_scratch__122__ * acegen_scratch__52__ +
                      acegen_scratch__123__ * acegen_scratch__53__ +
                      acegen_scratch__107__ * acegen_scratch__54__;
  gradientOut[1][0] = acegen_scratch__19__ * acegen_scratch__46__ +
                      acegen_scratch__20__ * acegen_scratch__47__ +
                      acegen_scratch__21__ * acegen_scratch__48__ +
                      acegen_scratch__109__ * acegen_scratch__55__ +
                      acegen_scratch__121__ * acegen_scratch__56__ +
                      acegen_scratch__122__ * acegen_scratch__57__;
  gradientOut[1][1] = acegen_scratch__19__ * acegen_scratch__47__ +
                      acegen_scratch__20__ * acegen_scratch__49__ +
                      acegen_scratch__21__ * acegen_scratch__50__ +
                      acegen_scratch__121__ * acegen_scratch__55__ +
                      acegen_scratch__108__ * acegen_scratch__56__ +
                      acegen_scratch__123__ * acegen_scratch__57__;
  gradientOut[1][2] = acegen_scratch__19__ * acegen_scratch__48__ +
                      acegen_scratch__20__ * acegen_scratch__50__ +
                      acegen_scratch__21__ * acegen_scratch__51__ +
                      acegen_scratch__122__ * acegen_scratch__55__ +
                      acegen_scratch__123__ * acegen_scratch__56__ +
                      acegen_scratch__107__ * acegen_scratch__57__;
  gradientOut[2][0] = acegen_scratch__22__ * acegen_scratch__46__ +
                      acegen_scratch__23__ * acegen_scratch__47__ +
                      acegen_scratch__24__ * acegen_scratch__48__ +
                      acegen_scratch__109__ * acegen_scratch__58__ +
                      1e0 * acegen_scratch__121__ * acegen_scratch__59__ +
                      acegen_scratch__122__ * acegen_scratch__60__;
  gradientOut[2][1] = acegen_scratch__22__ * acegen_scratch__47__ +
                      acegen_scratch__23__ * acegen_scratch__49__ +
                      acegen_scratch__24__ * acegen_scratch__50__ +
                      acegen_scratch__121__ * acegen_scratch__58__ +
                      acegen_scratch__108__ * acegen_scratch__59__ +
                      1e0 * acegen_scratch__123__ * acegen_scratch__60__;
  gradientOut[2][2] = acegen_scratch__22__ * acegen_scratch__48__ +
                      acegen_scratch__23__ * acegen_scratch__50__ +
                      acegen_scratch__24__ * acegen_scratch__51__ +
                      acegen_scratch__122__ * acegen_scratch__58__ +
                      acegen_scratch__123__ * acegen_scratch__59__ +
                      acegen_scratch__107__ * acegen_scratch__60__;
}



// =====
// dim=2
// =====

// template <>
// template <typename Number>
// inline void
// NeoHookeanCached<2>::Residual_local(
//   [[maybe_unused]] AlignedVector<Number> &       acegen_scratch,
//   [[maybe_unused]] const Tensor<1, dim, Number> &uIn,
//   const Tensor<2, dim, Number> &                 graduIn,
//   Tensor<1, dim, Number> &                       valueOut,
//   Tensor<2, dim, Number> &                       gradientOut)
// {}

template <>
template <typename Number>
inline void
NeoHookeanCached<2>::Tangent_local(
  [[maybe_unused]] const Tensor<2, dim, Number> &graduIn,
  const Tensor<2, dim, Number> &                 gradduIn,
  const ArrayView<const Number> &                cacheIn,
  Tensor<2, dim, Number> &                       gradientOut)
{
  Number acegen_scratch__10__, acegen_scratch__11__, acegen_scratch__12__,
    acegen_scratch__14__, acegen_scratch__15__, acegen_scratch__17__,
    acegen_scratch__19__, acegen_scratch__20__, acegen_scratch__21__,
    acegen_scratch__22__, acegen_scratch__23__, acegen_scratch__24__,
    acegen_scratch__25__, acegen_scratch__30__, acegen_scratch__33__,
    acegen_scratch__34__, acegen_scratch__35__, acegen_scratch__36__,
    acegen_scratch__38__, acegen_scratch__9__;
  acegen_scratch__9__  = gradduIn[0][0];
  acegen_scratch__10__ = gradduIn[0][1];
  acegen_scratch__11__ = gradduIn[1][0];
  acegen_scratch__12__ = gradduIn[1][1];
  acegen_scratch__14__ = cacheIn[1];
  acegen_scratch__15__ = cacheIn[2];
  acegen_scratch__17__ = cacheIn[4];
  acegen_scratch__19__ = cacheIn[6];
  acegen_scratch__20__ = cacheIn[7];
  acegen_scratch__21__ = cacheIn[8];
  acegen_scratch__22__ = 1e0 + graduIn[0][0];
  acegen_scratch__23__ = graduIn[0][1];
  acegen_scratch__24__ = graduIn[1][0];
  acegen_scratch__25__ = 1e0 + graduIn[1][1];
  acegen_scratch__30__ = acegen_scratch__11__ * acegen_scratch__24__ +
                         acegen_scratch__22__ * acegen_scratch__9__;
  acegen_scratch__33__ = acegen_scratch__10__ * acegen_scratch__23__ +
                         acegen_scratch__12__ * acegen_scratch__25__;
  acegen_scratch__34__ =
    0.7071067811865476e0 * (acegen_scratch__10__ * acegen_scratch__22__ +
                            acegen_scratch__12__ * acegen_scratch__24__ +
                            acegen_scratch__11__ * acegen_scratch__25__ +
                            acegen_scratch__23__ * acegen_scratch__9__);
  acegen_scratch__35__ = acegen_scratch__14__ * acegen_scratch__33__ +
                         acegen_scratch__15__ * acegen_scratch__34__ +
                         acegen_scratch__30__ * cacheIn[0];
  acegen_scratch__36__ = acegen_scratch__14__ * acegen_scratch__30__ +
                         acegen_scratch__17__ * acegen_scratch__34__ +
                         acegen_scratch__33__ * cacheIn[3];
  acegen_scratch__38__ =
    0.7071067811865476e0 * (acegen_scratch__15__ * acegen_scratch__30__ +
                            acegen_scratch__17__ * acegen_scratch__33__ +
                            acegen_scratch__34__ * cacheIn[5]);
  gradientOut[0][0] = acegen_scratch__10__ * acegen_scratch__20__ +
                      acegen_scratch__22__ * acegen_scratch__35__ +
                      acegen_scratch__23__ * acegen_scratch__38__ +
                      acegen_scratch__19__ * acegen_scratch__9__;
  gradientOut[0][1] = acegen_scratch__10__ * acegen_scratch__21__ +
                      acegen_scratch__23__ * acegen_scratch__36__ +
                      acegen_scratch__22__ * acegen_scratch__38__ +
                      acegen_scratch__20__ * acegen_scratch__9__;
  gradientOut[1][0] = acegen_scratch__11__ * acegen_scratch__19__ +
                      acegen_scratch__12__ * acegen_scratch__20__ +
                      acegen_scratch__24__ * acegen_scratch__35__ +
                      acegen_scratch__25__ * acegen_scratch__38__;
  gradientOut[1][1] = acegen_scratch__11__ * acegen_scratch__20__ +
                      acegen_scratch__12__ * acegen_scratch__21__ +
                      acegen_scratch__25__ * acegen_scratch__36__ +
                      acegen_scratch__24__ * acegen_scratch__38__;
}


template <>
template <typename Number>
inline void
NeoHookeanCached<2>::cache_local(const Tensor<2, dim, Number> &graduIn,
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
