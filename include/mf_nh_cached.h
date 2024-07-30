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
              ArrayView<Number> &           cacheOut);
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
                                 ArrayView<Number> &           cacheOut)
{
  Number acegen_scratch__103__, acegen_scratch__110__, acegen_scratch__111__,
    acegen_scratch__113__, acegen_scratch__114__, acegen_scratch__115__,
    acegen_scratch__116__, acegen_scratch__117__, acegen_scratch__13__,
    acegen_scratch__133__, acegen_scratch__137__, acegen_scratch__138__,
    acegen_scratch__14__, acegen_scratch__141__, acegen_scratch__15__,
    acegen_scratch__16__, acegen_scratch__167__, acegen_scratch__168__,
    acegen_scratch__169__, acegen_scratch__17__, acegen_scratch__170__,
    acegen_scratch__171__, acegen_scratch__172__, acegen_scratch__173__,
    acegen_scratch__174__, acegen_scratch__175__, acegen_scratch__176__,
    acegen_scratch__177__, acegen_scratch__178__, acegen_scratch__179__,
    acegen_scratch__18__, acegen_scratch__180__, acegen_scratch__181__,
    acegen_scratch__182__, acegen_scratch__183__, acegen_scratch__184__,
    acegen_scratch__185__, acegen_scratch__19__, acegen_scratch__20__,
    acegen_scratch__21__, acegen_scratch__28__, acegen_scratch__29__,
    acegen_scratch__30__, acegen_scratch__31__, acegen_scratch__32__,
    acegen_scratch__33__, acegen_scratch__40__, acegen_scratch__42__,
    acegen_scratch__43__, acegen_scratch__44__, acegen_scratch__45__,
    acegen_scratch__46__, acegen_scratch__47__, acegen_scratch__48__,
    acegen_scratch__50__, acegen_scratch__53__, acegen_scratch__56__,
    acegen_scratch__57__, acegen_scratch__58__, acegen_scratch__74__,
    acegen_scratch__75__, acegen_scratch__76__, acegen_scratch__77__,
    acegen_scratch__78__, acegen_scratch__79__, acegen_scratch__81__,
    acegen_scratch__84__;
  acegen_scratch__13__  = 1e0 + graduIn[0][0];
  acegen_scratch__14__  = graduIn[0][1];
  acegen_scratch__15__  = graduIn[0][2];
  acegen_scratch__16__  = graduIn[1][0];
  acegen_scratch__17__  = 1e0 + graduIn[1][1];
  acegen_scratch__18__  = graduIn[1][2];
  acegen_scratch__19__  = graduIn[2][0];
  acegen_scratch__20__  = graduIn[2][1];
  acegen_scratch__167__ = acegen_scratch__13__ * acegen_scratch__14__ +
                          acegen_scratch__16__ * acegen_scratch__17__ +
                          acegen_scratch__19__ * acegen_scratch__20__;
  acegen_scratch__21__  = 1e0 + graduIn[2][2];
  acegen_scratch__169__ = acegen_scratch__14__ * acegen_scratch__15__ +
                          acegen_scratch__17__ * acegen_scratch__18__ +
                          acegen_scratch__20__ * acegen_scratch__21__;
  acegen_scratch__168__ = acegen_scratch__13__ * acegen_scratch__15__ +
                          acegen_scratch__16__ * acegen_scratch__18__ +
                          acegen_scratch__19__ * acegen_scratch__21__;
  acegen_scratch__28__ = (acegen_scratch__13__ * acegen_scratch__13__) +
                         (acegen_scratch__16__ * acegen_scratch__16__) +
                         (acegen_scratch__19__ * acegen_scratch__19__);
  acegen_scratch__29__ = (acegen_scratch__14__ * acegen_scratch__14__) +
                         (acegen_scratch__17__ * acegen_scratch__17__) +
                         (acegen_scratch__20__ * acegen_scratch__20__);
  acegen_scratch__30__ = (acegen_scratch__15__ * acegen_scratch__15__) +
                         (acegen_scratch__18__ * acegen_scratch__18__) +
                         (acegen_scratch__21__ * acegen_scratch__21__);
  acegen_scratch__31__  = 0.1414213562373095e1 * acegen_scratch__169__;
  acegen_scratch__32__  = 0.1414213562373095e1 * acegen_scratch__168__;
  acegen_scratch__33__  = 0.1414213562373095e1 * acegen_scratch__167__;
  acegen_scratch__170__ = 2e0 * acegen_scratch__167__;
  acegen_scratch__47__  = (acegen_scratch__167__ * acegen_scratch__167__);
  acegen_scratch__76__ =
    acegen_scratch__28__ * acegen_scratch__29__ - acegen_scratch__47__;
  acegen_scratch__184__ = 4e0 * acegen_scratch__76__;
  acegen_scratch__171__ = 2e0 * acegen_scratch__168__;
  acegen_scratch__50__  = acegen_scratch__168__ * acegen_scratch__170__;
  acegen_scratch__77__  = -(acegen_scratch__28__ * acegen_scratch__31__) +
                         0.7071067811865476e0 * acegen_scratch__50__;
  acegen_scratch__45__ = (acegen_scratch__168__ * acegen_scratch__168__);
  acegen_scratch__75__ =
    acegen_scratch__28__ * acegen_scratch__30__ - acegen_scratch__45__;
  acegen_scratch__183__ = 4e0 * acegen_scratch__75__;
  acegen_scratch__141__ = acegen_scratch__169__ * acegen_scratch__171__ -
                          acegen_scratch__170__ * acegen_scratch__30__;
  acegen_scratch__138__ = acegen_scratch__169__ * acegen_scratch__170__ -
                          acegen_scratch__171__ * acegen_scratch__29__;
  acegen_scratch__133__ =
    -2e0 * acegen_scratch__169__ * acegen_scratch__28__ + acegen_scratch__50__;
  acegen_scratch__79__ = acegen_scratch__169__ * acegen_scratch__32__ -
                         acegen_scratch__30__ * acegen_scratch__33__;
  acegen_scratch__78__ = -(acegen_scratch__29__ * acegen_scratch__32__) +
                         acegen_scratch__169__ * acegen_scratch__33__;
  acegen_scratch__74__ = -(acegen_scratch__169__ * acegen_scratch__169__) +
                         acegen_scratch__29__ * acegen_scratch__30__;
  acegen_scratch__42__ = -(acegen_scratch__29__ * acegen_scratch__45__) -
                         acegen_scratch__30__ * acegen_scratch__47__ +
                         acegen_scratch__169__ * acegen_scratch__50__ +
                         acegen_scratch__28__ * acegen_scratch__74__;
  acegen_scratch__185__ = 0.7071067811865476e0 / acegen_scratch__42__;
  acegen_scratch__172__ = sqrt(acegen_scratch__42__);
  acegen_scratch__173__ =
    (175000e0 / 81e0) / (acegen_scratch__172__ * acegen_scratch__172__);
  acegen_scratch__176__ =
    (acegen_scratch__173__ * acegen_scratch__79__) / acegen_scratch__42__;
  acegen_scratch__178__ =
    (acegen_scratch__173__ * acegen_scratch__78__) / acegen_scratch__42__;
  acegen_scratch__103__ = acegen_scratch__173__ * acegen_scratch__77__;
  acegen_scratch__81__  = 1e0 / (acegen_scratch__42__ * acegen_scratch__42__);
  acegen_scratch__84__  = -(acegen_scratch__77__ * acegen_scratch__81__);
  acegen_scratch__40__ =
    (25000e0 / 81e0) * (-3e0 + 14e0 * log(acegen_scratch__172__));
  acegen_scratch__177__ =
    -(acegen_scratch__40__ * acegen_scratch__78__ * acegen_scratch__81__);
  acegen_scratch__175__ =
    -(acegen_scratch__40__ * acegen_scratch__79__ * acegen_scratch__81__);
  acegen_scratch__174__ =
    (175000e0 / ((acegen_scratch__172__ * acegen_scratch__172__) *
                 acegen_scratch__42__) -
     81e0 * acegen_scratch__40__ * acegen_scratch__81__) /
    162e0;
  acegen_scratch__117__ = (acegen_scratch__175__ + acegen_scratch__176__) / 2e0;
  acegen_scratch__182__ = 4e0 * acegen_scratch__117__;
  acegen_scratch__116__ = (acegen_scratch__177__ + acegen_scratch__178__) / 2e0;
  acegen_scratch__181__ = 4e0 * acegen_scratch__116__;
  acegen_scratch__115__ = (acegen_scratch__103__ / acegen_scratch__42__ +
                           acegen_scratch__40__ * acegen_scratch__84__) /
                          2e0;
  acegen_scratch__114__ = acegen_scratch__174__ * acegen_scratch__76__;
  acegen_scratch__113__ = acegen_scratch__174__ * acegen_scratch__75__;
  acegen_scratch__111__ =
    0.7071067811865476e0 * (acegen_scratch__175__ + acegen_scratch__176__);
  acegen_scratch__110__ =
    0.7071067811865476e0 * (acegen_scratch__177__ + acegen_scratch__178__);
  acegen_scratch__53__  = acegen_scratch__185__ * acegen_scratch__40__;
  acegen_scratch__137__ = -0.1414213562373095e1 * acegen_scratch__53__;
  acegen_scratch__44__  = 0.7071067811865476e0 * acegen_scratch__53__;
  acegen_scratch__180__ = 1e0 * acegen_scratch__44__;
  acegen_scratch__179__ = 2e0 * acegen_scratch__44__;
  acegen_scratch__43__ =
    (25000e0 / 27e0) + acegen_scratch__179__ * acegen_scratch__74__;
  acegen_scratch__46__ =
    (25000e0 / 27e0) + acegen_scratch__179__ * acegen_scratch__75__;
  acegen_scratch__48__ =
    (25000e0 / 27e0) + acegen_scratch__179__ * acegen_scratch__76__;
  acegen_scratch__56__ = acegen_scratch__141__ * acegen_scratch__180__;
  acegen_scratch__57__ = acegen_scratch__138__ * acegen_scratch__180__;
  acegen_scratch__58__ = 1e0 * acegen_scratch__133__ * acegen_scratch__180__;
  valueOut[0]          = 0e0;
  valueOut[1]          = 1000e0;
  valueOut[2]          = 0e0;
  gradientOut[0][0]    = acegen_scratch__13__ * acegen_scratch__43__ +
                      acegen_scratch__14__ * acegen_scratch__56__ +
                      acegen_scratch__15__ * acegen_scratch__57__;
  gradientOut[0][1] = acegen_scratch__14__ * acegen_scratch__46__ +
                      acegen_scratch__13__ * acegen_scratch__56__ +
                      acegen_scratch__15__ * acegen_scratch__58__;
  gradientOut[0][2] = acegen_scratch__15__ * acegen_scratch__48__ +
                      acegen_scratch__13__ * acegen_scratch__57__ +
                      acegen_scratch__14__ * acegen_scratch__58__;
  gradientOut[1][0] = acegen_scratch__16__ * acegen_scratch__43__ +
                      acegen_scratch__17__ * acegen_scratch__56__ +
                      acegen_scratch__18__ * acegen_scratch__57__;
  gradientOut[1][1] = acegen_scratch__17__ * acegen_scratch__46__ +
                      acegen_scratch__16__ * acegen_scratch__56__ +
                      acegen_scratch__18__ * acegen_scratch__58__;
  gradientOut[1][2] = acegen_scratch__18__ * acegen_scratch__48__ +
                      acegen_scratch__16__ * acegen_scratch__57__ +
                      acegen_scratch__17__ * acegen_scratch__58__;
  gradientOut[2][0] = acegen_scratch__19__ * acegen_scratch__43__ +
                      acegen_scratch__20__ * acegen_scratch__56__ +
                      acegen_scratch__21__ * acegen_scratch__57__;
  gradientOut[2][1] = acegen_scratch__20__ * acegen_scratch__46__ +
                      acegen_scratch__19__ * acegen_scratch__56__ +
                      acegen_scratch__21__ * acegen_scratch__58__;
  gradientOut[2][2] = acegen_scratch__21__ * acegen_scratch__48__ +
                      acegen_scratch__19__ * acegen_scratch__57__ +
                      acegen_scratch__20__ * acegen_scratch__58__;
  cacheOut[0] =
    4e0 * acegen_scratch__174__ * (acegen_scratch__74__ * acegen_scratch__74__);
  cacheOut[1]  = 4e0 * (acegen_scratch__30__ * acegen_scratch__44__ +
                       acegen_scratch__113__ * acegen_scratch__74__);
  cacheOut[2]  = 4e0 * (acegen_scratch__29__ * acegen_scratch__44__ +
                       acegen_scratch__114__ * acegen_scratch__74__);
  cacheOut[3]  = 4e0 * (-(acegen_scratch__31__ * acegen_scratch__44__) +
                       acegen_scratch__115__ * acegen_scratch__74__);
  cacheOut[4]  = acegen_scratch__181__ * acegen_scratch__74__;
  cacheOut[5]  = acegen_scratch__182__ * acegen_scratch__74__;
  cacheOut[6]  = acegen_scratch__113__ * acegen_scratch__183__;
  cacheOut[7]  = 4e0 * (acegen_scratch__28__ * acegen_scratch__44__ +
                       acegen_scratch__114__ * acegen_scratch__75__);
  cacheOut[8]  = acegen_scratch__115__ * acegen_scratch__183__;
  cacheOut[9]  = 4e0 * (-(acegen_scratch__32__ * acegen_scratch__44__) +
                       acegen_scratch__116__ * acegen_scratch__75__);
  cacheOut[10] = acegen_scratch__182__ * acegen_scratch__75__;
  cacheOut[11] = acegen_scratch__114__ * acegen_scratch__184__;
  cacheOut[12] = acegen_scratch__115__ * acegen_scratch__184__;
  cacheOut[13] = acegen_scratch__181__ * acegen_scratch__76__;
  cacheOut[14] = 4e0 * (-(acegen_scratch__33__ * acegen_scratch__44__) +
                        acegen_scratch__117__ * acegen_scratch__76__);
  cacheOut[15] =
    2e0 *
    (acegen_scratch__137__ * acegen_scratch__28__ +
     acegen_scratch__133__ *
       (acegen_scratch__103__ * acegen_scratch__185__ +
        acegen_scratch__42__ * acegen_scratch__53__ * acegen_scratch__84__));
  cacheOut[16] = 2e0 * (acegen_scratch__110__ * acegen_scratch__133__ +
                        acegen_scratch__33__ * acegen_scratch__53__);
  cacheOut[17] = 2e0 * (acegen_scratch__111__ * acegen_scratch__133__ +
                        acegen_scratch__32__ * acegen_scratch__53__);
  cacheOut[18] = 2e0 * (acegen_scratch__110__ * acegen_scratch__138__ +
                        acegen_scratch__137__ * acegen_scratch__29__);
  cacheOut[19] = 2e0 * (acegen_scratch__111__ * acegen_scratch__138__ -
                        acegen_scratch__137__ * acegen_scratch__169__);
  cacheOut[20] = 2e0 * (acegen_scratch__111__ * acegen_scratch__141__ +
                        acegen_scratch__137__ * acegen_scratch__30__);
  cacheOut[21] = acegen_scratch__43__;
  cacheOut[22] = acegen_scratch__56__;
  cacheOut[23] = acegen_scratch__57__;
  cacheOut[24] = acegen_scratch__46__;
  cacheOut[25] = acegen_scratch__58__;
  cacheOut[26] = acegen_scratch__48__;
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
  Number acegen_scratch__16__, acegen_scratch__17__, acegen_scratch__18__,
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
    acegen_scratch__70__, acegen_scratch__74__, acegen_scratch__78__,
    acegen_scratch__79__, acegen_scratch__80__, acegen_scratch__81__,
    acegen_scratch__82__, acegen_scratch__83__, acegen_scratch__84__,
    acegen_scratch__88__, acegen_scratch__89__, acegen_scratch__90__;
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
  acegen_scratch__70__ = acegen_scratch__16__ * acegen_scratch__52__ +
                         acegen_scratch__19__ * acegen_scratch__55__ +
                         acegen_scratch__22__ * acegen_scratch__58__;
  acegen_scratch__74__ = acegen_scratch__17__ * acegen_scratch__53__ +
                         acegen_scratch__20__ * acegen_scratch__56__ +
                         acegen_scratch__23__ * acegen_scratch__59__;
  acegen_scratch__78__ = acegen_scratch__18__ * acegen_scratch__54__ +
                         acegen_scratch__21__ * acegen_scratch__57__ +
                         acegen_scratch__24__ * acegen_scratch__60__;
  acegen_scratch__79__ =
    0.7071067811865476e0 * (acegen_scratch__18__ * acegen_scratch__53__ +
                            acegen_scratch__17__ * acegen_scratch__54__ +
                            acegen_scratch__21__ * acegen_scratch__56__ +
                            acegen_scratch__20__ * acegen_scratch__57__ +
                            acegen_scratch__24__ * acegen_scratch__59__ +
                            acegen_scratch__23__ * acegen_scratch__60__);
  acegen_scratch__80__ =
    0.7071067811865476e0 * (acegen_scratch__18__ * acegen_scratch__52__ +
                            acegen_scratch__16__ * acegen_scratch__54__ +
                            acegen_scratch__21__ * acegen_scratch__55__ +
                            acegen_scratch__19__ * acegen_scratch__57__ +
                            acegen_scratch__24__ * acegen_scratch__58__ +
                            acegen_scratch__22__ * acegen_scratch__60__);
  acegen_scratch__81__ =
    0.7071067811865476e0 * (acegen_scratch__17__ * acegen_scratch__52__ +
                            acegen_scratch__16__ * acegen_scratch__53__ +
                            acegen_scratch__20__ * acegen_scratch__55__ +
                            acegen_scratch__19__ * acegen_scratch__56__ +
                            acegen_scratch__23__ * acegen_scratch__58__ +
                            acegen_scratch__22__ * acegen_scratch__59__);
  acegen_scratch__82__ = acegen_scratch__26__ * acegen_scratch__74__ +
                         acegen_scratch__27__ * acegen_scratch__78__ +
                         acegen_scratch__28__ * acegen_scratch__79__ +
                         acegen_scratch__29__ * acegen_scratch__80__ +
                         acegen_scratch__30__ * acegen_scratch__81__ +
                         acegen_scratch__70__ * cacheIn[0];
  acegen_scratch__83__ = acegen_scratch__26__ * acegen_scratch__70__ +
                         acegen_scratch__32__ * acegen_scratch__78__ +
                         acegen_scratch__33__ * acegen_scratch__79__ +
                         acegen_scratch__34__ * acegen_scratch__80__ +
                         acegen_scratch__35__ * acegen_scratch__81__ +
                         acegen_scratch__74__ * cacheIn[6];
  acegen_scratch__84__ = acegen_scratch__27__ * acegen_scratch__70__ +
                         acegen_scratch__32__ * acegen_scratch__74__ +
                         acegen_scratch__37__ * acegen_scratch__79__ +
                         acegen_scratch__38__ * acegen_scratch__80__ +
                         acegen_scratch__39__ * acegen_scratch__81__ +
                         acegen_scratch__78__ * cacheIn[11];
  acegen_scratch__88__ =
    0.7071067811865476e0 * (acegen_scratch__30__ * acegen_scratch__70__ +
                            acegen_scratch__35__ * acegen_scratch__74__ +
                            acegen_scratch__39__ * acegen_scratch__78__ +
                            acegen_scratch__42__ * acegen_scratch__79__ +
                            acegen_scratch__44__ * acegen_scratch__80__ +
                            acegen_scratch__81__ * cacheIn[20]);
  acegen_scratch__89__ =
    0.7071067811865476e0 * (acegen_scratch__29__ * acegen_scratch__70__ +
                            acegen_scratch__34__ * acegen_scratch__74__ +
                            acegen_scratch__38__ * acegen_scratch__78__ +
                            acegen_scratch__41__ * acegen_scratch__79__ +
                            acegen_scratch__44__ * acegen_scratch__81__ +
                            acegen_scratch__80__ * cacheIn[18]);
  acegen_scratch__90__ =
    0.7071067811865476e0 * (acegen_scratch__28__ * acegen_scratch__70__ +
                            acegen_scratch__33__ * acegen_scratch__74__ +
                            acegen_scratch__37__ * acegen_scratch__78__ +
                            acegen_scratch__41__ * acegen_scratch__80__ +
                            acegen_scratch__42__ * acegen_scratch__81__ +
                            acegen_scratch__79__ * cacheIn[15]);
  gradientOut[0][0] = acegen_scratch__16__ * acegen_scratch__46__ +
                      acegen_scratch__17__ * acegen_scratch__47__ +
                      acegen_scratch__18__ * acegen_scratch__48__ +
                      acegen_scratch__52__ * acegen_scratch__82__ +
                      acegen_scratch__53__ * acegen_scratch__88__ +
                      acegen_scratch__54__ * acegen_scratch__89__;
  gradientOut[0][1] = acegen_scratch__16__ * acegen_scratch__47__ +
                      acegen_scratch__17__ * acegen_scratch__49__ +
                      acegen_scratch__18__ * acegen_scratch__50__ +
                      acegen_scratch__53__ * acegen_scratch__83__ +
                      acegen_scratch__52__ * acegen_scratch__88__ +
                      acegen_scratch__54__ * acegen_scratch__90__;
  gradientOut[0][2] = acegen_scratch__16__ * acegen_scratch__48__ +
                      acegen_scratch__17__ * acegen_scratch__50__ +
                      acegen_scratch__18__ * acegen_scratch__51__ +
                      acegen_scratch__54__ * acegen_scratch__84__ +
                      acegen_scratch__52__ * acegen_scratch__89__ +
                      acegen_scratch__53__ * acegen_scratch__90__;
  gradientOut[1][0] = acegen_scratch__19__ * acegen_scratch__46__ +
                      acegen_scratch__20__ * acegen_scratch__47__ +
                      acegen_scratch__21__ * acegen_scratch__48__ +
                      acegen_scratch__55__ * acegen_scratch__82__ +
                      acegen_scratch__56__ * acegen_scratch__88__ +
                      acegen_scratch__57__ * acegen_scratch__89__;
  gradientOut[1][1] = acegen_scratch__19__ * acegen_scratch__47__ +
                      acegen_scratch__20__ * acegen_scratch__49__ +
                      acegen_scratch__21__ * acegen_scratch__50__ +
                      acegen_scratch__56__ * acegen_scratch__83__ +
                      acegen_scratch__55__ * acegen_scratch__88__ +
                      acegen_scratch__57__ * acegen_scratch__90__;
  gradientOut[1][2] = acegen_scratch__19__ * acegen_scratch__48__ +
                      acegen_scratch__20__ * acegen_scratch__50__ +
                      acegen_scratch__21__ * acegen_scratch__51__ +
                      acegen_scratch__57__ * acegen_scratch__84__ +
                      acegen_scratch__55__ * acegen_scratch__89__ +
                      acegen_scratch__56__ * acegen_scratch__90__;
  gradientOut[2][0] = acegen_scratch__22__ * acegen_scratch__46__ +
                      acegen_scratch__23__ * acegen_scratch__47__ +
                      acegen_scratch__24__ * acegen_scratch__48__ +
                      acegen_scratch__58__ * acegen_scratch__82__ +
                      acegen_scratch__59__ * acegen_scratch__88__ +
                      acegen_scratch__60__ * acegen_scratch__89__;
  gradientOut[2][1] = acegen_scratch__22__ * acegen_scratch__47__ +
                      acegen_scratch__23__ * acegen_scratch__49__ +
                      acegen_scratch__24__ * acegen_scratch__50__ +
                      acegen_scratch__59__ * acegen_scratch__83__ +
                      acegen_scratch__58__ * acegen_scratch__88__ +
                      acegen_scratch__60__ * acegen_scratch__90__;
  gradientOut[2][2] = acegen_scratch__22__ * acegen_scratch__48__ +
                      acegen_scratch__23__ * acegen_scratch__50__ +
                      acegen_scratch__24__ * acegen_scratch__51__ +
                      acegen_scratch__60__ * acegen_scratch__84__ +
                      acegen_scratch__58__ * acegen_scratch__89__ +
                      acegen_scratch__59__ * acegen_scratch__90__;
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
{}


template <>
template <typename Number>
inline void
NeoHookeanCached<2>::cache_local(const Tensor<2, dim, Number> &graduIn,
                                 Tensor<1, dim, Number> &      valueOut,
                                 Tensor<2, dim, Number> &      gradientOut,
                                 ArrayView<Number> &           cacheOut)
{}
