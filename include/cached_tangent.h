
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>

using namespace dealii;

namespace Cached
{
  template <int DIMENSION>
  struct Tangent
  {
    const static unsigned int dim      = DIMENSION;
    const static unsigned int n_cached = dim == 3 ? 27 : 9;

    /**
     * Evaluates gradient contribution gradientOut
     *  graduIn gradient of linearization point,
     * gradduIn gradient of current iterate and
     * cacheIn  cached data
     */
    template <typename Number>
    static inline void
    evaluate_reference([[maybe_unused]] const Tensor<2, dim, Number> &graduIn,
                       const Tensor<2, dim, Number> &                 gradduIn,
                       const ArrayView<const Number> &                cacheIn,
                       Tensor<2, dim, Number> &gradientOut);

    /**
     * Evaluates gradient contribution gradientOut in deformend configuration
     * gradduIn gradient of current iterate and
     * cacheIn  cached data
     */
    template <typename Number>
    static inline void
    evaluate_deformed(const Tensor<2, dim, Number> & gradduIn,
                      const ArrayView<const Number> &cacheIn,
                      Tensor<2, dim, Number> &       gradientOut);
  };



  template <>
  template <typename Number>
  inline void
  Tangent<3>::evaluate_reference(
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



  template <>
  template <typename Number>
  inline void
  Tangent<3>::evaluate_deformed(const Tensor<2, dim, Number> & gradduIn,
                                const ArrayView<const Number> &cacheIn,
                                Tensor<2, dim, Number> &       gradientOut)
  {
    Number acegen_scratch__1__, acegen_scratch__11__, acegen_scratch__12__,
      acegen_scratch__13__, acegen_scratch__14__, acegen_scratch__15__,
      acegen_scratch__17__, acegen_scratch__18__, acegen_scratch__19__,
      acegen_scratch__2__, acegen_scratch__20__, acegen_scratch__22__,
      acegen_scratch__23__, acegen_scratch__24__, acegen_scratch__26__,
      acegen_scratch__27__, acegen_scratch__29__, acegen_scratch__3__,
      acegen_scratch__31__, acegen_scratch__32__, acegen_scratch__33__,
      acegen_scratch__34__, acegen_scratch__35__, acegen_scratch__36__,
      acegen_scratch__4__, acegen_scratch__5__, acegen_scratch__55__,
      acegen_scratch__56__, acegen_scratch__57__, acegen_scratch__6__,
      acegen_scratch__64__, acegen_scratch__65__, acegen_scratch__66__,
      acegen_scratch__7__, acegen_scratch__8__, acegen_scratch__9__;
    acegen_scratch__1__  = gradduIn[0][0];
    acegen_scratch__2__  = gradduIn[0][1];
    acegen_scratch__3__  = gradduIn[0][2];
    acegen_scratch__4__  = gradduIn[1][0];
    acegen_scratch__5__  = gradduIn[1][1];
    acegen_scratch__6__  = gradduIn[1][2];
    acegen_scratch__7__  = gradduIn[2][0];
    acegen_scratch__8__  = gradduIn[2][1];
    acegen_scratch__9__  = gradduIn[2][2];
    acegen_scratch__11__ = cacheIn[1];
    acegen_scratch__12__ = cacheIn[2];
    acegen_scratch__13__ = cacheIn[3];
    acegen_scratch__14__ = cacheIn[4];
    acegen_scratch__15__ = cacheIn[5];
    acegen_scratch__17__ = cacheIn[7];
    acegen_scratch__18__ = cacheIn[8];
    acegen_scratch__19__ = cacheIn[9];
    acegen_scratch__20__ = cacheIn[10];
    acegen_scratch__22__ = cacheIn[12];
    acegen_scratch__23__ = cacheIn[13];
    acegen_scratch__24__ = cacheIn[14];
    acegen_scratch__26__ = cacheIn[16];
    acegen_scratch__27__ = cacheIn[17];
    acegen_scratch__29__ = cacheIn[19];
    acegen_scratch__31__ = cacheIn[21];
    acegen_scratch__32__ = cacheIn[22];
    acegen_scratch__33__ = cacheIn[23];
    acegen_scratch__34__ = cacheIn[24];
    acegen_scratch__35__ = cacheIn[25];
    acegen_scratch__36__ = cacheIn[26];
    acegen_scratch__55__ =
      0.7071067811865476e0 * (acegen_scratch__6__ + acegen_scratch__8__);
    acegen_scratch__56__ =
      0.7071067811865476e0 * (acegen_scratch__3__ + acegen_scratch__7__);
    acegen_scratch__57__ =
      0.7071067811865476e0 * (acegen_scratch__2__ + acegen_scratch__4__);
    acegen_scratch__64__ =
      0.7071067811865476e0 * (acegen_scratch__1__ * acegen_scratch__15__ +
                              acegen_scratch__20__ * acegen_scratch__5__ +
                              acegen_scratch__27__ * acegen_scratch__55__ +
                              acegen_scratch__29__ * acegen_scratch__56__ +
                              acegen_scratch__24__ * acegen_scratch__9__ +
                              acegen_scratch__57__ * cacheIn[20]);
    acegen_scratch__65__ =
      0.7071067811865476e0 * (acegen_scratch__1__ * acegen_scratch__14__ +
                              acegen_scratch__19__ * acegen_scratch__5__ +
                              acegen_scratch__26__ * acegen_scratch__55__ +
                              acegen_scratch__29__ * acegen_scratch__57__ +
                              acegen_scratch__23__ * acegen_scratch__9__ +
                              acegen_scratch__56__ * cacheIn[18]);
    acegen_scratch__66__ =
      0.7071067811865476e0 * (acegen_scratch__1__ * acegen_scratch__13__ +
                              acegen_scratch__18__ * acegen_scratch__5__ +
                              acegen_scratch__26__ * acegen_scratch__56__ +
                              acegen_scratch__27__ * acegen_scratch__57__ +
                              acegen_scratch__22__ * acegen_scratch__9__ +
                              acegen_scratch__55__ * cacheIn[15]);
    gradientOut[0][0] =
      acegen_scratch__2__ * acegen_scratch__32__ +
      acegen_scratch__3__ * acegen_scratch__33__ +
      acegen_scratch__11__ * acegen_scratch__5__ +
      acegen_scratch__13__ * acegen_scratch__55__ +
      acegen_scratch__14__ * acegen_scratch__56__ +
      acegen_scratch__15__ * acegen_scratch__57__ +
      acegen_scratch__12__ * acegen_scratch__9__ +
      acegen_scratch__1__ * (acegen_scratch__31__ + cacheIn[0]);
    gradientOut[0][1] = acegen_scratch__1__ * acegen_scratch__32__ +
                        acegen_scratch__2__ * acegen_scratch__34__ +
                        acegen_scratch__3__ * acegen_scratch__35__ +
                        acegen_scratch__64__;
    gradientOut[0][2] = acegen_scratch__1__ * acegen_scratch__33__ +
                        acegen_scratch__2__ * acegen_scratch__35__ +
                        acegen_scratch__3__ * acegen_scratch__36__ +
                        acegen_scratch__65__;
    gradientOut[1][0] = acegen_scratch__31__ * acegen_scratch__4__ +
                        acegen_scratch__32__ * acegen_scratch__5__ +
                        acegen_scratch__33__ * acegen_scratch__6__ +
                        acegen_scratch__64__;
    gradientOut[1][1] =
      acegen_scratch__1__ * acegen_scratch__11__ +
      acegen_scratch__32__ * acegen_scratch__4__ +
      acegen_scratch__18__ * acegen_scratch__55__ +
      acegen_scratch__19__ * acegen_scratch__56__ +
      acegen_scratch__20__ * acegen_scratch__57__ +
      acegen_scratch__35__ * acegen_scratch__6__ +
      acegen_scratch__17__ * acegen_scratch__9__ +
      acegen_scratch__5__ * (acegen_scratch__34__ + cacheIn[6]);
    gradientOut[1][2] = acegen_scratch__33__ * acegen_scratch__4__ +
                        acegen_scratch__35__ * acegen_scratch__5__ +
                        acegen_scratch__36__ * acegen_scratch__6__ +
                        acegen_scratch__66__;
    gradientOut[2][0] = acegen_scratch__65__ +
                        acegen_scratch__31__ * acegen_scratch__7__ +
                        acegen_scratch__32__ * acegen_scratch__8__ +
                        acegen_scratch__33__ * acegen_scratch__9__;
    gradientOut[2][1] = acegen_scratch__66__ +
                        acegen_scratch__32__ * acegen_scratch__7__ +
                        acegen_scratch__34__ * acegen_scratch__8__ +
                        acegen_scratch__35__ * acegen_scratch__9__;
    gradientOut[2][2] =
      acegen_scratch__1__ * acegen_scratch__12__ +
      acegen_scratch__17__ * acegen_scratch__5__ +
      acegen_scratch__22__ * acegen_scratch__55__ +
      acegen_scratch__23__ * acegen_scratch__56__ +
      acegen_scratch__24__ * acegen_scratch__57__ +
      acegen_scratch__33__ * acegen_scratch__7__ +
      acegen_scratch__35__ * acegen_scratch__8__ +
      acegen_scratch__9__ * (acegen_scratch__36__ + cacheIn[11]);
  }

  // =====
  // dim=2
  // =====


  template <>
  template <typename Number>
  inline void
  Tangent<2>::evaluate_reference(
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
  Tangent<2>::evaluate_deformed(const Tensor<2, dim, Number> & gradduIn,
                                const ArrayView<const Number> &cacheIn,
                                Tensor<2, dim, Number> &       gradientOut)
  {
    Number acegen_scratch__1__, acegen_scratch__11__, acegen_scratch__12__,
      acegen_scratch__13__, acegen_scratch__2__, acegen_scratch__22__,
      acegen_scratch__26__, acegen_scratch__3__, acegen_scratch__4__,
      acegen_scratch__6__, acegen_scratch__7__, acegen_scratch__9__;
    acegen_scratch__1__  = gradduIn[0][0];
    acegen_scratch__2__  = gradduIn[0][1];
    acegen_scratch__3__  = gradduIn[1][0];
    acegen_scratch__4__  = gradduIn[1][1];
    acegen_scratch__6__  = cacheIn[1];
    acegen_scratch__7__  = cacheIn[2];
    acegen_scratch__9__  = cacheIn[4];
    acegen_scratch__11__ = cacheIn[6];
    acegen_scratch__12__ = cacheIn[7];
    acegen_scratch__13__ = cacheIn[8];
    acegen_scratch__22__ =
      0.7071067811865476e0 * (acegen_scratch__2__ + acegen_scratch__3__);
    acegen_scratch__26__ =
      0.7071067811865476e0 * (acegen_scratch__1__ * acegen_scratch__7__ +
                              acegen_scratch__4__ * acegen_scratch__9__ +
                              acegen_scratch__22__ * cacheIn[5]);
    gradientOut[0][0] =
      acegen_scratch__12__ * acegen_scratch__2__ +
      acegen_scratch__4__ * acegen_scratch__6__ +
      acegen_scratch__22__ * acegen_scratch__7__ +
      acegen_scratch__1__ * (acegen_scratch__11__ + cacheIn[0]);
    gradientOut[0][1] = acegen_scratch__1__ * acegen_scratch__12__ +
                        acegen_scratch__13__ * acegen_scratch__2__ +
                        acegen_scratch__26__;
    gradientOut[1][0] = acegen_scratch__26__ +
                        acegen_scratch__11__ * acegen_scratch__3__ +
                        acegen_scratch__12__ * acegen_scratch__4__;
    gradientOut[1][1] =
      acegen_scratch__12__ * acegen_scratch__3__ +
      acegen_scratch__1__ * acegen_scratch__6__ +
      acegen_scratch__22__ * acegen_scratch__9__ +
      acegen_scratch__4__ * (acegen_scratch__13__ + cacheIn[3]);
  }


} // namespace Cached