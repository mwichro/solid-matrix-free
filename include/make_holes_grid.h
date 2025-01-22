
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>


namespace MakeGrid
{
  using namespace dealii;

  void
  merge(parallel::distributed::Triangulation<3> &    dst,
        const std::vector<const Triangulation<2> *> &triangulations,
        const double                                 scale,
        const Point<3> &                             center_dim_1,
        const Point<3> &                             center_dim_2,
        const Point<3> &                             center_dim_3,
        const double                                 extrusion_height,
        const unsigned int                           extrusion_slices)
  {
    Triangulation<2> triangulation_2d;
    GridGenerator::merge_triangulations(triangulations,
                                        triangulation_2d,
                                        0.01 * scale,
                                        true);
    for (auto &id : triangulation_2d.get_manifold_ids())
      if (id != numbers::flat_manifold_id)
        triangulation_2d.set_manifold(id, FlatManifold<2>());
    GridGenerator::extrude_triangulation(
      triangulation_2d, extrusion_slices, extrusion_height * scale, dst, true);

    for (const auto &cell : dst.active_cell_iterators())
      for (unsigned int face_no = 0; face_no < GeometryInfo<3>::faces_per_cell;
           ++face_no)
        if (!cell->at_boundary(face_no))
          if (cell->manifold_id() != cell->neighbor(face_no)->manifold_id())
            cell->face(face_no)->set_all_manifold_ids(
              cell->face(face_no)->manifold_id());

    // we need to set manifolds:
    Tensor<1, 3> dir;
    dir[0] = 0.;
    dir[1] = 0.;
    dir[2] = 1.;
    CylindricalManifold<3> cylindrical_manifold_1(dir, center_dim_1);
    CylindricalManifold<3> cylindrical_manifold_2(dir, center_dim_2);
    CylindricalManifold<3> cylindrical_manifold_3(dir, center_dim_3);

    dst.set_manifold(1, cylindrical_manifold_1);
    dst.set_manifold(2, cylindrical_manifold_2);
    dst.set_manifold(3, cylindrical_manifold_3);
  }

  void
  merge(parallel::distributed::Triangulation<2> &    dst,
        const std::vector<const Triangulation<2> *> &triangulations,
        const double                                 scale,
        const Point<2> &                             center_dim_1,
        const Point<2> &                             center_dim_2,
        const Point<2> &                             center_dim_3,
        const double /*extrusion_height*/,
        const unsigned int /*extrusion_slices*/)

  {
    GridGenerator::merge_triangulations(triangulations,
                                        dst,
                                        0.01 * scale,
                                        true);

    // we need to set manifolds:
    PolarManifold<2> polar_manifold_1(center_dim_1);
    PolarManifold<2> polar_manifold_2(center_dim_2);
    PolarManifold<2> polar_manifold_3(center_dim_3);

    dst.set_manifold(1, polar_manifold_1);
    dst.set_manifold(2, polar_manifold_2);
    dst.set_manifold(3, polar_manifold_3);
  }


  template <int dim>
  void
  make_holes_grid(parallel::distributed::Triangulation<dim> &triangulation,
                  double                                     parameters_scale,
                  double                                     extrusion_height,
                  double                                     extrusion_slices)
  {
    // plate with a hole and 2 inclusions (geometry from Miehe 2007,
    // On multiscale FE analyses...)
    Point<2> center_1, center_2, center_3;
    center_1[0]      = -0.2 * parameters_scale;
    center_1[1]      = -0.2 * parameters_scale;
    center_2[0]      = -0.2 * parameters_scale;
    center_2[1]      = 0.2 * parameters_scale;
    center_3[0]      = 0.2 * parameters_scale;
    center_3[1]      = 0.0 * parameters_scale;
    const double R   = 0.15 * parameters_scale;
    const double R2  = 0.2 * parameters_scale;
    const double pLR = 0.1 * parameters_scale;
    const double pBT = 0.2 * parameters_scale;

    // inclusion:
    Triangulation<2> sphere_2, sphere_3;

    auto create_inclusion = [&](Triangulation<2> &       out,
                                const Point<2> &         center,
                                const double             radius,
                                const types::manifold_id tfi_manifold_id,
                                const types::manifold_id ball_id) -> void {
      Triangulation<2> sphere;
      GridGenerator::hyper_ball(sphere, center, radius);

      for (const auto &cell : sphere.active_cell_iterators())
        {
          if (cell->center().distance(center) < 1e-8 * parameters_scale)
            {
              cell->set_all_manifold_ids(numbers::flat_manifold_id);
            }
          else
            {
              cell->set_manifold_id(tfi_manifold_id);
            }
          cell->set_material_id(2);
        }

      sphere.set_manifold(tfi_manifold_id, FlatManifold<2>());
      sphere.refine_global(1);
      // at this point we have 8 faces across circumference
      GridGenerator::flatten_triangulation(sphere, out);
      out.set_all_manifold_ids_on_boundary(ball_id);
    };

    create_inclusion(sphere_2, center_2, R, 7, 2);
    create_inclusion(sphere_3, center_3, R, 8, 3);

    Triangulation<2> plate_1;
    GridGenerator::plate_with_a_hole(plate_1,
                                     R /*inner_radius*/,
                                     R2 /*outer_radius*/,
                                     0. /*pad_bottom*/,
                                     0. /*pad_top*/,
                                     pLR /*pad_left*/,
                                     0. /*pad_right*/,
                                     center_1 /*center*/,
                                     1 /*polar_manifold_id*/,
                                     4 /*tfi_manifold_id*/,
                                     1. /*L*/,
                                     1. /*n_slices*/,
                                     false /*colorize*/);
    for (const auto &cell : plate_1.active_cell_iterators())
      for (unsigned int face_no = 0; face_no < GeometryInfo<2>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no) &&
            cell->face(face_no)->manifold_id() != 1)
          cell->face(face_no)->set_all_manifold_ids(numbers::flat_manifold_id);

    Triangulation<2> plate_2;
    GridGenerator::plate_with_a_hole(plate_2,
                                     R /*inner_radius*/,
                                     R2 /*outer_radius*/,
                                     0. /*pad_bottom*/,
                                     0. /*pad_top*/,
                                     pLR /*pad_left*/,
                                     0. /*pad_right*/,
                                     center_2 /*center*/,
                                     2 /*polar_manifold_id*/,
                                     5 /*tfi_manifold_id*/,
                                     1. /*L*/,
                                     1. /*n_slices*/,
                                     false /*colorize*/);
    for (const auto &cell : plate_2.active_cell_iterators())
      for (unsigned int face_no = 0; face_no < GeometryInfo<2>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no) &&
            cell->face(face_no)->manifold_id() != 2)
          cell->face(face_no)->set_all_manifold_ids(numbers::flat_manifold_id);

    Triangulation<2> plate_3;
    GridGenerator::plate_with_a_hole(plate_3,
                                     R /*inner_radius*/,
                                     R2 /*outer_radius*/,
                                     pBT /*pad_bottom*/,
                                     pBT /*pad_top*/,
                                     0. /*pad_left*/,
                                     pLR /*pad_right*/,
                                     center_3 /*center*/,
                                     3 /*polar_manifold_id*/,
                                     6 /*tfi_manifold_id*/,
                                     1. /*L*/,
                                     1. /*n_slices*/,
                                     false /*colorize*/);
    for (const auto &cell : plate_3.active_cell_iterators())
      for (unsigned int face_no = 0; face_no < GeometryInfo<2>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no) &&
            cell->face(face_no)->manifold_id() != 3)
          cell->face(face_no)->set_all_manifold_ids(numbers::flat_manifold_id);

    Triangulation<2>                       top, bottom;
    const std::vector<std::vector<double>> step_sizes = {
      {0.1 * parameters_scale,
       0.2 * parameters_scale,
       0.2 * parameters_scale,
       0.2 * parameters_scale,
       0.2 * parameters_scale,
       0.1 * parameters_scale},
      {0.1 * parameters_scale}};
    Point<2> bl, tr;
    bl[0] = -0.5 * parameters_scale;
    bl[1] = 0.4 * parameters_scale;
    tr[0] = 0.5 * parameters_scale;
    tr[1] = 0.5 * parameters_scale;
    GridGenerator::subdivided_hyper_rectangle(top, step_sizes, bl, tr);

    bl[1] = -0.5 * parameters_scale;
    tr[1] = -0.4 * parameters_scale;
    GridGenerator::subdivided_hyper_rectangle(bottom, step_sizes, bl, tr);

    Point<dim> center_dim_1, center_dim_2, center_dim_3;
    for (unsigned int d = 0; d < 2; ++d)
      {
        center_dim_1[d] = center_1[d];
        center_dim_2[d] = center_2[d];
        center_dim_3[d] = center_3[d];
      }

    merge(triangulation,
          {&plate_1, &plate_2, &plate_3, &sphere_2, &sphere_3, &top, &bottom},
          parameters_scale,
          center_dim_1,
          center_dim_2,
          center_dim_3,
          extrusion_height,
          extrusion_slices);


    for (unsigned int i = 4; i <= 8; ++i)
      triangulation.set_manifold(i, FlatManifold<dim>());
    for (unsigned int i = 4; i <= 8; ++i)
      {
        TransfiniteInterpolationManifold<dim> transfinite_manifold;
        transfinite_manifold.initialize(triangulation);
        triangulation.set_manifold(i, transfinite_manifold);
      }

    const double tol_boundary = 1e-6 * parameters_scale;
    for (auto cell : triangulation.active_cell_iterators())
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary() == true)
          {
            if (std::abs(cell->face(face)->center()[1] -
                         (-0.5 * parameters_scale)) < tol_boundary)
              cell->face(face)->set_boundary_id(1); // -Y faces
            else if (std::abs(cell->face(face)->center()[1] -
                              0.5 * parameters_scale) < tol_boundary)
              cell->face(face)->set_boundary_id(11); // +Y faces
          }
  }

} // namespace MakeGrid