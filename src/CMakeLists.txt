add_executable(adaptive adaptive.cc driver.hh
               estimator.hh
               operator.hh
               bctype.hh
               coefficients.hh
               norm.hh
               adaptive.ini)


dune_symlink_to_source_files(FILES adaptive.ini ldomain.msh ldomain.geo)
