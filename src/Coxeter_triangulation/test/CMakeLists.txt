project(Coxeter_triangulation_test)

include(GUDHI_test_coverage)


add_executable ( Coxeter_triangulation_permutahedral_representation_test perm_rep_test.cpp )
target_link_libraries( Coxeter_triangulation_permutahedral_representation_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
gudhi_add_coverage_test(Coxeter_triangulation_permutahedral_representation_test)

add_executable ( Coxeter_triangulation_freudenthal_triangulation_test freud_triang_test.cpp )
target_link_libraries( Coxeter_triangulation_freudenthal_triangulation_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
gudhi_add_coverage_test(Coxeter_triangulation_freudenthal_triangulation_test)

add_executable ( Coxeter_triangulation_functions_test function_test.cpp )
target_link_libraries( Coxeter_triangulation_functions_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
gudhi_add_coverage_test(Coxeter_triangulation_functions_test)

add_executable ( Coxeter_triangulation_oracle_test oracle_test.cpp )
target_link_libraries( Coxeter_triangulation_oracle_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
gudhi_add_coverage_test(Coxeter_triangulation_oracle_test)

add_executable ( Coxeter_triangulation_manifold_tracing_test manifold_tracing_test.cpp )
target_link_libraries( Coxeter_triangulation_manifold_tracing_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
gudhi_add_coverage_test(Coxeter_triangulation_manifold_tracing_test)

add_executable ( Coxeter_triangulation_cell_complex_test cell_complex_test.cpp )
target_link_libraries( Coxeter_triangulation_cell_complex_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
gudhi_add_coverage_test(Coxeter_triangulation_cell_complex_test)

if (TBB_FOUND)
   target_link_libraries(Coxeter_triangulation_cell_complex_test ${TBB_LIBRARIES})
endif()

