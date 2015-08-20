cet_test(test_Range SOURCES test_Range.cxx
  LIBRARIES LbTupler
)

cet_test(test_geometry SOURCES test_geometry.cxx
  LIBRARIES lbne_Geometry LbTupler
)

cet_test(draw_geometry SOURCES draw_geometry.cxx
  LIBRARIES LbTupler
)

install_headers()
install_source()

