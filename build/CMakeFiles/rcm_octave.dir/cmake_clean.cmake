file(REMOVE_RECURSE
  "rcm_octave.o"
  "librcm_octave.pdb"
  "librcm_octave.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/rcm_octave.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
