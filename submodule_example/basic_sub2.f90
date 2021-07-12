submodule (my_math) quad_pi
contains
module procedure pi_mult
  print *,"Running in quad_pi"
  tau = 4*pi
end procedure pi_mult
end submodule quad_pi
