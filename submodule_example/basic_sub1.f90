submodule (my_math) double_pi
contains
module procedure pi_mult
  print *,"Running in double_pi"
  tau = 2*pi
end procedure pi_mult
end submodule double_pi
