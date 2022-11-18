gfortran nrutil.f90 nrtype.f90 broydn.f90 wal1.f90 fdjac.f90  qrdcmp.f90 fmin.f90 lnsrch.f90 rsolv.f90 qrupdt.f90 rotate.f90 pythag.f90 -fdebug -fmax-errors=1 -fbounds-check -fsanitize=address -O0 -o ex

#-fmax-errors=1  ===> para no primeiro erro
#-w ====> ignora os warnings
# para identificar erros: -fdebug
#https://gcc.gnu.org/onlinedocs/gfortran/gnu-fortran-command-options/options-for-debugging-your-program-or-gnu-fortran.html

