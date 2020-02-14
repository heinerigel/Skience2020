!--------------------------------------------------------------------------
program sem1D
!--------------------------------------------------------------------------
  use module_sem1d
!
  type(model), target  :: mod
  type(simu) :: sim
  logical :: h=.true.
!

  call load_sem_dat(sim)
  if (sim%horder<0) then
     call read_model(modname,mod) 
  else
     call read_model(modname,mod,h) 
  endif
  call mesh_model(mod,sim) 
  if (sim%horder>0) call init_correctors(sim)
  if (sim%scheme==ACC) then
     call newmark(sim)
  else
     call velocity_stress(sim)
  endif
!--------------------------------------------------------------------------
end program sem1D
!--------------------------------------------------------------------------
