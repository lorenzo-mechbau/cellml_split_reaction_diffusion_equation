PROGRAM CELLML_SPLIT_REACTION_DIFFUSION_EQUATION

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif








  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM CELLML_SPLIT_REACTION_DIFFUSION_EQUATION








