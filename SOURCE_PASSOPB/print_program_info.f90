!--------------------------------------------------------------------
! This subroutine print out the program title.
!
! Developer : M.K. Lee
! Last updated : 2013.01.03
!
! Bug report : mk@ices.utexas.edu
!--------------------------------------------------------------------
subroutine print_program_info
  implicit none
  call system('clear')
  write(*,'(a120)')" "
  write(*,'(a120)')" "
  write(*,'(a73)')" ######  ####### ####### #     #  #####  ######     #     #####  #    # "
  write(*,'(a73)')" #     # #     # #     # ##    # #     # #     #   # #   #     # #   #  "
  write(*,'(a73)')" #     # #     # #     # # #   # #       #     #  #   #  #       #  #   "
  write(*,'(a73)')" ######  #     # #     # #  #  # #  #### ######  #     # #       ###    "
  write(*,'(a73)')" #       #     # #     # #   # # #     # #     # ####### #       #  #   "
  write(*,'(a73)')" #       #     # #     # #    ## #     # #     # #     # #     # #   #  "
  write(*,'(a73)')" #       ####### ####### #     #  #####  ######  #     #  #####  #    # "

  write(*,'(a73)')""
  write(*,'(a73)')""
  write(*,'(a73)')"  ####  #     #    #    #     # #     # ####### #       "
  write(*,'(a73)')" #    # #     #   # #   ##    # ##    # #       #       "
  write(*,'(a73)')" #      #     #  #   #  # #   # # #   # #       #       "
  write(*,'(a73)')" #      ####### #     # #  #  # #  #  # #####   #       "
  write(*,'(a73)')" #      #     # ####### #   # # #   # # #       #       "
  write(*,'(a73)')" #    # #     # #     # #    ## #    ## #       #       "
  write(*,'(a73)')"  ####  #     # #     # #     # #     # ####### ####### "
  write(*,'(a73)')" "
  write(*,'(a73)')"                                                             Version 0.0"
  write(*,'(a73)')"                                          Bug report: mk@ices.utexas.edu"
  write(*,'(a73)')" "

  return
end subroutine print_program_info
