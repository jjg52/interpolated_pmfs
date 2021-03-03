program potential_mean_force
! include one or more water molecules in QM region
        use dynamo
        implicit none
        logical, dimension(:), allocatable      :: acs
        character( len=6 )                      :: si
        integer                                 :: a1, a2, a3, a4, i, my_random
        real(dp)                                :: eq
        call dynamo_header
        read(*,*) i
        read(*,*) eq
        call encode_integer ( i, si, '(i3)' )
        call mm_file_process ( 'borra', 'BAL.opls' )
        call mm_system_construct ( 'borra', trim(si) // '.seq' )
! edit the next line to specify the initial coordinates
        call coordinates_read ( trim(si) // ".crd" )
        allocate( acs(1:natoms) )
        acs = .false.
        acs = atom_selection( &
                                subsystem      = (/ "SOLUTE" /), &
                                residue_number = (/ 1 /), &
                                residue_name   = (/ "BAL" /) )
! change 999 to the number of the water molecule to include
!        acs = acs .or. atom_selection( &
!                                subsystem      = (/ "WAT" /), &
!                                residue_number = (/ 700 /), &
!                                residue_name   = (/ "HOH" /) )
! add more lines as necessary to include more waters in QM region
!        acs = acs .or. atom_selection( &
!                                subsystem      = (/ "WAT" /), &
!                                residue_number = (/ 195 /), &
!                                residue_name   = (/ "HOH" /) )
! add more lines as necessary to include more waters in QM region
!        acs = acs .or. atom_selection( &
!                                subsystem      = (/ "WAT" /), &
!                                residue_number = (/ 112 /), &
!                                residue_name   = (/ "HOH" /) )
! add more lines as necessary to include more waters in QM region
!        acs = acs .or. atom_selection( &
!                                subsystem      = (/ "WAT" /), &
!                                residue_number = (/ 402 /), &
!                                residue_name   = (/ "HOH" /) )
! add more lines as necessary to include more waters in QM region
!        acs = acs .or. atom_selection( &
!                                subsystem      = (/ "WAT" /), &
!                                residue_number = (/ 455 /), &
!                                residue_name   = (/ "HOH" /) )
        call mopac_setup ( &
                method    = "AM1", &
                charge    = 0, &
                selection = acs )
        call energy_initialize
        call energy_non_bonding_options ( &
                list_cutoff   = 15.5_dp, &
                outer_cutoff  = 13.0_dp, &
                inner_cutoff  = 12.5_dp, &
                minimum_image = .true. )
        a1 = atom_number ( &
                atom_name      = "C4", &
                residue_number = 1, &
                subsystem      = "SOLUTE" )
        a2 = atom_number ( &
                atom_name      = "BR6", &
                residue_number = 1, &
                subsystem      = "SOLUTE" )
        a3 = atom_number ( &
                atom_name      = "C4", &
                residue_number = 1, &
                subsystem      = "SOLUTE" )
        a4 = atom_number ( &
                atom_name      = "O9", &
                residue_number = 1, &
                subsystem      = "SOLUTE" )

                
                call constraint_initialize
                call constraint_point_define ( atom_selection ( atom_number = (/ a1 /) ) )
                call constraint_point_define ( atom_selection ( atom_number = (/ a2 /) ) )
                call constraint_point_define ( atom_selection ( atom_number = (/ a3 /) ) )
                call constraint_point_define ( atom_selection ( atom_number = (/ a4 /) ) )                
                call constraint_define ( &
                        type        = "MULTIPLE_DISTANCE", &
                        fc          = 2500.0_dp, &
                        eq          = eq, &
                        weights     = (/ 1.0_dp, -1.0_dp /), &
                        file        = "dat." // trim ( si ) )
                call gradient
! Relax MD
                call random_initialize ( my_random() + i )
                call velocity_assign ( 300._dp, .false. )
                call dynamics_options ( &
                        time_step       = .001_dp, &
                        print_frequency = 100, &
                        steps           = 25000 )
                call langevin_verlet_dynamics ( 300._dp, 100._dp )
! Production MD
                call dynamics_options ( &
                        time_step       = .001_dp, &
                        print_frequency = 500, &
                        steps           = 25000 )
                call constraint_writing_start
                call langevin_verlet_dynamics ( 300._dp, 100._dp )
                call constraint_writing_stop
                call coordinates_write ( "m_" // trim( si ) // ".crd" )
                call velocity_write ( "m_" // trim( si ) // ".vel" )
                call pdb_write( 'm_' // trim( si ) // '.pdb' )
        deallocate( acs )
                call dynamo_footer
end program
