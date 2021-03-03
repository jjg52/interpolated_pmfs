program const
	use dynamo
	implicit none

	real( kind=dp )                     :: dst
        character(12)                       :: sl
	integer                             :: i, ires, nsatoms, nsresid, o
	logical, allocatable, dimension(:)  :: flg 
	real( kind=dp ), dimension(1:3)     :: ref
        character(len = 20)                 :: opls, seq, crd

	call dynamo_header
        read(*,*) sl
!
! file "opls" must be in the same directory as this program
! and must contain parameters for all the atoms of the solute
!  
	call mm_file_process( "borra", 'BAL.opls' )
!
! file "THPA.seq" contains sequence information for the solute
!
	call mm_system_construct( "borra", 'BAL.seq' )
	call mm_system_write( "borra.prt" )
	nsatoms = natoms
	nsresid = nresid
!
! file "wat.seq" contains sequence information for the solvent
!
	call mm_system_construct( "borra", "wat.seq" )
	call mm_system_write( "borra.wat" )

	call mm_system_read( "borra.prt" )
	call mm_system_append( "borra.wat" )

	allocate( flg(1:natoms) )
	           flg = .false.
	flg(1:nsatoms) = .true.
!
! files "THPA.crd" and "wat.crd" contain cartesian coordinates for the
! solute and solvent, respectively
!
	call coordinates_read( 'cart_zma_'//trim(sl)//'.crd', selection = flg )
	call coordinates_read( "wat.crd", selection = .not. flg )

!	ref(1:3) = atmcrd(1:3,atom_number( subsystem = "X", residue_number = 3, atom_name = "C6" ) )
!	do i = 1, nsatoms
!		atmcrd(1:3,i) = atmcrd(1:3,i) - ref(1:3)
!	end do
!
! move the origin to the centre of mass
!
	call translate_to_center( atmcrd(1:3,nsatoms+1:natoms), atmmas(nsatoms+1:natoms) )
!
! remove solvent molecules with any atom less than or equal to 2.8 Å from any solute atom
!
	dst = 2.8_dp ** 2
	flg = .false.
	do ires = nsresid+1, nresid
		o = resind(ires)+1
		do i = 1, nsatoms
			if( atmnum(i) /= 1 ) then
				if( sum( ( atmcrd(1:3,i) - atmcrd(1:3,o) ) ** 2 ) <= dst ) then
					flg(resind(ires)+1:resind(ires+1)) = .true.
					exit
				end if
			end if
		end do
	end do
	write(*,*) ' >> DELETE', count( flg )
	call mm_system_delete( flg )
!
! write out the sequence and coordinates of the solute in solvent box
!
        call sequence_write( trim(sl)//".seq" )
	call coordinates_write( trim(sl)//".crd" )
        call pdb_write( trim(sl)//".pdb" )

	deallocate( flg )
end
