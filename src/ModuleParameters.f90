!###########################################################################################################################################################

module ModuleParameters

    use ModuleSireType

	implicit none

	public :: ReadSpec

    !###########################################################################################################################################################

	type Parameters

		! from specfile
		integer :: nSire, nSnp, nDam
		character(len=:), allocatable :: PedigreeFile, GenotypesFileOffspring, PhaseFileSire
        logical :: overwriteHalfSibPhase !< Phase half sibs
        
        ! from program
        integer :: nOffspringTotal
        type(Sire), allocatable, dimension(:) :: SireArray
        integer, allocatable, dimension(:,:) :: Offspring2SireArray
        


	end type Parameters

    !###########################################################################################################################################################

    contains 

        !###########################################################################################################################################################

        function initSpecDefaults(ped) result(Parameters)
            use PedigreeModule
            
            type(Parameters) :: Parameters


            params%nSire = ped%sireList%length
            params%nDam = ped%damList%length
            params%overwriteHalfSibPhase = .false.

            if (ped%nGenotyped > 0) then
                params%nsnp = ped%pedigree(ped%genotypeMap(1))%individualGenotype%length

            else 
                write(error_unit,*) "ERROR - Ferdosi given a pedigree with no genotyped animals"
                stop
            endif
            
            if 


        end function initSpecDefaults

        subroutine ReadSpec(AllParameters, specfile)

            use AlphaHouseMod, only : ToLower
            implicit none

            character(len=*), intent(in) :: specfile
            class(Parameters), intent(inout) :: AllParameters

            integer :: FileLength, stat, i, LineRead, sizeline, fileunit
            character(len=1) :: comma
            character(len=30) :: SpecParam, buffer
            character(len=:), allocatable :: TempStr

            print *, "Reading Specfile"
            open(newunit=fileunit, file=specfile, status="old")


            AllParameters%overwriteHalfSibPhase = .false.
            FileLength = 0
            buffer=""
            TempStr=""

            do
                read(fileunit, *, iostat=stat) SpecParam
                if (stat/=0) exit
                FileLength = FileLength + 1
            enddo  

            rewind(fileunit)

            do i=1, FileLength
                read(fileunit,'(a30,A)', advance='NO', iostat=stat) SpecParam
                read(fileunit,'(a1)', advance='NO', iostat=stat) comma

                select case(trim(ToLower(SpecParam)))
         
                    case('pedigreefile')
                        stat=0
                        LineRead=0
                        do while (stat==0)
                            read (fileunit, "(A)", advance='NO', size=sizeline, iostat=stat) buffer
                            LineRead=LineRead+1
                            AllParameters%PedigreeFile = AllParameters%PedigreeFile // buffer
                        enddo
                        if (LineRead == 0) then
                            print *, "PedigreeFile not set properly in spec file"
                            stop 1
                        endif

                    case('genotypesfileoffspring')
                        stat=0
                        LineRead=0
                        do while (stat==0)
                            read (fileunit, "(A)", advance='NO', size=sizeline, iostat=stat) buffer
                            LineRead=LineRead+1
                            AllParameters%GenotypesFileOffspring = AllParameters%GenotypesFileOffspring // buffer
                        enddo
                        if (LineRead == 0) then
                            print *, "GenotypesFileOffspring not set properly in spec file"
                            stop 2
                        endif

                    case('phasefilesire')
                        stat=0
                        LineRead=0
                        do while (stat==0)
                            read (fileunit, "(A)", advance='NO', size=sizeline, iostat=stat) buffer
                            LineRead=LineRead+1
                            AllParameters%phaseFileSire = AllParameters%phaseFileSire // buffer
                        enddo
                        if (LineRead == 0) then
                            print *, "PhaseFileSire not set properly in spec file"
                            stop 3
                        endif

                    case('numberofsires')
                        read(fileunit, *, iostat=stat) AllParameters%nSire
                        if (stat /= 0) then
                            print *, "NumberOfSires not set properly in spec file"
                            stop 4
                        endif   

                    case('numberofsnp')
                        read(fileunit, *, iostat=stat) AllParameters%nSnp
                        if (stat /= 0) then
                            print *, "NumberOfSnp not set properly in spec file"
                            stop 5
                        endif   

                    case default
                        print *, "Error in specfile, please check", SpecParam
                        stop 100
                end select
            
            enddo

            close(fileunit)


        end subroutine ReadSpec

        !###########################################################################################################################################################


end module ModuleParameters

!###########################################################################################################################################################
