!###########################################################################################################################################################

module ModuleRunFerdosi
    
    implicit none

    !###########################################################################################################################################################

    contains
        
        !###########################################################################################################################################################


        !---------------------------------------------------------------------------
        !> @brief Function to run Ferdosi from externall
        
        !> @author  David Wilson david.wilson@roslin.ed.ac.uk
        !> @date    October 26, 2016
        !---------------------------------------------------------------------------
        subroutine doFerdosi(AllParameters,ped, sireDamOpt)
            use ModuleParameters

            type(Parameters) :: AllParameters
            type(PedigreeHolder) :: ped
            integer, optional, intent(in) ::sireDamOpt


            ! TODO add in support for iterating through dams :P 

            
            ! This trims the array to only include animals with 5 or more offspring, and that they are genotyped
            call ped%sireList%removeIndividualsBasedOnThreshold(nOffsThresh=5, genotyped=.true.)
            allocate(AllParameters%SireArray(ped%sireList%length))
            AllParameters%nSire = ped%sireList%length
            block
                type(IndividualLinkedListNode),pointer :: tmpInd
                integer :: i
                tmpInd => ped%sireList%first
                do i=1, ped%sireList%length
                    call tmpInd%item%initPhaseArrays(AllParameters%nsnp)
                    AllParameters%SireArray(i)%ind => tmpInd%item
                    tmpInd => tmpInd%next
                enddo

            end block 
            call RunFerdosiPerSire(AllParameters, ped)



        end subroutine doFerdosi


        subroutine RunFerdosiPerSire(AllParameters, ped)

            use ModuleParameters
            use Ferdosi
            use PedigreeModule
            implicit none

            type(Parameters), intent(inout) :: AllParameters
            type(PedigreeHolder) :: ped
            type(IndividualLinkedListNode),pointer :: ind

            integer :: i
            
            print*, "Running Ferdosi for each sire..."
            ind => ped%sireList%first
            do i=1, ped%sireList%length
                call FerdosiRunner(ind%item, AllParameters%nSnp, AllParameters%overwriteHalfSibPhase)
                ind => ind%next
            enddo

            
        end subroutine RunFerdosiPerSire

        !###########################################################################################################################################################

 
 
end module ModuleRunFerdosi

!###########################################################################################################################################################
