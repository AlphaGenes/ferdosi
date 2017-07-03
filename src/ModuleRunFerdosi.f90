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
        subroutine doFerdosi(ped,AllParametersIn, sireDamOpt)
            use ModuleParameters


            class(PedigreeHolder), intent(inout) :: ped
            type(Parameters), optional :: AllParametersIn
            integer, optional, intent(in) ::sireDamOpt

            type(Parameters) :: AllParameters


            if (present(AllParametersIn)) then
                AllParameters = AllParametersIn
            else 
                AllParameters = initSpecDefaults(ped)
            endif
            ! This trims the array to only include animals with 5 or more offspring, and that they are genotyped
            call ped%sireList%removeIndividualsBasedOnThreshold(nOffsThresh=5)
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

            if (present(sireDamOpt))  then

                if (sireDamOpt == 2) then
                    call RunFerdosiPerDam(AllParameters,ped)
                else !< default, run per sire
                    call RunFerdosiPerSire(AllParameters, ped)
                endif
            else !< default, run per sire
                call RunFerdosiPerSire(AllParameters, ped)
            endif
            call convertPhaseInfo(ped)


        end subroutine doFerdosi


        subroutine convertPhaseInfo(ped)
            use HaplotypeModule
            use PedigreeModule

            type(PedigreeHolder), intent(inout) :: ped
            type(haplotype) :: hap1, hap2
            integer :: i
            do i=1, ped%pedigreeSize
                hap1 = ped%pedigree(i)%individualPhase(1)
                hap2 = ped%pedigree(i)%individualPhase(2)
                call ped%pedigree(i)%individualGenotype%setFromHaplotypesIfMissing(hap1,hap2)
                ! TODO check if the following is needed
                ! ped%pedigree(i)%genotyped = .true.
            enddo 



        end subroutine convertPhaseInfo


        subroutine RunFerdosiPerSire(AllParameters, ped)

            use ModuleParameters
            use Ferdosi
            use PedigreeModule
            implicit none

            type(Parameters), intent(inout) :: AllParameters
            type(PedigreeHolder) :: ped
            type(IndividualLinkedListNode),pointer :: ind

            integer :: i
            print *, "Ferdosi passed ", ped%sireList%length, " sires"
            print*, "Running Ferdosi for each sire..."
            ind => ped%sireList%first
            do i=1, ped%sireList%length
                call FerdosiRunner(ind%item, AllParameters%nSnp, AllParameters%overwriteHalfSibPhase)
                ind => ind%next
            enddo

            
        end subroutine RunFerdosiPerSire



        subroutine RunFerdosiPerDam(AllParameters, ped)

            use ModuleParameters
            use Ferdosi
            use PedigreeModule
            implicit none

            type(Parameters), intent(inout) :: AllParameters
            type(PedigreeHolder) :: ped
            type(IndividualLinkedListNode),pointer :: ind

            integer :: i
            print *, "Ferdosi passed ", ped%damList%length, " dams"
            print*, "Running Ferdosi for each dam..."
            ind => ped%damList%first
            do i=1, ped%damList%length
                call FerdosiRunner(ind%item, AllParameters%nSnp, AllParameters%overwriteHalfSibPhase)
                ind => ind%next
            enddo

            
        end subroutine RunFerdosiPerDam
        !###########################################################################################################################################################

 
 
end module ModuleRunFerdosi

!###########################################################################################################################################################
