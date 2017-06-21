!###########################################################################################################################################################

module ModuleWriteFiles
    
    implicit none

    !###########################################################################################################################################################

    contains
        
        !###########################################################################################################################################################


        subroutine WritePhase(AllParameters, ped)

            use ModuleParameters
            use PedigreeModule
            implicit none

            type(Parameters), intent(in) :: AllParameters
            type(PedigreeHolder), intent(in) :: ped
            type(IndividualLinkedListNode), pointer :: tmpSire
            integer :: i
            integer :: FileUnit
            character(len=30) :: StrSnp, OutFmt


            print*, "Writing imputed sire phase..."

            open(newunit=FileUnit, file="ImputedSirePhase.txt", status="replace")

            write(StrSnp,*) AllParameters%nSnp
            OutFmt='(a20,'//trim(adjustl(StrSnp))//'i2)'

            tmpSire => ped%sireList%first
            do i=1, ped%sireList%length
                write(FileUnit, OutFmt) tmpSire%item%originalId, tmpSire%item%individualPhase(1)%toIntegerArray()
                write(FileUnit, OutFmt) tmpSire%item%originalId, tmpSire%item%individualPhase(2)%toIntegerArray()
                tmpSire => tmpSire%next
            enddo

            close(FileUnit)
  
        end subroutine WritePhase

        !###########################################################################################################################################################


        subroutine WriteAccuracies(AllParameters)

            use ModuleParameters

            implicit none

            type(Parameters), intent(in) :: AllParameters

            integer :: i
            integer :: FileUnit


            print*, "Writing Sire Accuracies..."

            open(newunit=FileUnit, file="SirePhaseAccuracies.txt", status="replace")

            write(FileUnit,'(a86)') "SireID        Phase1Accuracy      Phase2Accuracy       Phase1Yield         Phase2Yield"

            do i=1, AllParameters%nSire
               write(FileUnit, '(a, 4f20.10)') AllParameters%SireArray(i)%ind%originalId, AllParameters%SireArray(i)%phaseAccP1, AllParameters%SireArray(i)%phaseAccP2, AllParameters%SireArray(i)%YieldP1, AllParameters%SireArray(i)%YieldP2
            enddo

            close(FileUnit)

        end subroutine WriteAccuracies

        !###########################################################################################################################################################


end module ModuleWriteFiles

!###########################################################################################################################################################
