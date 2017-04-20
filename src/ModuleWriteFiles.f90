!###########################################################################################################################################################

module ModuleWriteFiles
    
    implicit none

    !###########################################################################################################################################################

    contains
        
        !###########################################################################################################################################################


        subroutine WritePhase(AllParameters)

            use ModuleParameters

            implicit none

            type(Parameters), intent(in) :: AllParameters

            integer :: i
            integer :: FileUnit
            character(len=30) :: StrSnp, OutFmt


            print*, "Writing imputed sire phase..."

            open(newunit=FileUnit, file="ImputedSirePhase.txt", status="replace")

            write(StrSnp,*) AllParameters%nSnp
            OutFmt='(i20,'//trim(adjustl(StrSnp))//'f7.1)'

            do i=1, AllParameters%nSire
                write(FileUnit, OutFmt) AllParameters%SireArray(i)%ID, AllParameters%SireArray(i)%MyPhase(:,1)
                write(FileUnit, OutFmt) AllParameters%SireArray(i)%ID, AllParameters%SireArray(i)%MyPhase(:,2)
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
               write(FileUnit, '(1i20, 4f20.10)') AllParameters%SireArray(i)%ID, AllParameters%SireArray(i)%PhaseAccP1, AllParameters%SireArray(i)%PhaseAccP2, AllParameters%SireArray(i)%YieldP1, AllParameters%SireArray(i)%YieldP2
            enddo

            close(FileUnit)

        end subroutine WriteAccuracies

        !###########################################################################################################################################################


end module ModuleWriteFiles

!###########################################################################################################################################################
