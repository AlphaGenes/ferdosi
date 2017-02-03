#ifdef OS_UNIX

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"

#else
#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#endif

!###########################################################################################################################################################

module Ferdosi

	implicit none

	public :: FerdosiRunner

	contains

	!###########################################################################################################################################################
		subroutine FerdosiRunner(nHalfSib, nSnp, ParentPhase, HalfSibGeno, HalfSibPhase)

			implicit none

			integer, intent(in) :: nHalfSib, nSnp
			integer, intent(in), dimension (:,:) :: HalfSibGeno(nHalfSib, nSnp)
			integer, intent(in), dimension(:,:,:) :: ParentPhase(1,nSnp,2), HalfSibPhase(nHalfSib,nSnp,2) !Phase arrays are dimension(nIndividuals, nSnp, 2) for phase 1 (paternal) and phase 2 (maternal)

			integer, allocatable, dimension(:) :: OpposingHomoPositions

			allocate(OpposingHomoPositions(nSnp))

			OpposingHomoPositions = 9

			!Part 1 of Ferdosi: Paternal strand detection based on opposing homozygous markers in half-sib families
			call GetOpposingHomozygoteSnp(nHalfSib,nSnp,HalfSibGeno,OpposingHomoPositions)
			call DefineCommonParentPhase(nHalfSib,nSnp,ParentPhase,HalfSibPhase,OpposingHomoPositions)


			!Part 2 of Ferdosi: Phasing and imputation of sire genotypes



			!Part 3 of Ferdosi: Phasing of genotypes of half-sib families




			deallocate(OpposingHomoPositions)


		end subroutine FerdosiRunner

	!###########################################################################################################################################################

		subroutine GetOpposingHomozygoteSnp(nHalfSib,nSnp,HalfSibGeno,OpposingHomoPositions)
			
			implicit none

			integer, intent(in) :: nHalfSib, nSnp
			integer, intent(in) dimension(:) :: OpposingHomoPositions(nSnp)
			integer, intent(in) dimension(:,:,:) :: HalfSibPhase(nHalfSib,nSnp) 

			integer :: i
			integer, allocatable, dimension(:) :: TempGenos

			allocate(TempGenos(nHalfSib))

			TempGenos = 9

			do i=1, nSnp
				TempGenos = 9
				TempGenos = HalfSibPhase(:,i)
				if((ANY(TempGenos==0)) .and. (ANY(TempGenos==2))) then
					OpposingHomoPositions(i) = 1
				endif
			enddo

		end subroutine

	!###########################################################################################################################################################

end module Ferdosi

!###########################################################################################################################################################

program FerdosiTest

	implicit none

	use Ferdosi

	call FerdosiRunner(nHalfSib, nSnp, ParentPhase, HalfSibGeno, HalfSibPhase)

end program FerdosiTest

