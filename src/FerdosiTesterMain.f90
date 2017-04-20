! #ifdef OS_UNIX
! #DEFINE DASH "/"
! #else
! #DEFINE DASH "\"
! #endif

!###########################################################################################################################################################

! FerdosiRunnerProgram is a simple program for imputation in crop species
! FerdosiRunnerProgram was written by Serap Gonen, start date April 2017
! You are welcome to use / modify FerdosiRunnerProgram in any way that you wish, provided that you credit the original 
! authors and share any improvements you make with us
! The latest version will be maintained at http://www.alphagenes.roslin.ed.ac.uk/alphasuite/
! Please report any bugs you find to serap.gonen@roslin.ed.ac.uk

!###########################################################################################################################################################

program FerdosiTesterMain

	use Ferdosi
	use ModuleSireType
    use ModuleReadFiles
    use ModuleAccuracy
    use ModuleParameters
    use ModuleGet
    use ModuleTLC
    use ModuleRunFerdosi
    use ModuleWriteFiles

	implicit none

	type(Parameters) :: AllParameters

	call ReadSpec(AllParameters, "FerdosiSpec.txt")

	call ReadPedigree(AllParameters)

	call ReadGenotypes(AllParameters)

	call RunFerdosiPerSire(AllParameters)

	call CalculateAccuracyStats(AllParameters)

	call WritePhase(AllParameters)
	
	call WriteAccuracies(AllParameters)

	
	!TODO deallocate all arrays

end program FerdosiTesterMain