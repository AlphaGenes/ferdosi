!###########################################################################################################################################################

module ModuleReadFiles
    
    implicit none

    !###########################################################################################################################################################

    contains
        
        !###########################################################################################################################################################


        subroutine ReadPedigree(AllParameters)

            use ModuleParameters
            use ModuleGet

            implicit none

            type(Parameters), intent(inout) :: AllParameters

            integer :: i
            integer :: fileunit, FileLength, stat
            integer :: InID, SireID, SirePosID
            integer :: SireCounter, nSnpTemp, OffspringCount

            allocate(AllParameters%SireArray(AllParameters%nSire))
            SireCounter = 0

            print*, "Reading Pedigree..."

            open(newunit=fileunit, file=AllParameters%PedigreeFile, status="old")

            FileLength = 0

            do
                read(fileunit, *, iostat=stat) 
                if (stat/=0) exit
                FileLength = FileLength + 1
            enddo  
            AllParameters%nOffspringTotal = FileLength
            allocate(AllParameters%Offspring2SireArray(4,FileLength))
            AllParameters%Offspring2SireArray = -99
            rewind(fileunit)

            do i=1, FileLength
                InID = 0
                SireID = 0
                SirePosID = 0 
                
                read(fileunit, *, iostat=stat) InID, SireID
                AllParameters%Offspring2SireArray(1,i) = InID
                if (ANY(AllParameters%Offspring2SireArray(2,:)==SireID)) then 
                    call GetIDType(AllParameters%SireArray, AllParameters%nSire, SireID, SirePosID)
                    AllParameters%Offspring2SireArray(2,i) = SireID
                    AllParameters%Offspring2SireArray(3,i) = SirePosID

                    AllParameters%SireArray(SirePosID)%nOffSpring = AllParameters%SireArray(SirePosID)%nOffSpring + 1
                    AllParameters%Offspring2SireArray(4,i) = AllParameters%SireArray(SirePosID)%nOffSpring 
                else
                    SireCounter = SireCounter + 1
                    AllParameters%SireArray(SireCounter)%ID = SireID
                    AllParameters%SireArray(SireCounter)%PosID = SireCounter
                    AllParameters%Offspring2SireArray(2,i) = SireID
                    AllParameters%Offspring2SireArray(3,i) = SireCounter
                    AllParameters%SireArray(SireCounter)%nOffSpring = 1
                    AllParameters%Offspring2SireArray(4,i) = AllParameters%SireArray(SireCounter)%nOffSpring 
                endif 
            enddo

            nSnpTemp = AllParameters%nSnp
            do i=1, AllParameters%nSire 
                OffspringCount = AllParameters%SireArray(i)%nOffSpring
                allocate(AllParameters%SireArray(i)%OffspringGenotypes(OffspringCount, nSnpTemp))
                allocate(AllParameters%SireArray(i)%OffspringPhase(OffspringCount, nSnpTemp, 2))
                allocate(AllParameters%SireArray(i)%MyPhase(nSnpTemp,2))
                AllParameters%SireArray(i)%OffspringGenotypes = MissingPhaseCode
                AllParameters%SireArray(i)%OffspringPhase = MissingPhaseCode
                AllParameters%SireArray(i)%MyPhase = MissingPhaseCode
                allocate(AllParameters%SireArray(i)%OffspringIds(OffspringCount))
                AllParameters%SireArray(i)%OffspringIds = -99
            enddo

            close(fileunit)

        end subroutine ReadPedigree

        !###########################################################################################################################################################

        subroutine ReadGenotypes(AllParameters)

            use ModuleParameters
            use ModuleGet

            implicit none

            type(Parameters), intent(inout) :: AllParameters

            integer :: i
            integer :: fileunit, stat, IDIn, PosID, FileLength, MySire, MySireOffSpring
            real, allocatable, dimension(:) :: InGenotypes

            allocate(InGenotypes(AllParameters%nSnp))

            IDIn = 0
            PosID = 0
            InGenotypes = MissingPhaseCode
            FileLength = 0

            print*, "Reading Genotypes..."

            open(newunit=fileunit, file=AllParameters%GenotypesFileOffspring, status="old")

            do
                read(fileunit, *, iostat=stat) 
                if (stat/=0) exit
                FileLength = FileLength + 1
            enddo  

            rewind(fileunit)
            
            do i=1, FileLength
                IDIn = 0
                PosID = 0
                InGenotypes = MissingPhaseCode
                MySire = MissingPhaseCode
                MySireOffSpring =MissingPhaseCode
                read(fileunit, *, iostat=stat) IDIn, InGenotypes(:)
                call GetIDInt(AllParameters%Offspring2SireArray(1,:), AllParameters%nOffspringTotal, IDIn, PosID) 
                MySire=AllParameters%Offspring2SireArray(3,PosID)
                MySireOffSpring=AllParameters%Offspring2SireArray(4,PosID)
                AllParameters%SireArray(MySire)%OffspringIds(MySireOffSpring) = IDIn
                AllParameters%SireArray(MySire)%OffspringGenotypes(MySireOffSpring,:) = InGenotypes
            enddo 

            deallocate(InGenotypes)
            close(fileunit)

        end subroutine ReadGenotypes

        !###########################################################################################################################################################

 
end module ModuleReadFiles

!###########################################################################################################################################################
