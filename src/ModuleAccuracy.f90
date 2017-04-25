!###########################################################################################################################################################

module ModuleAccuracy
    
    implicit none

    !###########################################################################################################################################################

    contains
        
        !###########################################################################################################################################################


        subroutine CalculateAccuracyStats(AllParameters)

            use ModuleParameters
            use ModuleGet
            use pedigreeModule
            use ConstantModule
            
            implicit none

            type(Parameters), intent(inout) :: AllParameters
            integer :: i
            integer :: FileUnit, stat
            integer :: PosID, YieldP1, YieldP2, FirstHetPos
            character(len=IDLENGTH) :: InID
            double precision :: PhaseAccP1, PhaseAccP2
            real, allocatable, dimension(:,:) :: TruePhase, ImputedPhase, ImputedPhaseFormatted


            write(*,*), "Calculating sire phase accuracies..." 


            open(newunit=FileUnit, file=AllParameters%phaseFileSire, status="old")

            allocate(TruePhase(2,AllParameters%nSnp))
            allocate(ImputedPhaseFormatted(2,AllParameters%nSnp))
            allocate(ImputedPhase(2,AllParameters%nSnp))
            TruePhase = -99
            ImputedPhase = -99
            ImputedPhaseFormatted = -99
            
            do 
                TruePhase = -99
                ImputedPhase = -99
                ImputedPhaseFormatted = -99
                PosID = -99
                FirstHetPos = 1
                read(FileUnit, *, iostat=stat) InID, TruePhase(1,:)
                if (stat/=0) exit   
                read(FileUnit, *, iostat=stat) InID, TruePhase(2,:)
                if (stat/=0) exit
                
                call GetIDType(AllParameters%SireArray, AllParameters%nSire, InID, PosID)

                ImputedPhase(1,:) = AllParameters%SireArray(PosID)%ind%phaseInfo(:,1)
                ImputedPhase(2,:) = AllParameters%SireArray(PosID)%ind%phaseInfo(:,2)


                do i=1, AllParameters%nSnp
                    if ((SUM(TruePhase(:,i))==1) .and. (SUM(ImputedPhase(:,i))==1)) then
                        FirstHetPos = i 
                        exit 
                    endif
                enddo

                YieldP1 = 0
                PhaseAccP1 = 0.0
                YieldP2 = 0
                PhaseAccP2 = 0.0
                
               
                call DetermineSamePhase(ImputedPhase, TruePhase, ImputedPhaseFormatted, AllParameters%nSnp, FirstHetPos)

                YieldP1 = (AllParameters%nSnp - (count(ImputedPhaseFormatted(1,:)==-99)))
                call CalculateCorrelation(YieldP1, AllParameters%nSnp, TruePhase(1,:), ImputedPhaseFormatted(1,:), PhaseAccP1)
                YieldP2 = (AllParameters%nSnp - (count(ImputedPhaseFormatted(2,:)==-99)))
                call CalculateCorrelation(YieldP2, AllParameters%nSnp, TruePhase(2,:), ImputedPhaseFormatted(2,:), PhaseAccP2)
                AllParameters%SireArray(PosID)%YieldP1 = float(YieldP1) / float(AllParameters%nSnp)
                AllParameters%SireArray(PosID)%YieldP2 = float(YieldP2) / float(AllParameters%nSnp)
                AllParameters%SireArray(PosID)%phaseAccP1 = PhaseAccP1
                AllParameters%SireArray(PosID)%phaseAccP2 = PhaseAccP2

            enddo 

            deallocate(TruePhase)
            deallocate(ImputedPhase)
            deallocate(ImputedPhaseFormatted)

            close(FileUnit)

        end subroutine CalculateAccuracyStats



        !###########################################################################################################################################################


        subroutine CalculateCorrelation(Yield, nWantedSNp, TrueArray, ImputedArray, CorTrueImp)
            
            use AlphaStatMod

            implicit none
            
            
            integer,intent(in) :: nWantedSNp,Yield
            
            real,intent(in), dimension(:) :: ImputedArray(nWantedSNp), TrueArray(nWantedSNp)
            double precision, intent(inout) :: CorTrueImp
            
            integer :: i,p
            type(CorrelationReal32) :: CorTrueImpTemp
            real, allocatable,dimension(:) :: TrueTmp,ImpTmp
            
            
            if (Yield==nWantedSNp) then
                CorTrueImpTemp = Cor(TrueArray,ImputedArray)
            endif
            
            if (Yield.lt.nWantedSNp) then
                allocate(TrueTmp(Yield))
                allocate(ImpTmp(Yield))
                TrueTmp = -99
                ImpTmp = -99
                p=1
                do i=1,nWantedSNp
                    if (ImputedArray(i)/=-99) then
                        TrueTmp(p)=TrueArray(i)
                        ImpTmp(p)=ImputedArray(i)
                        p=p+1
                    endif
                enddo
                
                CorTrueImpTemp=Cor(TrueTmp,ImpTmp)

                if (allocated(TrueTmp)) deallocate(TrueTmp)
                if (allocated(ImpTmp)) deallocate(ImpTmp)
            endif
            
            CorTrueImp = DBLE(CorTrueImpTemp%Cor)

        end subroutine CalculateCorrelation

        !###########################################################################################################################################################

        subroutine DetermineSamePhase(ImputedPhaseIn, TruePhaseIn, ImputedPhaseOut, nSnpAll, FirstHeteroPos)

            implicit none

            integer, intent(in) :: nSnpAll, FirstHeteroPos
            real, intent(in), dimension(:,:) :: ImputedPhaseIn(2,nSnpAll), TruePhaseIn(2,nSnpAll)
            real, intent(inout), dimension(:,:) :: ImputedPhaseOut(2,nSnpAll)

            integer :: Phase1, Phase2

            Phase1 = 0
            Phase2 = 0
            ImputedPhaseOut = -99

            if (FirstHeteroPos /= 0 ) then
                if(ImputedPhaseIn(1,FirstHeteroPos)==TruePhaseIn(1,FirstHeteroPos)) then
                    Phase1 = 1
                    Phase2 = 2
                elseif (ImputedPhaseIn(1,FirstHeteroPos)==TruePhaseIn(2,FirstHeteroPos)) then
                    Phase1 = 2
                    Phase2 = 1
                else
                    Phase1 = 1
                    Phase2 = 2
                endif
            else 
                Phase1 = 1
                Phase2 = 2
            endif

            ImputedPhaseOut(Phase1,:) = ImputedPhaseIn(1,:)
            ImputedPhaseOut(Phase2,:) = ImputedPhaseIn(2,:)



        end subroutine DetermineSamePhase

        !###########################################################################################################################################################


end module ModuleAccuracy

!###########################################################################################################################################################
