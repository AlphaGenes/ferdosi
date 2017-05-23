!#############################################################################################################################################################################################################################
module ModuleGet

	implicit none

	!#############################################################################################################################################################################################################################
	contains

		!#############################################################################################################################################################################################################################
		subroutine GetIDType(InArray, InArrayDim, InID, InSeqID)

			use ModuleSireType
			use constantModule, only : dummyAnimalPrepre

			implicit none

			integer, intent(in) :: InArrayDim
			character(len=*), intent(in) :: InID
			integer, intent(inout) :: InSeqID
			type(Sire), dimension(:), intent(in) :: InArray(InArrayDim)
			
			integer :: i

			InSeqID = 0

			do i=1, InArrayDim
				if ((trim(InArray(i)%ind%originalId)) == trim(InID) .or. (trim(InArray(i)%ind%originalId)) == dummyAnimalPrepre//trim(InID) ) then
					InSeqID = i
				else
					continue
				endif 
				if (InSeqID /= 0) exit
			enddo

		end subroutine GetIDType

		!#############################################################################################################################################################################################################################
		
		subroutine GetIDInt(InArray, InArrayDim, InID, InSeqID)

			implicit none

			integer, intent(in) :: InArrayDim, InID
			integer, intent(inout) :: InSeqID
			integer, dimension(:), intent(in) :: InArray(InArrayDim)
			
			integer :: i

			InSeqID = 0

			do i=1, InArrayDim
				if (InArray(i) == InID) then
					InSeqID = i
				else
					continue
				endif 
				if (InSeqID /= 0) exit
			enddo

		end subroutine GetIDInt

		!#############################################################################################################################################################################################################################

end module ModuleGet

!#############################################################################################################################################################################################################################