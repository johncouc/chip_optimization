MODULE eigv_module
	IMPLICIT NONE
	SAVE
	
	PUBLIC eig
	INTERFACE eig
		MODULE PROCEDURE eig_sp
		MODULE PROCEDURE eig_dp
	END INTERFACE

	CONTAINS
		SUBROUTINE eig_sp(A, lambdar, lambdai)
			REAL(selected_real_kind(6,37)), DIMENSION(:,:), INTENT(in)				:: A
			REAL(selected_real_kind(6,37)), DIMENSION(SIZE(A,DIM=1), SIZE(A,DIM=2))	:: AA
			REAL(selected_real_kind(6,37)), DIMENSION(:),  INTENT(out)				:: lambdar, lambdai
			REAL(selected_real_kind(6,37)), DIMENSION(:,:), ALLOCATABLE				:: VLR
		  	REAL(selected_real_kind(6,37)), DIMENSION(:), ALLOCATABLE				:: work
			INTEGER info, l_work
			AA = A
			ALLOCATE(work(1))
			CALL SGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, -1, info)
			l_work = work(1)
			DEALLOCATE(work)
			ALLOCATE(work(l_work))
			CALL SGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, l_work, info)
			DEALLOCATE(work)			
		END SUBROUTINE

		SUBROUTINE eig_dp(A, lambdar, lambdai)
			REAL(selected_real_kind(15,307)), DIMENSION(:,:), INTENT(in)				:: A
			REAL(selected_real_kind(15,307)), DIMENSION(SIZE(A,DIM=1), SIZE(A,DIM=2))	:: AA
			REAL(selected_real_kind(15,307)), DIMENSION(:),  INTENT(out)				:: lambdar, lambdai
			REAL(selected_real_kind(15,307)), DIMENSION(:,:), ALLOCATABLE				:: VLR
		  	REAL(selected_real_kind(15,307)), DIMENSION(:), ALLOCATABLE					:: work
			INTEGER info, l_work
			AA = A
			ALLOCATE(work(1))
			CALL DGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, -1, info)
			l_work = work(1)
			DEALLOCATE(work)
			ALLOCATE(work(l_work))
			CALL DGEEV('N','N',SIZE(A, DIM=2), AA, SIZE(A, DIM=2), lambdar, lambdai, VLR, 1, VLR, 1, work, l_work, info)
			DEALLOCATE(work)			
		END SUBROUTINE
END MODULE
