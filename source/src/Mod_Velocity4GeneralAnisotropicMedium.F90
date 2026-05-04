#ifdef DOUBLE    
#define RP 8
#else
#define RP 4
#endif

#ifdef COMPLEX
#define DATATYPE complex(RP)
#else
#define DATATYPE real(RP)
#endif

  !> @file
  !> Contains the module Mod_Velocity4GeneralAnisotropicMedium.
  !

  !> Contains the routines for calculating phase and ray velocities.
  module Mod_Velocity4GeneralAnisotropicMedium
  implicit none
  private
#ifdef COMPLEX
  public:: solvePhaseDirection
#endif
  public:: qP, qSV, qSH, solveRayDirection, splitQsWaves, phaseAndRayVelocity
  public realKernal_a_q, realKernal_a_i, realKernal

  enum, bind(c)
    enumerator RAY_VELOCITY, RAY_SLOWNESS, RAY_ATTENUATION, IMAGINARY_VELOCITY
    enumerator RAY_QUALITY, PHASE_QUALITY, PHASE_ATTENUATION, PHASE_VELOCITY
  end enum

  enum, bind(c)
    !> qP wave
    enumerator:: qP  = 1
    !> qSV or qS1 wave
    enumerator:: qSV = 2
    !> qSH or qS2 wave
    enumerator:: qSH = 3
  end enum

  abstract interface
  pure function realKernalFun(n, a, kernal, cType) result(res)
  implicit none
  integer, intent(in):: n
  complex(RP), intent(in):: a(n)
  complex(RP), intent(in):: kernal(n)
  character(3), intent(in):: cType
  real(RP) res(n)
  end function
  end interface
  procedure(realKernalFun), pointer:: realKernal => null()

  contains
  !> Caculate the 6 upper triangular elements of Christoffel matrix Γ_jk = a_ijkl * n_i * n_l.
  pure function christoffelMatrix(a, n) result(res)
  implicit none
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> Phase velocity direction (normal vector).
  DATATYPE, intent(in):: n(3)
  !> @return 6 upper triangular elements.
  DATATYPE nn(6), res(6)
  !integer, parameter:: IND(3,3) = reshape([1,6,5,6,2,4,5,4,3],[3,3]) !convert Aijkl into Amn
  nn = [n(1)*n(1), n(1)*n(2), n(1)*n(3), n(2)*n(2), n(2)*n(3), n(3)*n(3)]
  res(1) = a(1,1)*nn(1) +   2.0_RP*a(1,6)*nn(2) +   2.0_RP*a(1,5)*nn(3) + a(6,6)*nn(4) +   2.0_RP*a(5,6)*nn(5) + a(5,5)*nn(6)
  res(2) = a(1,6)*nn(1) + (a(1,2)+a(6,6))*nn(2) + (a(1,4)+a(5,6))*nn(3) + a(2,6)*nn(4) + (a(2,5)+a(4,6))*nn(5) + a(4,5)*nn(6)
  res(3) = a(1,5)*nn(1) + (a(1,4)+a(5,6))*nn(2) + (a(1,3)+a(5,5))*nn(3) + a(4,6)*nn(4) + (a(3,6)+a(4,5))*nn(5) + a(3,5)*nn(6)
  res(4) = a(6,6)*nn(1) +   2.0_RP*a(2,6)*nn(2) +   2.0_RP*a(4,6)*nn(3) + a(2,2)*nn(4) +   2.0_RP*a(2,4)*nn(5) + a(4,4)*nn(6)
  res(5) = a(5,6)*nn(1) + (a(2,5)+a(4,6))*nn(2) + (a(3,6)+a(4,5))*nn(3) + a(2,4)*nn(4) + (a(2,3)+a(4,4))*nn(5) + a(3,4)*nn(6)
  res(6) = a(5,5)*nn(1) +   2.0_RP*a(4,5)*nn(2) +   2.0_RP*a(3,5)*nn(3) + a(4,4)*nn(4)  +  2.0_RP*a(3,4)*nn(5) + a(3,3)*nn(6)
  end function

  !> Calculate the partial derivative of 6 upper triangular elements of Christoffel matrix
  !> with respect to 3 components of phase velocity vector: ?Γ/?n.
  pure function christoffel_n(a, n) result(res)
  implicit none
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> Phase velocity direction (normal vector).
  DATATYPE, intent(in):: n(3)
  !> @return 3×6 partial derivatives.
  DATATYPE res(3,6), mat(6,3), B(6), n2(3)
  n2 = 2.0_RP * n
  mat(:,1) = [a(1,1), a(1,6), a(1,5), a(6,6), a(5,6), a(5,5)]
  mat(:,2) = [a(1,6), 0.5_RP*(a(1,2)+a(6,6)), 0.5_RP*(a(1,4)+a(5,6)), a(2,6), 0.5_RP*(a(2,5)+a(4,6)), a(4,5)]
  B(:)     = [a(1,5), 0.5_RP*(a(1,4)+a(5,6)), 0.5_RP*(a(1,3)+a(5,5)), a(4,6), 0.5_RP*(a(3,6)+a(4,5)), a(3,5)]
  mat(:,3) = B(:)
  res(1,:) = matmul(mat,n2) !偏n1
  mat(:,1) = mat(:,2)
  mat(:,2) = [a(6,6), a(2,6), a(4,6), a(2,2), a(2,4), a(4,4)]
  mat(:,3) = [a(5,6), 0.5_RP*(a(2,5)+a(4,6)), 0.5_RP*(a(3,6)+a(4,5)), a(2,4), 0.5_RP*(a(2,3)+a(4,4)), a(3,4)]
  res(2,:) = matmul(mat,n2) !偏n2
  mat(:,1) = B(:)
  mat(:,2) = mat(:,3)
  mat(:,3) = [a(5,5), a(4,5), a(3,5), a(4,4), a(3,4), a(3,3)]
  res(3,:) = matmul(mat,n2) !偏n3
  end function

  !> Calculate the partial derivative of 6 upper triangular elements of Christoffel matrix
  !> with respect to 6×6 density normalized elastic matrix: ?Γ/?a.
  pure function christoffel_a(n) result(res)
  implicit none
  !> Phase velocity direction (normal vector).
  DATATYPE, intent(in):: n(3)
  !> @return 6×6×6 partial derivatives.
  DATATYPE res(6,6,6), n11, n12, n13, n22, n23, n33
  integer i, j
  n11 = n(1) * n(1); n12 = n(1) * n(2); n13 = n(1) * n(3)
  n22 = n(2) * n(2); n23 = n(2) * n(3); n33 = n(3) * n(3)
  res = 0.0_RP
  !F11 的偏导
  res(1,1,1) = n11; res(1,6,1) = 2.0_RP*n12; res(1,5,1) = 2.0_RP*n13; res(6,6,1) = n22
  res(5,6,1) = 2.0_RP*n23; res(5,5,1) = n33
  !F12 的偏导
  res(1,6,2) = n11; res(1,2,2) = n12; res(6,6,2) = n12; res(1,4,2) = n13; res(5,6,2) = n13
  res(2,6,2) = n22; res(2,5,2) = n23; res(4,6,2) = n23; res(4,5,2) = n33
  !F13 的偏导
  res(1,5,3) = n11; res(1,4,3) = n12; res(5,6,3) = n12; res(1,3,3) = n13; res(5,5,3) = n13
  res(4,6,3) = n22; res(3,6,3) = n23; res(4,5,3) = n23; res(3,5,3) = n33
  !F22 的偏导
  res(6,6,4) = n11; res(2,6,4) = 2.0_RP*n12; res(4,6,4) = 2.0_RP*n13; res(2,2,4) = n22
  res(2,4,4) = 2.0_RP*n23; res(4,4,4) = n33
  !F23 的偏导
  res(5,6,5) = n11; res(2,5,5) = n12; res(4,6,5) = n12; res(3,6,5) = n13; res(4,5,5) = n13
  res(2,4,5) = n22; res(2,3,5) = n23; res(4,4,5) = n23; res(3,4,5) = n33
  !F33 的偏导
  res(5,5,6) = n11; res(4,5,6) = 2.0_RP*n12; res(3,5,6) = 2.0_RP*n13; res(4,4,6) = n22
  res(3,4,6) = 2.0_RP*n23; res(3,3,6) = n33
  !对称, 可删除
  do i = 1, 6
    do j = i+1, 6
      res(j,i,:) = res(i,j,:)
    end do
  end do
  end function

  !> Calculate the second-order partial derivative of 6 upper triangular elements of Christoffel matrix
  !> with respect to 3 components of phase velocity vector and to 6×6 density normalized elastic matrix: ?Γ/(?n ?a).
  pure function christoffel_n_a(n) result(res)
  implicit none
  !> Phase velocity direction (normal vector).
  DATATYPE, intent(in):: n(3)
  !> @return 6×6×3×6 partial derivatives.
  DATATYPE res(6,6,3,6), n1, n2, n3
  DATATYPE, parameter:: zero = 0.0_RP
  integer i, j
  n1 = n(1); n2 = n(2); n3 = n(3)
  res = 0.0_RP
  !F11 的偏导
  res(1,1,:,1) = [n1, zero, zero]; res(1,5,:,1) = [n3, zero, n1]; res(1,6,:,1) = [n2, n1, zero]
  res(5,5,:,1) = [zero, zero, n3]; res(5,6,:,1) = [zero, n3, n2]; res(6,6,:,1) = [zero, n2, zero]
  res(:,:,:,1) = res(:,:,:,1) * 2.0_RP
  !F12 的偏导
  res(1,2,:,2) = [n2, n1, zero]; res(1,4,:,2) = [n3, zero, n1];          res(1,6,:,2) = [n1*2.0_RP, zero, zero]
  res(2,5,:,2) = [zero, n3, n2]; res(2,6,:,2) = [zero, n2*2.0_RP, zero]; res(4,5,:,2) = [zero, zero, n3*2.0_RP]
  res(4,6,:,2) = [zero, n3, n2]; res(5,6,:,2) = [n3, zero, n1];          res(6,6,:,2) = [n2, n1, zero]
  !F13 的偏导
  res(1,3,:,3) = [n3, zero, n1];          res(1,4,:,3) = [n2, n1, zero]; res(1,5,:,3) = [n1*2.0_RP, zero, zero]
  res(3,5,:,3) = [zero, zero, n3*2.0_RP]; res(3,6,:,3) = [zero, n3, n2]; res(4,5,:,3) = [zero, n3, n2]
  res(4,6,:,3) = [zero, n2*2.0_RP, zero]; res(5,5,:,3) = [n3, zero, n1]; res(5,6,:,3) = [n2, n1, zero]
  !F22 的偏导
  res(2,2,:,4) = [zero, n2, zero]; res(2,4,:,4) = [zero, n3, n2]; res(2,6,:,4) = [n2, n1, zero]
  res(4,4,:,4) = [zero, zero, n3]; res(4,6,:,4) = [n3, zero, n1]; res(6,6,:,4) = [n1, zero, zero]
  res(:,:,:,4) = res(:,:,:,4) * 2.0_RP
  !F23 的偏导
  res(2,3,:,5) = [zero, n3, n2];          res(2,4,:,5) = [zero, n2*2.0_RP, zero]; res(2,5,:,5) = [n2, n1, zero]
  res(3,4,:,5) = [zero, zero, n3*2.0_RP]; res(3,6,:,5) = [n3, zero, n1];          res(4,4,:,5) = [zero, n3, n2]
  res(4,5,:,5) = [n3, zero, n1];          res(4,6,:,5) = [n2, n1, zero];          res(5,6,:,5) = [n1*2.0_RP, zero, zero]
  !F33 的偏导
  res(3,3,:,6) = [zero, zero, n3]; res(3,4,:,6) = [zero, n3, n2]; res(3,5,:,6) = [n3, zero, n1]
  res(4,4,:,6) = [zero, n2, zero]; res(4,5,:,6) = [n2, n1, zero]; res(5,5,:,6) = [n1, zero, zero]
  res(:,:,:,6) = res(:,:,:,6) * 2.0_RP
  !对称, 可删除
  do i = 1, 6
    do j = i+1, 6
      res(j,i,:,:) = res(i,j,:,:)
    end do
  end do
  end function

  !> Calculate the second-order partial derivative of 6 upper triangular elements of Christoffel matrix
  !> with respect to components of phase velocity vector n_i and n_j: ?Γ/(?n_i ?n_j).
  pure function christoffel_ni_nj(a) result(res)
  implicit none
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> @return 3×3×6 partial derivatives.
  DATATYPE res(3,3,6)
  !F11 的偏导
  res(:,1,1) = [a(1,1), a(1,6), a(1,5)]
  res(:,2,1) = [a(1,6), a(6,6), a(5,6)]
  res(:,3,1) = [a(1,5), a(5,6), a(5,5)]
  res(:,:,1) = res(:,:,1) * 2
  !F12 的偏导
  res(:,1,2) = [2*a(1,6), a(1,2)+a(6,6), a(1,4)+a(5,6)]
  res(:,2,2) = [a(1,2)+a(6,6), 2*a(2,6), a(2,5)+a(4,6)]
  res(:,3,2) = [a(1,4)+a(5,6), a(2,5)+a(4,6), 2*a(4,5)]
  !F13 的偏导
  res(:,1,3) = [2*a(1,5), a(1,4)+a(5,6), a(1,3)+a(5,5)]
  res(:,2,3) = [a(1,4)+a(5,6), 2*a(4,6), a(3,6)+a(4,5)]
  res(:,3,3) = [a(1,3)+a(5,5), a(3,6)+a(4,5), 2*a(3,5)]
  !F22 的偏导
  res(:,1,4) = [a(6,6), a(2,6), a(4,6)]
  res(:,2,4) = [a(2,6), a(2,2), a(2,4)]
  res(:,3,4) = [a(4,6), a(2,4), a(4,4)]
  res(:,:,4) = res(:,:,4) * 2
  !F23 的偏导
  res(:,1,5) = [2*a(5,6), a(2,5)+a(4,6), a(3,6)+a(4,5)]
  res(:,2,5) = [a(2,5)+a(4,6), 2*a(2,4), a(2,3)+a(4,4)]
  res(:,3,5) = [a(3,6)+a(4,5), a(2,3)+a(4,4), 2*a(3,4)]
  !F33 的偏导
  res(:,1,6) = [a(5,5), a(4,5), a(3,5)]
  res(:,2,6) = [a(4,5), a(4,4), a(3,4)]
  res(:,3,6) = [a(3,5), a(3,4), a(3,3)]
  res(:,:,6) = res(:,:,6) * 2
  end function

  !> Calculate the coefficients B, C, D of the cubic function of one variable:
  !> λ**3 + Bλ**2 + Cλ + D = 0.
  pure subroutine BCD(cstl, B, C, D)
  implicit none
  !> 6 upper triangular elements of Christoffel matrix.
  DATATYPE, intent(in) :: cstl(6)
  !> Coefficient B.
  DATATYPE, intent(out):: B
  !> Coefficient C.
  DATATYPE, intent(out):: C
  !> Coefficient D.
  DATATYPE, intent(out):: D
  B = -(cstl(1)+cstl(4)+cstl(6))
  C = (cstl(1)+cstl(4))*cstl(6) + cstl(1)*cstl(4) - (cstl(2)*cstl(2) + cstl(3)*cstl(3) + cstl(5)*cstl(5))
  D = cstl(1)*cstl(5)*cstl(5) + cstl(4)*cstl(3)*cstl(3) + cstl(6)*cstl(2)*cstl(2) - &
    (cstl(1)*cstl(4)*cstl(6) + 2.0_RP*cstl(2)*cstl(3)*cstl(5))
  end subroutine

  !> Caculate the partial derivative of coefficients B, C, D
  !> with respect to 3 components of phase velocity vector.
  pure subroutine BCD_n(cstl, cstl_n, B_n, C_n, D_n)
  implicit none
  !> 6 upper triangular elements of Christoffel matrix.
  DATATYPE, intent(in) :: cstl(6)
  !> ?Γ/?n.
  DATATYPE, intent(in) :: cstl_n(3,6)
  !> Partial derivative ?B/?n.
  DATATYPE, intent(out):: B_n(3)
  !> Partial derivative ?C/?n.
  DATATYPE, intent(out):: C_n(3)
  !> Partial derivative ?D/?n.
  DATATYPE, intent(out):: D_n(3)
  B_n(:) = -(cstl_n(:,1) + cstl_n(:,4) + cstl_n(:,6))
  C_n(:) = (cstl(4)+cstl(6))*cstl_n(:,1) + (cstl(1)+cstl(6))*cstl_n(:,4) + (cstl(1)+cstl(4))*cstl_n(:,6) - &
    2.0_RP*(cstl(2)*cstl_n(:,2) + cstl(3)*cstl_n(:,3) + cstl(5)*cstl_n(:,5))
  D_n(:) = (cstl(5)*cstl(5)-cstl(4)*cstl(6))*cstl_n(:,1) + 2.0_RP*(cstl(2)*cstl(6)-cstl(3)*cstl(5))*cstl_n(:,2) + &
    2.0_RP*(cstl(3)*cstl(4)-cstl(2)*cstl(5))*cstl_n(:,3) + (cstl(3)*cstl(3)-cstl(1)*cstl(6))*cstl_n(:,4) + &
    2.0_RP*(cstl(1)*cstl(5)-cstl(2)*cstl(3))*cstl_n(:,5) + (cstl(2)*cstl(2)-cstl(1)*cstl(4))*cstl_n(:,6)
  end subroutine

  !> Calculate the partial derivative of coefficients B, C, D
  !> with respect to density normalized elastic matrix: ?B/?a, ?C/?a, ?D/?a.
  pure subroutine BCD_a(cstl, cstl_a, B_a, C_a, D_a)
  implicit none
  !> 6 upper triangular elements of Christoffel matrix.
  DATATYPE, intent(in) :: cstl(6)
  !> ?Γ/?a.
  DATATYPE, intent(in) :: cstl_a(6,6,6)
  !> Partial derivative ?B/?a.
  DATATYPE, intent(out):: B_a(6,6)
  !> Partial derivative ?C/?a.
  DATATYPE, intent(out):: C_a(6,6)
  !> Partial derivative ?D/?a.
  DATATYPE, intent(out):: D_a(6,6)
  B_a(:,:) = -(cstl_a(:,:,1) + cstl_a(:,:,4) + cstl_a(:,:,6))
  C_a(:,:) = (cstl(4)+cstl(6))*cstl_a(:,:,1) + (cstl(1)+cstl(6))*cstl_a(:,:,4) + (cstl(1)+cstl(4))*cstl_a(:,:,6) - &
    2.0_RP*(cstl(2)*cstl_a(:,:,2) + cstl(3)*cstl_a(:,:,3) + cstl(5)*cstl_a(:,:,5))
  D_a(:,:) = (cstl(5)*cstl(5)-cstl(4)*cstl(6))*cstl_a(:,:,1) + 2.0_RP*(cstl(2)*cstl(6)-cstl(3)*cstl(5))*cstl_a(:,:,2) + &
    2.0_RP*(cstl(3)*cstl(4)-cstl(2)*cstl(5))*cstl_a(:,:,3) + (cstl(3)*cstl(3)-cstl(1)*cstl(6))*cstl_a(:,:,4) + &
    2.0_RP*(cstl(1)*cstl(5)-cstl(2)*cstl(3))*cstl_a(:,:,5) + (cstl(2)*cstl(2)-cstl(1)*cstl(4))*cstl_a(:,:,6)
  end subroutine

  !> Calculate the second-order partial derivative of coefficients B, C, D
  !> with respect to phase velocity vector n and density normalized elastic matrix a.
  pure subroutine BCD_n_a(cstl, cstl_n, cstl_a, cstl_n_a, B_n_a, C_n_a, D_n_a)
  implicit none
  !> 6 upper triangular elements of Christoffel matrix.
  DATATYPE, intent(in) :: cstl(6)
  !> ?Γ/?n.
  DATATYPE, intent(in) :: cstl_n(3,6)
  !> ?Γ/?a.
  DATATYPE, intent(in) :: cstl_a(6,6,6)
  !> Second-order partial derivative ?Γ/(?n ?a).
  DATATYPE, intent(in) :: cstl_n_a(6,6,3,6)
  !> Second-order partial derivative ?B/(?n ?a).
  DATATYPE, intent(out):: B_n_a(6,6,3)
  !> Second-order partial derivative ?C/(?n ?a).
  DATATYPE, intent(out):: C_n_a(6,6,3)
  !> Second-order partial derivative ?D/(?n ?a).
  DATATYPE, intent(out):: D_n_a(6,6,3)
  integer i
  B_n_a(:,:,:) = -(cstl_n_a(:,:,:,1) + cstl_n_a(:,:,:,4) + cstl_n_a(:,:,:,6))
  do concurrent(i=1:3)
    C_n_a(:,:,i) = (cstl_n(i,4)+cstl_n(i,6))*cstl_a(:,:,1) + (cstl(4)+cstl(6))*cstl_n_a(:,:,i,1) + &
      (cstl_n(i,1)+cstl_n(i,6))*cstl_a(:,:,4) + (cstl(1)+cstl(6))*cstl_n_a(:,:,i,4) + &
      (cstl_n(i,1)+cstl_n(i,4))*cstl_a(:,:,6) + (cstl(1)+cstl(4))*cstl_n_a(:,:,i,6) - &
      2.0_RP*( cstl_n(i,2)*cstl_a(:,:,2) + cstl(2)*cstl_n_a(:,:,i,2) + &
      cstl_n(i,3)*cstl_a(:,:,3) + cstl(3)*cstl_n_a(:,:,i,3) + &
      cstl_n(i,5)*cstl_a(:,:,5) + cstl(5)*cstl_n_a(:,:,i,5) )
    D_n_a(:,:,i) = (2.0_RP*cstl(5)*cstl_n(i,5)-cstl(6)*cstl_n(i,4)-cstl(4)*cstl_n(i,6))*cstl_a(:,:,1) + (cstl(5)*cstl(5)-cstl(4)*cstl(6))*cstl_n_a(:,:,i,1) + &
      (2.0_RP*cstl(3)*cstl_n(i,3)-cstl(6)*cstl_n(i,1)-cstl(1)*cstl_n(i,6))*cstl_a(:,:,4) + (cstl(3)*cstl(3)-cstl(1)*cstl(6))*cstl_n_a(:,:,i,4) + &
      (2.0_RP*cstl(2)*cstl_n(i,2)-cstl(4)*cstl_n(i,1)-cstl(1)*cstl_n(i,4))*cstl_a(:,:,6) + (cstl(2)*cstl(2)-cstl(1)*cstl(4))*cstl_n_a(:,:,i,6) + &
      2.0_RP*(cstl(6)*cstl_n(i,2)+cstl(2)*cstl_n(i,6)-cstl(5)*cstl_n(i,3)-cstl(3)*cstl_n(i,5))*cstl_a(:,:,2) + &
      2.0_RP*(cstl(2)*cstl(6)-cstl(3)*cstl(5))*cstl_n_a(:,:,i,2) + &
      2.0_RP*(cstl(4)*cstl_n(i,3)+cstl(3)*cstl_n(i,4)-cstl(5)*cstl_n(i,2)-cstl(2)*cstl_n(i,5))*cstl_a(:,:,3) + &
      2.0_RP*(cstl(3)*cstl(4)-cstl(2)*cstl(5))*cstl_n_a(:,:,i,3) + &
      2.0_RP*(cstl(5)*cstl_n(i,1)+cstl(1)*cstl_n(i,5)-cstl(3)*cstl_n(i,2)-cstl(2)*cstl_n(i,3))*cstl_a(:,:,5) + &
      2.0_RP*(cstl(1)*cstl(5)-cstl(2)*cstl(3))*cstl_n_a(:,:,i,5)
  end do
  end subroutine

  !> Calculate the second-order partial derivative: ?B/(?n_i ?n_j), ?C/(?n_i ?n_j), ?D/(?n_i ?n_j).
  pure subroutine BCD_ni_nj(cstl, cstl_n, cstl_ni_nj, B_ni_nj, C_ni_nj, D_ni_nj)
  implicit none
  !> 6 upper triangular elements of Christoffel matrix.
  DATATYPE, intent(in) :: cstl(6)
  !> ?Γ/?n.
  DATATYPE, intent(in) :: cstl_n(3,6)
  !> Second-order partial derivative ?Γ/(?n_i ?n_j).
  DATATYPE, intent(in) :: cstl_ni_nj(3,3,6)
  !> Second-order partial derivative ?B/(?n_i ?n_j).
  DATATYPE, intent(out):: B_ni_nj(3,3)
  !> Second-order partial derivative ?C/(?n_i ?n_j).
  DATATYPE, intent(out):: C_ni_nj(3,3)
  !> Second-order partial derivative ?D/(?n_i ?n_j).
  DATATYPE, intent(out):: D_ni_nj(3,3)
  integer i, j
  B_ni_nj(:,:) = -(cstl_ni_nj(:,:,1) + cstl_ni_nj(:,:,4) + cstl_ni_nj(:,:,6))
  do concurrent(i=1:3, j=1:3)
    C_ni_nj(j,i) = (cstl_n(i,4)+cstl_n(i,6))*cstl_n(j,1) + (cstl(4)+cstl(6))*cstl_ni_nj(j,i,1) + &
      (cstl_n(i,1)+cstl_n(i,6))*cstl_n(j,4) + (cstl(1)+cstl(6))*cstl_ni_nj(j,i,4) + &
      (cstl_n(i,1)+cstl_n(i,4))*cstl_n(j,6) + (cstl(1)+cstl(4))*cstl_ni_nj(j,i,6) - &
      2.0_RP*( cstl_n(i,2)*cstl_n(j,2) + cstl(2)*cstl_ni_nj(j,i,2) + &
      cstl_n(i,3)*cstl_n(j,3) + cstl(3)*cstl_ni_nj(j,i,3) + &
      cstl_n(i,5)*cstl_n(j,5) + cstl(5)*cstl_ni_nj(j,i,5) )
    D_ni_nj(j,i) = (2.0_RP*cstl(5)*cstl_n(i,5)-cstl(6)*cstl_n(i,4)-cstl(4)*cstl_n(i,6))*cstl_n(j,1) + (cstl(5)*cstl(5)-cstl(4)*cstl(6))*cstl_ni_nj(j,i,1) + &
      (2.0_RP*cstl(3)*cstl_n(i,3)-cstl(6)*cstl_n(i,1)-cstl(1)*cstl_n(i,6))*cstl_n(j,4) + (cstl(3)*cstl(3)-cstl(1)*cstl(6))*cstl_ni_nj(j,i,4) + &
      (2.0_RP*cstl(2)*cstl_n(i,2)-cstl(4)*cstl_n(i,1)-cstl(1)*cstl_n(i,4))*cstl_n(j,6) + (cstl(2)*cstl(2)-cstl(1)*cstl(4))*cstl_ni_nj(j,i,6) + &
      2.0_RP*(cstl(6)*cstl_n(i,2)+cstl(2)*cstl_n(i,6)-cstl(5)*cstl_n(i,3)-cstl(3)*cstl_n(i,5))*cstl_n(j,2) + &
      2.0_RP*(cstl(2)*cstl(6)-cstl(3)*cstl(5))*cstl_ni_nj(j,i,2) + &
      2.0_RP*(cstl(4)*cstl_n(i,3)+cstl(3)*cstl_n(i,4)-cstl(5)*cstl_n(i,2)-cstl(2)*cstl_n(i,5))*cstl_n(j,3) + &
      2.0_RP*(cstl(3)*cstl(4)-cstl(2)*cstl(5))*cstl_ni_nj(j,i,3) + &
      2.0_RP*(cstl(5)*cstl_n(i,1)+cstl(1)*cstl_n(i,5)-cstl(3)*cstl_n(i,2)-cstl(2)*cstl_n(i,3))*cstl_n(j,5) + &
      2.0_RP*(cstl(1)*cstl(5)-cstl(2)*cstl(3))*cstl_ni_nj(j,i,5)
  end do
  end subroutine

  !> Solve the cubic function of one variable:
  !> λ**3 + Bλ**2 + Cλ + D = 0.
  pure function unaryCubicEquation(B, C, D) result(res)
  implicit none
  complex(RP), parameter:: w = cmplx(-0.5_RP,sqrt(0.75_RP),RP), w2 = w*w
  !> Coefficients B, C, D.
  complex(RP), intent(in):: B, C, D
  !> @return 3 complex roots.
  complex(RP) res(3), u, v, P, Q
  integer sig
  sig = 1
  P = B*C/6.0_RP - B**3/27.0_RP - 0.5_RP*D
  Q = (3.0_RP*C-B*B) / 9.0_RP
  v = sqrt(P*P+Q**3)
  !判断正负, 计算 u
  u = P + v
  v = P - v
  if(abs(v)>abs(u)) then
    u = v
    sig = -1
  end if
  u = u**(1.0_RP/3.0_RP)
  !计算 v
  if(abs(u)>1.0e-15_RP) then
    v = -Q/u
  else
    v = 0.0_RP
    u = 0.0_RP
  end if
  !3个解
  res = -B / 3.0_RP
  res(1) = res(1) + u    + v
  res(2) = res(2) + u*w2 + v*w
  res(3) = res(3) + u*w  + v*w2
  if( any(abs(res**3 + B*res*res + C*res + D)>1.0e-10) ) error stop 'error causes in FUNCTION unaryCubicEquation'
  end function unaryCubicEquation

  !> Calculate the 3 eigenvalues of 3×3 Christoffel matrix by mkl.
  pure function eigenValue(cstl) result(res)
  use lapack95, only: syev, geev
  !> 6 upper triangular elements of Christoffel matrix.
  DATATYPE, intent(in):: cstl(6)
  !> 3 real/complex eigen values
  DATATYPE res(3), a(3,3), temp(3,3)
  integer i
  a(1,:) = cstl(1:3)
  a(2,:) = [cstl(2), cstl(4:5)]
  a(3,:) = [cstl(3), cstl(5), cstl(6)]
#ifdef COMPLEX
  call geev(a,w=res,vr=temp)
#else
  call syev(a, res, "V", "L")
#endif
  !排序,第一个最大
  select case(maxloc(abs(res),dim=1))
  case(2)
    res = cshift(res,1)
  case(3)
    res = cshift(res,-1)
  end select
  end function

  !> Calculate phase and ray velocities and their partial derivatives by analytic fomula.
  pure subroutine phaseAndRayVelocity(a, n, pv, rv, vi, pv_a, rv_a, vi_nj)
  implicit none
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> Phase velocity direction (normal vector).
  DATATYPE, intent(in):: n(3)
  !> 3 phase velocities.
  DATATYPE, intent(out):: pv(3)
  !> 3 ray velocities.
  DATATYPE, intent(out):: rv(3)
  !> 3 ray velocity vectors.
  DATATYPE, intent(out):: vi(3,3)
  !> \a optional, partial derivative: ?pv/?a.
  DATATYPE, intent(out), optional:: pv_a(6,6,3)
  !> \a optional, partial derivative: ?rv/?a.
  DATATYPE, intent(out), optional:: rv_a(6,6,3)
  !> \a optional, partial derivative: ?vi/?n_j.
  DATATYPE, intent(out), optional:: vi_nj(3,3,3)
  real(RP), parameter:: EPS = 1.0e-7_RP
  DATATYPE cstl(6), cstl_n(3,6), cstl_a(6,6,6), cstl_n_a(6,6,3,6), vi_a(6,6,3,3)
  DATATYPE B, C, D, B_n(3), C_n(3), D_n(3), B_a(6,6), C_a(6,6), D_a(6,6), B_n_a(6,6,3), C_n_a(6,6,3), D_n_a(6,6,3)
  DATATYPE cstl_ni_nj(3,3,6), B_ni_nj(3,3), C_ni_nj(3,3), D_ni_nj(3,3)
  DATATYPE g(3), c2(3), c4(3)
  integer nRoot
  integer i, j
  !计算christoffel矩阵、多项式系数
  cstl = christoffelMatrix(a, n)
  call BCD(cstl, B, C, D)
  !计算偏 n
  cstl_n = christoffel_n(a, n)
  call BCD_n(cstl, cstl_n, B_n, C_n, D_n)
  !计算偏 a
  if(present(pv_a) .or. present(rv_a)) then
    cstl_a = christoffel_a(n)
    call BCD_a(cstl, cstl_a, B_a, C_a, D_a)
    !计算偏 n 偏 a
    if(present(rv_a)) then
      cstl_n_a = christoffel_n_a(n)
      call BCD_n_a(cstl, cstl_n, cstl_a, cstl_n_a, B_n_a, C_n_a, D_n_a)
    end if
  end if

#ifdef MKL
  c2(:) = eigenValue(cstl)
#else 
  !求相速度
#ifdef COMPLEX
  c2(:) = unaryCubicEquation(B, C, D)
#else
  !强制数据类型转换为 complex
  c2(:) = unaryCubicEquation(cmplx(B,0.0_RP,RP), cmplx(C,0.0_RP,RP), cmplx(D,0.0_RP,RP))
#endif
#endif

  !区分qSV和qSH, 可删除
  if(abs(a(6,6)*(n(1)*n(1)+n(2)*n(2))+a(4,4)*n(3)*n(3)-c2(2))<EPS) then
    !if(c2(3)>1.6**2) then
    c4(1) = c2(2); c2(2) = c2(3); c2(3) = c4(1)
  end if
  c4(:) = c2(:) * c2(:) !相速度的高次
  pv = sqrt(c2) !相速度
  !判断根的个数
  if(abs(pv(2)-pv(3))<EPS) then
    nRoot = 2
    if(abs(pv(1)-pv(3))<EPS) nRoot = 1
  else
    nRoot = 3
  end if

  !计算速度
  select case(nRoot)
  case(1)
    g(:) = 6.0_RP*pv(:) !分母, 6c
    do concurrent(i=1:3)
      vi(:,i) = -B_n(:) / g(i)
    end do
  case(2)
    g(1) = 2.0_RP*pv(1)*(3.0_RP*c4(1) + 2.0_RP*B*c2(1) + C) !分母, 6c**5+4Bc**3+2Cc
    vi(:,1) = -(c4(1)*B_n(:) + c2(1)*C_n(:) + D_n(:)) / g(1)
    g(2) = 4.0_RP*pv(2)*(3.0_RP*c2(2) + B) !分母, 12c**3+4Bc
    vi(:,2) = -(2.0_RP*c2(2)*B_n(:) + C_n(:)) / g(2)
    g(3) = g(2); vi(:,3) = vi(:,2)
  case(3)
    g(:) = 2.0_RP*pv(:)*(3.0_RP*c4(:) + 2.0_RP*B*c2(:) + C) !分母, 6c**5+4Bc**3+2Cc
    do concurrent(i=1:3)
      vi(:,i) = -(c4(i)*B_n(:) + c2(i)*C_n(:) + D_n(:)) / g(i)
    end do
  end select
  do concurrent(i=1:3)
    rv(i) = sqrt(vi(1,i)*vi(1,i)+vi(2,i)*vi(2,i)+vi(3,i)*vi(3,i))
  end do
  !######### 速度计算完毕, 以下可选 #########

  !计算相速度偏 a
  if(present(pv_a) .or. present(rv_a)) then
    select case(nRoot)
    case(1)
      do concurrent(i=1:3)
        pv_a(:,:,i) = -B_a(:,:) / g(i)
      end do
    case(2)
      pv_a(:,:,1) = -(c4(1)*B_a(:,:) + c2(1)*C_a(:,:) + D_a(:,:)) / g(1)
      pv_a(:,:,2) = -(2.0_RP*c2(2)*B_a(:,:) + C_a(:,:)) / g(2)
      pv_a(:,:,3) = pv_a(:,:,2)
    case(3)
      do concurrent(i=1:3)
        pv_a(:,:,i) = -(c4(i)*B_a(:,:) + c2(i)*C_a(:,:) + D_a(:,:)) / g(i)
      end do
    end select
  end if

  !计算群速度分量偏 a
  if(present(rv_a)) then
    vi_a = 0.0_RP
    select case(nRoot)
    case(1)
      do j = 1, 3 !波类型
        do i = 1, 3 !分量
          vi_a(:,:,i,j) = -(6.0_RP*pv_a(:,:,j)*vi(i,j)+B_n_a(:,:,i)) / g(j)
        end do
      end do
    case(2)
      j = 1 !qP
      do i = 1, 3 !分量
        vi_a(:,:,i,j) = (30.0_RP*c4(j)+12.0_RP*B*c2(j)+2.0_RP*C)*pv_a(:,:,j) + 2.0_RP*pv(j)*(2.0_RP*c2(j)*B_a(:,:)+C_a(:,:))
        vi_a(:,:,i,j) = vi_a(:,:,i,j) * vi(i,j) !第一项
        vi_a(:,:,i,j) = vi_a(:,:,i,j) + 2.0_RP*pv(j)*pv_a(:,:,j)*(2.0_RP*c2(j)*B_n(i)+C_n(i)) + c4(j)*B_n_a(:,:,i) + &
          c2(j)*C_n_a(:,:,i) + D_n_a(:,:,i) !两项之和
      end do
      vi_a(:,:,:,j) = -vi_a(:,:,:,j) / g(j) !除以分母
      j = 2 !qS1
      do i = 1, 3 !分量
        vi_a(:,:,i,j) = (36.0_RP*c2(j)+4.0_RP*B)*pv_a(:,:,j) + 4.0_RP*pv(j)*B_a(:,:)
        vi_a(:,:,i,j) = vi_a(:,:,i,j) * vi(i,j) !第一项
        vi_a(:,:,i,j) = vi_a(:,:,i,j) + 4.0_RP*pv(j)*pv_a(:,:,j)*B_n(i) + 2.0_RP*c2(j)*B_n_a(:,:,i) + C_n_a(:,:,i) !两项之和
      end do
      vi_a(:,:,:,j) = -vi_a(:,:,:,j) / g(j) !除以分母
      vi_a(:,:,:,3) = vi_a(:,:,:,2) !qS2
    case(3)
      do concurrent( j = 1: 3) !波类型
        do i = 1, 3 !分量
          vi_a(:,:,i,j) = (30.0_RP*c4(j)+12.0_RP*B*c2(j)+2.0_RP*C)*pv_a(:,:,j) + 2.0_RP*pv(j)*(2.0_RP*c2(j)*B_a(:,:)+C_a(:,:))
          vi_a(:,:,i,j) = vi_a(:,:,i,j) * vi(i,j) !第一项
          vi_a(:,:,i,j) = vi_a(:,:,i,j) + 2.0_RP*pv(j)*pv_a(:,:,j)*(2.0_RP*c2(j)*B_n(i)+C_n(i)) + c4(j)*B_n_a(:,:,i) + &
            c2(j)*C_n_a(:,:,i) + D_n_a(:,:,i) !两项之和
        end do
        vi_a(:,:,:,j) = -vi_a(:,:,:,j) / g(j) !除以分母
      end do
    end select
    !群速度偏 a
    rv_a = 0.0_RP
    do concurrent( j = 1: 3) !波类型
      do i = 1, 3 !分量
        rv_a(:,:,j) = rv_a(:,:,j) + vi(i,j)*vi_a(:,:,i,j)
      end do
      rv_a(:,:,j) = rv_a(:,:,j) / rv(j)
    end do
  end if
  !######### 相、群速度偏 a 计算完毕, 以下可选 #########

  !计算群速度分量vi偏 nj
  if(present(vi_nj)) then
    cstl_ni_nj = christoffel_ni_nj(a)
    call BCD_ni_nj(cstl, cstl_n, cstl_ni_nj, B_ni_nj, C_ni_nj, D_ni_nj)
    select case(nRoot)
    case(1)
      do concurrent( i=1:3, j=1:3)
        vi_nj(j,i,:) = -(6.0_RP*vi(i,:)*vi(j,:)+B_ni_nj(j,i)) / g(:)
      end do
    case(2)
      do concurrent( i=1:3, j=1:3)
        vi_nj(j,i,:) = (36.0_RP*c2(:)+4.0_RP*B)*vi(j,:) + 4.0_RP*pv(:)*B_n(j)
        vi_nj(j,i,:) = vi_nj(j,i,:) * vi(i,:) !第一项
        vi_nj(j,i,:) = vi_nj(j,i,:) + 4.0_RP*pv(:)*vi(j,:)*B_n(i) + 2.0_RP*c2(:)*B_ni_nj(j,i) + C_ni_nj(j,i) !两项之和
        vi_nj(j,i,:) = -vi_nj(j,i,:) / g(:) !除以分母
      end do
    case(3)
      do concurrent( i=1:3, j=1:3)
        vi_nj(j,i,:) = (30.0_RP*c4(:)+12.0_RP*B*c2(:)+2.0_RP*C)*vi(j,:) + 2.0_RP*pv(:)*(2.0_RP*c2(:)*B_n(j)+C_n(j))
        vi_nj(j,i,:) = vi_nj(j,i,:) * vi(i,:) !第一项
        vi_nj(j,i,:) = vi_nj(j,i,:) + 2.0_RP*pv(:)*vi(j,:)*(2.0_RP*c2(:)*B_n(i)+C_n(i)) + c4(:)*B_ni_nj(j,i) + &
          c2(:)*C_ni_nj(j,i) + D_ni_nj(j,i) !两项之和
        vi_nj(j,i,:) = -vi_nj(j,i,:) / g(:) !除以分母
      end do
    end select
  end if
  end subroutine

  !> Calculate phase and ray velocities by MKL (eigenvector method)
  pure subroutine phaseAndRayVelocity_mkl(a, n, pv, rv, vi)
  use lapack95, only: syev, geev
  implicit none
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> Phase velocity direction (normal vector).
  DATATYPE, intent(in):: n(3)
  !> 3 phase velocities.
  DATATYPE, intent(out):: pv(3)
  !> 3 ray velocities.
  DATATYPE, intent(out):: rv(3)
  !> 3 ray velocity vectors.
  DATATYPE, intent(out):: vi(3,3)
  DATATYPE cstl(6), g(3,3), temp(3,3)
  integer i
  !计算christoffel矩阵
  cstl = christoffelMatrix(a, n)
  g(:,1) = cstl(1:3)
  g(2:3,2) = cstl(4:5)
  g(3,3) = cstl(6)
#ifdef COMPLEX
  g(1,2) = g(2,1); g(1:2,3) = g(3,1:2)
  call geev(g,w=pv,vr=temp)
  g = temp
  !排序
  i = maxloc(abs(pv),dim=1)
  if(i/=3) then
    cstl(1) = pv(3); pv(3) = pv(i); pv(i) = cstl(1)
    cstl(1:3) = g(:,3); g(:,3) = g(:,i); g(:,i) = cstl(1:3)
  end if
  i = maxloc(abs(pv(1:2)),dim=1)
  if(i/=2) then
    cstl(1) = pv(2); pv(2) = pv(i); pv(i) = cstl(1)
    cstl(1:3) = g(:,2); g(:,2) = g(:,i); g(:,i) = cstl(1:3)
  end if
  !归一化
  do concurrent (i=1:3)
    !g(:,i) = g(:,i) / sqrt(sum(g(:,i)*g(:,i)))
  end do
#else
  call syev(g, pv, "V", "L")
#endif
  pv = sqrt(pv)
  !计算速度
  do i=1,3
    vi(:,i) = vi_mkl(a, n(:)/pv(i), g(:,i), conjg(g(:,i)))
    rv(i) = sqrt(vi(1,i)*vi(1,i)+vi(2,i)*vi(2,i)+vi(3,i)*vi(3,i))
  end do
  end subroutine

  pure function vi_mkl(a, p, g1, g2) result(res)
  implicit none
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> Phase slowness vector.
  DATATYPE, intent(in):: p(3)
  !> eigen vector
  DATATYPE, intent(in):: g1(3), g2(3)
  !> @return vi.
  DATATYPE res(3), tmp(3,3)
  integer, parameter:: IND(3,3) = reshape([1,6,5,6,2,4,5,4,3],[3,3]) !convert Aijkl into Amn
  integer i, j, k, L, ii, jj
  tmp = 0d0
  do i = 1, 3
    do L = 1, 3
      do j = 1, 3
        do k = 1, 3
          ii = ind(i,j)
          jj = ind(k,L)
          tmp(i,L) = tmp(i,L) + a(ii,jj)*g1(j)*g2(k)
        end do
      end do
    end do
  end do
  !tmp = (tmp + transpose(tmp)) *0.5d0
  res = 0d0
  do i = 1, 3
    res(i) = sum(tmp(i,:)*p(:))
  end do
  end function

  !> Calculate the partial derivatives of phase and ray velocities with respect to two angles \a theta and \a phi.
  pure subroutine velocity_angles(theta, phi, vi_nj, vi, rv, pv_theta, pv_phi, vi_theta, vi_phi, rv_theta, rv_phi)
  implicit none
  real(RP), parameter:: DEG2RAD = acos(-1.0_RP)/180.0_RP !角度转弧度系数
  integer, parameter:: n = 3 !两个S波
  !> Angle \a theta (in rad).
  DATATYPE, intent(in):: theta
  !> Angle \a phi (in rad).
  DATATYPE, intent(in):: phi
  !> Partial derivative: ?vi/?n_j.
  DATATYPE, intent(in):: vi_nj(3,3,n)
  !> Ray velocity vectors.
  DATATYPE, intent(in):: vi(3,n)
  !> Ray velocities.
  DATATYPE, intent(in):: rv(n)
  !> Partial derivative: ?pv/?theta.
  DATATYPE, intent(out):: pv_theta(n)
  !> Partial derivative: ?pv/?phi.
  DATATYPE, intent(out):: pv_phi(n)
  !> Partial derivative: ?vi/?theta.
  DATATYPE, intent(out):: vi_theta(3,n)
  !> Partial derivative: ?vi/?phi.
  DATATYPE, intent(out):: vi_phi(3,n)
  !> Partial derivative: ?rv/?theta.
  DATATYPE, intent(out):: rv_theta(n)
  !> Partial derivative: ?rv/?phi.
  DATATYPE, intent(out):: rv_phi(n)
  DATATYPE sT, sP, cT, cP, n_theta(3), n_phi(3)
  integer i, j
  !计算三角函数
  !sT = theta * DEG2RAD; sP = phi * DEG2RAD
  !cT = cos(sT); cP = cos(sP); sT = sin(sT); sP = sin(sP)
  cT = cos(theta); cP = cos(phi); sT = sin(theta); sP = sin(phi)
  !计算n对theta和phi的偏导
  n_theta = [cT*cP, cT*sP, -sT] * DEG2RAD
  n_phi(1) = -sT*sP; n_phi(2) = sT*cP; n_phi(3) = 0.0_RP; n_phi = n_phi * DEG2RAD
  do concurrent (i=1:n)
    pv_theta(i) = dot_product(vi(:,i), n_theta)
    pv_phi(i)   = dot_product(vi(:,i), n_phi)
    do j = 1, 3
      vi_theta(j,i) = dot_product(vi_nj(:,j,i), n_theta)
      vi_phi(j,i)   = dot_product(vi_nj(:,j,i), n_phi)
    end do
    rv_theta(i) = sum(vi(:,i)*vi_theta(:,i)) / rv(i)
    rv_phi(i)   = sum(vi(:,i)*vi_phi(:,i))   / rv(i)
  end do
  end subroutine

  !> Slpit qS1 and qS2 waves.
  subroutine splitQsWaves(pv, rv, vi, pv_theta, pv_phi, vi_theta, vi_phi, rv_theta, rv_phi)
  implicit none
  integer, parameter:: nw = 2 !两种波
  integer, parameter:: nc = 1 !采样密度
  integer, parameter:: nTheta = 180*nc, nPhi = 360*nc
  real(RP),parameter:: dA = 0.5d0 / nc
  !> Distribution of phase velocities.
  DATATYPE, intent(inout):: pv(nw,0:nTheta,0:nPhi)
  !> Distribution of ray velocities.
  DATATYPE, intent(inout):: rv(nw,0:nTheta,0:nPhi)
  !> Distribution of ray velocity vectors.
  DATATYPE, intent(inout):: vi(3,nw,0:nTheta,0:nPhi)
  !> Partial derivative: ?pv/?theta.
  DATATYPE, intent(inout):: pv_theta(nw,0:nTheta,0:nPhi)
  !> Partial derivative: ?pv/?phi.
  DATATYPE, intent(inout):: pv_phi(nw,0:nTheta,0:nPhi)
  !> Partial derivative: ?vi/?theta.
  DATATYPE, intent(inout):: vi_theta(3,nw,0:nTheta,0:nPhi)
  !> Partial derivative: ?vi/?phi.
  DATATYPE, intent(inout):: vi_phi(3,nw,0:nTheta,0:nPhi)
  !> Partial derivative: ?rv/?theta.
  DATATYPE, intent(inout):: rv_theta(nw,0:nTheta,0:nPhi)
  !> Partial derivative: ?rv/?phi.
  DATATYPE, intent(inout):: rv_phi(nw,0:nTheta,0:nPhi)
  real(RP), parameter:: EPS = 1.0d-5
  DATATYPE pv0(nw), pv1(nw), rv0(nw), rv1(nw), vi0(3,nw), vi1(3,nw)
  logical, allocatable:: changed(:,:)
  integer, allocatable:: imark(:,:)
  real(RP), allocatable:: trust(:,:), error(:,:)
  integer n, loct(2), ix, iy, i
  real(RP) e1, e2

  !分配临时数组
  allocate(changed(0:nTheta,0:nPhi), imark(0:nTheta,0:nPhi))
  allocate(trust(0:nTheta,0:nPhi), error(0:nTheta,0:nPhi))

  changed = .false.
  !差异最大的点, 作为基准点
  imark = 0; imark(:,nPhi) = 3 ! phi=360 不参与计算
  trust = 0.0; error = huge(1.0)
  loct = maxloc(abs(pv(1,:,:)-pv(2,:,:)), mask=(imark/=3)) - 1
  trust(loct(1),loct(2)) = abs(pv(1,loct(1),loct(2))-pv(2,loct(1),loct(2)))
  error(loct(1),loct(2)) = 0.0
  imark(loct(1),loct(2)) = 2 !0未计算, 1不确定, 2确定, 3已使用

  n = count(imark<3)
  do while(n>0)
    n=n-1
    !确定性最大的点
    loct = maxloc(trust,mask=(imark<=2)) - 1
    imark(loct(1),loct(2)) = 3
    !周围需要处理的点 ix, iy
    do i = loct(2)-1, loct(2)+1 !phi
      if(i<0) then
        ix = nPhi-1
      elseif(i==nPhi) then
        ix = 0
      else
        ix = i
      end if
      do iy = loct(1)-1, loct(1)+1 !theta
        if(iy<0.or.iy>nTheta) cycle
        if(imark(iy,ix)==3) cycle !已经确定的点
        !计算误差
        call ERR(iy, ix, e1, e2)
        if(e2<e1) then
          if(e1>trust(iy,ix)) then
            call exchangeData(iy, ix) !交换
            error(iy,ix) = e2; trust(iy,ix) = e1
          end if
        else
          if(e2>trust(iy,ix)) then
            error(iy,ix) = e1; trust(iy,ix) = e2
          end if
        end if
        if(trust(iy,ix)>0.1) then
          imark(iy,ix) = 2 !确定的
        else
          imark(iy,ix) = 1 !不确定
        end if
      end do
    end do
  end do

  !处理边界
  imark(:,nPhi) = imark(:,0)
  trust(:,nPhi) = trust(:,0)
  error(:,nPhi) = error(:,0)

  pv(:,:,nPhi) = pv(:,:,0)
  pv_theta(:,:,nPhi) = pv_theta(:,:,0)
  pv_phi(:,:,nPhi) = pv_phi(:,:,0)
  rv(:,:,nPhi) = rv(:,:,0)
  vi(:,:,:,nPhi) = vi(:,:,:,0)
  vi_theta(:,:,:,nPhi) = vi_theta(:,:,:,0)
  vi_phi(:,:,:,nPhi) = vi_phi(:,:,:,0)
  rv(:,:,nPhi) = rv(:,:,0)
  rv_theta(:,:,nPhi) = rv_theta(:,:,0)

  !print*,count(imark==0),count(imark==1),count(imark==2),count(imark==3)
  write(*,*) '最大误差: ', maxval(error), maxloc(error)
  write(*,*) '最小信度: ', minval(trust)
  write(*,*) '最小T-E: ', minval(trust-error)
  write(*,*) '修改点数: ', count(changed)

  deallocate(changed, imark, trust, error)

  contains
  !> Calculate the difference of phase velocity between point (iTheta, iPhi) and its neighbors.
  !> @details if e1<e2, then the first one is qS1 and the second is qS2;
  !> else the first one is qS2 and the second is qS1.
  pure subroutine ERR(iTheta, iPhi, e1, e2)
  !> Point position.
  integer, intent(in):: iTheta, iPhi
  !> Errors.
  real(RP), intent(out):: e1, e2
  real(RP) e3, e4
  integer i, j, ix, iy, xn, yn, n
  DATATYPE pv0(nw), pv1(nw), p_theta(nw), p_phi(nw)
  !遍历周围的确定点
  n = 0 !数量, 不超过8
  e1 = 0; e2 = 0
  p_theta=0; p_phi=0
  do i = iPhi-1, iPhi+1 !phi
    if(i<0) then
      ix = nPhi-1
    elseif(i==nPhi) then
      ix = 0
    else
      ix = i
    end if
    xn = i - iPhi !梯度方向
    do iy = iTheta-1, iTheta+1 !theta
      if(iy<0 .or. iy>nTheta .or. (ix==iPhi.and.iy==iTheta) ) cycle
      if(imark(iy,ix)<2) cycle !跳过不确定的点
      yn = iy - iTheta
      n = n + 1
      !相速度连续
      pv0 = pv(1:nw,iTheta,iPhi) + pv_phi(1:nw,iTheta,iPhi)*dA*xn + pv_theta(1:nw,iTheta,iPhi)*dA*yn
      pv1 = pv(1:nw,iy,ix) - pv_phi(1:nw,iy,ix)*dA*xn - pv_theta(1:nw,iy,ix)*dA*yn
      e1 = e1 + sum( abs(pv0-pv1)**2 )
      e2 = e2 + sum( abs(pv0-pv1(nw:1:-1))**2 )
      !相速度光滑
      p_theta = p_theta  + pv_theta(1:nw,iy,ix)
      p_phi = p_phi + pv_phi(1:nw,iy,ix)

      !群速度连续
      !pv0 = rv(1:nw,iTheta,iPhi) + rv_phi(1:nw,iTheta,iPhi)*dA*xn + rv_theta(1:nw,iTheta,iPhi)*dA*yn
      !pv1 = rv(1:nw,iy,ix) - rv_phi(1:nw,iy,ix)*dA*xn - rv_theta(1:nw,iy,ix)*dA*yn
      !e1 = e1 + sum( abs(pv0-pv1) )
      !e2 = e2 + sum( abs(pv0-pv1(nw:1:-1)) )
      !速度分量
      !do j = 1, 3
      !  pv0 = vi(j,1:nw,iTheta,iPhi) + vi_phi(j,1:nw,iTheta,iPhi)*dA*xn + vi_theta(j,1:nw,iTheta,iPhi)*dA*yn
      !  pv1 = vi(j,1:nw,iy,ix) - vi_phi(j,1:nw,iy,ix)*dA*xn - vi_theta(j,1:nw,iy,ix)*dA*yn
      !  e1 = e1 + sum( abs(pv0-pv1) )
      !  e2 = e2 + sum( abs(pv0-pv1(nw:1:-1)) )
      !end do
    end do
  end do
  p_theta = p_theta / n
  p_phi = p_phi / n
  e3 = sum( abs(p_theta-pv_theta(1:nw,iTheta,iPhi))**2+abs(p_phi-pv_phi(1:nw,iTheta,iPhi))**2 )
  e4 = sum( abs(p_theta-pv_theta(nw:1:-1,iTheta,iPhi))**2+abs(p_phi-pv_phi(nw:1:-1,iTheta,iPhi))**2 )
  e1 = sqrt( (e1+e3)/(2*n+4) )
  e2 = sqrt( (e2+e4)/(2*n+4) )
  end subroutine

  !> Exchange the data of two qS waves
  subroutine exchangeData(iTheta, iPhi)
  implicit none
  !> Point position.
  integer, intent(in):: iTheta, iPhi
  DATATYPE vec(3), scl
  changed(iTheta,iPhi) = .not. changed(iTheta,iPhi)
  scl = pv(1,iTheta,iPhi); pv(1,iTheta,iPhi) = pv(2,iTheta,iPhi); pv(2,iTheta,iPhi) = scl
  scl = pv_theta(1,iTheta,iPhi); pv_theta(1,iTheta,iPhi) = pv_theta(2,iTheta,iPhi); pv_theta(2,iTheta,iPhi) = scl
  scl = pv_phi(1,iTheta,iPhi); pv_phi(1,iTheta,iPhi) = pv_phi(2,iTheta,iPhi); pv_phi(2,iTheta,iPhi) = scl
  scl = rv(1,iTheta,iPhi); rv(1,iTheta,iPhi) = rv(2,iTheta,iPhi); rv(2,iTheta,iPhi) = scl
  scl = rv_theta(1,iTheta,iPhi); rv_theta(1,iTheta,iPhi) = rv_theta(2,iTheta,iPhi); rv_theta(2,iTheta,iPhi) = scl
  scl = rv_phi(1,iTheta,iPhi); rv_phi(1,iTheta,iPhi) = rv_phi(2,iTheta,iPhi); rv_phi(2,iTheta,iPhi) = scl
  vec = vi(:,1,iTheta,iPhi); vi(:,1,iTheta,iPhi) = vi(:,2,iTheta,iPhi); vi(:,2,iTheta,iPhi) = vec
  vec = vi_theta(:,1,iTheta,iPhi); vi_theta(:,1,iTheta,iPhi) = vi_theta(:,2,iTheta,iPhi); vi_theta(:,2,iTheta,iPhi) = vec
  vec = vi_phi(:,1,iTheta,iPhi); vi_phi(:,1,iTheta,iPhi) = vi_phi(:,2,iTheta,iPhi); vi_phi(:,2,iTheta,iPhi) = vec
  end subroutine exchangeData
  end subroutine

  !> Calculate the phase and ray velocities, and their partial derivatives on specified ray direction.
  pure subroutine solveRayDirection(waveType, a, N, n0, c, v, c_a, err)
  implicit none
  !> Wave type.
  integer, intent(in):: waveType
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> Specified ray direction.
  DATATYPE, intent(in)::  N(3)
  !> in: initial phase direction; out: final phase direction.
  DATATYPE, intent(inout):: n0(3)
  !> Phase velocity.
  DATATYPE, intent(out):: c
  !> Ray velocity.
  DATATYPE, intent(out):: v
  !> Partial derivative: ∂c/∂a.
  DATATYPE, intent(out):: c_a(6,6)
  !> Errors between specified and calculated ray directions.
  real(RP), intent(out):: err
  DATATYPE pv(3), rv(3), vi(3,3), pv_a(6,6,3), vi_nj(3,3,3)
  DATATYPE matA(3,3), vecB(3), vecX(3)
  DATATYPE nray(3)
  integer its, i, j
  nray = N / sqrt(sum(N*N))
  do its = 1, 1000
    call phaseAndRayVelocity(a, n0, pv, rv, vi, pv_a=pv_a, vi_nj=vi_nj)
    ! 保证Q值为正
    !if(v%im>1.0d-10) then
    !  n0(1)%im = -n0(1)%im
    !  n0(2)%im = -n0(2)%im
    !  n0(3)%im = -n0(3)%im
    !  call phaseAndRayVelocity(a, n0, pv, rv, vi, pv_a=pv_a, vi_nj=vi_nj)
    !end  if
    v = rv(waveType)
    c = pv(waveType)
    vecB(:) = nray(:) - vi(:,waveType)/rv(waveType) !误差
    err = norm2(abs(vecB))
    if(err<1.0d-7) exit
    do concurrent (i=1:3, j=1:3)
      matA(i,j) = vi_nj(j,i,waveType)*rv(waveType) - vi(i,waveType)/rv(waveType)*dot_product(vi(:,waveType),vi_nj(j,:,waveType))
    end do
    matA= matA / (rv(waveType)*rv(waveType))
    !求解
    call SIRT(3, matA, vecB, vecX)
    n0 = n0 + vecX
    n0 = n0 / sqrt(sum(n0*n0))
  end do
  c_a = pv_a(:,:,waveType)
  !write(*,"('num of its: ',i0,' error = ', g0)")  its, err
  end subroutine

#ifdef COMPLEX  
  !> Calculate the phase and ray velocities on specified phase direction.
  subroutine solvePhaseDirection(waveType, a, theta, phi, N, n0, c, v, err)
  implicit none
  !> Wave type.
  integer, intent(in):: waveType
  !> Density normalized elastic matrix.
  DATATYPE, intent(in):: a(6,6)
  !> in: Specified real angle; out: complex angel that its real part is equal to the initial value.
  DATATYPE, intent(inout):: theta, phi
  !> Phase direction.
  DATATYPE, intent(out):: n0(3)
  !> Ray direction.
  DATATYPE, intent(out):: N(3)
  !> Phase velocity.
  DATATYPE, intent(out):: c
  !> Ray velocity.
  DATATYPE, intent(out):: v
  !> RMSE of the imaginary part o the complex ray direction.
  real(RP), intent(out):: err
  DATATYPE pv(3), rv(3), vi(3,3), pv_a(6,6,3), rv_a(6,6,3), vi_nj(3,3,3)
  DATATYPE matA(3,3), vecB(3), vecX(3), temp(3,3)
  DATATYPE n_theta(3), n_phi(3), nray(3)
  integer its, i, j
  temp = 0
  n_phi= 0
  do its = 1, 1000
    n0 = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    n_theta = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]
    n_phi(1:2) = [-sin(theta)*sin(phi), sin(theta)*cos(phi)]
    call phaseAndRayVelocity(a, n0, pv, rv, vi, pv_a, rv_a, vi_nj)
    v = rv(waveType)
    c = pv(waveType)
    N = vi(:,waveType)/rv(waveType)
    vecB(:) =  - imag(N) !误差，虚部
    err = norm2(abs(vecB))
    if(err<1.0d-7) exit
    do concurrent (i=1:3, j=1:3) !Ni_nj
      matA(i,j) = vi_nj(j,i,waveType)*rv(waveType) - vi(i,waveType)/rv(waveType)*dot_product(vi(:,waveType),vi_nj(j,:,waveType))
    end do
    matA= matA / (rv(waveType)*rv(waveType))
    ! Ni_angle
    do i=1,3
      temp(i,1) = dot_product(matA(:,i), n_theta) * cmplx(0d0,1d0,RP)
      temp(i,2) = dot_product(matA(:,i), n_phi) * cmplx(0d0,1d0,RP)
    end do
    matA = imag(temp)
    !最小二乘解
    !call LSQR(3, matA, vecB, vecX)
    call SIRT(3, matA, vecB, vecX)
    theta%im = theta%im + real(vecX(1),RP)
    phi%im = phi%im + real(vecX(2),RP)
  end do
  write(*,"('慢度方向迭代次数: ',i0,' error = ', g0)")  its, err
  end subroutine
#endif

  !> Calculate the real sensitive Kernal function
  pure function realKernal_a_q(n, a, kernal, cType) result(res)
  implicit none
  !> Number of independent elastic coefficients.
  integer, intent(in):: n
  !> elastic coefficients array
  complex(RP), intent(in):: a(n)
  !> complex sensitive kernal
  complex(RP), intent(in):: kernal(n)
  !> the type of partial derivative, Velocity/Attenuation to a/Q: V_a, A_a, V_q, A_q
  character(3), intent(in):: cType
  complex(RP), parameter:: IM = cmplx(0.0d0,1.0d0,RP)
  complex(RP) coef(n)
  real(RP) Q(n), res(n)
  where(abs(imag(a(:)))>1.0e-30)
    Q(:) = -real(a(:)) / imag(a(:))
  else where
    Q(:) = 1.0e30
  end where
  select case(cType)
  case('V_a')
    coef(:) = 1.0_RP - IM/Q(:) !偏 a
    res(:) = real(kernal(:) * coef(:))
  case('A_a')
    coef(:) = 1.0_RP - IM/Q(:) !偏 a
    res(:) = imag(kernal(:) * coef(:))
  case('V_q')
    coef(:) = real(a(:)) / (Q(:)*Q(:)) * IM !偏 Q
    res(:) = real(kernal(:) * coef(:))
  case('A_q')
    coef(:) = real(a(:)) / (Q(:)*Q(:)) * IM !偏 Q
    res(:) = imag(kernal(:) * coef(:))
  end select
  end function

  !> Calculate the real sensitive Kernal function
  pure function realKernal_a_i(n, a, kernal, cType) result(res)
  implicit none
  !> Number of independent elastic coefficients.
  integer, intent(in):: n
  !> elastic coefficients array
  complex(RP), intent(in):: a(n)
  !> complex sensitive kernal
  complex(RP), intent(in):: kernal(n)
  !> the type of partial derivative, Velocity/Attenuation to a/Q: V_a, A_a, V_q, A_q
  character(3), intent(in):: cType
  complex(RP), parameter:: IM = cmplx(0.0d0,1.0d0,RP)
  complex(RP) coef(n)
  real(RP) Q(n), res(n)
  where(abs(imag(a(:)))>1.0e-30)
    Q(:) = -real(a(:)) / imag(a(:))
  else where
    Q(:) = 1.0e30
  end where
  select case(cType)
  case('V_a')
    coef(:) = 1.0_RP!偏 a
    res(:) = real(kernal(:) * coef(:))
  case('A_a')
    coef(:) = 1.0_RP!偏 a
    res(:) = imag(kernal(:) * coef(:))
  case('V_q')
    coef(:) = IM !偏 img
    res(:) = real(kernal(:) * coef(:))
  case('A_q')
    coef(:) = IM !偏 img
    res(:) = imag(kernal(:) * coef(:))
  end select
  end function


  !> Solve linear equations Ax=B by SIRT.
  pure subroutine SIRT(n, matA, vecB, vecX)
  !迭代重建算法
  implicit none
  !> Order of equations.
  integer,  intent(in):: n
  !> Matrix A.
  DATATYPE, intent(in):: matA(n,n)
  !> Vector B.
  DATATYPE, intent(in):: vecB(n)
  !> Vector x
  DATATYPE, intent(out)::vecX(n)
  DATATYPE s(n)
  integer j
  s = sum(matA**2,dim=2)!行求和
  where(abs(s)>1.0d-6) s = vecB / s
  do concurrent (j=1:n)
    vecX(j) = sum(s*matA(:,j))
  end do
  vecX = vecX / n
  end subroutine

  end module

