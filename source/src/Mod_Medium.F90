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
  !> @brief Contains module Mod_Medium.
  !

  !> @brief Set the elastic coefficient matrix.
  module Mod_Medium
  implicit none
  !> Convert real/complex vector (\a x, \a y, \a z) to real/complex angles (\a theta, \a phi).
  !> - @ref mod_medium::normal2angle_real "normal2Angle_real" for real data
  !> - @ref mod_medium::normal2angle_complex "normal2Angle_complex" for complex data
  interface normal2Angle
  module procedure:: normal2Angle_real
  module procedure:: normal2Angle_complex
  end interface
  contains

  !> Calculate the attenuation angle of a complex direction vector.
  subroutine attenuationAngle(n, D)
  implicit none
  !> Complex direction vector.
  complex(RP), intent(in)::n(3)
  !> Attenuation angle.
  real(RP), intent(out):: D
  real(RP), parameter:: EPS = 1.0d-7
  complex(RP) c
  real(8) re(3), im(3), a, b
  re(1) = n(1)%re; im(1) = n(1)%im
  re(2) = n(2)%re; im(2) = n(2)%im
  re(3) = n(3)%re; im(3) = n(3)%im
  !D = atan2d(norm2(im),norm2(re))
  !return
  a = norm2(re); if(abs(a)>EPS) re = re / a
  b = dot_product(re,im)
  c = cmplx(a,b,RP) !均匀项系数
  im = im - b*re !非均匀部分
  a = norm2(im)
  !非均匀角度
  D = atan2d(a,b)
  end subroutine

  !> Convert real vector (\a x, \a y, \a z) to real angles (\a theta, \a phi).
  elemental subroutine normal2Angle_real(x,y,z,theta,phi)
  implicit none
  !> x component of the vector.
  real(RP), intent(in):: x
  !> y component of the vector.
  real(RP), intent(in):: y
  !> z component of the vector.
  real(RP), intent(in):: z
  !> Value range is [0, 180].
  real(RP), intent(out):: theta
  !> Value range is [0, 360].
  real(RP), intent(out):: phi
  theta = atan2d(sqrt(x*x+y*y),z)
  theta = abs(theta)
  phi   = atan2d(y,x)
  if(phi<0.0_RP) phi = 360.0_RP + phi
  end subroutine

  !> Convert complex vector (\a x, \a y, \a z) to complex angles (\a theta, \a phi).
  elemental subroutine normal2Angle_complex(x,y,z,theta,phi)
  implicit none
  !> x component of the vector.
  complex(RP), intent(in):: x
  !> y component of the vector.
  complex(RP), intent(in):: y
  !> z component of the vector.
  complex(RP), intent(in):: z
  !> Value range is [0, 180].
  complex(RP), intent(out):: theta
  !> Value range is [0, 360].
  complex(RP), intent(out):: phi
  complex(RP), parameter:: im = cmplx(0.0_RP,1.0_RP,RP)
  real(RP), parameter:: EPS = 1.0d-7, RAD2DEG = 180d0 / acos(-1d0)
  complex(RP) cz
  if(abs(z)<EPS) then
    theta = 90.0d0
  else !atan2d(H,z)
    cz = sqrt(x*x+y*y) / z * im
    cz = (1d0+cz)/(1d0-cz)
    theta = -0.5d0*im*log(cz) * RAD2DEG
    if(theta%re<0) theta = -theta
  end if

  if(abs(x)<EPS) then
    if(y%re>0) then
      phi = 90.0d0
    else
      phi = 270d0
    end if
  else !atan2d(y,x)
    cz = y/x * im
    cz = (1d0+cz)/(1d0-cz)
    phi = -0.5d0*im*log(cz) * RAD2DEG
    if(phi%re<0) phi = 360d0 - phi
  end if
  end subroutine

  !> @brief Convert the elastic coefficients array to a 6×6 matrix.
  pure function mediumArray2Matrix(n, a) result(res)
  implicit none
  !> Number of independent elastic coefficients.
  integer, intent(in):: n
  !> Elastic coefficients array.
  !> @verbatim
  !> isotropic:    2 elements - c11, c44
  !> VTI:          5 elements - c11, c13, c33, c44, c66
  !> Orthorhombic: 9 elements - c11, c12, c13, c22, c23, c33, c44, c55, c66
  !> Monoclinic:  13 elements - c11, c12, c13, c16, c22, c23, c26, c33, c36, c44, c45, c55, c66
  !> Triclinic:   21 elements - upper triangular matrix @endverbatim
  DATATYPE, intent(in):: a(n)
  !> @return 6×6 matrix
  DATATYPE res(6,6)
  integer i, j
  res = 0.0_RP
  !判断介质类型, 计算上三角元素
  select case (n)
  case (2) !isotropic
    res(1,1) = a(1); res(2,2) = a(1); res(3,3) = a(1); res(4,4) = a(2); res(5,5) = a(2); res(6,6) = a(2)
    res(1,2) = a(1) - 2.0*a(2); res(1,3) = res(1,2); res(2,3) = res(1,2)
  case (5) !VTI
    res(1,1) = a(1); res(1,3) = a(2); res(3,3) = a(3); res(4,4) = a(4); res(6,6) = a(5)
    !非独立参数
    res(1,2)=res(1,1)-2*res(6,6); res(2,2)=res(1,1); res(2,3)=res(1,3); res(5,5)=res(4,4)
    !res(3,1)=res(1,3); res(3,2)=res(1,3); res(2,1)=res(1,2)
  case (9) !orthorhombic
    res(1,1) = a(1); res(1,2) = a(2); res(1,3) = a(3); res(2,2) = a(4); res(2,3) = a(5)
    res(3,3) = a(6); res(4,4) = a(7); res(5,5) = a(8); res(6,6) = a(9)
  case (13) !Monoclinic, 对称面XOY
    res(1,1:3) = a(1:3); res(1,6) = a(4); res(2,2:3) = a(5:6); res(2,6) = a(7); res(3,3) = a(8)
    res(3,6) = a(9); res(4,4:5) = a(10:11); res(5,5) = a(12); res(6,6) = a(13)
  case (21) !Triclinic
    res(1,1:6) = a(1:6); res(2,2:6) = a(7:11); res(3,3:6) = a(12:15)
    res(4,4:6) = a(16:18); res(5,5:6) = a(19:20); res(6,6) = a(21)
  end select
  !对称性
  do i = 1, 6
    do j = i+1, 6
      res(j,i) = res(i,j)
    end do
  end do
  end function

  !> @brief Convert the 6×6 kernal matrix to a array.
  pure function kernalMatrix2Array(n, a) result(res)
  implicit none
  !> Number of independent elastic coefficients.
  integer, intent(in):: n
  !> matrix.
  !> @verbatim
  !> isotropic:    2 elements - c11, c44
  !> VTI:          5 elements - c11, c13, c33, c44, c66
  !> Orthorhombic: 9 elements - c11, c12, c13, c22, c23, c33, c44, c55, c66
  !> Monoclinic:  13 elements - c11, c12, c13, c16, c22, c23, c26, c33, c36, c44, c45, c55, c66
  !> Triclinic:   21 elements - upper triangular matrix @endverbatim
  DATATYPE, intent(in):: a(6,6)
  !> @return kernal array
  DATATYPE res(n)
  res(:) = 0.0_RP
  !判断介质类型, 计算上三角元素
  select case (n)
  case (2) !isotropic
    res(2) = a(2,1) + a(3,1) + a(3,2)
    res(1) = a(1,1) + a(2,2) + a(3,3) + res(2)
    res(2) = a(4,4) + a(5,5) + a(6,6) - 2*res(2)
  case (5) !VTI
    res(1) = a(1,1) + a(1,2) + a(2,2)
    res(2) = a(1,3) + a(2,3)
    res(3) = a(3,3)
    res(4) = a(4,4) + a(5,5)
    res(5) = a(6,6) - 2*a(1,2)
  case (9) !orthorhombic
    res(1:3) = a(1:3,1)
    res(4:5) = a(2:3,2)
    res(6) = a(3,3)
    res(7) = a(4,4)
    res(8) = a(5,5)
    res(9) = a(6,6)
  case (13) !Monoclinic, 对称面XOY
    res(:) = [a(1:3,1), a(6,1), a(2:3,2), a(6,2), a(3,3), a(6,3), a(4:5,4), a(5,5), a(6,6)]
  case (21) !Triclinic
    res(1:6)   = a(1:6,1)
    res(7:11)  = a(2:6,2)
    res(12:15) = a(3:6,3)
    res(16:18) = a(4:6,4)
    res(19:20) = a(5:6,5)
    res(21)    = a(6,6)
  end select
  end function
  end module


