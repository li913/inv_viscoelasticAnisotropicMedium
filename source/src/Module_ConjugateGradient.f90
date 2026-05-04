module ConjugateGradient
implicit none
private
public CG_RP
contains

pure subroutine CG_RP( damp, csr_values, csr_columns, csr_rowIndex, dt, dm)
implicit none
real(4),intent(in):: damp, csr_values(:)
integer,intent(in):: csr_columns(:), csr_rowIndex(:)
real(8),intent(in):: dt(:)
real(8),intent(out):: dm(:)
real rms_T, rms_M
real, allocatable:: solu(:), deltaT(:)

allocate(deltaT(size(dt)), solu(size(dm)))
deltaT = dt
call CG( csr_values, csr_columns, csr_rowIndex, deltaT, damp, rms_T, rms_M, solu)
dm = solu 
end subroutine

!--------------------------------------------------------------------
!      the conjugate gradient method is developed to solve the
!      minimum norm,least-squares and constrained problem in
!      non-linear seismic tomography

!      csr_values(:)   CSR 格式雅克比矩阵
!      csr_columns(:)
!      csr_rowIndex(:)
!      deltaT(:)       走时差, obs-calTime
!      damp            damping factor，阻尼因子
!      rms_T           均方根走时误差
!      rms_M           mean square root of the data-fit
!      solu(:)         速度增量dv
!--------------------------------------------------------------------
pure subroutine CG( csr_values, csr_columns, csr_rowIndex, deltaT, damp, rms_T, rms_M, solu)
implicit none
real,   intent(in):: csr_values(:), deltaT(:), damp
integer,intent(in):: csr_columns(:), csr_rowIndex(:)
real,  intent(out):: rms_T, rms_M, solu(:)
real,    parameter:: EPS = 1.0e-6
integer, allocatable:: icol0(:), icol(:), ictn(:), nPerRow(:)
real,allocatable:: q0(:), p0(:), r0(:)
real r00, sum0, alpha, r11
integer nEqs, nUnk, nEle
integer its, i, j, L

nUnk = size(solu)      !未知元素个数
nEqs = size(deltaT)    !方程数, 即射线数
nEle = size(csr_values)!非零值个数

!走时均方根误差
rms_T = 0.0; rms_M = 0.0; solu = 0.0
if(nEqs==0 .or. nUnk==0 .or. nEle==0) then
  return
else
  rms_T = sqrt( sum( deltaT*deltaT ) / size(deltaT) )
  if(rms_T<EPS) return
end if

! CSC 格式中每列非零值数、新序列中元素行号、在原序列中的编号
allocate( ictn(nUnk), icol0(nEle), icol(nEle) )
!每行非零元素个数
allocate(nPerRow(nEqs))

!此处储存格式稍异于CSR
do i = nEqs, 1, -1
  nPerRow(i) = csr_rowIndex(i+1) - csr_rowIndex(i)
end do
!call trans(nEqs,nUnk,csr_columns,nPerRow,icol0,icol,ictn)
!  将行存储(射线)的Jacob矩阵转为列存储(未知元素)
call trans_li(size(csr_values),nEqs,nUnk,csr_columns,nPerRow,icol0,icol,ictn)

allocate(q0(1:nEqs), r0(1:nUnk), p0(1:nUnk))
q0 = deltaT
L = 0
do i = 1, nUnk
  sum0 = 0.0
  do j = 1, ictn(i) !每列非零值个数
    L = L + 1
    sum0 = sum0 + csr_values(icol(L)) * q0(icol0(L))
  end do
  r0(i) = sum0
  p0(i) = sum0
end do

do its = 1, 10000
  L=0
  do i = 1, nEqs
    sum0 = 0.0
    do j = 1, nPerRow(i)
      L = L + 1
      sum0 = sum0 + csr_values(L) * p0(csr_columns(L))
    end do
    q0(i) = sum0
  end do
  r00 = sum( r0*r0 )  !平方和(1:nUnk)
  alpha = r00 / (damp*sum(p0*p0)+sum(q0*q0))
  solu = solu + alpha * p0  !(1:nUnk)
  L = 0
  r11=0.0
  do i = 1 , nUnk
    sum0 = 0.0
    do j = 1, ictn(i)
      L=L+1
      sum0 = sum0 + csr_values(icol(L))*q0(icol0(L))
    end do
    r0(i) = r0(i) - alpha*(damp*p0(i)+sum0)
    r11 = r11 + r0(i)*r0(i)
  end do
  rms_M = sqrt(r11/real(nUnk)) !更新量
  if(rms_M<=EPS) exit !!更新量过小，退出迭代
  p0 = r0 + ( r11/r00 ) * p0  !(1:nUnk)
end do
rms_M = sqrt( sum( solu * solu ) / real(nUnk) ) !模型更改量
deallocate(q0,p0,r0,icol0,icol,ictn,nPerRow)
end subroutine 

!--------------------------------------------------------------------c
!	将行存储(射线)的Jacob矩阵转为列存储(未知元素)
!	nelement......非零值个数
!	nEqs......射线数
!	nUnk......未知元素个数
!	csr_columns......非零值对应的列号(未知元素编号)
!	nPerRow......每条射线中非零值个数
!	icol0......Jacob元素重新排列后各元素对应的射线编号
!	icol......Jacob元素重新排列后各元素在原来序列中的编号
!	ictn......每个未知元素的非零值个数
!--------------------------------------------------------------------c
pure subroutine trans_li(nelement,nEqs,nUnk,csr_columns,nPerRow,icol0,icol,ictn)
implicit none
integer,intent(in):: nelement, nEqs, nUnk
integer,intent(in):: csr_columns(nelement),nPerRow(nEqs)
integer,intent(out)::icol0(nelement),icol(nelement),ictn(nUnk)
integer,allocatable::iloc(:), iray(:)
integer i, m, n
!  初始化变量
ictn = 0; icol = 0; icol0 = 0

!  每个未知元素的非零值个数ictn
do i = 1, nelement
  m=csr_columns(i)  !列号
  ictn(m) = ictn(m) + 1
end do

!  每个未知元素的第一个非零值在新序列中的位置iloc
allocate( iloc(nUnk) )
m = 0; iloc = 1
do i = 2, nUnk
  m = m + ictn(i-1)
  iloc(i) = m + 1
end do

!  原序列中每个元素的行号(属于第几条射线)iray
allocate( iray(nelement) )
n = 0
do i = 1, nEqs
  m = nPerRow(i)
  n = n + m
  iray(n-m+1:n) = i
end do

!  Jacob元素重新排列后各元素在原来序列中的编号icol
do i = 1, nelement
  m=csr_columns(i)  !列号
  n = iloc(m) !新序列中的位置
  icol(n) = i
  icol0(n) = iray(i)
  iloc(m) = iloc(m) + 1
end do

deallocate( iloc, iray )

end subroutine trans_li

!-----------------------------------------------------------------c
!     subroutine 'trans'is to transpose the a-matrix              c
!-----------------------------------------------------------------c
subroutine trans(nEqs,nUnk,csr_columns,nPerRow,icol0,icol,ictn)
implicit none
integer nEqs,nUnk
integer csr_columns(*),nPerRow(*),icol0(*),icol(*),ictn(*)
integer i, j, k, m
integer nt,nelea,ir1,ic

nt = 0
do k = 1, nUnk
  m = 0
  nelea = 0
  do i = 1, nEqs
    ir1 = nPerRow(i)
    do j = 1, ir1
      nelea = nelea + 1
      ic = csr_columns(nelea)
      if(k/=ic) cycle
      m = m+1
      nt = nt+1
      icol0(nt) = i
      icol(nt) = nelea
    end do
  end do
  ictn(k) = m
end do
end subroutine trans
end module 
