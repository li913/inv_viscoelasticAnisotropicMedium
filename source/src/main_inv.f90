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

  subroutine get_input_file(fln)
  implicit none
  character(*), intent(out):: fln
  integer i
  i = COMMAND_ARGUMENT_COUNT()
  write(*,"(A)", advance='no') 'input file: '
  if(i>0) then
    call GET_COMMAND_ARGUMENT(1,fln)
    write(*,"(A)") trim(fln)
  else
    read(*,"(A)") fln
  end if
  end subroutine

  program main_inv
  use Mod_Velocity4GeneralAnisotropicMedium
  use ConjugateGradient
  use Mod_Medium
  use lapack95
  implicit none
  real(RP), parameter:: DEG2RAD = acos(-1.0_RP)/180.0_RP !角度转弧度系数
  complex(RP), parameter:: IM = cmplx(0.0d0,1.0d0,8)
  !基础参数
  namelist /basepar/ radius, nPar, nRay, waveType, invPar, calDat, maxIts
  real(RP) radius   !球形样品半径
  integer  nPar     !独立弹性系数个数
  integer  nRay     !射线数量
  integer  waveType !使用的波类型: 1: qP, 2: qSV/qS1, 3: qSH/qS2
  integer  calDat   !计算哪些走时数据 1: real time, 2: image time; 3: re+im time
  integer  invPar   !需要反演哪些参数 1: a, 2: q, 3: a+q; -1: a, -2: img, -3: a+img;
  integer  maxIts   !最大迭代次数, 0表示仅正演

  !初始/参考模型, 参考(真实)模型仅用于计算反演误差,可忽略
  namelist /mediumPar/ E, Qi, refE, refQi
  real(RP), allocatable:: E(:), Qi(:), refE(:), refQi(:)
  !射线角
  namelist /raydirection/ rayAngle
  real(RP), allocatable:: rayAngle(:,:)
  !实走时
  namelist /realTime/ reTime
  real(RP), allocatable:: reTime(:)
  !虚走时
  namelist /imaginaryTime/ imTime
  real(RP), allocatable:: imTime(:)

  !inv
  real(RP), allocatable:: traveltime(:), obstime(:), timeErr(:)
  real(RP), allocatable:: jacob(:,:)
  real(RP) aveT, rmsT, rmsE, rmsQi
  real(RP), allocatable:: dm(:)
  real(4), allocatable:: csr_value(:)
  integer, allocatable:: csr_colnum(:), csr_index(:)

  !正演
  DATATYPE, allocatable:: aa(:), kernal(:)
  DATATYPE a(6,6), theta, phi, vi(3), n(3), rv, pv, pv_a(6,6), cTime
  real(RP), allocatable:: rayErr(:)

  integer row1, row2, col1, col2, nRow, nCol
  integer its, nEle, i
  real(RP) limit, lit1, lit2
  character(512):: fln = "GA_inp_0.1_0.1.txt"

  !读取文件
  call get_input_file(fln)
  write(*,*)

  open(10, file=fln,status='old')
  read(10,nml=basepar)
  close(10)

  allocate(E(nPar), Qi(nPar), refE(nPar), refQi(nPar))
  E = 0; Qi = 0; refE = -999; refQi = -999
  open(10, file=fln,status='old')
  read(10,nml=mediumpar)
  close(10)

  allocate(rayAngle(2,nRay))
  open(10, file=fln,status='old')
  read(10,nml=raydirection)
  close(10)

  if(maxIts>0) then
    !观测走时
    allocate( obstime(2*nRay), source=0.0_RP )
    if(calDat==1 .or. calDat==3) then
      allocate(reTime(nRay))
      open(10, file=fln,status='old')
      read(10,nml=realTime)
      close(10)
      obstime(1:nRay) = reTime(:)
      deallocate(reTime)
    end if
    if(calDat==2 .or. calDat==3) then
      allocate(imTime(nRay))
      open(10, file=fln,status='old')
      read(10,nml=imaginaryTime)
      close(10)
      obstime(nRay+1:) = imTime(:)
      deallocate(imTime)
    end if
  end if
  !文件读取完毕

  !实际数据量
  if(calDat==1) then
    row1 = 1; row2 = nRay
  else if(calDat==2) then
    row1 = nRay + 1; row2 = 2 * nRay
  else if(calDat==3) then
    row1 = 1; row2 = 2 * nRay
  end if
  nRow = row2 - row1 + 1

  if(abs(invPar)==1) then
    col1 = 1; col2 = nPar
  else if(abs(invPar)==2) then
    col1 = nPar + 1; col2 = 2 * nPar
  else if(abs(invPar)==3) then
    col1 = 1; col2 = 2 * nPar
  end if
  nCol = col2 - col1 + 1

  !用哪种核函数
  if(invPar>0) then
    realKernal => realKernal_a_q
  else
    realKernal => realKernal_a_i
  end if

  !分配空间
  allocate(aa(nPar), kernal(nPar))
  allocate(traveltime(2*nRay), rayErr(nRay))
  if(maxIts>0) then
    !走时相关
    allocate(timeErr(2*nRay), source=0.0_RP)
    allocate(jacob(2*nRay,2*nPar), dm(2*nPar))
    allocate(csr_value(nRow*nCol), csr_index(nRow+1), csr_colnum(nRow*nCol))
  end if

  open(10, file='traveltime.txt')
  if(maxIts>0) then
    open(11, file='rms.csv')
    write( *,"(A)")  'its        aveT          rmsT          rmsE          rmsQi'
    write(11,"(A)")  'its, aveT, rmsT, rmsE, rmsQi'
    open(12, file='result.csv')
  end if

  do its = 1, max(1,maxIts)
    !黏弹模型
    aa(:)%re = E(:)
    if(invPar>0) then
      aa(:)%im = -E(:) / Qi(:)
    else
      aa(:)%im = Qi(:)
    end if
    a = mediumArray2Matrix(nPar, aa)
    !计算走时及敏感核
    do i = 1, nRay
      theta = rayAngle(1,i) * DEG2RAD
      phi   = rayAngle(2,i) * DEG2RAD
      vi(:) = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
      n = vi(:) !初始相速度方向为射线方向
      call solveRayDirection(waveType, a, vi, n, pv, rv, pv_a(:,:), rayErr(i))
      !走时
      cTime = 2D0 * radius / rv
      traveltime(i) = cTime%re
      traveltime(nRay+i) = cTime%im
      !计算jacob
      if(maxIts>0) then
        kernal(:) = 2D0 * radius * kernalMatrix2Array(nPar,pv_a(:,:)) / (-pv*rv)
        jacob(i,:)  = [realKernal(nPar, aa, kernal, "V_a"), realKernal(nPar, aa, kernal, "V_q")]
        jacob(nRay+i,:) = [realKernal(nPar, aa, kernal, "A_a"), realKernal(nPar, aa, kernal, "A_q")]
      end if
    end do
    !正演完毕, 输出信息
    write(10,"(A,g0)") 'its=', its
    write(10,"(g0)") traveltime(row1:row2)

    if(maxIts<1) exit

    !误差
    timeErr(row1:row2) =  obstime(row1:row2) - traveltime(row1:row2)
    aveT = sum(abs(timeErr(row1:row2))) / nRow
    rmsT = norm2(timeErr) / sqrt(real(nRow))
    if(nPar==5) then !VTI介质特殊性
      if(waveType==qSH) then
        rmsE  = norm2(E(4:5)-refE(4:5))   / sqrt(2.0d0)
        rmsQi = norm2(Qi(4:5)-refQi(4:5)) / sqrt(2.0d0)
      else
        rmsE  = norm2(E(1:4)-refE(1:4))  / 2
        rmsQi = norm2(Qi(1:4)-refQi(1:4)) / 2
      end if
    else
      rmsE  = norm2(E(:)-refE(:))   / sqrt(real(nPar))
      rmsQi = norm2(Qi(:)-refQi(:)) / sqrt(real(nPar))
    end if
    if(refE(1)<-990)  rmsE  = -99 !无参考模型
    if(refQi(1)<-990) rmsQi = -99
    write( *,"(i3,1x,*(f14.7))")  its, aveT, rmsT, rmsE, rmsQi
    write(11,"(*(g0,:,','))")     its, aveT, rmsT, rmsE, rmsQi

    if(rmsT<1.0e-7) exit

    !jacob(:,1:nPar) = jacob(:,1:nPar) / maxval(abs(jacob(:,1:nPar)))
    !jacob(:,nPar+1:) = jacob(:,nPar+1:) / maxval(abs(jacob(:,nPar+1:)))

    !反演方法选择
    dm(:) = 0
    select case(3)
    case(1)
      !稀疏存储
      call csr(nRow, nCol, jacob(row1:row2,col1:col2), csr_value, csr_colnum, csr_index, nEle)
      call CG_RP( 0.1, csr_value(1:nEle), csr_colnum(1:nEle), csr_index, timeErr(row1:row2), dm(col1:col2))
    case(2)
      call SIRT(nRow, nCol, jacob(row1:row2,col1:col2), timeErr(row1:row2), dm(col1:col2))
    case(3)
      call gelss(jacob(row1:row2,col1:col2), timeErr(row1:row2))
      dm(col1:col2) = timeErr(row1:row1+nCol-1)
    end select

    !VTI介质特殊性
    if(nPar==5) then
      if(waveType==qSH) then
        dm(1:3) = 0
        dm(6:8) = 0
      else
        dm(5) = 0
        dm(10) = 0
      end if
    end if

    !约束结果
    !if(its<=5) then
    !  lit1 = 1
    !  lit2 = 0.01
    !else if(its<=10) then
    !  lit1 = 0.5
    !  lit2 = 0.005
    !else if(its<=15) then
    !  lit1 = 0.1
    !  lit2 = 0.001
    !else
    !  lit1 = 0.05
    !  lit2 = 0.0005
    !end if
    !limit = maxval( abs(dm(1:nPar)) )
    !if(lit1<limit) dm(1:nPar) = dm(1:nPar) / limit * lit1
    !limit = maxval( abs(dm(nPar+1:)) )
    !if(invPar==-2.or.invPar==-3) then
    !  if(lit2<limit) dm(nPar+1:) = dm(nPar+1:) / limit * lit2
    !else
    !  if(lit2<limit) dm(nPar+1:) = dm(nPar+1:) / limit * lit1
    !end if

    if(its<=20) then
      limit = 1
    else if(its<=40) then
      limit = 0.5
    else if(its<=60) then
      limit = 0.25
    else
      limit = 0.02
    end if
    limit = cosd(90.0*its/150.0)
    !limit = max(0.05,limit)
    if(maxval(abs(dm))>limit) dm = limit * dm / maxval(abs(dm))



    !更新
    select case( abs(invPar) )
    case(1)
      E(:) = E(:) + dm(1:nPar)
    case(2)
      Qi(:) = Qi(:) + dm(nPar+1:)
    case(3)
      E(:) = E(:) + dm(1:nPar)
      Qi(:) = Qi(:) + dm(nPar+1:)
    end select

    write(12,"(a,g0)") 'its=', its
    write(12,"(*(g0,:,','))")  E(:)
    write(12,"(*(g0,:,','))")  Qi(:)
  end do

  !输出提示
  if(maxIts>0) then
    write(*,*)
    write(*,"('result:')")
    write(*,"('E: ', *(f7.3,2x))")  E(:)
    if(invPar>0) then
      write(*,"('Q: ', *(f7.3,2x))")  Qi(:)
      Qi(:) = -E(:) / Qi(:)
      write(*,"('i: ', *(f7.3,2x))")  Qi(:)
    else
      write(*,"('i: ', *(f7.3,2x))")  Qi(:)
      where(abs(Qi)>1.0e-10)
        Qi(:) = -E(:) / Qi(:)
      else where
        Qi(:) = 999
      end where
      write(*,"('Q: ', *(f7.3,2x))")  Qi(:)
    end if
  else
    write(*,*) 'forward finished.'
  end if
  read(*,*)
  contains

  subroutine csr(nRow, nCol, array, A, colnum, index, nEle)
  implicit none
  integer, intent(in):: nRow, nCol
  real(RP), intent(in):: array(nRow,nCol)
  real(4), intent(out):: A(:)
  integer, intent(out):: colnum(:), index(nRow+1), nEle
  real(RP) v
  integer i, j, n
  A = 0
  colnum = 0
  index(1) = 1
  nEle = 0
  do j = 1, nRow
    n = 0
    do i = 1, nCol
      v = array(j,i)
      if(abs(v)<1.0d-7) cycle
      nEle = nEle + 1
      n = n + 1
      A(nEle) = v
      colnum(nEle) = i
    end do
    index(j+1) = index(j) + n
  end do
  end subroutine

  !//************************************************************************
  !//计算规则点的速度，近似提取最近点的值
  ! 输入:
  ! theta, phi    ... 规则点的位置
  ! vArray(3,:,:) ... 所有非规则点的速度分量
  ! 输出:
  ! v             ... 规则点的速度
  !//************************************************************************
  !attributes(device) &
  pure subroutine get_regularVelocityArray(theta, phi, vArray, v)
  implicit none
  real(RP),      value:: theta, phi
  real(RP), intent(in):: vArray(:,:,:)
  real(RP),intent(out):: v
  real(RP) maxDot, x, y, z, dot
  real(RP), parameter:: deg = 10.0 !搜索范围
  real(RP), parameter:: dTheta0 = 1.0, dPhi0 = 1.0, deg2Rad = acos(-1.0)/180.0
  real(RP) v0
  integer nD, nU, mD, mU, i, j, iP
  !计算theta的循环范围（扩大10度）, 缩减计算量
  nD = (theta-deg) / dTheta0; nD = max(1,nD)
  nU = (theta+deg) / dTheta0; nU = min(size(vArray,2),nU)
  !phi的搜索范围
  mD = 1; mU = size(vArray,3)
  if(theta>2.0*deg.AND.theta<180.0-2.0*deg) then
    mD = (phi-2.0*deg) / dPhi0
    mu = (phi+2.0*deg) / dPhi0
  end if

  !角度转为弧度
  theta = theta * deg2Rad
  phi   = phi   * deg2rad
  !角度转为向量
  z = sin(theta)
  x = cos(phi) * z
  y = sin(phi) * z
  z = cos(theta)
  !查找最近的速度值
  maxDot = -1.0e38
  !do j = 1, size(vArray,2)
  do j = nD, nU
    do i = mD, mU !phi
      !iP = i - size(vArray,3)
      !if(iP<1) iP = iP + size(vArray,3)
      ip = modulo(i,size(vArray,3)) + 1
      v0 = norm2( vArray(:,j,ip) )
      dot = (vArray(1,j,iP)*x + vArray(2,j,iP)*y + vArray(3,j,iP)*z)/v0 !点乘dot值越大夹角越小
      if(dot>maxDot) then
        maxDot = dot; v = v0
      end if
    end do
  end do
  end subroutine

  pure subroutine SIRT(m, n, matA, vecB, vecX)
  !迭代重建算法
  implicit none
  integer,  intent(in):: m, n
  real(RP), intent(in):: matA(m,n), vecB(m)
  real(RP), intent(out)::vecX(n)
  real(RP) s(m)
  integer j
  s = sum(matA**2,dim=2)!行求和
  where(abs(s)>1.0d-12) s = vecB / s
  do concurrent (j=1:n)
    vecX(j) = sum(s*matA(:,j))
  end do
  vecX = vecX / m
  end subroutine

  end program

