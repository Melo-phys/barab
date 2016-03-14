module globals
  
  !n -> number of nodes
  !ne -> number of linked nodes
  !ktot -> number of links (ones on the conections matrix)
  !ll -> loops counter
  !flag -> flag
  !pseed -> pseudo seed
  !m -> conections matrix
  !v -> neighbors array
  !k -> conectivity array (number of neighbors) 
  !o -> origin  matrix (needs to be updated)
  !b -> betweenes array
  !dis -> distances matrix
  !nc -> number of minimum paths matrix
  !cod -> program code for the origins matrix
  !cl -> clustering array
  !del -> deleted nodes
  !nodes -> existing nodes
  
  integer,parameter :: nl=1000000,nn=10000
  integer :: n,ne,ktot,ll,pseed,cod,jj,flag
  integer,allocatable :: v(:),k(:),o(:,:),del(:),p2(:),p2k(:),tt(:),ttk(:),pk(:)
  real*8,allocatable :: cl(:),clhist(:,:),lk(:,:),dkk(:,:,:)
  real*8 :: frac,ddd,dddd
  real*8,parameter ::a=1, b=0, c=0
  
  integer,allocatable :: kold(:),oold(:,:)
  real*8,allocatable :: clold(:)
  
end module globals

program evol
  !*********************************************************!
  !                            TWANG                        !
  !*********************************************************!
  !
  !
  !Dupliation follows (1-cl)*k
  !Horizontal transfer now is forced.
  !It has one variable to control the horizontal 
  !duplication ratio. Keeps the origins of each gene. 
  !Mutation Added
  !Old way to pick nodes
  !Clhist and Dist added      ***WORK***
  ! 
  !VARIABLES:
  !nn -> number of nodes of the final network
  !i,j,in -> loop indexes
  !pk -> conectivity histogram 
  !z -> random variable
  !hd -> horizontal-duplication ratio
  

  
  use globals
  implicit none
  integer :: i,j,l,in,no,kmax,soma,sk,sk2,sk3,no1,no2,numi,numf,num
  integer,allocatable :: P(:,:),lll(:),d(:,:,:),t(:,:,:),nk(:)
  real*8 :: z,hd,km,km2,km3
  character(14) :: clfile
  character(14) :: pkfile
  character(14) :: lfile
  character(8) :: ddfile
  character(5) :: aaaa
  character(5) :: numberi
  character(5) :: numberf
  character(5) :: number
  character(5) :: randnum
!  open(1,file='pk.dat')
  open(2,file='comoliga.dat')
  open(3,file='nlinks.dat')
  open(4,file='del.dat')
!  open(55,file='cl.dat')
!  open(6,file='arq.dat')
!  read(6,*)frac, hd
  hd=1
!  open(8,file='dist.dat')
  !frac=0.5
  !hd=0.9d0
  call getarg(1,numberi)
  read(numberi,'(I5)')numi
  call getarg(2,numberf)
  read(numberf,'(I5)')numf
  call getarg(3,randnum)
  read(randnum,'(I5)')pseed
    
  do num=numi,numf
    number=''
    write(number,'(I5.5)')num
    call inicia
    ddfile="dd_"//number
    open(123123,file=ddfile)
    ddd=0
    dddd=0
    allocate(tt(nl),p2(nl))
    
    do in=2,4
      clfile(1:9)='cl_'//number//"_"
      pkfile(1:9)='pk_'//number//"_"
      lfile(1:9)='ll_'//number//"_"
!    open(111,file='intermednome')
!    write(111,'(I5)')in
      aaaa=''
      write(aaaa,'(I5.5)')n
!    read(111,'(A)')aaaa
      close(1111)
      clfile(10:14)=adjustl(aaaa)
      pkfile(10:14)=adjustl(aaaa)
      lfile(10:14)=adjustl(aaaa)
      write(103,*)clfile,' ',pkfile,' ',aaaa
!    open(101,file=clfile)
!    open(102,file=pkfile)
!    open(104,file=lfile)
      
!    do while((n<in*(nn/10)).and.(ktot.lt.nl))
      do while((n<(10**in)).and.(ktot.lt.nl))
        ll=ll+1
        call random_number(z)
        if(z.le.hd)then        !chooses horizontal transfer or dupication
          flag=0
          do while (flag.ne.1)
            call horizontal
          end do
        else 
          call duplica
        end if
!      write(111,*)v(1:ktot+1),'&'
        write(3,*)n,ktot/2
        
!      allocate(tt(n),p2(n))
        kmax=maxval(k)
        allocate(ttk(kmax+1))
        ttk=0
        ddd=0
        do i=1,n
          call cluster2(i)
        end do
        ddd=sum(ttk)*1./6
        write(123123,*)n,ddd,ddd-dddd
        dddd=ddd
        deallocate(ttk)
!      deallocate(tt,p2)
        
      end do
      
      clfile(1:9)='cl_'//number//"_"
      pkfile(1:9)='pk_'//number//"_"
      lfile(1:9)='ll_'//number//"_"
!    open(111,file='intermednome')
!    write(111,'(I5)')in
      aaaa=''
      write(aaaa,'(I5.5)')n
!    read(111,'(A)')aaaa
      close(1111)
      clfile(10:14)=adjustl(aaaa)
      pkfile(10:14)=adjustl(aaaa)
      lfile(10:14)=adjustl(aaaa)
      write(103,*)clfile,' ',pkfile,' ',aaaa
      open(101,file=clfile,form='UNFORMATTED')
      open(102,file=pkfile,form='UNFORMATTED')
      open(104,file=lfile,form='UNFORMATTED')
      
      kmax=maxval(k)+1
!    allocate(pk(kmax),p2(n),p2k(kmax),tt(n),ttk(kmax))
      allocate(pk(kmax),p2k(kmax),lll(kmax),ttk(kmax),lk(kmax,kmax),dkk(kmax,kmax,kmax))
      kmax=maxval(k)
      pk=0
      p2=0
      p2k=0
      tt=0
      ttk=0
      dkk=0
      
      do i=1,n
        pk(k(i))=pk(k(i))+1
      end do
      
      do i=1,n
        call cluster3(i)
      end do
      do i=1,n
!    write(999,*)tt(i),k(i)
        ttk(k(i))=ttk(k(i))+tt(i)
        p2k(k(i))=p2k(k(i))+p2(i)
      end do
!  do i=1,kmax
!    if(p2k(i).gt.0)write(45,*)i,real(ttk(i))/(p2k(i)*pk(i))
!  end do
      allocate(clhist(kmax,2))
      !clhist=0
      !do i=1,n
      !clhist(k(i),1)=clhist(k(i),1)+1
      !clhist(k(i),2)=clhist(k(i),2)+cl(i)
      !end do
      
      do i=1,kmax
        write(102)i,real(pk(i))
      end do
      !do i=1,kmax
      !if((clhist(i,1).ne.0).and.(clhist(i,2).ne.0))then
      !write(101,*)i,clhist(i,2)/clhist(i,1),real(ttk(i)),real(p2k(i)),real(ttk(i))/(p2k(i)*pk(i)),real(ttk(i))/(i*(i-1)*pk(i))
      !end if
      !end do
      
      lll=0
      sk=0
      i=1
!    do j=1,k(i)
!      lll(k(i))=lll(k(i))+real(k(v(sk+j)))/k(i)
!    end do
!    do i=2,n!sum(k)
!      sk=sum(k(1:i-1))
!      !write(111,*)i,k(i),sum(k(1:i-1)),kmax
!      do j=1,k(i)
!        lll(k(i))=lll(k(i))+real(k(v(sk+j)))/k(i)
!      end do
!    end do
      
      
      do i=1,kmax
        do j=1,kmax
          write(104)i,j,lk(j,i)
        end do
      end do
      do l=1,kmax
        do i=1,kmax
          do j=1,kmax
            if(dkk(j,i,l).gt.0)write(101)i,l,j,dkk(j,i,l)*.5
          end do
        end do
      end do
      
      close(101)
      close(102)
      close(104)
      
      
!    deallocate(pk,p2,p2k,tt,ttk,clhist)
      deallocate(pk,p2k,ttk,clhist,lll,lk,dkk)
      
      
      !write(100+in,*)'##########################'
!    write(1000+in)n,pseed,cod,ll,ktot
      !do i=1,n
!    write(1000+in)k
!    write(1000+in)v
      !end do
!    write(1000+in)o
!    close(1000+in)
      
      
!    allocate(pk(n))   !makes the conectivity histogram
!    pk=0
!    do i=1,n
!      j=sum(m(i,:))
!      pk(j)=pk(j)+1
!    end do
!    do i=1,n
!      if(pk(i).ne.0)write(100+in,*)real(i)*2/real(ktot),real(pk(i))/real(n)
!    end do
!    deallocate(pk)
!    close(100)
!  end do
      
      
!  kmax=maxval(k)
!  allocate(clhist(kmax,2))
!  clhist=0
!  do i=1,n
!    !    write(*,*)'9'
!    clhist(k(i),1)=clhist(k(i),1)+1
!    clhist(k(i),2)=clhist(k(i),2)+cl(i)
!  end do
!  do i=1,kmax
!    !    write(*,*)'10'
!    !write(33,*)clhist(i,1),clhist(i,2)/clhist(i,1)
!    if((clhist(i,1).ne.0).and.(clhist(i,2).ne.0))write(55,*)real(i),clhist(i,2)/clhist(i,1)
!  end do
!  deallocate(clhist)
      
!  allocate(P(kmax,kmax))
!  p=0
!  do i=1,n
!    do j=1,k(i)
!      !      write(*,*)'11'
!      p(k(i),k(v(i,j)))=p(k(i),k(v(i,j)))+1
!      p(k(v(i,j)),k(i))=p(k(v(i,j)),k(i))+1
!    end do
!  end do
!  !  write(*,*)kmax
!  do i=1,kmax
!    soma=0
!    do j=1,kmax
!      !      write(*,*)'12'
!      soma=soma+j*p(j,i)
!    end do
!    write(8,*)i,soma/i
!    !    write(*,*)'oi'
    end do
    
    
    
    write(4,*)del
    deallocate(v,k,cl,o,del,tt,p2)
  end do
    
end program evol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inicia
  use globals
  implicit none
  integer,allocatable :: m(:,:)
  integer :: i, j ,l !loop indexes
  
  !open(12,file='inicia.dat',status='old')
  
  
  
  !read(12,*)pseed
  !allocate(m(n,n),o(2,n))
  !read(12,*)m
  !read(12,*)o
  !deallocate(m,o)
  n=5
  allocate(m(n,n),v(nl),k(n),cl(n),o(2,n),del(1))
  ktot=0
  k=0
  v=0
  m=0
  cl=0.d0
  
!  jj=n/(nn/10)
!  jj=n/(nn)
  
  
  call setseed
  !close(12)
  
  do i=1,n
    do j=i+1,n
      m(i,j)=1
      m(j,i)=1
      call liga(i,j)
      !call liga(j,i)
    end do
  end do
!  do i=1,n-1
!    m(i,i+1)=1
!    m(i+1,i)=1
!    call liga(i,i+1)
!  end do
!  do i=2,5
!    m(1,i)=1
!    m(i,1)=1
!    call liga(1,i)
!  end do
    
  o=0
  do i=1,n
    call cluster(i)
  end do
  
  
end subroutine inicia

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine liga(n1,n2) !links the two nodes
  use globals
  implicit none
  integer :: n1,n2,sk1,sk2      !nodes to be linked

  ktot=ktot+2		!number of links increase by two

  k(n1)=k(n1)+1		!updates the conectivity array
  sk1=sum(k(1:n1))
  v(sk1:ktot)=eoshift(v(sk1:ktot),-1)
  v(sk1)=n2	!updates the neighbors array

  k(n2)=k(n2)+1
  sk2=sum(k(1:n2))
  v(sk2:ktot)=eoshift(v(sk2:ktot),-1)
  v(sk2)=n1
end subroutine liga

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cluster(n1) !compute the clustering for the node
  
  !n1 -> node to compute the clustering of
  !n2 -> neighbors of n1
  !i,j,l -> loop indexes
  !combina -> normalization
  
  use globals
  implicit none
  integer :: n1,n2,i,j,l,sk1,sk2
  real*8 :: combina
  
  sk1=sum(k(1:n1-1))
  do i=1,k(n1)
    n2=v(sk1+i)		!n2 is one neighbor of n1
!    write(333,*)n1,k(n1),sk1,sk1+i,n2,sum(k(1:5)),v(1:16)
!    read(*,*)sk2
    sk2=sum(k(1:n2-1))
    Do j=1,k(n2)		!run trough all neighbors of n2
      Do l=1,k(n1)	!run trough all neighbors of n1
        if(v(sk2+j).eq.v(sk1+l))cl(n1)=cl(n1)+1	!if the j-neighbor of n2 is neighbor of the l-neighbor of n1
      enddo
    enddo
  enddo
  
  if(k(n1)>=2)then
    combina=k(n1)*(k(n1)-1.d0)/2.d0  !DUAS DIVISÕES POR DOIS?
    cl(n1)=cl(n1)/2.d0/combina
  else
    cl(n1)=0.d0
  endif
  
end subroutine cluster

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine horizontal
  
  !Subroutine for the horizontal transfer
  !link -> nodes wich the new node is linked to
  !i,j -> loop indexes
  !pl -> probability of linking
  !p -> parameter of probability
  !z -> random number
  
  use globals
  implicit none
  integer :: link(n),i,j,sk
  real*8,parameter :: p=1
  real*8 :: pl,z,norm
  
  link=0
  j=0
  i=0
  norm=(sum(k**(0.5)))
  pl=0
  call random_number(z)
  do while(j==0)
    i=i+1
    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
    j=1
    link(j)=i  !keeps all the nodes that links
    if (z>pl) then
      j=0
    endif
  end do

  i=0
  pl=0
  call random_number(z)
  do while(j==1)
    i=i+1
    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
    j=2
    link(j)=i  !keeps all the nodes that links
    if (z>pl) then
      j=1
    endif
  end do

  i=0
  pl=0
  call random_number(z)
  do while(j==2)
    i=i+1
    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
    j=3
    link(j)=i  !keeps all the nodes that links
    if (z>pl) then
      j=2
    endif
  end do
!
!  i=0
!  pl=0
!  call random_number(z)
!  do while(j==3)
!    i=i+1
!    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
!    j=4
!    link(j)=i  !keeps all the nodes that links
!    if (z>pl) then
!      j=3
!    endif
!  end do
!
!  i=0
!  pl=0
!  call random_number(z)
!  do while(j==4)
!    i=i+1
!    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
!    j=5
!    link(j)=i  !keeps all the nodes that links
!    if (z>pl) then
!      j=4
!    endif
!  end do
!  
!  i=0
!  pl=0
!  call random_number(z)
!  do while(j==5)
!    i=i+1
!    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
!    j=6
!    link(j)=i  !keeps all the nodes that links
!    if (z>pl) then
!      j=5
!    endif
!  end do
!  
!  i=0
!  pl=0
!  call random_number(z)
!  do while(j==6)
!    i=i+1
!    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
!    j=7
!    link(j)=i  !keeps all the nodes that links
!    if (z>pl) then
!      j=6
!    endif
!  end do
!  
!  i=0
!  pl=0
!  call random_number(z)
!  do while(j==7)
!    i=i+1
!    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
!    j=8
!    link(j)=i  !keeps all the nodes that links
!    if (z>pl) then
!      j=7
!    endif
!  end do
!  
!  i=0
!  pl=0
!  call random_number(z)
!  do while(j==8)
!    i=i+1
!    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
!    j=9
!    link(j)=i  !keeps all the nodes that links
!    if (z>pl) then
!      j=8
!    endif
!  end do
!  
!  i=0
!  pl=0
!  call random_number(z)
!  do while(j==9)
!    i=i+1
!    pl=pl+(k(i)**(0.5))/norm
!    write(666,*)i,k(i),z,pl,j
!    j=10
!    link(j)=i  !keeps all the nodes that links
!    if (z>pl) then
!      j=9
!    endif
!  end do
  
  if(j>0) then
    call atumatriz
    flag=1
    n=n+1 !adds just one node
    do i=1,j
      call liga(n,link(i))
    end do
    
    call cluster(n)
    sk=sum(k(1:n-1))
    do i=1,k(n)
      call cluster(v(sk+i))
    end do
    o(2,n)=0
    write(2,*)n,ll,'h',k(n),cl(n)
  end if

  
end subroutine horizontal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine duplica
  
  !d -> father node
  !h -> neighbor of the new node
  !i,ii,j,l -> loop index
  !duplic -> duplicated nodes array
  !norm -> normalization
  !z -> random variable
  !pm -> mutation probability
  !zb -> duplication probability
  !p -> paramater
  
  use globals
  implicit none
  integer :: d,h,i,ii,j,l,duplic(n),km,sk
  real*8 :: norm,z,zb,p
  
  p=1
  zb=0
  j=0
  duplic=0
  norm=sum((1-cl)/k)
  !do i=1,n
  !   zb(i)=zb(i-1)+(1.0-cl(i))/k(i)/norm
  !enddo
  do i=1,n
    zb=p*(1-cl(i))/k(i)/norm
    call random_number(z)
    if (z<zb)then
      j=j+1
      duplic(j)=i !keeps all the father genes
    end if
  end do
  
  if (j>0) then
    call atumatriz
    n=n+1
    call random_number(z)
    d=ceiling(z*j)
    d=duplic(d)
    call ligadup(d)
    
    !   call random_number(z) !mutation
    !   ii=ceiling(z*n)
    !   if(m(n,ii)==1)then
    !      call desliga(i,n)
    !   else
    !      call liga(ii,n)
    !   end if
    
    
    call cluster(n)
    sk=sum(k(1:n-1))
    do i=1,k(n)
      h=v(sk+i)
      call cluster(h)
    end do
    call cluster(d)
    sk=sum(k(1:d-1))
    do i=1,k(d)
      h=v(sk+i)
      call cluster(h)
    end do
    
    
    write(2,*)n,ll,'d',k(n),k(d),d,cl(n),cl(d)    
    o(2,n)=duplic(1)
    
  end if
  
end subroutine duplica

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine atumatriz
  use globals
  implicit none
  
  
  allocate(kold(n),oold(2,n),clold(n))
  
  kold=k
  clold=cl
  oold=o
  deallocate(k,cl,o)
  allocate(k(n+1),cl(n+1),o(2,n+1))
  k=0
  cl=0.d0
  o=0
  k(1:n)=kold
  cl(1:n)=clold
  o(1:2,1:n)=oold
  o(1,n+1)=cod
  deallocate(kold,oold,clold)
  
end subroutine atumatriz

!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setseed
  use globals
  implicit none
  
  integer :: I,j
  integer,allocatable :: seed(:)
  real :: z
  if(pseed.ne.0)then
    call random_seed(size=I)
    allocate(seed(I))
    seed=pseed
    call random_seed(put =seed)
    
    do j=1,I
      call random_number(z)
      seed(j)=nint(999999999*z)
    end do
    call random_seed(put =seed)
  end if
end subroutine setseed

!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ligadup(d)
  use globals
  implicit none

  integer :: d,i,ii,j,numv,kk,dd,sk,lig(k(d)),des(k(d))
  integer,allocatable :: viz(:)
  real*8 :: z
  
  call liga(d,n)
  sk=sum(k(1:d-1))
  kk=k(d)
  allocate(viz(kk))
  viz(:)=v(sk+1:sk+kk)
  
  do ii=1,kk-1 !link the new node with all neighbors of the father
    i=viz(ii)
    call liga(i,n)
  end do

!  write(444,*)viz,'&',v(sk+1:sk+kk)
  do ii=1,kk-1
    i=viz(ii)
    call random_number(z)
    if(z<frac/2.d0)then
      call desliga(d,i)
    else if(z<frac)then
      call desliga(n,i)
    end if
  
  !do ii=1,kk-1
  !  i=viz(ii)
  !  call random_number(z)
  !  if(z<frac/2.d0)then
  !    call desliga(d,i)
  !  else if(z<frac)then
  !    call desliga(n,i)
  !  end if
    
    
  end do


end subroutine ligadup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine desliga(n1,n2)
  use globals
  implicit none
  
  integer :: i,n1,n2,i1,i2,sk1,sk2
  i1=-1000000000
!  write(*,*)n1,n2
  sk1=sum(k(1:n1-1))
  do i=1,k(n1)
    if(v(sk1+i)==n2) i1=i
  end do
!  write(333,*)n1,n2,sk1,k(n1),i1,ktot,n
  k(n1)=k(n1)-1		!updates the conectivity array
  v(sk1+i1:ktot)=eoshift(v(sk1+i1:ktot),1)

  sk2=sum(k(1:n2-1))
  do i=1,k(n2)
    if(v(sk2+i)==n1) i2=i
  end do
  
  k(n2)=k(n2)-1
  v(sk2+i2:ktot)=eoshift(v(sk2+i2:ktot),1)
  
!  v(n1,i1:k(n1)+1)=v(n1,i1+1:k(n1)+2)	!updates the neighbors array
!  v(n2,i2:k(n2)+1)=v(n2,i2+1:k(n2)+2)
  ktot=ktot-2	

end subroutine desliga

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine cluster2(n1) !compute the clustering for the node

  !n1 -> node to compute the clustering of
  !n2 -> neighbors of n1
  !i,j,l -> loop indexes
  !combina -> normalization

  use globals
  implicit none
  integer :: n1,n2,i,j,l,sk1,sk2
  real*8 :: combina

  sk1=sum(k(1:n1-1))
  tt(n1)=0
  do i=1,k(n1)
    n2=v(sk1+i)         !n2 is one neighbor of n1
    !      if(n2.eq.0)write(888,*)n1,i,k(n1),v(n1,i)
    sk2=sum(k(1:n2-1))
    Do j=1,k(n2)               !run trough all neighbors of n2
      Do l=1,k(n1)    !run trough all neighbors of n1
        if(v(sk2+j).eq.v(sk1+l))cl(n1)=cl(n1)+1 !if the j-neighbor of n2 is neighbor of the l-neighbor of n1
        if(v(sk2+j).eq.v(sk1+l))then !if the j-neighbor of n2 is neighbor of the l-neighbor of n1
          ttk(k(n1))=ttk(k(n1))+1
          ttk(k(v(sk2+j)))=ttk(k(v(sk2+j)))+1 
          ttk(k(v(sk1+l)))=ttk(k(v(sk1+l)))+1
        end if
      enddo
    enddo
  enddo
  p2(n1)=k(n1)*(k(n1)-1)/2
  if(k(n1)>=2)then
    combina=real(k(n1))*(real(k(n1))-1.d0)/2.d0  !DUAS DIVIS<D5>ES POR DOIS?
    cl(n1)=cl(n1)/(2.d0*combina)
  else
    cl(n1)=0.d0
  endif

end subroutine cluster2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cluster3(n1) !compute the clustering for the node

  !n1 -> node to compute the clustering of
  !n2 -> neighbors of n1
  !i,j,l -> loop indexes
  !combina -> normalization

  use globals
  implicit none
  integer :: n1,n2,i,j,l,sk1,sk2
  real*8 :: combina

  sk1=sum(k(1:n1-1))
  tt(n1)=0
  do i=1,k(n1)
    n2=v(sk1+i)         !n2 is one neighbor of n1
    !      if(n2.eq.0)write(888,*)n1,i,k(n1),v(n1,i)
    lk(k(n1),k(n2))=lk(k(n1),k(n2))+0.5
    sk2=sum(k(1:n2-1))
    Do j=1,k(n2)               !run trough all neighbors of n2
      Do l=1,k(n1)    !run trough all neighbors of n1
        if(v(sk2+j).eq.v(sk1+l))cl(n1)=cl(n1)+1 !if the j-neighbor of n2 is neighbor of the l-neighbor of n1
        if(v(sk2+j).eq.v(sk1+l))then !if the j-neighbor of n2 is neighbor of the l-neighbor of n1
          ttk(k(n1))=ttk(k(n1))+1
          ttk(k(v(sk2+j)))=ttk(k(v(sk2+j)))+1
          ttk(k(v(sk1+l)))=ttk(k(v(sk1+l)))+1
          dkk(k(n1),k(n2),k(v(sk2+j)))=dkk(k(n1),k(n2),k(v(sk2+j)))+1
        end if
      enddo
    enddo
  enddo
  p2(n1)=k(n1)*(k(n1)-1)/2
  if(k(n1)>=2)then
    combina=real(k(n1))*(real(k(n1))-1.d0)/2.d0  !DUAS DIVIS<D5>ES POR DOIS?
    cl(n1)=cl(n1)/(2.d0*combina)
  else
    cl(n1)=0.d0
  endif
  
end subroutine cluster3
