
   PROGRAM PGremod
   use declare
   use mods_remod
   IMPLICIT NONE

!  SETTING INFORMATION:

   CALL SETINFO(JOBID,NCYL,NUNIT,LAVE,GMODE,LEFTCELL,RITECELL)


   WRITE(*,*)'=========================================================='
   WRITE(*,*)'Remodeler version 1'

   IF(JOBID==1)THEN
      PRINT*,'To generate the start sacculus ...'

   ELSEIF(JOBID==2)THEN
      PRINT*,'Growth of sacculus ....'
   ELSE
      PRINT*,'Division of sacculus ....'
   END IF

!--------------------------------------------
!  Must comment out these for serial jobs:
   nthread=4
   call omp_set_num_threads(nthread)
!--------------------------------------------

   WRITE(*,*)'--------------------------------------'

!  COMMON PARAMETERS:
   PI=3.141592653589793239D0

!  Spring constant and relaxed length for glycan:
   K_G=577.0D0 ! unit=1e-20 J/nm**2
   L_G=2.005D0 ! unit=nm

!  Worm-like chain model force constant for peptide crosslink:
   K_P=1.49D0

!  Effective contour length for crosslink:
   L_P=3.8D0

!  Bending stiffness for glycan:
   KTHETA=8.36D0 ! unit=1e-20 J
   THETA_0=PI

!  Turgor pressure:
   PRES=0.03D0   ! unit=1e-20 J/nm**3

!  Generate random seed:
   CALL INIT_RANDOM_SEED()

!  To estimate running time of the job:
   CALL WHATTIME(TIMESTART)

52 FORMAT(A4,2X,I9,4X,A15,2X,F6.3,2x,F6.3)
91 FORMAT(A14,2X,I2,1X,A7,I2,1X,A8,1X,I2,1X,A4)

   IF(JOBID==2)THEN
      GOTO 72
   ELSEIF(JOBID==3)THEN
      GOTO 73
   END IF

!===================================================================

!  CREATION OF A SACCULUS

!   PRINT*,'NUMBER OF HOOPS ON THE CYLINDER'
!   READ(*,*)NCYL
!  NCYL=64

!   PRINT*,'NUMBER OF UNITS PER HOOP'
!   READ(*,*)NUNIT
!  NUNIT=100

!   PRINT*,'AVERAGE STRAND LENGTH'
!   READ(*,*)LAVE
!   LAVE=10

   NCLOSEDHOOP=2

!  Set crosslink length:
   L_PEP=2.0D0

   IF(NUNIT>1000)THEN
      WRITE(*,*)'WANRNING: THERE ARE MORE THAN 1000 UNITS PER HOOP '
      WRITE(*,*)'PICK A SMALLER SIZE TO RUN A JOB PROPERLY'
      STOP
   END IF

!=====================================================================================
!  number of hoops per each polar cap:
   NCAP=NUNIT*L_G/L_PEP/4

!  We estimate the total number of beads in the system to be less than:
   NESTI=(NCYL+2*NCAP)*NUNIT

!  Coordinates of beads:
   ALLOCATE(X(NESTI),Y(NESTI),Z(NESTI),APOST(NESTI),ANGLE(NESTI))   

!  Direction of peptides on beads:
   ALLOCATE(PEPDIR(NESTI))

!  This is used to distinguish types of beads:
   ALLOCATE(ATYP(NESTI))

!  Bonds of glycan and peptide:
   ALLOCATE(BONDGLY(2,NESTI),BONDPEP(2,NESTI))

!  Types of peptide bonds:
   ALLOCATE(BONTYP(NESTI))

!  Bonds in general:
   ALLOCATE(BOND(2,2*NESTI))

!  Partners of beads, connected by bonds:
   ALLOCATE(NPART(NESTI),PART(3,NESTI))

!  PG strands are represented as PGID:
   ALLOCATE(PGID(NUNIT,2*NESTI/LAVE))
!  Strand length are PGLEN:
   ALLOCATE(PGLEN(2*NESTI/LAVE)); PGLEN=0
!  Strand types are PGTYP:
   ALLOCATE(PGTYP(2*NESTI/LAVE),XPOST(2*NESTI/LAVE))
!  Strand directions are PGDIR:
   ALLOCATE(PGDIR(2*NESTI/LAVE))
   PGDIR=1

!  Pore periphery is called LOOP, length=LOOPLEN, type=LOOPTYP
   ALLOCATE(LOOP(1000,NESTI),LOOPLEN(NESTI))
   ALLOCATE(LOOPTYP(NESTI))
   LOOPTYP=1

!  To mark beads as 0 or 1 for convenience
   ALLOCATE(MARK(2*NESTI))


!  If there are NUNIT glycan subunits per hoop, then the radius
!  of the cylinder is:

   RADIUS=NUNIT*L_G/2.0D0/PI

!  If there are NCYL hoops and the peptide cross-link length is 
!  L_PEP then the length of the cylinder is:

   LCYL=(NCYL-1)*L_PEP

!  Let choose the center of mass at (0,0,0) then the center of 
!  hoop Nth is at (-L_CYN/2+(N-1)*L_PEP,0,0).

!  number of atoms (beads) is NATOM, NPG is
!  the number of residues (glycan strands).

!========================================================================
!************** NOW WE SET UP COORDINATES FOR THE LEFT CAP:

   NREPEAT=0

20 NATOM=0; NPG=0


   DO N=1,NCAP


      NPG=NPG+1

      PGTYP(NPG)=0

      PHI=N*PI/NCAP/2

! CENTER OF THE HOOP IS AT 

      X0=-LCYL/2.0D0-L_PEP-RADIUS*COS(PHI)


! RADIUS OF THE HOOP IS:

      RAD=RADIUS*SIN(PHI)

      IF(RAD<0.0D0)THEN
         PRINT*,'ERROR: NEGATIVE RADIUS OF CAP HOOPS'
         STOP
      END IF

! NUMBER OF GLYCAN UNITS ON THIS HOOP IS:

      PGLEN(NPG)=INT((2.0D0*PI*RAD+0.000001D0)/L_G)


      IF(PGLEN(NPG)<5)THEN
         PRINT*,'ERROR: HOOP IS TOO SMALL'
         STOP
      END IF

      IF(MOD(PGLEN(NPG),2)==1)THEN
         PGLEN(NPG)=PGLEN(NPG)-1
      END IF


! In the y-z plane, position of atoms can be defined by angle
! Phi with an interval:

      DPHI=2*PI/PGLEN(NPG); 

      CALL RANDOM_NUMBER(R); PHI0=2*PI*R

      DO I=0,PGLEN(NPG)-1

         NATOM=NATOM+1; PGID(I+1,NPG)=NATOM
         PHI=I*DPHI+PHI0

         IF(PHI>=2.0D0*PI)THEN
            PHI=PHI-2.0D0*PI
         END IF

         X(NATOM)=X0
         Y(NATOM)=RAD*COS(PHI)
         Z(NATOM)=RAD*SIN(PHI)

         ANGLE(NATOM)=PHI


      END DO

   END DO


!==================================================================
! NOW WE SET UP COORDINATES FOR THE CYLINDER:

! In the y-z plane, position of atoms can be defined by angle
! Phi with an interval:

   DPHI=2*PI/NUNIT

! Each glycan strand (or residue) will be chosen to have its length 
! between NDOWN and NUP and average length is LAVE:

   NUP=LAVE+LAVE/2; NDOWN=LAVE-LAVE/2

   DO N=1,NCYL


      X0=-LCYL/2+(N-1)*L_PEP

! We start the first atom of the hoop at a random position:

      CALL RANDOM_NUMBER(R); NPHI0=NUNIT*R

      IF(MOD(NPHI0,2)==1)THEN
         NPHI0=NPHI0+1
      END IF

      IF(NPHI0>=NUNIT)THEN
         NPHI0=NPHI0-NUNIT
      END IF

! But not the first and last hoops of the cylinder:

      IF(N<=NCLOSEDHOOP.OR.N>=NCYL-NCLOSEDHOOP+1)THEN
         NPHI0=0
      END IF

      NCOUNT=NUNIT

      NPOST=NPHI0

      DO N0=1,1000

         IF(NCOUNT<=0)THEN
            EXIT
         END IF

         CALL RANDOM_NUMBER(R)
         LENGTH=NDOWN+R*(NUP-NDOWN)


         NLEFT=NCOUNT-LENGTH

         IF(NLEFT<NDOWN)THEN
            LENGTH=NCOUNT
         END IF

         NPG=NPG+1
         XPOST(NPG)=N

         IF(MOD(N,2)==1)THEN
            PGTYP(NPG)=1
         ELSE
            PGTYP(NPG)=2
         END IF

         IF(N<=NCLOSEDHOOP.OR.N>=NCYL-NCLOSEDHOOP+1)THEN
            PGTYP(NPG)=0
            LENGTH=NUNIT
         END IF

         PGLEN(NPG)=LENGTH

         DO J=1,LENGTH

            NATOM=NATOM+1


            PGID(J,NPG)=NATOM

            X(NATOM)=X0

            PHI=NPOST*DPHI

            ANGLE(NATOM)=PHI

            Y(NATOM)=RADIUS*COS(PHI)

            Z(NATOM)=RADIUS*SIN(PHI)

            APOST(NATOM)=NPOST

            NPOST=NPOST+1

            IF(NPOST>=NUNIT)THEN
               NPOST=NPOST-NUNIT
            END IF

         END DO

         NCOUNT=NCOUNT-LENGTH


         IF(NPOST>=NUNIT)THEN
            NPOST=NPOST-NUNIT
         END IF

      END DO

      IF(N==NCLOSEDHOOP)THEN
         NATOMCAP=2*NATOM
         PGCAP1=NPG
      END IF

      IF(N==NCYL-NCLOSEDHOOP+1)THEN
         PGCAP2=NPG
      END IF

   END DO

!==========================================================================

! NOW WE SET UP COORDINATES FOR THE RIGHT CAP:

   DO N=1,NCAP


      NPG=NPG+1

      PGTYP(NPG)=0

      PHI=(N-1)*PI/NCAP/2

! CENTER OF THE HOOP IS AT 

      X0=LCYL/2.0D0+L_PEP+RADIUS*SIN(PHI)

! RADIUS OF THE HOOP IS:

      RAD=RADIUS*COS(PHI)

      IF(RAD<0.0D0)THEN
          PRINT*,'ERROR: NEGATIVE RADIUS OF CAP HOOPS'
         STOP
      END IF

! NUMBER OF GLYCAN UNITS ON THIS HOOP IS:

      PGLEN(NPG)=INT((2.0D0*PI*RAD+0.000001D0)/L_G)

      IF(MOD(PGLEN(NPG),2)==1)THEN
         PGLEN(NPG)=PGLEN(NPG)-1
      END IF

! In the y-z plane, position of atoms can be defined by angle
! Phi with an interval:

      DPHI=2*PI/PGLEN(NPG)
      CALL RANDOM_NUMBER(R); PHI0=2*PI*R

      DO I=0,PGLEN(NPG)-1

         NATOM=NATOM+1; PGID(I+1,NPG)=NATOM
         PHI=I*DPHI+PHI0

         IF(PHI>=2.0D0*PI)THEN
            PHI=PHI-2.0D0*PI
         END IF

         X(NATOM)=X0
         Y(NATOM)=RAD*COS(PHI)
         Z(NATOM)=RAD*SIN(PHI)

         ANGLE(NATOM)=PHI


      END DO

   END DO



!===============================================================================


! NOW ASSIGN BONDS TO THE SYSTEM:
    

   ATYP=3


!  We start with the GLY-GLy bonds. This is straight forward:

   NBONDGLY=0
      
   DO N=1,NPG

      LENGTH=PGLEN(N)-1

!       For cap hoops at the ends, they are closed so there is a bond
!       between the first and the last atoms in the hoops:

      IF(PGTYP(N)==0)THEN
         LENGTH=LENGTH+1
      END IF

      DO IPG=1,LENGTH
 
         I=PGID(IPG,N); J=I+1
         IF(IPG==PGLEN(N))THEN
            J=PGID(1,N)
         END IF

         NBONDGLY=NBONDGLY+1; BONDGLY(1,NBONDGLY)=I; BONDGLY(2,NBONDGLY)=J

      END DO


   END DO

!===================

!  Peptide bonds:

   PEPDIR=0

   NBONDPEP=0

! -- START WITH THE CYLINDER


   DO N=NCAP+1,NPG-NCAP-1

      DO J=1,PGLEN(N)

         N1=PGID(J,N)

         IF(ATYP(N1)==1)THEN
            CYCLE
         END IF

         IF(N==NCAP+1.AND.MOD(J-1,2)/=0)THEN
            CYCLE
         END IF

         PEPDIR(N1)=1

         STOPA=0


         DO N0=N+1,NPG-NCAP

            IF(XPOST(N0)>XPOST(N)+1)THEN
               CYCLE
            END IF

            DO J0=1,PGLEN(N0)

               N2=PGID(J0,N0)

               IF(APOST(N2)==APOST(N1))THEN

                  NBONDPEP=NBONDPEP+1

                  BONDPEP(1,NBONDPEP)=N1

                  BONDPEP(2,NBONDPEP)=N2

                  ATYP(N1)=1

                  ATYP(N2)=1

                  PEPDIR(N2)=-1

                  STOPA=1

                  EXIT
               END IF

            END DO

            IF(STOPA==1)THEN
               EXIT
            END IF

         END DO

      END DO

   END DO

!-----------------------------

! -- CONTINUE WITH THE FISRT CAP

   ALLOCATE(JMARK(NUNIT),JMARKMIRROR(NUNIT))

   DO NC1=NCAP,1,-1

      JMARK=0
      JMARKMIRROR=0

      JB1=0
      JB2=0

      NC2=NC1+1

      CHECK=0


      DO J1=1,PGLEN(NC1)/2

         IF(CHECK==1.AND.NC1/=1)THEN
            CHECK=0
            CYCLE
         END IF


         N1=PGID(J1,NC1)


         PHI0=PI


         DO J2=1,PGLEN(NC2)

            IF(JMARK(J2)==1)THEN
               CYCLE
            END IF

            N2=PGID(J2,NC2)


            IF(ATYP(N2)==1)THEN
               CYCLE
            END IF

            DPHI=ABS(ANGLE(N1)-ANGLE(N2))

            IF(DPHI>PI)THEN
               DPHI=2.0D0*PI-DPHI
            END IF

            IF(DPHI<PHI0)THEN
               PHI0=DPHI
               NPICK=N2
               JPICK=J2
            END IF

         END DO


         IF(PHI0>2*PI/PGLEN(NC2))THEN
            CYCLE
         END IF

         NBONDPEP=NBONDPEP+1

         BONDPEP(1,NBONDPEP)=N1

         BONDPEP(2,NBONDPEP)=NPICK

         ATYP(N1)=1

         ATYP(NPICK)=1

         PEPDIR(N1)=1

         PEPDIR(NPICK)=-1

         CHECK=1

         IF(JB1==0)THEN
            JB1=JPICK
         END IF

!        MIRROR SYMMETRIC PEPTIDE BOND:

         J=J1+PGLEN(NC1)/2


         N1SYM=PGID(J,NC1)

         IF(ATYP(N1SYM)==1.AND.J<PGLEN(NC1))THEN

            J=J+1


            N1SYM=PGID(J,NC1)

         END IF

         IF(ATYP(N1SYM)==1)THEN
            CYCLE
         END IF

         J=JPICK+PGLEN(NC2)/2

         IF(J>PGLEN(NC2))THEN
            J=J-PGLEN(NC2)
         END IF

         N2SYM=PGID(J,NC2)

         IF(ATYP(N2SYM)==1)THEN

            J=J+1

            IF(J>PGLEN(NC2))THEN
               J=J-PGLEN(NC2)
            END IF

            N2SYM=PGID(J,NC2)

         END IF

         IF(JMARKMIRROR(J)==1)THEN
            CYCLE
         END IF

         IF(ATYP(N2SYM)==1)THEN
            CYCLE
         END IF

         NBONDPEP=NBONDPEP+1

         BONDPEP(1,NBONDPEP)=N1SYM

         BONDPEP(2,NBONDPEP)=N2SYM

         ATYP(N1SYM)=1

         ATYP(N2SYM)=1

         PEPDIR(N1SYM)=1

         PEPDIR(N2SYM)=-1

         IF(JB2==0)THEN
            JB2=J

!           MARK NON-BOND ZONES:

            DJ=JB2-JB1

            IF(DJ<0)THEN
               DJ=DJ+PGLEN(NC2)
            END IF

            DO J=JB1,JB1+DJ

               IF(J<=PGLEN(NC2))THEN
                  JMARKMIRROR(J)=1
               ELSE
                  JMARKMIRROR(J-PGLEN(NC2))=1
               END IF

            END DO

            DJ=JB1-JB2

            IF(DJ<0)THEN
               DJ=DJ+PGLEN(NC2)
            END IF

            DO J=JB2,JB2+DJ

               IF(J<=PGLEN(NC2))THEN
                  JMARK(J)=1
               ELSE
                  JMARK(J-PGLEN(NC2))=1
               END IF

            END DO

         END IF

      END DO

   END DO


!--------------------------------------------

! -- CONTINUE WITH THE SECOND CAP


   DO NC1=NPG-NCAP+1,NPG

      JMARK=0
      JMARKMIRROR=0

      JB1=0
      JB2=0

      NC2=NC1-1

      CHECK=0

      DO J1=1,PGLEN(NC1)/2

         IF(CHECK==1.AND.NC1/=NPG)THEN
            CHECK=0
            CYCLE
         END IF


         N1=PGID(J1,NC1)

         PHI0=PI


         DO J2=1,PGLEN(NC2)

            IF(JMARK(J2)==1)THEN
               CYCLE
            END IF

            N2=PGID(J2,NC2)

            IF(ATYP(N2)==1)THEN
               CYCLE
            END IF

            DPHI=ABS(ANGLE(N1)-ANGLE(N2))

            IF(DPHI>PI)THEN
               DPHI=2.0D0*PI-DPHI
            END IF

            IF(PHI0>DPHI)THEN
               PHI0=DPHI
               NPICK=N2
               JPICK=J2
            END IF

         END DO

         IF(PHI0>2*PI/PGLEN(NC2))THEN
            CYCLE
         END IF

         NBONDPEP=NBONDPEP+1

         BONDPEP(1,NBONDPEP)=N1
         BONDPEP(2,NBONDPEP)=NPICK

         ATYP(N1)=1
         ATYP(NPICK)=1

         PEPDIR(N1)=-1

         PEPDIR(NPICK)=1

         CHECK=CHECK+1

         IF(JB1==0)THEN
            JB1=JPICK
         END IF

!        MIRROR SYMMETRIC PEPTIDE BOND:

         J=J1+PGLEN(NC1)/2


         N1SYM=PGID(J,NC1)

         IF(ATYP(N1SYM)==1.AND.J<PGLEN(NC1))THEN

            J=J+1


            N1SYM=PGID(J,NC1)

         END IF

         IF(ATYP(N1SYM)==1)THEN
            CYCLE
         END IF

         J=JPICK+PGLEN(NC2)/2

         IF(J>PGLEN(NC2))THEN
            J=J-PGLEN(NC2)
         END IF

         N2SYM=PGID(J,NC2)

         IF(ATYP(N2SYM)==1)THEN

            J=J+1

            IF(J>PGLEN(NC2))THEN
               J=J-PGLEN(NC2)
            END IF

            N2SYM=PGID(J,NC2)

         END IF

         IF(ATYP(N2SYM)==1)THEN
            CYCLE
         END IF

         IF(JMARKMIRROR(J)==1)THEN
            CYCLE
         END IF

         NBONDPEP=NBONDPEP+1

         BONDPEP(1,NBONDPEP)=N1SYM
         BONDPEP(2,NBONDPEP)=N2SYM

         ATYP(N1SYM)=1
         ATYP(N2SYM)=1

         PEPDIR(N1SYM)=-1

         PEPDIR(N2SYM)=1

         IF(JB2==0)THEN
            JB2=J

!           MARK NON-BOND ZONES:

            DJ=JB2-JB1

            IF(DJ<0)THEN
               DJ=DJ+PGLEN(NC2)
            END IF

            DO J=JB1,JB1+DJ

               IF(J<=PGLEN(NC2))THEN
                  JMARKMIRROR(J)=1
               ELSE
                  JMARKMIRROR(J-PGLEN(NC2))=1
               END IF

            END DO

            DJ=JB1-JB2

            IF(DJ<0)THEN
               DJ=DJ+PGLEN(NC2)
            END IF

            DO J=JB2,JB2+DJ

               IF(J<=PGLEN(NC2))THEN
                  JMARK(J)=1
               ELSE
                  JMARK(J-PGLEN(NC2))=1
               END IF

            END DO

         END IF

      END DO

   END DO

   BONTYP(1:NBONDPEP)=1

   DEALLOCATE(JMARK,JMARKMIRROR)
!-------------------------------------------


   NBOND=NBONDGLY+NBONDPEP


   BOND(1,1:NBONDGLY)=BONDGLY(1,1:NBONDGLY)
   BOND(2,1:NBONDGLY)=BONDGLY(2,1:NBONDGLY)

   BOND(1,NBONDGLY+1:NBOND)=BONDPEP(1,1:NBONDPEP)
   BOND(2,NBONDGLY+1:NBOND)=BONDPEP(2,1:NBONDPEP)

!------------------------------------------------
   NPART=0
   PART=0

   DO N=1,NBOND
      N1=BOND(1,N)
      N2=BOND(2,N)

      NPART(N1)=NPART(N1)+1
      PART(NPART(N1),N1)=N2

      NPART(N2)=NPART(N2)+1

      PART(NPART(N2),N2)=N1
   END DO

!=======================================================================================

!  DEFINE LOOPS:

   LOOPLEN=0

! THE FIRST LOOP IS THE FIRST CLOSED HOOP:

   NLOOP=1; LOOPLEN(1)=PGLEN(1)
   LENGTH=LOOPLEN(1)
 
   DO J=1,LENGTH
      LOOP(J,1)=PGID(LENGTH-J+1,1)
   END DO


! WE USE A PEPTIDE BOND TO START A LOOP. EACH PEPTIDE BOND IS ASSOCIATED WITH
! TWO LOOPS. WE USE MARKING TO TELL IF A PEPTIDE BOND HAS BEEN ASSIGNED TO A
! LOOP OR TWO.

   MARK=0

! FIND CENTER OF THE CELL:

   XCEN0=SUM(X(1:NATOM))/NATOM
   YCEN0=SUM(Y(1:NATOM))/NATOM
   ZCEN0=SUM(Z(1:NATOM))/NATOM



   DO N0=NBONDGLY+1,NBOND
      IF(MARK(N0)<2)THEN

         NLOOP=NLOOP+1

         ILOOP=2; LOOP(1,NLOOP)=BOND(2,N0); LOOP(2,NLOOP)=BOND(1,N0)

         DO NRUN=1,1000
            N=LOOP(ILOOP,NLOOP)
            IF(NPART(N)==2)THEN
               I=PART(1,N); J=PART(2,N)
               IF(I==LOOP(ILOOP-1,NLOOP))THEN
                  NEXT=J
               ELSE
                  NEXT=I
               END IF
            END IF

            IF(NPART(N)==3)THEN
               I=PART(1,N); J=PART(2,N); K=PART(3,N)
               IF(I==LOOP(ILOOP-1,NLOOP))THEN
                  NEXT1=J; NEXT2=K
               ELSEIF(J==LOOP(ILOOP-1,NLOOP))THEN
                  NEXT1=K; NEXT2=I
               ELSE
                  NEXT1=J; NEXT2=I
               END IF

               X1=X(N)-X(LOOP(ILOOP-1,NLOOP))
               Y1=Y(N)-Y(LOOP(ILOOP-1,NLOOP))
               Z1=Z(N)-Z(LOOP(ILOOP-1,NLOOP))

               X2=X(NEXT1)-X(N); Y2=Y(NEXT1)-Y(N); Z2=Z(NEXT1)-Z(N)

               NX=Y1*Z2-Y2*Z1; NY=Z1*X2-Z2*X1; NZ=X1*Y2-X2*Y1

               W1=NX*(X(N)-XCEN0)+NY*(Y(N)-YCEN0)+NZ*(Z(N)-ZCEN0)

               X2=X(NEXT2)-X(N); Y2=Y(NEXT2)-Y(N); Z2=Z(NEXT2)-Z(N)

               NX=Y1*Z2-Y2*Z1; NY=Z1*X2-Z2*X1; NZ=X1*Y2-X2*Y1

               W2=NX*(X(N)-XCEN0)+NY*(Y(N)-YCEN0)+NZ*(Z(N)-ZCEN0)

               IF(W1>W2)THEN
                  NEXT=NEXT1
               ELSE
                  NEXT=NEXT2
               END IF
            END IF

            IF(NEXT==LOOP(1,NLOOP))THEN
               EXIT
            END IF

            ILOOP=ILOOP+1; LOOP(ILOOP,NLOOP)=NEXT

            DO IBOND=1,NBOND
               IF(BOND(1,IBOND)==N.AND.BOND(2,IBOND)==NEXT)THEN
                  MARK(IBOND)=MARK(IBOND)+1
                  EXIT
               ELSE IF(BOND(1,IBOND)==NEXT.AND.BOND(2,IBOND)==N)THEN
                  MARK(IBOND)=MARK(IBOND)+1
                  EXIT
               END IF
            END DO

         END DO

         LOOPLEN(NLOOP)=ILOOP

      END IF

   END DO


! THE LAST LOOP IS THE LAST CLOSED HOOP:

   NLOOP=NLOOP+1; LOOPLEN(NLOOP)=PGLEN(NPG)
   DO I=1,PGLEN(NPG)
      LOOP(I,NLOOP)=PGID(I,NPG)
   END DO

!  NUMBER OF LOOPS = NBONDPEP - NUMBER OF SHORT STRANDS + 2 LOOPS AT THE ENDS

   IF(NLOOP/=NBONDPEP-(NPG-2*(NCAP+NCLOSEDHOOP))+2)THEN
      PRINT*,'ERROR IN FINDING LOOPS! STOP NOW.'
      STOP
   END IF

   JREPEAT=0

   DO N=1,NLOOP
      IF(LOOPLEN(N)>14)THEN
         JREPEAT=JREPEAT+(LOOPLEN(N)-14)/4
      END IF

   END DO

   IF(JREPEAT>15)THEN
      NREPEAT=NREPEAT+1
      PRINT*,'BIG LOOPS =',JREPEAT,'REPEAT SETTING UP #',NREPEAT
      GOTO 20
   END IF

!----------------------------------------------------------
   PRINT*,'NUMBER OF BEADS =',NATOM
   PRINT*,'NUMBER OF STRANDS =',NPG
   PRINT*,'NUMBER OF BONDS =',NBOND
   PRINT*,'NUMBER OF GLYCAN BONDS =',NBONDGLY
   PRINT*,'NUMBER OF PEPTIDE BONDS =',NBONDPEP
   PRINT*,'NUMBER OF LOOPS =',NLOOP

!====================================================

!  ADD MATURE FEATURE:

!   ALLOCATE(MATURE(NATOM))
!   MATURE=1

!  ASSIGN DONOR AND ACCETOR PROPERTIES FOR PEPTIDES:

   ALLOCATE(DNOR(NATOM),ATOR(NATOM))

   DNOR=-1

   ATOR=1

   DNOR(BONDPEP(1,1:NBONDPEP))=0

   ATOR(BONDPEP(2,1:NBONDPEP))=0

!====================================================
! NOW WRITE OUT THE PSF:

50 CALL WRITEPSF(NATOM,NPG,PGLEN,PGTYP,NBOND,BOND)

   PRINT*,'========================================================'
   PRINT*,'NOW RELAXING THE SACCULUS'


!-----------------------------------
!------------------------------------------

!   IF(PARAL==1)THEN
!      call omp_set_num_threads(nthread)
!   END IF

!=====================================================================

!  Create trajactory file:

   open(10,file='PG0000.dcd',form='unformatted')

   coor='CORD'; nframe=10000; ifirst=0; nfreq=100; ntotal=100
   zeros5=0; jdelta=1; peroff=0; zeros7=0; twentyfour=24

   write(10)coor,nframe+1,ifirst,nfreq,ntotal,zeros5,jdelta,peroff,zeros7,twentyfour

   two=2; string1='welcome to the heaven!'; string2='go to hell!'

   write(10)two,string1,string2

   write(10)natom

!  The first frame is what created:

   allocate(xw(natom),yw(natom),zw(natom))
   xw=x; yw=y; zw=z
   write(10)(xw(i),i=1,natom)
   write(10)(yw(i),i=1,natom)
   write(10)(zw(i),i=1,natom)
   jprint=1

!  number of MD steps:
   NSTEP=1000*NATOM

   NRATIO=MAX(1,INT(SQRT(K_G/KTHETA)))

!  frequency of printing:
   NPRINT=1000

   NPEP=1
   NTHE=1

!  a small number:
   DELTA=0.005D0

!----------------

!  forces on beads:
      ALLOCATE(FX(NATOM),FY(NATOM),FZ(NATOM),F(NATOM))
!  forces due to glycan bonds:
      ALLOCATE(FXGLY(NATOM),FYGLY(NATOM),FZGLY(NATOM))
!  forces due to crosslinks:
      ALLOCATE(FXPEP(NATOM),FYPEP(NATOM),FZPEP(NATOM))
!  forces due to strand bending:
      ALLOCATE(FXTHETA(NATOM),FYTHETA(NATOM),FZTHETA(NATOM))
!  forces due to turgor pressure:
      ALLOCATE(FXPRES(NATOM),FYPRES(NATOM),FZPRES(NATOM))
!  to define new coordinates:

      IF(JOBID==1)THEN
         ALLOCATE(XNEW(NATOM),YNEW(NATOM),ZNEW(NATOM))
      END IF
!---------------

!  Initial calculation of energies:

!  ENGGLY, ENGTHETA are energies due to glycan bonds and bending

   ENGGLY=0.0D0; ENGTHETA=0.0D0
   FXGLY=0.0D0; FYGLY=0.0D0; FZGLY=0.0D0

   JFORCE=1

   CALL ESTRAND(ENGGLY,FXGLY,FYGLY,FZGLY,ENGTHETA,FXTHETA,FYTHETA,FZTHETA, &
                  NPG,PGID,PGLEN,PGTYP,X,Y,Z,L_G,K_G,THETA_0,KTHETA,JFORCE)

!-------------
!  ENGPEP is energy due to crosslinks

   ENGPEP=0.0D0;

   CALL EPEPBONDS(JFORCE,ENGPEP,FXPEP,FYPEP,FZPEP,X,Y,Z,NBONDPEP,BONDPEP,L_P,K_P)

!------------------------
!  And the last type is the contribution from pressure:


   ENGPRES=0.0D0; VOLM=0.0D0

   CALL EPRES(JFORCE,NATOM,ENGPRES,FXPRES,FYPRES,FZPRES,X,Y,Z,NLOOP,LOOP,LOOPLEN,PRES,VOLM,NTHREAD)

!----------------------------------------------------------------

   ENG=ENGPRES+ENGGLY+ENGPEP+ENGTHETA
   ENGBOND=ENGGLY+ENGPEP

   FX=FXGLY+FXPEP+FXTHETA+FXPRES
   FY=FYGLY+FYPEP+FYTHETA+FYPRES
   FZ=FZGLY+FZPEP+FZTHETA+FZPRES

   F=SQRT(FX**2+FY**2+FZ**2); FMAX=MAXVAL(F); FAVE=SUM(F/NATOM)

   N=0
   WRITE(*,11)'STEP',N,'AVE FORCE (pN) =',10*FAVE,'MAX FORCE (pN) =',10*MAXVAL(F)
   WRITE(*,12)'E_BONDS','E_ANGLES','E_PRESSURE','E_TOTAL (1e-20 J)','VOLUME (nm**3)'
   WRITE(*,13)ENGBOND,ENGTHETA,ENGPRES,ENG,VOLM

11     FORMAT(A4,1X,I8,4X,A18,1X,E10.3,4X,A18,1X,E10.3)
12     FORMAT(4X,A7,4X,A8,2X,A10,4X,A18,2X,A14)
13     FORMAT(F11.2,1X,F11.1,1X,F11.1,1X,F11.1,3X,F15.1)
14     FORMAT(A4,1X,I8,4X,A18,1X,E10.3,4X,A19,1X,E10.3)

!------------------------------

!--- CHECKING THE RANGE OF ENERGY CHANGES:

   TEM=0.0D0 ! This is to define temperature

   JFORCE=0

   DO I=1,500

      ENGGLYNEW=0.0D0; ENGTHETANEW=0.0D0

      DO N=1,NATOM

         CALL RANDOM_NUMBER(R)
         XNEW(N)=X(N)+DELTA*R

         CALL RANDOM_NUMBER(R)
         YNEW(N)=Y(N)+DELTA*R

         CALL RANDOM_NUMBER(R)
         ZNEW(N)=Z(N)+DELTA*R

      END DO

      FXGLY=0.0D0; FYGLY=0.0D0; FZGLY=0.0D0

      ENGGLYNEW=0.0D0; ENGPEPNEW=0.0D0; ENGTHETANEW=0.0D0; ENGPRESNEW=0.0D0; VOLMNEW=0.0D0

      CALL ESTRAND(ENGGLYNEW,FXGLY,FYGLY,FZGLY,ENGTHETANEW,FXTHETA,FYTHETA,FZTHETA, &
                     NPG,PGID,PGLEN,PGTYP,XNEW,YNEW,ZNEW,L_G,K_G,THETA_0,KTHETA,JFORCE)


      CALL EPEPBONDS(JFORCE,ENGPEPNEW,FXPEP,FYPEP,FZPEP,XNEW,YNEW,ZNEW,NBONDPEP,BONDPEP,L_P,K_P)

      CALL EPRES(JFORCE,NATOM,ENGPRESNEW,FXPRES,FYPRES,FZPRES,XNEW,YNEW,ZNEW, &
                 NLOOP,LOOP,LOOPLEN,PRES,VOLMNEW,NTHREAD)

      ENEW=ENGGLYNEW+ENGTHETANEW+ENGPEPNEW+ENGPRESNEW


      IF(ENEW>TEM)THEN
         TEM=ENEW
      END IF

   END DO

   TEM=(TEM-ENG)*5.0D0

!------------------------------------------------

!  Now run dynamics:

   CHECK=0

   STOPP=0

   DO N=1,NSTEP

      IF(MOD(N,10000)==0)THEN
         IF(NPEP<NRATIO)THEN
            NPEP=NPEP+1
            NTHE=NPEP
         END IF
      END IF

!  Here NPEP and NTHE are frequencies to calculate forces from crosslinks and bending

!  Now we have the forces on the atoms. Let move the atoms a small
!  distance along the directions of forces:

      F=SQRT(FX**2+FY**2+FZ**2); FMAX=MAXVAL(F); FAVE=SUM(F/NATOM)


      IF(CHECK==100)THEN
         DELTA=MIN(0.005D0,DELTA*2.0D0)
         CHECK=0
      END IF

      DO I=1,NATOM

         XNEW(I)=X(I)+DELTA*FX(I)
         YNEW(I)=Y(I)+DELTA*FY(I)
         ZNEW(I)=Z(I)+DELTA*FZ(I)

      END DO

!-------------------------------------------------------------------
!    Then repeat calculating the energy and forces of the system


      ENGGLYNEW=0.0D0; ENGTHETANEW=0.0D0
      FXGLY=0.0D0; FYGLY=0.0D0; FZGLY=0.0D0

      IF(MOD(N,NTHE)==0)THEN
         JFORCE=1
      ELSE
         JFORCE=0
      ENDIF

      CALL ESTRAND(ENGGLYNEW,FXGLY,FYGLY,FZGLY,ENGTHETANEW,FXTHETA,FYTHETA,FZTHETA, &
                     NPG,PGID,PGLEN,PGTYP,XNEW,YNEW,ZNEW,L_G,K_G,THETA_0,KTHETA,JFORCE)

!-------------

      IF(MOD(N,NPEP)==0)THEN
         JFORCE=1
      ELSE
         JFORCE=0
      ENDIF

      ENGPEPNEW=0.0D0;

      CALL EPEPBONDS(JFORCE,ENGPEPNEW,FXPEP,FYPEP,FZPEP,XNEW,YNEW,ZNEW,NBONDPEP,BONDPEP,L_P,K_P)

!------------------------
!       And the last type is the contribution from pressure:

      IF(MOD(N,50)==0)THEN
         JFORCE=1
      ELSE
         JFORCE=0
      ENDIF

      ENGPRESNEW=0.0D0; VOLMNEW=0.0D0

      CALL EPRES(JFORCE,NATOM,ENGPRESNEW,FXPRES,FYPRES,FZPRES,XNEW,YNEW,ZNEW, &
                 NLOOP,LOOP,LOOPLEN,PRES,VOLMNEW,NTHREAD)

!-------------------------------------------------------------------
      ENGNEW=ENGPRESNEW+ENGGLYNEW+ENGPEPNEW+ENGTHETANEW

!       If the energy decreases, move from X,Y,Z to XNEW,YNEW,ZNEW

      IF(ENG>ENGNEW)THEN

         X=XNEW; Y=YNEW; Z=ZNEW

!       And update energy and forces:

         ENG=ENGNEW; VOLM=VOLMNEW; ENGBOND=ENGGLYNEW+ENGPEPNEW
         ENGTHETA=ENGTHETANEW; ENGPRES=ENGPRESNEW

         FX=FXGLY+FXPEP+FXTHETA+FXPRES
         FY=FYGLY+FYPEP+FYTHETA+FYPRES
         FZ=FZGLY+FZPEP+FZTHETA+FZPRES

         CHECK=CHECK+1

!       If the energy does not decrease, move to XNEW,YNEW,ZNEW with a 
!       probability:

      ELSE

         CHECK=0

         PROB=EXP(-(ENGNEW-ENG)/TEM)

         DELTA=MAX(0.0001D0,DELTA*MAX(PROB,0.5D0))

         CALL RANDOM_NUMBER(R)


         IF(R<=PROB)THEN
            CHECK=CHECK+1

            X=XNEW; Y=YNEW; Z=ZNEW
            ENG=ENGNEW; VOLM=VOLMNEW
            ENGBOND=ENGGLYNEW+ENGPEPNEW; ENGTHETA=ENGTHETANEW; ENGPRES=ENGPRESNEW
            FX=FXGLY+FXPEP+FXTHETA+FXPRES
            FY=FYGLY+FYPEP+FYTHETA+FYPRES
            FZ=FZGLY+FZPEP+FZTHETA+FZPRES

         END IF

      END IF

      IF(MOD(N,NPRINT)==0)THEN

         F=SQRT(FX**2+FY**2+FZ**2); FMAX=MAXVAL(F); FAVE=SUM(F/NATOM)

         WRITE(*,*)  !'# OF DENIED STEPS = ',DENY
         WRITE(*,14)'STEP',N,'AVE FORCE (pN) =',10*FAVE,'AVE DISPLACE (nm) =',GAM*FAVE
         WRITE(*,12)'E_BONDS','E_ANGLES','E_PRESSURE','E_TOTAL (1e-20 J)','VOLUME (nm**3)'
         WRITE(*,13)ENGBOND,ENGTHETA,ENGPRES,ENG,VOLM

      END IF

!-------------------------------------------------------------------
      IF(MOD(N,100*jprint*jprint)==0)THEN
         xw=x; yw=y; zw=z
         write(10)(xw(i),i=1,natom)
         write(10)(yw(i),i=1,natom)
         write(10)(zw(i),i=1,natom)
         jprint=jprint+1
      END IF

!-------------------------

      IF(FAVE<0.001D0)THEN
         STOPP=1+STOPP
      END IF

      IF(STOPP>=100)THEN

         xw=x; yw=y; zw=z
         write(10)(xw(i),i=1,natom)
         write(10)(yw(i),i=1,natom)
         write(10)(zw(i),i=1,natom)

         DEALLOCATE(FX,FY,FZ,F)
         DEALLOCATE(FXGLY,FYGLY,FZGLY)
         DEALLOCATE(FXPEP,FYPEP,FZPEP)
         DEALLOCATE(FXTHETA,FYTHETA,FZTHETA)
         DEALLOCATE(FXPRES,FYPRES,FZPRES)
         DEALLOCATE(XNEW,YNEW,ZNEW)

         PRINT*,'SYSTEM REACHES EQUILIBRIUM AT STEP',N

         rewind(10)
         close(10)
         EXIT
      END IF

   END DO

!--------------------------------------------------
!  CALCULATE THE SACCULUS RADIUS:

   RADIUS=0.0D0
   NCOUNT=0

   DO N=1,NATOM

      IF(ABS(X(N))<20.0D0)THEN
         RADIUS=RADIUS+SQRT(Y(N)*Y(N)+Z(N)*Z(N))
         NCOUNT=NCOUNT+1
      END IF

   END DO

   RADIUS=RADIUS/NCOUNT

   IF(JOBID==3)THEN
      RADIUS=ORAD
   END IF

   CALL OUTPUTS(NATOM,X,Y,Z,DNOR,ATOR,PEPDIR,NPG,PGID,PGLEN,PGTYP,PGDIR, &
                NBOND,NBONDGLY,NBONDPEP,BOND,BONDGLY,BONDPEP,BONTYP, &
                NLOOP,LOOP,LOOPLEN,LOOPTYP,NATOMCAP,PGCAP1,PGCAP2,RADIUS)     


   WRITE(*,*)
   WRITE(*,*)'-------------------'
   print*,'Relaxation of the system is finished'
   WRITE(*,*)'****************************************************'

   OPEN(10,FILE='restart.inp')
   N=0
   WRITE(10,*)N,N
   WRITE(10,*)N
   CLOSE(10)

   CALL WHATTIME(TIMERUN)

   TIMERUN=TIMERUN-TIMESTART

   days=timerun/1440
   hours=(timerun-1440*days)/60
   mins=timerun-1440*days-60*hours

   write(*,91)'RUNNING TIME =',DAYS,'days : ',HOURS,'hours : ',MINS,'mins'


!  In case of Division job:

   IF(JOBID==3)THEN
      IF(LEFTCELL==1)THEN

         COMMAND='mkdir -p leftcell'
         CALL SYSTEM(COMMAND)

         COMMAND='mv PGstart.psf PG0000.dcd coor0000.inp restart.inp PG0000.pdb config0000.inp leftcell'
         CALL SYSTEM(COMMAND)

         COMMAND='cp -p *.info leftcell'

         CALL SYSTEM(COMMAND)

         WRITE(*,*)
         WRITE(*,*)'LEFT DAUGHTER CELL CREATED'
         WRITE(*,*)
         WRITE(*,*)'================================================'

!         DEALLOCATE(BOND)

      ELSE

         COMMAND='mkdir -p ritecell'
         CALL SYSTEM(COMMAND)

         COMMAND='mv PGstart.psf PG0000.dcd coor0000.inp restart.inp PG0000.pdb config0000.inp ritecell'
         CALL SYSTEM(COMMAND)

         COMMAND='cp -p *.info ritecell'

         CALL SYSTEM(COMMAND)

         WRITE(*,*)
         WRITE(*,*)'RIGHT DAUGHTER CELL CREATED'
         WRITE(*,*)
         WRITE(*,*)'================================================'

!         DEALLOCATE(BOND)

      END IF

   END IF


   IF(JOBID==1.OR.JOBID==3)THEN
      STOP
   END IF

!============================================================================
!============================================================================
!============================================================================
!============================================================================

72 IF(JOBID==2)THEN
!      PRINT*,'SELECT GROWTH MODE'
!      PRINT*,'0 for uncomplexed enzymes model'
!      PRINT*,'1 for multi-enzyme complex model'
!      PRINT*,'2 for addition of "cleaved crosslink capture"'
!      PRINT*,'3 for addition of "bend-induced termination"'
!      PRINT*,'4 for addition of "fixed transglycosylase orientation"'
!      PRINT*,'5 for addition of "crosslink before terminate"'
!      PRINT*,'6 for addition of "tail hydrolysis"'
!      PRINT*,'7 for addition of "paired transpeptidases"'
!      PRINT*,'8 for addition of "peptide maturation"'
!      PRINT*,'9 for addition of "crosslink-dependent processivity"'
!      PRINT*,'10 for addition of "first crosslink on the same side"'
!      PRINT*,'11 for addition of "hole-dependent processivity"'
!      PRINT*,'12 for addition of "paired-strand insertion"'
!      PRINT*,'13 for final model"'

!      READ(*,*)GMODE
!gmode=1
      PRINT*,'============================================'

      IF(GMODE<0.OR.GMODE>13)THEN
         PRINT*,'WRONG CHOICE. STOP NOW'
         STOP

      ELSEIF(GMODE==0)THEN
         PRINT*,'SACCULUS GROWTH IN UNCOMPLEXED ENZYMES MODE, GMODE=',GMODE
      ELSEIF(GMODE==1)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "MULTIENZYME COMPLEX" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==2)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "CLEAVED CROSSLINK CAPTURE" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==3)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "BEND-INDUCED TERMINATION" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==4)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "FIXED TRANSGLYCOSYLASE ORIENTATION" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==5)THEN
         PRINT*,'GROWTH AFTER ADDING "CROSSLINK BEFORE TERMINATE" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==6)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "TAIL HYDROLYSIS" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==7)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "PAIRED TRANSPEPTIDASES" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==8)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "PEPTIDE MATURATION" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==9)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "CROSSLINK-DEPENDENT PROCESSIVITY" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==10)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "FIRST CROSSLINK ON THE SAME SIDE" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==11)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "HOLE-DEPENDENT PROCESSIVITY" HYPOTHESIS, GMODE=',GMODE
      ELSEIF(GMODE==12)THEN
         PRINT*,'SACCULUS GROWTH AFTER ADDING "PAIRED-STRAND INSERTION MODE" HYPOTHESIS, GMODE=',GMODE

      ELSEIF(GMODE==13)THEN
         PRINT*,'SACCULUS GROWTH IN THE FINAL MODEL, GMODE=',GMODE

      END IF

   END IF

   OPEN(10,FILE='restart.inp')
   READ(10,*)NSTART,JFILE
   READ(10,*)LASTTIME
   CLOSE(10)

!  NSTART tells which configuration # to start
!  JFILE tells trajactory file # will start from JFILE+1
!  LASTTIME is the running time of the previously stopped job

!  Start with a number of complexes:
   NSYN=6

!  Tell how much would the sacculus grow:
   GROWTH=1.0D0
   IF(GMODE==0)THEN
      GROWTH=0.1D0
   ELSEIF(GMODE==1)THEN
      GROWTH=0.2D0
   ELSEIF(GMODE<4)THEN
      GROWTH=0.5D0
   ELSEIF(GMODE==4)THEN
      GROWTH=0.7D0
   END IF

!  Number of complexes increases with growth until:
   NSYNMAX=NSYN*(1.0D0+GROWTH)

!  Frequency of print out trajectory:
   IF(GMODE<2)THEN
      JPRINT1=500
      JPRINT2=1000000
   ELSEIF(GMODE<4)THEN
      JPRINT1=1000
      JPRINT2=2000000
   ELSEIF(GMODE<6)THEN
      JPRINT1=2000
      JPRINT2=5000000
   ELSE
      JPRINT1=10000
      JPRINT2=10000000
   END IF

!  Set parameters:
   CALL SETPARA(L_P,LSTART,LREP,LSWITCH,K_P,KSWITCH,MSWITCH,KGTASE,KTPASE,LTPASE,KSUR,KEDASE,LEDASE, &
        KGTTP,LGTTP,KGTED,LGTED,KSIDE,LSIDE,KWALL,KLEAD,KPAIR,LPAIR,DELTA,INVDELTA,BETA,INVMPBP,GMODE)



! -- SETUP GROWTH CONDITIONS:

!  Probability to activate a complex:
   PSTART=1.0D0/10000

!  To deactivate a complex:
   PSTOP=1.0D0/50000

!  To terminate transglycosylation:

   IF(GMODE<5)THEN
      PTERM=1.0D0/1000000
   ELSEIF(GMODE<9)THEN
      PTERM=1.0D0/200000 ! Because there is crosslink-before-terminate to prevent termination
   ELSEIF(GMODE<12)THEN
      PTERM=1.0D0/1000000  ! Because termination increases with crosslinkage
   ELSE
      PTERM=1.0D0/10000000 ! Crosslink-before-terminate is removed
   END IF


!  For a transglycosylase to translocate:
   P_TRANSFAST=1.0D0/20000
   IF(GMODE>=12)THEN
      P_TRANSFAST=1.0D0/30000
   END IF

   P_TRANSSLOW=1.0D0/2000000
   IF(GMODE<13)THEN
      P_TRANSSLOW=P_TRANSFAST
   END IF

!  For a transglycosylase to load a precursor:
   P_GLYIN=1.0D0/1000
   IF(GMODE>=12)THEN
      P_GLYIN=1.0D0/10000
   END IF

!  To release the precursor:
   P_GLYOUT=1.0D0/10000

!  For a transpeptidase to capture a donor peptide:

   DREACT=2.0D0
   DREACT2=DREACT*DREACT

!  To release the donor peptide:
   P_PEPOUT=1.0D0/10000

!  For an endopeptidase to release a captured peptide
   IF(GMODE<2)THEN
      P_EDHOLD=1.0D0/1000 ! Edase releases peptides right away
   ELSEIF(GMODE<13)THEN
      P_EDHOLD=1.0D0/10000000 ! Edase release peptides very slowly
   ELSE
      P_EDHOLD=1.0D0/10000  ! Edase release peptides quickly
   END IF

!  Probability of tail hydrolysis every 10000 steps:
   PLYT=1.0D0/10

!  Probability of maturation of peptides every 10000 steps:
   PMATURE=1.0D0/10

   CALL GETINFO(NSTART,NATOM,NPG,PGLENMAX,NBONDGLY,NBONDPEP,NLOOP,LOOPLENMAX)

!  Because the number of beads increase, need to allocate with an excess:
   NATOMMAX=NATOM+1000

!  Same with number of strands:
   NPGMAX=NPG+1000

!  Number of glycan bonds:
   NGLYMAX=NBONDGLY+1000

!  Number of peptide crosslinks:
   NPEPMAX=NBONDPEP+1000

!  Number of loops:
   NLOOPMAX=NLOOP+1000
   LOOPLENMAX=MAX(LOOPLENMAX*2,60)

   ALLOCATE(X(NATOMMAX),Y(NATOMMAX),Z(NATOMMAX),DNOR(NATOMMAX),ATOR(NATOMMAX),PEPDIR(NATOMMAX))
   ALLOCATE(PGID(PGLENMAX,NPGMAX),PGLEN(NPGMAX),PGTYP(NPGMAX),PGDIR(NPGMAX))
   ALLOCATE(BONDGLY(2,NGLYMAX),BONDPEP(2,NPEPMAX),BONTYP(NPEPMAX))
   ALLOCATE(LOOP(LOOPLENMAX,NLOOPMAX),LOOPLEN(NLOOPMAX),LOOPTYP(NLOOPMAX))

   ALLOCATE(SYNDIR(NSYNMAX),SYNTHESIS(2,NSYNMAX),GTLOAD(2,NSYNMAX),SYNPG(2,NSYNMAX),SYNLOOP(NSYNMAX))
   ALLOCATE(OLDSYNPG(2,NSYNMAX))
   SYNTHESIS=0; GTLOAD=0; SYNPG=0; SYNLOOP=0; OLDSYNPG=0

!  SYNDIR tells direction of a complex
!  SYNTHESIS tells if a complex is active
!  GTLOAD tells if a GTASE is loaded with a precursor
!  SYNPG tells the strand IDs being synthesized
!  SYNLOOP tells the loop ID the complex is in

   ALLOCATE(XGTASE(2,NSYNMAX),YGTASE(2,NSYNMAX),ZGTASE(2,NSYNMAX))
   ALLOCATE(XTPASE(3,NSYNMAX),YTPASE(3,NSYNMAX),ZTPASE(3,NSYNMAX))
   ALLOCATE(XEDASE(NSYNMAX),YEDASE(NSYNMAX),ZEDASE(NSYNMAX))
   ALLOCATE(XEDASEOLD(NSYNMAX),YEDASEOLD(NSYNMAX),ZEDASEOLD(NSYNMAX))

!  This assumes there are 2 GTASES (transglycosylases), 3 TPASES (transpeptidases) and 1 EDASE (endopeptidase)

   ALLOCATE(GLYTIP(2,NSYNMAX),GLYSEC(2,NSYNMAX),TPPEP(2,3,NSYNMAX),EDPEP(2,NSYNMAX),SYNRAD(NSYNMAX))
   GLYTIP=0; GLYSEC=0; TPPEP=0; EDPEP=0

!  This is to define local radii at GTASE and TPASE when EDASE is not in the complex:
   ALLOCATE(GTRAD(NSYNMAX))

!  GLYTIP tells bead ID of the strand tip which is being synthesized
!  GLYSEC tells the second bead, next to the tip
!  TPPEP tells bead ID of peptides held by TPASE
!  EDPEP tells bead ID of peptides held by EDASE
!  SYNRAD tells the local radius of the sacculus where the complex is at

!  To tell if GTASE is already translocate:
   ALLOCATE(GTATRANS(2,NSYNMAX))
   GTATRANS=0

!  to tell whose turn it is to wait before adding a new PG in the pair mode:
   ALLOCATE(GTAWAIT(NSYNMAX))
   GTAWAIT=0

!  to mark crosslinking between paired strands:
   ALLOCATE(JPAIR(NSYNMAX))

   DO NS=1,NSYNMAX
      CALL RANDOM_NUMBER(R)
      IF(R>0.5D0)THEN
         JPAIR(NS)=1
      ELSE
         JPAIR(NS)=-1
      END IF
   END DO

!  To force deactivate a complex:
   ALLOCATE(JDEACT(NSYNMAX))
   JDEACT=0

!  To re-initiate a complex:
   ALLOCATE(JREINI(NSYNMAX))
   JREINI=0

!  Signal for crosslinking:
   ALLOCATE(SIGCROSS(3,NSYNMAX))
   SIGCROSS=0

!  Signal for cleaving a crosslink:
   ALLOCATE(SIGCLEAVE(NSYNMAX))
   SIGCLEAVE=0

!  Tells if EDASE is locked in a crosslink, ready to cleave:
   ALLOCATE(EDLOCKIN(NSYNMAX))
   EDLOCKIN=0

!  Tells if EDASE captures a crosslink:
   ALLOCATE(EDCAP(NSYNMAX))
   EDCAP=0

!  For signaling crosslinks which are part of trimeric crosslinks:
   ALLOCATE(SIGBOND(2,200))
   SIGBOND=0

!  THIS IS USED TO TRACK PG UNITS AND BONDS IN ASSOCIATION WITH ENDOPEPTIDASE:
   ALLOCATE(EDHOLD(3,NSYNMAX))
   EDHOLD=0

!  TRACK CROSS-LINKAGE OF STRAND ASSOCIATED WITH THE COMPLEXES:
   ALLOCATE(CRLKAGE(2,NSYNMAX))
   CRLKAGE=0

!  TO TELL IF A COMPLEX IS TRAPPED:

   ALLOCATE(JTRAP(NSYNMAX))
   JTRAP=0
   NTRAP=1000000
   IF(GMODE>=12)THEN
      NTRAP=5000000
   END IF
!------------------------------------------------
!  Inputs from the original configuration:

   CALL PG_INPUT(GMODE,NSTART,NATOM,NATOMDEL,OLDNATOMDEL,NATOMCAP,NATOMSTART,DNOR,ATOR,PEPDIR,X,Y,Z,PGCAP1,PGCAP2,ORAD, &
                 NPG,NPGSTART,PGID,PGLEN,PGTYP,PGDIR,NBONDGLY,NBONDGLYSTART,NBONDPEP,NBONDDEL, &
                 BONDGLY,BONDPEP,BONTYP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                 NSYN,SYNDIR,SYNTHESIS,GTLOAD,SYNPG,SYNLOOP,GLYTIP,GLYSEC,TPPEP,EDPEP, &
                 GTATRANS,JDEACT,SIGCROSS,SIGCLEAVE,EDLOCKIN,EDCAP,EDHOLD,CRLKAGE,SYNRATIO, &
                 XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD)

!  NATOMCAP is the number of atoms on two caps
!  NATOMDEL is the number of atoms removed from the sacculus by lytic transglycosylases
!  PGCAP1, PGCAP2 are indices of 2 residues (hoops) connecting 2 caps with the cylinder
!  NATOMSTART is the number of atoms at NSTART=0, this is used for visualization purpose
!  NPGSTART is the number of residues (strands) at NSTART=0

!------------------------------------------------

!  For visualization, write PSF:
   CALL VISUALPSF(GROWTH,NSTART,NATOMSTART,NATOMCAP,NBONDGLYSTART,NPG,PGID,PGLEN,PGTYP, &
                  NSYNMAX,NTOTAL,NGLYNEW,NPEPTOTAL,BOND,BONDGLY)

!  NTOTAL is the max number of atoms in visualization including glycan atoms, peptide atoms,
!  NPEPTOTAL is the max number of peptides in visualization

!  Start a trajactory file:
   CALL RANDOM_NUMBER(R)
   JUNIT=80*R+11
   JFILE=JFILE+1

   CALL DCDHEADER(JUNIT,JFILE,NTOTAL)

   NFRAME=0

!------------------------------------------------

!  FORCES ON BEADS AND PBPS:

   ALLOCATE(FX(NATOMMAX),FY(NATOMMAX),FZ(NATOMMAX),FAMAG(NATOMMAX))
   ALLOCATE(FXGLY(NATOMMAX),FYGLY(NATOMMAX),FZGLY(NATOMMAX))
   ALLOCATE(FXPEP(NATOMMAX),FYPEP(NATOMMAX),FZPEP(NATOMMAX))
   ALLOCATE(FXTHETA(NATOMMAX),FYTHETA(NATOMMAX),FZTHETA(NATOMMAX))
   ALLOCATE(FXPRES(NATOMMAX),FYPRES(NATOMMAX),FZPRES(NATOMMAX))
   ALLOCATE(FXSYN(NATOMMAX),FYSYN(NATOMMAX),FZSYN(NATOMMAX))
   ALLOCATE(FYSUR(NATOMMAX),FZSUR(NATOMMAX),ATOMRAD(NATOMMAX))

   ALLOCATE(FXGTASE(2,NSYNMAX),FYGTASE(2,NSYNMAX),FZGTASE(2,NSYNMAX))
   ALLOCATE(FXTPASE(3,NSYNMAX),FYTPASE(3,NSYNMAX),FZTPASE(3,NSYNMAX))
   ALLOCATE(FXEDASE(NSYNMAX),FYEDASE(NSYNMAX),FZEDASE(NSYNMAX),FMAGPBP(NSYNMAX))
   ALLOCATE(FXCAPGT(NSYNMAX),FXCAPEDASE(NSYNMAX))

!  To take into account diffusion, use random forces:

   NRAND=1000000

   JRESET=NRAND/5

   NSYNREP=MIN(NSYN+2,NSYNMAX)
   ALLOCATE(RXGTASE(2*NSYNREP,NRAND),RYGTASE(2*NSYNREP,NRAND),RZGTASE(2*NSYNREP,NRAND))
   ALLOCATE(RXTPASE(3*NSYNREP,NRAND),RYTPASE(3*NSYNREP,NRAND),RZTPASE(3*NSYNREP,NRAND))
   ALLOCATE(RXEDASE(NSYNREP,NRAND),RYEDASE(NSYNREP,NRAND),RZEDASE(NSYNREP,NRAND))

!  SETS OF RANDOM NUMBERS FOR CONVENIENCE:

   ALLOCATE(RANDS(NRAND,NSYNMAX),JRANDS(NSYNMAX))
   JRANDS=NRAND

!------------------------------------------------

!  MISCELLANEOUS:
!  LGTASE denotes the relaxed distance from the end of new strand to GTASE
!  XLEAD,YLEAD,ZLEAD denote the prefered orientation of GTASE
!  UX,UY,UZ are components of the axis of sacculus

   ALLOCATE(LGTASE(2,NSYNMAX),XLEAD(NSYNMAX),YLEAD(NSYNMAX),ZLEAD(NSYNMAX))
   ALLOCATE(XCEN(NSYNMAX),YCEN(NSYNMAX),ZCEN(NSYNMAX),UX(NSYNMAX),UY(NSYNMAX),UZ(NSYNMAX))

   LGTASE=0.5D0

!  NEIGHBOR BONDS FOR TOPOLOGICAL CONSTRAINT ON GTASE:
!  For each complex, count number of neighbor bonds to GTASE as NEINUM. NEIBOND denote the 
!  atoms on those bonds

   ALLOCATE(NEINUM(NSYNMAX),NEIBOND(2,100,NSYNMAX))
   ALLOCATE(GLYNUM(NSYNMAX),NEIGLY(2,50,NSYNMAX))

!  Signaling bonds in the pair insertion mode:
   CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

!=============================================================================

!  DYNAMICS OF GROWTH:

   DO JSTEP=1,1000000000

!  BEGINNING EACH STEP, CHECK AND UPDATE CONDITIONS

      JFORCE=0

      IF(MOD(JSTEP-1,10)==0)THEN
         JFORCE=1
      END IF

!     UPDATE THE SYSTEM:

      IF(MOD(JSTEP-1,NRAND)==0)THEN


!        ADDING A NEW COMPLEX:

         IF((NATOM-NATOMCAP/2)/SYNRATIO>NSYN.AND.NSYN<NSYNMAX)THEN

            XCAP1=SUM(X(PGID(1:PGLEN(PGCAP1),PGCAP1)))/PGLEN(PGCAP1)
            XCAP2=SUM(X(PGID(1:PGLEN(PGCAP2),PGCAP2)))/PGLEN(PGCAP2)

            NSYN=NSYN+1

            CALL ADDSYNCOM(NATOM,X,Y,Z,XCAP1,XCAP2,NSYN,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE, &
                 XEDASE,YEDASE,ZEDASE,SYNDIR,GMODE,JREINI)


            NSYNREP=MIN(NSYN+2,NSYNMAX)

            DEALLOCATE(RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE)

            ALLOCATE(RXGTASE(2*NSYNREP,NRAND),RYGTASE(2*NSYNREP,NRAND),RZGTASE(2*NSYNREP,NRAND))
            ALLOCATE(RXTPASE(3*NSYNREP,NRAND),RYTPASE(3*NSYNREP,NRAND),RZTPASE(3*NSYNREP,NRAND))
            ALLOCATE(RXEDASE(NSYNREP,NRAND),RYEDASE(NSYNREP,NRAND),RZEDASE(NSYNREP,NRAND))

         END IF

!        UPDATE RANDOM FORCES:

         CALL SETRAND(RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE,JRFORCE,NRAND,NSYNREP,PI)


         CALL REFRAME(NATOM,X,Y,Z,NSYN,XEDASE,YEDASE,ZEDASE,XCEN,YCEN,ZCEN,UX,UY,UZ)

      END IF


      IF(MOD(JSTEP-1,10000)==0)THEN

         XCAP1=SUM(X(PGID(1:PGLEN(PGCAP1),PGCAP1)))/PGLEN(PGCAP1)
         XCAP2=SUM(X(PGID(1:PGLEN(PGCAP2),PGCAP2)))/PGLEN(PGCAP2)

         CALL LEADIR(NSYN,SYNDIR,XGTASE,YGTASE,ZGTASE,XCEN,YCEN,ZCEN,UX,UY,UZ,XLEAD,YLEAD,ZLEAD)


         IF(GMODE>5)THEN
            CALL LYTGTASE(NATOM,NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,PLYT,NATOMDEL,NSYN,SYNPG,SYNLOOP,LOOPLEN,LOOP)

            NDEL=NATOMDEL-OLDNATOMDEL

            IF(NDEL>0)THEN

               CALL UPDATESYS(NATOM,NATOMSTART,NDEL,DNOR,ATOR,PEPDIR,ATOMRAD,X,Y,Z,NPG,NPGSTART,PGID, &
                      PGLEN,PGTYP,PGDIR,NSYN,SYNPG,OLDSYNPG,GLYTIP,GLYSEC,TPPEP,EDPEP,EDHOLD, &
                        NLOOP,LOOP,LOOPLEN,LOOPTYP,NBONDGLY,NBONDGLYSTART,NBONDPEP,BONDGLY,BONDPEP,BONTYP)

               OLDNATOMDEL=NATOMDEL

               CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

            END IF

         END IF

!        MATURATION OF PEPTIDES:

         IF(GMODE>7.AND.GMODE<13)THEN
            CALL MATURATION(NATOM,DNOR,SYNPG,PGID,PGLEN,SYNLOOP,LOOP,LOOPLEN,NSYN,PMATURE)
         END IF


         DO NS=1,NSYN

            IF(JREINI(NS)==1)THEN

               CALL ADDSYNCOM(NATOM,X,Y,Z,XCAP1,XCAP2,NS,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE, &
                    XEDASE,YEDASE,ZEDASE,SYNDIR,GMODE,JREINI)

            END IF

            CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                      NEIBOND,NEINUM,NEIGLY,GLYNUM)

         END DO

      END IF

!-------------------------------------------------------------


!     NOW WE CHECK AND UPDATE CONDITIONS OF ENZYMES

!========================================== PBP1 TRANSGLYCOSYLATION RULES:


      DO NS=1,NSYN

         IF(JDEACT(NS)==1)THEN
            CYCLE
         ENDIF

!        reset random numbers:

         IF(JRANDS(NS)>NRAND-100)THEN
            CALL RANDOM_NUMBER(RANDS(:,NS))
            JRANDS(NS)=0
         ENDIF

!        activating a complex:

         IF((SYNTHESIS(1,NS)==0.OR.SYNTHESIS(2,NS)==0).AND.JREINI(NS)==0)THEN

            CALL ACTIVATE(NS,SYNTHESIS,SYNLOOP,XGTASE,YGTASE,ZGTASE,X,Y,Z,NEIGLY,GLYNUM, &
                          NLOOP,LOOP,LOOPLEN,LOOPTYP,PSTART,GMODE)

         END IF


!        LOADING PRECURSORS:

         ILOOP=SYNLOOP(NS)

         IF(SYNTHESIS(1,NS)==1.AND.GTLOAD(1,NS)==0.AND.GTATRANS(1,NS)==1)THEN


!           HOLE-DEPENDENT PROCESSIVITY:
            AREA=1.0D0

            IF(GMODE>10)THEN
               AREA=MAX(1.0D0,0.01D0*LOOPLEN(ILOOP)**2)
            END IF

            JRANDS(NS)=JRANDS(NS)+1

            IF(P_GLYIN*AREA>RANDS(JRANDS(NS),NS))THEN
               GTLOAD(1,NS)=1

            END IF

         END IF

         IF(SYNTHESIS(2,NS)==1.AND.GTLOAD(2,NS)==0.AND.GTATRANS(2,NS)==1.AND.GMODE>=12)THEN

            AREA=MAX(1.0D0,0.01D0*LOOPLEN(ILOOP)**2)

            JRANDS(NS)=JRANDS(NS)+1

            IF(P_GLYIN*AREA>RANDS(JRANDS(NS),NS))THEN
               GTLOAD(2,NS)=1

            END IF

         END IF

!        UNLOADING PRECURSORS:

         IF(GTLOAD(1,NS)==1)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_GLYOUT>RANDS(JRANDS(NS),NS))THEN
               GTLOAD(1,NS)=0
            END IF

         END IF

         IF(GTLOAD(2,NS)==1)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_GLYOUT>RANDS(JRANDS(NS),NS))THEN
               GTLOAD(2,NS)=0
            END IF

         END IF

!        TO PREVENT UNCROSSLINKED STRAND FROM ELONGATING FOREVER
         IF(GTLOAD(1,NS)==1.AND.SYNPG(1,NS)>0.AND.CRLKAGE(1,NS)==0)THEN

            IPG=SYNPG(1,NS)

            IF(PGLEN(IPG)>25)THEN

               GTLOAD(1,NS)=0
               JDEACT(NS)=1


            END IF

         END IF

         IF(GTLOAD(2,NS)==1.AND.SYNPG(2,NS)>0.AND.CRLKAGE(2,NS)==0)THEN

            IPG=SYNPG(2,NS)

            IF(PGLEN(IPG)>25)THEN


               GTLOAD(2,NS)=0
               JDEACT(NS)=1


            END IF

         END IF

!--------

!        TRANSLOCATION OF GTASES ON STRANDS:

         IF(GTATRANS(1,NS)==0.AND.SYNPG(1,NS)>0.AND.GTAWAIT(NS)/=1)THEN

            IPG=SYNPG(1,NS)
            NTIP=PGID(PGLEN(IPG),IPG)

            IF(CRLKAGE(1,NS)==0.OR.DNOR(NTIP)==0.OR.ATOR(NTIP)==0)THEN

               JRANDS(NS)=JRANDS(NS)+1

               IF(P_TRANSFAST>RANDS(JRANDS(NS),NS))THEN
                  GTATRANS(1,NS)=1
               END IF

            ELSE

               JRANDS(NS)=JRANDS(NS)+1

               IF(P_TRANSSLOW>RANDS(JRANDS(NS),NS))THEN
                  GTATRANS(1,NS)=1
               END IF

            END IF

         ELSEIF(GTATRANS(1,NS)==0.AND.SYNPG(1,NS)==0)THEN

            GTATRANS(1,NS)=1

         END IF

         IF(GTATRANS(2,NS)==0.AND.SYNPG(2,NS)>0.AND.GTAWAIT(NS)/=2)THEN

            IPG=SYNPG(2,NS)
            NTIP=PGID(PGLEN(IPG),IPG)

            IF(CRLKAGE(2,NS)==0.OR.DNOR(NTIP)==0.OR.ATOR(NTIP)==0)THEN

               JRANDS(NS)=JRANDS(NS)+1

               IF(P_TRANSFAST>RANDS(JRANDS(NS),NS))THEN
                  GTATRANS(2,NS)=1
               END IF

            ELSE

               JRANDS(NS)=JRANDS(NS)+1

               IF(P_TRANSSLOW>RANDS(JRANDS(NS),NS))THEN
                  GTATRANS(2,NS)=1
               END IF

            END IF

         ELSEIF(GTATRANS(2,NS)==0.AND.SYNPG(2,NS)==0)THEN

            GTATRANS(2,NS)=1

         END IF


!        RELAXED DISTANCE FROM GTASE TO TIP

         LGTASE(1:2,NS)=0.5D0

         IF(SYNPG(1,NS)>0.AND.GTLOAD(1,NS)==1)THEN
            LGTASE(1,NS)=0.5D0+L_G
         END IF

         IF(SYNPG(2,NS)>0.AND.GTLOAD(2,NS)==1)THEN
            LGTASE(2,NS)=0.5D0+L_G
         END IF

!        ELONGATION IS DELAYED FOR SHORT UNCROSSLINKED STRANDS IN PAIRED STRAND MODE:

         IF(GMODE>=12)THEN

            IF(CRLKAGE(1,NS)==1)THEN

               IPG=SYNPG(1,NS)

               NTIP=PGID(PGLEN(IPG),IPG)

               IF(DNOR(NTIP)/=0.AND.ATOR(NTIP)/=0)THEN
                  LGTASE(1,NS)=0.5D0
               END IF

            END IF

            IF(CRLKAGE(2,NS)==1)THEN

               IPG=SYNPG(2,NS)

               NTIP=PGID(PGLEN(IPG),IPG)

               IF(DNOR(NTIP)/=0.AND.ATOR(NTIP)/=0)THEN
                  LGTASE(2,NS)=0.5D0
               END IF

            END IF

         END IF

!        ELONGATION IS DELAYED AS EDASE IS ENGAGING A CROSSLINK:

         IF(EDLOCKIN(NS)==1.AND.GMODE==13)THEN
            LGTASE(1:2,NS)=0.5D0
         END IF



!-------------------------------------------


!        INITIATE NEW STRANDS:

         IF(GTLOAD(1,NS)==1.AND.SYNPG(1,NS)==0)THEN

            JS=1

            CALL FIRSTPG(JS,NS,NPG,PGID,PGLEN,PGTYP,PGDIR,NATOM,DNOR,ATOR,PEPDIR,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                 SYNPG,GTLOAD,SYNDIR,GLYTIP,JFORCE,GTATRANS,ATOMRAD,SYNRAD,XLEAD,YLEAD,ZLEAD,GMODE,JPAIR,GTAWAIT)

            JTRAP(NS)=JSTEP

         END IF

         IF(GTLOAD(2,NS)==1.AND.SYNPG(2,NS)==0)THEN

            JS=2

            CALL FIRSTPG(JS,NS,NPG,PGID,PGLEN,PGTYP,PGDIR,NATOM,DNOR,ATOR,PEPDIR,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                 SYNPG,GTLOAD,SYNDIR,GLYTIP,JFORCE,GTATRANS,ATOMRAD,SYNRAD,XLEAD,YLEAD,ZLEAD,GMODE,JPAIR,GTAWAIT)

            JTRAP(NS)=JSTEP

         END IF


!        ELONGATING STRANDS:

         IF(GTLOAD(1,NS)==1.AND.SYNPG(1,NS)/=0.AND.MOD(JSTEP,100)==0)THEN

            JS=1

            CALL ELONGATE(JS,NS,SYNPG,PGID,PGLEN,NBONDGLY,BONDGLY,X,Y,Z,XGTASE,YGTASE,ZGTASE,L_G, &
                 GLYTIP,GLYSEC,GTLOAD,NATOM,DNOR,ATOR,PEPDIR,JFORCE,GTATRANS,ATOMRAD,SYNRAD,GTAWAIT,GMODE,JPAIR)

            IF(GTLOAD(1,NS)==0)THEN

               CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                      NEIBOND,NEINUM,NEIGLY,GLYNUM)

               JTRAP(NS)=JSTEP

               IF(PGLEN(SYNPG(1,NS))>75)THEN
                  JDEACT(NS)=1
               END IF

!              adding bend-induced termination hypothesis:
               CALL BENDING(NS,SYNPG,PGID,PGLEN,CRLKAGE,DNOR,ATOR,X,Y,Z,XGTASE,YGTASE,ZGTASE,BETA,PI,JBEND)

               IF(JBEND==1)THEN
                  JDEACT(NS)=1
               END IF

            END IF

         END IF

         IF(GTLOAD(2,NS)==1.AND.SYNPG(2,NS)/=0.AND.MOD(JSTEP,100)==0)THEN

            JS=2

            CALL ELONGATE(JS,NS,SYNPG,PGID,PGLEN,NBONDGLY,BONDGLY,X,Y,Z,XGTASE,YGTASE,ZGTASE,L_G, &
                 GLYTIP,GLYSEC,GTLOAD,NATOM,DNOR,ATOR,PEPDIR,JFORCE,GTATRANS,ATOMRAD,SYNRAD,GTAWAIT,GMODE,JPAIR)

            IF(GTLOAD(2,NS)==0)THEN

               CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                      NEIBOND,NEINUM,NEIGLY,GLYNUM)

               JTRAP(NS)=JSTEP

               IF(PGLEN(SYNPG(2,NS))>75)THEN
                  JDEACT(NS)=1
               END IF

               CALL BENDING(NS,SYNPG,PGID,PGLEN,CRLKAGE,DNOR,ATOR,X,Y,Z,XGTASE,YGTASE,ZGTASE,BETA,PI,JBEND)

               IF(JBEND==1)THEN
                  JDEACT(NS)=1
               END IF

            END IF

         END IF



!==================  TRANSPEPTIDATION RULES:



         IF(SYNTHESIS(1,NS)==0.OR.SYNTHESIS(2,NS)==0)THEN
            CYCLE
         END IF

!        THE FIRST TPASE LOADING DONOR SITE:

         IF(SYNPG(1,NS)/=0.AND.TPPEP(1,1,NS)==0.AND.(CRLKAGE(1,NS)+CRLKAGE(2,NS)>0.OR.GMODE<12))THEN

            IPG=SYNPG(1,NS)
            NTIP=PGID(PGLEN(IPG),IPG)

            IF(DNOR(NTIP)==1.AND.(PEPDIR(NTIP)==SYNDIR(NS).OR.GMODE<7))THEN

               DX=XTPASE(1,NS)-X(NTIP)
               DY=YTPASE(1,NS)-Y(NTIP)
               DZ=ZTPASE(1,NS)-Z(NTIP)

               DIST=SQRT(DX**2+DY**2+DZ**2)

               IF(DIST<DREACT)THEN

                  P_PEPIN=(1.0-DIST/DREACT)**2

                  JRANDS(NS)=JRANDS(NS)+1

                  IF(P_PEPIN>RANDS(JRANDS(NS),NS))THEN
                     TPPEP(1,1,NS)=NTIP
                  END IF

               END IF

            END IF

!  from the mode with two TPASES, but not the final model, crosslinking can happen to the first two beads:

            BACKCRLK=1

            IF(PGLEN(IPG)==1.OR.GMODE<7.OR.GMODE>=12)THEN
               BACKCRLK=0
            END IF

            IF(PGLEN(IPG)>1)THEN

               NPREV=PGID(PGLEN(IPG)-1,IPG)

               IF(DNOR(NPREV)/=1.OR.PEPDIR(NPREV)==-SYNDIR(NS))THEN
                  BACKCRLK=0
               END IF

            END IF


            IF(BACKCRLK==1)THEN

               DX=XTPASE(1,NS)-X(NPREV)
               DY=YTPASE(1,NS)-Y(NPREV)
               DZ=ZTPASE(1,NS)-Z(NPREV)

               DIST=SQRT(DX**2+DY**2+DZ**2)

               IF(DIST<DREACT)THEN

                  P_PEPIN=(1.0-DIST/DREACT)**2

                  JRANDS(NS)=JRANDS(NS)+1

                  IF(P_PEPIN>RANDS(JRANDS(NS),NS))THEN
                     TPPEP(1,1,NS)=NPREV
                  END IF

               END IF

            END IF

!        peptide could be released:

         ELSEIF(TPPEP(1,1,NS)>0)THEN

            NA=TPPEP(1,1,NS)

            DX=XTPASE(1,NS)-X(NA)
            DY=YTPASE(1,NS)-Y(NA)
            DZ=ZTPASE(1,NS)-Z(NA)

            IF(DX*DX+DY*DY+DZ*DZ>DREACT2)THEN
               TPPEP(1:2,1,NS)=0
            END IF

         END IF

!        THE SECOND TPASE LOADING DONOR PEPTIDE, ONLY IF GMODE>6:

         IF(SYNPG(1,NS)/=0.AND.TPPEP(1,2,NS)==0.AND.GMODE>6)THEN

            IPG=SYNPG(1,NS)

            NTIP=PGID(PGLEN(IPG),IPG)

!           Note that, "first crosslink on the same side" applied if GMODE = 10,11:

            IF(DNOR(NTIP)==1.AND.PEPDIR(NTIP)==-SYNDIR(NS).AND.(CRLKAGE(1,NS)>0.OR.GMODE<10.OR.GMODE>=12))THEN

               DX=XTPASE(2,NS)-X(NTIP)
               DY=YTPASE(2,NS)-Y(NTIP)
               DZ=ZTPASE(2,NS)-Z(NTIP)

               DIST=SQRT(DX**2+DY**2+DZ**2)

               IF(DIST<DREACT)THEN

                  P_PEPIN=(1.0-DIST/DREACT)**2

                  JRANDS(NS)=JRANDS(NS)+1

                  IF(P_PEPIN>RANDS(JRANDS(NS),NS))THEN
                     TPPEP(1,2,NS)=NTIP
                  END IF

               END IF

            END IF

!  from the mode with two TPASES, but not the final model, crosslinking can happen to the first two beads:

            BACKCRLK=1

            IF(PGLEN(IPG)==1.OR.GMODE<7.OR.GMODE>=12)THEN
               BACKCRLK=0
            END IF

            IF(PGLEN(IPG)>1)THEN

               NPREV=PGID(PGLEN(IPG)-1,IPG)

               IF(DNOR(NPREV)/=1.OR.PEPDIR(NPREV)==SYNDIR(NS))THEN
                  BACKCRLK=0
               END IF

            END IF

!           Note that, "first crosslink on the same side" applied if GMODE = 10,11:

            IF((GMODE==10.OR.GMODE==11).AND.CRLKAGE(1,NS)<2)THEN
               BACKCRLK=0
            END IF

               
            IF(BACKCRLK==1)THEN

               DX=XTPASE(2,NS)-X(NPREV)
               DY=YTPASE(2,NS)-Y(NPREV)
               DZ=ZTPASE(2,NS)-Z(NPREV)

               DIST=SQRT(DX**2+DY**2+DZ**2)

               IF(DIST<DREACT)THEN

                  P_PEPIN=(1.0-DIST/DREACT)**2

                  JRANDS(NS)=JRANDS(NS)+1

                  IF(P_PEPIN>RANDS(JRANDS(NS),NS))THEN
                     TPPEP(1,2,NS)=NPREV
                  END IF

               END IF

            END IF

!        peptide could be released:

         ELSEIF(TPPEP(1,2,NS)>0)THEN

            NA=TPPEP(1,2,NS)

            DX=XTPASE(2,NS)-X(NA)
            DY=YTPASE(2,NS)-Y(NA)
            DZ=ZTPASE(2,NS)-Z(NA)

            IF(DX*DX+DY*DY+DZ*DZ>DREACT2)THEN
               TPPEP(1:2,2,NS)=0
            END IF

         END IF

!        THE THIRD TPASE LOADING DONOR SITE:

         IF(SYNPG(2,NS)/=0.AND.TPPEP(1,3,NS)==0.AND.CRLKAGE(1,NS)+CRLKAGE(2,NS)>0)THEN

            IPG=SYNPG(2,NS)
            NTIP=PGID(PGLEN(IPG),IPG)

            IF(DNOR(NTIP)==1.AND.PEPDIR(NTIP)==-SYNDIR(NS))THEN

               DX=XTPASE(3,NS)-X(NTIP)
               DY=YTPASE(3,NS)-Y(NTIP)
               DZ=ZTPASE(3,NS)-Z(NTIP)

               DIST=SQRT(DX**2+DY**2+DZ**2)

               IF(DIST<DREACT)THEN

                  P_PEPIN=(1.0-DIST/DREACT)**2

                  JRANDS(NS)=JRANDS(NS)+1

                  IF(P_PEPIN>RANDS(JRANDS(NS),NS))THEN
                     TPPEP(1,3,NS)=NTIP
                  END IF

               END IF

            END IF

!        peptide could be released:

         ELSEIF(TPPEP(1,3,NS)>0)THEN

            NA=TPPEP(1,3,NS)

            DX=XTPASE(3,NS)-X(NA)
            DY=YTPASE(3,NS)-Y(NA)
            DZ=ZTPASE(3,NS)-Z(NA)

            IF(DX*DX+DY*DY+DZ*DZ>DREACT2)THEN
               TPPEP(1:2,3,NS)=0
            END IF

         END IF

!------  THE FIRST TPASE WOULD RELEASE PEPTIDE HOOKS:

         IF(TPPEP(1,1,NS)/=0.AND.TPPEP(2,1,NS)==0)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_PEPOUT>RANDS(JRANDS(NS),NS))THEN
               TPPEP(1,1,NS)=0
            END IF
         END IF

         IF(TPPEP(2,1,NS)/=0)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_PEPOUT>RANDS(JRANDS(NS),NS))THEN
               TPPEP(2,1,NS)=0
            END IF
         END IF

!------  THE SECOND TPASE WOULD RELEASE PEPTIDE HOOKS:

         IF(TPPEP(1,2,NS)/=0.AND.TPPEP(2,2,NS)==0)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_PEPOUT>RANDS(JRANDS(NS),NS))THEN
               TPPEP(1,2,NS)=0
            END IF
         END IF

         IF(TPPEP(2,2,NS)/=0)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_PEPOUT>RANDS(JRANDS(NS),NS))THEN
               TPPEP(2,2,NS)=0
            END IF
         END IF

!------  THE THIRD TPASE WOULD RELEASE PEPTIDE HOOKS:

         IF(TPPEP(1,3,NS)/=0.AND.TPPEP(2,3,NS)==0)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_PEPOUT>RANDS(JRANDS(NS),NS))THEN
               TPPEP(1,3,NS)=0
            END IF
         END IF

         IF(TPPEP(2,3,NS)/=0)THEN
            JRANDS(NS)=JRANDS(NS)+1

            IF(P_PEPOUT>RANDS(JRANDS(NS),NS))THEN
               TPPEP(2,3,NS)=0
            END IF
         END IF


!        THE FIRST TRANSPEPTIDASE CROSSLINKING:

         IF(TPPEP(1,1,NS)>0.AND.TPPEP(2,1,NS)==0)THEN

            JS=1

            CALL CRLK2SAC(JS,NS,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,BONDPEP,BONTYP, &
                 PEPDIR,X,Y,Z,SYNLOOP,LOOP,LOOPLEN,DNOR,ATOR,EDHOLD,EDLOCKIN,GMODE,DREACT)

         END IF

!        THE SECOND TRANSPEPTIDASE CROSSLINKING:

         IF(TPPEP(1,2,NS)>0.AND.TPPEP(2,2,NS)==0)THEN

            IF(GMODE<12)THEN

               JS=2

               CALL CRLK2SAC(JS,NS,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,BONDPEP,BONTYP, &
                    PEPDIR,X,Y,Z,SYNLOOP,LOOP,LOOPLEN,DNOR,ATOR,EDHOLD,EDLOCKIN,GMODE,DREACT)

            ELSE

               CALL PAIRCRLK(NS,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,PEPDIR, &
                    X,Y,Z,ATOR,CRLKAGE,DREACT)

            END IF

         END IF

!        THE THIRD TPASE CROSSLINKING:

         IF(TPPEP(1,3,NS)>0.AND.TPPEP(2,3,NS)==0)THEN

            JS=3

            CALL CRLK2SAC(JS,NS,XTPASE,YTPASE,ZTPASE,TPPEP,SYNPG,PGID,PGLEN,BONDPEP,BONTYP, &
                 PEPDIR,X,Y,Z,SYNLOOP,LOOP,LOOPLEN,DNOR,ATOR,EDHOLD,EDLOCKIN,GMODE,DREACT)

         END IF


!        SIGNAL FOR THE FIRST TPASE TO CROSSLINK:

         IF(TPPEP(2,1,NS)>0)THEN

            IF(GMODE<12)THEN
               SIGCROSS(1,NS)=1
            ELSE
               JS=1
               CALL CROSSSIGNAL(NS,JS,TPPEP,X,Y,Z,XGTASE,YGTASE,ZGTASE,SIGCROSS,DELTA)
            END IF

         END IF

!        SIGNAL FOR THE SECOND TPASE TO CROSSLINK:

         IF(TPPEP(2,2,NS)>0)THEN

            IF(GMODE<12)THEN
               SIGCROSS(2,NS)=1
            ELSE
               JS=2
               CALL CROSSSIGNAL(NS,JS,TPPEP,X,Y,Z,XGTASE,YGTASE,ZGTASE,SIGCROSS,DELTA)

            END IF

         END IF

!        SIGNAL FOR THE THIRD TPASE TO CROSSLINK:

         IF(TPPEP(2,3,NS)>0)THEN

            JS=3
            CALL CROSSSIGNAL(NS,JS,TPPEP,X,Y,Z,XGTASE,YGTASE,ZGTASE,SIGCROSS,DELTA)

         END IF



!----------------------------

!     POST CROSSLINKING:


         JCHANGE=0

         IF(SIGCROSS(1,NS)==1)THEN

            JS=1

            CALL POSTCRLK(JS,NS,TPPEP,CRLKAGE,SIGCROSS,SYNLOOP,NATOM,DNOR,ATOR,X,Y,Z, &
                  XGTASE,YGTASE,ZGTASE,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                   NBONDGLY,BONDGLY,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,JFORCE,EDPEP,JDEACT,GMODE)


            JCHANGE=1

         END IF

         IF(SIGCROSS(2,NS)==1)THEN

            JS=2

            CALL POSTCRLK(JS,NS,TPPEP,CRLKAGE,SIGCROSS,SYNLOOP,NATOM,DNOR,ATOR,X,Y,Z, &
                  XGTASE,YGTASE,ZGTASE,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                   NBONDGLY,BONDGLY,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,JFORCE,EDPEP,JDEACT,GMODE)


            JCHANGE=1

         END IF

         IF(SIGCROSS(3,NS)==1)THEN

            JS=3

            CALL POSTCRLK(JS,NS,TPPEP,CRLKAGE,SIGCROSS,SYNLOOP,NATOM,DNOR,ATOR,X,Y,Z, &
                  XGTASE,YGTASE,ZGTASE,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                   NBONDGLY,BONDGLY,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,JFORCE,EDPEP,JDEACT,GMODE)


            JCHANGE=1

         END IF

         IF(JCHANGE==1)THEN

            CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                      NEIBOND,NEINUM,NEIGLY,GLYNUM)

            CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

            JTRAP(NS)=JSTEP


         END IF


!============================  ENDOPEPTIDATION RULES:




         IF(EDHOLD(1,NS)>0.AND.EDLOCKIN(NS)==0)THEN

            JRANDS(NS)=JRANDS(NS)+1

            IF(P_EDHOLD>RANDS(JRANDS(NS),NS))THEN
               EDHOLD(1:3,NS)=0
            END IF

            IF((TPPEP(1,1,NS)==0.OR.TPPEP(1,3,NS))==0.AND.GMODE>=12)THEN
               EDHOLD(1:3,NS)=0
            END IF

         END IF

!        TURN ON CAPTURING MODE:

         IF(EDHOLD(1,NS)==0.AND.EDPEP(1,NS)==0.AND.EDPEP(2,NS)==0.AND.EDCAP(NS)==0)THEN
            IF(GMODE==0)THEN
               EDCAP(NS)=1
            ELSEIF(SYNLOOP(NS)>0)THEN
               EDCAP(NS)=1
            END IF

            IF(GMODE==13.AND.(TPPEP(1,1,NS)==0.OR.TPPEP(1,3,NS)==0))THEN
               EDCAP(NS)=0
            END IF
         END IF

!        edase is in capture mode, then capturing crosslinks:

         IF(EDHOLD(1,NS)==0.AND.EDPEP(1,NS)==0.AND.EDPEP(2,NS)==0.AND.EDCAP(NS)==1)THEN

            CALL ENDOPEP(NS,SYNLOOP,LOOP,LOOPLEN,NBONDPEP,BONDPEP,BONTYP,SYNPG,OLDSYNPG,PGID,PGLEN,EDHOLD, &
                      DNOR,X,Y,Z,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD,GMODE)

         END IF

!        to signal edase to cleave if it is holding a crosslink:

         IF(EDHOLD(1,NS)>0.AND.MOD(JSTEP,10)==0)THEN

            IF(GMODE<7)THEN
               SIGCLEAVE(NS)=1
            ELSE
               CALL CLEAVESIGNAL(NS,EDHOLD,X,Y,Z,SIGCLEAVE,EDLOCKIN,BONTYP,GMODE)
            END IF

         END IF

!        ENDOPEPTIDASE RELEASING PEPTIDE:

         IF((EDPEP(1,NS)>0.OR.EDPEP(2,NS)>0).AND.MOD(JSTEP-1,10)==0.AND.GMODE<13)THEN

            ILOOP=SYNLOOP(NS)

            JCHECK1=0

            JCHECK2=0

            DO JL=1,LOOPLEN(ILOOP)

               IF(EDPEP(1,NS)==LOOP(JL,ILOOP))THEN
                  JCHECK1=1
               END IF

               IF(EDPEP(2,NS)==LOOP(JL,ILOOP))THEN
                  JCHECK2=1
               END IF

            END DO

            IF(JCHECK1==0)THEN
               EDPEP(1,NS)=0
            END IF

            IF(JCHECK2==0)THEN
               EDPEP(2,NS)=0
            END IF

         END IF

!        IMMEDIATE RELEASE DUE TO PULLING FORCE:

         IF((EDPEP(1,NS)>0.OR.EDPEP(2,NS)>0).AND.GMODE==12)THEN

            DX=XEDASE(NS)-0.5D0*(XGTASE(1,NS)+XGTASE(2,NS))
            DY=YEDASE(NS)-0.5D0*(YGTASE(1,NS)+YGTASE(2,NS))
            DZ=ZEDASE(NS)-0.5D0*(ZGTASE(1,NS)+ZGTASE(2,NS))

            IF(DX*XLEAD(NS)+DY*YLEAD(NS)+DZ*ZLEAD(NS)<-1.0D0)THEN
               EDPEP(1:2,NS)=0
            END IF

         END IF

!        PROBABILISTIC RELEASE:

         IF(EDPEP(1,NS)>0)THEN

            JRANDS(NS)=JRANDS(NS)+1

            IF(P_EDHOLD>RANDS(JRANDS(NS),NS))THEN
               EDPEP(1,NS)=0
            END IF

         END IF

         IF(EDPEP(2,NS)>0)THEN

            JRANDS(NS)=JRANDS(NS)+1

            IF(P_EDHOLD>RANDS(JRANDS(NS),NS))THEN
               EDPEP(2,NS)=0
            END IF

         END IF



!----------------------------------

!     POST CLEAVING:

         IF(SIGCLEAVE(NS)==1)THEN

            CALL POSTCLEAVE(NS,NSYN,SYNLOOP,LOOP,LOOPLEN,LOOPTYP,NLOOP,NLOOPDEL,EDHOLD,SIGCLEAVE,EDLOCKIN, &
                            JDEACT,JFORCE,X,Y,Z,NBONDDEL,NBONDPEP,BONDPEP,BONTYP,DNOR,ATOR,EDPEP,EDCAP,GMODE)

            CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                      NEIBOND,NEINUM,NEIGLY,GLYNUM)

            CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

            JTRAP(NS)=JSTEP

         END IF

      END DO

!----------------------------------
!     TO REINITIATE COMPLEXES AT POLES:

      IF(GMODE==10.OR.GMODE==11)THEN

         IF(MOD(JSTEP-1,10000)==0)THEN


            DO NS=1,NSYN
               IF(ABS(XCAP1-XEDASE(NS))<2.0D0.OR.ABS(XCAP2-XEDASE(NS))<2.0D0)THEN
                  JDEACT(NS)=1
                  JREINI(NS)=1
               END IF
            END DO



         END IF

      END IF


!======================= TERMINATION AND DEACTIVATION


      DO NS=1,NSYN

         IF(SYNPG(1,NS)==0.AND.SYNPG(2,NS)==0.AND.SYNTHESIS(1,NS)==1.AND.SYNTHESIS(2,NS)==1)THEN

            JRANDS(NS)=JRANDS(NS)+1

            IF(PSTOP>RANDS(JRANDS(NS),NS))THEN
               JDEACT(NS)=1
            END IF

         END IF

!        if there is deactivating signal:

         IF(JDEACT(NS)==1.OR.JSTEP-JTRAP(NS)>NTRAP)THEN

            JDEACT(NS)=0

            SYNTHESIS(1:2,NS)=0

            OLDSYNPG(1:2,NS)=SYNPG(1:2,NS)

            SYNPG(1:2,NS)=0

            GLYTIP(1:2,NS)=0

            GLYSEC(1:2,NS)=0

            GTLOAD(1:2,NS)=0

            CRLKAGE(1:2,NS)=0

            IF(SYNLOOP(NS)>0)THEN

               LOOPTYP(SYNLOOP(NS))=1

               DO NS1=1,NSYN

                  IF(NS1/=NS.AND.SYNLOOP(NS1)==SYNLOOP(NS))THEN
                     LOOPTYP(SYNLOOP(NS))=2
                     EXIT
                  END IF

               END DO

            END IF

            SYNLOOP(NS)=0

            TPPEP(1:2,1,NS)=0

            TPPEP(1:2,2,NS)=0

            TPPEP(1:2,3,NS)=0

            SIGCROSS(1:3,NS)=0

            EDHOLD(1:3,NS)=0

            EDPEP(1:2,NS)=0

            SIGCLEAVE(NS)=0

            EDLOCKIN(NS)=0

            EDCAP(NS)=0

            JTRAP(NS)=JSTEP

            JFORCE=1

            GTAWAIT(NS)=0


            CYCLE

         END IF

!----------------------

!        FIRST STRAND TERMINATION

         IF(GMODE>=12.AND.(CRLKAGE(1,NS)<2.OR.CRLKAGE(2,NS)<2.OR.TPPEP(1,1,NS)/=0.OR.TPPEP(1,2,NS)/=0))THEN
            GOTO 22
         END IF

!        AVOID TERMINATION WHILE TPASE IS CROSSLINKING IN THE SINGLE-STRAND MODE:

         IF(GMODE<12.AND.(SYNPG(1,NS)==0.OR.TPPEP(1,1,NS)/=0))THEN
            CYCLE
         END IF

!        WHEN PAIRED TPASES HYPOTHESIS APPLIED:

         IF(GMODE>6.AND.GMODE<12.AND.TPPEP(1,2,NS)/=0)THEN
            CYCLE
         END IF

!        CROSSLINK-BEFORE-TERMINATE:

         IF(GMODE>4.AND.GMODE<13.AND.(EDPEP(1,NS)/=0.OR.EDPEP(2,NS)/=0))THEN
            CYCLE
         END IF

         IPG=SYNPG(1,NS)

         ILOOP=SYNLOOP(NS)

!        HOLE-DEPENDENT PROCESSIVITY:

         IF(GMODE<11)THEN
            AREA=1.0D0
         ELSE
            AREA=MAX(1.0D0,0.01D0*LOOPLEN(ILOOP)**2)
         END IF

!        CROSSLINK-DEPENDENT PROCESSIVITY:

         IF(GMODE<9)THEN
            JCRLK=1
         ELSEIF(GMODE<12)THEN
            JCRLK=CRLKAGE(1,NS)-2
         ELSE
            JCRLK=CRLKAGE(1,NS)
         END IF


         JLONG=0

!        LIMIT GLYCAN LENGTH:

         IF(GMODE>2.AND.PGLEN(IPG)>=25)THEN
            JLONG=1
         END IF

         JRANDS(NS)=JRANDS(NS)+1

         IF(PTERM*JCRLK>RANDS(JRANDS(NS),NS)*AREA.OR.JLONG==1)THEN

            OLDSYNPG(1,NS)=SYNPG(1,NS)

            SYNPG(1,NS)=0

            GTLOAD(1,NS)=0

            GLYTIP(1,NS)=0

            GLYSEC(1,NS)=0

            SIGCROSS(1:2,NS)=0

            CRLKAGE(1,NS)=0

            JFORCE=1


         END IF


!        SECOND STRAND TERMINATION

         IF(GMODE<12)THEN
            CYCLE
         END IF

22       IF(CRLKAGE(2,NS)>1.AND.CRLKAGE(1,NS)>1.AND.TPPEP(1,3,NS)==0.AND.TPPEP(1,2,NS)==0)THEN

            IPG=SYNPG(2,NS)

            ILOOP=SYNLOOP(NS)

            AREA=MAX(1.0D0,0.01D0*LOOPLEN(ILOOP)**2)

            JRANDS(NS)=JRANDS(NS)+1

            IF(PTERM*CRLKAGE(2,NS)>RANDS(JRANDS(NS),NS)*AREA.OR.PGLEN(IPG)>=25)THEN

               OLDSYNPG(2,NS)=SYNPG(2,NS)

               SYNPG(2,NS)=0

               GTLOAD(2,NS)=0

               GLYTIP(2,NS)=0

               GLYSEC(2,NS)=0

               SIGCROSS(2:3,NS)=0

               CRLKAGE(2,NS)=0

               JFORCE=1

            END IF

         END IF


      END DO


!================================================

!     UPDATE SYSTEM:

      IF(NATOM>NATOMMAX-10.OR.NPG>NPGMAX-10.OR.NBONDPEP>NPEPMAX-10.OR.NLOOP>NLOOPMAX-10)THEN

         ALLOCATE(XNEW(NATOM),YNEW(NATOM),ZNEW(NATOM),DNORNEW(NATOM),ATORNEW(NATOM),PEPDIRNEW(NATOM),ATOMRADNEW(NATOM))

         XNEW(1:NATOM)=X(1:NATOM)
         YNEW(1:NATOM)=Y(1:NATOM)
         ZNEW(1:NATOM)=Z(1:NATOM)

         DNORNEW(1:NATOM)=DNOR(1:NATOM)
         ATORNEW(1:NATOM)=ATOR(1:NATOM)

         PEPDIRNEW(1:NATOM)=PEPDIR(1:NATOM)
         ATOMRADNEW(1:NATOM)=ATOMRAD(1:NATOM)

         DEALLOCATE(X,Y,Z,DNOR,ATOR,PEPDIR,ATOMRAD)

         NATOMMAX=NATOM+1000

         ALLOCATE(X(NATOMMAX),Y(NATOMMAX),Z(NATOMMAX),DNOR(NATOMMAX),ATOR(NATOMMAX),PEPDIR(NATOMMAX),ATOMRAD(NATOMMAX))

         X(1:NATOM)=XNEW(1:NATOM)
         Y(1:NATOM)=YNEW(1:NATOM)
         Z(1:NATOM)=ZNEW(1:NATOM)

         DNOR(1:NATOM)=DNORNEW(1:NATOM)
         ATOR(1:NATOM)=ATORNEW(1:NATOM)

         PEPDIR(1:NATOM)=PEPDIRNEW(1:NATOM)
         ATOMRAD(1:NATOM)=ATOMRADNEW(1:NATOM)

         DEALLOCATE(XNEW,YNEW,ZNEW,DNORNEW,ATORNEW,PEPDIRNEW,ATOMRADNEW)

!        ---------------------------------------------------------------------

         ALLOCATE(NEWPGID(PGLENMAX,NPG),NEWPGLEN(NPG),NEWPGTYP(NPG),NEWPGDIR(NPG))

         DO NR=1,NPG
            NEWPGID(1:PGLEN(NR),NR)=PGID(1:PGLEN(NR),NR)
         END DO

         NEWPGLEN(1:NPG)=PGLEN(1:NPG)
         NEWPGTYP(1:NPG)=PGTYP(1:NPG)
         NEWPGDIR(1:NPG)=PGDIR(1:NPG)

         DEALLOCATE(PGID,PGLEN,PGTYP,PGDIR)

         NPGMAX=NPG+1000

         ALLOCATE(PGID(PGLENMAX,NPGMAX),PGLEN(NPGMAX),PGTYP(NPGMAX),PGDIR(NPGMAX))

         DO NR=1,NPG
            PGID(1:NEWPGLEN(NR),NR)=NEWPGID(1:NEWPGLEN(NR),NR)
         END DO

         PGLEN(1:NPG)=NEWPGLEN(1:NPG)
         PGTYP(1:NPG)=NEWPGTYP(1:NPG)
         PGDIR(1:NPG)=NEWPGDIR(1:NPG)

         DEALLOCATE(NEWPGID,NEWPGLEN,NEWPGTYP,NEWPGDIR)

!        ---------------------------------------------------------------------
         ALLOCATE(NEWBONDGLY(2,NBONDGLY))

         NEWBONDGLY(1,1:NBONDGLY)=BONDGLY(1,1:NBONDGLY)
         NEWBONDGLY(2,1:NBONDGLY)=BONDGLY(2,1:NBONDGLY)

         DEALLOCATE(BONDGLY)

         NGLYMAX=NBONDGLY+1000

         ALLOCATE(BONDGLY(2,NGLYMAX))

         BONDGLY(1,1:NBONDGLY)=NEWBONDGLY(1,1:NBONDGLY)
         BONDGLY(2,1:NBONDGLY)=NEWBONDGLY(2,1:NBONDGLY)

         DEALLOCATE(NEWBONDGLY)

!        ---------------------------------------------------------------------
         ALLOCATE(NEWBONDPEP(2,NBONDPEP-NBONDDEL),NEWBONTYP(NBONDPEP-NBONDDEL))

         NNEW=0

         DO NB=1,NBONDPEP
            IF(BONTYP(NB)==0)THEN
               CYCLE
            END IF

            NNEW=NNEW+1
            NEWBONTYP(NNEW)=BONTYP(NB)
            NEWBONDPEP(1:2,NNEW)=BONDPEP(1:2,NB)

            DO NS=1,NSYN

               IF(EDHOLD(3,NS)==NB)THEN
                  EDHOLD(3,NS)=NNEW
                  EXIT
               END IF

            END DO

         END DO

         IF(NNEW/=NBONDPEP-NBONDDEL)THEN
            PRINT*,'ERROR IN COUNTING DELETED BONDS'
            STOP
         END IF

         NBONDPEP=NNEW
         NBONDDEL=0

         DEALLOCATE(BONDPEP,BONTYP)

         NPEPMAX=NBONDPEP+1000

         ALLOCATE(BONDPEP(2,NPEPMAX),BONTYP(NPEPMAX))

         BONDPEP(1,1:NBONDPEP)=NEWBONDPEP(1,1:NBONDPEP)
         BONDPEP(2,1:NBONDPEP)=NEWBONDPEP(2,1:NBONDPEP)
         BONTYP(1:NBONDPEP)=NEWBONTYP(1:NBONDPEP)

         CALL SETSIGBOND(NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND)

         DEALLOCATE(NEWBONDPEP,NEWBONTYP)

!        ---------------------------------------------------------------------
         LOOPLENMAX=MAXVAL(LOOPLEN(1:NLOOP))

         ALLOCATE(NEWLOOP(LOOPLENMAX,NLOOP-NLOOPDEL),NEWLOOPLEN(NLOOP-NLOOPDEL),NEWLOOPTYP(NLOOP-NLOOPDEL))

         NNEW=0

         DO NL=1,NLOOP

            IF(LOOPTYP(NL)==0)THEN
               CYCLE
            END IF

            NNEW=NNEW+1

            NEWLOOP(1:LOOPLEN(NL),NNEW)=LOOP(1:LOOPLEN(NL),NL)
            NEWLOOPLEN(NNEW)=LOOPLEN(NL)
            NEWLOOPTYP(NNEW)=LOOPTYP(NL)

            IF(LOOPTYP(NL)==2)THEN
               DO NS=1,NSYN
                  IF(SYNLOOP(NS)==NL)THEN
                     SYNLOOP(NS)=NNEW
                     EXIT
                  END IF
               END DO
            END IF

         END DO

         IF(NNEW/=NLOOP-NLOOPDEL)THEN
            PRINT*,'ERROR IN COUNTING DELETED LOOPS'
            STOP
         END IF

         NLOOP=NNEW
         NLOOPDEL=0

         DEALLOCATE(LOOP,LOOPLEN,LOOPTYP)

         NLOOPMAX=NLOOP+1000
         LOOPLENMAX=MAX(LOOPLENMAX*2,60)

         ALLOCATE(LOOP(LOOPLENMAX,NLOOPMAX),LOOPLEN(NLOOPMAX),LOOPTYP(NLOOPMAX))

         DO NL=1,NLOOP
            LOOP(1:NEWLOOPLEN(NL),NL)=NEWLOOP(1:NEWLOOPLEN(NL),NL)
         END DO

         LOOPLEN(1:NLOOP)=NEWLOOPLEN(1:NLOOP); LOOPLEN(NLOOP+1:NLOOPMAX)=0

         LOOPTYP(1:NLOOP)=NEWLOOPTYP(1:NLOOP)

         DEALLOCATE(NEWLOOP,NEWLOOPLEN,NEWLOOPTYP)

!------------------------

         DEALLOCATE(FX,FY,FZ,FAMAG)
         DEALLOCATE(FXGLY,FYGLY,FZGLY)
         DEALLOCATE(FXPEP,FYPEP,FZPEP)
         DEALLOCATE(FXTHETA,FYTHETA,FZTHETA)
         DEALLOCATE(FXPRES,FYPRES,FZPRES)
         DEALLOCATE(FXSYN,FYSYN,FZSYN)
         DEALLOCATE(FYSUR,FZSUR)

         ALLOCATE(FX(NATOMMAX),FY(NATOMMAX),FZ(NATOMMAX),FAMAG(NATOMMAX))
         ALLOCATE(FXGLY(NATOMMAX),FYGLY(NATOMMAX),FZGLY(NATOMMAX))
         ALLOCATE(FXPEP(NATOMMAX),FYPEP(NATOMMAX),FZPEP(NATOMMAX))
         ALLOCATE(FXTHETA(NATOMMAX),FYTHETA(NATOMMAX),FZTHETA(NATOMMAX))
         ALLOCATE(FXPRES(NATOMMAX),FYPRES(NATOMMAX),FZPRES(NATOMMAX))
         ALLOCATE(FXSYN(NATOMMAX),FYSYN(NATOMMAX),FZSYN(NATOMMAX))
         ALLOCATE(FYSUR(NATOMMAX),FZSUR(NATOMMAX))

         JFORCE=1

         DO NS=1,NSYN
            IF(SYNLOOP(NS)/=0)THEN

               CALL SETNEIBOND(NS,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP,X,Y,Z,XGTASE,YGTASE,ZGTASE, &
                      NEIBOND,NEINUM,NEIGLY,GLYNUM)


            END IF

         END DO

      END IF

!     ------------------------------------------------------------------------
!     this is to cope with large loop

      IF(MOD(JSTEP-1,1000)==0)THEN

         IF(MAXVAL(LOOPLEN(1:NLOOP))*2>LOOPLENMAX)THEN

            LOOPLENMAX=MAXVAL(LOOPLEN(1:NLOOP))*2

            ALLOCATE(NEWLOOP(MAXVAL(LOOPLEN(1:NLOOP)),NLOOP))

            DO NL=1,NLOOP
               NEWLOOP(1:LOOPLEN(NL),NL)=LOOP(1:LOOPLEN(NL),NL)
            END DO

            DEALLOCATE(LOOP)

            ALLOCATE(LOOP(LOOPLENMAX,NLOOPMAX))

            DO NL=1,NLOOP
               LOOP(1:LOOPLEN(NL),NL)=NEWLOOP(1:LOOPLEN(NL),NL)
            END DO

            DEALLOCATE(NEWLOOP)

         END IF



      END IF

!======================================================================

!     FORCE CALCULATION:

      CALL FGLYCANS(FXGLY,FYGLY,FZGLY,FXTHETA,FYTHETA,FZTHETA, &
                     NPG,PGID,PGLEN,PGTYP,X,Y,Z,L_G,K_G,THETA_0,KTHETA,JFORCE,DELTA,INVDELTA,BETA,PI)


      CALL FPEPBONDS(JFORCE,FXPEP,FYPEP,FZPEP,X,Y,Z,NBONDPEP,BONDPEP,BONTYP,NSIG,SIGBOND, &
                     L_P,K_P,LSTART,LREP,LSWITCH,KSWITCH,MSWITCH)

      CALL FPRES(JFORCE,NATOM,FXPRES,FYPRES,FZPRES,X,Y,Z,NLOOP,LOOP,LOOPLEN,LOOPTYP,PRES,NTHREAD,DELTA,INVDELTA)

!     -----------------------------

      CALL RANDFORCES(NSYN,FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE,FXEDASE,FYEDASE,FZEDASE, &
                      RXGTASE,RYGTASE,RZGTASE,RXTPASE,RYTPASE,RZTPASE,RXEDASE,RYEDASE,RZEDASE,JRFORCE, &
                      YGTASE,ZGTASE,YTPASE,ZTPASE,YEDASE,ZEDASE)


!     -----------------------------

      CALL SYNCOMP(NSYN,SYNDIR,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE, &
                      KPAIR,LPAIR,KGTTP,LGTTP,KGTED,LGTED,KSIDE,LSIDE,DELTA,INVDELTA,BETA, &
                      FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE,FXEDASE,FYEDASE,FZEDASE, &
                      XLEAD,YLEAD,ZLEAD,KLEAD,GMODE)


!     -----------------------------

      CALL GTAHOLD(NSYN,GLYTIP,GLYSEC,KGTASE,LGTASE,KLEAD,KTHETA,DELTA,INVDELTA,BETA,X,Y,Z, &
                      XLEAD,YLEAD,ZLEAD,FXSYN,FYSYN,FZSYN,XGTASE,YGTASE,ZGTASE,FXGTASE,FYGTASE,FZGTASE)


!     -----------------------------

      CALL TPAHOLD(NSYN,X,Y,Z,FXSYN,FYSYN,FZSYN,TPPEP,XTPASE,YTPASE,ZTPASE,FXTPASE,FYTPASE,FZTPASE,KTPASE,LTPASE)

!     -----------------------------

      CALL EDAHOLD(NSYN,EDHOLD,EDPEP,KEDASE,LEDASE,X,Y,Z,FXSYN,FYSYN,FZSYN, &
                   XEDASE,YEDASE,ZEDASE,FXEDASE,FYEDASE,FZEDASE,PEPDIR,GMODE)

!     -----------------------------------------

      CALL STERIC(NSYN,SYNTHESIS,NEINUM,NEIBOND,GLYNUM,NEIGLY,GLYTIP,KWALL, &
                     X,Y,Z,FXSYN,FYSYN,FZSYN,FXGTASE,FYGTASE,FZGTASE,XGTASE,YGTASE,ZGTASE)

!     -----------------------------------------
!     CONSTRAIN PBPS AND UNLINKED PG TO THE SURFACE:

      IF(MOD(JSTEP-1,200)==0)THEN
         CALL CALSURF(NSYN,SYNPG,PGID,PGLEN,XEDASE,YEDASE,ZEDASE,NATOM,DNOR,ATOR,X,Y,Z, &
                      SYNRAD,ATOMRAD,GTRAD,GMODE,XGTASE,YGTASE,ZGTASE)
      END IF

      CALL SURFLINK(NSYN,SYNRAD,GTRAD,YGTASE,ZGTASE,YTPASE,ZTPASE,YEDASE,ZEDASE,FYGTASE,FZGTASE,FYTPASE,FZTPASE, &
                    FYEDASE,FZEDASE,NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,ATOMRAD,Y,Z,FYSUR,FZSUR,KSUR)

!     CONSTRAIN ENZYMES WITHIN TWO CAPS:

      IF(MOD(JSTEP-1,100)==0)THEN

         FXCAPEDASE(1:NSYN)=0.0D0
         FXCAPGT(1:NSYN)=0.0D0


         DO NS=1,NSYN
            IF(XEDASE(NS)<XCAP1)THEN
               FXCAPEDASE(NS)=10*(XCAP1-XEDASE(NS))
            ELSEIF(XEDASE(NS)>XCAP2)THEN
               FXCAPEDASE(NS)=10*(XCAP2-XEDASE(NS))
            END IF

            IF(GMODE>0)THEN
               FXCAPGT(NS)=FXCAPEDASE(NS)
               CYCLE
            END IF

            IF(XGTASE(1,NS)<XCAP1)THEN
               FXCAPGT(NS)=10*(XCAP1-XGTASE(1,NS))
            ELSEIF(XGTASE(1,NS)>XCAP2)THEN
               FXCAPGT(NS)=10*(XCAP2-XGTASE(1,NS))
            END IF

         END DO


      END IF

      FXGTASE(1,1:NSYN)=FXGTASE(1,1:NSYN)+FXCAPGT(1:NSYN)
      FXGTASE(2,1:NSYN)=FXGTASE(2,1:NSYN)+FXCAPGT(1:NSYN)

      FXTPASE(1,1:NSYN)=FXTPASE(1,1:NSYN)+FXCAPGT(1:NSYN)
      FXTPASE(2,1:NSYN)=FXTPASE(2,1:NSYN)+FXCAPGT(1:NSYN)
      FXTPASE(3,1:NSYN)=FXTPASE(3,1:NSYN)+FXCAPGT(1:NSYN)

      FXEDASE(1:NSYN)=FXEDASE(1:NSYN)+FXCAPEDASE(1:NSYN)

!======================================================================

      FX=FXGLY+FXPEP+FXTHETA+FXPRES+FXSYN
      FY=FYGLY+FYPEP+FYTHETA+FYPRES+FYSYN+FYSUR
      FZ=FZGLY+FZPEP+FZTHETA+FZPRES+FZSYN+FZSUR

      FAMAG(1:NATOM)=ABS(FX(1:NATOM))+ABS(FY(1:NATOM))+ABS(FZ(1:NATOM))

      FAMAX=MAXVAL(FAMAG(1:NATOM))


      FMAGPBP(1:NSYN)=ABS(FXGTASE(1,1:NSYN))+ABS(FYGTASE(1,1:NSYN))+ABS(FZGTASE(1,1:NSYN))

      FPMAX1=MAXVAL(FMAGPBP(1:NSYN))

      IF(GMODE>=12)THEN
         FMAGPBP(1:NSYN)=ABS(FXGTASE(2,1:NSYN))+ABS(FYGTASE(2,1:NSYN))+ABS(FZGTASE(2,1:NSYN))
         FPMAX2=MAXVAL(FMAGPBP(1:NSYN))
      ELSE
         FPMAX2=0.0D0
      END IF


      FMAGPBP(1:NSYN)=ABS(FXTPASE(1,1:NSYN))+ABS(FYTPASE(1,1:NSYN))+ABS(FZTPASE(1,1:NSYN))

      FPMAX3=MAXVAL(FMAGPBP(1:NSYN))

      IF(GMODE>6)THEN
         FMAGPBP(1:NSYN)=ABS(FXTPASE(2,1:NSYN))+ABS(FYTPASE(2,1:NSYN))+ABS(FZTPASE(2,1:NSYN))
         FPMAX4=MAXVAL(FMAGPBP(1:NSYN))
      ELSE
         FPMAX4=0.0D0
      END IF

      IF(GMODE>=12)THEN
         FMAGPBP(1:NSYN)=ABS(FXTPASE(3,1:NSYN))+ABS(FYTPASE(3,1:NSYN))+ABS(FZTPASE(3,1:NSYN))
         FPMAX5=MAXVAL(FMAGPBP(1:NSYN))
      ELSE
         FPMAX5=0.0D0
      END IF

      FMAGPBP(1:NSYN)=ABS(FXEDASE(1:NSYN))+ABS(FYEDASE(1:NSYN))+ABS(FZEDASE(1:NSYN))

      FPMAX6=MAXVAL(FMAGPBP(1:NSYN))

      FPBPMAX=MAX(FPMAX1,FPMAX2,FPMAX3,FPMAX4,FPMAX5,FPMAX6)

      GAM=MIN(0.005D0/FAMAX,0.1D0/FPBPMAX)

!---------------------------------
!-------------------------------------------------------

      CALL NEWCOOR(NATOM,X,Y,Z,FX,FY,FZ,NSYN,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE, &
                   XEDASEOLD,YEDASEOLD,ZEDASEOLD,FXGTASE,FYGTASE,FZGTASE,FXTPASE,FYTPASE,FZTPASE, &
                   FXEDASE,FYEDASE,FZEDASE,INVMPBP,GAM,GMODE)

!---------------------------------
!---------------------------------

      IF(MOD(JSTEP-1,JPRINT1)==0)THEN

         WRITE(*,52)'STEP',JSTEP,'MAX FORCE (nN)',0.01*FAMAX,0.01*FPBPMAX

         CALL WRITEDCD(JUNIT,NTOTAL,NATOMSTART,NATOM,NPGSTART,NPG,DNOR,ATOR,PEPDIR,X,Y,Z,SYNPG, &
              SYNTHESIS,PGID,PGLEN,NGLYNEW,NBONDPEP,BONDPEP,BONTYP,NPEPTOTAL,NSYN,NSYNMAX,TPPEP, &
              EDHOLD,EDPEP,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,NFRAME)

         IF(NFRAME>=500)THEN
            CLOSE(JUNIT)
            JFILE=JFILE+1
            CALL DCDHEADER(JUNIT,JFILE,NTOTAL)
            NFRAME=0
         END IF


      END IF

!---------------------------------
!---------------------------------

      IF(MOD(JSTEP,JPRINT2)==0)THEN


         CALL PG_OUT(NSTART,NATOM,NATOMDEL,OLDNATOMDEL,NATOMCAP,NATOMSTART,DNOR,ATOR,PEPDIR,X,Y,Z,PGCAP1,PGCAP2,ORAD, &
                 NPG,NPGSTART,PGID,PGLEN,PGTYP,PGDIR,NBONDGLY,NBONDGLYSTART,NBONDPEP,NBONDDEL, &
                 BONDGLY,BONDPEP,BONTYP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                 NSYN,SYNDIR,SYNTHESIS,GTLOAD,SYNPG,SYNLOOP,GLYTIP,GLYSEC,TPPEP,EDPEP, &
                 GTATRANS,JDEACT,SIGCROSS,SIGCLEAVE,EDLOCKIN,EDCAP,EDHOLD,CRLKAGE,SYNRATIO, &
                 XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD)


         CALL WHATTIME(TIMERUN)

         OPEN(10,FILE='restart.inp')
         WRITE(10,*)NSTART,JFILE
         WRITE(10,*)TIMERUN-TIMESTART+LASTTIME
         CLOSE(10)

      END IF

!---------------------------------

      IF(NATOM>=NATOMSTART+NGLYNEW-50)THEN

         IF(MOD(JSTEP,10)==0)THEN

            CALL WRITEDCD(JUNIT,NTOTAL,NATOMSTART,NATOM,NPGSTART,NPG,DNOR,ATOR,PEPDIR,X,Y,Z,SYNPG, &
              SYNTHESIS,PGID,PGLEN,NGLYNEW,NBONDPEP,BONDPEP,BONTYP,NPEPTOTAL,NSYN,NSYNMAX,TPPEP, &
              EDHOLD,EDPEP,XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,NFRAME)

            CALL PG_OUT(NSTART,NATOM,NATOMDEL,OLDNATOMDEL,NATOMCAP,NATOMSTART,DNOR,ATOR,PEPDIR,X,Y,Z,PGCAP1,PGCAP2,ORAD, &
                 NPG,NPGSTART,PGID,PGLEN,PGTYP,PGDIR,NBONDGLY,NBONDGLYSTART,NBONDPEP,NBONDDEL, &
                 BONDGLY,BONDPEP,BONTYP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                 NSYN,SYNDIR,SYNTHESIS,GTLOAD,SYNPG,SYNLOOP,GLYTIP,GLYSEC,TPPEP,EDPEP, &
                 GTATRANS,JDEACT,SIGCROSS,SIGCLEAVE,EDLOCKIN,EDCAP,EDHOLD,CRLKAGE,SYNRATIO, &
                 XGTASE,YGTASE,ZGTASE,XTPASE,YTPASE,ZTPASE,XEDASE,YEDASE,ZEDASE,XEDASEOLD,YEDASEOLD,ZEDASEOLD)


            CALL WHATTIME(TIMERUN)

            OPEN(10,FILE='restart.inp')
            WRITE(10,*)NSTART,JFILE
            WRITE(10,*)TIMERUN-TIMESTART+LASTTIME
            CLOSE(10)

            EXIT

         END IF

      END IF

   END DO  ! END DO INSERTION

!================================================

   TIMERUN=TIMERUN-TIMESTART+LASTTIME

   days=timerun/1440
   hours=(timerun-1440*days)/60
   mins=timerun-1440*days-60*hours

   write(*,91)'RUNNING TIME =',DAYS,'days : ',HOURS,'hours : ',MINS,'mins'


   IF(JOBID==2)THEN
      STOP
   END IF

!============================================================================
!============================================================================
!============================================================================
!============================================================================

73 IF(JOBID==3)THEN
      PRINT*,'Check mother cell presence: configmother.inp and coormother.inp'
   END IF

!  Create the restart file:

!   OPEN(10,FILE='restart.inp')
!   N=0
!   WRITE(10,*)N,N
!   WRITE(10,*)N
!   CLOSE(10)

!  Pick daughter cells:

!   PRINT*,'CREATE LEFT CELL?'
!   READ(*,*)ANSWER
!   IF(ANSWER(1:1)=='Y'.OR.ANSWER(1:1)=='y')THEN
!      LEFTCELL=1
!   ELSE
!      LEFTCELL=0
!   END IF

!   PRINT*,'CREATE RIGHT CELL?'
!   READ(*,*)ANSWER
!   IF(ANSWER(1:1)=='Y'.OR.ANSWER(1:1)=='y')THEN
!      RITECELL=1
!   ELSE
!      RITECELL=0
!   END IF

   IF(LEFTCELL==1.AND.RITECELL==1)THEN
      PRINT*,'CAN NOT CREATE BOTH DAUGHTER CELLS AT THE SAME TIME'
      STOP
   ELSEIF(LEFTCELL==0.AND.RITECELL==0)THEN
      STOP
   END IF

!   PRINT*,'NUMBER OF CYLINDRICAL PG WANTED?'
!   READ(*,*)NSIZE

!  To further stabilize caps:
   NCLOSEDHOOP=1

!  Read info from mother cell:
   NSTART=1000000
   CALL GETINFO(NSTART,NATOM,NPG,PGLENMAX,NBONDGLY,NBONDPEP,NLOOP,LOOPLENMAX)

   NATOMMAX=NATOM+1000
   NPGMAX=NPG+1000
   NGLYMAX=NBONDGLY+1000
   NPEPMAX=NBONDPEP+1000
   NLOOPMAX=NLOOP+1000
   LOOPLENMAX=MAX(LOOPLENMAX*2,60)
   PGLENMAX=PGLENMAX+50

   ALLOCATE(X(NATOMMAX),Y(NATOMMAX),Z(NATOMMAX),DNOR(NATOMMAX),ATOR(NATOMMAX),PEPDIR(NATOMMAX))
   ALLOCATE(PGID(PGLENMAX,NPGMAX),PGLEN(NPGMAX),PGTYP(NPGMAX),PGDIR(NPGMAX))
   ALLOCATE(BONDGLY(2,NGLYMAX),BONDPEP(2,NPEPMAX),BONTYP(NPEPMAX))
   ALLOCATE(LOOP(LOOPLENMAX,NLOOPMAX),LOOPLEN(NLOOPMAX),LOOPTYP(NLOOPMAX))

   CALL PGMOTHER(NATOMSTART,NATOM,NATOMDEL,OLDNATOMDEL,DNOR,ATOR,PEPDIR,X,Y,Z, &
                 NPG,PGID,PGLEN,PGTYP,PGDIR, &
                 NBONDGLY,NBONDPEP,NBONDDEL,BONDGLY,BONDPEP,BONTYP, &
                 NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP,NSIZE,ORAD)



!  Make all PG mature:
   DO NA=1,NATOM

      IF(DNOR(NA)==1)THEN
         DNOR(NA)=-1
      END IF
   END DO

!--------------------------------------------------------
!  Clean up the system:

10   CALL LOOPBREAK(NPG,PGID,PGLEN,PGTYP,NATOM,DNOR,ATOR,NBONDPEP,NBONDDEL,BONDPEP,BONTYP, &
                  NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP)

   CALL PEPCLEAN(NBONDPEP,NBONDDEL,BONDPEP,BONTYP,NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,NATOM, &
                 NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP)

   CALL LYTGLY(NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,NATOMDEL)

   NDEL=NATOMDEL-OLDNATOMDEL

   CALL PG_UPDATE(NATOM,NDEL,DNOR,ATOR,PEPDIR,X,Y,Z,NPG,PGID, &
                  PGLEN,PGTYP,PGDIR, &
                  NLOOP,LOOP,LOOPLEN,LOOPTYP,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP)

   IF(NDEL>0)THEN
      OLDNATOMDEL=NATOMDEL
      GOTO 10
   END IF

!  Mark strands as old PG:

   DO NR=1,NPG
      IF(PGTYP(NR)==3)THEN
         PGTYP(NR)=1
      END IF
   END DO

   NATOMSTART=NATOM
   NATOMDEL=0
   OLDNATOMDEL=0

   NPGSTART=NPG

!  Updates bonds:

   ALLOCATE(NEWBONDPEP(2,NBONDPEP-NBONDDEL),NEWBONTYP(NBONDPEP-NBONDDEL))

   NNEW=0

   DO NB=1,NBONDPEP
      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      NNEW=NNEW+1
      NEWBONTYP(NNEW)=BONTYP(NB)
      NEWBONDPEP(1:2,NNEW)=BONDPEP(1:2,NB)


   END DO

   IF(NNEW/=NBONDPEP-NBONDDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED BONDS'
      STOP
   END IF

   NBONDPEP=NNEW
   NBONDDEL=0

   DEALLOCATE(BONDPEP,BONTYP)

   NPEPMAX=NBONDPEP+1000

   ALLOCATE(BONDPEP(2,NPEPMAX),BONTYP(NPEPMAX))

   BONDPEP(1,1:NBONDPEP)=NEWBONDPEP(1,1:NBONDPEP)
   BONDPEP(2,1:NBONDPEP)=NEWBONDPEP(2,1:NBONDPEP)
   BONTYP(1:NBONDPEP)=NEWBONTYP(1:NBONDPEP)

   DEALLOCATE(NEWBONDPEP,NEWBONTYP)

!-----------------------------------
!  Update LOOPS:

   LOOPLENMAX=MAXVAL(LOOPLEN(1:NLOOP))

   ALLOCATE(NEWLOOP(LOOPLENMAX,NLOOP-NLOOPDEL),NEWLOOPLEN(NLOOP-NLOOPDEL))

   NNEW=0

   DO NL=1,NLOOP

      IF(LOOPTYP(NL)==0)THEN
         CYCLE
      END IF

      NNEW=NNEW+1

      NEWLOOP(1:LOOPLEN(NL),NNEW)=LOOP(1:LOOPLEN(NL),NL)
      NEWLOOPLEN(NNEW)=LOOPLEN(NL)

   END DO

   IF(NNEW/=NLOOP-NLOOPDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED LOOPS'
      STOP
   END IF

   NLOOP=NNEW
   NLOOPDEL=0

   DEALLOCATE(LOOP,LOOPLEN,LOOPTYP)

   NLOOPMAX=NLOOP+1000
   LOOPLENMAX=MAX(LOOPLENMAX*2,200)

   ALLOCATE(LOOP(LOOPLENMAX,NLOOPMAX),LOOPLEN(NLOOPMAX),LOOPTYP(NLOOPMAX))

   DO NL=1,NLOOP
      LOOP(1:NEWLOOPLEN(NL),NL)=NEWLOOP(1:NEWLOOPLEN(NL),NL)
   END DO

   LOOPLEN(1:NLOOP)=NEWLOOPLEN(1:NLOOP); LOOPLEN(NLOOP+1:NLOOPMAX)=0

   LOOPTYP=1

   DEALLOCATE(NEWLOOP,NEWLOOPLEN)

!============================

!  DEFINE THE DIVISION PLANE OF THE SACCULUS:

   XCEN0=SUM(X(1:NATOM))/NATOM

   XCAP1=XCEN0-1000.0D0
   XCAP2=XCEN0+1000.0D0

   DO NR=1,NPG

      IF(PGTYP(NR)/=0)THEN
         CYCLE
      END IF

      XPG=SUM(X(PGID(1:PGLEN(NR),NR)))/PGLEN(NR)

      IF(XPG>XCAP1.AND.XPG<XCEN0)THEN
         PGCAP1=NR
         XCAP1=XPG
      END IF

      IF(XPG<XCAP2.AND.XPG>XCEN0)THEN
         PGCAP2=NR
         XCAP2=XPG
      END IF

   END DO


   DELTAX=0.0D0

   XSHIFT=0.0D0

30 XCEN0=XCEN0+DELTAX

   XSHIFT=XSHIFT+DELTAX

   NCOUNT=0

   DO N=1,NATOM

      IF(X(N)>XCAP1.AND.X(N)<XCEN0.AND.LEFTCELL==1)THEN
         NCOUNT=NCOUNT+1
      END IF

      IF(X(N)<XCAP2.AND.X(N)>XCEN0.AND.RITECELL==1)THEN
         NCOUNT=NCOUNT+1
      END IF

   END DO

   IF(ABS(NCOUNT-NSIZE)>100)THEN

      IF(NCOUNT>NSIZE)THEN
         DELTAX=0.1D0*(RITECELL-LEFTCELL)
      ELSE
         DELTAX=0.1D0*(LEFTCELL-RITECELL)
      END IF

      GOTO 30

   END IF

   CALL DIVISIONPLANE(NTHREAD,NATOM,X,Y,Z,NPG,PGID,PGLEN,PGTYP,PGDIR,NBONDPEP,BONDPEP,BONTYP, &
                      NLOOP,LOOP,LOOPLEN,LOOPTYP,DNOR,ATOR,VX,VY,VZ,XCEN0,YCEN0,ZCEN0,RADIUS)

   VX0=VX; VY0=VY;VZ0=VZ

   PRINT*,'SACCULUS RADIUS AT THE DIVISION PLANE =',RADIUS
   PRINT*,'AXIS DIRECTION',VX0,VY0,VZ0

   VX0=XCEN0; VY0=YCEN0;VZ0=ZCEN0
   PRINT*,'CENTER AT',VX0,VY0,VZ0


   LENHOOP=2*PI*RADIUS/L_G


   IF(MOD(LENHOOP,2)==1)THEN
      LENHOOP=LENHOOP-1
   END IF


   PRINT*,'DIVISION HOOP LENGTH =',LENHOOP

   PRINT*,'-----------------------------------------------'
   WRITE(*,*)

   IF(LENHOOP>PGLENMAX)THEN

      PGLENMAX=LENHOOP+1

      ALLOCATE(NEWPGID(PGLENMAX,NPG))

      DO NR=1,NPG

         NEWPGID(1:PGLEN(NR),NR)=PGID(1:PGLEN(NR),NR)

      END DO

      DEALLOCATE(PGID)

      ALLOCATE(PGID(PGLENMAX,NPGMAX))

      DO NR=1,NPG

         PGID(1:PGLEN(NR),NR)=NEWPGID(1:PGLEN(NR),NR)

      END DO

      DEALLOCATE(NEWPGID)

   ENDIF

!---------------------------------------------------

!  MARK STRANDS WHICH CROSS THE CENTRAL PLANE:

   ALLOCATE(MARK(NPG+1000))
   MARK=0

   NDIV=0

   DO NR=1,NPG

      XPG=SUM(X(PGID(1:PGLEN(NR),NR)))/PGLEN(NR)

      IF(ABS(XPG-XCEN0)>40.0D0)THEN
         CYCLE
      END IF

      DO JR=1,PGLEN(NR)-1

         NA1=PGID(JR,NR)
         NA2=PGID(JR+1,NR)

         DX=X(NA1)-XCEN0
         DY=Y(NA1)-YCEN0
         DZ=Z(NA1)-ZCEN0

         ARG1=DX*VX+DY*VY+DZ*VZ

         DX=X(NA2)-XCEN0
         DY=Y(NA2)-YCEN0
         DZ=Z(NA2)-ZCEN0

         ARG2=DX*VX+DY*VY+DZ*VZ


         IF(ARG1*ARG2<0.0D0)THEN

            NDIV=NDIV+1
            MARK(NR)=1
         END IF

      END DO

   END DO

!================================================
!  DIVIDE MARKED STRANDS AT ANY POINT THEY CROSS THE CENTRAL PLANE:

   DO N=1,NDIV

      CALL DIVSTRAND(X,Y,Z,XCEN0,YCEN0,ZCEN0,VX,VY,VZ,NPG,PGID,PGLEN,PGTYP,PGDIR,MARK, &
           NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP,DNOR,ATOR,NBONDPEP,NBONDDEL,BONDPEP,BONTYP)

      CALL LYTGLY(NPG,PGID,PGLEN,MARK,DNOR,ATOR,NATOMDEL)

      NDEL=NATOMDEL-OLDNATOMDEL

      CALL CLEANUP(NATOM,NDEL,DNOR,ATOR,PEPDIR,X,Y,Z,NPG,PGID,PGLEN,PGTYP,PGDIR,MARK, &
                   NLOOP,LOOP,LOOPLEN,LOOPTYP,NBONDPEP,BONDPEP,BONTYP)

      OLDNATOMDEL=NATOMDEL

      CALL RMPEPBOND(NPG,PGID,PGLEN,MARK,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,DNOR,ATOR,NLOOP,LOOP,LOOPLEN,LOOPTYP)
   END DO

   DEALLOCATE(MARK)

!-------------------------------------------------
!  CLEAN UP TAILS AGAIN:

   CALL PEPCLEAN(NBONDPEP,NBONDDEL,BONDPEP,BONTYP,NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,NATOM, &
                 NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP)

   CALL LYTGLY(NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,NATOMDEL)

   NDEL=NATOMDEL-OLDNATOMDEL

   CALL PG_UPDATE(NATOM,NDEL,DNOR,ATOR,PEPDIR,X,Y,Z,NPG,PGID, &
                  PGLEN,PGTYP,PGDIR, &
                  NLOOP,LOOP,LOOPLEN,LOOPTYP,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP)

   OLDNATOMDEL=NATOMDEL

!----------------------------------------------
!  DEFINE GLYCAN BONDS AFTER DIVIDING STRANDS:

   NBONDGLY=0

   DO NR=1,NPG

      IF(PGLEN(NR)==1)THEN
         CYCLE
      END IF

      DO JR=1,PGLEN(NR)-1
         NBONDGLY=NBONDGLY+1
         BONDGLY(1,NBONDGLY)=PGID(JR,NR)
         BONDGLY(2,NBONDGLY)=PGID(JR+1,NR)
      END DO

      IF(PGTYP(NR)==0)THEN
         NBONDGLY=NBONDGLY+1
         BONDGLY(1,NBONDGLY)=PGID(PGLEN(NR),NR)
         BONDGLY(2,NBONDGLY)=PGID(1,NR)
      END IF

   END DO

!  Clean up pep-bonds:

   ALLOCATE(NEWBONDPEP(2,NBONDPEP-NBONDDEL),NEWBONTYP(NBONDPEP-NBONDDEL))

   NNEW=0

   DO NB=1,NBONDPEP
      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      NNEW=NNEW+1
      NEWBONTYP(NNEW)=BONTYP(NB)
      NEWBONDPEP(1:2,NNEW)=BONDPEP(1:2,NB)


   END DO

   IF(NNEW/=NBONDPEP-NBONDDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED BONDS'
      STOP
   END IF

   NBONDPEP=NNEW
   NBONDDEL=0

   DEALLOCATE(BONDPEP,BONTYP)

   NPEPMAX=NBONDPEP+1000

   ALLOCATE(BONDPEP(2,NPEPMAX),BONTYP(NPEPMAX))

   BONDPEP(1,1:NBONDPEP)=NEWBONDPEP(1,1:NBONDPEP)
   BONDPEP(2,1:NBONDPEP)=NEWBONDPEP(2,1:NBONDPEP)
   BONTYP(1:NBONDPEP)=NEWBONTYP(1:NBONDPEP)

   DEALLOCATE(NEWBONDPEP,NEWBONTYP)

!---------------------------------------------------------

!  SET UP THE CONDITIONS FOR DIVISION HOOP AND ROTATION MATRIX:

   ALLOCATE(XNEW(LENHOOP),YNEW(LENHOOP),ZNEW(LENHOOP))

   X0=XCEN0; Y0=YCEN0+1.0D0; Z0=ZCEN0+1.0D0

   ARG=VY+VZ

   XPRJ=X0-ARG*VX
   YPRJ=Y0-ARG*VY
   ZPRJ=Z0-ARG*VZ

   DX=XPRJ-XCEN0
   DY=YPRJ-YCEN0
   DZ=ZPRJ-ZCEN0

   DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

   XNEW(1)=XCEN0+DX*RADIUS/DIST
   YNEW(1)=YCEN0+DY*RADIUS/DIST
   ZNEW(1)=ZCEN0+DZ*RADIUS/DIST

!  ROTATIONAL MATRIX FROM AXIS

   DPHI=2*PI/LENHOOP

   SINP=SIN(DPHI)
   COSP=COS(DPHI)

   RXX=COSP+VX*VX*(1.0D0-COSP)
   RXY=VX*VY*(1.0D0-COSP)-VZ*SINP
   RXZ=VX*VZ*(1.0D0-COSP)+VY*SINP

   RYX=VY*VX*(1.0D0-COSP)+VZ*SINP
   RYY=COSP+VY*VY*(1.0D0-COSP)
   RYZ=VY*VZ*(1.0D0-COSP)-VX*SINP

   RZX=VZ*VX*(1.0D0-COSP)-VY*SINP
   RZY=VZ*VY*(1.0D0-COSP)+VX*SINP
   RZZ=COSP+VZ*VZ*(1.0D0-COSP)

   DO N=2,LENHOOP

      DX=XNEW(N-1)-XCEN0
      DY=YNEW(N-1)-YCEN0
      DZ=ZNEW(N-1)-ZCEN0

      DXNEW=RXX*DX+RXY*DY+RXZ*DZ
      DYNEW=RYX*DX+RYY*DY+RYZ*DZ
      DZNEW=RZX*DX+RZY*DY+RZZ*DZ

      XNEW(N)=XCEN0+DXNEW
      YNEW(N)=YCEN0+DYNEW
      ZNEW(N)=ZCEN0+DZNEW

   END DO

!------------------------------------------------------------

!  LOCATE THE FIRST LOOP TO START

   X0=XNEW(1); Y0=YNEW(1); Z0=ZNEW(1)

   CALL GETLOOP(X0,Y0,Z0,XCEN0,YCEN0,ZCEN0,X,Y,Z,NLOOP,LOOP,LOOPLEN,LOOPTYP,ILOOP)


   NPG=NPG+1

   PGDIR(NPG)=1

   PGTYP(NPG)=0

!-----------------------------------------------

   DO N=1,LENHOOP+1


      NATOM=NATOM+1


      IF(N<=LENHOOP)THEN
         X(NATOM)=XNEW(N)
         Y(NATOM)=YNEW(N)
         Z(NATOM)=ZNEW(N)

      ELSE
         X(NATOM)=XNEW(1)
         Y(NATOM)=YNEW(1)
         Z(NATOM)=ZNEW(1)
      END IF

      DNOR(NATOM)=1
      ATOR(NATOM)=1

      IF(N==1)THEN
         CALL RANDOM_NUMBER(R)

         IF(R>0.5D0)THEN
            PEPDIR(NATOM)=1
         ELSE
            PEPDIR(NATOM)=-1
         END IF

         PRINT*,'FIRST PEPTIDE DIRECTION',PEPDIR(NATOM)

      ELSE
         PEPDIR(NATOM)=-PEPDIR(NATOM-1)
      END IF

      IF(N<=LENHOOP)THEN
         PGLEN(NPG)=N

         PGID(N,NPG)=NATOM

         IF(N>1)THEN
            NBONDGLY=NBONDGLY+1

            BONDGLY(1,NBONDGLY)=NATOM-1
            BONDGLY(2,NBONDGLY)=NATOM

         END IF


      END IF

!------------------------------------------------------
!-----------------------------------------------
!     CLEAVING PEPTIDE BOND:

      IF(N>1)THEN

         CALL CLEAVEPEP(NATOM,X,Y,Z,DNOR,ATOR,XCEN0,YCEN0,ZCEN0,ILOOP,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP, &
                        NBONDPEP,NBONDDEL,BONDPEP,BONTYP)

      END IF

!---------------------------------------------------------
!     CROSSLINKING:

      ND=NATOM
      JD=N

      IF(N==LENHOOP+1)THEN
         ND=NATOM-LENHOOP
         JD=1
      END IF

      CALL LINKING(ND,JD,NLINK,DNOR,ATOR,PEPDIR,VX,VY,VZ,XCEN0,YCEN0,ZCEN0,RADIUS,X,Y,Z,ILOOP,LOOP,LOOPLEN)

!---------------------------
!     CROSSLINKING THE CURRENT BEAD:

      IF(NLINK>0)THEN

         CALL POSTLINK(NATOM,ND,NLINK,DNOR,ATOR,NLOOP,NLOOPDEL,ILOOP,LOOP,LOOPLEN,LOOPTYP, &
                       NPG,PGID,PGLEN,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,PEPDIR,NBONDGLY,BONDGLY)

      END IF
!--------------------------------------------
!     CROSSLINKING THE PREVIOUS BEAD:

      IF(N>1)THEN

         ND=NATOM-1

         IF(DNOR(ND)/=0.AND.ATOR(ND)/=0)THEN
            CALL LINKING(ND,JD,NLINK,DNOR,ATOR,PEPDIR,VX,VY,VZ,XCEN0,YCEN0,ZCEN0,RADIUS,X,Y,Z,ILOOP,LOOP,LOOPLEN)

            IF(NLINK>0)THEN

               CALL POSTLINK(NATOM,ND,NLINK,DNOR,ATOR,NLOOP,NLOOPDEL,ILOOP,LOOP,LOOPLEN,LOOPTYP, &
                             NPG,PGID,PGLEN,NBONDPEP,NBONDDEL,BONDPEP,BONTYP,PEPDIR,NBONDGLY,BONDGLY)


            END IF

         END IF

      END IF


   END DO

!---------------------------------

   NATOM=NATOM-1

   DO NB=1,NBONDPEP

      IF(BONTYP(NB)==2)THEN
         BONTYP(NB)=1
      END IF

   END DO

!  Clean up pep-bonds:

   ALLOCATE(NEWBONDPEP(2,NBONDPEP-NBONDDEL),NEWBONTYP(NBONDPEP-NBONDDEL))

   NNEW=0

   DO NB=1,NBONDPEP
      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      NNEW=NNEW+1
      NEWBONTYP(NNEW)=BONTYP(NB)
      NEWBONDPEP(1:2,NNEW)=BONDPEP(1:2,NB)


   END DO

   IF(NNEW/=NBONDPEP-NBONDDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED BONDS'
      STOP
   END IF

   NBONDPEP=NNEW
   NBONDDEL=0

   DEALLOCATE(BONDPEP,BONTYP)

   NPEPMAX=NBONDPEP+1000

   ALLOCATE(BONDPEP(2,NPEPMAX),BONTYP(NPEPMAX))

   BONDPEP(1,1:NBONDPEP)=NEWBONDPEP(1,1:NBONDPEP)
   BONDPEP(2,1:NBONDPEP)=NEWBONDPEP(2,1:NBONDPEP)
   BONTYP(1:NBONDPEP)=NEWBONTYP(1:NBONDPEP)

   DEALLOCATE(NEWBONDPEP,NEWBONTYP)

!-------------

   NA1=NATOM

   NA2=NATOM-LENHOOP+1

   NBONDGLY=NBONDGLY+1

   BONDGLY(1,NBONDGLY)=NA1
   BONDGLY(2,NBONDGLY)=NA2

!   MAKE THE HOOP AND DIVIDE THE ILOOP:

   CALL DIVIDELASTLOOP(NATOM,NA1,NA2,ILOOP,NLOOP,LOOP,LOOPLEN,LOOPTYP,NLOOPDEL, &
                       NBONDPEP,NBONDDEL,BONDPEP,BONTYP,NBONDGLY,BONDGLY)

!---------------------------------

!  SEAL BIG LOOPS ON THE CENTRAL PLANE:

   DO JR=1,PGLEN(NPG)

      NA=PGID(JR,NPG)

      JCOUNT=0

      DO NL=1,NLOOP
         IF(LOOPTYP(NL)==0)THEN
            CYCLE
         END IF

         DO JL=1,LOOPLEN(NL)

            IF(LOOP(JL,NL)==NA)THEN
               JCOUNT=JCOUNT+1

               LIST(JCOUNT)=NL
               EXIT
            END IF

         END DO

      END DO

      IF(JCOUNT==0)THEN
         CYCLE
      END IF

      DO J=1,JCOUNT

         ILOOP=LIST(J)

         J1=0

         J2=0

         DO JL=1,LOOPLEN(ILOOP)

            IF(DNOR(LOOP(JL,ILOOP))/=0.AND.ATOR(LOOP(JL,ILOOP))/=0.AND.PEPDIR(LOOP(JL,ILOOP))==1)THEN
               J1=J1+1
               DNORTEMP(J1)=LOOP(JL,ILOOP)
            END IF

            IF(DNOR(LOOP(JL,ILOOP))/=0.AND.ATOR(LOOP(JL,ILOOP))/=0.AND.PEPDIR(LOOP(JL,ILOOP))==-1)THEN
               J2=J2+1
               ATORTEMP(J2)=LOOP(JL,ILOOP)
            END IF
         END DO

         IF(J1==0.OR.J2==0)THEN

            CYCLE
         END IF

         DIST=100.0D0

         DO JJ1=1,J1

            N1=DNORTEMP(JJ1)

            DO JJ2=1,J2

               N2=ATORTEMP(JJ2)

               DX=X(N2)-X(N1)
               DY=Y(N2)-Y(N1)
               DZ=Z(N2)-Z(N1)

!              THEY HAVE TO POINT TO EACH OTHER:

               IF(DX*VX+DY*VY+DZ*VZ<0.00001D0)THEN

                  CYCLE
               END IF

!              THEY HAVE TO BE ON DIFFERENT STRANDS:

               JEXIT=0

               DO NR=1,NPG

                  DO JR0=1,PGLEN(NR)

                     IF(PGID(JR0,NR)==N1)THEN
                        NR0=NR
                        JEXIT=1
                        EXIT
                     END IF

                  END DO

                  IF(JEXIT==1)THEN
                     EXIT
                  END IF

               END DO

               JCYCLE=0

               DO JR0=1,PGLEN(NR0)

                  IF(PGID(JR0,NR0)==N2)THEN
                     JCYCLE=1
                     EXIT
                  END IF

               END DO

               IF(JCYCLE==1)THEN
                  CYCLE
               END IF

               IF(DIST>DX*DX+DY*DY+DZ*DZ)THEN
                  DIST=DX*DX+DY*DY+DZ*DZ
                  NL1=N1
                  NL2=N2
               END IF

            END DO

         END DO

         IF(DIST>16.0D0)THEN
            CYCLE
         END IF

         CALL DIVIDELOOP(ILOOP,NL1,NL2,NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP,NBONDPEP,BONDPEP,BONTYP,DNOR,ATOR)

      END DO

   END DO

!---------------------------------
!  CLEAN UP TAILS AGAIN:

   CALL PEPCLEAN(NBONDPEP,NBONDDEL,BONDPEP,BONTYP,NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,NATOM, &
                 NLOOP,NLOOPDEL,LOOP,LOOPLEN,LOOPTYP)

   CALL LYTGLY(NPG,PGID,PGLEN,PGTYP,DNOR,ATOR,NATOMDEL)

   NDEL=NATOMDEL-OLDNATOMDEL

   CALL PG_UPDATE(NATOM,NDEL,DNOR,ATOR,PEPDIR,X,Y,Z,NPG,PGID, &
                  PGLEN,PGTYP,PGDIR, &
                  NLOOP,LOOP,LOOPLEN,LOOPTYP,NBONDGLY,NBONDPEP,BONDGLY,BONDPEP,BONTYP)

   OLDNATOMDEL=NATOMDEL

!--------------------------------------
!  CLEAN PEPTIDE BOND LIST:

   NNEW=0

   DO NB=1,NBONDPEP
      IF(BONTYP(NB)==0)THEN
         CYCLE
      END IF

      NNEW=NNEW+1
      BONDPEP(1:2,NNEW)=BONDPEP(1:2,NB)
   END DO

   IF(NNEW/=NBONDPEP-NBONDDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED BONDS'
      STOP
   END IF

   NBONDPEP=NNEW
   NBONDDEL=0

!-----------------------------
!  CLEAN LOOP LIST:

   NNEW=0

   DO NL=1,NLOOP

      IF(LOOPTYP(NL)==0)THEN
         CYCLE
      END IF

      NNEW=NNEW+1

      LOOP(1:LOOPLEN(NL),NNEW)=LOOP(1:LOOPLEN(NL),NL)
      LOOPLEN(NNEW)=LOOPLEN(NL)
      LOOPTYP(NNEW)=LOOPTYP(NL)

   END DO

   IF(NNEW/=NLOOP-NLOOPDEL)THEN
      PRINT*,'ERROR IN COUNTING DELETED LOOPS',NNEW,NLOOP,NLOOPDEL
      STOP
   END IF

   NLOOP=NNEW
   NLOOPDEL=0

   DEALLOCATE(XNEW,YNEW,ZNEW)

!  CLEAN UP THE DIVISION HOOP

   DO JR=1,LENHOOP
      NA=PGID(JR,NPG)

      IF(DNOR(NA)==1)THEN
         DNOR(NA)=-1
      END IF

   END DO

!--------------------------------------------------
!=================================================
!  CREATE THE LEFT DAUGHTER SACCULUS:

   ALLOCATE(XNEW(NATOM),YNEW(NATOM),ZNEW(NATOM),DNORNEW(NATOM),ATORNEW(NATOM),PEPDIRNEW(NATOM),MAP(NATOM))
   ALLOCATE(NEWPGID(PGLENMAX,NPG),NEWPGLEN(NPG),NEWPGTYP(NPG),NEWPGDIR(NPG))
   ALLOCATE(NEWBONDGLY(2,NBONDGLY),NEWBONDPEP(2,NBONDPEP),NEWBONTYP(NBONDPEP))
   NEWBONTYP=1

   ALLOCATE(NEWLOOP(LOOPLENMAX,NLOOP),NEWLOOPLEN(NLOOP),NEWLOOPTYP(NLOOP))
   NEWLOOPTYP=1

   ALLOCATE(JMARK(LENHOOP+1),PARTNER(LENHOOP+1))

!  Set crosslink length:
   L_PEP=2.0D0

   NCAP=PI*RADIUS/L_PEP/2

   IF(LEFTCELL==0)THEN
      GOTO 40
   END IF

!  ACQUIRE PG IN THE LEFT HALF:

   NANEW=0

   NPGNEW=0

   NBGLYNEW=0

   NBPEPNEW=0

   NLPNEW=0

   DO NR=1,NPG

      LENGTH=PGLEN(NR)

      XPG=SUM(X(PGID(1:LENGTH,NR)))/LENGTH
      YPG=SUM(Y(PGID(1:LENGTH,NR)))/LENGTH
      ZPG=SUM(Z(PGID(1:LENGTH,NR)))/LENGTH

      IF((XPG-XCEN0)*VX+(YPG-YCEN0)*VY+(ZPG-ZCEN0)*VZ>0.0D0.AND.NR<NPG)THEN
         CYCLE
      END IF

      NPGNEW=NPGNEW+1

      NEWPGLEN(NPGNEW)=LENGTH

      NEWPGTYP(NPGNEW)=PGTYP(NR)

      NEWPGDIR(NPGNEW)=PGDIR(NR)

      DO JR=1,LENGTH

         NANEW=NANEW+1

         PEPDIRNEW(NANEW)=PEPDIR(PGID(JR,NR))

         DNORNEW(NANEW)=DNOR(PGID(JR,NR))

!        FOR THE DIVISION HOOP, ONLY BOND DIRECTION TO THE LEFT IS KEPT:
         IF(NR==NPG)THEN
            IF(PEPDIRNEW(NANEW)==1)THEN
               DNORNEW(NANEW)=-1
            END IF
         END IF

         ATORNEW(NANEW)=ATOR(PGID(JR,NR))

         MAP(PGID(JR,NR))=NANEW

         NEWPGID(JR,NPGNEW)=NANEW

         XNEW(NANEW)=X(PGID(JR,NR))

         YNEW(NANEW)=Y(PGID(JR,NR))

         ZNEW(NANEW)=Z(PGID(JR,NR))

      END DO

      IF(LENGTH<2)THEN
         CYCLE
      END IF

      NEWBONDGLY(1,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(1:LENGTH-1,NPGNEW)
      NEWBONDGLY(2,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(2:LENGTH,NPGNEW)

      NBGLYNEW=NBGLYNEW+LENGTH-1

      IF(PGTYP(NR)==0)THEN

         NBGLYNEW=NBGLYNEW+1

         NEWBONDGLY(1,NBGLYNEW)=NEWPGID(LENGTH,NPGNEW)
         NEWBONDGLY(2,NBGLYNEW)=NEWPGID(1,NPGNEW)

      END IF

   END DO

!  ACQUIRE PEPTIDE BONDS IN THE LEFT HALF:

   DO NB=1,NBONDPEP

      XB=0.5D0*(X(BONDPEP(1,NB))+X(BONDPEP(2,NB)))
      YB=0.5D0*(Y(BONDPEP(1,NB))+Y(BONDPEP(2,NB)))
      ZB=0.5D0*(Z(BONDPEP(1,NB))+Z(BONDPEP(2,NB)))

      IF((XB-XCEN0)*VX+(YB-YCEN0)*VY+(ZB-ZCEN0)*VZ>0.0D0)THEN
         CYCLE
      END IF

      NBPEPNEW=NBPEPNEW+1

      NEWBONDPEP(1:2,NBPEPNEW)=MAP(BONDPEP(1:2,NB))

   END DO


!  ACQUIRE LOOPS IN THE LEFT HALF:

   DO NL=1,NLOOP

      LENGTH=LOOPLEN(NL)

      XLOOP=SUM(X(LOOP(1:LENGTH,NL)))/LENGTH
      YLOOP=SUM(Y(LOOP(1:LENGTH,NL)))/LENGTH
      ZLOOP=SUM(Z(LOOP(1:LENGTH,NL)))/LENGTH

      IF((XLOOP-XCEN0)*VX+(YLOOP-YCEN0)*VY+(ZLOOP-ZCEN0)*VZ>0.0D0)THEN
         CYCLE
      END IF

      NLPNEW=NLPNEW+1

      NEWLOOP(1:LENGTH,NLPNEW)=MAP(LOOP(1:LENGTH,NL))

      NEWLOOPLEN(NLPNEW)=LENGTH

   END DO

!print*,'#loop on left',nlpnew

!   LENGTH=MAXVAL(NEWPGLEN(1:NPGNEW-1))

!================================================================
!  BUILD THE LARGEST CLOSED HOOP OF THE RIGHT CAP

   X0=XCEN0; Y0=YCEN0; Z0=ZCEN0

   LENGTH=LENHOOP

         NPGNEW=NPGNEW+1

         DO JR=1,LENGTH

            NANEW=NANEW+1

            NEWPGID(JR,NPGNEW)=NANEW

            DNORNEW(NANEW)=1
            ATORNEW(NANEW)=1

         END DO

         NEWPGLEN(NPGNEW)=LENGTH

         NEWPGTYP(NPGNEW)=0

         NEWPGDIR(NPGNEW)=1

         PEPDIRNEW(NEWPGID(1:LENGTH,NPGNEW))=-PEPDIRNEW(NEWPGID(1:LENGTH,NPGNEW-1))

         XNEW(NEWPGID(1:LENGTH,NPGNEW))=XNEW(NEWPGID(1:LENGTH,NPGNEW-1))+VX*L_PEP

         YNEW(NEWPGID(1:LENGTH,NPGNEW))=YNEW(NEWPGID(1:LENGTH,NPGNEW-1))+VY*L_PEP

         ZNEW(NEWPGID(1:LENGTH,NPGNEW))=ZNEW(NEWPGID(1:LENGTH,NPGNEW-1))+VZ*L_PEP

!  GLYCAN BONDS:

         NEWBONDGLY(1,NBGLYNEW+1:NBGLYNEW+LENGTH)=NEWPGID(1:LENGTH,NPGNEW)

         NEWBONDGLY(2,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(2:LENGTH,NPGNEW)

         NEWBONDGLY(2,NBGLYNEW+LENGTH)=NEWPGID(1,NPGNEW)

         NBGLYNEW=NBGLYNEW+LENGTH

!  PEP BONDS AND LOOPS:

         DO JR=1,LENGTH

            IF(PEPDIRNEW(NEWPGID(JR,NPGNEW))==1)THEN
               CYCLE
            END IF

            NBPEPNEW=NBPEPNEW+1

            NEWBONDPEP(1,NBPEPNEW)=NEWPGID(JR,NPGNEW-1)
            NEWBONDPEP(2,NBPEPNEW)=NEWPGID(JR,NPGNEW)


            DNORNEW(NEWPGID(JR,NPGNEW-1:NPGNEW))=0
            ATORNEW(NEWPGID(JR,NPGNEW-1:NPGNEW))=0

            NLPNEW=NLPNEW+1

            NEWLOOP(1,NLPNEW)=NEWPGID(JR,NPGNEW-1)

            IF(JR<LENGTH)THEN
               NEWLOOP(2,NLPNEW)=NEWPGID(JR+1,NPGNEW-1)
               NEWLOOP(5,NLPNEW)=NEWPGID(JR+1,NPGNEW)
            ELSE
               NEWLOOP(2,NLPNEW)=NEWPGID(1,NPGNEW-1)
               NEWLOOP(5,NLPNEW)=NEWPGID(1,NPGNEW)
            END IF

            IF(JR<LENGTH-1)THEN
               NEWLOOP(3,NLPNEW)=NEWPGID(JR+2,NPGNEW-1)
               NEWLOOP(4,NLPNEW)=NEWPGID(JR+2,NPGNEW)
            ELSE IF(JR==LENGTH-1)THEN
               NEWLOOP(3,NLPNEW)=NEWPGID(1,NPGNEW-1)
               NEWLOOP(4,NLPNEW)=NEWPGID(1,NPGNEW)
            ELSE
               NEWLOOP(3,NLPNEW)=NEWPGID(2,NPGNEW-1)
               NEWLOOP(4,NLPNEW)=NEWPGID(2,NPGNEW)
            END IF

            NEWLOOP(6,NLPNEW)=NEWPGID(JR,NPGNEW)

            NEWLOOPLEN(NLPNEW)=6

         END DO

!-- UPDATE CENTER OF HOOP:

         X0=X0+VX*L_PEP
         Y0=Y0+VY*L_PEP
         Z0=Z0+VZ*L_PEP

!=============================================================
!  BUILD THE REST OF THE RIGHT CAP AS HALF OF A HEMISPHERE


!  CENTER OF THE HEMISPHERE:

   X00=X0!-VX*RADIUS/2
   Y00=Y0!-VY*RADIUS/2
   Z00=Z0!-VZ*RADIUS/2

   DO NC=1,NCAP-1

      PHI=NC*PI/NCAP/2

      DIST=RADIUS*SIN(PHI)

! NEW HOOP CENTER:

      X0=X00+DIST*VX
      Y0=Y00+DIST*VY
      Z0=Z00+DIST*VZ

! NEW HOOP RADIUS:

      RAD=RADIUS*COS(PHI)

! ADD NEW HOOP:

      NPGNEW=NPGNEW+1

      NEWPGTYP(NPGNEW)=0

      NEWPGDIR(NPGNEW)=1

      LENGTH=INT((2.0D0*PI*RAD+0.000001D0)/L_G)

      IF(LENGTH<5)THEN
         PRINT*,'ERROR: HOOP IS TOO SMALL'
         STOP
      END IF

      IF(MOD(LENGTH,2)==1)THEN
         LENGTH=LENGTH-1
      END IF

      NEWPGLEN(NPGNEW)=LENGTH

! SETUP COORDINATE FOR THE FIRST UNIT ON HOOP:

      CALL RANDOM_NUMBER(R)

      JR=(LENGTH-1)*R+2


      DX=XNEW(NEWPGID(JR,NPGNEW-1))-XNEW(NEWPGID(1,NPGNEW-1))
      DY=YNEW(NEWPGID(JR,NPGNEW-1))-YNEW(NEWPGID(1,NPGNEW-1))
      DZ=ZNEW(NEWPGID(JR,NPGNEW-1))-ZNEW(NEWPGID(1,NPGNEW-1))

      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

      NANEW=NANEW+1

      XNEW(NANEW)=X0+DX*RAD/DIST
      YNEW(NANEW)=Y0+DY*RAD/DIST
      ZNEW(NANEW)=Z0+DZ*RAD/DIST

      ATORNEW(NANEW)=1
      DNORNEW(NANEW)=1

      PEPDIRNEW(NANEW)=0

      NEWPGID(1,NPGNEW)=NANEW

! SETUP ROTATION MATRIX FOR THIS HOOP:

      PHI=2*PI/LENGTH

      SINP=SIN(PHI)
      COSP=COS(PHI)

      RXX=COSP+VX*VX*(1.0D0-COSP)
      RXY=VX*VY*(1.0D0-COSP)-VZ*SINP
      RXZ=VX*VZ*(1.0D0-COSP)+VY*SINP

      RYX=VY*VX*(1.0D0-COSP)+VZ*SINP
      RYY=COSP+VY*VY*(1.0D0-COSP)
      RYZ=VY*VZ*(1.0D0-COSP)-VX*SINP

      RZX=VZ*VX*(1.0D0-COSP)-VY*SINP
      RZY=VZ*VY*(1.0D0-COSP)+VX*SINP
      RZZ=COSP+VZ*VZ*(1.0D0-COSP)

! NOW SETUP COORDINATE FOR THE REST OF THE HOOP:

      DO JR=2,LENGTH

         NANEW=NANEW+1

         ATORNEW(NANEW)=1
         DNORNEW(NANEW)=1

         PEPDIRNEW(NANEW)=0

         NEWPGID(JR,NPGNEW)=NANEW

         DX=XNEW(NANEW-1)-X0
         DY=YNEW(NANEW-1)-Y0
         DZ=ZNEW(NANEW-1)-Z0

         DXNEW=RXX*DX+RXY*DY+RXZ*DZ
         DYNEW=RYX*DX+RYY*DY+RYZ*DZ
         DZNEW=RZX*DX+RZY*DY+RZZ*DZ

         XNEW(NANEW)=X0+DXNEW
         YNEW(NANEW)=Y0+DYNEW
         ZNEW(NANEW)=Z0+DZNEW

      END DO

! ADD NEW GLY BONDS:

      NEWBONDGLY(1,NBGLYNEW+1:NBGLYNEW+LENGTH)=NEWPGID(1:LENGTH,NPGNEW)

      NEWBONDGLY(2,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(2:LENGTH,NPGNEW)

      NEWBONDGLY(2,NBGLYNEW+LENGTH)=NEWPGID(1,NPGNEW)

      NBGLYNEW=NBGLYNEW+LENGTH

! ADD NEW PEP BONDS FROM HOOP A ON THE LEFT TO HOOP B ON THE RIGHT THEN MAKE LOOPS:

      JMARK=0
      PARTNER=0

      JA1=0

      JSKIP=0

      LENB=LENGTH

      LENA=NEWPGLEN(NPGNEW-1)

      DO JB=1,LENB

         IF(JSKIP==1.AND.NC<NCAP-1)THEN
            JSKIP=0
            CYCLE
         END IF

         NB=NEWPGID(JB,NPGNEW)


         DIST=16.0D0

         JGET=0

         DO JA=1,LENA

            IF(JMARK(JA)==1)THEN
                CYCLE
            END IF

            NA=NEWPGID(JA,NPGNEW-1)

            IF(DNORNEW(NA)==0.OR.ATORNEW(NA)==0)THEN
               CYCLE
            END IF

            DX=XNEW(NA)-XNEW(NB)
            DY=YNEW(NA)-YNEW(NB)
            DZ=ZNEW(NA)-ZNEW(NB)

            IF(DIST>DX*DX+DY*DY+DZ*DZ)THEN
               DIST=DX*DX+DY*DY+DZ*DZ
               NGET=NA
               JGET=JA
            END IF

         END DO

         IF(JGET==0)THEN
            CYCLE
         END IF

         NBPEPNEW=NBPEPNEW+1

         NEWBONDPEP(1,NBPEPNEW)=NGET

         NEWBONDPEP(2,NBPEPNEW)=NB

         ATORNEW(NB)=0

         DNORNEW(NGET)=0

         PEPDIRNEW(NB)=-1
         PEPDIRNEW(NGET)=1

         PARTNER(JGET)=JB

!        MARK THE AREA ON HOOP A AS ALREADY SEARCHED:

         IF(JA1==0)THEN
            JMARK(JGET)=1
         ELSE
            DJ=JGET-JA1

            IF(DJ<0)THEN
               DJ=DJ+LENA
            END IF

            DO J=JA1,JA1+DJ

               IF(J<=LENA)THEN
                  JMARK(J)=1
               ELSE
                  JMARK(J-LENA)=1
               END IF
            END DO

         END IF

         JA1=JGET

         JSKIP=1

      END DO


! ADD LOOPS:


      DO JA=1,LENA

         IF(PARTNER(JA)==0)THEN
            CYCLE
         END IF

         DO J=JA+1,2*LENA

            IF(J>LENA)THEN
               JA1=J-LENA
            ELSE
               JA1=J
            END IF

            IF(PARTNER(JA1)/=0)THEN
               EXIT
            END IF

         END DO

         JB=PARTNER(JA)

         JB1=PARTNER(JA1)

         NLPNEW=NLPNEW+1

         LOLEN=0

         DJ=JA1-JA

         IF(DJ<0)THEN
            DJ=DJ+LENA
         END IF

         DO J=JA,JA+DJ

            LOLEN=LOLEN+1

            IF(J<=LENA)THEN
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J,NPGNEW-1)
            ELSE
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J-LENA,NPGNEW-1)
            END IF

         END DO

         DJ=JB1-JB

         IF(DJ<0)THEN
            DJ=DJ+LENB
         END IF

         DO J=JB1,JB1-DJ,-1

            LOLEN=LOLEN+1

            IF(J>0)THEN
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J,NPGNEW)
            ELSE
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J+LENB,NPGNEW)
            END IF

         END DO

         NEWLOOPLEN(NLPNEW)=LOLEN

      END DO


   END DO

!--------------------------------------------------
!  THE LAST LOOP = LAST HOOP

   NLPNEW=NLPNEW+1

   NEWLOOP(1:LENGTH,NLPNEW)=NEWPGID(1:LENGTH,NPGNEW)

   NEWLOOPLEN(NLPNEW)=LENGTH

!----------------------------------------------

!  CLEAN UP DNORNEW, ATORNEW:

   DNORNEW=-1
   ATORNEW=1

   DNORNEW(NEWBONDPEP(1,1:NBPEPNEW))=0
   ATORNEW(NEWBONDPEP(2,1:NBPEPNEW))=0

   NATOMDEL=0
   CALL LYTGLY(NPGNEW,NEWPGID,NEWPGLEN,NEWPGTYP,DNORNEW,ATORNEW,NATOMDEL)

   NDEL=NATOMDEL

   CALL PG_UPDATE(NANEW,NDEL,DNORNEW,ATORNEW,PEPDIRNEW,XNEW,YNEW,ZNEW,NPGNEW,NEWPGID, &
                  NEWPGLEN,NEWPGTYP,NEWPGDIR, &
                  NLPNEW,NEWLOOP,NEWLOOPLEN,NEWLOOPTYP,NBGLYNEW,NBPEPNEW,NEWBONDGLY,NEWBONDPEP,NEWBONTYP)

   NBOND=NBGLYNEW+NBPEPNEW

   ALLOCATE(BOND(2,NBOND))

   BOND(1,1:NBGLYNEW)=NEWBONDGLY(1,1:NBGLYNEW)
   BOND(2,1:NBGLYNEW)=NEWBONDGLY(2,1:NBGLYNEW)

   BOND(1,NBGLYNEW+1:NBOND)=NEWBONDPEP(1,1:NBPEPNEW)
   BOND(2,NBGLYNEW+1:NBOND)=NEWBONDPEP(2,1:NBPEPNEW)

!  DETERMINE 2 HOOPS SEPARATING 2 CAPS FROM THE CYLINDER:

   XCEN0=SUM(XNEW(1:NANEW))/NANEW

   XCAP1=XCEN0-1000.0D0
   XCAP2=XCEN0+1000.0D0

   DO NR=1,NPGNEW

      IF(NEWPGTYP(NR)/=0)THEN
         CYCLE
      END IF

      XPG=SUM(XNEW(NEWPGID(1:NEWPGLEN(NR),NR)))/NEWPGLEN(NR)

      IF(XPG>XCAP1.AND.XPG<XCEN0)THEN
         NEWPGCAP1=NR
         XCAP1=XPG
      END IF

      IF(XPG<XCAP2.AND.XPG>XCEN0)THEN
         NEWPGCAP2=NR
         XCAP2=XPG
      END IF

   END DO

!  CALCULATE NUMBER OF PG UNITS ON TWO CAPS:

   NEWNATOMCAP=0

   DO NR=1,NPGNEW
      IF(NEWPGTYP(NR)==0)THEN
         NEWNATOMCAP=NEWNATOMCAP+NEWPGLEN(NR)
      END IF
   END DO

!-----------------------
!  Change representation to finalize:

   NATOM=NANEW
   X(1:NATOM)=XNEW(1:NANEW)
   Y(1:NATOM)=YNEW(1:NANEW)
   Z(1:NATOM)=ZNEW(1:NANEW)

   DNOR(1:NATOM)=DNORNEW(1:NANEW)
   ATOR(1:NATOM)=ATORNEW(1:NANEW)
   PEPDIR(1:NATOM)=PEPDIRNEW(1:NANEW)

   NPG=NPGNEW

   PGLEN(1:NPG)=NEWPGLEN(1:NPG)
   PGTYP(1:NPG)=NEWPGTYP(1:NPG)
   PGDIR(1:NPG)=NEWPGDIR(1:NPG)

   DO NR=1,NPG
      PGID(1:PGLEN(NR),NR)=NEWPGID(1:PGLEN(NR),NR)
   END DO

   NBONDGLY=NBGLYNEW

   NBONDPEP=NBPEPNEW

   BONDGLY(1,1:NBONDGLY)=NEWBONDGLY(1,1:NBONDGLY)
   BONDGLY(2,1:NBONDGLY)=NEWBONDGLY(2,1:NBONDGLY)

   BONDPEP(1,1:NBONDPEP)=NEWBONDPEP(1,1:NBONDPEP)
   BONDPEP(2,1:NBONDPEP)=NEWBONDPEP(2,1:NBONDPEP)
   BONTYP(1:NBONDPEP)=NEWBONTYP(1:NBONDPEP)

   NLOOP=NLPNEW
   LOOPLEN(1:NLOOP)=NEWLOOPLEN(1:NLOOP)
   LOOPTYP(1:NLOOP)=NEWLOOPTYP(1:NLOOP)

   DO NL=1,NLOOP
      LOOP(1:LOOPLEN(NL),NL)=NEWLOOP(1:LOOPLEN(NL),NL)
   END DO

   NATOMCAP=NEWNATOMCAP

   PGCAP1=NEWPGCAP1

   PGCAP2=NEWPGCAP2

   GOTO 50

!   CALL OUTPUTS(NANEW,XNEW,YNEW,ZNEW,DNORNEW,ATORNEW,PEPDIRNEW,NPGNEW,NEWPGID,NEWPGLEN,NEWPGTYP,NEWPGDIR, &
!                NBOND,NBGLYNEW,NBPEPNEW,BOND,NEWBONDGLY,NEWBONDPEP,NEWBONTYP, &
!                NLPNEW,NEWLOOP,NEWLOOPLEN,NEWLOOPTYP,NEWNATOMCAP,NEWPGCAP1,NEWPGCAP2)


!   COMMAND='mkdir -p leftcell'
!   CALL SYSTEM(COMMAND)

!   COMMAND='mv restart.inp coorstart.inp config0000.inp PGout* leftcell'
!   CALL SYSTEM(COMMAND)

!   WRITE(*,*)
!   WRITE(*,*)'CREATED LEFT DAUGHTER CELL'
!   WRITE(*,*)
!   WRITE(*,*)'================================================'

!   DEALLOCATE(BOND)

!==============================================================================

! NOW CREATE THE RIGHT DAUGHTER CELL:

40 PRINT*,'Creating RIGHT daughter cell ...'

!  ACQUIRE PG IN THE LEFT HALF:

   NANEW=0

   NPGNEW=0

   NBGLYNEW=0

   NBPEPNEW=0

   NLPNEW=0

   DO NR=1,NPG

      LENGTH=PGLEN(NR)

      XPG=SUM(X(PGID(1:LENGTH,NR)))/LENGTH
      YPG=SUM(Y(PGID(1:LENGTH,NR)))/LENGTH
      ZPG=SUM(Z(PGID(1:LENGTH,NR)))/LENGTH

      IF((XPG-XCEN0)*VX+(YPG-YCEN0)*VY+(ZPG-ZCEN0)*VZ<0.0D0.AND.NR<NPG)THEN
         CYCLE
      END IF

      NPGNEW=NPGNEW+1

      NEWPGLEN(NPGNEW)=LENGTH

      NEWPGTYP(NPGNEW)=PGTYP(NR)

      NEWPGDIR(NPGNEW)=PGDIR(NR)

      DO JR=1,LENGTH

         NANEW=NANEW+1

         PEPDIRNEW(NANEW)=PEPDIR(PGID(JR,NR))

         DNORNEW(NANEW)=DNOR(PGID(JR,NR))
         IF(NR==NPG.AND.PEPDIRNEW(NANEW)==-1)THEN
            DNORNEW(NANEW)=-1
         END IF

         ATORNEW(NANEW)=ATOR(PGID(JR,NR))

         MAP(PGID(JR,NR))=NANEW

         NEWPGID(JR,NPGNEW)=NANEW

         XNEW(NANEW)=X(PGID(JR,NR))

         YNEW(NANEW)=Y(PGID(JR,NR))

         ZNEW(NANEW)=Z(PGID(JR,NR))

      END DO

      IF(LENGTH<2)THEN
         CYCLE
      END IF

      NEWBONDGLY(1,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(1:LENGTH-1,NPGNEW)
      NEWBONDGLY(2,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(2:LENGTH,NPGNEW)

      NBGLYNEW=NBGLYNEW+LENGTH-1

      IF(PGTYP(NR)==0)THEN

         NBGLYNEW=NBGLYNEW+1

         NEWBONDGLY(1,NBGLYNEW)=NEWPGID(LENGTH,NPGNEW)
         NEWBONDGLY(2,NBGLYNEW)=NEWPGID(1,NPGNEW)

      END IF

   END DO

!  ACQUIRE PEPTIDE BONDS IN THE RIGHT HALF:

   DO NB=1,NBONDPEP

      XB=0.5D0*(X(BONDPEP(1,NB))+X(BONDPEP(2,NB)))
      YB=0.5D0*(Y(BONDPEP(1,NB))+Y(BONDPEP(2,NB)))
      ZB=0.5D0*(Z(BONDPEP(1,NB))+Z(BONDPEP(2,NB)))

      IF((XB-XCEN0)*VX+(YB-YCEN0)*VY+(ZB-ZCEN0)*VZ<0.0D0)THEN
         CYCLE
      END IF

      NBPEPNEW=NBPEPNEW+1

      NEWBONDPEP(1:2,NBPEPNEW)=MAP(BONDPEP(1:2,NB))

   END DO

!  ACQUIRE LOOPS IN THE RIGHT HALF:

   DO NL=1,NLOOP

      LENGTH=LOOPLEN(NL)

      XLOOP=SUM(X(LOOP(1:LENGTH,NL)))/LENGTH
      YLOOP=SUM(Y(LOOP(1:LENGTH,NL)))/LENGTH
      ZLOOP=SUM(Z(LOOP(1:LENGTH,NL)))/LENGTH

      IF((XLOOP-XCEN0)*VX+(YLOOP-YCEN0)*VY+(ZLOOP-ZCEN0)*VZ<0.0D0)THEN
         CYCLE
      END IF

      NLPNEW=NLPNEW+1

      NEWLOOP(1:LENGTH,NLPNEW)=MAP(LOOP(1:LENGTH,NL))

      NEWLOOPLEN(NLPNEW)=LENGTH

   END DO

!================================================================
!  BUILD THE LARGEST CLOSED HOOPS OF THE LEFT CAP

   X0=XCEN0; Y0=YCEN0; Z0=ZCEN0

   LENGTH=NEWPGLEN(NPGNEW)

         NPGNEW=NPGNEW+1

         DO JR=1,LENGTH

            NANEW=NANEW+1

            NEWPGID(JR,NPGNEW)=NANEW

            ATORNEW(NANEW)=1
            DNORNEW(NANEW)=1

         END DO

         NEWPGLEN(NPGNEW)=LENGTH

         NEWPGTYP(NPGNEW)=0

         NEWPGDIR(NPGNEW)=1

         PEPDIRNEW(NEWPGID(1:LENGTH,NPGNEW))=-PEPDIRNEW(NEWPGID(1:LENGTH,NPGNEW-1))

         XNEW(NEWPGID(1:LENGTH,NPGNEW))=XNEW(NEWPGID(1:LENGTH,NPGNEW-1))-VX*L_PEP

         YNEW(NEWPGID(1:LENGTH,NPGNEW))=YNEW(NEWPGID(1:LENGTH,NPGNEW-1))-VY*L_PEP

         ZNEW(NEWPGID(1:LENGTH,NPGNEW))=ZNEW(NEWPGID(1:LENGTH,NPGNEW-1))-VZ*L_PEP

!  GLYCAN BONDS:

         NEWBONDGLY(1,NBGLYNEW+1:NBGLYNEW+LENGTH)=NEWPGID(1:LENGTH,NPGNEW)

         NEWBONDGLY(2,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(2:LENGTH,NPGNEW)

         NEWBONDGLY(2,NBGLYNEW+LENGTH)=NEWPGID(1,NPGNEW)

         NBGLYNEW=NBGLYNEW+LENGTH

!  PEP BONDS AND LOOPS:

         DO JR=1,LENGTH

            IF(PEPDIRNEW(NEWPGID(JR,NPGNEW))==-1)THEN
               CYCLE
            END IF

            NBPEPNEW=NBPEPNEW+1

            NEWBONDPEP(1,NBPEPNEW)=NEWPGID(JR,NPGNEW)
            NEWBONDPEP(2,NBPEPNEW)=NEWPGID(JR,NPGNEW-1)


            ATORNEW(NEWPGID(JR,NPGNEW-1:NPGNEW))=0
            DNORNEW(NEWPGID(JR,NPGNEW-1:NPGNEW))=0


            NLPNEW=NLPNEW+1

            NEWLOOP(1,NLPNEW)=NEWPGID(JR,NPGNEW)

            IF(JR<LENGTH)THEN
               NEWLOOP(2,NLPNEW)=NEWPGID(JR+1,NPGNEW)
               NEWLOOP(5,NLPNEW)=NEWPGID(JR+1,NPGNEW-1)
            ELSE
               NEWLOOP(2,NLPNEW)=NEWPGID(1,NPGNEW)
               NEWLOOP(5,NLPNEW)=NEWPGID(1,NPGNEW-1)
            END IF

            IF(JR<LENGTH-1)THEN
               NEWLOOP(3,NLPNEW)=NEWPGID(JR+2,NPGNEW)
               NEWLOOP(4,NLPNEW)=NEWPGID(JR+2,NPGNEW-1)
            ELSE IF(JR==LENGTH-1)THEN
               NEWLOOP(3,NLPNEW)=NEWPGID(1,NPGNEW)
               NEWLOOP(4,NLPNEW)=NEWPGID(1,NPGNEW-1)
            ELSE
               NEWLOOP(3,NLPNEW)=NEWPGID(2,NPGNEW)
               NEWLOOP(4,NLPNEW)=NEWPGID(2,NPGNEW-1)
            END IF

            NEWLOOP(6,NLPNEW)=NEWPGID(JR,NPGNEW-1)

            NEWLOOPLEN(NLPNEW)=6

         END DO

!-- UPDATE CENTER OF HOOP:

         X0=X0-VX*L_PEP
         Y0=Y0-VY*L_PEP
         Z0=Z0-VZ*L_PEP

!================================================================
!  BUILD THE REST OF THE LEFT CAP AS HALF OF A HEMISPHERE:


   X00=X0!+VX*RADIUS/2
   Y00=Y0!+VY*RADIUS/2
   Z00=Z0!+VZ*RADIUS/2

   DO NC=1,NCAP-1

      PHI=NC*PI/NCAP/2

      DIST=RADIUS*SIN(PHI)

! NEW HOOP CENTER:

      X0=X00-DIST*VX
      Y0=Y00-DIST*VY
      Z0=Z00-DIST*VZ


! NEW HOOP RADIUS:

      RAD=RADIUS*COS(PHI)


! ADD NEW HOOP:

      NPGNEW=NPGNEW+1

      NEWPGTYP(NPGNEW)=0

      NEWPGDIR(NPGNEW)=1

      LENGTH=INT((2.0D0*PI*RAD+0.000001D0)/L_G)

      IF(LENGTH<5)THEN
         PRINT*,'ERROR: HOOP IS TOO SMALL'
         STOP
      END IF

      IF(MOD(LENGTH,2)==1)THEN
         LENGTH=LENGTH-1
      END IF

      NEWPGLEN(NPGNEW)=LENGTH

! SETUP COORDINATE FOR THE FIRST UNIT ON HOOP:

      CALL RANDOM_NUMBER(R)

      JR=(LENGTH-1)*R+2

      DX=XNEW(NEWPGID(JR,NPGNEW-1))-XNEW(NEWPGID(1,NPGNEW-1))
      DY=YNEW(NEWPGID(JR,NPGNEW-1))-YNEW(NEWPGID(1,NPGNEW-1))
      DZ=ZNEW(NEWPGID(JR,NPGNEW-1))-ZNEW(NEWPGID(1,NPGNEW-1))

      DIST=SQRT(DX*DX+DY*DY+DZ*DZ)

      NANEW=NANEW+1

      XNEW(NANEW)=X0+DX*RAD/DIST
      YNEW(NANEW)=Y0+DY*RAD/DIST
      ZNEW(NANEW)=Z0+DZ*RAD/DIST

      DNORNEW(NANEW)=1
      ATORNEW(NANEW)=1

      PEPDIRNEW(NANEW)=0

      NEWPGID(1,NPGNEW)=NANEW

! SETUP ROTATION MATRIX FOR THIS HOOP:

      PHI=2*PI/LENGTH

      SINP=SIN(PHI)
      COSP=COS(PHI)

      RXX=COSP+VX*VX*(1.0D0-COSP)
      RXY=VX*VY*(1.0D0-COSP)-VZ*SINP
      RXZ=VX*VZ*(1.0D0-COSP)+VY*SINP

      RYX=VY*VX*(1.0D0-COSP)+VZ*SINP
      RYY=COSP+VY*VY*(1.0D0-COSP)
      RYZ=VY*VZ*(1.0D0-COSP)-VX*SINP

      RZX=VZ*VX*(1.0D0-COSP)-VY*SINP
      RZY=VZ*VY*(1.0D0-COSP)+VX*SINP
      RZZ=COSP+VZ*VZ*(1.0D0-COSP)

! NOW SETUP COORDINATE FOR THE REST OF THE HOOP:

      DO JR=2,LENGTH

         NANEW=NANEW+1

         DNORNEW(NANEW)=1
         ATORNEW(NANEW)=1

         PEPDIRNEW(NANEW)=0

         NEWPGID(JR,NPGNEW)=NANEW

         DX=XNEW(NANEW-1)-X0
         DY=YNEW(NANEW-1)-Y0
         DZ=ZNEW(NANEW-1)-Z0

         DXNEW=RXX*DX+RXY*DY+RXZ*DZ
         DYNEW=RYX*DX+RYY*DY+RYZ*DZ
         DZNEW=RZX*DX+RZY*DY+RZZ*DZ

         XNEW(NANEW)=X0+DXNEW
         YNEW(NANEW)=Y0+DYNEW
         ZNEW(NANEW)=Z0+DZNEW

      END DO
! ADD NEW GLY BONDS:

      NEWBONDGLY(1,NBGLYNEW+1:NBGLYNEW+LENGTH)=NEWPGID(1:LENGTH,NPGNEW)

      NEWBONDGLY(2,NBGLYNEW+1:NBGLYNEW+LENGTH-1)=NEWPGID(2:LENGTH,NPGNEW)

      NEWBONDGLY(2,NBGLYNEW+LENGTH)=NEWPGID(1,NPGNEW)

      NBGLYNEW=NBGLYNEW+LENGTH

! ADD NEW PEP BONDS FROM HOOP A ON THE LEFT TO HOOP B ON THE RIGHT THEN MAKE LOOPS:

      JMARK=0
      PARTNER=0

      JB1=0

      JSKIP=0

      LENA=LENGTH

      LENB=NEWPGLEN(NPGNEW-1)

      DO JA=1,LENA

         IF(JSKIP==1.AND.NC<NCAP-1)THEN
            JSKIP=0
            CYCLE
         END IF

         NA=NEWPGID(JA,NPGNEW)

         DIST=16.0D0

         JGET=0

         DO JB=1,LENB

            IF(JMARK(JB)==1)THEN
               CYCLE
            END IF

            NB=NEWPGID(JB,NPGNEW-1)

            IF(DNORNEW(NB)==0.OR.ATORNEW(NB)==0)THEN
               CYCLE
            END IF

            DX=XNEW(NA)-XNEW(NB)
            DY=YNEW(NA)-YNEW(NB)
            DZ=ZNEW(NA)-ZNEW(NB)


            IF(DIST>DX*DX+DY*DY+DZ*DZ)THEN
               DIST=DX*DX+DY*DY+DZ*DZ
               NGET=NB
               JGET=JB
            END IF

         END DO

         IF(JGET==0)THEN
            CYCLE
         END IF

         NBPEPNEW=NBPEPNEW+1

         NEWBONDPEP(1,NBPEPNEW)=NA

         NEWBONDPEP(2,NBPEPNEW)=NGET

         DNORNEW(NA)=0

         ATORNEW(NGET)=0

         PEPDIRNEW(NA)=1
         PEPDIRNEW(NGET)=-1

         PARTNER(JGET)=JA

!        MARK ALREADY SEARCHED AEAR ON HOOP B:

         IF(JB1==0)THEN
            JMARK(JGET)=1
         ELSE

            DJ=JGET-JB1

            IF(DJ<0)THEN
               DJ=DJ+LENB
            END IF

            DO J=JB1,JB1+DJ

               IF(J<=LENB)THEN
                  JMARK(J)=1
               ELSE
                  JMARK(J-LENB)=1
               END IF

            END DO

         END IF

         JB1=JGET

         JSKIP=1

      END DO

! ADD LOOPS:


      DO JB=1,LENB

         IF(PARTNER(JB)==0)THEN
            CYCLE
         END IF

         DO J=JB+1,2*LENB

            IF(J<=LENB)THEN
               JB1=J
            ELSE
               JB1=J-LENB
            END IF

            IF(PARTNER(JB1)/=0)THEN
               EXIT
            END IF

         END DO

         JA=PARTNER(JB)

         JA1=PARTNER(JB1)

         NLPNEW=NLPNEW+1

         LOLEN=0

         DJ=JA1-JA

         IF(DJ<0)THEN
            DJ=DJ+LENA
         END IF

         DO J=JA,JA+DJ

            LOLEN=LOLEN+1

            IF(J<=LENA)THEN
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J,NPGNEW)
            ELSE
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J-LENA,NPGNEW)
            END IF

         END DO

         DJ=JB1-JB

         IF(DJ<0)THEN
            DJ=DJ+LENB
         END IF

         DO J=JB1,JB1-DJ,-1

            LOLEN=LOLEN+1

            IF(J>0)THEN
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J,NPGNEW-1)
            ELSE
               NEWLOOP(LOLEN,NLPNEW)=NEWPGID(J+LENB,NPGNEW-1)
            END IF

         END DO

         NEWLOOPLEN(NLPNEW)=LOLEN


      END DO


   END DO

!--------------------------------------------------------
!  THE LAST LOOP = LAST HOOP

   NLPNEW=NLPNEW+1

   DO J=LENGTH,1,-1

      NEWLOOP(LENGTH-J+1,NLPNEW)=NEWPGID(J,NPGNEW)

   END DO

   NEWLOOPLEN(NLPNEW)=LENGTH

!----------------------------------------------

!  CLEAN UP DNORNEW,ATORNEW:

   DNORNEW=-1
   ATORNEW=1

   DNORNEW(NEWBONDPEP(1,1:NBPEPNEW))=0
   ATORNEW(NEWBONDPEP(2,1:NBPEPNEW))=0

   NATOMDEL=0
   CALL LYTGLY(NPGNEW,NEWPGID,NEWPGLEN,NEWPGTYP,DNORNEW,ATORNEW,NATOMDEL)

   NDEL=NATOMDEL

   CALL PG_UPDATE(NANEW,NDEL,DNORNEW,ATORNEW,PEPDIRNEW,XNEW,YNEW,ZNEW,NPGNEW,NEWPGID, &
                  NEWPGLEN,NEWPGTYP,NEWPGDIR, &
                  NLPNEW,NEWLOOP,NEWLOOPLEN,NEWLOOPTYP,NBGLYNEW,NBPEPNEW,NEWBONDGLY,NEWBONDPEP,NEWBONTYP)


   NBOND=NBGLYNEW+NBPEPNEW

   ALLOCATE(BOND(2,NBOND))

   BOND(1,1:NBGLYNEW)=NEWBONDGLY(1,1:NBGLYNEW)
   BOND(2,1:NBGLYNEW)=NEWBONDGLY(2,1:NBGLYNEW)

   BOND(1,NBGLYNEW+1:NBOND)=NEWBONDPEP(1,1:NBPEPNEW)
   BOND(2,NBGLYNEW+1:NBOND)=NEWBONDPEP(2,1:NBPEPNEW)

!  DETERMINE 2 HOOPS SEPARATING 2 CAPS FROM THE CYLINDER:

   XCEN0=SUM(XNEW(1:NANEW))/NANEW

   XCAP1=XCEN0-1000.0D0
   XCAP2=XCEN0+1000.0D0

   DO NR=1,NPGNEW

      IF(NEWPGTYP(NR)/=0)THEN
         CYCLE
      END IF

      XPG=SUM(XNEW(NEWPGID(1:NEWPGLEN(NR),NR)))/NEWPGLEN(NR)

      IF(XPG>XCAP1.AND.XPG<XCEN0)THEN
         NEWPGCAP1=NR
         XCAP1=XPG
      END IF

      IF(XPG<XCAP2.AND.XPG>XCEN0)THEN
         NEWPGCAP2=NR
         XCAP2=XPG
      END IF

   END DO

!-----------------------------------
!  CALCULATE NUMBER OF PG UNITS ON TWO CAPS:

   NEWNATOMCAP=0

   DO NR=1,NPGNEW
      IF(NEWPGTYP(NR)==0)THEN
         NEWNATOMCAP=NEWNATOMCAP+NEWPGLEN(NR)
      END IF
   END DO

!-----------------------------------

!  Change representation to finalize:

   NATOM=NANEW
   X(1:NATOM)=XNEW(1:NANEW)
   Y(1:NATOM)=YNEW(1:NANEW)
   Z(1:NATOM)=ZNEW(1:NANEW)

   DNOR(1:NATOM)=DNORNEW(1:NANEW)
   ATOR(1:NATOM)=ATORNEW(1:NANEW)
   PEPDIR(1:NATOM)=PEPDIRNEW(1:NANEW)

   NPG=NPGNEW

   PGLEN(1:NPG)=NEWPGLEN(1:NPG)
   PGTYP(1:NPG)=NEWPGTYP(1:NPG)
   PGDIR(1:NPG)=NEWPGDIR(1:NPG)

   DO NR=1,NPG
      PGID(1:PGLEN(NR),NR)=NEWPGID(1:PGLEN(NR),NR)
   END DO

   NBONDGLY=NBGLYNEW

   NBONDPEP=NBPEPNEW

   BONDGLY(1,1:NBONDGLY)=NEWBONDGLY(1,1:NBONDGLY)
   BONDGLY(2,1:NBONDGLY)=NEWBONDGLY(2,1:NBONDGLY)

   BONDPEP(1,1:NBONDPEP)=NEWBONDPEP(1,1:NBONDPEP)
   BONDPEP(2,1:NBONDPEP)=NEWBONDPEP(2,1:NBONDPEP)
   BONTYP(1:NBONDPEP)=NEWBONTYP(1:NBONDPEP)

   NLOOP=NLPNEW
   LOOPLEN(1:NLOOP)=NEWLOOPLEN(1:NLOOP)
   LOOPTYP(1:NLOOP)=NEWLOOPTYP(1:NLOOP)

   DO NL=1,NLOOP
      LOOP(1:LOOPLEN(NL),NL)=NEWLOOP(1:LOOPLEN(NL),NL)
   END DO

   NATOMCAP=NEWNATOMCAP

   PGCAP1=NEWPGCAP1

   PGCAP2=NEWPGCAP2

   GOTO 50



!==========================================
!==========================================
!==========================================
!==========================================

END
