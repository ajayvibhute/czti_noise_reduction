c radec2CZTI.f v1.0
c
c Code to compute Astrosat CZT imager camera coordinates
c (theta_x,theta_y) or (theta,phi)
c for a source of known RA, Dec, given the attitude of
c the spacecraft, specified as one of the following
c options:
c 1. Roll_RA, Roll_Dec and ROLL_ROT
c 2. Q_SAT (4 components of attitude quaternion)
c 3. RA, Dec of all three principal axes: Roll, Pitch, Yaw
c
c compile as: f77 -o radec2CZTI radec2CZTI.f
c
c Dipankar Bhattacharya 7 March 2016
c based on RPY2CZTItxty.f (DB, 7 Apr 2015)
c
c $Id: radec2CZTI.f,v 1.1 2016/03/07 16:06:42 dipankar Exp $
c $Log: radec2CZTI.f,v $
c Revision 1.1  2016/03/07 16:06:42  dipankar
c Initial revision
c
c
	subroutine fcompute_tx_ty(RAroll,Decroll,roll_rot,RAs,Decs,tx,ty)

	implicit none
	real*8  cztipcsx(3),cztipcsy(3),cztipcsz(3)
	real*8  roll(3),pitch(3),yaw(3)
	real*8  radec_r(2),radec_p(2),radec_y(2)
	real*8  src(3),srcmid(3),srcfin(3)
        real*8  d2r,dotpnorm,anorm
	real*8  RAroll,Decroll,RApitch,Decpitch,RAyaw,Decyaw
	real*8  RAs,Decs,tx,ty
c	real*8  RAsrc,Decsrc
c	real*8  thetax,thetay
	real*8  theta,phi
	real*8  roll_rot,q_sat(4)
	integer i
	
	data cztipcsx/0.0d0,0.0d0,-1.0d0/
	data cztipcsy/0.0d0,1.0d0,0.0d0/
	data cztipcsz/1.0d0,0.0d0,0.0d0/

	d2r=acos(-1.d0)/180.d0

	RAs=d2r*RAs
	Decs=d2r*Decs

	call roll_rot2rpy(RAroll,Decroll,roll_rot,radec_r,radec_p,radec_y)

	RAroll=radec_r(1)
	Decroll=radec_r(2)
	RApitch=radec_p(1)
	Decpitch=radec_p(2)
	RAyaw=radec_y(1)
	Decyaw=radec_y(2)

	RAroll=d2r*RAroll
	Decroll=d2r*Decroll
	RApitch=d2r*RApitch
	Decpitch=d2r*Decpitch
	RAyaw=d2r*RAyaw
	Decyaw=d2r*Decyaw

	roll(1)=cos(Decroll)*cos(RAroll)
	roll(2)=cos(Decroll)*sin(RAroll)
	roll(3)=sin(Decroll)

	pitch(1)=cos(Decpitch)*cos(RApitch)
	pitch(2)=cos(Decpitch)*sin(RApitch)
	pitch(3)=sin(Decpitch)

	yaw(1)=cos(Decyaw)*cos(RAyaw)
	yaw(2)=cos(Decyaw)*sin(RAyaw)
	yaw(3)=sin(Decyaw)


	src(1)=cos(Decs)*cos(RAs)
	src(2)=cos(Decs)*sin(RAs)
	src(3)=sin(Decs)

	srcmid(1)=dotpnorm(roll,src)
	srcmid(2)=dotpnorm(pitch,src)
	srcmid(3)=dotpnorm(yaw,src)

	srcfin(1)=dotpnorm(cztipcsx,srcmid)
	srcfin(2)=dotpnorm(cztipcsy,srcmid)
	srcfin(3)=dotpnorm(cztipcsz,srcmid)

	tx=atan2(srcfin(1),srcfin(3))/d2r
	ty=atan2(srcfin(2),srcfin(3))/d2r

	if(abs(tx).lt.1.d-10) tx=0.0
	if(abs(ty).lt.1.d-10) ty=0.0

	anorm=0.d0

	do i=1,3
	   anorm=anorm+srcfin(i)*srcfin(i)
	end do

	anorm=sqrt(anorm)
	theta=acos(srcfin(3)/anorm)/d2r
	phi=atan2(srcfin(2),srcfin(1))/d2r
	if (phi.lt.0.d0) phi=phi+360.d0

	return
	end


	function dotpnorm(v1,v2)
	implicit none
	real*8 dotpnorm,v1(3),v2(3)
	real*8 v1mag,v2mag

	v1mag=sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
	v2mag=sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))

	if((v1mag.le.0.d0).or.(v2mag.le.0.d0))then
           print *,"null vector?"
           dotpnorm=0.0d0
	else
	   dotpnorm=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
	   dotpnorm=dotpnorm/(v1mag*v2mag)
	endif

	return
	end

c conversions.f
c
c A collection of subroutines to transform the representation
c of Astrosat attitude between different systems:
c (Quaternion), (RA, Dec, rot of roll axis) or
c (RA, Dec of the three axes: roll, pitch, yaw) 
c
c quat2rpy
c roll_rot2rpy
c rpy2roll_rot
c rpy2quat
c
c roll_rot2rpy and rpy2quat can be called in sequence to
c achieve roll_rot to quaternion transformation.  The inverse
c can be done by combining quat2rpy and rpy2roll_rot.
c
c Dipankar Bhattacharya 7 March 2016 
c
	subroutine quat2rpy(Q,roll,pitch,yaw)
	implicit none
	real*8 Q(4),roll(2),pitch(2),yaw(2)
	real*8 q11,q12,q13,q14
	real*8 q22,q23,q24
	real*8 q33,q34,q44
	real*8 rmat(3,3)
	real*8 d2r
c	integer i,j

	d2r=acos(-1.0d0)/180.0

	q44=Q(4)*Q(4)
	q11=Q(1)*Q(1)
	q22=Q(2)*Q(2)
	q33=Q(3)*Q(3)
	q12=Q(1)*Q(2)
	q13=Q(1)*Q(3)
	q23=Q(2)*Q(3)
	q34=Q(3)*Q(4)
	q24=Q(2)*Q(4)
	q14=Q(1)*Q(4)

	rmat(1,1)=q44+q11-q22-q33
	rmat(1,2)=2.0d0*(q12+q34)
	rmat(1,3)=2.0d0*(q13-q24)
	rmat(2,1)=2.0d0*(q12-q34)
	rmat(2,2)=q44-q11+q22-q33
	rmat(2,3)=2.0d0*(q23+q14)
	rmat(3,1)=2.0d0*(q13+q24)
	rmat(3,2)=2.0d0*(q23-q14)
	rmat(3,3)=q44-q11-q22+q33

c	do i=1,3
c	   print *,(rmat(i,j),j=1,3)
c	end do

        yaw(2)=asin(rmat(1,3))/d2r
        if (abs(rmat(1,1)).lt.1.d-15.and.
     c            abs(rmat(1,2)).lt.1.d-15) then
           yaw(1)=0.d0
        else
           yaw(1)=atan2(rmat(1,2),rmat(1,1))/d2r
           if (yaw(1).lt.0.d0) yaw(1)=yaw(1)+360.0d0
        endif

        roll(2)=asin(rmat(2,3))/d2r
        if (abs(rmat(2,1)).lt.1.d-15.and.
     c            abs(rmat(2,2)).lt.1.d-15) then
           roll(1)=0.d0
        else
           roll(1)=atan2(rmat(2,2),rmat(2,1))/d2r
           if (roll(1).lt.0.d0) roll(1)=roll(1)+360.0d0
        endif

        pitch(2)=asin(rmat(3,3))/d2r
        if (abs(rmat(3,1)).lt.1.d-15.and.
     c            abs(rmat(3,2)).lt.1.d-15) then
           pitch(1)=0.d0
        else
           pitch(1)=atan2(rmat(3,2),rmat(3,1))/d2r
           if (pitch(1).lt.0.d0) pitch(1)=pitch(1)+360.0d0
        endif

	return
	end

	subroutine roll_rot2rpy
     c     (roll_RA,roll_Dec,roll_Rot,roll,pitch,yaw)
	implicit none
	real*8 roll_RA,roll_Dec,roll_Rot
	real*8 roll(2),pitch(2),yaw(2)
	real*8 d2r,ra,dec,rot
	real*8 alp,amp,anp,temp
	real*8 alr,amr,anr
	real*8 aly,amy,any

	d2r=acos(-1.d0)/180.d0

	ra=roll_RA*d2r
	dec=roll_Dec*d2r
	rot=roll_Rot*d2r

	roll(1)=roll_RA
	roll(2)=roll_Dec
	alr=cos(ra)*cos(dec)
	amr=sin(ra)*cos(dec)
	anr=sin(dec)

	amp=-cos(rot)*sin(dec)*sin(ra)-sin(rot)*cos(ra)
	alp=-cos(rot)*sin(dec)*cos(ra)+sin(rot)*sin(ra)
	anp=sqrt(1.d0-alp*alp-amp*amp)
	temp=(alp*alr+amp*amr)
	anp=-sign(anp,temp)
	if (dec.lt.0.d0) anp=-anp

        pitch(2)=asin(anp)/d2r
        if (abs(alp).lt.1.d-15.and.
     c            abs(amp).lt.1.d-15) then
           pitch(1)=0.d0
        else
           pitch(1)=atan2(amp,alp)/d2r
           if (pitch(1).lt.0.d0) pitch(1)=pitch(1)+360.0d0
        endif

	aly=amr*anp-anr*amp
	amy=anr*alp-alr*anp
	any=alr*amp-amr*alp

        yaw(2)=asin(any)/d2r
        if (abs(aly).lt.1.d-15.and.
     c            abs(amy).lt.1.d-15) then
           yaw(1)=0.d0
        else
           yaw(1)=atan2(amy,aly)/d2r
           if (yaw(1).lt.0.d0) yaw(1)=yaw(1)+360.0d0
        endif

	return
	end

	subroutine rpy2roll_rot(roll,pitch,yaw,roll_rot)
	implicit none
	real*8 roll(2),pitch(2),yaw(2),roll_rot
	real*8 d2r,ra,dec,ct,st
	real*8 rmat(3,3)

	d2r=acos(-1.d0)/180.d0
	
	ra=yaw(1)*d2r
	dec=yaw(2)*d2r
	rmat(1,1)=cos(ra)*cos(dec)
	rmat(1,2)=sin(ra)*cos(dec)
	rmat(1,3)=sin(dec)

	ra=roll(1)*d2r
	dec=roll(2)*d2r
	rmat(2,1)=cos(ra)*cos(dec)
	rmat(2,2)=sin(ra)*cos(dec)
	rmat(2,3)=sin(dec)

	ra=pitch(1)*d2r
	dec=pitch(2)*d2r
	rmat(3,1)=cos(ra)*cos(dec)
	rmat(3,2)=sin(ra)*cos(dec)
	rmat(3,3)=sin(dec)

	ct=rmat(3,3)*cos(roll(2)*d2r)-sin(roll(2)*d2r)*
     c     (rmat(3,1)*cos(roll(1)*d2r)+rmat(3,2)*sin(roll(1)*d2r))
	st=rmat(3,2)*cos(roll(1)*d2r)-rmat(3,1)*sin(roll(1)*d2r)
	roll_rot=-atan2(st,ct)/d2r

	return
	end

	subroutine rpy2quat(roll,pitch,yaw,Q)
        implicit none
	real*8 roll(2),pitch(2),yaw(2),Q(4)
        real*8 R(3,3)
	real*8 d2r,ra,dec
        real*8 w2_1,w2_2,w2_3,w2_4,w

	d2r=acos(-1.d0)/180.d0
	
	ra=yaw(1)*d2r
	dec=yaw(2)*d2r
	R(1,1)=cos(ra)*cos(dec)
	R(1,2)=sin(ra)*cos(dec)
	R(1,3)=sin(dec)

	ra=roll(1)*d2r
	dec=roll(2)*d2r
	R(2,1)=cos(ra)*cos(dec)
	R(2,2)=sin(ra)*cos(dec)
	R(2,3)=sin(dec)

	ra=pitch(1)*d2r
	dec=pitch(2)*d2r
	R(3,1)=cos(ra)*cos(dec)
	R(3,2)=sin(ra)*cos(dec)
	R(3,3)=sin(dec)

        w2_1=1.0d0+R(1,1)+R(2,2)+R(3,3)
        w2_2=1.0d0+R(1,1)-R(2,2)-R(3,3)
        w2_3=1.0d0-R(1,1)+R(2,2)-R(3,3)
        w2_4=1.0d0-R(1,1)-R(2,2)+R(3,3)

        if (w2_1.gt.0.0d0) then
           w=sqrt(w2_1)
           Q(1)=(R(2,3)-R(3,2))/(2.0d0*w)
           Q(2)=(R(3,1)-R(1,3))/(2.0d0*w)
           Q(3)=(R(1,2)-R(2,1))/(2.0d0*w)
           Q(4)=w/2.0d0
           return
        endif
        if (w2_2.gt.0.0d0) then
           w=sqrt(w2_2)
           Q(1)=w/2.0d0
           Q(2)=(R(1,2)+R(2,1))/(2.0d0*w)
           Q(3)=(R(3,1)+R(1,3))/(2.0d0*w)
           Q(4)=(R(2,3)-R(3,2))/(2.0d0*w)
           return
        endif
        if (w2_3.gt.0.0d0) then
           w=sqrt(w2_3)
           Q(1)=(R(1,2)+R(2,1))/(2.0d0*w)
           Q(2)=w/2.0d0
           Q(3)=(R(2,3)+R(3,2))/(2.0d0*w)
           Q(4)=(R(3,1)-R(1,3))/(2.0d0*w)
           return
        endif
        if (w2_4.gt.0.0d0) then
           w=sqrt(w2_4)
           Q(1)=(R(3,1)+R(1,3))/(2.0d0*w)
           Q(2)=(R(2,3)+R(3,2))/(2.0d0*w)
           Q(3)=w/2.0d0
           Q(4)=(R(1,2)-R(2,1))/(2.0d0*w)
           return
        endif

        print *,"quaternion conversion error"
        return
        end
