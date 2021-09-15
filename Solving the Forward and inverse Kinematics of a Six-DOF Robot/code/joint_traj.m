function Joint_Space = joint_traj(q0, qf, qdot0, qdotf, tf)
q00=inverse_position_kinematics(5, q0, [2;0;1]);
qff=inverse_position_kinematics(5, qf, [2;0;2]);



  for i=1:tf
      r = 1
     c0=q00(r);
 c1=qdot0(r);
 c3=( qdotf(r)-(3*c1)-((2*qff(r))/tf)-(2*c0))/(5*tf^2);
 c2= (qff(r)-(c3*tf^3)-c0-(c1*tf))/tf^2;
 c = [ c0 ;c1 ;c2 ; c3 ];
     
    t=i;
 
Joint_Spacer= c(1)+(c(2)*t)+ (c(3)*t^2) + (c(4)*t^3)
Joint_Space(r,i) = Joint_Spacer
  end

  for i=1:tf
      l = 2
     c0=q00(l);
 c1=qdot0(l);
 c3=( qdotf(l)-(3*c1)-((2*qff(l))/tf)-(2*c0))/(5*tf^2);
 c2= (qff(l)-(c3*tf^3)-c0-(c1*tf))/tf^2;
 c = [ c0 ;c1 ;c2 ; c3 ];
     
    t=i;
 
Joint_Spacel= c(1)+(c(2)*t)+ (c(3)*t^2) + (c(4)*t^3)
Joint_Space(l,i) = Joint_Spacel
  end
   for i=1:tf
      e = 3
     c0=q00(e);
 c1=qdot0(e);
 c3=( qdotf(e)-(3*c1)-((2*qff(e))/tf)-(2*c0))/(5*tf^2);
 c2= (qff(e)-(c3*tf^3)-c0-(c1*tf))/tf^2;
 c = [ c0 ;c1 ;c2 ; c3 ];
     
    t=i;
 
Joint_Spacee= c(1)+(c(2)*t)+ (c(3)*t^2) + (c(4)*t^3)
Joint_Space(e,i) = Joint_Spacee
   end
   for i=1:tf
      w = 4
     c0=q00(w);
 c1=qdot0(w);
 c3=( qdotf(w)-(3*c1)-((2*qff(w))/tf)-(2*c0))/(5*tf^2);
 c2= (qff(w)-(c3*tf^3)-c0-(c1*tf))/tf^2;
 c = [ c0 ;c1 ;c2 ; c3 ];
     
    t=i;
 
Joint_Spacew= c(1)+(c(2)*t)+ (c(3)*t^2) + (c(4)*t^3)
Joint_Space(w,i) = Joint_Spacew
   end
   for i=1:tf
      o = 5
     c0=q00(o);
 c1=qdot0(o);
 c3=( qdotf(o)-(3*c1)-((2*qff(o))/tf)-(2*c0))/(5*tf^2);
 c2= (qff(o)-(c3*tf^3)-c0-(c1*tf))/tf^2;
 c = [ c0 ;c1 ;c2 ; c3 ];
     
    t=i;
 
Joint_Spaceo= c(1)+(c(2)*t)+ (c(3)*t^2) + (c(4)*t^3)
Joint_Space(o,i) = Joint_Spaceo
  end
 for i=1:tf
      z= 6
     c0=q00(z);
 c1=qdot0(z);
 c3=( qdotf(z)-(3*c1)-((2*qff(z))/tf)-(2*c0))/(5*tf^2);
 c2= (qff(z)-(c3*tf^3)-c0-(c1*tf))/tf^2;
 c = [ c0 ;c1 ;c2 ; c3 ];
     
    t=i;
 
Joint_Spacez= c(1)+(c(2)*t)+ (c(3)*t^2) + (c(4)*t^3)
Joint_Space(z,i) = Joint_Spacez
  end
 
Joint_Space= [ Joint_Space(1,1) Joint_Space(1,2) Joint_Space(1,3) Joint_Space(1,4) Joint_Space(1,5) Joint_Space(1,6) Joint_Space(1,7) Joint_Space(1,8) Joint_Space(1,9) Joint_Space(1,10);
               Joint_Space(2,1) Joint_Space(2,2) Joint_Space(2,3) Joint_Space(2,4) Joint_Space(2,5) Joint_Space(2,6) Joint_Space(2,7) Joint_Space(2,8) Joint_Space(2,9) Joint_Space(2,10);
               Joint_Space(3,1) Joint_Space(3,2) Joint_Space(3,3) Joint_Space(3,4) Joint_Space(3,5) Joint_Space(3,6) Joint_Space(3,7) Joint_Space(3,8) Joint_Space(3,9) Joint_Space(3,10);
               Joint_Space(4,1) Joint_Space(4,2) Joint_Space(4,3) Joint_Space(4,4) Joint_Space(4,5) Joint_Space(4,6) Joint_Space(4,7) Joint_Space(4,8) Joint_Space(4,9) Joint_Space(4,10);
               Joint_Space(5,1) Joint_Space(5,2) Joint_Space(5,3) Joint_Space(5,4) Joint_Space(5,5) Joint_Space(5,6) Joint_Space(5,7) Joint_Space(5,8) Joint_Space(5,9) Joint_Space(5,10);
               Joint_Space(6,1) Joint_Space(6,2) Joint_Space(6,3) Joint_Space(6,4) Joint_Space(6,5) Joint_Space(6,6) Joint_Space(6,7) Joint_Space(6,8) Joint_Space(6,9) Joint_Space(6,10)]


end
 