1 -  potential_energy()
* Here at first we called the transformation functions that was created before to get the postion of centermass of each joint 
* Then to get the potential energy of each joint we use the following equation  U1 = - 1/2 * m(1) *[ 0   0  -g] * P_cm1 substituiting the mass "m" from the input 
* and define g  = 9.8 inside the function 
* Then we get the sum of all the U's and assign it to U_i 
 function U_i = potential_energy(q, m)
   
   T12 = transformation_func(q(1) , 0 , 5 , 90) ;
   T23 = transformation_func(q(2) , 0 , 5 , 0) ;
   T34 = transformation_func(q(3) , 0 , 0 , 90) ;
   T45 = transformation_func(q(4) , 5+5 , 0 , -90) ;
   T56 = transformation_func(q(5) , 0 , 0 , 90) ;
   T67 = transformation_func(q(6) , 5+5 , 0 , 0) ;
   
  g = 9.8 ;
  
   P_cm1 = [(T12(1,4))/2 ; (T12(2,4))/2 ; (T12(3,4))/2]  
    U1 = - 1/2 * m(1) *[ 0   0  -g] * P_cm1;  
   P_cm2 = [T12(1,4)+(T23(1,4))/2 ; T12(2,4)+(T23(2,4))/2 ; T12(3,4)+(T23(3,4))/2] 
    U2 = - 1/2 * m(2) *[ 0   0  -g] * P_cm2;  
   P_cm3 = [T12(1,4)+T23(1,4)+(T34(1,4))/2 ; T12(2,4)+T23(2,4)+(T34(2,4))/2  ; T12(3,4)+T23(3,4)+(T34(3,4))/2 ] 
    U3 = - 1/2 * m(3) *[ 0   0  -g] * P_cm3;
   P_cm4 = [T12(1,4)+T23(1,4)+T34(1,4)+(T45(1,4))/2 ; T12(2,4)+T23(2,4)+T34(2,4)+(T45(2,4))/2 ; T12(3,4)+T23(3,4)+T34(3,4)+(T45(3,4))/2] 
    U4 = - 1/2 * m(4) *[ 0   0  -g] * P_cm4;
   P_cm5 =  [T12(1,4)+T23(1,4)+T34(1,4)+T45(1,4)+(T56(1,4))/2 ; T12(2,4)+T23(2,4)+T34(2,4)+T45(2,4)+(T56(2,4))/2 ; T12(3,4)+T23(3,4)+T34(3,4)+T45(3,4)+(T56(3,4))/2]   
    U5 = - 1/2 * m(5) *[ 0   0  -g] * P_cm5;
   P_cm6 =[T12(1,4)+T23(1,4)+T34(1,4)+T45(1,4)+T56(1,4)+(T67(1,4))/2 ; T12(2,4)+T23(2,4)+T34(2,4)+T45(2,4)+T56(2,4)+(T67(2,4))/2 ; T12(3,4)+T23(3,4)+T34(3,4)+T45(3,4)+T56(3,4)+(T67(3,4))/2]   
   U6 = - 1/2 * m(6) *[ 0   0  -g] * P_cm6;
    
    U = [U1 ; U2; U3 ; U4 ; U5 ; U6]
    
    U_i = U1+U2+U3+U4+U5+U6
 
  end

2 - kinemtic_energy() 
* first we get the P_cm1 with same method optained in the previuos function 
* After that to get the time derivative we replace every q with q(t) by sub in the P_cm1 
* Then to get the velocity we diff with respect to time then to get the qdot and double_dot 
  we replace every  diff(q1_dot(t),t)  with double_dot and  the same method used earlier   
* and the same for diff(q4_dot(t),t) 
* then we use  T1 = (1/2 * m(1) * (V_cm1') * V_cm1) + (1/2 * (W1') * I(1)* W1) to get the 
 KINMETIC ENERGY THEN WE SUMN ALL OF THEM TO GET TH FINAL VALUE OF T_i


function T_i = kinemtic_energy(q, qdot, m, I)
   T12 = transformation_func(q(1) , 0 , 5 , 90) ;
   T23 = transformation_func(q(2) , 0 , 5 , 0) ;
   T34 = transformation_func(q(3) , 0 , 0 , 90) ;
   T45 = transformation_func(q(4) , 5+5 , 0 , -90) ;
   T56 = transformation_func(q(5) , 0 , 0 , 90) ;
   T67 = transformation_func(q(6) , 5+5 , 0 , 0) ;
   
 syms q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t) t q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot
 
   P_cm1 = [(T12(1,4))/2 ; (T12(2,4))/2 ; (T12(3,4))/2]  
   p1 = subs(P_cm1,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  
   V_cm1 = (diff(p1,t)) 
     subs(V_cm1,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
   
   W1 = [ 0 ; 0 ; qdot(1)]
   T1 = (1/2 * m(1) * (V_cm1') * V_cm1) + (1/2 * (W1') * I(1)* W1)
   P_cm2 =  [T12(1,4)+(T23(1,4))/2 ; T12(2,4)+(T23(2,4))/2 ; T12(3,4)+(T23(3,4))/2] 
   p2 = subs(P_cm2,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  
   V_cm2 = (diff(p2,t))
    subs(V_cm2,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
   
   
   W2 = W1 + (T12(1:3 , 1:3)*[0 ; 0 ; qdot(2)])
   T2 = (1/2 * m(2) * (V_cm2') * V_cm2) + (1/2 * (W2') * I(2)* W2)
   P_cm3 = [T12(1,4)+T23(1,4)+(T34(1,4))/2 ; T12(2,4)+T23(2,4)+(T34(2,4))/2  ; T12(3,4)+T23(3,4)+(T34(3,4))/2 ] 
  p3= subs(P_cm3,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  
   V_cm3 = (diff(p3,t))
     subs(V_cm3,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
   
   W3 = W2 + (T23(1:3 , 1:3)*[0 ; 0 ; qdot(3)])
   T3 = (1/2 * m(3) * (V_cm3') * V_cm3) + (1/2 * (W3') * I(3)* W3)
   P_cm4 = [T12(1,4)+T23(1,4)+T34(1,4)+(T45(1,4))/2 ; T12(2,4)+T23(2,4)+T34(2,4)+(T45(2,4))/2 ; T12(3,4)+T23(3,4)+T34(3,4)+(T45(3,4))/2] 
   p4= subs(P_cm4,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  
   V_cm4 = (diff(p4,t))
    subs(V_cm4,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
   W4 = W3 + (T34(1:3 , 1:3)*[0 ; 0 ; qdot(4)])
   T4 = (1/2 * m(4) * (V_cm4') * V_cm4) + (1/2 * (W4') * I(4)* W4)
   P_cm5 =  [T12(1,4)+T23(1,4)+T34(1,4)+T45(1,4)+(T56(1,4))/2 ; T12(2,4)+T23(2,4)+T34(2,4)+T45(2,4)+(T56(2,4))/2 ; T12(3,4)+T23(3,4)+T34(3,4)+T45(3,4)+(T56(3,4))/2]   
    p5= subs(P_cm5,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
 
   V_cm5 = (diff(p5,t))
    subs(V_cm5,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
   
   W5 = W4 + (T45(1:3 , 1:3)*[0 ; 0 ; qdot(5)])
   T5 = (1/2 * m(5) * (V_cm5') * V_cm5) + (1/2 * (W5') * I(5)* W5)
   P_cm6 =  [T12(1,4)+T23(1,4)+T34(1,4)+T45(1,4)+T56(1,4)+(T67(1,4))/2 ; T12(2,4)+T23(2,4)+T34(2,4)+T45(2,4)+T56(2,4)+(T67(2,4))/2 ; T12(3,4)+T23(3,4)+T34(3,4)+T45(3,4)+T56(3,4)+(T67(3,4))/2]
   p6= subs(P_cm6,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
 
   V_cm6 = (diff(p6,t))
     subs(V_cm6,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
   
   W6 = W5 + (T56(1:3 , 1:3)*[0 ; 0 ; qdot(6)])
   T6 = (1/2 * m(6) *(V_cm6') * V_cm6) + (1/2 * (W6') * I(6)* W6)
   
   T = [ T1 ; T2 ; T3 ; T4 ; T5 ; T6]
   T_i = T1 + T2 + T3 + T4 + T5 + T6
    
  
  end


3 - Lagrange()  
* In this function we defined the I metrcis earlier in the top side of the code and we took it as inputs in here 
* then we called the T function and the U function created earlier the we sub T from U to get L  

  function L = Lagrange() 
  syms q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot  I1 I2 I3 I4 I5 I6  m_1 m_2 m_3 m_4 m_5 m_6
  T =potential_energy([q1;q2;q3;q4;q5;q6], [m_1;m_2;m_3;m_4;m_5;m_6]);
  U = kinemtic_energy([q1;q2;q3;q4;q5;q6], [q1_dot ;q2_dot ;q3_dot ;q4_dot ;q5_dot ;q6_dot],[m_1;m_2;m_3;m_4;m_5;m_6],[I1;I2;I3;I4;I5;I6]);
  L = T - U
  
  end

4 -  Torques()
* we take the L = Lagrange()  as an input and the qs and qdot 
* then we defined them from its vector form to get the scaler values 
we use them to diff L with each q and qdot 
* then we get the diff with each qdot and use the sub method we used earlier to diff with respect to time 
* and subtract the from the diff wrt to q 
* then we get th coeeftint of each qdot and qdouble_dot with coeff() for each Taw to get th M and C and B 
and subistitute them to creat  the equation of motion 
* How ever The equation of motion did not appear her because of some reasons which are the matrix exceeded or a non linear problem or non consentent matrix 

function T = Torques(L, q, qdot)
  q1_dot = qdot(1)
  q2_dot = qdot(2)
  q3_dot = qdot(3)
  q4_dot = qdot(4)
  q5_dot = qdot(5)
  q6_dot = qdot(6)
   q1 = q(1)
  q2 = q(2)
  q3 = q(3)
  q4 = q(4)
  q5 = q(5)
  q6 = q(6)
 t1 = diff(L,q1_dot) 
 t2 = diff(L,q2_dot)
 t3 = diff(L,q3_dot)
 t4 = diff(L,q4_dot)
 t5 = diff(L,q5_dot)
 t6 = diff(L,q6_dot)
 t7 = diff(L,q1) 
 t8 = diff(L,q2)
 t9 = diff(L,q3)
 t10= diff(L,q4)
 t11= diff(L,q5)
 t12= diff(L,q6)
syms  q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t) t q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot
 tt1 = subs(t1,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  tt2 = subs(t2,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot },{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  tt3 = subs(t3,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot},{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  tt4 =subs(t4,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot},{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  tt5 = subs(t5,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot},{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})
  tt6 = subs(t6,{q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot},{q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q1_dot(t) q2_dot(t) q3_dot(t) q4_dot(t) q5_dot(t) q6_dot(t)})

 Taw1 = (diff(tt1)./diff(t)) - t7
 Taw2 = (diff(tt2)./diff(t)) - t8 
 Taw3 = (diff(tt3)./diff(t)) - t9
 Taw4 = (diff(tt4)./diff(t)) - t10
 Taw5 = (diff(tt5)./diff(t)) - t11
 Taw6 = (diff(tt6)./diff(t)) - t12 
 
 Taw11 = subs(Taw1,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
 Taw22 = subs(Taw2,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t) diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
 Taw33 = subs(Taw3,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t) diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
 Taw44 = subs(Taw4,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t) diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
 Taw55 = subs(Taw5,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t) diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
 Taw66 = subs(Taw6,{ diff(q1_dot(t),t)  diff(q2_dot(t),t)  diff(q3_dot(t),t)  diff(q4_dot(t),t)  diff(q5_dot(t),t)  diff(q6_dot(t),t)  diff(q1(t),t)  diff(q2(t),t)  diff(q3(t),t)  diff(q4(t),t)  diff(q5(t),t)  diff(q6(t),t)},{q1_double_dot q2_double_dot q3_double_dot q4_double_dot q5_double_dot q6_double_dot q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot
})
taw1 = subs(Taw11,{q1_dot^2 q2_dot^2 q3_dot^2 q4_dot^2 q5_dot^2 q6_dot^2},{str2sym('q1_dot^2') str2sym('q2_dot^2') str2sym('q3_dot^2') str2sym('q4_dot^2') str2sym('q5_dot^2') str2sym('q6_dot^2')})
taw2 = subs(Taw22,{q1_dot^2 q2_dot^2 q3_dot^2 q4_dot^2 q5_dot^2 q6_dot^2},{str2sym('q1_dot^2') str2sym('q2_dot^2') str2sym('q3_dot^2') str2sym('q4_dot^2') str2sym('q5_dot^2') str2sym('q6_dot^2')})
taw3 = subs(Taw33,{q1_dot^2 q2_dot^2 q3_dot^2 q4_dot^2 q5_dot^2 q6_dot^2},{str2sym('q1_dot^2') str2sym('q2_dot^2') str2sym('q3_dot^2') str2sym('q4_dot^2') str2sym('q5_dot^2') str2sym('q6_dot^2')})
taw4 = subs(Taw44,{q1_dot^2 q2_dot^2 q3_dot^2 q4_dot^2 q5_dot^2 q6_dot^2},{str2sym('q1_dot^2') str2sym('q2_dot^2') str2sym('q3_dot^2') str2sym('q4_dot^2') str2sym('q5_dot^2') str2sym('q6_dot^2')})
taw5 = subs(Taw55,{q1_dot^2 q2_dot^2 q3_dot^2 q4_dot^2 q5_dot^2 q6_dot^2},{str2sym('q1_dot^2') str2sym('q2_dot^2') str2sym('q3_dot^2') str2sym('q4_dot^2') str2sym('q5_dot^2') str2sym('q6_dot^2')})
taw6 = subs(Taw66,{q1_dot^2 q2_dot^2 q3_dot^2 q4_dot^2 q5_dot^2 q6_dot^2},{str2sym('q1_dot^2') str2sym('q2_dot^2') str2sym('q3_dot^2') str2sym('q4_dot^2') str2sym('q5_dot^2') str2sym('q6_dot^2')})


M1 = [ coeffs(taw1,q1_double_dot) coeffs(taw2,q1_double_dot)  coeffs(taw3,q1_double_dot)  coeffs(taw4,q1_double_dot)  coeffs(taw5,q1_double_dot) coeffs(taw6,q1_double_dot)]
M2 = [ coeffs(taw1,q2_double_dot) coeffs(taw2,q2_double_dot)  coeffs(taw3,q2_double_dot)  coeffs(taw4,q2_double_dot)  coeffs(taw5,q2_double_dot) coeffs(taw6,q2_double_dot)]
M3 = [ coeffs(taw1,q3_double_dot) coeffs(taw2,q3_double_dot)  coeffs(taw3,q3_double_dot)  coeffs(taw4,q3_double_dot)  coeffs(taw5,q3_double_dot) coeffs(taw6,q3_double_dot)]
M4 = [ coeffs(taw1,q4_double_dot) coeffs(taw2,q4_double_dot)  coeffs(taw3,q4_double_dot)  coeffs(taw4,q4_double_dot)  coeffs(taw5,q4_double_dot) coeffs(taw6,q4_double_dot)]
M5 = [ coeffs(taw1,q5_double_dot) coeffs(taw2,q5_double_dot)  coeffs(taw3,q5_double_dot)  coeffs(taw4,q5_double_dot)  coeffs(taw5,q5_double_dot) coeffs(taw6,q5_double_dot)]
M6 = [ coeffs(taw1,q6_double_dot) coeffs(taw2,q6_double_dot)  coeffs(taw3,q6_double_dot)  coeffs(taw4,q6_double_dot)  coeffs(taw5,q6_double_dot) coeffs(taw6,q6_double_dot)]
B1 = [ coeffs(taw1,q1_dot) coeffs(taw2,q1_dot)  coeffs(taw3,q1_dot)  coeffs(taw4,q1_dot)  coeffs(taw5,q1_dot) coeffs(taw6,q1_dot)]
B2 = [ coeffs(taw1,q2_dot) coeffs(taw2,q2_dot)  coeffs(taw3,q2_dot)  coeffs(taw4,q2_dot)  coeffs(taw5,q2_dot) coeffs(taw6,q2_dot)]
B3 = [ coeffs(taw1,q3_dot) coeffs(taw2,q3_dot)  coeffs(taw3,q3_dot)  coeffs(taw4,q3_dot)  coeffs(taw5,q3_dot) coeffs(taw6,q3_dot)]
B4 = [ coeffs(taw1,q4_dot) coeffs(taw2,q4_dot)  coeffs(taw3,q4_dot)  coeffs(taw4,q4_dot)  coeffs(taw5,q4_dot) coeffs(taw6,q4_dot)] 
B5 = [ coeffs(taw1,q5_dot) coeffs(taw2,q5_dot)  coeffs(taw3,q5_dot)  coeffs(taw4,q5_dot)  coeffs(taw5,q5_dot) coeffs(taw6,q5_dot)]
B6 = [ coeffs(taw1,q6_dot) coeffs(taw2,q6_dot)  coeffs(taw3,q6_dot)  coeffs(taw4,q6_dot)  coeffs(taw5,q6_dot) coeffs(taw6,q6_dot)] 
% C1 = [ coeffs(taw1,str2sym('q1_dot^2')) coeffs(taw2,str2sym('q1_dot^2'))  coeffs(taw3,str2sym('q1_dot^2'))  coeffs(taw4,str2sym('q1_dot^2'))  coeffs(taw5,str2sym('q1_dot^2')) coeffs(taw6,str2sym('q1_dot^2'))]
% C2 = [ coeffs(taw1,str2sym('q2_dot^2')) coeffs(taw2,str2sym('q2_dot^2'))  coeffs(taw3,str2sym('q2_dot^2'))  coeffs(taw4,str2sym('q2_dot^2'))  coeffs(taw5,str2sym('q2_dot^2')) coeffs(taw6,str2sym('q2_dot^2'))]
% C3 = [ coeffs(taw1,str2sym('q3_dot^2')) coeffs(taw2,str2sym('q3_dot^2'))  coeffs(taw3,str2sym('q3_dot^2'))  coeffs(taw4,str2sym('q3_dot^2'))  coeffs(taw5,str2sym('q3_dot^2')) coeffs(taw6,str2sym('q3_dot^2'))]
% C4 = [ coeffs(taw1,str2sym('q4_dot^2')) coeffs(taw2,str2sym('q4_dot^2'))  coeffs(taw3,str2sym('q4_dot^2'))  coeffs(taw4,str2sym('q4_dot^2'))  coeffs(taw5,str2sym('q4_dot^2')) coeffs(taw6,str2sym('q4_dot^2'))]
% C5 =[ coeffs(taw1,str2sym('q5_dot^2')) coeffs(taw2,str2sym('q5_dot^2'))  coeffs(taw3,str2sym('q5_dot^2'))  coeffs(taw4,str2sym('q5_dot^2'))  coeffs(taw5,str2sym('q5_dot^2')) coeffs(taw6,str2sym('q5_dot^2'))]
% C6 =[ coeffs(taw1,str2sym('q6_dot^2')) coeffs(taw2,str2sym('q6_dot^2'))  coeffs(taw3,str2sym('q6_dot^2'))  coeffs(taw4,str2sym('q6_dot^2'))  coeffs(taw5,str2sym('q6_dot^2')) coeffs(taw6,str2sym('q6_dot^2'))]  
T= [ M1  M2  M3  M4  M5  M6]
  end
  