syms t Iyz1 Iyz2 Iyz3 Iyz4 Iyz5 Iyz6 Izx6  Izy6 Iz6 Izx5  Izy5 Iz5 Izx4  Izy4 Iz4 Izx2 Izx3  Izy3 Iz3 Izy2 Iz2 Izx1  Izy1 Iz1 Iyx6 Iy6 Iyx5 Iy5 Ixz5 Iyx4 Iy4 Iyx3 Iy3 Iyx2 Iy2 q1 q2 q3 q4 q5 q6 L0 L1 L2 L3 L4 L5 L6 T alpha theta m1 m2 m3 m4 m5 m6 Ix1 Ixy1 Ixz1 Ix2 Ixy2 Ixz2 Ix3 Ixy3 Ixz3 Ix4 Ixy4 Ixz4 Ix5 Ixy5 Ixz5   Ix6 Ixy6 Ixz6 Iyx1 Iy1 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot

I1 = [ Ix1 Ixy1 Ixz1 ; Iyx1 Iy1 Iyz1 ; Izx1  Izy1 Iz1]
I2 = [ Ix2 Ixy2 Ixz2 ; Iyx2 Iy2 Iyz2 ; Izx2  Izy2 Iz2]
I3 = [ Ix3 Ixy3 Ixz3 ; Iyx3 Iy3 Iyz3 ; Izx3  Izy3 Iz3]
I4 = [ Ix4 Ixy4 Ixz4 ; Iyx4 Iy4 Iyz4 ; Izx4  Izy4 Iz4]
I5 = [ Ix5 Ixy5 Ixz5 ; Iyx5 Iy5 Iyz5 ; Izx5  Izy5 Iz5]
I6 = [ Ix6 Ixy6 Ixz6 ; Iyx6 Iy6 Iyz6 ; Izx6  Izy6 Iz6]

% potential_energy([q1;q2;q3;q4;q5;q6], [m1;m2;m3;m4;m5;m6])
%  kinemtic_energy([q1;q2;q3;q4;q5;q6], [q1_dot ;q2_dot ;q3_dot ;q4_dot ;q5_dot ;q6_dot],[m1;m2;m3;m4;m5;m6],[I1;I2;I3;I4;I5;I6])
L = Lagrange();

T = Torques(L,[q1;q2;q3;q4;q5;q6], [q1_dot ;q2_dot ;q3_dot ;q4_dot ;q5_dot ;q6_dot])
  function T = transformation_func(theta , d , a , alpha) 
           T = [cosd(theta)  -sind(theta)*cosd(alpha) sind(theta)*sind(alpha) a*cosd(theta) ;
       sind(theta) cosd(theta)*cosd(alpha) -cosd(theta)*sind(alpha) a*sind(theta);
       0 sind(alpha) cosd(alpha) d;
       0 0 0 1]
  % X= T01 .* T12 .* T23 .* T34 .* T45 .* T56 .* T67


  end
  function X = end_effector_position(q)
  T01 = transformation_func(0 , 5 , 0 , 0) ;
 T12 = transformation_func(q(1) , 0 , 5 , 90) ;
 T23 = transformation_func(q(2) , 0 , 5 , 0) ;
 T34 = transformation_func(q(3) , 0 , 0 , 90) ;
 T45 = transformation_func(q(4) , 5+5 , 0 , -90) ;
 T56 = transformation_func(q(5) , 0 , 0 , 90) ;
 T67 = transformation_func(q(6) , 5+5 , 0 , 0) ;

       x = T01*T12*T23*T34*T45*T56*T67 



 
  
  end
  function J_inv = inverse_jacobian_matrix(q1,q2,q3,q4,q5,q6)
         x = simplify( 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180))
    y = simplify(5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))
    z = simplify( 5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5)                                                                                                                                                                   

         J_inv = [ diff(x,q1)  diff(x,q2)  diff(x,q3) diff(x,q4) diff(x,q5) diff(x,q6) ;
                   diff(y,q1)  diff(y,q2)  diff(y,q3) diff(y,q4) diff(y,q5) diff(y,q6);
                   diff(z,q1)  diff(z,q2)  diff(z,q3) diff(z,q4) diff(z,q5) diff(z,q6)]
                   
          
  end
  function q = inverse_position_kinematics(max_iterations, q_0, X)
  
      q_0=[q_0(1);q_0(2); q_0(3); q_0(4);q_0(5);q_0(6)];
      X = [X(1) ; X(2); X(3)]
      q_i=[q_0(1);q_0(2); q_0(3); q_0(4);q_0(5);q_0(6)];
      q1 = q_i(1);
      q2 = q_i(2);
      q3 = q_i(3);
      q4 = q_i(4);
      q5 = q_i(5);
      q6 = q_i(6);
      
      
       H=[10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*sin((pi*q4)/180))/180 - (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180))/180 + (pi*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/180) - (pi*sin((pi*q1)/180))/36 - (pi*cos((pi*q2)/180)*sin((pi*q1)/180))/36 - (pi*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))/18 - (pi*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180))/18 - (pi*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180))/18 - (pi*cos((pi*q1)/180)*sin((pi*q2)/180))/36 - 10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q1)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180))/18 - (pi*cos((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180))/18 - 10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q1)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180))/18 - (pi*cos((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18,  10*sin((pi*q5)/180)*((pi*cos((pi*q4)/180)*sin((pi*q1)/180))/180 - (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q4)/180))/180 + (pi*cos((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)*sin((pi*q4)/180))/180),   (pi*cos((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)))/18 - (pi*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*sin((pi*q5)/180))/18, 0
(pi*cos((pi*q1)/180))/36 + 10*sin((pi*q5)/180)*((pi*sin((pi*q1)/180)*sin((pi*q4)/180))/180 + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180))/180 - (pi*cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/180) + (pi*cos((pi*q1)/180)*cos((pi*q2)/180))/36 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180))/18 + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180))/18 + (pi*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))/18 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180))/36 - 10*sin((pi*q5)/180)*((pi*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q1)/180))/18 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))/18 - 10*sin((pi*q5)/180)*((pi*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q1)/180))/18 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18, -10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*cos((pi*q4)/180))/180 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)*sin((pi*q4)/180))/180 + (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q4)/180))/180), - (pi*cos((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)))/18 - (pi*sin((pi*(q2 + q3))/180)*sin((pi*q1)/180)*sin((pi*q5)/180))/18, 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0,                                                                                                                                                                                           (pi*sin((pi*(q2 + q3))/180))/18 + (pi*cos((pi*q2)/180))/36 + (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180))/36 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180))/18 - (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180))/36,                                                                                                                                                                          (pi*sin((pi*(q2 + q3))/180))/18 + (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180))/36 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180))/18 - (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180))/36,                                                                                                               (pi*cos((pi*(q4 + q5))/180)*sin((pi*(q2 + q3))/180))/36 - (pi*sin((pi*(q2 + q3))/180)*cos((pi*(q4 - q5))/180))/36,                                                                                                           (pi*cos((pi*(q4 + q5))/180)*sin((pi*(q2 + q3))/180))/36 + (pi*cos((pi*(q2 + q3))/180)*sin((pi*q5)/180))/18 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*(q4 - q5))/180))/36, 0];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    
     
%      
%     L = [ 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180);
%                   5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180);
%                   5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5 ]                                                                                                                                                                  
%   
%               
    q2= pinv(H)*X;
     q_i = q_i - q2 ;
      
      
   for   i = 0 : max_iterations 
     
  
   
  
     q1 = q_i(1);
      q2 = q_i(2);
      q3 = q_i(3);
      q4 = q_i(4);
      q5 = q_i(5);
      q6 = q_i(6);
      
      
       H=[10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*sin((pi*q4)/180))/180 - (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180))/180 + (pi*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/180) - (pi*sin((pi*q1)/180))/36 - (pi*cos((pi*q2)/180)*sin((pi*q1)/180))/36 - (pi*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))/18 - (pi*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180))/18 - (pi*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180))/18 - (pi*cos((pi*q1)/180)*sin((pi*q2)/180))/36 - 10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q1)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180))/18 - (pi*cos((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180))/18 - 10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q1)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180))/18 - (pi*cos((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18,  10*sin((pi*q5)/180)*((pi*cos((pi*q4)/180)*sin((pi*q1)/180))/180 - (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q4)/180))/180 + (pi*cos((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)*sin((pi*q4)/180))/180),   (pi*cos((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)))/18 - (pi*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*sin((pi*q5)/180))/18, 0
(pi*cos((pi*q1)/180))/36 + 10*sin((pi*q5)/180)*((pi*sin((pi*q1)/180)*sin((pi*q4)/180))/180 + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180))/180 - (pi*cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/180) + (pi*cos((pi*q1)/180)*cos((pi*q2)/180))/36 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180))/18 + (pi*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180))/18 + (pi*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))/18 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180))/36 - 10*sin((pi*q5)/180)*((pi*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q1)/180))/18 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18, (pi*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))/18 - 10*sin((pi*q5)/180)*((pi*cos((pi*q2)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q3)/180))/180 + (pi*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180))/180) + (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q1)/180))/18 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180))/18, -10*sin((pi*q5)/180)*((pi*cos((pi*q1)/180)*cos((pi*q4)/180))/180 - (pi*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)*sin((pi*q4)/180))/180 + (pi*cos((pi*q2)/180)*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q4)/180))/180), - (pi*cos((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)))/18 - (pi*sin((pi*(q2 + q3))/180)*sin((pi*q1)/180)*sin((pi*q5)/180))/18, 0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0,                                                                                                                                                                                           (pi*sin((pi*(q2 + q3))/180))/18 + (pi*cos((pi*q2)/180))/36 + (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180))/36 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180))/18 - (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180))/36,                                                                                                                                                                          (pi*sin((pi*(q2 + q3))/180))/18 + (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180))/36 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180))/18 - (pi*cos((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180))/36,                                                                                                               (pi*cos((pi*(q4 + q5))/180)*sin((pi*(q2 + q3))/180))/36 - (pi*sin((pi*(q2 + q3))/180)*cos((pi*(q4 - q5))/180))/36,                                                                                                           (pi*cos((pi*(q4 + q5))/180)*sin((pi*(q2 + q3))/180))/36 + (pi*cos((pi*(q2 + q3))/180)*sin((pi*q5)/180))/18 + (pi*sin((pi*(q2 + q3))/180)*cos((pi*(q4 - q5))/180))/36, 0];
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    
     
     
%     L = [ 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180);
%                   5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180);
%                   5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5 ] ;                                                                                                                                                                 

              
    q2= pinv(H)*X ;
     q_i = q_i - q2 ;
       q = q_i;
               
   end
     
           q = q_i

  
   
  end
  function J = jacobian_matrix(q1,q2,q3,q4,q5,q6) 
      
    T01 = transformation_func(0 , 5 , 0 , 0) ;
 T12 = transformation_func(q1 , 0 , 5 , 90) ;
 T23 = transformation_func(q2 , 0 , 5 , 0) ;
 T34 = transformation_func(q3 , 0 , 0 , 90) ;
 T45 = transformation_func(q4 , 5+5 , 0 , -90) ;
 T56 = transformation_func(q5 , 0 , 0 , 90) ;
 T67 = transformation_func(q6 , 5+5 , 0 , 0) 
 T02= T01*T12
 T03 = T01*T12*T23
 T04 = T01*T12*T23*T34
 T05 = T01*T12*T23*T34*T45
 T06 = T01*T12*T23*T34*T45*T56
 T07 = T01*T12*T23*T34*T45*T56*T67

   z1 = T01(1:3,3);
   z2 = T02(1:3,3);
   z3 = T03(1:3,3);
   z4 = T04(1:3,3);
   z5 = T05(1:3,3);
   z6 = T06(1:3,3);
   z7 = T07(1:3,3);

   
   O1= T01(1:3,4);
   O2= T02(1:3,4);
   O3= T03(1:3,4);
   O4= T04(1:3,4);
   O5= T05(1:3,4);
   O6= T06(1:3,4);
   O7= T07(1:3,4);
   
   jw2=z1;
   jw3=z2;
   jw4=z3;
   jw5=z4;
   jw6=z5;
   jw7=z6;
   on1 = (O7-O1);
    on2 = (O7-O2);
     on3 = (O7-O3);
      on4 = (O7-O4);
       on5 = (O7-O5);
        on6 = (O7-O6);
   
   jv2= cross(z1,on1);
   jv3=cross(z2,on2);
   jv4=cross(z3,on3) ;
   jv5=cross(z4,on4);
   jv6=cross(z5,on5);
   jv7=cross(z6,on6);
   
   J = [ jv2 jv3 jv4 jv5 jv6 jv7 ; jw2 jw3 jw4 jw5 jw6 jw7 ]
   
  
  
     
   
  end
  function V_F = forward_velocity_kinematics(q1,q2,q3,q4,q5,q6, q_dot)
    JM = jacobian_matrix(q1,q2,q3,q4,q5,q6) ;
    q_dot = [ q_dot(1) ;q_dot(2); q_dot(3); q_dot(4) ;q_dot(5); q_dot(6)];
    V_F = JM * q_dot
  end
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
  function Task_Space01 = task_traj_01(X0, Xf, tf)
% helical
x0= X0(1);
y0= X0(2);
z0= X0(3);
xf= Xf(1);
yf= Xf(2);
zf= Xf(3);

r= x0/cosd(0)
w= 4*180 
L0=z0
 c= (zf-L0)/tf
 t= 0 
 
 xt=r*cosd(w*t);
    yt=r*sind(w*t);
    zt=L0+(c*t);
    X= [ xt ; yt  ; zt ] 
    q =inverse_position_kinematics(3, [160,120,60,180,90,75], X0);
    Q =inverse_position_kinematics(3, q , X);
%  q =inverse_position_kinematics(3, [160,120,60,180,90,75], X0);
for i = 0.0:0.1:tf
    
    t = i 
    xt=r*cosd(w*t);
    yt=r*sind(w*t);
    zt=L0+(c*t);
    X= [ xt ; yt  ; zt ] 
    e =inverse_position_kinematics(3, Q , X)
    
    
    
end

 Task_Space01= [ xt ; yt  ; zt ] 
 
 
 

 end
  function Task_Space02 = task_traj_02(X0, Xf, tf)
% cone
x0= X0(1);
y0= X0(2);
z0= X0(3);
xf= Xf(1);
yf= Xf(2);
zf= Xf(3);

r= x0/cosd(0)
k= r-xf
w= 4*180 
L0=z0
c= (zf-L0)/tf
t=0
xt=(r-(k*t))*cosd(w*t);
yt=(r-(k*t))*sind(w*t);
 zt=L0+(c*t);
X= [ xt ; yt  ; zt ] 
 q =inverse_position_kinematics(3, [160,120,60,180,90,75], X0);
    Q =inverse_position_kinematics(3, q , X);

for i = 0.0:0.1:tf
    
    t = i 
    xt=(r-(k*t))*cosd(w*t);
    yt=(r-(k*t))*sind(w*t);
    zt=L0+(c*t);
    X= [ xt ; yt  ; zt ] 
    e =inverse_position_kinematics(3, Q , X)
    
    
end

 Task_Space02= [ xt ; yt  ; zt ] 

  end
  function Task_Space03 = task_traj_03(X0, Xf, tf)
% linear
x0= X0(1);
y0= X0(2);
z0= X0(3);
xf= Xf(1);
yf= Xf(2);
zf= Xf(3);

alphax= (xf-x0)/tf
alphay= (yf-y0)/tf
alphaz= (zf-z0)/tf
t=0
   xt= x0+(alphax*t);
    yt= y0+(alphay*t);
    zt= z0+(alphaz*t);
X= [ xt ; yt  ; zt ] 
 q =inverse_position_kinematics(3, [160,120,60,180,90,75], X0);
    Q =inverse_position_kinematics(3, q , X);

for i = 0.0:0.1:tf
    
    t = i 
    xt= x0+(alphax*t);
    yt= y0+(alphay*t);
    zt= z0+(alphaz*t);
      X= [ xt ; yt  ; zt ] 
    e =inverse_position_kinematics(3, Q , X)
    
    
end
Task_Space03 = [xt;yt;zt]
  end
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
  function L = Lagrange() 
  syms q1 q2 q3 q4 q5 q6 q1_dot q2_dot q3_dot q4_dot q5_dot q6_dot  I1 I2 I3 I4 I5 I6  m_1 m_2 m_3 m_4 m_5 m_6
  T =potential_energy([q1;q2;q3;q4;q5;q6], [m_1;m_2;m_3;m_4;m_5;m_6]);
  U = kinemtic_energy([q1;q2;q3;q4;q5;q6], [q1_dot ;q2_dot ;q3_dot ;q4_dot ;q5_dot ;q6_dot],[m_1;m_2;m_3;m_4;m_5;m_6],[I1;I2;I3;I4;I5;I6]);
  L = T - U
  
  end
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
  