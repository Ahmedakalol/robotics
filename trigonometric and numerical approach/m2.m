syms q1 q2 q3 q4 q5 q6 L0 L1 L2 L3 L4 L5 L6 T alpha theta q_double_dot1 q_double_dot2 q_double_dot3 q_double_dot4 q_double_dot5 q_double_dot6 q_dot1 q_dot2 q_dot3 q_dot4 q_dot5 q_dot6 
  
  L = [ 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180);
                  5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180);
                  5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5 ]                                                                                                                                                                  
  
   J_inv = inverse_jacobian_matrix(q1,q2,q3,q4,q5,q6)
   q = inverse_position_kinematics(8, [160,120,60,180,90,75], L)
   J= jacobian_matrix(q1,q2,q3,q4,q5,q6) 
   V_F = forward_velocity_kinematics(q1,q2,q3,q4,q5,q6, [q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6])
   q_dot = inverse_velocity_kinematics(q1,q2,q3,q4,q5,q6, V_F)
   j_dot = jacobian_derivative(q1,q2,q3,q4,q5,q6,[q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6])
   A_F =  forward_acceleration_kinematics(q1,q2,q3,q4,q5,q6, [q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6], [q_double_dot1;q_double_dot2;q_double_dot3;q_double_dot4;q_double_dot5;q_double_dot6])
   q_double_dot = simplify( inverse_acceleration_kinematics(q1,q2,q3,q4,q5,q6, [q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6], A_F))
 

  function T = transformation_func(theta , d , a , alpha) 
           T = [cosd(theta)  -sind(theta)*cosd(alpha) sind(theta)*sind(alpha) a*cosd(theta) ;
       sind(theta) cosd(theta)*cosd(alpha) -cosd(theta)*sind(alpha) a*sind(theta);
       0 sind(alpha) cosd(alpha) d;
       0 0 0 1]
  % X= T01 .* T12 .* T23 .* T34 .* T45 .* T56 .* T67


  end
  function X = end_effector_position(q1,q2,q3,q4,q5,q6)
  T01 = transformation_func(0 , 5 , 0 , 0) ;
 T12 = transformation_func(q1 , 0 , 5 , 90) ;
 T23 = transformation_func(q2 , 0 , 5 , 0) ;
 T34 = transformation_func(q3 , 0 , 0 , 90) ;
 T45 = transformation_func(q4 , 5+5 , 0 , -90) ;
 T56 = transformation_func(q5 , 0 , 0 , 90) ;
 T67 = transformation_func(q6 , 5+5 , 0 , 0) ;

       x = simplify(T01*T12*T23*T34*T45*T56*T67 )



 
  
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    
     
     
    L = [ 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180);
                  5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180);
                  5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5 ]                                                                                                                                                                  
  
              
    q2= pinv(H)*L ;
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    
     
     
    L = [ 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180);
                  5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180);
                  5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5 ] ;                                                                                                                                                                 

              
    q2= pinv(H)*L ;
     q_i = q_i - q2 ;
       q = q_i
               
   end
     
    
  
   
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
 T03 = simplify(T01*T12*T23)
 T04 = simplify(T01*T12*T23*T34)
 T05 = simplify (T01*T12*T23*T34*T45)
 T06 = simplify(T01*T12*T23*T34*T45*T56)
 T07 = simplify(T01*T12*T23*T34*T45*T56*T67)

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
  function j_dot =jacobian_derivative(q1,q2,q3,q4,q5,q6,q_dot)
  J = jacobian_matrix(q1,q2,q3,q4,q5,q6)
   q_dot1=q_dot(1);
     q_dot2=q_dot(2);
     q_dot3=q_dot(3);
     q_dot4=q_dot(4);
     q_dot5=q_dot(5);
     q_dot6=q_dot(6); 
    jv11 = diff(J(1,1),q1)*q_dot1
    jv21 = diff(J(2,1),q1)*q_dot1
    jv31 = diff(J(3,1),q1)*q_dot1
    jv41 = diff(J(4,1),q1)*q_dot1
    jv51 = diff(J(5,1),q1)*q_dot1
    jv61 = diff(J(6,1),q1)*q_dot1
    jv12 = diff(J(1,2),q2)*q_dot2
    jv22 = diff(J(2,2),q2)*q_dot2
    jv32 = diff(J(3,2),q2)*q_dot2
    jv42 = diff(J(4,2),q2)*q_dot2
    jv52 = diff(J(5,2),q2)*q_dot2
    jv62 = diff(J(6,2),q2)*q_dot2
    jv13 = diff(J(1,3),q3)*q_dot3
    jv23 = diff(J(2,3),q3)*q_dot3;
    jv33 = diff(J(3,3),q3)*q_dot3
    jv43 = diff(J(4,3),q3)*q_dot3
    jv53 = diff(J(5,3),q3)*q_dot3
    jv63 = diff(J(6,3),q3)*q_dot3
    jv14 = diff(J(1,4),q4)*q_dot4
    jv24 = diff(J(2,4),q4)*q_dot4
    jv34 = diff(J(3,4),q4)*q_dot4
    jv44 = diff(J(4,4),q4)*q_dot4
    jv54 = diff(J(5,4),q4)*q_dot4
    jv64 = diff(J(6,4),q4)*q_dot4
    jv15 = diff(J(1,5),q5)*q_dot5
    jv25 = diff(J(2,5),q5)*q_dot5
    jv35 = diff(J(3,5),q5)*q_dot5
    jv45 = diff(J(4,5),q5)*q_dot5
    jv55 = diff(J(5,5),q5)*q_dot5
    jv65 = diff(J(6,5),q5)*q_dot5
    jv16 = diff(J(1,6),q6)*q_dot6
    jv26 = diff(J(2,6),q6)*q_dot6
    jv36 = diff(J(3,6),q6)*q_dot6
    jv46 = diff(J(4,6),q6)*q_dot6
    jv56 = diff(J(5,6),q6)*q_dot6
    jv66 = diff(J(6,6),q6)*q_dot6
    
    
    j_dot = [ jv11 jv12  jv13 jv14 jv15 jv16 ;
        jv21 jv22 jv23 jv24 jv25 jv26 ; 
        jv31 jv32 jv33 jv34 jv35 jv36 ;
        jv41 jv42 jv43 jv44 jv45 jv46 ;  
        jv51 jv52 jv53 jv54 jv55 jv56 ; 
        jv61 jv62 jv63 jv64 jv65 jv66 ]
  end
  function A_F = forward_acceleration_kinematics(q1,q2,q3,q4,q5,q6, q_dot, q_double_dot)
   
   j_dot =jacobian_derivative(q1,q2,q3,q4,q5,q6,q_dot)                                              
  J = jacobian_matrix(q1,q2,q3,q4,q5,q6)
    A = j_dot * q_dot 
    B = J * q_double_dot 
    A_F = A + B
  end
  function q_dot = inverse_velocity_kinematics(q1,q2,q3,q4,q5,q6, V_F)
     
    
      J = jacobian_matrix(q1,q2,q3,q4,q5,q6)
  
      q_dot = inv(J) * V_F 
  
  end
  function q_double_dot = inverse_acceleration_kinematics(q1,q2,q3,q4,q5,q6, q_dot, A_F)
      
      J = jacobian_matrix(q1,q2,q3,q4,q5,q6)
      j_dot =jacobian_derivative(q1,q2,q3,q4,q5,q6,q_dot)
      A = j_dot * q_dot
      B = A_F - A
      q_double_dot = inv(J) * B 
  end
  