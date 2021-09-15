syms q1 q2 q3 q4 q5 q6 L0 L1 L2 L3 L4 L5 L6 T alpha theta
%  jacobian_matrix(q1,q2,q3,q4,q5,q6) 
%  
%  L = [ 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180);
%                   5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180);
%                   5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5 ]                                                                                                                                                                  


% 
%     inverse_jacobian_matrix(q1,q2,q3,q4,q5,q6)

    
    forward_velocity_kinematics(279.8473,106.3596 ,76.5489, 304.6278,76.5047,0, [2;2;2;2;2;2])

% %    inverse_position_kinematics(7, [40;90;45;150;130;40],L)
% task_traj_01([2;0;1], [2;0;2], 1.0)
% task_traj_02([2;0;1], [2;0;2], 1.0)
% task_traj_03([2;2;1], [-2;2.5;3], 1.0)
% q= joint_traj([160,120,60,180,90,75], [160,120,60,180,90,75], [0,0,0,0,0,0], [0,0,0,0,0,0], 10)
% % for i = 1: 10 
%   end_effector_position(q(1:6,i))
% end
% X3 = end_effector_position3([q1,q2,q3])
% J_inv3 = inverse_jacobian_matrix3(q1,q2,q3)
Task_Space3 = task_traj_3([2;2;1], [-2;2.5;3], 1.0)


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
