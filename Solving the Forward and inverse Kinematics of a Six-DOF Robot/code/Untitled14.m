syms q1 q2 q3 q4 q5 q6 L0 L1 L2 L3 L4 L5 L6 T alpha theta q_double_dot1 q_double_dot2 q_double_dot3 q_double_dot4 q_double_dot5 q_double_dot6 q_dot1 q_dot2 q_dot3 q_dot4 q_dot5 q_dot6 X
  
  L = [ 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180);
                  5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180);
                  5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5 ]                                                                                                                                                                  
   
%    q = inverse_position_kinematics(8, [160,120,60,180,90,75], L)
% %    J_inv = inverse_jacobian_matrix(q1,q2,q3,q4,q5,q6)
%    V_F = forward_velocity_kinematics(q1,q2,q3,q4,q5,q6, [q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6])
%    q_dot = inverse_velocity_kinematics(q1,q2,q3,q4,q5,q6, V_F)
%    j_dot = jacobian_derivative(q1,q2,q3,q4,q5,q6,[q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6])
%    A_F =  forward_acceleration_kinematics(q1,q2,q3,q4,q5,q6, [q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6], [q_double_dot1;q_double_dot2;q_double_dot3;q_double_dot4;q_double_dot5;q_double_dot6])
%    q_double_dot = simplify( inverse_acceleration_kinematics(q1,q2,q3,q4,q5,q6, [q_dot1;q_dot2;q_dot3;q_dot4;q_dot5;q_dot6], A_F))
% %  
task_traj_01([2;0;1], [2;0;2], 1.0)
task_traj_02([2;0;1], [2;0;2], 1.0)
task_traj_03([2;2;1], [-2;2.5;3], 1.0)
q= joint_traj([160,120,60,180,90,75], [160,120,60,180,90,75], [0,0,0,0,0,0], [0,0,0,0,0,0], 10)
for i = 1: 10 
  end_effector_position(q(1:6,i))
end
