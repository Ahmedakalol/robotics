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