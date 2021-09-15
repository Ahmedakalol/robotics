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




