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