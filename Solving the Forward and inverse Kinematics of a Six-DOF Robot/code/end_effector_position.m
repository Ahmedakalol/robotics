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
 