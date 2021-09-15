 function J_inv = inverse_jacobian_matrix(q1,q2,q3,q4,q5,q6)
         x = simplify( 5*cos((pi*q1)/180) + 5*cos((pi*q1)/180)*cos((pi*q2)/180) + 10*sin((pi*q5)/180)*(sin((pi*q1)/180)*sin((pi*q4)/180) + cos((pi*q1)/180)*cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180) - cos((pi*q1)/180)*cos((pi*q4)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q1)/180)*cos((pi*q5)/180) + 10*cos((pi*q1)/180)*cos((pi*q2)/180)*sin((pi*q3)/180) + 10*cos((pi*q1)/180)*cos((pi*q3)/180)*sin((pi*q2)/180))
    y = simplify(5*sin((pi*q1)/180) + 5*cos((pi*q2)/180)*sin((pi*q1)/180) - 10*sin((pi*q5)/180)*(cos((pi*q1)/180)*sin((pi*q4)/180) - cos((pi*q2)/180)*cos((pi*q3)/180)*cos((pi*q4)/180)*sin((pi*q1)/180) + cos((pi*q4)/180)*sin((pi*q1)/180)*sin((pi*q2)/180)*sin((pi*q3)/180)) + 10*cos((pi*q2)/180)*sin((pi*q1)/180)*sin((pi*q3)/180) + 10*cos((pi*q3)/180)*sin((pi*q1)/180)*sin((pi*q2)/180) + 10*sin((pi*(q2 + q3))/180)*cos((pi*q5)/180)*sin((pi*q1)/180))
    z = simplify( 5*sin((pi*q2)/180) - 10*cos((pi*(q2 + q3))/180) + 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 + q5))/180) - 10*cos((pi*(q2 + q3))/180)*cos((pi*q5)/180) - 5*sin((pi*(q2 + q3))/180)*sin((pi*(q4 - q5))/180) + 5)                                                                                                                                                                   

         J_inv = [ diff(x,q1)  diff(x,q2)  diff(x,q3) diff(x,q4) diff(x,q5) diff(x,q6) ;
                   diff(y,q1)  diff(y,q2)  diff(y,q3) diff(y,q4) diff(y,q5) diff(y,q6);
                   diff(z,q1)  diff(z,q2)  diff(z,q3) diff(z,q4) diff(z,q5) diff(z,q6)]
                   
          
  end