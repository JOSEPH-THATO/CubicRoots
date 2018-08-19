function x = Cubic(a,b,c,d)
    if a == 0 && b == 0
       x = -d/c;
    elseif a == 0
       delta = c.^2 - 4.00*b*d;
       delta = sqrt(delta);
       x1 = (-c + delta)/(2*b);
       x2 = (-c - delta)/(2*b);
       x = [x1,x2];
    else
       F = ((3.0 * c / a) - (b.^2/a.^2)) / 3.0;
       G = (((2.0 * b.^3.0) / a.^3.0) - ((9.0 * b * c) / a.^2.0) + (27.0 * d / a)) /27.0;
       H = G.^2.0 / 4.0 + F.^3.0 / 27.0; 
       if F == 0 && G == 0 && H == 0
           if (d / a) >= 0
               sol = ((d / a).^(1/3.0)) * (-1);
           else 
               sol = (-d / a).^(1/3.0);
           end
           x = [sol,sol,sol];
       elseif( H <= 0) 
            i = sqrt((G.^2.0 / 4.0) - H);
            j = i.^(1 / 3.0);
            k = acos(-(G / (2 * i)));
            L = j * -1;
            M = cos(k / 3.0);
            N = sqrt(3) * sin(k / 3.0);
            P = (b / (3.0 * a)) * (-1) ;
            x1 = 2 * j * cos(k / 3.0) - (b / (3.0 * a));
            x2 = L * (M + N) + P;
            x3 = L * (M - N) + P;
            x = [x1,x2,x3];
         
       else
            U = 0;
            S = 0;
            R = -(G / 2.0) + sqrt(H);
            if R >= 0
                S = R.^(1 / 3.0);
            else
                S = (-R).^(1 / 3.0) * -1;
                T = -(G / 2.0) - sqrt(H)
                if T >= 0
                    U = (T.^(1 / 3.0));
                else
                    U = ((-T).^(1 / 3.0)) * -1;
                end
            end     
            x1 = (S + U) - (b / (3.0 * a))
            x2 = -(S + U) / 2 - (b / (3.0 * a)) + (S - U) * sqrt(3) * 0.5j
            x3 = -(S + U) / 2 - (b / (3.0 * a)) - (S - U) * sqrt(3) * 0.5j
            x = [x1,x2,x3];
       end
           
    end
end