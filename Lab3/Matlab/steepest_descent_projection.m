% Defining the function
clear;
syms x y
f(x,y) = (1/3)*x^2 + 3*y^2; 
fhandle = matlabFunction(f);

gradf = gradient(f , [x, y]); 

%Plot function f(x,y)

[X,Y] = meshgrid(-10:10);

figure;
surf(X,Y,fhandle(X,Y)); 
shading interp;  
title('Function f(x,y)', 'FontSize', 16);
grid on;

%First case

[x1, y1, f1_values, k1] = steepest_descent_projection(f, gradf, 5, -5, 0.01, 0.5, 5);
disp(k1);
disp(['Size of f1_values: ', num2str(length(f1_values))]);
figure;
plot(1:k1, f1_values, 'ob', 'MarkerSize', 5);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

ylim([0, 300]);
xlim([0, 60]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g with initial (x0,y0) = (5,-5) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%Second case

[x2, y2, f2_values, k2] = steepest_descent_projection(f, gradf, -5, 10, 0.01, 0.1, 15);
disp(k2);
disp(['Size of f2_values: ', num2str(length(f2_values))]);
figure;
plot(1:k2, f2_values, 'ob', 'MarkerSize', 5);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

xlim([0, 60]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g with initial (x0,y0) = (-5,10) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%Third case

[x3, y3, f3_values, k3] = steepest_descent_projection(f, gradf, 8, -10, 0.01, 0.2, 0.1);
disp(k3);
disp(['Size of f2_values: ', num2str(length(f3_values))]);
figure;
plot(1:k3, f3_values, 'ob', 'MarkerSize', 5);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

xlim([0, 60]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g with initial (x0,y0) = (8,-10) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%--------------------------------------------------Steepest Descent with Projection---------------------------------------------------------------------------%
function [x, y, f_values, k] = steepest_descent_projection(f, gradf, initial_x, initial_y, e, g, s)

   x = []; % stores all x values
   y = []; % stores all y values
   f_values = []; % stores all g(k) values
   dk_values = zeros(2, 20);
   xbar = [];
   ybar = [];
   xpbar = [];
   ypbar = [];

   k = 1;  %Number of iterations

   x(1) = initial_x;
   y(1) = initial_y;

   x(k) = projection(x(k), 'x');
   y(k) = projection(y(k), 'y');

   dk_values(:,1) = -gradf(x(1), y(1));

   f_values(1) =  f(x(1), y(1));
   disp(['f(5,-5):', num2str(f_values(1))]);
    
   while (norm(gradf(x(k),y(k))) > e)
          
       xbar(k) = x(k) + s * dk_values(1,k);
       ybar(k) = y(k) + s * dk_values(2,k);
       xpbar(k) = projection(xbar(k), 'x');
       ypbar(k) = projection(ybar(k), 'y');

       x(k + 1) = (x(k) + g * (xpbar(k) - x(k)));
       y(k + 1) = (y(k) + g * (ypbar(k) - y(k)));

       f_values(k + 1) = f(x(k + 1), y(k + 1));
       dk_values(:, k + 1) = -gradf(x(k+1),y(k+1));    

       if(k > 60)
            k = k + 1;
            break
       end 

       k = k + 1;
   end
end

%--------------------------------------------------Projection---------------------------------------------------------------------------%

function [varbar] = projection(var, var_choice)

    if (strcmp(var_choice,'x'))
        if(var < -10)
            a = -10;
        elseif (var > 5)
            a = 5;
        else
            a = var;
        end
    end
  
    if (strcmp(var_choice,'y'))
        if(var < -8)
            a = -8;
        elseif (var > 12)
            a = 12;
        else
            a = var;
        end
    end 

     varbar = a;

end