% Defining the function
clear;
syms x y
f(x,y) = (1/3)*x^2 + 3*y^2; 
fhandle = matlabFunction(f);

gradf = gradient(f , [x, y]); 


%Plot function f(x,y)

[X,Y] = meshgrid(-6:6);

figure;
surf(X,Y,fhandle(X,Y)); % first type of plot
shading interp;  % Smooth the surface
title('Function f(x,y)', 'FontSize', 16);
grid on;

%First case

[x1, y1, f1_values, k1] = steepest_descent(f, gradf, 2, 2, 0.0001, 0.1);
disp(k1)
disp(f1_values(k1));
figure
for i=1:1:k1
    plot(i, f1_values(i),'ob','MarkerSize', 5)
    hold on
end

ylim([0, 14]);
xlim([0, 200]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g = 0.1 with initial (x0,y0) = (2,2) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%Second case

[x2, y2, f2_values, k2] = steepest_descent(f, gradf, 2, 2, 0.0001, 0.3);
disp(k2)
disp(f2_values(k2));
figure
for i=1:1:k2
    plot(i, f2_values(i),'ob','MarkerSize', 5)
    hold on
end

ylim([0, 14]);
xlim([0, 60]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g = 0.3 with initial (x0,y0) = (2,2) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%Thord case

[x3, y3, f3_values, k3] = steepest_descent(f, gradf, 2, 2, 0.0001, 3);
disp(k3)
disp(f3_values(k3));

figure
for i=1:1:k3
    plot(i, f3_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0, 200]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g = 3 with initial (x0,y0) = (2,2) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%Plot in logarithmic scale

figure
for i=1:1:k3
    semilogy(i, f3_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0, 200]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('Logarithmic scale: f(x, y)', 'FontSize', 14);
title('Constant g = 3 with initial (x0,y0) = (2,2) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%Fourth case

[x4, y4, f4_values, k4] = steepest_descent(f, gradf, 2, 2, 0.0001, 5);
disp(k4)
disp(f4_values(k4));
figure
for i=1:1:k4
    plot(i, f4_values(i),'ob','MarkerSize', 5)
    hold on
end

xlim([0, 200]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g  = 5 with initial (x0,y0) = (2,2) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%Plot in logarithmic scale

figure
for i=1:1:k4
    semilogy(i, f4_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0, 200]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('Logarithmic scale: f(x, y)', 'FontSize', 14);
title('Constant g = 5 with initial (x0,y0) = (2,2) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 14);
grid on;

%--------------------------------------------------Steepest Descent---------------------------------------------------------------------------%
function [x, y, f_values, k] = steepest_descent(f, gradf, initial_x, initial_y, e, g)

   x = []; % stores all x values
   y = []; % stores all y values
   f_values = []; % stores all g(k) values
   dk_values = zeros(2, 20);
   k = 1;  %Number of iterations

   x(1) = initial_x;
   y(1) = initial_y;
   
   disp(['inx and iny', num2str(x(1)), num2str(y(1))]);
   dk_values(:,1) = -gradf(x(1), y(1));

   f_values(1) =  f(x(1), y(1));
   disp(['f(2,2):', num2str(f_values(1))]);
    
   while (norm(gradf(x(k),y(k))) >= e)
          
       x(k + 1) = (x(k) + g * dk_values(1, k));
       y(k + 1) = (y(k) + g * dk_values(2, k));
       f_values(k + 1) = f(x(k + 1), y(k + 1));
       dk_values(:, k + 1) = -gradf(x(k+1),y(k+1));    

       if(k > 200)
        break
       end 

       k = k + 1;
   end
end

