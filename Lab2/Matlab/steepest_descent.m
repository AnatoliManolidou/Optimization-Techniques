% Defining the function
clear;
syms x y
f(x,y) = x^5 * exp(-(x^2 + y^2)); 
fhandle = matlabFunction(f);

gradf = gradient(f , [x, y]); 


%Plot function f(x,y)

[X,Y] = meshgrid(-5:0.1:5);
Z = (X.^5).*(exp(-(X.^2) - (Y.^2)));

surf(X,Y,Z); % first type of plot
shading interp;  % Smooth the surface
title('Function f(x,y)', 'FontSize', 16);
grid on;
figure; 

contourf(X,Y,Z); % second type of plot
title('Contour lines of f(x,y)', 'FontSize', 16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Constant g%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x1, y1, f1_values, k1] = steepest_descent(f, gradf,'constant', 0, 0, 0.0001);

figure
for i=1:1:k1
    plot(i, f1_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k1+2])
ylim([0, 0.01]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g with initial (x0,y0) = (0,0) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;


[x2, y2, f2_values, k2] = steepest_descent(f, gradf,'constant', -1, 1, 0.0001);
figure
for i=1:1:k2
    plot(i, f2_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k2+2])
ylim([-1, -0.15]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g with initial (x0,y0) = (-1,1) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;

[x3, y3, f3_values, k3] = steepest_descent(f, gradf,'constant', 1, -1, 0.0001);
figure
for i=1:1:k3
    plot(i, f3_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k3+2])
ylim([0, 0.15]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Constant g with initial(x0,y0) = (1,-1) ','Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Armijo%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x4, y4, f4_values, k4] = steepest_descent(f, gradf,'armijo', 0, 0, 0.0001);

figure
for i=1:1:k4
    plot(i, f4_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k4+2])
ylim([0, 0.01]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Armijo with initial (x0,y0) = (0,0) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;


[x5, y5, f5_values, k5] = steepest_descent(f, gradf,'armijo', -1, 1, 0.0001);
figure
for i=1:1:k5
    plot(i, f5_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k5+2])
ylim([-1, -0.15]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Armijo with initial (x0,y0) = (-1,1) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;

[x6, y6, f6_values, k6] = steepest_descent(f, gradf,'armijo', 1, -1, 0.0001);
figure
for i=1:1:k6
    plot(i, f6_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k6+2])
ylim([0, 0.15]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Armijo with initial (x0,y0) = (1,-1) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Minimization of g using Bisection with Derivatives%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x7, y7, f7_values, k7] = steepest_descent(f, gradf,'minimization', 0, 0, 0.0001);

figure
for i=1:1:k7
    plot(i, f7_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k7+2])
ylim([0, 0.01]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Minimization of g with initial (x0,y0) = (0,0) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;


[x8, y8, f8_values, k8] = steepest_descent(f, gradf,'minimization', -1, 1, 0.0001);
figure
for i=1:1:k8
    plot(i, f5_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k8+2])
ylim([-1, -0.15]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Minimization of g with initial (x0,y0) = (-1,1) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;

[x9, y9, f9_values, k9] = steepest_descent(f, gradf,'minimization', 1, -1, 0.0001);
figure
for i=1:1:k9
    plot(i, f9_values(i),'ob','MarkerSize', 5)
    hold on
end
xlim([0 k9+2])
ylim([0, 0.15]);
xlabel('k Iterations', 'FontSize', 14);
ylabel('f(x, y)', 'FontSize', 14);
title('Minimization of g with initial (x0,y0) = (1,-1) ',' Convergence of the Objective Function over k Iterations', 'FontSize', 16);
grid on;

%--------------------------------------------------Steepest Descent---------------------------------------------------------------------------%
function [x, y, f_values, k] = steepest_descent(f, gradf, g_choice, initial_x, initial_y, e)

   syms x y g

   x = []; % stores all x values
   y = []; % stores all y values
   g = []; % stores all g(k) values
   f_values = []; % stores all g(k) values
   dk_values = zeros(2, 20);
   m = []; % stores all m(k) values for armijo
   k = 1;  %Number of iterations

   x(1) = initial_x;
   y(1) = initial_y;
   
   disp(['inx and iny', num2str(x(1)), num2str(y(1))]);
   dk_values(:,1) = -gradf(x(1), y(1));

   f_values(1) =  f(x(1), y(1));
   disp(['f(1,1):', num2str(f_values(1))]);
    
   while (norm(gradf(x(k),y(k))) >= e)
      
       if (strcmp(g_choice,'constant'))
           g(k) = 0.6;   
       elseif (strcmp(g_choice,'armijo'))
           
           a1 = 0.001;
           b1 = 0.3;
           s = 0.5;

           temp_mk = 1;
           temp_g = s*(b1^temp_mk);
         
           while (f(x(k),y(k)) - f(x(k) + g * dk_values(1, k) , y(k) + g * dk_values(2, k)) < (-a1 * (b1^temp_mk) * s * (dk_values(: , k)') * gradf(x(k),y(k))))
               temp_mk = temp_mk + 1;
               temp_g = s*(b1^temp_mk);
           end

           m(k) = temp_mk;
           g(k) = temp_g;
       elseif (strcmp(g_choice,'minimization'))
           h(g) = f(x(k) + g * dk_values(1,k) , y(k) + g * dk_values(2,k));
           [a, b, c] = bisection_derivatives(h, 0, 5, 0.005);
   
           g(k) = (a(c) + b(c))/ 2;
       end
       x(k + 1) = (x(k) + g(k) * dk_values(1, k));
       y(k + 1) = (y(k) + g(k) * dk_values(2, k));
       f_values(k + 1) = f(x(k + 1), y(k + 1));
       dk_values(:, k + 1) = -gradf(x(k+1),y(k+1));    

       if(k > 200)
        break
       end 

       k = k + 1;
   end
end

%-----------------------------------------------BISECTION WITH DERIVATIVES FUNCTION------------------------------------------------------------------------------%
function [a, b, k] = bisection_derivatives(h, in_a, in_b, l)
    % Initializing the variables

    syms g;
    a = [];
    b = [];
    x = [];
    k = 1;

    a(1) = in_a; 
    b(1) = in_b;                   

    dh = diff(h,g);
    n = ceil(log(l/(b(k)-a(k)))/log(1/2))

    for k=1:n

        x(k) = (a(k) + b(k))/2;
        dhxk = subs(dh, x(k));

        if(dhxk > 0)
            b(k) = x(k);
        elseif(dhxk < 0)
            a(k) = x(k);
        else
            break;
        end
    end
end