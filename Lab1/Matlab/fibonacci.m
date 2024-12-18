% Defining the functions 
f1 = @(x) (((x-2)^2) + (x * (log(x+3) / log(exp(1)))));
f2 = @(x) (exp(-2 * x) + (x-2)^2);
f3 = @(x) (exp(x) * (x^3 - 1) + (x-1) * sin(x));

% Defining the initial interval
in_a = -1;
in_b = 3;

l_values = [];
e = 0.001;
step = 0.01;
k = 1;

for(l = 0.003:step: 0.1) % Testing for different values of l
    [a1, b1, f1_count] = Fibonacci(f1, in_a, in_b, l, e);    
    [a2, b2, f2_count] = Fibonacci(f2, in_a, in_b, l, e);
    [a3, b3, f3_count] = Fibonacci(f3, in_a, in_b, l, e);
    if(k == 1)
        interv1a = a1;
        interv1b = b1;
        interv2a = a2;
        interv2b = b2;
        interv3a = a3;
        interv3b = b3;
    end
    f1count(k) = f1_count;
    f2count(k) = f2_count;
    f3count(k) = f3_count;
    l_values(k) = l;
    k = k + 1;
end  

% Plot results for f1(x)
figure;
subplot(3,1,1);
plot(l_values, f1count, '-o'); 
ylabel('Function Evaluations');
xlabel('l values');
title('Fibonacci Method for f1(x): evaluations of f1(x) against different values of l');
grid on;

% Plot results for f2(x)
subplot(3,1,2);
plot(l_values, f2count, '-o'); 
ylabel('Function Evaluations');
xlabel('l values');
title('Fibonacci Method for f2(x): evaluations of f2(x) against different values of l');
grid on;

% Plot results for f3(x)
subplot(3,1,3);
plot(l_values, f3count,'-o'); 
ylabel('Function Evaluations');
xlabel('l values');
title('Fibonacci Method for f3(x): evaluations of f3(x) against different values of l');
grid on;

% Plot ak vs k and bk vs k for different l values
figure;
plot(1:length(interv1a), interv1a, '-o', 'DisplayName', 'a_k(l=0.003)');
hold on;
plot(1:length(interv1b), interv1b, '-o', 'DisplayName', 'b_k(l=0.003)');
plot(1:length(a1), a1, '-o', 'DisplayName', 'a_k(l=0.1)');
plot(1:length(b1), b1, '-o', 'DisplayName', 'b_k(l=0.1)');

xlabel('Iteration k');
ylabel('Interval Boundaries');
legend;
title('Fibonacci Method for f1(x): ak and bk against iteration k');
grid on;

% Plot ak vs k and bk vs k for different l values
figure;
plot(1:length(interv2a), interv2a, '-o', 'DisplayName', 'a_k(l=0.003)');
hold on;
plot(1:length(interv2b), interv2b, '-o', 'DisplayName', 'b_k(l=0.003)');
plot(1:length(a2), a2, '-o', 'DisplayName', 'a_k(l=0.1)');
plot(1:length(b2), b2, '-o', 'DisplayName', 'b_k(l=0.1)');

xlabel('Iteration k');
ylabel('Interval Boundaries');
legend;
title('Fibonacci Method for f2(x): ak and bk against iteration k');
grid on;

% Plot ak vs k and bk vs k for different l values
figure;
plot(1:length(interv3a), interv3a, '-o', 'DisplayName', 'a_k(l=0.003)');
hold on;
plot(1:length(interv3b), interv3b, '-o', 'DisplayName', 'b_k(l=0.003)');
plot(1:length(a3), a3, '-o', 'DisplayName', 'a_k(l=0.1)');
plot(1:length(b3), b3, '-o', 'DisplayName', 'b_k(l=0.1)');

xlabel('Iteration k');
ylabel('Interval Boundaries');
legend;
title('Fibonacci Method for f3(x): ak and bk against iteration k');
grid on;
%---------------------------------------FIBONACCI----------------------------------------------------%
function [a, b, f_count] = Fibonacci(f, in_a, in_b, l, e)
    % Initializing interval boundaries
    a(1) = in_a; 
    b(1) = in_b;
    f_count = 0;                           
    k = 1;     

    % Initializing Fibonacci sequence with the first two values
    F = [1, 1];  
    n = 2;

    % Calculating Fibonacci sequence until it meets the required interval length
    while F(n) <= (b(1) - a(1)) / l
        F(n+1) = F(n) + F(n-1);
        n = n + 1;
    end

    % Seting initial x1 and x2 based on the Fibonacci ratios
    x1(k) = a(k) + (F(n-k-1) / F(n-k+1)) * (b(k) - a(k)); 
    x2(k) = a(k) + (F(n-k) / F(n-k+1)) * (b(k) - a(k)); 

    if (b(k) - a(k) < l)
        return
    end

    while (k < n - 2)
       
        if f(x1(k)) > f(x2(k))
            a(k+1) = x1(k);
            b(k+1) = b(k);
            x1(k + 1) = x2(k);  
            x2(k+1) = a(k+1) + (F(n-k-1) / F(n-k)) * (b(k+1) - a(k+1));
        else
            a(k+1) = a(k);
            b(k+1) = x2(k);
            x1(k + 1) = a(k+1) + (F(n-k-2) / F(n-k)) * (b(k+1) - a(k+1)); 
            x2(k+1) = x1(k + 1);
        end 

        k = k + 1;
        if(k == 1)
            f_count = f_count + 2;
        else
            f_count = f_count + 1;
        end

    end

    % Adjustments if k == n - 2
    if (k == n - 2)
        if f(x1(k-1)) > f(x2(k-1))
            x1(k) = x2(k-1);
            x2(k) = x1(k) + e;
        else
            x1(k) = x2(k) - e;
            x2(k) = x1(k-1);
        end

        if(k == 1)
            f_count = f_count + 2;
        else
            f_count = f_count + 1;
        end

        if (f(x1(k)) > f(x2(k)))
            a(k+1) = x1(k);
            b(k+1) = b(k);
        else
            a(k+1) = a(k);
            b(k+1) = x2(k);
        end
    end
end