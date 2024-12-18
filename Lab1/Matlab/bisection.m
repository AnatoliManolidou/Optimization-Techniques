% Defining the functions
f1 = @(x) ((x - 2)^2 + x * (log(x + 3) / log(exp(1))));
f2 = @(x) (exp(-2 * x) + (x - 2)^2);
f3 = @(x) (exp(x) * (x^3 - 1) + (x - 1) * sin(x));

% Defining the initial interval
in_a = -1;
in_b = 3;

e_values = [];          
l_values = [];

% First Case: Constant l and varying e
l1 = 0.01;
step = 0.0005;
k = 1;

for e1 = 0.001:step:(l1 / 2) - 0.0005 % Testing for different values of e
    [a11, b11, f11_count] = bisection_(f1, in_a, in_b, l1, e1);
    [a12, b12, f12_count] = bisection_(f2, in_a, in_b, l1, e1);
    [a13, b13, f13_count] = bisection_(f3, in_a, in_b, l1, e1);
    
    e_values(k) = e1;
    f11count(k) = f11_count;
    f12count(k) = f12_count;
    f13count(k) = f13_count;

    k = k + 1;
end    

% Plot results for f1(x)
figure;
subplot(3,1,1);
plot(e_values, f11count, '-o'); 
ylabel('Function Evaluations');
xlabel('e values');
title('Bisection Method for f1(x): Evaluations of f1(x) against different values of e');
grid on;

% Plot results for f2(x)
subplot(3,1,2);
plot(e_values, f12count, '-o'); 
ylabel('Function Evaluations');
xlabel('e values');
title('Bisection Method for f2(x): Evaluations of f2(x) against different values of e');
grid on;

% Plot results for f3(x)
subplot(3,1,3);
plot(e_values, f13count,'-o'); 
ylabel('Function Evaluations');
xlabel('e values');
title('Bisection Method for f3(x): Evaluations of f3(x) against different values of e');
grid on;

% Second Case: Constant e and varying l
e2 = 0.0001;
step = 0.01;
k = 1;

interv1a = [];
interv1b = [];
interv2a = [];
interv2b = [];
interv3a = [];
interv3b = [];

for l2 = 0.003:step:0.1  % Testing for different values of l
    [a21, b21, f21_count] = bisection_(f1, in_a, in_b, l2, e2);
    [a22, b22, f22_count] = bisection_(f2, in_a, in_b, l2, e2);
    [a23, b23, f23_count] = bisection_(f3, in_a, in_b, l2, e2);

    if(k == 1)
        interv1a = a21;
        interv1b = b21;
        interv2a = a22;
        interv2b = b22;
        interv3a = a23;
        interv3b = b23;
    end
    l_values(k) = l2;
    f21count(k) = f21_count;
    f22count(k) = f22_count;
    f23count(k) = f23_count;

    k = k + 1;
end  

% Plot results for f1(x)
figure;
subplot(3,1,1);
plot(l_values, f21count, '-o'); 
ylabel('Function Evaluations');
xlabel('l values');
title('Bisection Method for f1(x): evaluations of f1(x) against different values of l');
grid on;

% Plot results for f2(x)
subplot(3,1,2);
plot(l_values, f22count, '-o'); 
ylabel('Function Evaluations');
xlabel('l values');
title('Bisection Method for f2(x): evaluations of f2(x) against different values of l');
grid on;

% Plot results for f3(x)
subplot(3,1,3);
plot(l_values, f23count, '-o'); 
ylabel('Function Evaluations');
xlabel('l values');
title('Bisection Method for f3(x): evaluations of f3(x) against different values of l');
grid on;

% Plot ak vs k and bk vs k for different l values
figure;
plot(1:length(interv1a), interv1a, '-o', 'DisplayName', 'a_k(l=0.003)');
hold on;
plot(1:length(interv1b), interv1b, '-o', 'DisplayName', 'b_k(l=0.003)');
plot(1:length(a21), a21, '-o', 'DisplayName', 'a_k(l=0.1)');
plot(1:length(b21), b21, '-o', 'DisplayName', 'b_k(l=0.1)');

xlabel('Iteration k');
ylabel('Interval Boundaries');
legend;
title('Bisection Method for f1(x): ak and bk against iteration k');
grid on;

% Plot ak vs k and bk vs k for different l values
figure;
plot(1:length(interv2a), interv2a, '-o', 'DisplayName', 'a_k(l=0.003)');
hold on;
plot(1:length(interv2b), interv2b, '-o', 'DisplayName', 'b_k(l=0.003)');
plot(1:length(a22), a22, '-o', 'DisplayName', 'a_k(l=0.1)');
plot(1:length(b22), b22, '-o', 'DisplayName', 'b_k(l=0.1)');

xlabel('Iteration k');
ylabel('Interval Boundaries');
legend;
title('Bisection Method for f2(x): ak and bk against iteration k');
grid on;

% Plot ak vs k and bk vs k for different l values
figure;
plot(1:length(interv3a), interv3a, '-o', 'DisplayName', 'a_k(l=0.003)');
hold on;
plot(1:length(interv3b), interv3b, '-o', 'DisplayName', 'b_k(l=0.003)');
plot(1:length(a23), a23, '-o', 'DisplayName', 'a_k(l=0.1)');
plot(1:length(b23), b23, '-o', 'DisplayName', 'b_k(l=0.1)');

xlabel('Iteration k');
ylabel('Interval Boundaries');
legend;
title('Bisection Method for f3(x): ak and bk against iteration k');
grid on;

%-----------------------------------------------BISECTION FUNCTION------------------------------------------------------------------------------%
function [a, b, f_count] = bisection_(f, in_a, in_b, l, e)
    % Initializing the variables
    a(1) = in_a; 
    b(1) = in_b;
    f_count = 0;                    
    k = 1;    
    x1 = [];
    x2 = [];

    while (((b(k) - a(k)) > l)&& (e <= (l / 2) - 0.0005))
        % Calculating the points x1 and x2 
        x1(k) = ((a(k) + b(k)) / 2 - e);
        x2(k) = ((a(k) + b(k)) / 2 + e);

        % Updating interval [a, b] based on function values
        if f(x1(k)) < f(x2(k))
            a(k+1) = a(k);
            b(k+1) = x2(k);
        else
            a(k+1) = x1(k);
            b(k+1) = b(k);
        end

        f_count = f_count + 2; % Increment counters
        k = k + 1; % Move to next iteration
    end
end