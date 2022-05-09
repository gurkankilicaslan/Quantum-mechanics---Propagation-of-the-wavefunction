clc;
clear;
clear all;

L = 100;
k = 1/L^2;
m = 1; 
hbar = 1; 
N = 200;
dt = 0.25; 
dx = L/(N+1);
f0 = sin(pi*(1:N)'/(N+1)) + sin(2*pi*(1:N)'/(N+1)) + sin(3*pi*(1:N)'/(N+1));
f0 = f0/norm(f0)/sqrt(dx);
f = f0;
V = zeros(N,N); 
X = zeros(N,N); 
P = zeros(N,N); 
T = zeros(N,N);
X(1,1) = 1/(N+1)*L;
P(1,2) = hbar/2i/dx; 
T(1,1:2) = [2/dx^2*hbar^2/2/m -1/dx^2*hbar^2/2/m];
for n = 2 : N - 1
    X(n,n) = n/(N+1)*L;
    P(n,n+1) = hbar/2i/dx;
    P(n,n-1) = -hbar/2i/dx;
    T(n,n-1) = -1/dx^2*hbar^2/2/m;
    T(n,n) = 2/dx^2*hbar^2/2/m;
    T(n,n+1) = -1/dx^2*hbar^2/2/m;
end
X(N,N) = N/(N+1)*L;
P(N,N-1) = -hbar/2i/dx;
T(N,N-1:N) = [-1/dx^2*hbar^2/2/m 2/dx^2*hbar^2/2/m];
V(1,1) = k*(X(1,1)-X(ceil(N/2),ceil(N/2)))^2/2;
for n = 2 : N - 1
    V(n,n) = k*(X(n,n)-X(ceil(N/2),ceil(N/2)))^2/2;
end
V(N,N) = k*(X(N,N)-X(ceil(N/2),ceil(N/2)))^2/2;
H = T + V; 
U = H/hbar/1i; 

fig = figure;
set(fig,'KeyPressFcn',@dummy);
for t = 1 : 1000000
    
    k1 = dt*U*f;
    k2 = dt*U*(f + k1/2);
    k3 = dt*U*(f + k2/2);
    k4 = dt*U*(f + k3);
    f = f + (k1 + 2*k2 + 2*k3 + k4)/6;

    if (mod(t,100) == 1)
        if (fig.CurrentCharacter == 'h') 
            fig.CurrentCharacter = '~';
            [f,value] = measure(H,f,dx);
            display(['Measured energy is ' num2str(value/(hbar*sqrt(k/m)))]);
        elseif (fig.CurrentCharacter == 'x')
            fig.CurrentCharacter = '~';
            [f,value] = measure(X,f,dx);
            display(['Measured position is ' num2str(value)]);
        elseif (fig.CurrentCharacter == 'p')
            fig.CurrentCharacter = '~';
            [f,value] = measure(P,f,dx);
            display(['Measured momentum is ' num2str(value)]);
        elseif (fig.CurrentCharacter == 'r')
            fig.CurrentCharacter = '~';
            f = f0;
            disp('Wavefunction was reset'); 
        end
        plot(L*(0:N+1)/(N+1),[0; conj(f).*f; 0],"k",'LineWidth',3);
        axis([0 L 0 0.09]);
        set(gca,'FontSize',20);
        xlabel('Position');
        ylabel('|\Psi|^2');
        title("Analytical solution-propagation of the wavefunction")
        grid on
        drawnow;
    end
end

function [f,value] = measure(Operator,f,dx) 
[evec,eval] = eig(Operator);
c = evec'*f;
c2 = conj(c).*c/(c'*c);
r = rand; 
s = 0;
for n = 1 : size(f,1)
    s = s + c2(n);
    if (s > r)
        break
    end
end
f = evec(:,n);
f = f/norm(f)/sqrt(dx);
value = eval(n,n);
end

function dummy(~,~,~)
end
