%% Vortex-Lattice Method for Incompressible Aerodynamics Analysis

% Authors: Alonso Lucas, Sara; Errasti Odriozola, Jon;
%          Sarabia Vargas, Alejandro; Terreros Sanchez, Carlos

% Description:
% This MATLAB code implements the Vortex-Lattice method to analyze the
% incompressible steady-state aerodynamics of a given wing geometry. The
% code calculates and visualizes aerodynamic coefficients, lift and drag
% distributions, and other relevant parameters for a specified range of
% angles of attack (alpha) or sweep angles (psi). Users can customize the
% angle sweep and analyze the aerodynamic behavior of different wing
% geometries.

% Usage:
% 1. Configure the parameters and ranges for analysis in the 'main' function.
% 2. Run the script to perform aerodynamic calculations and generate results.

% Note:
% This code is part of a collaborative project for [Course/Research Name].
% Please refer to the respective documentation for detailed insights and
% considerations for the Vortex-Lattice method and its application.

% Dependencies:
% - LeyBiotSavart3D: A custom function for vortex panel integration.

% Last Updated: June 2020

clear

%% INPUT (MENU)

b = 14;     % Wingspan
cr = 5;     % Root chord
E = 0.5;    % Wing taper ratio (0-1)
psi = 40*pi/180;    % Wing sweep angle (radians)
alphag = linspace(0,10,11)*pi/180;    % Angle of attack (radians)
U_inf = 200;    % Flight velocity (m/s)
Ny = 20;    % Number of divisions along y-axis
Nx = 10;    % Number of divisions along x-axis
h = 0;      % Flight altitude (m)

%% Air Parameters - ISA Conditions

Tamb = 288.15 - 0.0065 * h;    % Ambient temperature (K)
Pamb = 101325 * (Tamb / 288.15)^(9.81 / (287.075 * 0.0065));    % Ambient pressure (Pa)
ro = Pamb / (287.074 * Tamb);    % Air density (kg/m^3)

%% PRELIMINARY CALCULATIONS

ct = E * cr;    % Tip chord
b2 = b / 2;     % Half wingspan
S = cr * b * (1 + E) / 2;    % Wing area
AR = b^2 / S;    % Aspect ratio
cp = S / b;     % Mean chord

%% NACA 2412

m = 2 / 100;
p = 4 / 10;
e = 12 / 100;
%% DESARROLLO

for k=1:length(alphag)
        alpha=alphag(k);
       %  for l=1:length(psig)                  
       % psi=psig(k);
        
    
% Cuerda en fcion. de y

for i=1:Ny+1
    yt(i)=-b/2+(i-1)*b/Ny;
    c_y(i)=((ct-cr)/b2)*abs(yt(i))+cr;
end

%% Leading and trailing edge calculation
x_ataque = tan(psi) * abs(yt);
x_salida = x_ataque + c_y;

%% Primary panel meshing
for j = 1:Ny+1
    for i = 1:Nx+1
        x_v(i, j) = x_ataque(j) + (x_salida(j) - x_ataque(j)) * (i - 1) / Nx;
        y_v(i, j) = yt(j);
        c_v(i, j) = c_y(j);
    end
end

%% Vertices of panels 1234 calculation
G = 1;
for i = 1:Nx
    for j = 1:Ny
        x1(G) = x_v(i, j); y1(G) = y_v(i, j);
        x2(G) = x_v(i, j + 1); y2(G) = y_v(i, j + 1);
        x3(G) = x_v(i + 1, j + 1); y3(G) = y_v(i + 1, j + 1);
        x4(G) = x_v(i + 1, j); y4(G) = y_v(i + 1, j);
        G = G + 1;
    end
end

% Calculate average chord and placement point for each panel
N = Nx * Ny;

for i = 1:N
    cm(i) = ((x3(i) + x4(i)) - (x1(i) + x2(i))) / 2;
    xp(i) = (x2(i) + x1(i)) / 2;
    xp_2(i) = (x3(i) + x4(i)) / 2;
    xptotal(i) = xp(i) + 0.75 * (xp_2(i) - xp(i));
    yp(i) = (y2(i) + y1(i)) / 2;
    yp_2(i) = (y3(i) + y4(i)) / 2;
    yptotal(i) = yp(i) + 0.75 * (yp_2(i) - yp(i));
    cmedia(i) = sqrt((xp(i) - xp_2(i)).^2 + (yp(i) - yp_2(i)).^2);
    c_yp(i) = ((ct - cr) / b2) * abs(yptotal(i)) + cr;
end

for i = 1:Nx
    for j = 1:Ny
        G = (i - 1) * Ny + j;
        Xptotal(i, j) = xptotal(G);
        Yptotal(i, j) = yptotal(G);
        C_yp(i, j) = c_yp(G);
    end
end

% Calculate profile for placement points
Xp_total = Xptotal - abs(Yptotal) .* tan(psi);
Xcurvatura = Xp_total ./ C_yp;

for j = 1:Ny
    for i = 1:Nx
        if Xcurvatura(i, j) <= p
            Dx(i, j) = (2 * m) / p^2 * (p - Xcurvatura(i, j));
        else
            Dx(i, j) = ((2 * m) / (1 - p)^2) * (p - Xcurvatura(i, j));
        end
        theta1(i, j) = Dx(i, j);
    end 
end

for i = 1:Nx
    for j = 1:Ny
        G = (i - 1) * Ny + j;
        thetaC(G) = theta1(i, j);
    end
end


%% Calculate x and y coordinates for A, B, C, D points
xAtotal = x1 + (x4 - x1) / 4;
yAtotal = y1;

xBtotal = x2 + (x3 - x2) / 4;
yBtotal = y2;

xCtotal = xBtotal(Ny+1:Nx*Ny);
yCtotal = yBtotal;

xCtotal(Nx*Ny-Ny+1:Nx*Ny) = x_salida(2:Ny+1) + 1/4;

xDtotal = xAtotal(Ny+1:Nx*Ny);
yDtotal = yAtotal;

xDtotal(Nx*Ny-Ny+1:Nx*Ny) = x_salida(1:Ny) + 1/4;

xCtotal(Nx*Ny+1:Nx*Ny+Ny) = 100 * cr;
xDtotal(Nx*Ny+1:Nx*Ny+Ny) = 100 * cr;

xAtotal(Nx*Ny+1:Nx*Ny+Ny) = xDtotal(Nx*Ny-Ny+1:Nx*Ny);
xBtotal(Nx*Ny+1:Nx*Ny+Ny) = xCtotal(Nx*Ny-Ny+1:Nx*Ny);

yAtotal(Nx*Ny+1:Nx*Ny+Ny) = yAtotal(1:Ny);
yBtotal(Nx*Ny+1:Nx*Ny+Ny) = yBtotal(1:Ny);
yCtotal(Nx*Ny+1:Nx*Ny+Ny) = yCtotal(1:Ny);
yDtotal(Nx*Ny+1:Nx*Ny+Ny) = yDtotal(1:Ny);


%% Calculations of normals
normal_P_z = ((xptotal(1:Nx*Ny) - x1(1:Nx*Ny)) .* (yptotal(1:Nx*Ny) - y2(1:Nx*Ny))) ...
           - ((xptotal(1:Nx*Ny) - x2(1:Nx*Ny)) .* (yptotal(1:Nx*Ny) - y1(1:Nx*Ny)));

normal_z = normal_P_z ./ normal_P_z;

%% Calculations of influence factors for segments ABCD

N = Nx * Ny;
Velocidad_total = zeros(1, N);
z = zeros(1, N);

for i = 1:N
    xcoloc = xptotal(i);
    ycoloc = yptotal(i);
    zcoloc = z(i);
    
    [Velocidad_z_AB(1:Nx*Ny + Ny)] = LeyBiotSavart3D(xcoloc, ycoloc, zcoloc, xAtotal(1:Nx*Ny + Ny), yAtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny), xBtotal(1:Nx*Ny + Ny), yBtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny));
    [Velocidad_z_BC(1:Nx*Ny + Ny)] = LeyBiotSavart3D(xcoloc, ycoloc, zcoloc, xBtotal(1:Nx*Ny + Ny), yBtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny), xCtotal(1:Nx*Ny + Ny), yCtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny));
    [Velocidad_z_CD(1:Nx*Ny + Ny)] = LeyBiotSavart3D(xcoloc, ycoloc, zcoloc, xCtotal(1:Nx*Ny + Ny), yCtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny), xDtotal(1:Nx*Ny + Ny), yDtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny));
    [Velocidad_z_DA(1:Nx*Ny + Ny)] = LeyBiotSavart3D(xcoloc, ycoloc, zcoloc, xDtotal(1:Nx*Ny + Ny), yDtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny), xAtotal(1:Nx*Ny + Ny), yAtotal(1:Nx*Ny + Ny), zeros(1, Nx*Ny + Ny));
    
    Velocidad_z_total = Velocidad_z_AB + Velocidad_z_BC + Velocidad_z_CD + Velocidad_z_DA;
    VelocidadBCDA = Velocidad_z_BC + Velocidad_z_DA;
    
    Velocidad_total(i, 1:Nx*Ny + Ny) = Velocidad_z_total .* ones(1, Nx*Ny + Ny);
    Matriz_bprima1(i, 1:Nx*Ny + Ny) = VelocidadBCDA .* ones(1, Nx*Ny + Ny);
end


Matriz_A_aux = Velocidad_total;

Matriz_A(:, 1:Nx*Ny) = Matriz_A_aux(:, 1:Nx*Ny);
Matriz_A(:, Nx*Ny - Ny + 1:Nx*Ny) = Matriz_A(:, Nx*Ny - Ny + 1:Nx*Ny) + Matriz_A_aux(:, Nx*Ny + 1:Nx*Ny + Ny);

Matriz_b = -(U_inf * alpha) * normal_z(1:N) + (U_inf * thetaC .* normal_z(1:N));

Matriz_Gamma = Matriz_A \ Matriz_b'; % Circulation density matrix

Matriz_bprima(:, 1:Nx*Ny) = Matriz_bprima1(:, 1:Nx*Ny);

for i = 1:Nx
    for j = 1:Ny
        I = (i - 1) * Ny + j;
        
        % Calculate area of panels
        Area(i, j) = 0.5 * (x4(I) - x1(I) + x3(I) - x2(I)) * (y2(I) - y1(I));
        Areavector(I) = Area(i, j);
        
        % Calculate average chord and related quantities
        h(i, j) = cm(I);
        hyp(i, j) = c_yp(I);
        xcentrogamma(i, j) = 0.5 * (xAtotal(I) + xBtotal(I));
        ycentrogamma(i, j) = 0.5 * (yAtotal(I) + yBtotal(I));
    end
end

% Calculate variables for each bar
ymed(1:Ny) = 0.5 * (y_v(1, 2:(Ny+1)) + y_v(1, 1:Ny));
cy(1:Ny) = 0.5 * (x_v(Nx+1, 2:(Ny+1)) + x_v(Nx+1, 1:Ny) - x_v(1, (2:(Ny+1))) - x_v(1, 1:Ny));

xab(1:Ny) = 0.5 * (x_ataque(2:Ny+1) + x_ataque(1:Ny));
xsb(1:Ny) = 0.5 * (x_salida(2:Ny+1) + x_salida(1:Ny));
cb = xsb - xab;


% Calculate circulation densities
for i = 1:Nx
    for j = 1:Ny
        G = (i - 1) * Ny + j;
        Gamma(i, j) = Matriz_Gamma(G, 1);
    end
end

for j = 1:Ny
    Gamma_total(1, j) = Gamma(1, j) / h(1, j);
    
    for i = 2:Nx
        Gamma_total(i, j) = (Gamma(i, j) - Gamma(i - 1, j)) / h(i, j);
    end
end

for j = 1:Ny
    Circulacion(j) = Gamma_total(1:Nx, j)' * h(1:Nx, j);
end

% Calculate induced velocities
Velocidad_inducida = Matriz_bprima * Matriz_Gamma;

for i = 1:Nx
    for j = 1:Ny
        I = (i - 1) * Ny + j;
        Velocidad_inducidap(i, j) = Velocidad_inducida(I);
    end
end

% Calculate incDprima
for i = 1:Nx
    for j = 1:Ny
        incDprima(i, j) = Gamma_total(i, j) * Velocidad_inducidap(i, j);
    end
end


%% Calculos de la sustentacion

%Distribucion de presiones
P=ro*U_inf*Gamma_total;

%Incremento de Sustentacion
incL(1:Nx,1:Ny)=ro*U_inf*(Gamma_total.*Area(1:Nx,1:Ny));
%Incremento de Resistencia
incD=-ro*incDprima(1:Nx,1:Ny).*(Area(1:Nx,1:Ny));
%Distribucion de Sustentacion
distr_L=ro*U_inf*Circulacion;

%Coeficiente de presiones
Cp=-2*Gamma_total/U_inf;
%Sustentacion
L1=sum(sum(incL));
D=sum(sum(incD));
L=sum(sum(Area(1:Nx,1:Ny).*P));

AreaTotal=sum(sum(Area(1:Nx,1:Ny))); %Superficie del ala completa

%Coeficientes aerodinamicos
Cl(k)=L/(1/2*ro*U_inf^2*AreaTotal);  %(k)
Cd(k)=D/(1/2*ro*U_inf^2*AreaTotal);  %(k)

distr_Cl=distr_L./(0.5*ro*U_inf^2*cb);

cly=2*Circulacion./(U_inf*cy);

%Momento de cabeceo

   for i=1:Nx
         for j=1:Ny
             I=(i-1)*Ny+j;
             xptotalp(i,j)=xptotal(I);
             yptotalp(i,j)=yptotal(I);
         end
     end
  
     for i=1:Nx
            for j=1:Ny
                I=(i-1)*Ny+j;
                xAtaque(I)=x_ataque(j)+(i-1)*(x_salida(j)-x_ataque(j))/Nx;
          
                Mya(I)=ro*U_inf*sum(sum(Gamma_total.*Area.*(xptotalp-xAtaque(I))));
                
            end
     end
    
  
      incM(1:Nx,1:Ny)=incL(1:Nx,1:Ny).*(xptotalp(1:Nx,1:Ny));
      Momya=-sum(sum(incM));
      CM(k)=Momya/(1/2*ro*U_inf^2*AreaTotal*cp);  %(k)
%% Representacion grafica

figure(14)
plot(xAtaque/cr,Mya)
%legend({'Momentos de cabeceo'},'Location','Northeast')
title('	Momentos de cabeceo')
hold on
legend ('\alpha = 0^o','\alpha = 2^o','\alpha = 4^o','\alpha = 5^o','\alpha = 6^o','\alpha = 8^o','\alpha = 10^o')
xlabel('Borde de ataque'),ylabel('Momento de cabeceo')
grid on

figure(15)
plot(ymed,distr_Cl)
legend({'Cl'},'Location','Southeast')
title('Coeficiente de sustentacion a lo largo de la envergadura')
hold on
xlabel('Envergadura'),ylabel('Coeficiente de sustentacion')
legend ('\alpha = 0^o','\alpha = 2^o','\alpha = 4^o','\alpha = 5^o','\alpha = 6^o','\alpha = 8^o','\alpha = 10^o')
grid on
    end  
   % end

% Define common style for all plots
lineStyle = '-o';
legendLocation = 'Best';
gridStyle = 'on';

% Plot CM vs Cl
figure(1)
plot(CM, Cl, lineStyle)
title('Cl vs CM')
xlabel('CM')
ylabel('Cl')
grid(gridStyle)
legend('Cl vs CM', 'Location', legendLocation)

% Plot Cl vs Cd
figure(2)
plot(Cd, Cl, lineStyle)
title('Cl vs Cd')
xlabel('Cd')
ylabel('Cl')
grid(gridStyle)
legend('Cl vs Cd', 'Location', legendLocation)

% Plot Cl vs Angle of Attack
figure(3)
plot(alphag * 180 / pi, Cl, lineStyle)
title('Cl vs Angle of Attack')
xlabel('\alpha (^o)')
ylabel('Cl')
grid(gridStyle)
legend('Cl vs Angle of Attack', 'Location', legendLocation)

% Plot CM vs Angle of Attack
figure(5)
plot(alphag * 180 / pi, CM, lineStyle)
title('CM vs Angle of Attack')
xlabel('\alpha (^o)')
ylabel('CM')
grid(gridStyle)
legend('CM vs Angle of Attack', 'Location', legendLocation)

% Plot Cd vs Angle of Attack
figure(7)
plot(alphag * 180 / pi, Cd, lineStyle)
title('Cd vs Angle of Attack')
xlabel('\alpha (^o)')
ylabel('Cd')
grid(gridStyle)
legend('Cd vs Angle of Attack', 'Location', legendLocation)

% Plot Moments of Pitch vs Chord
figure(9)
plot(xAtaque / cr, Mya, lineStyle)
title('Moments of Pitch')
xlabel('Normalized Chord Length (x/c)')
ylabel('Moment of Pitch')
grid(gridStyle)
legend('\alpha = 0^o', '\alpha = 2^o', '\alpha = 4^o', '\alpha = 5^o', '\alpha = 6^o', '\alpha = 8^o', '\alpha = 10^o', 'Location', legendLocation)

% Plot Lift Coefficient Distribution
figure(10)
plot(ymed, distr_Cl, lineStyle)
title('Lift Coefficient Distribution')
xlabel('Spanwise Position (y/b)')
ylabel('Lift Coefficient (Cl)')
grid(gridStyle)
legend('\alpha = 0^o', '\alpha = 2^o', '\alpha = 4^o', '\alpha = 5^o', '\alpha = 6^o', '\alpha = 8^o', '\alpha = 10^o', 'Location', legendLocation)

% Plot Lift Distribution
figure(11)
plot(ymed, distr_L, lineStyle)
title('Lift Distribution Along Span')
xlabel('Spanwise Position (y/b)')
ylabel('Lift (L)')
grid(gridStyle)
legend('Lift Distribution', 'Location', legendLocation)

     
    
% Plot for Cp Distribution as Mesh
figure(12)
mesh(xptotalp / cr, yptotalp / b2, Cp)
title('Cp Distribution Along the Wing')
xlabel('x/c')
ylabel('y/b')
zlabel('Cp')
legend('Cp', 'Location', 'southeast')
grid on

% Plot for Cp Distribution as Individual Surfaces
figure(13)
rotate3d on
colormap(hot);

for i = 1:Nx
    for j = 1:Ny
        h1 = surf(x_v(i:i+1, j:j+1), y_v(i:i+1, j:j+1), zeros(size(x_v(i:i+1, j:j+1))), Cp(i, j));
        c = colorbar;
        c.Label.String = 'Cp';
        title('Cp Distribution')
        xlabel('x')
        ylabel('y')
        zlabel('z')
        grid on
        axis equal
        hold on
    end
end

hold off
