%% GRUPO 13 Alonso Lucas, Sara; Errasti Odriozola, Jon; 
%%          Sarabia Vargas, Alejandro; Terreros Sanchez, Carlos

% Este codigo calcula la aerodinamica estacionaria incompresible para un
% ala con una geometria dada mediante el metodo Vortex-Lattice. El codigo
% esta predeterminado para calcular las graficas de los coeficientes
% aerodinamicos respecto al angulo de ataque, alpha. No obstante, el usuario
% puede cambiar el bucle de resolucion de las graficas al inicio del codigo
% y obtener las mismas graficas dependientes del angulo de flecha, psi.
% (linea 51-52 y Representacion a partir de la linea 333)
 

clear

%% INPUT (MENU)

b = 14;  %Envergadura alar
cr = 5;  %Cuerda en raiz del ala
E = 0.5;  %Valor estrechamiento del ala (0-1)
psi=40*pi/180; %Angulo flecha del ala (rad)  %psig=[0:5:40]*pi/180
alphag=[2]*pi/180;   %Angulo de ataque(rad)  %alpha=5*pi/180  adaptar este dato invariable para el estudio de variacion de psi
U_inf=200; %Velocidad de vuelo (m/s)
Ny = 20;  %Num. divisiones eje y en ala
Nx = 10;  %Num. divisiones eje x en ala
h=0;      %Altura de vuelo (m)

%% Parametros del aire - Condiciones ISA

Tamb = 288.15-0.0065*h; %Temperatura ambiental (K)
Pamb = 101325*(Tamb/288.15)^(9.81/(287.075*0.0065)); %Presion ambiental (Pa)
ro = Pamb/(287.074*Tamb); %Densidad del aire (kg/m^3)

%% CALCULOS PREVIOS

ct = E*cr;  %Cuerda en punta
b2 = b/2;  %Mitad de envergadura
S=cr*b*(1+E)/2;  %Superficie alar
AR=b^2/S;  %AR 
cp = S/b;  %Cuerda promedio

%% NACA 2412

m = 2/100;
p = 4/10;
e = 12/100;

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

%% Borde de ataque y borde salida

for i=1:Ny+1
    x_ataque(i)=tan(psi)*abs(yt(i));
    x_salida(i)=x_ataque(i)+c_y(i);
end

%% Mallado primario de paneles 1234

for j=1:Ny+1
    for i=1:Nx+1
        x_v(i,j)=x_ataque(j)+(x_salida(j)-x_ataque(j))*(i-1)/Nx;
        y_v(i,j)=yt(j);
        c_v(i,j)=c_y(j);
    end
    
end
    
%% Vertices de los paneles 1234

for i = 1:Nx
    for j = 1:Ny
  G=(i-1)*Ny+j;
 x1(G) = x_v(i,j); y1(G)=y_v(i,j);
 x2(G) = x_v(i,j+1); y2(G)=y_v(i,j+1); 
 x3(G) = x_v(i+1,j+1); y3(G)=y_v(i+1,j+1); 
 x4(G) = x_v(i+1,j); y4(G)=y_v(i+1,j); 
     end
end

%Cuerda media de cada panel y del punto de colocacion para cada panel

N=Nx*Ny;

for i=1:N
cm(i) = ((x3(i)+x4(i))-(x1(i)+x2(i)))/2;
xp(i) = (x2(i)+x1(i))./2;
xp_2(i) = (x3(i)+x4(i))./2;
xptotal(i) = xp(i)+0.75.*(xp_2(i)-xp(i));
yp(i) = (y2(i)+y1(i))./2;
yp_2(i) = (y3(i)+y4(i))./2;
yptotal(i)=yp(i)+0.75.*(yp_2(i)-yp(i));
cmedia(i)=sqrt((xp(i)-xp_2(i)).^2+(yp(i)-yp_2(i)).^2);
c_yp(i)=((ct-cr)/b2)*abs(yptotal(i))+cr;
end

for i=1:Nx
    for j=1:Ny
        G=(i-1)*Ny+j;
        Xptotal(i,j)=xptotal(G);
        Yptotal(i,j)=yptotal(G);
        C_yp(i,j)=c_yp(G);
        
    end
end
%Perfil para puntos de colocacion

Xp_total=Xptotal-abs(Yptotal).*tan(psi);
Xcurvatura=Xp_total./C_yp;

for j=1:Ny
    for i=1:Nx
        if Xcurvatura(i,j) <= p
           Dx(i,j)=(2*m)/p^2*(p-Xcurvatura(i,j));
        elseif Xcurvatura(i,j) > p
           Dx(i,j)=((2*m)/(1-p)^2)*(p-Xcurvatura(i,j));
        end
           theta1(i,j)=Dx(i,j);       
    end 
end

for i=1:Nx
    for j=1:Ny
        G=(i-1)*Ny+j;
        thetaC(G)=theta1(i,j);
    end
end
%% Vertices de los paneles ABCD

xAtotal(1:Nx*Ny) = x1(1:Nx*Ny)+(x4(1:Nx*Ny)-x1(1:Nx*Ny))./4; yAtotal(1:Nx*Ny)=y1(1:Nx*Ny); 

xBtotal(1:Nx*Ny) = x2(1:Nx*Ny)+(x3(1:Nx*Ny)-x2(1:Nx*Ny))./4; yBtotal(1:Nx*Ny)=y2(1:Nx*Ny); 

xCtotal(1:Nx*Ny-Ny) = xBtotal((Ny+1):Nx*Ny); yCtotal(1:Nx*Ny)=yBtotal(1:Nx*Ny);
xCtotal((Nx*Ny-Ny+1):Nx*Ny)=x_salida(2:Ny+1)+1/4;

xDtotal(1:Nx*Ny-Ny) = xAtotal((Ny+1):Nx*Ny); yDtotal(1:Nx*Ny)=yAtotal(1:Nx*Ny);
xDtotal((Nx*Ny-Ny+1):Nx*Ny)=x_salida(1:Ny)+1/4;

xCtotal(Nx*Ny+1:Nx*Ny+Ny)=100*cr;
xDtotal(Nx*Ny+1:Nx*Ny+Ny)=100*cr;

xAtotal((Nx*Ny+1:Nx*Ny+Ny))=xDtotal((Nx*Ny-Ny+1):Nx*Ny);
xBtotal((Nx*Ny+1:Nx*Ny+Ny))=xCtotal((Nx*Ny-Ny+1):Nx*Ny);

yAtotal((Nx*Ny+1:Nx*Ny+Ny))=yAtotal(1:Ny);
yBtotal((Nx*Ny+1:Nx*Ny+Ny))=yBtotal(1:Ny);
yCtotal((Nx*Ny+1:Nx*Ny+Ny))=yCtotal(1:Ny);
yDtotal((Nx*Ny+1:Nx*Ny+Ny))=yDtotal(1:Ny);

%%  Calculos de normales

normal_P_z = -(xptotal(1:Nx*Ny)-x1(1:Nx*Ny)).*(yptotal(1:Nx*Ny)-y2(1:Nx*Ny))...
+(xptotal(1:Nx*Ny)-x2(1:Nx*Ny)).*(yptotal(1:Nx*Ny)-y1(1:Nx*Ny));

normal_z=normal_P_z./normal_P_z;

%% Calculos del factor de influencia de los segmentos ABCD

N=Nx*Ny;
Velocidad_total=zeros(1,N);

%Particularizamos en z=0
z=zeros(1,N);

for i=1:N
    xcoloc=xptotal(i); ycoloc=yptotal(i); zcoloc=z(i);
    
    [Velocidad_z_AB(1:Nx*Ny+Ny)]=LeyBiotSavart3D(xcoloc,ycoloc,zcoloc,xAtotal(1:Nx*Ny+Ny),yAtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny),xBtotal(1:Nx*Ny+Ny),yBtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny));
    [Velocidad_z_BC(1:Nx*Ny+Ny)]=LeyBiotSavart3D(xcoloc,ycoloc,zcoloc,xBtotal(1:Nx*Ny+Ny),yBtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny),xCtotal(1:Nx*Ny+Ny),yCtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny));
    [Velocidad_z_CD(1:Nx*Ny+Ny)]=LeyBiotSavart3D(xcoloc,ycoloc,zcoloc,xCtotal(1:Nx*Ny+Ny),yCtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny),xDtotal(1:Nx*Ny+Ny),yDtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny));
    [Velocidad_z_DA(1:Nx*Ny+Ny)]=LeyBiotSavart3D(xcoloc,ycoloc,zcoloc,xDtotal(1:Nx*Ny+Ny),yDtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny),xAtotal(1:Nx*Ny+Ny),yAtotal(1:Nx*Ny+Ny),zeros(1,Nx*Ny+Ny));
    
    Velocidad_z_total=Velocidad_z_AB+Velocidad_z_BC+Velocidad_z_CD+Velocidad_z_DA;
    
    VelocidadBCDA=Velocidad_z_BC+Velocidad_z_DA;
    
    Velocidad_total(i,1:Nx*Ny+Ny)=Velocidad_z_total.*ones(1,Nx*Ny+Ny);
    Matriz_bprima1(i,1:Nx*Ny+Ny)=VelocidadBCDA.*ones(1,Nx*Ny+Ny);
end

Matriz_A_aux=Velocidad_total;

Matriz_A(:,1:Nx*Ny)=Matriz_A_aux(:,1:Nx*Ny);
Matriz_A(:,Nx*Ny-Ny+1:Nx*Ny)=Matriz_A(:,Nx*Ny-Ny+1:Nx*Ny)+Matriz_A_aux(:,Nx*Ny+1:Nx*Ny+Ny);


Matriz_b=-(U_inf*(alpha)*normal_z(1:N))+(U_inf*thetaC.*normal_z(1:N));

Matriz_Gamma=(Matriz_A\Matriz_b'); %Matriz de densidades de circulacion

Matriz_bprima(:,1:Nx*Ny)=Matriz_bprima1(:,1:Nx*Ny);

%% Calculos necesarios para la obtencionn de la circulacion del panel como Matriz (Nx,Ny)

for i=1:Nx
    for j=1:Ny
        I=(i-1)*Ny+j;
        %Area de los paneles
        Area(i,j)=0.5*(x4(I)-x1(I)+x3(I)-x2(I))*(y2(I)-y1(I));
        Areavector(I)=Area(i,j);
        %Matriz de cuerdas medias de los paneles
        h(i,j)=cm(I);
        hyp(i,j)=c_yp(I);
        xcentrogamma(i,j)=1/2*(xAtotal(I)+xBtotal(I));
        ycentrogamma(i,j)=1/2*(yAtotal(I)+yBtotal(I));
    end
end

%Calculos de variables de cada barra

ymed(1:Ny)=0.5*(y_v(1,2:(Ny+1))+y_v(1,1:Ny));
cy(1:Ny)=0.5*(x_v(Nx+1,2:(Ny+1))+x_v(Nx+1,1:Ny)-x_v(1,(2:(Ny+1)))-x_v(1,1:Ny));

xab(1:Ny)=(x_ataque(2:Ny+1)+x_ataque(1:Ny))*0.5;
xsb(1:Ny)=(x_salida(2:Ny+1)+x_salida(1:Ny))*0.5;
cb=(xsb-xab);

%% Calculo de las densidades de circulacion

for i=1:Nx
    for j=1:Ny
        G=(i-1)*Ny+j;
        Gamma(i,j)=Matriz_Gamma(G,1);
    end
end

for j=1:Ny
    
      Gamma_total(1,j)=Gamma(1,j)/h(1,j);
      
     for i=2:Nx
   
        Gamma_total(i,j)=(Gamma(i,j)-Gamma(i-1,j))/h(i,j);
    end
end


for j=1:Ny
Circulacion(j)=Gamma_total(1:Nx,j)'*h(1:Nx,j);
end



%Bucle para calculo del incD

Velocidad_inducida= Matriz_bprima*Matriz_Gamma;

for i=1:Nx
    for j=1:Ny
        I=(i-1)*Ny+j;
 Velocidad_inducidap(i,j)=Velocidad_inducida(I);
    end
end

for i=1:Nx
    for j=1:Ny
        
        incDprima(i,j)=Gamma_total(i,j)*Velocidad_inducidap(i,j);
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

figure(1)
plot(CM,Cl,'-o')
title('	Cl vs CM')
hold on
xlabel('CM'),ylabel('Cl')
grid on

figure(2)
plot(Cd,Cl,'-o')
title('	Cl vs Cd')
hold on
xlabel('Cd'),ylabel('Cl')
grid on


figure(3)
plot(alphag*180/pi,Cl,'-o')
title('	Cl vs Angulo de ataque')
hold on
xlabel('\alpha(^o)'),ylabel('Cl')
grid on

%figure(4)
%plot(psig*180/pi,Cl,'-o')
%title('	Cl vs Angulo de flecha')
%hold on
%xlabel('\psi(^o)'),ylabel('Cl')
%grid on

figure(5)
plot(alphag*180/pi,CM,'-o')
title('	CM vs Angulo de ataque')
hold on
xlabel('\alpha(^o)'),ylabel('CM')
grid on

%figure(6)
%plot(psig*180/pi,CM,'-o')
%title('	CM vs Angulo de flecha')
%hold on
%xlabel('\psi(^o)'),ylabel('CM')
%grid on

figure(7)
plot(alphag*180/pi,Cd,'-o')
title('	Cd vs Angulo de ataque')
hold on
xlabel('\alpha(^o)'),ylabel('Cd')
grid on

%figure(8)
%plot(psig*180/pi,Cd,'-o')
%title('	Cd vs Angulo de flecha')
%hold on
%xlabel('\psi(^o)'),ylabel('Cd')
%grid on

figure(9)
plot(xAtaque/cr,Mya)
%legend({'Momentos de cabeceo'},'Location','Northeast')
title('	Momentos de cabeceo')
hold on
legend ('\alpha = 0^o','\alpha = 2^o','\alpha = 4^o','\alpha = 5^o','\alpha = 6^o','\alpha = 8^o','\alpha = 10^o')
xlabel('Borde de ataque'),ylabel('Momento de cabeceo')
grid on

figure(10)
plot(ymed,distr_Cl)
legend({'Cl'},'Location','Southeast')
title('Coeficiente de sustentacion a lo largo de la envergadura')
hold on
xlabel('Envergadura'),ylabel('Coeficiente de sustentacion')
legend ('\alpha = 0^o','\alpha = 2^o','\alpha = 4^o','\alpha = 5^o','\alpha = 6^o','\alpha = 8^o','\alpha = 10^o')
grid on

  
figure(11)

plot(ymed,distr_L)
legend({'L'},'Location','Southeast')
title('Sustentacion a lo largo de la envergadura')
hold on
xlabel('Envergadura'),ylabel('Sustentacion')
grid on
     
    
 figure(12)
mesh(xptotalp/cr,yptotalp/b2,Cp)
legend({'Cp'},'Location','Southeast')
title('Distribucion de Cp a lo largo del ala')
xlabel('x'),ylabel('y'),zlabel('z')
grid on

figure(13)

rotate3d on
colormap(hot);

for i=1:Nx
    for j=1:Ny
h1=surf(x_v(i:i+1,j:j+1),y_v(i:i+1,j:j+1),zeros(size(x_v(i:i+1,j:j+1))),Cp(i,j));
c=colorbar;
c.Label.String='Cp';
title('Distribucion de Cp')
xlabel('x'),ylabel('y'),zlabel('z')
grid on
axis equal
hold on
    end
    hold on
end