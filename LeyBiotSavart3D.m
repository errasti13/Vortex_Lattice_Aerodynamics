%% GRUPO 13 Alonso Lucas, Sara; Errasti Odriozola, Jon; 
%%          Sarabia Vargas, Alejandro; Terreros Sanchez, Carlos
%           Función para la aplicacion de la Ley de Biot-Savart


function [Velocidad_z,Velocidad_y,Velocidad_x]=LeyBiotSavart3D(xr,yr,zr,xP,yP,zP,xQ,yQ,zQ)

% Esta funcion calcula el factor de influencia 'w' que un hilo de torbellinos
% formado por un segmento 'PQ' genera en un punto arbitrario 'r'.

% Una vez enunciada la Ley de Biot-Savart y establecer las distintas
% relaciones para llegar a la formula final, definimos los modulos de los
% vectores que se utilizaran para calcular el factor de influencia:
% vector rrP y vector rrQ.

modulo_rrP=sqrt((xr-xP).^2 +(yr-yP).^2+(zr-zP).^2);
modulo_rrQ=sqrt((xr-xQ).^2 +(yr-yQ).^2+(zr-zQ).^2);

% Se procedera a calcular a continuacion los productos escalares
% necesarios para la obtencion del factor de influencia: rPQ . rrP y rPQ .
% rrQ.

prod_escalar_rPQrrP= (xQ-xP).*(xr-xP)+(yQ-yP).*(yr-yP)+(zQ-zP).*(zr-zP);
prod_escalar_rPQrrQ= (xQ-xP).*(xr-xQ)+(yQ-yP).*(yr-yQ)+(zQ-zP).*(zr-zQ);

% El siguiente paso sera calcular el producto vectorial entre rrP y rrQ y
% su respectivo modulo; como z=0, solo habra componente en z del producto
% vectorial.

rrPxrrQ_z=(xr-xP).*(yr-yQ)-(xr-xQ).*(yr-yP); %Componente z
rrPxrrQ_x=(zr-zP).*(yr-yQ)-(zr-zQ).*(yr-yP); %Componente x
rrPxrrQ_y=(xr-xP).*(zr-zQ)-(xr-xQ).*(zr-zP); %Componente y

% Por ultimo, se tiene todo lo necesario para obtener el factor de
% influencia.

modulo_rrPxrrQ_2=(rrPxrrQ_z.^2+rrPxrrQ_x.^2+rrPxrrQ_y.^2);

Velocidad_z=1/4/pi.*rrPxrrQ_z./(modulo_rrPxrrQ_2+1e-12).*...
(prod_escalar_rPQrrP./modulo_rrP-prod_escalar_rPQrrQ./modulo_rrQ);
Velocidad_y=1/4/pi.*rrPxrrQ_y./(modulo_rrPxrrQ_2+1e-12).*...
(prod_escalar_rPQrrP./modulo_rrP-prod_escalar_rPQrrQ./modulo_rrQ);
Velocidad_x=1/4/pi.*rrPxrrQ_x./(modulo_rrPxrrQ_2+1e-12).*...
(prod_escalar_rPQrrP./modulo_rrP-prod_escalar_rPQrrQ./modulo_rrQ);

end