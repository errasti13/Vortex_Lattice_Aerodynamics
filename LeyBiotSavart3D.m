% calculateVortexInfluence - Calculate the velocity components induced by a vortex segment
%
% Syntax:
%   [Velocidad_z, Velocidad_y, Velocidad_x] = calculateVortexInfluence(xr, yr, zr, xP, yP, zP, xQ, yQ, zQ)
%
% Inputs:
%   xr, yr, zr - Coordinates of the point where influence velocity is calculated
%   xP, yP, zP - Coordinates of the vortex segment start point
%   xQ, yQ, zQ - Coordinates of the vortex segment end point
%
% Outputs:
%   Velocidad_z - Influence velocity in the z direction at the specified point
%   Velocidad_y - Influence velocity in the y direction at the specified point
%   Velocidad_x - Influence velocity in the x direction at the specified point
%
% Description:
%   This function calculates the velocity components induced by a vortex
%   segment formed by points P and Q on a given point (xr, yr, zr) using
%   the Biot-Savart law. The resulting velocities are provided in the z, y,
%   and x directions.
%
% Example:
%   % Calculate the influence velocity components at a point due to a vortex segment
%   [Velocidad_z, Velocidad_y, Velocidad_x] = calculateVortexInfluence(xr, yr, zr, xP, yP, zP, xQ, yQ, zQ);
%
% Authors: Alonso Lucas, Sara; Errasti Odriozola, Jon;
%          Sarabia Vargas, Alejandro; Terreros Sanchez, Carlos


function [Velocidad_z, Velocidad_y, Velocidad_x] = calculateVortexInfluence(xr, yr, zr, xP, yP, zP, xQ, yQ, zQ)
    % Calculate vortex influence velocity components on a point
    
    % Calculate distances
    modulo_rrP = sqrt((xr - xP).^2 + (yr - yP).^2 + (zr - zP).^2);
    modulo_rrQ = sqrt((xr - xQ).^2 + (yr - yQ).^2 + (zr - zQ).^2);
    
    % Calculate scalar products
    prod_escalar_rPQrrP = (xQ - xP) .* (xr - xP) + (yQ - yP) .* (yr - yP) + (zQ - zP) .* (zr - zP);
    prod_escalar_rPQrrQ = (xQ - xP) .* (xr - xQ) + (yQ - yP) .* (yr - yQ) + (zQ - zP) .* (zr - zQ);
    
    % Calculate vector cross product components
    rrPxrrQ_z = (xr - xP) .* (yr - yQ) - (xr - xQ) .* (yr - yP); % Component z
    rrPxrrQ_x = (zr - zP) .* (yr - yQ) - (zr - zQ) .* (yr - yP); % Component x
    rrPxrrQ_y = (xr - xP) .* (zr - zQ) - (xr - xQ) .* (zr - zP); % Component y
    
    % Calculate the square of the vector cross product magnitude
    modulo_rrPxrrQ_2 = (rrPxrrQ_z.^2 + rrPxrrQ_x.^2 + rrPxrrQ_y.^2);
    
    % Calculate influence velocities
    small_value = 1e-12;
    Velocidad_z = 1 / (4 * pi) * rrPxrrQ_z ./ (modulo_rrPxrrQ_2 + small_value) .* ...
        (prod_escalar_rPQrrP ./ modulo_rrP - prod_escalar_rPQrrQ ./ modulo_rrQ);
    Velocidad_y = 1 / (4 * pi) * rrPxrrQ_y ./ (modulo_rrPxrrQ_2 + small_value) .* ...
        (prod_escalar_rPQrrP ./ modulo_rrP - prod_escalar_rPQrrQ ./ modulo_rrQ);
    Velocidad_x = 1 / (4 * pi) * rrPxrrQ_x ./ (modulo_rrPxrrQ_2 + small_value) .* ...
        (prod_escalar_rPQrrP ./ modulo_rrP - prod_escalar_rPQrrQ ./ modulo_rrQ);
end
