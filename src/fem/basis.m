%***********************************************
% FEMAT — Finite Element Analysis in MATLAB
%***********************************************

function [quad_pts, weights, shape_func, shape_grads_ref] = basis(eletype)
% Gauss quadrature and shape functions for QUAD4, TRI3, HEX8, TET4

if eletype == "QUAD4"
    g = 1 / sqrt(3);
    quad_pts = [-g -g; -g g; g -g; g g];
    weights = [1; 1; 1; 1];

    shape_func = zeros(4, 4);
    shape_grads_ref = zeros(4, 2, 4);
    for i = 1:4
        xi  = quad_pts(i, 1);
        eta = quad_pts(i, 2);
        shape_func(i, :) = [(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), ...
                            (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)] / 4;
        shape_grads_ref(:, :, i) = [-(1 - eta), -(1 - xi);
                                     (1 - eta), -(1 + xi);
                                     (1 + eta),  (1 + xi);
                                    -(1 + eta),  (1 - xi)] / 4;
    end

elseif eletype == "TRI3"
    quad_pts = [1/3, 1/3];
    weights = 1/2;

    shape_func = zeros(1, 3);
    shape_grads_ref = zeros(3, 2, 1);
    xi  = quad_pts(1, 1);
    eta = quad_pts(1, 2);
    shape_func(1, :) = [1 - xi - eta, xi, eta];
    shape_grads_ref(:, :, 1) = [-1, -1; 1, 0; 0, 1];

elseif eletype == "HEX8"
    g = 1 / sqrt(3);
    quad_pts = [-g -g -g; -g  g -g;  g -g -g;  g  g -g; ...
                -g -g  g; -g  g  g;  g -g  g;  g  g  g];
    weights = [1; 1; 1; 1; 1; 1; 1; 1];

    shape_func = zeros(8, 8);
    shape_grads_ref = zeros(8, 3, 8);
    for i = 1:8
        xi   = quad_pts(i, 1);
        eta  = quad_pts(i, 2);
        zeta = quad_pts(i, 3);
        shape_func(i, :) = [(1 - xi) * (1 - eta) * (1 - zeta), ...
                            (1 + xi) * (1 - eta) * (1 - zeta), ...
                            (1 + xi) * (1 + eta) * (1 - zeta), ...
                            (1 - xi) * (1 + eta) * (1 - zeta), ...
                            (1 - xi) * (1 - eta) * (1 + zeta), ...
                            (1 + xi) * (1 - eta) * (1 + zeta), ...
                            (1 + xi) * (1 + eta) * (1 + zeta), ...
                            (1 - xi) * (1 + eta) * (1 + zeta)] / 8;
        shape_grads_ref(:, :, i) = ...
            [-(1 - eta) * (1 - zeta), -(1 - xi) * (1 - zeta), -(1 - xi) * (1 - eta);
              (1 - eta) * (1 - zeta), -(1 + xi) * (1 - zeta), -(1 + xi) * (1 - eta);
              (1 + eta) * (1 - zeta),  (1 + xi) * (1 - zeta), -(1 + xi) * (1 + eta);
             -(1 + eta) * (1 - zeta),  (1 - xi) * (1 - zeta), -(1 - xi) * (1 + eta);
             -(1 - eta) * (1 + zeta), -(1 - xi) * (1 + zeta),  (1 - xi) * (1 - eta);
              (1 - eta) * (1 + zeta), -(1 + xi) * (1 + zeta),  (1 + xi) * (1 - eta);
              (1 + eta) * (1 + zeta),  (1 + xi) * (1 + zeta),  (1 + xi) * (1 + eta);
             -(1 + eta) * (1 + zeta),  (1 - xi) * (1 + zeta),  (1 - xi) * (1 + eta)] / 8;
    end

elseif eletype == "TET4"
    quad_pts = [1/4, 1/4, 1/4];
    weights = 1/6;

    shape_func = zeros(1, 4);
    shape_grads_ref = zeros(4, 3, 1);
    xi   = quad_pts(1, 1);
    eta  = quad_pts(1, 2);
    zeta = quad_pts(1, 3);
    shape_func(1, :) = [1 - xi - eta - zeta, xi, eta, zeta];
    shape_grads_ref(:, :, 1) = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1];

else
    error('Unknown element type: %s', eletype)
end

end
