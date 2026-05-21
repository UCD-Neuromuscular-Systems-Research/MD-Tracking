function [velocity, acceleration, curvature, unitTangent, unitNormal] = lineParameters(line,dl)
%lineParameters - Computes several features of a line
%
% Syntax:  [velocity, acceleration, curvature, unitTangent, unitNormal] = lineParameters(linePoints,spatialIncrement)
%
% Inputs:
%	linePoints - coordinates of points defining the line
%	spatialIncrement - length of line segments [m]
%
% Outputs:
%	velocity - gradient of line (vectors composing the line)
%	acceleration
%	curvature
%	unitTangentVector
%	unitNormalVector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% External library required: none
%
% See also: parameters.m

% Author: Diego P. Botelho, Ph.D.
% email address: diegobotelho@gmail.com 
% Last revision: 30 September 2018


velocity = gradient(line');
gradGradLine = gradient(velocity);

unitTangent = velocity./sqrt(sum(velocity.^2,1));
acceleration = gradGradLine./dl.^2;

% curvature = sqrt(sum(velocity.^2,1) .* sum(gradGradLine.^2,1)...
%     - sum(velocity.*gradGradLine,1).^2)./sum(velocity.^2,1).^(3/2);

% curvature = sqrt(sum((velocity.*gradGradLine).^2,1))./sqrt(sum(velocity.^2,1)).^3 ;

% curvature = sqrt(sum(gradGradLine.^2./(1+velocity.^2).^(3/2),1)) ;

curvature = 1./sqrt(sum(velocity.^2,1)).*sqrt(sum(gradient(unitTangent).^2,1)) ;

unitNormal = acceleration./curvature;