function new_X = pz_rescale(X, new_min, new_max)
    % Rescales a matrix of numbers to a new minimum and maximum
    %
    % Inputs:
    % X - numbers to rescale
    % new_min, new_max - rescaled range
    %
    % Outputs:
    % new_X - rescaled version of X
    %
    % ---------------------------------------------------------------------
    % Copyright (C) 2012 Peter Zeidman
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <http://www.gnu.org/licenses/>.   
    % ---------------------------------------------------------------------
    
    current_min = min(X(:));
    current_max = max(X(:));
    
    scale_factor = (current_max - current_min) / (new_max - new_min);
    new_X = new_min + (X - current_min) / scale_factor;    
end