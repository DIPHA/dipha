%  Copyright 2014 IST Austria
%
%  Contributed by: Jan Reininghaus
%
%  This file is part of DIPHA.
%
%  DIPHA is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Lesser General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  DIPHA is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License
%  along with DIPHA.  If not, see <http://www.gnu.org/licenses/>.

function plot_persistence_diagram_density( filename, dimension, resolution )
    %% load data
    [dimensions, birth_values, death_values] = load_persistence_diagram( filename );
    selection = (dimensions == dimension);
    diagram_points = [ birth_values( selection ), death_values( selection ) ];
    
    %% compute density
    min_value = min( birth_values );
    max_value = max( death_values );
    normalized_diagram_points = ( diagram_points - min_value) ./ ( max_value - min_value );
    raster = round( normalized_diagram_points * (resolution - 1) + 1 );
    I = raster(:,1);
    J = raster(:,2);
    S = ones( size(I,1), 1 );
    diagram_density = full(sparse( I, J, S, resolution, resolution ));

    %% plot density
    imagesc( [min_value, max_value ], [min_value, max_value], diagram_density' );  
    axis xy
    axis vis3d  
end