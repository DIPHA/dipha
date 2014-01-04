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

function plot_persistence_diagram( filename )
    %% load data
    [dimensions, birth_values, death_values] = load_persistence_diagram( filename );

    %% select a subset of points in the diagram based on lifetime
    persistence_threshold = 0.0;
    selection = ( death_values - birth_values ) > persistence_threshold;

    %% draw points
    colormap( jet( 13 ) );
    scatter( birth_values( selection ), death_values( selection ), 50, dimensions( selection ), 'filled', 'MarkerEdgeColor','k' );

    %% draw diagonal
    if( all( ~selection ) )
        min_value = -1;
        max_value = 1;
    else
        min_value = min( birth_values( selection ) );
        max_value = max( death_values( selection ) );
    end
    line( [min_value, max_value], [min_value, max_value], 'Color','k' );

    %% restrict axis
    safe_min = min_value - eps( min_value );
    safe_max = max_value + eps( max_value );
    axis( [safe_min, safe_max, safe_min, safe_max] );
    axis vis3d  
end