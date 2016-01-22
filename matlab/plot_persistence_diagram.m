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

function plot_persistence_diagram( filename, persistence_threshold )
    %% default values
    if nargin < 2
        persistence_threshold = 0;
    end
    
    %% load data
    [dims, births, deaths] = load_persistence_diagram( filename );

    %% apply persistence_threshold
    thresholded_points = deaths - births > persistence_threshold;
    dims = dims( thresholded_points );
    births = births( thresholded_points );
    deaths = deaths( thresholded_points );

    %% classify points in diagram     
    essentials = (dims < 0);
    ordinaries = ~essentials;
    
    %% draw parameters
    max_dims = max(floor(abs(dims+0.5)))
    cur_colormap = lines( max_dims + 1 );
    colormap( cur_colormap );
    marker_size = 50;
    
    %% draw essential points
    scatter( births( essentials ), deaths( essentials ), marker_size, cur_colormap(-dims( essentials ) - 1 + 1,:), ...
             'filled', 'MarkerEdgeColor','k', 'Marker', 's' );
         
    %% draw non essential points
    hold on
    scatter( births( ordinaries ), deaths( ordinaries ), marker_size, cur_colormap(dims( ordinaries ) + 1,:), ...
             'filled', 'MarkerEdgeColor','k', 'Marker', 'o' );
    hold off
    
    %% draw diagonal
    min_value = min( births );
    max_value = max( deaths );
    line( [min_value, max_value], [min_value, max_value], 'Color','k' );

    %% draw axis
    safe_min = min_value - 10 * eps( min_value );
    safe_max = max_value + 10 * eps( max_value );
    axis square  
    axis( [safe_min, safe_max, safe_min, safe_max] );
    xticks = get(gca,'XTick');
    set(gca,'YTick',xticks);

    
    %% draw legend
    legend('Essential','Ordinary','Location','SouthEast');
    
    %% draw axis labels
    xlabel('Birth');
    ylabel('Death');

    %% draw colobar
    caxis([0 max_dims])
    ylabel(colorbar('YTick',0:1:max_dims), 'Dimension');
    
    %% draw title
    set( title(['Data: ' filename]), 'Interpreter', 'none' );
    
    %% draw grid
    grid on;
end
