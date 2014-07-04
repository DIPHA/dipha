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

function plot_midlife_persistence_diagram( filename, persistence_threshold )
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
    cur_colormap = lines( max( dims ) + 1 );
    colormap( cur_colormap );
    marker_size = 50;
    
    %% transform points
    midlife = (births + deaths) / 2;
    persistence = deaths - births;
    
    %% draw essential points
    scatter( midlife( essentials ), persistence( essentials ), marker_size, cur_colormap(-dims( essentials ) - 1 + 1,:), ...
             'filled', 'MarkerEdgeColor','k', 'Marker', 's' );
         
    %% draw non essential points
    hold on
    scatter( midlife( ordinaries ), persistence( ordinaries ), marker_size, cur_colormap(dims( ordinaries ) + 1,:), ...
             'filled', 'MarkerEdgeColor','k', 'Marker', 'o' );
    hold off
    
    %% draw axis
    min_value = min( midlife );
    max_value = max( persistence );
    safe_min = min_value - 10 * eps( min_value );
    safe_max = max_value + 10 * eps( max_value );
    axis square  
    axis( [safe_min, safe_max, safe_min, safe_max] );
    xticks = get(gca,'XTick');
    set(gca,'YTick',xticks);

    %% draw legend
    legend('Essential','Ordinary','Location','NorthEastOutside');
    
    %% draw axis labels
    xlabel('Midlife: (birth + death) / 2');
    ylabel('Age: death - birth');

    %% draw colobar
    caxis([0 max(dims)])
    ylabel(colorbar('YTick',0:1:max(dims)), 'Dimension');
    
    %% draw title
    set( title(['Data: ' filename]), 'Interpreter', 'none' );
    
    %% draw grid
    grid on;
end