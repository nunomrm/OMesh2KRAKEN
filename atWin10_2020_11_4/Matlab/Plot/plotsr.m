function plotsr( Pos )

% Plot the source and receiver positions
% usage: plotsr( Pos )
% where Pos is the standard structure containing positions
% e.g. plotray( Pos )
%
% MBP Nov. 2008


% set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down
% 
% xlabel( 'Range (km)' )
% xlabel( 'Range (m)' )
% ylabel( 'Depth (m)' )
% title( TITLE )
%figure

hold on

Pos.Nrd = length( Pos.r.depth );
Pos.Nrr = length( Pos.r.range );

plot( 1000 * Pos.r.range(  1  ) * ones( Pos.Nsd, 1 ), Pos.s.depth, 'ko', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
        'MarkerSize', 20 )    % source depths
plot( 1000 * Pos.r.range( end ) * ones( Pos.Nrd, 1 ), Pos.r.depth, 'k<', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
        'MarkerSize', 5 )    % receiver depths
plot( 1000 * Pos.r.range, Pos.r.depth( 1 ) * ones( Pos.Nrr, 1 ),   'kv', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
        'MarkerSize', 5 )    % receiver ranges

hold off
zoom on
