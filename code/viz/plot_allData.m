function [] = plot_allData( data, ii )

T = data.Ts(ii);
M = 10;

xs = 1:T;
minVal=min(min(data.Xdata,[],2)); % Compute value for the rows
maxVal=max(max(data.Xdata,[],2));
ys = linspace( minVal, maxVal, M);

hold all;

% Plot colors with xs as the x-axis time series, ys determines the
% top&bottom levels.
hIM = imagesc( xs, ys, repmat(data.zTrue(ii), M, 1), [1 max( data.zTrueAll)] ); % Indicate which states correspond to which timesteps.
set( hIM, 'AlphaData', 0.65 );

% seq: function which
X = data.seq(ii);

for i=1:data.D
    plot( xs, X(i,:), 'color',[(rand(1,3))],'LineStyle','-.' );
    title( ['Sequence ' num2str(ii)], 'FontSize', 20 );
end
axis( [1 T ys(1) ys(end)] );

end