%-------------------------------------------------------------------------------
% Updright  : SPACE Research Centre, RMIT
% Project   : SPACE
% Author    : Changyong HE
% Language  : matlab
% CreateTime: 23 Nov 2017
% Function  : This function is used to make density plot
%-------------------------------------------------------------------------------

function H = densityplot(x,y,varargin)
% Parameters as 'name',value pairs:
% - 'nbins': Array in the form of [nxbins nybins] to set the
% number of bins in each dimension
% - 'edges': Cell in the form of {[x__edges] [y_edges]} to set
% custom bin edges for each dimension
%
% More details refer to hist3
if isrow(x)
    x = x(:);
end
if isrow(y)
    y = y(:);
end
% combine
X = [x,y];

% sumarisation
% N is a matrix containing the number of elements of X that fall in 
%   each bin of the grid.
% C returns the positions of the bin centers in a 1-by-2 cell array
%   of numeric vectors.
[N,C] = hist3(X, varargin{:});

if ~isempty(find(strcmpi(varargin,'edges'), 1))
  %Put values on the upper edges as if they were in the last
  %bin
  N(:,end-1)=N(:,end-1)+N(:,end);
  N(end-1,:)=N(end-1,:)+N(end,:);
  %Remove upper edge
  N(:,end)=[];
  N(end,:)=[];
  C{1}(end) = [];
  C{2}(end) = [];
end

%Get polygon half widths
wx=C{1}(:);
wy=C{2}(:);

% display
% modification
myColorMap=parula(64); %choose the colormap scheme 
R = linspace(1, myColorMap(1, 1), 15);  %// Red from 1 to 0
G = linspace(1, myColorMap(1, 2), 15);  %// Blue from 0 to 1
B = linspace(1, myColorMap(1, 3), 15);   %// Green all zero
myColorMap = [R', G', B'; myColorMap];

figure
H = pcolor(wx, wy, N');
box on
grid off
shading interp
set(H,'edgecolor', 'none');
colorbar
colormap(myColorMap)
set(gca, 'fontsize', 20, 'FontWeight','bold','linewidth',4)

%colormap jet
