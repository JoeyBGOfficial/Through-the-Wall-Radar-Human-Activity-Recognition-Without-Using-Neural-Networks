%% Function for Creating Contour Plots
% Author: The MathWorks, Inc.
% Date: 1984-2022 (Last modified: 2022).
%
% Information:
% This function generates a contour plot of a 2D matrix Z, with optional specification of
% x- and y-coordinates, contour levels, and line styles. It supports multiple syntaxes to
% create contour lines automatically or at user-specified levels, with optional axes and
% line specification inputs. The function handles plotting in a specified axes object and
% returns a contour matrix and graphics handle for further customization (e.g., labeling
% with CLABEL). It uses internal MATLAB graphics utilities to parse inputs, manage axes
% properties, and create the contour plot.
%
% The function supports automatic contour level selection, custom levels via a vector,
% and additional name-value pair properties for plot customization. It ensures proper axes
% setup, including tight axis limits and 2D view, and applies Neumann-like boundary handling
% implicitly through graphics settings.
%
% Inputs:
%   varargin - Variable input arguments, supporting the following forms:
%              - Z: 2D matrix to plot contours for (x- and y-coordinates derived from indices).
%              - X, Y, Z: Vectors or matrices defining the mesh and values to plot.
%              - Z, N or X, Y, Z, N: Integer N specifies the number of contour lines.
%              - Z, V or X, Y, Z, V: Vector V specifies contour levels.
%              - AX, ...: Axes handle followed by any of the above syntaxes.
%              - ..., LineSpec: Line type and color specification.
%              - ..., Name-Value pairs: Additional contour properties.
%
% Outputs:
%   cout - Contour matrix containing level and coordinate information (for CLABEL).
%   hand - Handle to the contour graphics object.
%
% Dependencies:
%   - matlab.graphics.chart.internal.contourobjHelper: Internal function for parsing arguments.
%   - matlab.graphics.chart.primitive.Contour: Graphics object for contour rendering.
%   - matlab.graphics.internal.InteractionInfoPanel: Utility for displaying interaction info.
%
% Reference:
%   - MATLAB Documentation: Contour Plot.
%
% See also: CONTOUR3, CONTOURF, CLABEL.

%% Function Body
function [cout, hand] = contour(varargin)
    %CONTOUR Contour plot.
    %   CONTOUR(Z) draws a contour plot of matrix Z in the x-y plane, with
    %   the x-coordinates of the vertices corresponding to column indices
    %   of Z and the y-coordinates corresponding to row indices of Z. The
    %   contour levels are chosen automatically.
    %
    %   CONTOUR(X,Y,Z) draws a contour plot of Z using vertices from the
    %   mesh defined by X and Y. X and Y can be vectors or matrices.
    %
    %   CONTOUR(Z,N) and CONTOUR(X,Y,Z,N) draw N contour lines, choosing
    %   the levels automatically.
    %
    %   CONTOUR(Z,V) and CONTOUR(X,Y,Z,V) draw a contour line for each
    %   level specified in vector V.  Use CONTOUR(Z,[v v]) or
    %   CONTOUR(X,Y,Z,[v v]) to draw contours for the single level v.
    %
    %   CONTOUR(AX, ...) plots into the axes AX.
    %
    %   [C,H] = CONTOUR(...) returns contour matrix C and a handle, H, to
    %   a contour object. These can be used as inputs to CLABEL. The
    %   structure of a contour matrix is described in the contour
    %   documentation.
    %
    %   CONTOUR(..., LineSpec) draws the contours using the line type and
    %   color specified by LineSpec (ignoring marker symbols).
    %
    %   To specify additional contour properties, you can follow the
    %   arguments in any of the syntaxes described above with name-value
    %   pairs.
    %
    %   Example:
    %      [c,h] = contour(peaks);
    %      clabel(c,h)
    %
    %   See also CONTOUR3, CONTOURF, CLABEL.
        
    % Copyright 1984-2022 The MathWorks, Inc.
    
    % Determine the number of outputs
    nout = nargout;
    
    [~, cax, args] = parseplotapi(varargin{:},'-mfilename',mfilename);

    narginchk(1, Inf);
    [pvpairs, ~, ~, errmsg, warnmsg] = matlab.graphics.chart.internal.contourobjHelper('parseargs', false, args{:});
    if ~isempty(errmsg)
        error(errmsg);
    end
    if ~isempty(warnmsg)
        warning(warnmsg);
    end
    
    % Prepend pvpairs specific to contour
    pvpairs = [{ ...
        'EdgeColor_I', 'flat', ...
        'FaceColor_I', 'none', ...
        'ShowText_I', 'off', ...
        'Is3D_I', 'off'} ...
        pvpairs];
    
    if isempty(cax) || ishghandle(cax, 'axes')
        showInteractionInfoPanel = isempty(cax) && isempty(get(groot,'CurrentFigure'));
        cax = newplot(cax);
        parax = cax;
        nextPlot = cax.NextPlot;
        if showInteractionInfoPanel
            % Maybe open the Interaction Info Panel
            matlab.graphics.internal.InteractionInfoPanel.maybeShow(cax);
        end            
    else
        parax = cax;
        cax = ancestor(cax, 'axes');
        nextPlot = 'add';
    end

    % Create contour with parent so defacultCreateFcn gets called.
    h = matlab.graphics.chart.primitive.Contour(pvpairs{:}, 'Parent', parax);
    
    if ismember(nextPlot, {'replace','replaceall'})
        view(cax, 2);
        cax.Box = 'on';
        cax.BoxStyle = 'full';
        cax.Layer = 'top';
        grid(cax,'off');
    end
    
    set(cax, 'XLimSpec', 'tight');
    set(cax, 'YLimSpec', 'tight');
    set(cax, 'ZLimSpec', 'tight');
    
    if nout > 0
        cout = h.ContourMatrix;
    end
    if nout > 1
        hand = h;
    end
end
