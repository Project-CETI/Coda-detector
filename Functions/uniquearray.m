function [outArray, oldInd, newInd, repCell] = uniquearray(...
  inArray, dim, eqnan, wb, gc)
% [B, ai, bi, rep] = UNIQUEARRAY(A, dim, eqnan, wb, gc)
% Finds the unique parts in an array that can be of any type, including cell
% arrays and structures with mixed data, and have any number of dimensions.
%
% INPUT:
% A     - an array. For example, a string, vector, matrix, structure array or a
%         cell array.
% dim   - optional integer, the considered dimension. For a matrix dim = 1 gives
%         the unique rows and dim = 2 gives the unique columns. Default is the
%         first non-singleton dimension.
% eqnan - optional logical (0/1). If 0 (default) NaN's are not considered equal,
%         if 1 they are considered equal.
% wb    - optional logical (0/1). Adds a wait bar if 1, default is 0.
% gc    - optional logical (0/1). Default (0) is to uses MATLAB's faster unique
%         function for supported cases. If 1 always uses the generic code (only
%         for testing!).
%
% OUTPUT:
% B   - an array with the unique parts in the considered dimension. The
%       parts will occur in the same order as in A.
% ai  - an index vector such that A(ai) is B.
% bi  - an index vector such that B(bi) is A. 
% rep - a cell array with three columns. The first column contains an index in
%       A to the first occurrence of a part, the second column how many time
%       this part is repeated in A and the third column indices in A to
%       repeated occurrences of the part (if any). 
% All indices above are for the considered dimension.
%
% NOTE:
% Can be slow for large arrays.
%
% EXAMPLES:
% clear
%
% unique([1 NaN NaN 2], 'stable')
% UNIQUEARRAY([1 NaN NaN 2], [], 1)  % here we can consider NaN's equal
%
% tmp = rand(10); A(:, :, 1) = tmp; A(:, :, 2) = tmp'; A(:, :, 3) = tmp;
% UNIQUEARRAY(A, 3)  % gives the unique sheets in A
%
% s(1).test = 2; s(2).test = 2; s(3).test = 'a'
% unique(s, 'stable')  % try to get the unique elements in s...
% UNIQUEARRAY(s)
%
% str = {'a', 'b', 'c'; 'a', 'b', 'c'}
% unique(str, 'rows', 'stable')
% UNIQUEARRAY(str)
%
% B = {{}, A, [1 2 3], {}, {1 2}, A; {}, s, 'b', {}, A, s}
% unique(B', 'rows', 'stable') % try to get the unique columns...
% UNIQUEARRAY(B, 2)
%
% see also UNIQUE
 
% VERSION 1.1, 2021
% - Updated slice function after suggestion by Stephen
% - Fixed calls to unique to handle two-dimensional string cells (unique does
%   not work correctly for these) 
%
% By Patrik Forssén (SatStar Ltd & Karlstad University)
 
%----------------------------------INPUT----------------------------------------
% Default output
outArray = inArray;
oldInd   = 0; newInd = 0;
repCell  = {};
% First input
dimSize = size(inArray);
if (isempty(inArray))
  % Empty array
  return
end
if (max(dimSize) == 1)
  % Array with just one item
  oldInd = 1; newInd = 1; repCell = {1, 0, []};
  return
end
 
% Secound input (optional)
if (nargin < 2 || isempty(dim))
  dim = find(dimSize > 1, 1, 'first');
end
% Check
nDim = ndims(inArray);
if (~isnumeric(dim) || ~isscalar(dim) || ~isreal(dim) || ...
    dim - round(dim) ~= 0)
  error('uniquearray:WrongType', ['Second input (dimension) must be an ', ...
    'integer scalar'])
end
if (dim < 0 || dim > nDim)
  error('uniquearray:WrongValue', ['Second input (dimension) must be ', ...
    'between 1 and ', num2str(nDim)])
end
if (dimSize(dim) == 1)
  % Requested singleton dimension
  oldInd = 1; newInd = 1; repCell = {1, 0, []};
  return
end
 
% Third input (optional)
chkLogical = @(x) ~isscalar(x) || (~islogical(x) && (~isnumeric(x) || ...
  ~isreal(x) || x - round(x) ~= 0 || x < 0 || x > 1));
if (nargin < 3 || isempty(eqnan))
  % Default
  eqnan = false;
end
if (chkLogical(eqnan))
  error('uniquearray:WrongLogical', ...
    'Third input (treat NaN''s as equal) must be an scalar logical or 0/1')
end
% Fourth input (optional)
if (nargin < 4 || isempty(wb))
  % Default
  wb = false;
end
if (chkLogical(wb))
  error('uniquearray:WrongLogical', ...
    'Fourth input (waitbar) must be an scalar logical or 0/1')
end
% Fifth input (optional)
if (nargin < 5 || isempty(gc))
  % Default
  gc = false;
end
if (chkLogical(gc))
  error('uniquearray:WrongLogical', ...
    'Fifth input (use generic code) must be an scalar logical or 0/1')
end
 
%--------------------------------CALCULATION------------------------------------
if (~gc && ~eqnan && nDim <= 2)
  % Try "unique", much faster when it works!
  nsDim = sum(dimSize > 1);  % Number of non-singleton dimensions
  if (dim == 1)
    try
      if (nsDim == 1)
        [outArray, oldInd, newInd] = unique(inArray, 'stable');
        unFlag = 1;
      else
        warning off
        [outArray, oldInd, newInd] = unique(inArray, 'stable', 'rows');
        warning on
        if (size(outArray, 2) ~= dimSize(2))
          % Did not work!
          outArray = inArray;
          unFlag   = 0;
        else
          unFlag   = 1;
        end
      end
    catch
      % Failed, try generic code
      unFlag = 0;
    end
  else
    try
      % Here we can try transpose (might not be defined for all classes!)
      warning off
      [outArray, oldInd, newInd] = unique(transpose(inArray), ...
        'stable', 'rows');
      warning on
      if (size(outArray, 2) ~= dimSize(1))
        % Did not work!
        outArray = inArray;
        unFlag   = 0;
      else
        outArray = transpose(outArray);
        unFlag = 1;
      end
    catch
      % Failed, try generic code
      unFlag = 0;
    end
  end
else
  % Have to use the generic code!
  unFlag = 0;
end
 
 
% Genric code
if (~unFlag)
  % Index vector
  indVec = (1 : size(inArray, dim))';
  oldInd = indVec;
  newInd = indVec;
  
  % Comparision function
  if (~eqnan)
    % NaN's are not equal
    compFun = @(x, y) isequal(x, y);
  else
    % NaN's are equal
    compFun = @(x, y) isequaln(x, y);
  end
  
  % Iterate
  refInd = 0;
  if (wb), wh = waitbar(0, 'Calculating unique parts...'); end
  while (1)
    refInd  = refInd + 1;
    nSlices = size(outArray, dim);
    if (wb), waitbar(refInd/nSlices, wh), end
    if (refInd > nSlices), break, end
    newInd(oldInd(refInd)) = refInd;
    
    refSlice = getSlice(outArray, dim, refInd);
    rmSlice  = [];
    for compInd = refInd+1 : nSlices
      compSlice = getSlice(outArray, dim, compInd);
      if (compFun(refSlice, compSlice))
        % Remove this slice
        rmSlice = [rmSlice, compInd];
        newInd(oldInd(compInd)) = refInd;
      end
    end
    keepSlice = setdiff(1:nSlices, rmSlice);
    outArray  = getSlice(outArray, dim, keepSlice);
    oldInd    = oldInd(keepSlice);
  end
  if (wb), delete(wh), end
end
 
  
% Make a cell array with repetitions
repCell = cell(length(oldInd), 3);
for itemNo = 1 : length(oldInd)
  rowInd = find(newInd == itemNo);
  repCell(itemNo, :) = {rowInd(1), length(rowInd)-1, rowInd(2:end)};
end
 
end
 
 
 
function slice = getSlice(A, dim, sliceInd)
% Get slice(s) of A in a specified dimension
% As suggested by Stephen
C      = repmat({':'}, 1, max(ndims(A), dim)); 
C{dim} = sliceInd; 
slice  = A(C{:});
 
end