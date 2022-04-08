# Enumerate permuations of single lines 
# to reconnect isolated buses.
function recursiveEnumPoints(path, isolated, lab)
  global points;
  if length(path) == 0
    _clab = lab(isolated(1));
    for l = _clab
      recursiveEnumPoints([l], isolated, lab);
    end
  elseif (length(path) + 1) != length(isolated)
    _lvl = length(path) + 1;
    _clab = lab(isolated(_lvl));    
    for l = _clab
      _npath = path;
      _npath(_lvl) = l;
      recursiveEnumPoints(_npath, isolated, lab);
    end
  else
    _lvl = length(path) + 1;
    _clab = lab(isolated(_lvl));
    for l = _clab
      path = [path, l];     
      _np = length(points(1,:));
      points(:, _np + 1) = path;
    end
  end
end