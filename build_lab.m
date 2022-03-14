## lab := data structure containing lines at bus
################################################
function _lab = build_lab(c)
  _lab = Map('KeyType', 'double', 'ValueType', 'any')
  _bi = build_bi(c);
  _nLines = length(c.branch(:,1));
  _numLab = ones(_nLines); # keeps track of how many lines are in b -> []
  
  ## init empty arrays
  ####################
  for b = 1:length(c.bus(:,1))
    _lab(b) = zeros(0);
  end
  
  ## push line ids into arrays
  ############################
  for l = 1:length(c.branch(:,1))
    _current_fbi = _bi(c.branch(l,1));
    _current_tbi = _bi(c.branch(l,2));
    
    _lines = _lab(_current_fbi);
    _lines(_numLab(_current_fbi)) = l;
    _lab(_current_fbi) = _lines;
    _numLab(_current_fbi) += 1;
    
    _lines = _lab(_current_tbi);
    _lines(_numLab(_current_tbi)) = l;
    _lab(_current_tbi) = _lines;
    _numLab(_current_tbi) += 1;
  end
end