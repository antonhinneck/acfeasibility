## bi := data structure containing bus ids
##########################################
function _bi = build_bi(c)
  _bi = Map('KeyType', 'double', 'ValueType', 'any')
  # init empty arrays
  for b = 1:length(c.bus(:,1))
    #_bi = typeinfo(c.bus(b, 1))(1);
    _idx = c.bus(b, 1);
    _bi(_idx) = b;
  end
end