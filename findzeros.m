function sol = findzeros(f,range,err)
if nargin < 3
  err = 1e4*eps;
end
if range(1)>=range(2)
    sol = [];
    return
end
% digits(4);
sol = double(vpasolve(f,symvar(f,1),range));

if(isempty(sol))
  return
else
  lowLimit = sol-err;
  highLimit = sol+err;
  temp = findzeros(f,[range(1) lowLimit],err);
  if ~isempty(temp)
    sol = sort([sol temp]);
  end
  temp = findzeros(f,[highLimit range(2)],err);
  if ~isempty(temp)
    sol = sort([sol temp]);
  end
  return
end
end