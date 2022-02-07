% Jan. 31, 2017, Yunfeng Chen, normalize the receiver function to unit
% amplitude
function seisout = normalize_RFs(seis)
% input: seis cell array contains RFs
% output: normalized RFs
seisout = cell(size(seis));
for n = 1:length(seis)
   x = seis{n};
   xout = x/max(x);
   seisout{n} = xout;
end