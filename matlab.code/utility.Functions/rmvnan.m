function Xno_nan= rmvnan(x)
% F: removes nan ROWS.

I = anynan(x);

if sum(I)==0;
	Xno_nan = x;
else
	Xno_nan = x(~I,:);
end;

