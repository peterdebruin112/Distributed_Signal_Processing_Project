function n = find_n(r)
    fun = @(n) r - sqrt(2*log(n)/n);
    options = optimoptions('fsolve', 'Display', 'off');
    n = fsolve(fun, 10, options);
end