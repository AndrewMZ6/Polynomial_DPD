function rounded_values = round2(input_values, round_power)

    d = 10^round_power;
    a = input_values*d;
    b = round(a);
    rounded_values = b./d;
    

endfunction