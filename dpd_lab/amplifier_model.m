function amplified_values = amplifier_model(input_values)

    # Polynomial amplifier model coefficients
    a0 = 1086.8;
    a1 = -1320.6;
    a2 = 509.41;
    a3 = -77.974;
    a4 = 15.419;
    #a5 = 0.023673;
    a5 = 0.0;


    # Truncate input values
    input_values(input_values >= 0.3832) = 0.3832;
    input_values(input_values <= -0.3832) = -0.3832;


    # Array of signs (+1 or -1). Since the polynomial affects positive and negative input values differently
    # signs of the values are stored separately, whilst values modules are sent to model. For loop is made
    # because of zero input values. Without it NaNs are created

    signs = zeros(1, length(input_values));
    for i = 1:length(input_values)
        if (!input_values(i))
            signs(i) = 1;
        else
            signs(i) = input_values(i)./abs(input_values(i));
        endif
    endfor


    # Уравнение полинома через который проходит входной сигнал. За подробностями откуда взялись цифры https://github.com/AndrewMZ6/Job/tree/main/DSP/DPD
    amplifier_function = @(x) a5 + a4*x + a3*(x.^2) + a2*(x.^3) + a1*(x.^4) + a0*(x.^5);

    # Применение полинома
    amplified_values_abs = amplifier_function(abs(input_values));

    # Расстановка знаков "по местам"
    amplified_values = amplified_values_abs.*signs;

endfunction
