function converted = units_conversion (to_convert, option)

% units_conversion performs the conversion from one unit to another.
% Further development has to be done.
%
% OUTPUT:
% - converted     [double]    Converted measure
%
% INPUT:
% - to_convert    [double]    Value to be converted
% - option        [char]      Option specifying the conversion
%
% - 'lbs_to_kg'   Pounds to kilograms
% - 'kg_to_lbs'   Kilograms to pounds
%
% - 'lbs_to_newton'   Pounds to Newton
% - 'newton_to_lbs'   Newton to pounds
% 
% - 'in_to_m'     Inches to metres
% - 'm_to_in'     Metres to inches

switch option
    case 'lbs_to_kg'
        converted = 0.45359237 * to_convert;

    case 'kg_to_lbs'
        converted = 2.2046226218 * to_convert;

    case 'lbs_to_newton'
        converted = 4.4482216153 * to_convert;

    case 'newton_to_lbs'
        converted = 0.2248089431 * to_convert;

    case 'in_to_m'
        converted = 0.0254 * to_convert;
    case 'm_to_in'
        converted = 39.3700787402 * to_convert;
end

end