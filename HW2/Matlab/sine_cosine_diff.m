function [sine_value, cosine_value] = sine_cosine_diff (sine_i, sine_j, cosine_i, cosine_j) 
    sine_value = sine_i*cosine_j -cosine_i*sine_j;
    cosine_value = cosine_i*cosine_j +sine_i*sine_j;
end