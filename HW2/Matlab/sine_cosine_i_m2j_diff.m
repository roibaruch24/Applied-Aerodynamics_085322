function [sin_i_m_2j, cos_i_m_2j] = sine_cosine_i_m2j_diff(sine_i, sine_j,cosine_i, cosine_j)
    sin_2j = 2 * sine_j .* cosine_j;
    cos_2j = cosine_j.^2 - sine_j.^2;
    sin_i_m_2j = sine_i .* cos_2j - cosine_i .* sin_2j;
    cos_i_m_2j = cosine_i .* cos_2j + sine_i .* sin_2j;
end