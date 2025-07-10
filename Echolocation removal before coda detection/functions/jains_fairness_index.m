function fairness = jains_fairness_index(x)
    % Computes Jain's Fairness Index for a given allocation vector x
    % Input: x - a vector of resource allocations (non-negative values)
    % Output: fairness - Jain's Fairness Index (value between 0 and 1)
    
    if nargin == 0 || isempty(x)
        error('Input vector x cannot be empty.');
    end
    
    if any(x < 0)
        error('All elements of x must be non-negative.');
    end

    numerator = sum(x)^2;
    denominator = length(x) * sum(x.^2);
    
    if denominator == 0
        fairness = 0;  % Avoid division by zero; fairness is 0 if all allocations are zero.
    else
        fairness = numerator / denominator;
    end
end