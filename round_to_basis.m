function [ basis_indicies ] = round_to_basis( A, indicies_order )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [m,n] = size(A);
    basis_indicies = [];
    current_rank = 0;
    for i = 1:length(indicies_order)
        col_index = indicies_order(i);
        candidate_basis = [basis_indicies; col_index];
        candidate_rank = spnrank(A(:,candidate_basis));
        if candidate_rank > current_rank
            basis_indicies = candidate_basis;
            current_rank = candidate_rank;
        end
        if current_rank == m
           break
        end
    end
end

