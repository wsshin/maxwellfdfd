function D = create_spdiag(vec)

D = spdiags(vec(:), 0, numel(vec), numel(vec));

