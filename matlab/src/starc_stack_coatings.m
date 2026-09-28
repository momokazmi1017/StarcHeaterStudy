function coat = starc_stack_coatings(p, names)
%STARC_STACK_COATINGS  Combine named coatings into one struct of column vectors.
%
%   coat = starc_stack_coatings(p, {'FEP','OSR','STARC'})
%
%   The result can be passed to STARC_CONFIG so all coatings are simulated
%   together in one vectorised run.

names = cellstr(names);
fields = {'alpha', 'eps_c', 'eps_h', 'Tt_C', 'w'};
for k = 1:numel(fields)
    f = fields{k};
    coat.(f) = cellfun(@(n) p.coat.(n).(f), names(:));
end
coat.name = cellfun(@(n) p.coat.(n).name, names(:), 'UniformOutput', false);
coat.key  = names(:);
end
