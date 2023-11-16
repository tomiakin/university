function indexes = Rep_cleaner(input)
% input   - vector containing the elements that are scrutinised for repetitions
% indexes - rows to keep in the original catalog

[~,indexes,~]=unique(input,'rows');
end