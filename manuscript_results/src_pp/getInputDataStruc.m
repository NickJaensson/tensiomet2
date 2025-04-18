function [InputDataStruc] = getInputDataStruc(DataStruc, InputVarList)

for i = 1:size(InputVarList, 2)
    InputDataStruc.(InputVarList{1,i}) = DataStruc.(InputVarList{1,i});
end

end