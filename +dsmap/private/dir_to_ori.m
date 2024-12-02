function oriList = dir_to_ori(dirList)
oriList = dirList;
for i = 1:length(dirList)
    if oriList(i) >= 180
        oriList(i) = oriList(i) - 180;
    end
end
end