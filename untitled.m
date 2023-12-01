%% 从文本文件中导入数据
% 用于从以下文本文件中导入数据的脚本:
%
%    filename: E:\work\work_project\Document\uhf_data\data_00.txt
%
% 由 MATLAB 于 2023-10-11 09:10:15 自动生成

%% 设置导入选项并导入数据
clear ;
opts = delimitedTextImportOptions("NumVariables", 1);

% 指定范围和分隔符
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% 指定列名称和类型
opts.VariableNames = "VarName1";
opts.VariableTypes = "double";

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 导入数据
data00 = readtable("200.txt", opts);
% 
x = table2array(data00);
x = transpose(x);
plot(x);
fileID = fopen("raw_data.txt","w");
writematrix(x,"raw_data.txt",'delimiter',' ');

fclose(fileID);

%% 清除临时变量
clear opts