function new_file_list = fastdepfun(paths)
% new_file_list = fastdepfun(paths)
% paths = same input as you use with depfun

[file_list] = depfun(paths,'-toponly','-quiet');

% Remove builtins (implement this part however you like)
mroot = matlabroot;
file_list = file_list(~strncmp(file_list,mroot,length(mroot)));

% Remove files already inspected (otherwise we get stuck in an infinite loop)
new_file_list = setdiff(file_list,paths);

if ~isempty(new_file_list)
    new_file_list = fastdepfun(new_file_list);
end
new_file_list = unique([file_list; new_file_list]);

