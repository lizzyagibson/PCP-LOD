function flag = is_same(SAME_THRESH, varargin)

flag = true;
for i = 1 : nargin-1
    for j = i+1 : nargin-1
        if max(max(abs(varargin{i}-varargin{j}))) > SAME_THRESH
            flag = false;
            return
        end
    end
end
