function srcdir = get_data_dir(loc)

cutoff=1;
srcdir='';
if ispc
    tag='\';
else
    tag='/';
end
for n=1:numel(loc)-1

    if strcmpi(loc(end-n),tag)
        cutoff=numel(loc)-n;
        break;
    end
end

srcdir=loc(1:cutoff);

% END
end