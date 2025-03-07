;;; SF_STR - The way STRING() should behave by default
function sf_str, stringin, format=format
return, strcompress(string(stringin, format=format), /rem)
end
