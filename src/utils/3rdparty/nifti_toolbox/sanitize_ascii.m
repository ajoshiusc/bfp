% Clean strings of non-printable characters or characters > 127
% characters > 127 get converted to utf-16 during fwrite, causing data positions
% to be shifted. This will prevent that.
% authored by David Shattuck
%
% simple test:
%  x=char(rand(1,80)*255)
%  fprintf(1,'%s\n%s\n',x,sanitize_ascii(x))

function str = sanitize_ascii(str)
str((str>0&str<32)|str>=127)=32;