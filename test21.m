percent = 12.1;
disp(['testing'])
w=50;
k=10;
perc = sprintf('%3.1f%%', percent); % 4 characters wide, percentage
disp([repmat(char(8), 1, (w+k)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
disp([repmat(char(8), 1, (w+k)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
disp([repmat(char(8), 1, (w+k)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
disp([repmat(char(8), 1, (w+k)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
