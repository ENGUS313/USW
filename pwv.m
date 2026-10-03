function out = pwv(image, peaks, x_m, z_m, fs)
% image: [depth,lateral,time] intensity; peaks: two [row,column] wall points.
% x_m/z_m: image coordinates in metres; fs: measured image frame rate (Hz).
% One continuous segment; requires Signal Processing Toolbox. PWV is signed m/s.
[nz, nx, nt] = size(image);
assert(isreal(image) && nt > 13 && nz >= 31, 'Invalid intensity image stack.');
assert(isequal(size(peaks), [2 2]) && all(isfinite(peaks(:))) && ...
    all(peaks(:) == round(peaks(:))), 'Use two integer [row,column] points.');
assert(all(peaks(:,1) >= 16 & peaks(:,1) <= nz-15) && ...
    all(peaks(:,2) >= 1 & peaks(:,2) <= nx), 'Wall points are outside valid bounds.');
assert(numel(x_m) == nx && numel(z_m) == nz && all(isfinite(x_m(:))) && ...
    all(isfinite(z_m(:))) && all(diff(z_m(:)) > 0), 'Invalid image coordinates.');
assert(isscalar(fs) && isfinite(fs) && fs > 40, 'Frame rate must exceed 40 Hz.');
distance = abs(x_m(peaks(2,2)) - x_m(peaks(1,2)));
assert(distance > 0, 'Wall points must have different lateral positions.');
wall = zeros(nt, 2);
for site = 1:2
    row = peaks(site,1); col = peaks(site,2);
    template = double(image(row-15:row+15,col,1));
    template = template - mean(template);
    assert(all(isfinite(template)) && norm(template) > 0, 'Invalid wall template.');
    template = template / (norm(template) + eps);
    for t = 1:nt
        scan = max(double(image(:,max(1,col-1):min(nx,col+1),t)), [], 2);
        assert(all(isfinite(scan)), 'Nonfinite wall samples.');
        candidates = max(16,round(row)-30):min(nz-15,round(row)+30);
        score = zeros(size(candidates));
        for j = 1:numel(candidates)
            v = scan(candidates(j)-15:candidates(j)+15);
            v = v - mean(v);
            score(j) = sum(template .* (v / (norm(v) + eps)));
        end
        [~, j] = max(score);
        row = candidates(j);
        if j > 1 && j < numel(score)
            p = polyfit(candidates(j-1:j+1), score(j-1:j+1), 2);
            if p(1) < -1e-8 && all(isfinite(p))
                delta = -p(2)/(2*p(1)) - row;
                if abs(delta) < 1, row = row + delta; end
            end
        end
        wall(t,site) = row;
    end
end
out = estimate_pwv(wall, z_m, distance, fs);
out.wall_px = wall;
end

function out = estimate_pwv(wall, z_m, distance, fs)
trace = interp1((1:numel(z_m))', z_m(:), wall, 'linear');
[b, a] = butter(4, 20/(fs/2), 'low');
trace = filtfilt(b, a, trace);
trace = trace - mean(trace, 1);
assert(all(vecnorm(trace) > 1e-12), 'No measurable wall motion.');
out.inverted = sum(trace(:,1).*trace(:,2)) < 0;
if out.inverted, trace(:,2) = -trace(:,2); end
upsample = 10;
t = (0:size(trace,1)-1)' / fs;
t_up = (0:(size(trace,1)-1)*upsample)' / (fs*upsample);
v = interp1(t, trace, t_up, 'spline');
v = v - mean(v, 1);
[r, lag] = xcorr(v(:,2), v(:,1), 'coeff');
[out.correlation, k] = max(r);
out.delay_s = lag(k) / (fs*upsample);
out.distance_m = distance;
out.pwv_m_s = NaN;
if out.delay_s ~= 0, out.pwv_m_s = distance / out.delay_s; end
out.wall_m = trace;
end
