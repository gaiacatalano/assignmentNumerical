function varargout = experiment_utils(action, varargin)
% Central utility dispatcher for experiments

    switch action
        case 'make_row'
            varargout{1} = make_row(varargin{:});
    
        case 'print_table'
            print_table_like_pdf(varargin{:});
    
        case 'update_counts'
            [varargout{1}, varargout{2}] = update_counts(varargin{:});

         % ===== NEW: NELDER-MEAD =====
        case 'make_row_nm'
            varargout{1} = make_row_nm(varargin{:});

        case 'print_table_nelder'
            print_table_nelder(varargin{:});
    
    end
end

function r = make_row(startLabel, methodName, it, xnorm, fval, gnorm, timeSec)
    r.start  = char(startLabel);
    r.method = char(methodName);
    r.iter   = it;
    r.xnorm  = xnorm;
    r.fval   = fval;
    r.gnorm  = gnorm;
    r.time   = timeSec;
end

function print_table_like_pdf(fid, n, hstep, rows)
    if nargin < 1 || isempty(fid), fid = 1; end

    fprintf(fid, '\n');
    fprintf(fid, ...
        'Table - Performance of Modified Newton variants (n = %d, h = %.1e)\n', ...
        n, hstep);

    wStart  = 14; wMethod = 26; wIter = 10;
    wX = 14; wF = 14; wG = 10; wT = 10;

    lineLen = wStart + wMethod + wIter + wX + wF + wG + wT + 6*3;
    fprintf(fid, '%s\n', repmat('-',1,lineLen));

    fprintf(fid, ['%-' num2str(wStart) 's | ' ...
                  '%-' num2str(wMethod) 's | ' ...
                  '%'  num2str(wIter)   's | ' ...
                  '%'  num2str(wX)      's | ' ...
                  '%'  num2str(wF)      's | ' ...
                  '%'  num2str(wG)      's | ' ...
                  '%'  num2str(wT)      's\n'], ...
            'Starting point','Method','Iterations','||x*||_2','Final f(x)','||∇f||','Time(s)');

    fprintf(fid, '%s\n', repmat('-',1,lineLen));

    prevStart = '';
    for i = 1:numel(rows)
        if ~strcmp(prevStart, rows(i).start) && ~isempty(prevStart)
            fprintf(fid, '%s\n', repmat('-',1,lineLen));
        end
        prevStart = rows(i).start;

        fprintf(fid, ['%-' num2str(wStart) 's | ' ...
                      '%-' num2str(wMethod) 's | ' ...
                      '%'  num2str(wIter)   'd | ' ...
                      '%'  num2str(wX)      '.4e | ' ...
                      '%'  num2str(wF)      '.4e | ' ...
                      '%'  num2str(wG)      '.2e | ' ...
                      '%'  num2str(wT)      '.4f\n'], ...
                rows(i).start, rows(i).method, rows(i).iter, ...
                rows(i).xnorm, rows(i).fval, rows(i).gnorm, rows(i).time);
    end

    fprintf(fid, '%s\n', repmat('-',1,lineLen));
end

function [succ, fail] = update_counts(gnorm, epsilon, succ, fail)
    if gnorm < epsilon
        succ = succ + 1;
    else
        fail = fail + 1;
    end
end


function r = make_row_nm(startLabel, methodName, it, xnorm, fval, D, timeSec)
    r.start  = char(startLabel);
    r.method = char(methodName);
    r.iter   = it;        % NaN if not available
    r.xnorm  = xnorm;
    r.fval   = fval;
    r.D      = D;
    r.time   = timeSec;
end


function print_table_nelder(fid, n, rows)
    if nargin < 1 || isempty(fid), fid = 1; end

    fprintf(fid, '\n');
    fprintf(fid, 'Table - Nelder-Mead results (n = %d)\n', n);

    wStart  = 14;
    wMethod = 14;
    wIter   = 10;
    wX      = 14;
    wF      = 14;
    wD      = 14;
    wT      = 10;

    lineLen = wStart + wMethod + wIter + wX + wF + wD + wT + 6*3;
    fprintf(fid, '%s\n', repmat('-',1,lineLen));

    fprintf(fid, ['%-' num2str(wStart) 's | ' ...
                  '%-' num2str(wMethod) 's | ' ...
                  '%'  num2str(wIter)   's | ' ...
                  '%'  num2str(wX)      's | ' ...
                  '%'  num2str(wF)      's | ' ...
                  '%'  num2str(wD)      's | ' ...
                  '%'  num2str(wT)      's\n'], ...
            'Starting point','Method','Iterations','||x*||_2','Final f(x)','D(simplex)','Time(s)');

    fprintf(fid, '%s\n', repmat('-',1,lineLen));

    prevStart = '';
    for i = 1:numel(rows)

        if ~strcmp(prevStart, rows(i).start) && ~isempty(prevStart)
            fprintf(fid, '%s\n', repmat('-',1,lineLen));
        end
        prevStart = rows(i).start;

        % Iterations may be NaN -> print '-'
        if isnan(rows(i).iter)
            iterStr = '-';
            fprintf(fid, ['%-' num2str(wStart) 's | ' ...
                          '%-' num2str(wMethod) 's | ' ...
                          '%'  num2str(wIter)   's | ' ...
                          '%'  num2str(wX)      '.4e | ' ...
                          '%'  num2str(wF)      '.4e | ' ...
                          '%'  num2str(wD)      '.4e | ' ...
                          '%'  num2str(wT)      '.4f\n'], ...
                    rows(i).start, rows(i).method, iterStr, ...
                    rows(i).xnorm, rows(i).fval, rows(i).D, rows(i).time);
        else
            fprintf(fid, ['%-' num2str(wStart) 's | ' ...
                          '%-' num2str(wMethod) 's | ' ...
                          '%'  num2str(wIter)   'd | ' ...
                          '%'  num2str(wX)      '.4e | ' ...
                          '%'  num2str(wF)      '.4e | ' ...
                          '%'  num2str(wD)      '.4e | ' ...
                          '%'  num2str(wT)      '.4f\n'], ...
                    rows(i).start, rows(i).method, rows(i).iter, ...
                    rows(i).xnorm, rows(i).fval, rows(i).D, rows(i).time);
        end
    end

    fprintf(fid, '%s\n', repmat('-',1,lineLen));
end




