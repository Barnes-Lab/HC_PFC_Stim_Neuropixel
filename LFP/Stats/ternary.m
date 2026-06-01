% Helper for inline ternary (or just replace with if/else)
function out = ternary(cond, a, b)
    if cond, out = a; else out = b; end
end
